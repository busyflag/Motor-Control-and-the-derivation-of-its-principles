#include <Stdlib.h>
#include <String.h>
#include <math.h>
#include <Psim.h>


#define PI 3.14159265358979323846
#define SQRT3 1.7320508075688772

/************** 电机参数 **************/
static double Rs = 0.294;          // 定子电阻，单位 Ohm，需要按实际电机修改
static double pole_pairs = 3.0;    // 极对数，例如 4 极电机 p=2

/************** DTC 控制参数 **************/
static double Flux_ref = 0.9;      // 定子磁链给定，单位 Wb，需要按电机调试
static double Flux_band = 0.01;    // 磁链滞环宽度
static double Te_band = 0.5;       // 转矩滞环宽度，单位 N.m

/************** 转速 PI 参数 **************/
static double Kp_speed = 200.0;
static double Ki_speed = 1;
static double Te_max = 200.0;       // 最大转矩给定
static double Te_min = -200.0;      // 最小转矩给定

/************** 状态变量 **************/
static double Fsalfa = 0.0;
static double Fsbeta = 0.0;

static double speed_int = 0.0;

static int HF_last = 1;

static int Sa = 0;
static int Sb = 0;
static int Sc = 0;


/************** 符号函数：大于等于 0 返回 1，小于 0 返回 0 **************/
int i_sign(double x)
{
    if (x >= 0.0)
        return 1;
    else
        return 0;
}


/************** 限幅函数 **************/
double limit_value(double x, double max, double min)
{
    if (x > max)
        return max;
    if (x < min)
        return min;
    return x;
}


/************** 扇区判断函数 **************/
int GetSector(double Falfa, double Fbeta)
{
    double A, B, C;
    int N;
    int Sector;

    int section[8] = {5, -1, 4, 3, 6, 1, -1, 2};

    A = Falfa;
    B = SQRT3 * Fbeta - Falfa;
    C = SQRT3 * Fbeta + Falfa;

    N = 4 * i_sign(A) + 2 * i_sign(B) + i_sign(C);

    Sector = section[N];

    /*
     * 理论上 Sector 不应该等于 -1。
     * 但仿真初始磁链很小时，可能出现异常。
     * 这里做保护。
     */
    if (Sector < 1 || Sector > 6)
    {
        double angle;

        angle = atan2(Fbeta, Falfa);

        if (angle < 0.0)
            angle += 2.0 * PI;

        Sector = (int)(angle / (PI / 3.0)) + 1;

        if (Sector < 1)
            Sector = 1;
        if (Sector > 6)
            Sector = 6;
    }

    return Sector;
}


/************** 根据电压矢量编号得到 Sa Sb Sc **************/
void VectorToSwitch(int Vector, int *pSa, int *pSb, int *pSc)
{
    /*
     * u0 = 000
     * u1 = 100
     * u2 = 110
     * u3 = 010
     * u4 = 011
     * u5 = 001
     * u6 = 101
     * u7 = 111
     */

    int SaTable[8] = {0, 1, 1, 0, 0, 0, 1, 1};
    int SbTable[8] = {0, 0, 1, 1, 1, 0, 0, 1};
    int ScTable[8] = {0, 0, 0, 0, 1, 1, 1, 1};

    if (Vector < 0)
        Vector = 0;

    if (Vector > 7)
        Vector = 7;

    *pSa = SaTable[Vector];
    *pSb = SbTable[Vector];
    *pSc = ScTable[Vector];
}


/************** PSIM 主仿真函数 **************/
void SimulationStep(
        double t, double delt, double *in, double *out,
        int *pnError, char * szErrorMsg,
        void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    /******************** 1. 读取输入 ********************/
    double speed_ref_rpm = 600;
	double delt=0.0001;
    double speed_rpm     = in[0];

    double ia = in[1];
    double ib = in[2];
    double ic = in[3];

    double Udc = 600;

    /******************** 2. 三相电流 Clarke 变换 ********************/
    double isalfa, isbeta;

    isalfa = (2.0 * ia - ib - ic) / 3.0;
    isbeta = SQRT3 * (ib - ic) / 3.0;

    /******************** 3. 根据上一拍开关状态计算定子电压 ********************/
    double usalfa, usbeta;

    usalfa = Udc * (2.0 * Sa - Sb - Sc) / 3.0;
    usbeta = Udc * SQRT3 * (Sb - Sc) / 3.0;

    /******************** 4. 定子磁链观测器 ********************/
    Fsalfa = Fsalfa + (usalfa - Rs * isalfa) * delt;
    Fsbeta = Fsbeta + (usbeta - Rs * isbeta) * delt;

    double Flux;

    Flux = sqrt(Fsalfa * Fsalfa + Fsbeta * Fsbeta);

    /******************** 5. 电磁转矩计算 ********************/
    double Te;

    Te = 1.5 * pole_pairs * (Fsalfa * isbeta - Fsbeta * isalfa);

    /******************** 6. 转速外环 PI，输出转矩给定 ********************/
    double speed_ref_rad;
    double speed_rad;
    double speed_err;
    double Te_ref_unsat;
    double Te_ref;

    speed_ref_rad = speed_ref_rpm * 2.0 * PI / 60.0;
    speed_rad     = speed_rpm     * 2.0 * PI / 60.0;

    speed_err = speed_ref_rad - speed_rad;

    speed_int = speed_int + speed_err * delt;

    Te_ref_unsat = Kp_speed * speed_err + Ki_speed * speed_int;

    Te_ref = limit_value(Te_ref_unsat, Te_max, Te_min);

    /*
     * 简单抗积分饱和：
     * 如果 PI 输出饱和，则回退本步积分。
     */
    if (Te_ref != Te_ref_unsat)
    {
        speed_int = speed_int - speed_err * delt;
    }

    /******************** 7. 磁链滞环比较器 ********************/
    int HF;

    if (Flux < Flux_ref - Flux_band)
    {
        HF = 1;       // 磁链偏小，需要增加磁链
    }
    else if (Flux > Flux_ref + Flux_band)
    {
        HF = -1;      // 磁链偏大，需要减小磁链
    }
    else
    {
        HF = HF_last; // 滞环内部保持
    }

    HF_last = HF;

    /******************** 8. 转矩滞环比较器 ********************/
    int HTe;
    double Te_err;

    Te_err = Te_ref - Te;

    if (Te_err > Te_band)
    {
        HTe = 1;       // 转矩偏小，需要增加转矩
    }
    else if (Te_err < -Te_band)
    {
        HTe = -1;      // 转矩偏大，需要减小转矩
    }
    else
    {
        HTe = 0;       // 转矩基本合适
    }

    /******************** 9. 扇区判断 ********************/
    int Sector;

    Sector = GetSector(Fsalfa, Fsbeta);

    /******************** 10. 空间电压矢量选择表 ********************/
    int DTCSVTable[6][6] =
    {
        {2, 3, 4, 5, 6, 1},
        {0, 7, 0, 7, 0, 7},
        {6, 1, 2, 3, 4, 5},
        {3, 4, 5, 6, 1, 2},
        {7, 0, 7, 0, 7, 0},
        {5, 6, 1, 2, 3, 4}
    };

    int Row;
    int Col;
    int Vector;

    Row = (1 - HTe) + ((1 - HF) / 2) * 3;
    Col = Sector - 1;

    if (Row < 0)
        Row = 0;
    if (Row > 5)
        Row = 5;

    if (Col < 0)
        Col = 0;
    if (Col > 5)
        Col = 5;

    Vector = DTCSVTable[Row][Col];

    /******************** 11. 电压矢量转换为开关状态 ********************/
    VectorToSwitch(Vector, &Sa, &Sb, &Sc);

    /******************** 12. 输出到 PSIM ********************/
    out[0] = (double)Sa;
    out[1] = (double)Sb;
    out[2] = (double)Sc;

    /*
     * 如果使用普通三相桥，则下桥臂一般取反。
     * 如果你的 PSIM 逆变器模块已经自动互补，则不要使用 out[3]~out[5]。
     */
    out[3] = 1.0 - (double)Fsalfa;
    out[4] = 1.0 - (double)Fsbeta;
    out[5] = 1.0 - (double)Sc;

    /*
     * 下面这些是观测量，可以接示波器观察。
     */
    out[6]  = Flux;
    out[7]  = Te;
    out[8]  = (double)Sector;
    out[9]  = (double)Vector;
    out[10] = Te_ref;
}


/************** PSIM 初始化函数 **************/
void SimulationBegin(
        const char *szId, int nInputCount, int nOutputCount,
        int nParameterCount, const char ** pszParameters,
        int *pnError, char * szErrorMsg,
        void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    Fsalfa = 0.0;
    Fsbeta = 0.0;

    speed_int = 0.0;

    HF_last = 1;

    Sa = 0;
    Sb = 0;
    Sc = 0;
}


/************** PSIM 结束函数 **************/
void SimulationEnd(
        const char *szId, void ** reserved_UserData,
        int reserved_ThreadIndex, void * reserved_AppPtr)
{
}
