#include <Stdlib.h>
#include <String.h>
#include <math.h>
#include <Psim.h>

// ===== 仿真步长 =====
double mydelt = 0.00001; 
#define PI 3.141592653589793f

// ================= 1. 外部输入与输出变量 =================
float n0 = 600.0;     // 目标转速 (RPM)
float n;               // 实际反馈转速 (RPM，来自传感器)
float Frm0 = 0.8;      // 目标转子磁链幅值 (Wb)
float ia, ib, ic;      // 三相反馈电流 (A，来自传感器)
float Udc = 600.0;     // 直流母线电压 (V)

// ================= 2. 电机铭牌物理参数  =================
double Np = 3.0;          // 极对数 = 极数6 / 2
double Rs = 0.294;        // 定子电阻 (欧姆)
double Rr = 0.156;        // 转子电阻 (欧姆)
double Lm = 0.041;        // 互感 (H)
double Ls = 0.04239;      // 定子总电感 = 漏感1.39m + 互感41m (H)
double Lr = 0.04174;      // 转子总电感 = 漏感0.74m + 互感41m (H)

double Tr = 0.26756;      // 转子时间常数 Tr = Lr / Rr
double sigma = 0.04994;   // 漏磁系数 sigma = 1 - (Lm*Lm)/(Ls*Lr)



// ================= 3. 控制器与系统参数 =================

float Immax = 25.0;    // 最大允许电流限制 (A)


float kpw = 20.0, kiw = 0.002;   // 转速环 PI 参数
float kpf = 10000.0, kif = 15.0;  // 磁链环 PI 参数
float kpi = 5000.0, kii = 10.1; // 电流环 PI 参数 

// ================= 4. 内部运行状态变量  =================
float wr0, wr;         // 目标角速度，实际角速度 (rad/s)
float isalfa, isbeta;  // 静态坐标系定子电流
float ism, ist;        // 旋转坐标系反馈电流 (d/q轴电流)
float ism0, ist0;      // 旋转坐标系目标电流 (PI输出)
float usm, ust;        // 旋转坐标系目标电压
float sm, st;          // 归一化后的电压
float salfa, sbeta;    // 静态坐标系目标电压
float Sa, Sb, Sc;      // SVPWM占空比

float delt_w, sumdelt_w = 0.0;     // 转速环误差及积分
float delt_Frm, sumdelt_Frm = 0.0; // 磁链环误差及积分
float delt_ism, sumdelt_ism = 0.0; // 励磁电流环误差及积分
float delt_ist, sumdelt_ist = 0.0; // 转矩电流环误差及积分

float ws_r, ws;        // 滑差角频率，同步角频率
float ftmp;            // 临时计算变量

// 磁链观测器状态量
float Fralfa = 0.001;  // 转子磁链 alpha 分量 
float Frbeta = 0.0;    // 转子磁链 beta 分量
float Frm = 0.001;     // 磁链幅值
float cos1 = 1.0, sin1 = 0.0; // 旋转变换角度的三角函数值

int SV_isign(float x)
{
	if(x>=0.0)return 1;
	else return 0;
}

void SVPWM_Generator(float salfa,float sbeta,float * sa,float * sb,float * sc)
{
//------------SVPWM   variable definition --------------------
	static int SV_section[6]={2,6,1,4,3,5};
	static int SV_N,SV_iSection;
	static float SV_A1,SV_B1,SV_C1;
	static float SV_T0,SV_T1,SV_T2;
	static float SV_uu,SV_uv,SV_uw;

	//---------------------------------------------------------------------------------
	SV_A1=sbeta;
	SV_B1=1.7320508*salfa-sbeta;
	SV_C1=-1.7320508*salfa-sbeta;
	SV_N=SV_isign(SV_A1)+2*SV_isign(SV_B1)+4*SV_isign(SV_C1);
	if(SV_N==0||SV_N==7)
	{
		*sa=0.0;
		*sb=0.0;
		*sc=0.0;
		return;
	}
	SV_iSection=SV_section[SV_N-1];
	//-----                --------------------
	salfa=1.5*salfa;
	sbeta=1.5*sbeta;
	switch(SV_iSection)
	{
	case 1:
            SV_T1=salfa-sbeta/1.7320508;
            SV_T2= 1.154700538379252*sbeta;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1.0)
            {
		SV_T1=SV_T1/(SV_T1+SV_T2);
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0+SV_T1+SV_T2;
            SV_uv=SV_T0/2.0+SV_T2;
            SV_uw=SV_T0/2.0;
		break;
	case 2:
            SV_T1=salfa+sbeta/1.7320508;
            SV_T2=-salfa+sbeta/1.7320508;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1)
	   {
                SV_T1=SV_T1/(SV_T1+SV_T2);
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0+SV_T1;
            SV_uv=SV_T0/2.0+SV_T1+SV_T2;
            SV_uw=SV_T0/2.0;
		break;
	case 3:
            SV_T1= 1.154700538379252*sbeta;
            SV_T2=-salfa-sbeta/1.7320508;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1.0)
	   {
                SV_T1=SV_T1/(SV_T1+SV_T2);
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0;
            SV_uv=SV_T0/2.0+SV_T1+SV_T2;
            SV_uw=SV_T0/2.0+SV_T2;
		break;
	case 4:
            SV_T1=-salfa+sbeta/1.7320508;
            SV_T2=- 1.154700538379252*sbeta;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1.0)
	   {
                SV_T1=SV_T1/(SV_T1+SV_T2);
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0;
            SV_uv=SV_T0/2.0+SV_T1;
            SV_uw=SV_T0/2.0+SV_T1+SV_T2;
		break;
	case 5:
            SV_T1=-salfa-sbeta/1.7320508;
            SV_T2=salfa-sbeta/1.7320508;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1.0)
	    {
                SV_T1=SV_T1/(SV_T1+SV_T2);
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0+SV_T2;
            SV_uv=SV_T0/2.0;
            SV_uw=SV_T0/2.0+SV_T1+SV_T2;
		break;
	case 6:
            SV_T1=- 1.154700538379252*sbeta;
            SV_T2=salfa+sbeta/1.7320508;
            SV_T0=1.0-SV_T1-SV_T2;
            if((SV_T1+SV_T2)>1.0)
	    {
                SV_T1=SV_T1/(SV_T1+SV_T2);  
                SV_T2=SV_T2/(SV_T1+SV_T2);
                SV_T0=0.0;
            }
            SV_uu=SV_T0/2.0+SV_T1+SV_T2;
            SV_uv=SV_T0/2.0;
            SV_uw=SV_T0/2.0+SV_T1;
		break;
	default:
		break;
	}
	*sa=SV_uu;
	*sb=SV_uv;
	*sc=SV_uw;
}
void SimulationStep(
		double t, double delt, double *in, double *out,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    
    n  = in[0];  // 第1个输入：实际转速
    ia = in[1];  // 第2个输入：A相电流
    ib = in[2];  // 第3个输入：B相电流
    ic = in[3];  // 第4个输入：C相电流

    
    wr0 = Np * 2.0 * PI * n0 / 60.0; 
    wr  = Np * 2.0 * PI * n / 60.0; 
    isalfa = 2.0 * (ia - 0.5 * ib - 0.5 * ic) / 3.0; 
    isbeta = 2.0 * (sqrt(3.0) / 2.0 * (ib - ic)) / 3.0; 
    ism = cos1 * isalfa + sin1 * isbeta; 
    ist = -sin1 * isalfa + cos1 * isbeta; 
    
    delt_w = wr0 - wr; 
    delt_Frm = Frm0 - Frm; 
    
    sumdelt_w = sumdelt_w + kiw * delt_w * mydelt; 
    sumdelt_Frm = sumdelt_Frm + kif * delt_Frm * mydelt; 
    
    ftmp = 30.0; 
    if(sumdelt_w > ftmp) sumdelt_w = ftmp; 
    if(sumdelt_w < -ftmp) sumdelt_w = -ftmp; 
    
    ftmp = 40.0; 
    if(sumdelt_Frm > ftmp) sumdelt_Frm = ftmp; 
    if(sumdelt_Frm < -ftmp) sumdelt_Frm = -ftmp; 
    
    ist0 = kpw * delt_w + sumdelt_w; 
    ism0 = kpf * delt_Frm + sumdelt_Frm; 
    
    ftmp = sqrt(ist0 * ist0 + ism0 * ism0); 
    if(ftmp > Immax)  
    { 
        ist0 = Immax * ist0 / ftmp; 
        //ism0 = Immax * ism0 / ftmp; 
    } 
    
    ws_r = Lm * ist / (Tr * Frm); 
    ws = wr + ws_r; 
    
    delt_ism = ism0 - ism; 
    delt_ist = ist0 - ist; 
    
    sumdelt_ism = sumdelt_ism + kii * delt_ism * mydelt; 
    sumdelt_ist = sumdelt_ist + kii * delt_ist * mydelt; 
    
    ftmp = 2000.0; 
    if(sumdelt_ism > ftmp) sumdelt_ism = ftmp; 
    if(sumdelt_ism < -ftmp) sumdelt_ism = -ftmp; 
    if(sumdelt_ist > ftmp) sumdelt_ist = ftmp; 
    if(sumdelt_ist < -ftmp) sumdelt_ist = -ftmp; 
    
    //-------------电流调节器-------------------------------------- 
    usm = (kpi * delt_ism + sumdelt_ism - ws * ist - Frm * Lm / (sigma * Ls * Lr) + (Rs * Lr * Lr + Rr * Lm * Lm) / (sigma * Ls * Lr * Lr) * ism) * sigma * Ls; 
    ust = (kpi * delt_ist + sumdelt_ist + ws * ism + Frm * Lm * ws / (sigma * Ls * Lr) + (Rs * Lr * Lr + Rr * Lm * Lm) / (sigma * Ls * Lr * Lr) * ist) * sigma * Ls; 
    
    sm = usm / Udc; 
    st = ust / Udc; 
    
    ftmp = sqrt(sm * sm + st * st) + 0.0001; 
    if(ftmp > (sqrt(3.0) / 3.0)) 
    { 
        sm = sqrt(3.0) / 3.0 * sm / ftmp; 
        st = sqrt(3.0) / 3.0 * st / ftmp; 
    } 
    
    salfa = cos1 * sm - sin1 * st; 
    sbeta = sin1 * sm + cos1 * st; 
    
    SVPWM_Generator(salfa, sbeta, &Sa, &Sb, &Sc); 
    
    Fralfa = Fralfa + (-wr * Frbeta - (Fralfa - Lm * isalfa) / Tr) * mydelt; 
    Frbeta = Frbeta + (wr * Fralfa - (Frbeta - Lm * isbeta) / Tr) * mydelt; 
    Frm = sqrt(Fralfa * Fralfa + Frbeta * Frbeta + 0.001); 
    cos1 = Fralfa / Frm; 
    sin1 = Frbeta / Frm; 
	float Tem_cal = 1.5 * Np * (Lm / Lr) * Frm * ist;

    
    out[0] = Sa; 
    out[1] = Sb; 
    out[2] = Sc;
	out[12]=isalfa;
	out[13]=isbeta;
	out[14]=Fralfa;
	out[15]=Frbeta;
	out[16]=Tem_cal;
}

void SimulationBegin(
		const char *szId, int nInputCount, int nOutputCount,
		 int nParameterCount, const char ** pszParameters,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    
    sumdelt_w = 0.0;
    sumdelt_Frm = 0.0;
    sumdelt_ism = 0.0;
    sumdelt_ist = 0.0;
    
    Fralfa = 0.001;
    Frbeta = 0.0;
    Frm = 0.001;
    cos1 = 1.0;
    sin1 = 0.0;
}

void SimulationEnd(const char *szId, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
}
