#include <Stdlib.h>
#include <String.h>
#include <math.h>
#include <Psim.h>

// ===== 仿真步长 =====
double mydelt = 5e-7; // 采样时间步长 (s)，必须与 PSIM 仿真步长一致

// ===== 控制目标与限制 (非常重要) =====
double setn = 10.0;   // 设定目标转速 (rpm)
double I_max = 20.0;  // 电机最大允许电流 (A)，外环输出的电流指令不能超过这个值，起到限幅保护作用

// ===== 外环：转速环 PI 参数 =====
// 注意：双闭环中，转速环的输出是“电流(A)”。
// 如果你以前的 knp 是 40，这里绝对不能用 40，否则 1rpm 的误差会要 40A 的电流！
double knp = 20.0;      // 速度环比例系数 (A/rpm) - 需要根据你的电机重新调参
double kni = 0.05;     // 速度环积分系数 (A/(rpm*s)) - 需要根据你的电机重新调参

// ===== 内环：电流环 PI 参数 =====
// 注意：电流环的输出是“电压(V)”。
double kip = 20.0;     // 电流环比例系数 (V/A) - 需要根据你的电机重新调参
double kii = 5.0;      // 电流环积分系数 (V/(A*s)) - 需要根据你的电机重新调参

// ===== 状态变量 =====
double Ia = 0;         // 实际电流 (A)
double n = 0;          // 实际转速 (rpm)
double Uin = 0;        // 母线电压 (V)

double deltn = 0;      // 转速误差 (rpm)
double sumdeltn = 0;   // 转速积分累加值 (A)
double I_ref = 0;      // 速度环输出的期望电流 (A)

double delti = 0;      // 电流误差 (A)
double sumdelti = 0;   // 电流积分累加值 (V)
double Ua = 0;         // 电流环输出的期望电压 (V)

double duty = 0.5;     // PWM 占空比 (0~1)

void SimulationStep(
		double t, double delt, double *in, double *out,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    // 1. 读取输入 (确保你的 C-Block 接线顺序正确)
    Ia  = in[0];    // 第1个引脚：实际电枢电流
    n   = in[1];    // 第2个引脚：实际转速
    // Iin = in[2]; // 第3个引脚：输入端电流 (控制中暂不使用)
    Uin = in[3];    // 第4个引脚：母线直流电压

    // ==========================================
    // 2. 外环：转速环 控制 (Speed -> Current)
    // ==========================================
    deltn = setn - n;
    sumdeltn = sumdeltn + kni * deltn * mydelt;
    
    // 转速环积分抗饱和限制 (限幅到最大电流)
    if(sumdeltn > I_max)  sumdeltn = I_max;
    if(sumdeltn < -I_max) sumdeltn = -I_max;
    
    I_ref = knp * deltn + sumdeltn;
    
    // 转速环输出限幅 (限制目标电流，保护电机)
    if(I_ref > I_max)  I_ref = I_max;
    if(I_ref < -I_max) I_ref = -I_max;


    // ==========================================
    // 3. 内环：电流环 控制 (Current -> Voltage)
    // ==========================================
    delti = I_ref - Ia;
    sumdelti = sumdelti + kii * delti * mydelt;
    
    // 电流环积分抗饱和限制 (限幅到母线电压)
    if(sumdelti > Uin)  sumdelti = Uin;
    if(sumdelti < -Uin) sumdelti = -Uin;
    
    Ua = kip * delti + sumdelti;
    
    // 电流环输出限幅 (限制目标电压不超过母线电压)
    if(Ua > Uin)  Ua = Uin;
    if(Ua < -Uin) Ua = -Uin;


    // ==========================================
    // 4. 计算 PWM 占空比 (Voltage -> Duty)
    // ==========================================
    if(Uin > 0.01)  
    {
        // 适用于 H桥双极性调制的公式
        duty = 0.5 * Ua / Uin + 0.5;
    }
    else
    {
        duty = 0.5;  // 母线无电压时，占空比为0.5（电机两端电压为0，停转）
    }
    
    // 占空比硬限幅，防止计算溢出
    if(duty > 1.0) duty = 1.0;
    if(duty < 0.0) duty = 0.0;
    
    // 5. 输出占空比到引脚
    out[0] = duty;
}

void SimulationBegin(
		const char *szId, int nInputCount, int nOutputCount,
		 int nParameterCount, const char ** pszParameters,
		 int *pnError, char * szErrorMsg,
		 void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
    // 初始化时可将积分清零
    sumdeltn = 0;
    sumdelti = 0;
}

void SimulationEnd(const char *szId, void ** reserved_UserData, int reserved_ThreadIndex, void * reserved_AppPtr)
{
}
