//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_capi.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2406
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sun Aug  6 01:50:56 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#include "rtw_capi.h"
#ifdef HOST_CAPI_BUILD
#include "SupervisoryController_capi_host.h"
#define sizeof(s)                      ((size_t)(0xFFFF))
#undef rt_offsetof
#define rt_offsetof(s,el)              ((uint16_T)(0xFFFF))
#define TARGET_CONST
#define TARGET_STRING(s)               (s)
#ifndef SS_UINT64
#define SS_UINT64                      24
#endif

#ifndef SS_INT64
#define SS_INT64                       25
#endif

#else                                  // HOST_CAPI_BUILD
#include "builtin_typeid_types.h"
#include "SupervisoryController.h"
#include "SupervisoryController_capi.h"
#ifdef LIGHT_WEIGHT_CAPI
#define TARGET_CONST
#define TARGET_STRING(s)               ((nullptr))
#else
#define TARGET_CONST                   const
#define TARGET_STRING(s)               (s)
#endif
#endif                                 // HOST_CAPI_BUILD

static rtwCAPI_BlockParameters rtBlockParameters[]{
  // addrMapIndex, blockPath,
  //  paramName, dataTypeIndex, dimIndex, fixPtIdx

  { 0, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 1, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/ywt"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 2, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 3, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/r_"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 4, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 5, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant12"),
    TARGET_STRING("Value"), 0, 1, 0 },

  { 6, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant13"),
    TARGET_STRING("Value"), 0, 2, 0 },

  { 7, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant2"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 8, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Gain2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 9, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Gain3"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 10, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("Mean"), 0, 4, 0 },

  { 11, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("StdDev"), 0, 4, 0 },

  { 12, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("Seed"), 0, 4, 0 },

  { 13, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 14, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 15, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Discrete Filter1"),
    TARGET_STRING("InitialStates"), 0, 0, 0 },

  { 16, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 17, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 18, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/r_"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 19, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/ywt"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 20, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 21, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 22, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant12"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 23, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant13"),
    TARGET_STRING("Value"), 0, 6, 0 },

  { 24, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant2"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 25, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant3"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 26, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant4"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 27, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Discrete-Time Integrator"),
    TARGET_STRING("gainval"), 0, 0, 0 },

  { 28, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Discrete-Time Integrator"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 29, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/theta"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 30, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/P"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 31, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/err"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 32, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/theta"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 33, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/P"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 34, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/err"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 35, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 36, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/F_zero"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 37, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 38, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 39, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 40, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 41, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 42, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 43, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 44, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 45, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 46, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 8, 0 },

  { 47, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 48, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/F_zero"),
    TARGET_STRING("Value"), 0, 9, 0 },

  { 49, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 50, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 51, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 52, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 53, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 54, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 55, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 56, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 57, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/u.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 58, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 59, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/y.wt_zero"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 60, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ymax_zero"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 61, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ymin_zero"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 62, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 63, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/Param Estimator (RLS)/Unit Delay3"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 64, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/Param Estimator (RLS)/Unit Delay3"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 65, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 66, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 67, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/u_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 68, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 69, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/uref_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 70, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 7, 0 },

  { 71, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 72, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/LastPcov"),
    TARGET_STRING("InitialCondition"), 0, 11, 0 },

  { 73, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 12, 0 },

  { 74, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 75, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 11, 0 },

  { 76, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 13, 0 },

  { 77, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 11, 0 },

  { 78, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 14, 0 },

  { 79, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 80, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 81, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ext.mv_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 82, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/umin_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 83, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 84, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 9, 0 },

  { 85, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 86, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 15, 0 },

  { 87, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 88, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 89, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 17, 0 },

  { 90, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 91, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 92, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer"),
    TARGET_STRING("Ndis"), 2, 0, 0 },

  { 93, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 94, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 95, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 96, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  {
    0, (nullptr), (nullptr), 0, 0, 0
  }
};

// Tunable variable parameters
static rtwCAPI_ModelParameters rtModelParameters[]{
  // addrMapIndex, varName, dataTypeIndex, dimIndex, fixPtIndex
  { 97, TARGET_STRING("nullEv"), 3, 0, 0 },

  { 98, TARGET_STRING("Aod"), 0, 20, 0 },

  { 99, TARGET_STRING("Aod1"), 0, 6, 0 },

  { 100, TARGET_STRING("Bod"), 0, 21, 0 },

  { 101, TARGET_STRING("Bod1"), 0, 6, 0 },

  { 102, TARGET_STRING("Cod"), 0, 22, 0 },

  { 103, TARGET_STRING("Cod1"), 0, 6, 0 },

  { 104, TARGET_STRING("Dmn"), 0, 1, 0 },

  { 105, TARGET_STRING("Dmn1"), 0, 6, 0 },

  { 106, TARGET_STRING("Dod"), 0, 1, 0 },

  { 107, TARGET_STRING("Dod1"), 0, 6, 0 },

  { 108, TARGET_STRING("beta"), 0, 0, 0 },

  { 109, TARGET_STRING("dt"), 0, 0, 0 },

  { 110, TARGET_STRING("lpfDen"), 0, 0, 0 },

  { 111, TARGET_STRING("lpfNum"), 0, 23, 0 },

  { 112, TARGET_STRING("mdlNum"), 0, 0, 0 },

  { 113, TARGET_STRING("uwt0"), 0, 4, 0 },

  { 114, TARGET_STRING("ywt0"), 0, 7, 0 },

  { 0, (nullptr), 0, 0, 0 }
};

#ifndef HOST_CAPI_BUILD

// Initialize Data Address
static void InitializeDataAddr(void* dataAddr[], SupervisoryController::P *rtP)
{
  dataAddr[0] = (void*) (&rtP->u_Y0);
  dataAddr[1] = (void*) (&rtP->ywt_Y0);
  dataAddr[2] = (void*) (&rtP->yhat_Y0);
  dataAddr[3] = (void*) (&rtP->r_Y0);
  dataAddr[4] = (void*) (&rtP->Constant1_Value_c);
  dataAddr[5] = (void*) (&rtP->Constant12_Value[0]);
  dataAddr[6] = (void*) (&rtP->Constant13_Value[0]);
  dataAddr[7] = (void*) (&rtP->Constant2_Value[0]);
  dataAddr[8] = (void*) (&rtP->Gain2_Gain);
  dataAddr[9] = (void*) (&rtP->Gain3_Gain);
  dataAddr[10] = (void*) (&rtP->excitation_Mean[0]);
  dataAddr[11] = (void*) (&rtP->excitation_StdDev[0]);
  dataAddr[12] = (void*) (&rtP->excitation_Seed[0]);
  dataAddr[13] = (void*) (&rtP->Saturation_UpperSat);
  dataAddr[14] = (void*) (&rtP->Saturation_LowerSat);
  dataAddr[15] = (void*) (&rtP->DiscreteFilter1_InitialStates);
  dataAddr[16] = (void*) (&rtP->u_Y0_b);
  dataAddr[17] = (void*) (&rtP->yhat_Y0_d);
  dataAddr[18] = (void*) (&rtP->r_Y0_m);
  dataAddr[19] = (void*) (&rtP->ywt_Y0_m);
  dataAddr[20] = (void*) (&rtP->Constant_Value);
  dataAddr[21] = (void*) (&rtP->Constant1_Value_e);
  dataAddr[22] = (void*) (&rtP->Constant12_Value_e[0]);
  dataAddr[23] = (void*) (&rtP->Constant13_Value_c[0]);
  dataAddr[24] = (void*) (&rtP->Constant2_Value_a);
  dataAddr[25] = (void*) (&rtP->Constant3_Value);
  dataAddr[26] = (void*) (&rtP->Constant4_Value[0]);
  dataAddr[27] = (void*) (&rtP->DiscreteTimeIntegrator_gainval);
  dataAddr[28] = (void*) (&rtP->DiscreteTimeIntegrator_IC);
  dataAddr[29] = (void*) (&rtP->paramEst1_o.theta_Y0);
  dataAddr[30] = (void*) (&rtP->paramEst1_o.P_Y0);
  dataAddr[31] = (void*) (&rtP->paramEst1_o.err_Y0);
  dataAddr[32] = (void*) (&rtP->paramEst2.theta_Y0);
  dataAddr[33] = (void*) (&rtP->paramEst2.P_Y0);
  dataAddr[34] = (void*) (&rtP->paramEst2.err_Y0);
  dataAddr[35] = (void*) (&rtP->E_zero_Value[0]);
  dataAddr[36] = (void*) (&rtP->F_zero_Value[0]);
  dataAddr[37] = (void*) (&rtP->G_zero_Value);
  dataAddr[38] = (void*) (&rtP->S_zero_Value);
  dataAddr[39] = (void*) (&rtP->duwt_zero_Value[0]);
  dataAddr[40] = (void*) (&rtP->ecrwt_zero_Value);
  dataAddr[41] = (void*) (&rtP->extmv_zero_Value[0]);
  dataAddr[42] = (void*) (&rtP->md_zero_Value);
  dataAddr[43] = (void*) (&rtP->mvtarget_zero_Value[0]);
  dataAddr[44] = (void*) (&rtP->switch_zero_Value);
  dataAddr[45] = (void*) (&rtP->umin_zero_Value[0]);
  dataAddr[46] = (void*) (&rtP->Constant1_Value[0]);
  dataAddr[47] = (void*) (&rtP->E_zero_Value_a[0]);
  dataAddr[48] = (void*) (&rtP->F_zero_Value_g[0]);
  dataAddr[49] = (void*) (&rtP->G_zero_Value_m);
  dataAddr[50] = (void*) (&rtP->S_zero_Value_g);
  dataAddr[51] = (void*) (&rtP->duwt_zero_Value_p[0]);
  dataAddr[52] = (void*) (&rtP->ecrwt_zero_Value_o);
  dataAddr[53] = (void*) (&rtP->extmv_zero_Value_e[0]);
  dataAddr[54] = (void*) (&rtP->md_zero_Value_p);
  dataAddr[55] = (void*) (&rtP->mvtarget_zero_Value_k[0]);
  dataAddr[56] = (void*) (&rtP->switch_zero_Value_e);
  dataAddr[57] = (void*) (&rtP->uwt_zero_Value[0]);
  dataAddr[58] = (void*) (&rtP->umin_zero_Value_d[0]);
  dataAddr[59] = (void*) (&rtP->ywt_zero_Value[0]);
  dataAddr[60] = (void*) (&rtP->ymax_zero_Value[0]);
  dataAddr[61] = (void*) (&rtP->ymin_zero_Value[0]);
  dataAddr[62] = (void*) (&rtP->Constant1_Value_j[0]);
  dataAddr[63] = (void*) (&rtP->paramEst1_o.UnitDelay3_InitialCondition);
  dataAddr[64] = (void*) (&rtP->paramEst2.UnitDelay3_InitialCondition);
  dataAddr[65] = (void*) (&rtP->ym_zero_Value[0]);
  dataAddr[66] = (void*) (&rtP->extmv_scale_Gain[0]);
  dataAddr[67] = (void*) (&rtP->u_scale_Gain[0]);
  dataAddr[68] = (void*) (&rtP->umin_scale4_Gain[0]);
  dataAddr[69] = (void*) (&rtP->uref_scale_Gain[0]);
  dataAddr[70] = (void*) (&rtP->ymin_scale1_Gain[0]);
  dataAddr[71] = (void*) (&rtP->ymin_scale2_Gain);
  dataAddr[72] = (void*) (&rtP->LastPcov_InitialCondition[0]);
  dataAddr[73] = (void*) (&rtP->Memory_InitialCondition[0]);
  dataAddr[74] = (void*) (&rtP->last_mv_InitialCondition[0]);
  dataAddr[75] = (void*) (&rtP->G_Value[0]);
  dataAddr[76] = (void*) (&rtP->H_Value[0]);
  dataAddr[77] = (void*) (&rtP->P0_Value[0]);
  dataAddr[78] = (void*) (&rtP->X0_Value[0]);
  dataAddr[79] = (void*) (&rtP->ym_zero_Value_c[0]);
  dataAddr[80] = (void*) (&rtP->extmv_scale_Gain_e[0]);
  dataAddr[81] = (void*) (&rtP->extmv_scale1_Gain[0]);
  dataAddr[82] = (void*) (&rtP->umin_scale1_Gain[0]);
  dataAddr[83] = (void*) (&rtP->umin_scale4_Gain_p[0]);
  dataAddr[84] = (void*) (&rtP->ymin_scale1_Gain_j[0]);
  dataAddr[85] = (void*) (&rtP->ymin_scale2_Gain_f);
  dataAddr[86] = (void*) (&rtP->Memory_InitialCondition_f[0]);
  dataAddr[87] = (void*) (&rtP->last_mv_InitialCondition_f[0]);
  dataAddr[88] = (void*) (&rtP->G_Value_a[0]);
  dataAddr[89] = (void*) (&rtP->H_Value_o[0]);
  dataAddr[90] = (void*) (&rtP->P0_Value_a[0]);
  dataAddr[91] = (void*) (&rtP->X0_Value_f[0]);
  dataAddr[92] = (void*) (&rtP->FixedHorizonOptimizer_Ndis);
  dataAddr[93] = (void*) (&rtP->isSqrtUsed_Value);
  dataAddr[94] = (void*) (&rtP->Lykyhatkk1_Y0);
  dataAddr[95] = (void*) (&rtP->isSqrtUsed_Value_d);
  dataAddr[96] = (void*) (&rtP->Lykyhatkk1_Y0_c);
  dataAddr[97] = (void*) (&rtP->nullEv);
  dataAddr[98] = (void*) (&rtP->Aod[0]);
  dataAddr[99] = (void*) (&rtP->Aod1[0]);
  dataAddr[100] = (void*) (&rtP->Bod[0]);
  dataAddr[101] = (void*) (&rtP->Bod1[0]);
  dataAddr[102] = (void*) (&rtP->Cod[0]);
  dataAddr[103] = (void*) (&rtP->Cod1[0]);
  dataAddr[104] = (void*) (&rtP->Dmn[0]);
  dataAddr[105] = (void*) (&rtP->Dmn1[0]);
  dataAddr[106] = (void*) (&rtP->Dod[0]);
  dataAddr[107] = (void*) (&rtP->Dod1[0]);
  dataAddr[108] = (void*) (&rtP->beta);
  dataAddr[109] = (void*) (&rtP->dt);
  dataAddr[110] = (void*) (&rtP->lpfDen);
  dataAddr[111] = (void*) (&rtP->lpfNum[0]);
  dataAddr[112] = (void*) (&rtP->mdlNum);
  dataAddr[113] = (void*) (&rtP->uwt0[0]);
  dataAddr[114] = (void*) (&rtP->ywt0[0]);
}

#endif

// Initialize Data Run-Time Dimension Buffer Address
#ifndef HOST_CAPI_BUILD

static void InitializeVarDimsAddr(int32_T* vardimsAddr[])
{
  vardimsAddr[0] = (nullptr);
}

#endif

#ifndef HOST_CAPI_BUILD

// Initialize logging function pointers
static void InitializeLoggingFunctions(RTWLoggingFcnPtr loggingPtrs[])
{
  loggingPtrs[0] = (nullptr);
  loggingPtrs[1] = (nullptr);
  loggingPtrs[2] = (nullptr);
  loggingPtrs[3] = (nullptr);
  loggingPtrs[4] = (nullptr);
  loggingPtrs[5] = (nullptr);
  loggingPtrs[6] = (nullptr);
  loggingPtrs[7] = (nullptr);
  loggingPtrs[8] = (nullptr);
  loggingPtrs[9] = (nullptr);
  loggingPtrs[10] = (nullptr);
  loggingPtrs[11] = (nullptr);
  loggingPtrs[12] = (nullptr);
  loggingPtrs[13] = (nullptr);
  loggingPtrs[14] = (nullptr);
  loggingPtrs[15] = (nullptr);
  loggingPtrs[16] = (nullptr);
  loggingPtrs[17] = (nullptr);
  loggingPtrs[18] = (nullptr);
  loggingPtrs[19] = (nullptr);
  loggingPtrs[20] = (nullptr);
  loggingPtrs[21] = (nullptr);
  loggingPtrs[22] = (nullptr);
  loggingPtrs[23] = (nullptr);
  loggingPtrs[24] = (nullptr);
  loggingPtrs[25] = (nullptr);
  loggingPtrs[26] = (nullptr);
  loggingPtrs[27] = (nullptr);
  loggingPtrs[28] = (nullptr);
  loggingPtrs[29] = (nullptr);
  loggingPtrs[30] = (nullptr);
  loggingPtrs[31] = (nullptr);
  loggingPtrs[32] = (nullptr);
  loggingPtrs[33] = (nullptr);
  loggingPtrs[34] = (nullptr);
  loggingPtrs[35] = (nullptr);
  loggingPtrs[36] = (nullptr);
  loggingPtrs[37] = (nullptr);
  loggingPtrs[38] = (nullptr);
  loggingPtrs[39] = (nullptr);
  loggingPtrs[40] = (nullptr);
  loggingPtrs[41] = (nullptr);
  loggingPtrs[42] = (nullptr);
  loggingPtrs[43] = (nullptr);
  loggingPtrs[44] = (nullptr);
  loggingPtrs[45] = (nullptr);
  loggingPtrs[46] = (nullptr);
  loggingPtrs[47] = (nullptr);
  loggingPtrs[48] = (nullptr);
  loggingPtrs[49] = (nullptr);
  loggingPtrs[50] = (nullptr);
  loggingPtrs[51] = (nullptr);
  loggingPtrs[52] = (nullptr);
  loggingPtrs[53] = (nullptr);
  loggingPtrs[54] = (nullptr);
  loggingPtrs[55] = (nullptr);
  loggingPtrs[56] = (nullptr);
  loggingPtrs[57] = (nullptr);
  loggingPtrs[58] = (nullptr);
  loggingPtrs[59] = (nullptr);
  loggingPtrs[60] = (nullptr);
  loggingPtrs[61] = (nullptr);
  loggingPtrs[62] = (nullptr);
  loggingPtrs[63] = (nullptr);
  loggingPtrs[64] = (nullptr);
  loggingPtrs[65] = (nullptr);
  loggingPtrs[66] = (nullptr);
  loggingPtrs[67] = (nullptr);
  loggingPtrs[68] = (nullptr);
  loggingPtrs[69] = (nullptr);
  loggingPtrs[70] = (nullptr);
  loggingPtrs[71] = (nullptr);
  loggingPtrs[72] = (nullptr);
  loggingPtrs[73] = (nullptr);
  loggingPtrs[74] = (nullptr);
  loggingPtrs[75] = (nullptr);
  loggingPtrs[76] = (nullptr);
  loggingPtrs[77] = (nullptr);
  loggingPtrs[78] = (nullptr);
  loggingPtrs[79] = (nullptr);
  loggingPtrs[80] = (nullptr);
  loggingPtrs[81] = (nullptr);
  loggingPtrs[82] = (nullptr);
  loggingPtrs[83] = (nullptr);
  loggingPtrs[84] = (nullptr);
  loggingPtrs[85] = (nullptr);
  loggingPtrs[86] = (nullptr);
  loggingPtrs[87] = (nullptr);
  loggingPtrs[88] = (nullptr);
  loggingPtrs[89] = (nullptr);
  loggingPtrs[90] = (nullptr);
  loggingPtrs[91] = (nullptr);
  loggingPtrs[92] = (nullptr);
  loggingPtrs[93] = (nullptr);
  loggingPtrs[94] = (nullptr);
  loggingPtrs[95] = (nullptr);
  loggingPtrs[96] = (nullptr);
  loggingPtrs[97] = (nullptr);
  loggingPtrs[98] = (nullptr);
  loggingPtrs[99] = (nullptr);
  loggingPtrs[100] = (nullptr);
  loggingPtrs[101] = (nullptr);
  loggingPtrs[102] = (nullptr);
  loggingPtrs[103] = (nullptr);
  loggingPtrs[104] = (nullptr);
  loggingPtrs[105] = (nullptr);
  loggingPtrs[106] = (nullptr);
  loggingPtrs[107] = (nullptr);
  loggingPtrs[108] = (nullptr);
  loggingPtrs[109] = (nullptr);
  loggingPtrs[110] = (nullptr);
  loggingPtrs[111] = (nullptr);
  loggingPtrs[112] = (nullptr);
  loggingPtrs[113] = (nullptr);
  loggingPtrs[114] = (nullptr);
}

#endif

// Data Type Map - use dataTypeMapIndex to access this structure
static TARGET_CONST rtwCAPI_DataTypeMap rtDataTypeMap[]{
  // cName, mwName, numElements, elemMapIndex, dataSize, slDataId, *
  //  isComplex, isPointer, enumStorageType
  { "double", "real_T", 0, 0, sizeof(real_T), (uint8_T)SS_DOUBLE, 0, 0, 0 },

  { "unsigned char", "boolean_T", 0, 0, sizeof(boolean_T), (uint8_T)SS_BOOLEAN,
    0, 0, 0 },

  { "int", "int32_T", 0, 0, sizeof(int32_T), (uint8_T)SS_INT32, 0, 0, 0 },

  { "struct", "event_bus", 4, 1, sizeof(event_bus), (uint8_T)SS_STRUCT, 0, 0, 0
  }
};

#ifdef HOST_CAPI_BUILD
#undef sizeof
#endif

// Structure Element Map - use elemMapIndex to access this structure
static TARGET_CONST rtwCAPI_ElementMap rtElementMap[]{
  // elementName, elementOffset, dataTypeIndex, dimIndex, fxpIndex
  { (nullptr), 0, 0, 0, 0 },

  { "r", rt_offsetof(event_bus, r), 0, 18, 0 },

  { "preT", rt_offsetof(event_bus, preT), 0, 19, 0 },

  { "moveT", rt_offsetof(event_bus, moveT), 0, 19, 0 },

  { "postT", rt_offsetof(event_bus, postT), 0, 19, 0 }
};

// Dimension Map - use dimensionMapIndex to access elements of ths structure
static rtwCAPI_DimensionMap rtDimensionMap[]{
  // dataOrientation, dimArrayIndex, numDims, vardimsIndex
  { rtwCAPI_SCALAR, 0, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 2, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 4, 2, 0 },

  { rtwCAPI_VECTOR, 6, 2, 0 },

  { rtwCAPI_VECTOR, 8, 2, 0 },

  { rtwCAPI_VECTOR, 10, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 12, 2, 0 },

  { rtwCAPI_VECTOR, 14, 2, 0 },

  { rtwCAPI_VECTOR, 16, 2, 0 },

  { rtwCAPI_VECTOR, 18, 2, 0 },

  { rtwCAPI_VECTOR, 20, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 22, 2, 0 },

  { rtwCAPI_VECTOR, 24, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 26, 2, 0 },

  { rtwCAPI_VECTOR, 28, 2, 0 },

  { rtwCAPI_VECTOR, 30, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 32, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 34, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 6, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 0, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 36, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 38, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 40, 2, 0 },

  { rtwCAPI_VECTOR, 42, 2, 0 }
};

// Dimension Array- use dimArrayIndex to access elements of this array
static uint_T rtDimensionArray[]{
  1,                                   // 0
  1,                                   // 1
  6,                                   // 2
  6,                                   // 3
  6,                                   // 4
  3,                                   // 5
  6,                                   // 6
  1,                                   // 7
  1,                                   // 8
  3,                                   // 9
  3,                                   // 10
  1,                                   // 11
  3,                                   // 12
  3,                                   // 13
  1,                                   // 14
  6,                                   // 15
  12,                                  // 16
  1,                                   // 17
  1,                                   // 18
  4,                                   // 19
  4,                                   // 20
  1,                                   // 21
  18,                                  // 22
  18,                                  // 23
  246,                                 // 24
  1,                                   // 25
  6,                                   // 26
  18,                                  // 27
  18,                                  // 28
  1,                                   // 29
  166,                                 // 30
  1,                                   // 31
  4,                                   // 32
  4,                                   // 33
  3,                                   // 34
  4,                                   // 35
  12,                                  // 36
  12,                                  // 37
  12,                                  // 38
  6,                                   // 39
  6,                                   // 40
  12,                                  // 41
  1,                                   // 42
  60                                   // 43
};

// Fixed Point Map
static rtwCAPI_FixPtMap rtFixPtMap[]{
  // fracSlopePtr, biasPtr, scaleType, wordLength, exponent, isSigned
  { (nullptr), (nullptr), rtwCAPI_FIX_RESERVED, 0, 0, (boolean_T)0 },
};

// Sample Time Map - use sTimeIndex to access elements of ths structure
static rtwCAPI_SampleTimeMap rtSampleTimeMap[]{
  // samplePeriodPtr, sampleOffsetPtr, tid, samplingMode
  {
    (nullptr), (nullptr), 0, 0
  }
};

static rtwCAPI_ModelMappingStaticInfo mmiStatic{
  // Signals:{signals, numSignals,
  //            rootInputs, numRootInputs,
  //            rootOutputs, numRootOutputs},
  //  Params: {blockParameters, numBlockParameters,
  //           modelParameters, numModelParameters},
  //  States: {states, numStates},
  //  Maps:   {dataTypeMap, dimensionMap, fixPtMap,
  //           elementMap, sampleTimeMap, dimensionArray},
  //  TargetType: targetType

  { (nullptr), 0,
    (nullptr), 0,
    (nullptr), 0 },

  { rtBlockParameters, 97,
    rtModelParameters, 18 },

  { (nullptr), 0 },

  { rtDataTypeMap, rtDimensionMap, rtFixPtMap,
    rtElementMap, rtSampleTimeMap, rtDimensionArray },
  "float",

  { 2318526301U,
    3658225509U,
    3904677418U,
    4250369410U },
  (nullptr), 0,
  (boolean_T)0
};

// Function to get C API Model Mapping Static Info
const rtwCAPI_ModelMappingStaticInfo*
  SupervisoryController_GetCAPIStaticMap(void)
{
  return &mmiStatic;
}

// Cache pointers into DataMapInfo substructure of RTModel
#ifndef HOST_CAPI_BUILD

void SupervisoryController_InitializeDataMapInfo(SupervisoryController::RT_MODEL
  *const rtM, SupervisoryController::P *rtP)
{
  // Set C-API version
  rtwCAPI_SetVersion(rtM->DataMapInfo.mmi, 1);

  // Cache static C-API data into the Real-time Model Data structure
  rtwCAPI_SetStaticMap(rtM->DataMapInfo.mmi, &mmiStatic);

  // Cache static C-API logging data into the Real-time Model Data structure
  rtwCAPI_SetLoggingStaticMap(rtM->DataMapInfo.mmi, (nullptr));

  // Cache C-API Data Addresses into the Real-Time Model Data structure
  InitializeDataAddr(rtM->DataMapInfo.dataAddress, rtP);
  rtwCAPI_SetDataAddressMap(rtM->DataMapInfo.mmi, rtM->DataMapInfo.dataAddress);

  // Cache C-API Data Run-Time Dimension Buffer Addresses into the Real-Time Model Data structure 
  InitializeVarDimsAddr(rtM->DataMapInfo.vardimsAddress);
  rtwCAPI_SetVarDimsAddressMap(rtM->DataMapInfo.mmi,
    rtM->DataMapInfo.vardimsAddress);

  // Set Instance specific path
  rtwCAPI_SetPath(rtM->DataMapInfo.mmi, (nullptr));
  rtwCAPI_SetFullPath(rtM->DataMapInfo.mmi, (nullptr));

  // Cache C-API logging function pointers into the Real-Time Model Data structure 
  InitializeLoggingFunctions(rtM->DataMapInfo.loggingPtrs);
  rtwCAPI_SetLoggingPtrs(rtM->DataMapInfo.mmi, rtM->DataMapInfo.loggingPtrs);

  // Cache the instance C-API logging pointer
  rtwCAPI_SetInstanceLoggingInfo(rtM->DataMapInfo.mmi, (nullptr));

  // Set reference to submodels
  rtwCAPI_SetChildMMIArray(rtM->DataMapInfo.mmi, (nullptr));
  rtwCAPI_SetChildMMIArrayLen(rtM->DataMapInfo.mmi, 0);
}

#else                                  // HOST_CAPI_BUILD
#ifdef __cplusplus

extern "C"
{

#endif

  void SupervisoryController_host_InitializeDataMapInfo
    (SupervisoryController_host_DataMapInfo_T *dataMap, const char *path)
  {
    // Set C-API version
    rtwCAPI_SetVersion(dataMap->mmi, 1);

    // Cache static C-API data into the Real-time Model Data structure
    rtwCAPI_SetStaticMap(dataMap->mmi, &mmiStatic);

    // host data address map is NULL
    rtwCAPI_SetDataAddressMap(dataMap->mmi, (nullptr));

    // host vardims address map is NULL
    rtwCAPI_SetVarDimsAddressMap(dataMap->mmi, (nullptr));

    // Set Instance specific path
    rtwCAPI_SetPath(dataMap->mmi, path);
    rtwCAPI_SetFullPath(dataMap->mmi, (nullptr));

    // Set reference to submodels
    rtwCAPI_SetChildMMIArray(dataMap->mmi, (nullptr));
    rtwCAPI_SetChildMMIArrayLen(dataMap->mmi, 0);
  }

#ifdef __cplusplus

}

#endif
#endif                                 // HOST_CAPI_BUILD

//
// File trailer for generated code.
//
// [EOF]
//
