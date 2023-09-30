//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_capi.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2507
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sat Sep 30 07:51:38 2023
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
#define SS_UINT64                      26
#endif

#ifndef SS_INT64
#define SS_INT64                       27
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

  { 1, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 2, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 3, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant12"),
    TARGET_STRING("Value"), 0, 1, 0 },

  { 4, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant13"),
    TARGET_STRING("Value"), 0, 2, 0 },

  { 5, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Constant2"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 6, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Gain2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 7, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Gain3"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 8, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("Mean"), 0, 4, 0 },

  { 9, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("StdDev"), 0, 4, 0 },

  { 10, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/excitation"),
    TARGET_STRING("Seed"), 0, 4, 0 },

  { 11, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 12, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 13, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Discrete Filter1"),
    TARGET_STRING("InitialStates"), 0, 0, 0 },

  { 14, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 15, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 16, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 17, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 18, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant12"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 19, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant13"),
    TARGET_STRING("Value"), 0, 6, 0 },

  { 20, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant2"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 21, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant3"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 22, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Constant4"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 23, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Discrete-Time Integrator"),
    TARGET_STRING("gainval"), 0, 0, 0 },

  { 24, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Discrete-Time Integrator"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 25, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 26, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 27, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 28, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 29, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 30, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 31, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant12"),
    TARGET_STRING("Value"), 0, 8, 0 },

  { 32, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant13"),
    TARGET_STRING("Value"), 0, 6, 0 },

  { 33, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant2"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 34, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant3"),
    TARGET_STRING("Value"), 0, 9, 0 },

  { 35, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Constant4"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 36, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Discrete-Time Integrator"),
    TARGET_STRING("gainval"), 0, 0, 0 },

  { 37, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Discrete-Time Integrator"),
    TARGET_STRING("InitialCondition"), 0, 7, 0 },

  { 38, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 39, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 40, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 41, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/yhat"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 42, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 43, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant1"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 44, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant12"),
    TARGET_STRING("Value"), 0, 8, 0 },

  { 45, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant13"),
    TARGET_STRING("Value"), 0, 6, 0 },

  { 46, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant2"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 47, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant3"),
    TARGET_STRING("Value"), 0, 9, 0 },

  { 48, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Constant4"),
    TARGET_STRING("Value"), 0, 10, 0 },

  { 49, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Discrete-Time Integrator"),
    TARGET_STRING("gainval"), 0, 0, 0 },

  { 50, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Discrete-Time Integrator"),
    TARGET_STRING("InitialCondition"), 0, 7, 0 },

  { 51, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 52, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 53, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/u"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 54, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/Constant"),
    TARGET_STRING("Value"), 0, 7, 0 },

  { 55, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/Discrete-Time Integrator"),
    TARGET_STRING("gainval"), 0, 0, 0 },

  { 56, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/Discrete-Time Integrator"),
    TARGET_STRING("InitialCondition"), 0, 7, 0 },

  { 57, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/Saturation"),
    TARGET_STRING("UpperLimit"), 0, 0, 0 },

  { 58, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/Saturation"),
    TARGET_STRING("LowerLimit"), 0, 0, 0 },

  { 59, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/theta"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 60, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/P"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 61, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/err"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 62, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/theta"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 63, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/P"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 64, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/err"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 65, TARGET_STRING("SupervisoryController/SupervisoryController/wtMod/ywt"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 66, TARGET_STRING("SupervisoryController/SupervisoryController/wtMod/y_"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 67, TARGET_STRING("SupervisoryController/SupervisoryController/wtMod/r_"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 68, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 69, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/F_zero"),
    TARGET_STRING("Value"), 0, 11, 0 },

  { 70, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 71, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 72, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 73, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 74, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 75, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 76, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 77, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 78, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 79, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 12, 0 },

  { 80, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 81, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/F_zero"),
    TARGET_STRING("Value"), 0, 13, 0 },

  { 82, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 83, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 84, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 85, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 86, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 87, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 88, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 89, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 90, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 91, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/y.wt_zero"),
    TARGET_STRING("Value"), 0, 14, 0 },

  { 92, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ymax_zero"),
    TARGET_STRING("Value"), 0, 14, 0 },

  { 93, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/ymin_zero"),
    TARGET_STRING("Value"), 0, 14, 0 },

  { 94, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 95, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 96, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/F_zero"),
    TARGET_STRING("Value"), 0, 15, 0 },

  { 97, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 98, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 99, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 100, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 101, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 102, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 103, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 104, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 105, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 106, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/y.wt_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 107, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/ymax_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 108, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/ymin_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 109, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 110, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 111, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/F_zero"),
    TARGET_STRING("Value"), 0, 15, 0 },

  { 112, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 113, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 114, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 115, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 116, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 117, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 118, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 119, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 120, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 121, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/y.wt_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 122, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/ymax_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 123, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/ymin_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 124, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Constant1"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 125, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/E_zero"),
    TARGET_STRING("Value"), 0, 4, 0 },

  { 126, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/F_zero"),
    TARGET_STRING("Value"), 0, 15, 0 },

  { 127, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/G_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 128, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/S_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 129, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/du.wt_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 130, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/ecr.wt_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 131, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/ext.mv_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 132, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/md_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 133, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/mv.target_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 134, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/switch_zero"),
    TARGET_STRING("Value"), 0, 0, 0 },

  { 135, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/umin_zero"),
    TARGET_STRING("Value"), 0, 5, 0 },

  { 136, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/y.wt_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 137, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/ymax_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 138, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/ymin_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 139, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst1/Param Estimator (RLS)/Unit Delay3"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 140, TARGET_STRING("SupervisoryController/SupervisoryController/paramEst2/Param Estimator (RLS)/Unit Delay3"),
    TARGET_STRING("InitialCondition"), 0, 0, 0 },

  { 141, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 3, 0 },

  { 142, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 143, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/u_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 144, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 145, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/uref_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 146, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 11, 0 },

  { 147, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 148, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/LastPcov"),
    TARGET_STRING("InitialCondition"), 0, 17, 0 },

  { 149, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 18, 0 },

  { 150, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 151, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 17, 0 },

  { 152, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 19, 0 },

  { 153, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 17, 0 },

  { 154, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 20, 0 },

  { 155, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 14, 0 },

  { 156, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 157, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ext.mv_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 158, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/umin_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 159, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 160, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 13, 0 },

  { 161, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 162, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 21, 0 },

  { 163, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/MPC Controller1/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 164, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 22, 0 },

  { 165, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 23, 0 },

  { 166, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 22, 0 },

  { 167, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 24, 0 },

  { 168, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 169, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 170, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/ext.mv_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 171, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/umin_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 172, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 173, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 15, 0 },

  { 174, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 175, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 25, 0 },

  { 176, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/MPC Controller1/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 177, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 26, 0 },

  { 178, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 27, 0 },

  { 179, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 26, 0 },

  { 180, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 28, 0 },

  { 181, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/ym_zero"),
    TARGET_STRING("Value"), 0, 16, 0 },

  { 182, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 183, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/ext.mv_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 184, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/umin_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 185, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 186, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 15, 0 },

  { 187, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 188, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 25, 0 },

  { 189, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/MPC Controller1/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 190, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/G"),
    TARGET_STRING("Value"), 0, 26, 0 },

  { 191, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/H"),
    TARGET_STRING("Value"), 0, 27, 0 },

  { 192, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/P0"),
    TARGET_STRING("Value"), 0, 26, 0 },

  { 193, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/X0"),
    TARGET_STRING("Value"), 0, 28, 0 },

  { 194, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/ext.mv_scale"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 195, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/ext.mv_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 196, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/umin_scale1"),
    TARGET_STRING("Gain"), 0, 5, 0 },

  { 197, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/umin_scale4"),
    TARGET_STRING("Gain"), 0, 4, 0 },

  { 198, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/ymin_scale1"),
    TARGET_STRING("Gain"), 0, 15, 0 },

  { 199, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/ymin_scale2"),
    TARGET_STRING("Gain"), 0, 0, 0 },

  { 200, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/Memory"),
    TARGET_STRING("InitialCondition"), 1, 25, 0 },

  { 201, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/last_x"),
    TARGET_STRING("InitialCondition"), 0, 29, 0 },

  { 202, TARGET_STRING("SupervisoryController/SupervisoryController/mpcg/MPC Controller1/MPC/last_mv"),
    TARGET_STRING("InitialCondition"), 0, 5, 0 },

  { 203, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer"),
    TARGET_STRING("Ndis"), 2, 0, 0 },

  { 204, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 205, TARGET_STRING("SupervisoryController/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 206, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 207, TARGET_STRING("SupervisoryController/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 208, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 209, TARGET_STRING("SupervisoryController/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  { 210, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/isSqrtUsed"),
    TARGET_STRING("Value"), 1, 0, 0 },

  { 211, TARGET_STRING("SupervisoryController/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate/L*(y[k]-yhat[k|k-1])"),
    TARGET_STRING("InitialOutput"), 0, 0, 0 },

  {
    0, (nullptr), (nullptr), 0, 0, 0
  }
};

// Tunable variable parameters
static rtwCAPI_ModelParameters rtModelParameters[]{
  // addrMapIndex, varName, dataTypeIndex, dimIndex, fixPtIndex
  { 212, TARGET_STRING("nullEv"), 3, 0, 0 },

  { 213, TARGET_STRING("Aod"), 0, 32, 0 },

  { 214, TARGET_STRING("Aod1"), 0, 1, 0 },

  { 215, TARGET_STRING("Aod2"), 0, 1, 0 },

  { 216, TARGET_STRING("Aod3"), 0, 1, 0 },

  { 217, TARGET_STRING("Bod"), 0, 33, 0 },

  { 218, TARGET_STRING("Bod1"), 0, 2, 0 },

  { 219, TARGET_STRING("Bod2"), 0, 2, 0 },

  { 220, TARGET_STRING("Bod3"), 0, 2, 0 },

  { 221, TARGET_STRING("Cod"), 0, 34, 0 },

  { 222, TARGET_STRING("Cod1"), 0, 35, 0 },

  { 223, TARGET_STRING("Cod2"), 0, 36, 0 },

  { 224, TARGET_STRING("Cod3"), 0, 36, 0 },

  { 225, TARGET_STRING("Dmn"), 0, 1, 0 },

  { 226, TARGET_STRING("Dmn1"), 0, 6, 0 },

  { 227, TARGET_STRING("Dod"), 0, 1, 0 },

  { 228, TARGET_STRING("Dod1"), 0, 37, 0 },

  { 229, TARGET_STRING("Dod2"), 0, 38, 0 },

  { 230, TARGET_STRING("Dod3"), 0, 38, 0 },

  { 231, TARGET_STRING("beta"), 0, 0, 0 },

  { 232, TARGET_STRING("dt"), 0, 0, 0 },

  { 233, TARGET_STRING("lpfDen"), 0, 0, 0 },

  { 234, TARGET_STRING("lpfNum"), 0, 39, 0 },

  { 235, TARGET_STRING("mdlNum"), 0, 0, 0 },

  { 236, TARGET_STRING("uwt0"), 0, 4, 0 },

  { 237, TARGET_STRING("ywt0"), 0, 11, 0 },

  { 0, (nullptr), 0, 0, 0 }
};

#ifndef HOST_CAPI_BUILD

// Initialize Data Address
static void InitializeDataAddr(void* dataAddr[], SupervisoryController::P *rtP)
{
  dataAddr[0] = (void*) (&rtP->u_Y0);
  dataAddr[1] = (void*) (&rtP->yhat_Y0);
  dataAddr[2] = (void*) (&rtP->Constant1_Value_c);
  dataAddr[3] = (void*) (&rtP->Constant12_Value[0]);
  dataAddr[4] = (void*) (&rtP->Constant13_Value[0]);
  dataAddr[5] = (void*) (&rtP->Constant2_Value[0]);
  dataAddr[6] = (void*) (&rtP->Gain2_Gain);
  dataAddr[7] = (void*) (&rtP->Gain3_Gain);
  dataAddr[8] = (void*) (&rtP->excitation_Mean[0]);
  dataAddr[9] = (void*) (&rtP->excitation_StdDev[0]);
  dataAddr[10] = (void*) (&rtP->excitation_Seed[0]);
  dataAddr[11] = (void*) (&rtP->Saturation_UpperSat);
  dataAddr[12] = (void*) (&rtP->Saturation_LowerSat);
  dataAddr[13] = (void*) (&rtP->DiscreteFilter1_InitialStates);
  dataAddr[14] = (void*) (&rtP->u_Y0_b);
  dataAddr[15] = (void*) (&rtP->yhat_Y0_d);
  dataAddr[16] = (void*) (&rtP->Constant_Value);
  dataAddr[17] = (void*) (&rtP->Constant1_Value_e);
  dataAddr[18] = (void*) (&rtP->Constant12_Value_e[0]);
  dataAddr[19] = (void*) (&rtP->Constant13_Value_c[0]);
  dataAddr[20] = (void*) (&rtP->Constant2_Value_a);
  dataAddr[21] = (void*) (&rtP->Constant3_Value);
  dataAddr[22] = (void*) (&rtP->Constant4_Value[0]);
  dataAddr[23] = (void*) (&rtP->DiscreteTimeIntegrator_gainval);
  dataAddr[24] = (void*) (&rtP->DiscreteTimeIntegrator_IC);
  dataAddr[25] = (void*) (&rtP->Saturation_UpperSat_k);
  dataAddr[26] = (void*) (&rtP->Saturation_LowerSat_h);
  dataAddr[27] = (void*) (&rtP->u_Y0_n);
  dataAddr[28] = (void*) (&rtP->yhat_Y0_k);
  dataAddr[29] = (void*) (&rtP->Constant_Value_p[0]);
  dataAddr[30] = (void*) (&rtP->Constant1_Value_pe);
  dataAddr[31] = (void*) (&rtP->Constant12_Value_i[0]);
  dataAddr[32] = (void*) (&rtP->Constant13_Value_g[0]);
  dataAddr[33] = (void*) (&rtP->Constant2_Value_c[0]);
  dataAddr[34] = (void*) (&rtP->Constant3_Value_d[0]);
  dataAddr[35] = (void*) (&rtP->Constant4_Value_n[0]);
  dataAddr[36] = (void*) (&rtP->DiscreteTimeIntegrator_gainva_b);
  dataAddr[37] = (void*) (&rtP->DiscreteTimeIntegrator_IC_n[0]);
  dataAddr[38] = (void*) (&rtP->Saturation_UpperSat_h);
  dataAddr[39] = (void*) (&rtP->Saturation_LowerSat_o);
  dataAddr[40] = (void*) (&rtP->u_Y0_h);
  dataAddr[41] = (void*) (&rtP->yhat_Y0_f);
  dataAddr[42] = (void*) (&rtP->Constant_Value_e[0]);
  dataAddr[43] = (void*) (&rtP->Constant1_Value_n);
  dataAddr[44] = (void*) (&rtP->Constant12_Value_f[0]);
  dataAddr[45] = (void*) (&rtP->Constant13_Value_a[0]);
  dataAddr[46] = (void*) (&rtP->Constant2_Value_m[0]);
  dataAddr[47] = (void*) (&rtP->Constant3_Value_g[0]);
  dataAddr[48] = (void*) (&rtP->Constant4_Value_f[0]);
  dataAddr[49] = (void*) (&rtP->DiscreteTimeIntegrator_gainva_k);
  dataAddr[50] = (void*) (&rtP->DiscreteTimeIntegrator_IC_c[0]);
  dataAddr[51] = (void*) (&rtP->Saturation_UpperSat_c);
  dataAddr[52] = (void*) (&rtP->Saturation_LowerSat_b);
  dataAddr[53] = (void*) (&rtP->u_Y0_j);
  dataAddr[54] = (void*) (&rtP->Constant_Value_f[0]);
  dataAddr[55] = (void*) (&rtP->DiscreteTimeIntegrator_gainv_kx);
  dataAddr[56] = (void*) (&rtP->DiscreteTimeIntegrator_IC_o[0]);
  dataAddr[57] = (void*) (&rtP->Saturation_UpperSat_cb);
  dataAddr[58] = (void*) (&rtP->Saturation_LowerSat_p);
  dataAddr[59] = (void*) (&rtP->paramEst1_o.theta_Y0);
  dataAddr[60] = (void*) (&rtP->paramEst1_o.P_Y0);
  dataAddr[61] = (void*) (&rtP->paramEst1_o.err_Y0);
  dataAddr[62] = (void*) (&rtP->paramEst2.theta_Y0);
  dataAddr[63] = (void*) (&rtP->paramEst2.P_Y0);
  dataAddr[64] = (void*) (&rtP->paramEst2.err_Y0);
  dataAddr[65] = (void*) (&rtP->ywt_Y0);
  dataAddr[66] = (void*) (&rtP->y_Y0);
  dataAddr[67] = (void*) (&rtP->r_Y0);
  dataAddr[68] = (void*) (&rtP->E_zero_Value[0]);
  dataAddr[69] = (void*) (&rtP->F_zero_Value[0]);
  dataAddr[70] = (void*) (&rtP->G_zero_Value);
  dataAddr[71] = (void*) (&rtP->S_zero_Value);
  dataAddr[72] = (void*) (&rtP->duwt_zero_Value[0]);
  dataAddr[73] = (void*) (&rtP->ecrwt_zero_Value);
  dataAddr[74] = (void*) (&rtP->extmv_zero_Value[0]);
  dataAddr[75] = (void*) (&rtP->md_zero_Value);
  dataAddr[76] = (void*) (&rtP->mvtarget_zero_Value[0]);
  dataAddr[77] = (void*) (&rtP->switch_zero_Value);
  dataAddr[78] = (void*) (&rtP->umin_zero_Value[0]);
  dataAddr[79] = (void*) (&rtP->Constant1_Value[0]);
  dataAddr[80] = (void*) (&rtP->E_zero_Value_a[0]);
  dataAddr[81] = (void*) (&rtP->F_zero_Value_g[0]);
  dataAddr[82] = (void*) (&rtP->G_zero_Value_m);
  dataAddr[83] = (void*) (&rtP->S_zero_Value_g);
  dataAddr[84] = (void*) (&rtP->duwt_zero_Value_p[0]);
  dataAddr[85] = (void*) (&rtP->ecrwt_zero_Value_o);
  dataAddr[86] = (void*) (&rtP->extmv_zero_Value_e[0]);
  dataAddr[87] = (void*) (&rtP->md_zero_Value_p);
  dataAddr[88] = (void*) (&rtP->mvtarget_zero_Value_k[0]);
  dataAddr[89] = (void*) (&rtP->switch_zero_Value_e);
  dataAddr[90] = (void*) (&rtP->umin_zero_Value_d[0]);
  dataAddr[91] = (void*) (&rtP->ywt_zero_Value[0]);
  dataAddr[92] = (void*) (&rtP->ymax_zero_Value[0]);
  dataAddr[93] = (void*) (&rtP->ymin_zero_Value[0]);
  dataAddr[94] = (void*) (&rtP->Constant1_Value_j[0]);
  dataAddr[95] = (void*) (&rtP->E_zero_Value_b[0]);
  dataAddr[96] = (void*) (&rtP->F_zero_Value_o[0]);
  dataAddr[97] = (void*) (&rtP->G_zero_Value_n);
  dataAddr[98] = (void*) (&rtP->S_zero_Value_m);
  dataAddr[99] = (void*) (&rtP->duwt_zero_Value_l[0]);
  dataAddr[100] = (void*) (&rtP->ecrwt_zero_Value_e);
  dataAddr[101] = (void*) (&rtP->extmv_zero_Value_c[0]);
  dataAddr[102] = (void*) (&rtP->md_zero_Value_pu);
  dataAddr[103] = (void*) (&rtP->mvtarget_zero_Value_e[0]);
  dataAddr[104] = (void*) (&rtP->switch_zero_Value_i);
  dataAddr[105] = (void*) (&rtP->umin_zero_Value_e[0]);
  dataAddr[106] = (void*) (&rtP->ywt_zero_Value_n[0]);
  dataAddr[107] = (void*) (&rtP->ymax_zero_Value_g[0]);
  dataAddr[108] = (void*) (&rtP->ymin_zero_Value_g[0]);
  dataAddr[109] = (void*) (&rtP->Constant1_Value_p[0]);
  dataAddr[110] = (void*) (&rtP->E_zero_Value_j[0]);
  dataAddr[111] = (void*) (&rtP->F_zero_Value_n[0]);
  dataAddr[112] = (void*) (&rtP->G_zero_Value_j);
  dataAddr[113] = (void*) (&rtP->S_zero_Value_i);
  dataAddr[114] = (void*) (&rtP->duwt_zero_Value_a[0]);
  dataAddr[115] = (void*) (&rtP->ecrwt_zero_Value_j);
  dataAddr[116] = (void*) (&rtP->extmv_zero_Value_m[0]);
  dataAddr[117] = (void*) (&rtP->md_zero_Value_o);
  dataAddr[118] = (void*) (&rtP->mvtarget_zero_Value_d[0]);
  dataAddr[119] = (void*) (&rtP->switch_zero_Value_k);
  dataAddr[120] = (void*) (&rtP->umin_zero_Value_b[0]);
  dataAddr[121] = (void*) (&rtP->ywt_zero_Value_g[0]);
  dataAddr[122] = (void*) (&rtP->ymax_zero_Value_d[0]);
  dataAddr[123] = (void*) (&rtP->ymin_zero_Value_e[0]);
  dataAddr[124] = (void*) (&rtP->Constant1_Value_h[0]);
  dataAddr[125] = (void*) (&rtP->E_zero_Value_l[0]);
  dataAddr[126] = (void*) (&rtP->F_zero_Value_l[0]);
  dataAddr[127] = (void*) (&rtP->G_zero_Value_a);
  dataAddr[128] = (void*) (&rtP->S_zero_Value_h);
  dataAddr[129] = (void*) (&rtP->duwt_zero_Value_lq[0]);
  dataAddr[130] = (void*) (&rtP->ecrwt_zero_Value_o0);
  dataAddr[131] = (void*) (&rtP->extmv_zero_Value_j[0]);
  dataAddr[132] = (void*) (&rtP->md_zero_Value_c);
  dataAddr[133] = (void*) (&rtP->mvtarget_zero_Value_n[0]);
  dataAddr[134] = (void*) (&rtP->switch_zero_Value_j);
  dataAddr[135] = (void*) (&rtP->umin_zero_Value_k[0]);
  dataAddr[136] = (void*) (&rtP->ywt_zero_Value_d[0]);
  dataAddr[137] = (void*) (&rtP->ymax_zero_Value_k[0]);
  dataAddr[138] = (void*) (&rtP->ymin_zero_Value_a[0]);
  dataAddr[139] = (void*) (&rtP->paramEst1_o.UnitDelay3_InitialCondition);
  dataAddr[140] = (void*) (&rtP->paramEst2.UnitDelay3_InitialCondition);
  dataAddr[141] = (void*) (&rtP->ym_zero_Value[0]);
  dataAddr[142] = (void*) (&rtP->extmv_scale_Gain[0]);
  dataAddr[143] = (void*) (&rtP->u_scale_Gain[0]);
  dataAddr[144] = (void*) (&rtP->umin_scale4_Gain[0]);
  dataAddr[145] = (void*) (&rtP->uref_scale_Gain[0]);
  dataAddr[146] = (void*) (&rtP->ymin_scale1_Gain[0]);
  dataAddr[147] = (void*) (&rtP->ymin_scale2_Gain);
  dataAddr[148] = (void*) (&rtP->LastPcov_InitialCondition[0]);
  dataAddr[149] = (void*) (&rtP->Memory_InitialCondition[0]);
  dataAddr[150] = (void*) (&rtP->last_mv_InitialCondition[0]);
  dataAddr[151] = (void*) (&rtP->G_Value[0]);
  dataAddr[152] = (void*) (&rtP->H_Value[0]);
  dataAddr[153] = (void*) (&rtP->P0_Value[0]);
  dataAddr[154] = (void*) (&rtP->X0_Value[0]);
  dataAddr[155] = (void*) (&rtP->ym_zero_Value_c[0]);
  dataAddr[156] = (void*) (&rtP->extmv_scale_Gain_e[0]);
  dataAddr[157] = (void*) (&rtP->extmv_scale1_Gain[0]);
  dataAddr[158] = (void*) (&rtP->umin_scale1_Gain[0]);
  dataAddr[159] = (void*) (&rtP->umin_scale4_Gain_p[0]);
  dataAddr[160] = (void*) (&rtP->ymin_scale1_Gain_j[0]);
  dataAddr[161] = (void*) (&rtP->ymin_scale2_Gain_f);
  dataAddr[162] = (void*) (&rtP->Memory_InitialCondition_f[0]);
  dataAddr[163] = (void*) (&rtP->last_mv_InitialCondition_f[0]);
  dataAddr[164] = (void*) (&rtP->G_Value_a[0]);
  dataAddr[165] = (void*) (&rtP->H_Value_o[0]);
  dataAddr[166] = (void*) (&rtP->P0_Value_a[0]);
  dataAddr[167] = (void*) (&rtP->X0_Value_f[0]);
  dataAddr[168] = (void*) (&rtP->ym_zero_Value_l[0]);
  dataAddr[169] = (void*) (&rtP->extmv_scale_Gain_g[0]);
  dataAddr[170] = (void*) (&rtP->extmv_scale1_Gain_b[0]);
  dataAddr[171] = (void*) (&rtP->umin_scale1_Gain_p[0]);
  dataAddr[172] = (void*) (&rtP->umin_scale4_Gain_g[0]);
  dataAddr[173] = (void*) (&rtP->ymin_scale1_Gain_f[0]);
  dataAddr[174] = (void*) (&rtP->ymin_scale2_Gain_g);
  dataAddr[175] = (void*) (&rtP->Memory_InitialCondition_j[0]);
  dataAddr[176] = (void*) (&rtP->last_mv_InitialCondition_b[0]);
  dataAddr[177] = (void*) (&rtP->G_Value_g[0]);
  dataAddr[178] = (void*) (&rtP->H_Value_k[0]);
  dataAddr[179] = (void*) (&rtP->P0_Value_c[0]);
  dataAddr[180] = (void*) (&rtP->X0_Value_k[0]);
  dataAddr[181] = (void*) (&rtP->ym_zero_Value_d[0]);
  dataAddr[182] = (void*) (&rtP->extmv_scale_Gain_h[0]);
  dataAddr[183] = (void*) (&rtP->extmv_scale1_Gain_e[0]);
  dataAddr[184] = (void*) (&rtP->umin_scale1_Gain_g[0]);
  dataAddr[185] = (void*) (&rtP->umin_scale4_Gain_f[0]);
  dataAddr[186] = (void*) (&rtP->ymin_scale1_Gain_e[0]);
  dataAddr[187] = (void*) (&rtP->ymin_scale2_Gain_e);
  dataAddr[188] = (void*) (&rtP->Memory_InitialCondition_b[0]);
  dataAddr[189] = (void*) (&rtP->last_mv_InitialCondition_i[0]);
  dataAddr[190] = (void*) (&rtP->G_Value_h[0]);
  dataAddr[191] = (void*) (&rtP->H_Value_oa[0]);
  dataAddr[192] = (void*) (&rtP->P0_Value_m[0]);
  dataAddr[193] = (void*) (&rtP->X0_Value_a[0]);
  dataAddr[194] = (void*) (&rtP->extmv_scale_Gain_b[0]);
  dataAddr[195] = (void*) (&rtP->extmv_scale1_Gain_j[0]);
  dataAddr[196] = (void*) (&rtP->umin_scale1_Gain_l[0]);
  dataAddr[197] = (void*) (&rtP->umin_scale4_Gain_m[0]);
  dataAddr[198] = (void*) (&rtP->ymin_scale1_Gain_o[0]);
  dataAddr[199] = (void*) (&rtP->ymin_scale2_Gain_n);
  dataAddr[200] = (void*) (&rtP->Memory_InitialCondition_o[0]);
  dataAddr[201] = (void*) (&rtP->last_x_InitialCondition[0]);
  dataAddr[202] = (void*) (&rtP->last_mv_InitialCondition_g[0]);
  dataAddr[203] = (void*) (&rtP->FixedHorizonOptimizer_Ndis);
  dataAddr[204] = (void*) (&rtP->isSqrtUsed_Value);
  dataAddr[205] = (void*) (&rtP->Lykyhatkk1_Y0);
  dataAddr[206] = (void*) (&rtP->isSqrtUsed_Value_d);
  dataAddr[207] = (void*) (&rtP->Lykyhatkk1_Y0_c);
  dataAddr[208] = (void*) (&rtP->isSqrtUsed_Value_a);
  dataAddr[209] = (void*) (&rtP->MeasurementUpdate_j.Lykyhatkk1_Y0);
  dataAddr[210] = (void*) (&rtP->isSqrtUsed_Value_p);
  dataAddr[211] = (void*) (&rtP->MeasurementUpdate_c.Lykyhatkk1_Y0);
  dataAddr[212] = (void*) (&rtP->nullEv);
  dataAddr[213] = (void*) (&rtP->Aod[0]);
  dataAddr[214] = (void*) (&rtP->Aod1[0]);
  dataAddr[215] = (void*) (&rtP->Aod2[0]);
  dataAddr[216] = (void*) (&rtP->Aod3[0]);
  dataAddr[217] = (void*) (&rtP->Bod[0]);
  dataAddr[218] = (void*) (&rtP->Bod1[0]);
  dataAddr[219] = (void*) (&rtP->Bod2[0]);
  dataAddr[220] = (void*) (&rtP->Bod3[0]);
  dataAddr[221] = (void*) (&rtP->Cod[0]);
  dataAddr[222] = (void*) (&rtP->Cod1[0]);
  dataAddr[223] = (void*) (&rtP->Cod2[0]);
  dataAddr[224] = (void*) (&rtP->Cod3[0]);
  dataAddr[225] = (void*) (&rtP->Dmn[0]);
  dataAddr[226] = (void*) (&rtP->Dmn1[0]);
  dataAddr[227] = (void*) (&rtP->Dod[0]);
  dataAddr[228] = (void*) (&rtP->Dod1[0]);
  dataAddr[229] = (void*) (&rtP->Dod2[0]);
  dataAddr[230] = (void*) (&rtP->Dod3[0]);
  dataAddr[231] = (void*) (&rtP->beta);
  dataAddr[232] = (void*) (&rtP->dt);
  dataAddr[233] = (void*) (&rtP->lpfDen);
  dataAddr[234] = (void*) (&rtP->lpfNum[0]);
  dataAddr[235] = (void*) (&rtP->mdlNum);
  dataAddr[236] = (void*) (&rtP->uwt0[0]);
  dataAddr[237] = (void*) (&rtP->ywt0[0]);
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
  loggingPtrs[115] = (nullptr);
  loggingPtrs[116] = (nullptr);
  loggingPtrs[117] = (nullptr);
  loggingPtrs[118] = (nullptr);
  loggingPtrs[119] = (nullptr);
  loggingPtrs[120] = (nullptr);
  loggingPtrs[121] = (nullptr);
  loggingPtrs[122] = (nullptr);
  loggingPtrs[123] = (nullptr);
  loggingPtrs[124] = (nullptr);
  loggingPtrs[125] = (nullptr);
  loggingPtrs[126] = (nullptr);
  loggingPtrs[127] = (nullptr);
  loggingPtrs[128] = (nullptr);
  loggingPtrs[129] = (nullptr);
  loggingPtrs[130] = (nullptr);
  loggingPtrs[131] = (nullptr);
  loggingPtrs[132] = (nullptr);
  loggingPtrs[133] = (nullptr);
  loggingPtrs[134] = (nullptr);
  loggingPtrs[135] = (nullptr);
  loggingPtrs[136] = (nullptr);
  loggingPtrs[137] = (nullptr);
  loggingPtrs[138] = (nullptr);
  loggingPtrs[139] = (nullptr);
  loggingPtrs[140] = (nullptr);
  loggingPtrs[141] = (nullptr);
  loggingPtrs[142] = (nullptr);
  loggingPtrs[143] = (nullptr);
  loggingPtrs[144] = (nullptr);
  loggingPtrs[145] = (nullptr);
  loggingPtrs[146] = (nullptr);
  loggingPtrs[147] = (nullptr);
  loggingPtrs[148] = (nullptr);
  loggingPtrs[149] = (nullptr);
  loggingPtrs[150] = (nullptr);
  loggingPtrs[151] = (nullptr);
  loggingPtrs[152] = (nullptr);
  loggingPtrs[153] = (nullptr);
  loggingPtrs[154] = (nullptr);
  loggingPtrs[155] = (nullptr);
  loggingPtrs[156] = (nullptr);
  loggingPtrs[157] = (nullptr);
  loggingPtrs[158] = (nullptr);
  loggingPtrs[159] = (nullptr);
  loggingPtrs[160] = (nullptr);
  loggingPtrs[161] = (nullptr);
  loggingPtrs[162] = (nullptr);
  loggingPtrs[163] = (nullptr);
  loggingPtrs[164] = (nullptr);
  loggingPtrs[165] = (nullptr);
  loggingPtrs[166] = (nullptr);
  loggingPtrs[167] = (nullptr);
  loggingPtrs[168] = (nullptr);
  loggingPtrs[169] = (nullptr);
  loggingPtrs[170] = (nullptr);
  loggingPtrs[171] = (nullptr);
  loggingPtrs[172] = (nullptr);
  loggingPtrs[173] = (nullptr);
  loggingPtrs[174] = (nullptr);
  loggingPtrs[175] = (nullptr);
  loggingPtrs[176] = (nullptr);
  loggingPtrs[177] = (nullptr);
  loggingPtrs[178] = (nullptr);
  loggingPtrs[179] = (nullptr);
  loggingPtrs[180] = (nullptr);
  loggingPtrs[181] = (nullptr);
  loggingPtrs[182] = (nullptr);
  loggingPtrs[183] = (nullptr);
  loggingPtrs[184] = (nullptr);
  loggingPtrs[185] = (nullptr);
  loggingPtrs[186] = (nullptr);
  loggingPtrs[187] = (nullptr);
  loggingPtrs[188] = (nullptr);
  loggingPtrs[189] = (nullptr);
  loggingPtrs[190] = (nullptr);
  loggingPtrs[191] = (nullptr);
  loggingPtrs[192] = (nullptr);
  loggingPtrs[193] = (nullptr);
  loggingPtrs[194] = (nullptr);
  loggingPtrs[195] = (nullptr);
  loggingPtrs[196] = (nullptr);
  loggingPtrs[197] = (nullptr);
  loggingPtrs[198] = (nullptr);
  loggingPtrs[199] = (nullptr);
  loggingPtrs[200] = (nullptr);
  loggingPtrs[201] = (nullptr);
  loggingPtrs[202] = (nullptr);
  loggingPtrs[203] = (nullptr);
  loggingPtrs[204] = (nullptr);
  loggingPtrs[205] = (nullptr);
  loggingPtrs[206] = (nullptr);
  loggingPtrs[207] = (nullptr);
  loggingPtrs[208] = (nullptr);
  loggingPtrs[209] = (nullptr);
  loggingPtrs[210] = (nullptr);
  loggingPtrs[211] = (nullptr);
  loggingPtrs[212] = (nullptr);
  loggingPtrs[213] = (nullptr);
  loggingPtrs[214] = (nullptr);
  loggingPtrs[215] = (nullptr);
  loggingPtrs[216] = (nullptr);
  loggingPtrs[217] = (nullptr);
  loggingPtrs[218] = (nullptr);
  loggingPtrs[219] = (nullptr);
  loggingPtrs[220] = (nullptr);
  loggingPtrs[221] = (nullptr);
  loggingPtrs[222] = (nullptr);
  loggingPtrs[223] = (nullptr);
  loggingPtrs[224] = (nullptr);
  loggingPtrs[225] = (nullptr);
  loggingPtrs[226] = (nullptr);
  loggingPtrs[227] = (nullptr);
  loggingPtrs[228] = (nullptr);
  loggingPtrs[229] = (nullptr);
  loggingPtrs[230] = (nullptr);
  loggingPtrs[231] = (nullptr);
  loggingPtrs[232] = (nullptr);
  loggingPtrs[233] = (nullptr);
  loggingPtrs[234] = (nullptr);
  loggingPtrs[235] = (nullptr);
  loggingPtrs[236] = (nullptr);
  loggingPtrs[237] = (nullptr);
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

  { "r", rt_offsetof(event_bus, r), 0, 30, 0 },

  { "preT", rt_offsetof(event_bus, preT), 0, 31, 0 },

  { "moveT", rt_offsetof(event_bus, moveT), 0, 31, 0 },

  { "postT", rt_offsetof(event_bus, postT), 0, 31, 0 }
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

  { rtwCAPI_MATRIX_COL_MAJOR, 16, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 18, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 20, 2, 0 },

  { rtwCAPI_VECTOR, 22, 2, 0 },

  { rtwCAPI_VECTOR, 24, 2, 0 },

  { rtwCAPI_VECTOR, 26, 2, 0 },

  { rtwCAPI_VECTOR, 28, 2, 0 },

  { rtwCAPI_VECTOR, 30, 2, 0 },

  { rtwCAPI_VECTOR, 32, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 34, 2, 0 },

  { rtwCAPI_VECTOR, 36, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 38, 2, 0 },

  { rtwCAPI_VECTOR, 40, 2, 0 },

  { rtwCAPI_VECTOR, 42, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 44, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 46, 2, 0 },

  { rtwCAPI_VECTOR, 48, 2, 0 },

  { rtwCAPI_VECTOR, 50, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 52, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 54, 2, 0 },

  { rtwCAPI_VECTOR, 56, 2, 0 },

  { rtwCAPI_VECTOR, 58, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 6, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 0, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 60, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 62, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 64, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 66, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 68, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 70, 2, 0 },

  { rtwCAPI_MATRIX_COL_MAJOR, 72, 2, 0 },

  { rtwCAPI_VECTOR, 74, 2, 0 }
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
  2,                                   // 14
  1,                                   // 15
  3,                                   // 16
  2,                                   // 17
  2,                                   // 18
  2,                                   // 19
  2,                                   // 20
  3,                                   // 21
  1,                                   // 22
  6,                                   // 23
  12,                                  // 24
  1,                                   // 25
  1,                                   // 26
  4,                                   // 27
  4,                                   // 28
  1,                                   // 29
  1,                                   // 30
  5,                                   // 31
  5,                                   // 32
  1,                                   // 33
  18,                                  // 34
  18,                                  // 35
  246,                                 // 36
  1,                                   // 37
  6,                                   // 38
  18,                                  // 39
  18,                                  // 40
  1,                                   // 41
  166,                                 // 42
  1,                                   // 43
  7,                                   // 44
  7,                                   // 45
  3,                                   // 46
  7,                                   // 47
  7,                                   // 48
  1,                                   // 49
  206,                                 // 50
  1,                                   // 51
  8,                                   // 52
  8,                                   // 53
  3,                                   // 54
  8,                                   // 55
  8,                                   // 56
  1,                                   // 57
  11,                                  // 58
  1,                                   // 59
  12,                                  // 60
  12,                                  // 61
  12,                                  // 62
  6,                                   // 63
  6,                                   // 64
  12,                                  // 65
  4,                                   // 66
  6,                                   // 67
  5,                                   // 68
  6,                                   // 69
  4,                                   // 70
  3,                                   // 71
  5,                                   // 72
  3,                                   // 73
  1,                                   // 74
  60                                   // 75
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

  { rtBlockParameters, 212,
    rtModelParameters, 26 },

  { (nullptr), 0 },

  { rtDataTypeMap, rtDimensionMap, rtFixPtMap,
    rtElementMap, rtSampleTimeMap, rtDimensionArray },
  "float",

  { 3487440580U,
    566627409U,
    1180415088U,
    223588936U },
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
