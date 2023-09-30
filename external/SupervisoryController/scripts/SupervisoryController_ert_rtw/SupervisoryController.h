//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2506
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sat Sep 30 06:00:49 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#ifndef RTW_HEADER_SupervisoryController_h_
#define RTW_HEADER_SupervisoryController_h_
#include <stdio.h>
#include "rtwtypes.h"
#include "rtw_modelmap.h"
#include <stddef.h>
#include "zero_crossing_types.h"

// Macros for accessing real-time model data structure
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm)         ((rtm)->DataMapInfo)
#endif

#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val)    ((rtm)->DataMapInfo = (val))
#endif

#ifndef DEFINED_TYPEDEF_FOR_event_bus_
#define DEFINED_TYPEDEF_FOR_event_bus_

struct event_bus
{
  real_T r[6];
  real_T preT;
  real_T moveT;
  real_T postT;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_mdl_bus_
#define DEFINED_TYPEDEF_FOR_mdl_bus_

struct mdl_bus
{
  real_T A[36];
  real_T B[18];
  real_T C[36];
  real_T D[18];
  real_T U[3];
  real_T Y[6];
  real_T X[6];
  real_T DX[6];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_mdl0_bus_
#define DEFINED_TYPEDEF_FOR_mdl0_bus_

struct mdl0_bus
{
  real_T A;
  real_T B[3];
  real_T C[3];
  real_T D[9];
  real_T U[3];
  real_T Y[3];
  real_T X;
  real_T DX;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_mdl1_bus_
#define DEFINED_TYPEDEF_FOR_mdl1_bus_

struct mdl1_bus
{
  real_T A[4];
  real_T B[6];
  real_T C[6];
  real_T D[9];
  real_T U[3];
  real_T Y[3];
  real_T X[2];
  real_T DX[2];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_mdl2_bus_
#define DEFINED_TYPEDEF_FOR_mdl2_bus_

struct mdl2_bus
{
  real_T A[4];
  real_T B[6];
  real_T C[6];
  real_T D[9];
  real_T U[3];
  real_T Y[3];
  real_T X[2];
  real_T DX[2];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_
#define DEFINED_TYPEDEF_FOR_struct_WTmPWsEMvOzNnnAVv5fQNC_

struct struct_WTmPWsEMvOzNnnAVv5fQNC
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  boolean_T UseWarmStart;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_
#define DEFINED_TYPEDEF_FOR_struct_WHjMt45Sk148iktWsfFxl_

struct struct_WHjMt45Sk148iktWsfFxl
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  real_T OptimalityTolerance;
  real_T ComplementarityTolerance;
  real_T StepTolerance;
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_
#define DEFINED_TYPEDEF_FOR_struct_lnQ9KXdSZFplhcBp5LBCc_

struct struct_lnQ9KXdSZFplhcBp5LBCc
{
  int32_T MaxIterations;
  real_T ConstraintTolerance;
  real_T DiscreteConstraintTolerance;
  boolean_T RoundingAtRootNode;
  int32_T MaxPendingNodes;
};

#endif

// Custom Type definition for Chart: '<Root>/SupervisoryController'
#ifndef struct_siswYcTR8LLamuD4YWmtXHC
#define struct_siswYcTR8LLamuD4YWmtXHC

struct siswYcTR8LLamuD4YWmtXHC
{
  real_T breaks[6];
  real_T coefs[15];
};

#endif                                 // struct_siswYcTR8LLamuD4YWmtXHC

#ifndef struct_cell_wrap_9
#define struct_cell_wrap_9

struct cell_wrap_9
{
  real_T f1[9];
};

#endif                                 // struct_cell_wrap_9

#ifndef SS_UINT64
#define SS_UINT64                      26
#endif

#ifndef SS_INT64
#define SS_INT64                       27
#endif

// Function to get C API Model Mapping Static Info
extern const rtwCAPI_ModelMappingStaticInfo*
  SupervisoryController_GetCAPIStaticMap(void);
extern "C"
{
  static real_T rtGetInf(void);
  static real32_T rtGetInfF(void);
  static real_T rtGetMinusInf(void);
  static real32_T rtGetMinusInfF(void);
}                                      // extern "C"

extern "C"
{
  static real_T rtGetNaN(void);
  static real32_T rtGetNaNF(void);
}                                      // extern "C"

#define NOT_USING_NONFINITE_LITERALS   1

extern "C"
{
  extern real_T rtInf;
  extern real_T rtMinusInf;
  extern real_T rtNaN;
  extern real32_T rtInfF;
  extern real32_T rtMinusInfF;
  extern real32_T rtNaNF;
  static void rt_InitInfAndNaN(size_t realSize);
  static boolean_T rtIsInf(real_T value);
  static boolean_T rtIsInfF(real32_T value);
  static boolean_T rtIsNaN(real_T value);
  static boolean_T rtIsNaNF(real32_T value);
  struct BigEndianIEEEDouble {
    struct {
      uint32_T wordH;
      uint32_T wordL;
    } words;
  };

  struct LittleEndianIEEEDouble {
    struct {
      uint32_T wordL;
      uint32_T wordH;
    } words;
  };

  struct IEEESingle {
    union {
      real32_T wordLreal;
      uint32_T wordLuint;
    } wordL;
  };
}                                      // extern "C"

// Class declaration for model SupervisoryController
class SupervisoryController final
{
  // public data and function members
 public:
  // Block signals and states (default storage) for system '<S1>/paramEst1'
  struct DW_paramEst1 {
    real_T Delay1_DSTATE[144];         // '<S324>/Delay1'
    real_T UnitDelay3_DSTATE[3];       // '<S322>/Unit Delay3'
    real_T Delay_DSTATE[12];           // '<S324>/Delay'
    boolean_T icLoad;                  // '<S324>/Delay1'
    boolean_T icLoad_n;                // '<S324>/Delay'
  };

  // Zero-crossing (trigger) state for system '<S1>/paramEst1'
  struct ZCE_paramEst1 {
    ZCSigState Delay1_Reset_ZCE;       // '<S324>/Delay1'
    ZCSigState Delay_Reset_ZCE;        // '<S324>/Delay'
  };

  // Block signals and states (default storage) for system '<S202>/MeasurementUpdate' 
  struct DW_MeasurementUpdate {
    boolean_T MeasurementUpdate_MODE;  // '<S202>/MeasurementUpdate'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    DW_MeasurementUpdate MeasurementUpdate_c;// '<S272>/MeasurementUpdate'
    DW_MeasurementUpdate MeasurementUpdate_j;// '<S202>/MeasurementUpdate'
    DW_paramEst1 paramEst2;            // '<S1>/paramEst2'
    DW_paramEst1 paramEst1_o;          // '<S1>/paramEst1'
    real_T r[6];                       // '<Root>/SupervisoryController'
    real_T y_l[6];                     // '<Root>/SupervisoryController'
    real_T ymax[6];                    // '<Root>/SupervisoryController'
    real_T ywt_o[6];                   // '<Root>/SupervisoryController'
    real_T y0_p[6];                    // '<Root>/SupervisoryController'
    real_T x0[6];                      // '<Root>/SupervisoryController'
    real_T u0_e[3];                    // '<Root>/SupervisoryController'
    real_T umax[3];                    // '<Root>/SupervisoryController'
    real_T uwt[3];                     // '<Root>/SupervisoryController'
    real_T theta_m[24];                // '<Root>/SupervisoryController'
    real_T thetaSgn_k[24];             // '<Root>/SupervisoryController'
    real_T r__p[3];                    // '<Root>/SupervisoryController'
    real_T y__g[3];                    // '<Root>/SupervisoryController'
    real_T y0_j[3];                    // '<Root>/SupervisoryController'
    real_T x0_l[2];                    // '<Root>/SupervisoryController'
    real_T u0_f[3];                    // '<Root>/SupervisoryController'
    real_T umax_i[3];                  // '<Root>/SupervisoryController'
    real_T uwt_m[3];                   // '<Root>/SupervisoryController'
    real_T Product3_a[8];              // '<S226>/Product3'
    real_T Product3_c[7];              // '<S156>/Product3'
    real_T DiscreteTimeIntegrator_DSTATE[2];// '<S6>/Discrete-Time Integrator'
    real_T last_mv_DSTATE[3];          // '<S300>/last_mv'
    real_T Delay_DSTATE[6];            // '<S9>/Delay'
    real_T MemoryX_DSTATE[8];          // '<S253>/MemoryX'
    real_T MemoryP_DSTATE[64];         // '<S253>/MemoryP'
    real_T DiscreteTimeIntegrator_DSTATE_m[2];// '<S4>/Discrete-Time Integrator' 
    real_T last_mv_DSTATE_i[3];        // '<S161>/last_mv'
    real_T MemoryX_DSTATE_c[8];        // '<S183>/MemoryX'
    real_T MemoryP_DSTATE_h[64];       // '<S183>/MemoryP'
    real_T last_mv_DSTATE_n[3];        // '<S91>/last_mv'
    real_T MemoryX_DSTATE_l[7];        // '<S113>/MemoryX'
    real_T MemoryP_DSTATE_e[49];       // '<S113>/MemoryP'
    real_T MemoryX_DSTATE_d[18];       // '<S43>/MemoryX'
    real_T MemoryP_DSTATE_h4[324];     // '<S43>/MemoryP'
    real_T traj[14400];                // '<Root>/SupervisoryController'
    real_T P0_1[144];                  // '<Root>/SupervisoryController'
    real_T P0_2[144];                  // '<Root>/SupervisoryController'
    real_T theta0_1[12];               // '<Root>/SupervisoryController'
    real_T theta0_2[12];               // '<Root>/SupervisoryController'
    real_T thetaSgn[24];               // '<Root>/SupervisoryController'
    real_T last_x_PreviousInput[11];   // '<S300>/last_x'
    real_T DiscreteFilter1_tmp[3];     // '<S2>/Discrete Filter1'
    real_T e_data[1199];
    real_T t_data[1199];
    real_T tmp_data[1199];
    real_T exVal;                      // '<Root>/SupervisoryController'
    real_T SFunction_o50;              // '<Root>/SupervisoryController'
    real_T SFunction_o51;              // '<Root>/SupervisoryController'
    real_T SFunction_o52;              // '<Root>/SupervisoryController'
    real_T SFunction_o53;              // '<Root>/SupervisoryController'
    real_T SFunction_o54;              // '<Root>/SupervisoryController'
    real_T SFunction_o55;              // '<Root>/SupervisoryController'
    real_T SFunction_o56;              // '<Root>/SupervisoryController'
    real_T DiscreteTimeIntegrator_DSTATE_j;// '<S3>/Discrete-Time Integrator'
    real_T sigPrev;                    // '<Root>/SupervisoryController'
    uint16_T waypt;                    // '<Root>/SupervisoryController'
    uint16_T trajSize;                 // '<Root>/SupervisoryController'
    int8_T DiscreteTimeIntegrator_PrevRese;// '<S6>/Discrete-Time Integrator'
    int8_T DiscreteTimeIntegrator_PrevRe_f;// '<S4>/Discrete-Time Integrator'
    int8_T DiscreteTimeIntegrator_PrevRe_b;// '<S3>/Discrete-Time Integrator'
    uint8_T is_EventHandler;           // '<Root>/SupervisoryController'
    uint8_T is_active_c6_SupervisoryControl;// '<Root>/SupervisoryController'
    boolean_T enAdapt_[6];             // '<Root>/SupervisoryController'
    boolean_T Memory_PreviousInput[206];// '<S300>/Memory'
    boolean_T Memory_PreviousInput_c[206];// '<S161>/Memory'
    boolean_T Memory_PreviousInput_d[166];// '<S91>/Memory'
    boolean_T integratorRST_oj;        // '<Root>/SupervisoryController'
    boolean_T evDone;                  // '<Root>/SupervisoryController'
    boolean_T icLoad_a;                // '<S183>/MemoryX'
    boolean_T icLoad_p;                // '<S183>/MemoryP'
    boolean_T icLoad_n;                // '<S113>/MemoryX'
    boolean_T icLoad_h;                // '<S113>/MemoryP'
    boolean_T MeasurementUpdate_MODE;  // '<S132>/MeasurementUpdate'
    boolean_T MeasurementUpdate_MODE_b;// '<S62>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state
  struct PrevZCX {
    ZCSigState SupervisoryController_Trig_ZCE;// '<Root>/SupervisoryController'
    ZCSigState MemoryX_Reset_ZCE;      // '<S253>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S253>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_g;    // '<S183>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_i;    // '<S183>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_l;    // '<S113>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_n;    // '<S113>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_o;    // '<S43>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_b;    // '<S43>/MemoryP'
    ZCE_paramEst1 paramEst2;           // '<S1>/paramEst2'
    ZCE_paramEst1 paramEst1_o;         // '<S1>/paramEst1'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T y[6];                       // '<Root>/y'
    real_T ymax[6];                    // '<Root>/ymax'
    real_T y0[6];                      // '<Root>/y0'
    real_T x0[6];                      // '<Root>/x0'
    real_T u0[3];                      // '<Root>/u0'
    real_T umax[3];                    // '<Root>/umax'
    real_T uwt[3];                     // '<Root>/uwt'
    boolean_T yo[6];                   // '<Root>/zeroCross'
    boolean_T enAdapt[6];              // '<Root>/enAdapt'
    real_T excitation;                 // '<Root>/excitation'
    real_T dPmod_;                     // '<Root>/dPmod_'
    real_T p_;                         // '<Root>/p_'
    real_T lambda;                     // '<Root>/lambda'
    real_T k_2;                        // '<Root>/k_2'
    event_bus nextEv;                  // '<Root>/nextEv'
    boolean_T iRST;                    // '<Root>/iRST'
    real_T measAvail;                  // '<Root>/measAvail'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T u[3];                       // '<Root>/u'
    real_T ywt[6];                     // '<Root>/ywt'
    real_T yhat[6];                    // '<Root>/yhat'
    real_T currTraj[6];                // '<Root>/currTraj'
    event_bus currEv;                  // '<Root>/currEv'
    real_T theta[24];                  // '<Root>/theta'
    real_T prmErr[6];                  // '<Root>/prmErr'
    real_T P_p[576];                   // '<Root>/P'
    boolean_T requestEvent;            // '<Root>/requestEvent'
    real_T sig;                        // '<Root>/sig'
  };

  // Parameters for system: '<S1>/paramEst1'
  struct P_paramEst1 {
    real_T theta_Y0;                   // Computed Parameter: theta_Y0
                                          //  Referenced by: '<S7>/theta'

    real_T P_Y0;                       // Computed Parameter: P_Y0
                                          //  Referenced by: '<S7>/P'

    real_T err_Y0;                     // Computed Parameter: err_Y0
                                          //  Referenced by: '<S7>/err'

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S322>/Unit Delay3'

  };

  // Parameters for system: '<S202>/MeasurementUpdate'
  struct P_MeasurementUpdate {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S226>/L*(y[k]-yhat[k|k-1])'

  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod[144];                   // Variable: Aod
                                          //  Referenced by: '<S12>/MATLAB Function'

    real_T Aod1[36];                   // Variable: Aod1
                                          //  Referenced by: '<S90>/MATLAB Function'

    real_T Aod2[36];                   // Variable: Aod2
                                          //  Referenced by: '<S160>/MATLAB Function'

    real_T Aod3[36];                   // Variable: Aod3
                                          //  Referenced by: '<S230>/MATLAB Function'

    real_T Bod[72];                    // Variable: Bod
                                          //  Referenced by: '<S12>/MATLAB Function'

    real_T Bod1[18];                   // Variable: Bod1
                                          //  Referenced by: '<S90>/MATLAB Function'

    real_T Bod2[18];                   // Variable: Bod2
                                          //  Referenced by: '<S160>/MATLAB Function'

    real_T Bod3[18];                   // Variable: Bod3
                                          //  Referenced by: '<S230>/MATLAB Function'

    real_T Cod[72];                    // Variable: Cod
                                          //  Referenced by: '<S12>/MATLAB Function'

    real_T Cod1[24];                   // Variable: Cod1
                                          //  Referenced by: '<S90>/MATLAB Function'

    real_T Cod2[30];                   // Variable: Cod2
                                          //  Referenced by: '<S160>/MATLAB Function'

    real_T Cod3[30];                   // Variable: Cod3
                                          //  Referenced by: '<S230>/MATLAB Function'

    real_T Dmn[36];                    // Variable: Dmn
                                          //  Referenced by: '<S12>/MATLAB Function'

    real_T Dmn1[9];                    // Variable: Dmn1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S160>/MATLAB Function'
                                          //    '<S230>/MATLAB Function'

    real_T Dod[36];                    // Variable: Dod
                                          //  Referenced by: '<S12>/MATLAB Function'

    real_T Dod1[12];                   // Variable: Dod1
                                          //  Referenced by: '<S90>/MATLAB Function'

    real_T Dod2[15];                   // Variable: Dod2
                                          //  Referenced by: '<S160>/MATLAB Function'

    real_T Dod3[15];                   // Variable: Dod3
                                          //  Referenced by: '<S230>/MATLAB Function'

    real_T beta;                       // Variable: beta
                                          //  Referenced by:
                                          //    '<Root>/SupervisoryController'
                                          //    '<S2>/Gain1'
                                          //    '<S3>/Gain1'
                                          //    '<S4>/Gain1'
                                          //    '<S5>/Gain1'
                                          //    '<S6>/Gain1'
                                          //    '<S9>/Gain'

    real_T dt;                         // Variable: dt
                                          //  Referenced by:
                                          //    '<Root>/SupervisoryController'
                                          //    '<S2>/MATLAB Function2'
                                          //    '<S3>/Gain2'
                                          //    '<S4>/Gain2'
                                          //    '<S5>/Gain2'
                                          //    '<S6>/Gain2'
                                          //    '<S9>/MATLAB Function'
                                          //    '<S322>/MATLAB Function1'
                                          //    '<S326>/MATLAB Function1'

    real_T lpfDen;                     // Variable: lpfDen
                                          //  Referenced by: '<S2>/Discrete Filter1'

    real_T lpfNum[60];                 // Variable: lpfNum
                                          //  Referenced by: '<S2>/Discrete Filter1'

    real_T mdlNum;                     // Variable: mdlNum
                                          //  Referenced by:
                                          //    '<S2>/MATLAB Function2'
                                          //    '<S322>/MATLAB Function1'
                                          //    '<S326>/MATLAB Function1'

    real_T uwt0[3];                    // Variable: uwt0
                                          //  Referenced by: '<S9>/Delay1'

    real_T ywt0[6];                    // Variable: ywt0
                                          //  Referenced by: '<S9>/Delay'

    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S86>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S10>/G_zero'

    real_T LastPcov_InitialCondition[324];// Expression: lastPcov
                                             //  Referenced by: '<S13>/LastPcov'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S10>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S10>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S13>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S13>/last_mv'

    real_T Constant12_Value[36];       // Expression: G.C
                                          //  Referenced by: '<S2>/Constant12'

    real_T Constant13_Value[18];       // Expression: G.D
                                          //  Referenced by: '<S2>/Constant13'

    real_T Constant2_Value[6];         // Expression: zeros(2*ns, 1)
                                          //  Referenced by: '<S2>/Constant2'

    real_T Constant1_Value[12];        // Expression: zeros(size(God.A,1),1)
                                          //  Referenced by: '<S12>/Constant1'

    real_T X0_Value[18];               // Expression: pInitialization.X0
                                          //  Referenced by: '<S43>/X0'

    real_T ym_zero_Value[6];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S13>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S10>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S10>/umin_zero'

    real_T Gain2_Gain;                 // Expression: -2
                                          //  Referenced by: '<S2>/Gain2'

    real_T Gain3_Gain;                 // Expression: 2
                                          //  Referenced by: '<S2>/Gain3'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S10>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S13>/umin_scale4'

    real_T F_zero_Value[6];            // Expression: zeros(1,6)
                                          //  Referenced by: '<S10>/F_zero'

    real_T ymin_scale1_Gain[6];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S13>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S10>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S13>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S10>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S10>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S13>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S10>/ecr.wt_zero'

    real_T P0_Value[324];              // Expression: pInitialization.P0
                                          //  Referenced by: '<S43>/P0'

    real_T Constant1_Value_c;          // Expression: 1
                                          //  Referenced by: '<S2>/Constant1'

    real_T H_Value[108];               // Expression: pInitialization.H
                                          //  Referenced by: '<S43>/H'

    real_T G_Value[324];               // Expression: pInitialization.G
                                          //  Referenced by: '<S43>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S13>/u_scale'

    real_T excitation_Mean[3];         // Expression: [0 0 0]
                                          //  Referenced by: '<S2>/excitation'

    real_T excitation_StdDev[3];       // Computed Parameter: excitation_StdDev
                                          //  Referenced by: '<S2>/excitation'

    real_T excitation_Seed[3];         // Expression: [12345 12346 12347]
                                          //  Referenced by: '<S2>/excitation'

    real_T DiscreteFilter1_InitialStates;// Expression: 0
                                            //  Referenced by: '<S2>/Discrete Filter1'

    real_T Saturation_UpperSat;        // Expression: 1000
                                          //  Referenced by: '<S2>/Saturation'

    real_T Saturation_LowerSat;        // Expression: 0
                                          //  Referenced by: '<S2>/Saturation'

    real_T Lykyhatkk1_Y0_c;            // Expression: 0
                                          //  Referenced by: '<S156>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0_b;                     // Computed Parameter: u_Y0_b
                                          //  Referenced by: '<S3>/u'

    real_T yhat_Y0_d;                  // Computed Parameter: yhat_Y0_d
                                          //  Referenced by: '<S3>/yhat'

    real_T Constant_Value;             // Expression: 0
                                          //  Referenced by: '<S3>/Constant'

    real_T DiscreteTimeIntegrator_gainval;
                           // Computed Parameter: DiscreteTimeIntegrator_gainval
                              //  Referenced by: '<S3>/Discrete-Time Integrator'

    real_T DiscreteTimeIntegrator_IC;  // Expression: 0
                                          //  Referenced by: '<S3>/Discrete-Time Integrator'

    real_T G_zero_Value_m;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S89>/G_zero'

    real_T ywt_zero_Value[4];          // Expression: zeros(4,1)
                                          //  Referenced by: '<S89>/y.wt_zero'

    real_T duwt_zero_Value_p[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S89>/du.wt_zero'

    real_T extmv_zero_Value_e[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S89>/ext.mv_zero'

    real_T extmv_scale_Gain_e[3];      // Expression: RMVscale
                                          //  Referenced by: '<S91>/ext.mv_scale'

    real_T mvtarget_zero_Value_k[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S89>/mv.target_zero'

    real_T extmv_scale1_Gain[3];       // Expression: RMVscale
                                          //  Referenced by: '<S91>/ext.mv_scale1'

    real_T last_mv_InitialCondition_f[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S91>/last_mv'

    real_T Constant4_Value[3];         // Expression: G0.B
                                          //  Referenced by: '<S3>/Constant4'

    real_T Constant12_Value_e[3];      // Expression: G0.C
                                          //  Referenced by: '<S3>/Constant12'

    real_T Constant13_Value_c[9];      // Expression: G0.D
                                          //  Referenced by: '<S3>/Constant13'

    real_T Constant3_Value;            // Expression: G0.A
                                          //  Referenced by: '<S3>/Constant3'

    real_T Constant2_Value_a;          // Expression: zeros(size(G0.A, 1), 1)
                                          //  Referenced by: '<S3>/Constant2'

    real_T Constant1_Value_j[6];       // Expression: zeros(size(God1.A,1),1)
                                          //  Referenced by: '<S90>/Constant1'

    real_T X0_Value_f[7];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S113>/X0'

    real_T ym_zero_Value_c[4];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S91>/ym_zero'

    real_T md_zero_Value_p;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S89>/md_zero'

    real_T umin_zero_Value_d[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S89>/umin_zero'

    real_T ymin_zero_Value[4];         // Expression: zeros(4,1)
                                          //  Referenced by: '<S89>/ymin_zero'

    real_T ymax_zero_Value[4];         // Expression: zeros(4,1)
                                          //  Referenced by: '<S89>/ymax_zero'

    real_T E_zero_Value_a[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S89>/E_zero'

    real_T umin_scale4_Gain_p[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S91>/umin_scale4'

    real_T F_zero_Value_g[4];          // Expression: zeros(1,4)
                                          //  Referenced by: '<S89>/F_zero'

    real_T ymin_scale1_Gain_j[4];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S91>/ymin_scale1'

    real_T S_zero_Value_g;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S89>/S_zero'

    real_T ymin_scale2_Gain_f;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S91>/ymin_scale2'

    real_T switch_zero_Value_e;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S89>/switch_zero'

    real_T ecrwt_zero_Value_o;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S89>/ecr.wt_zero'

    real_T P0_Value_a[49];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S113>/P0'

    real_T Constant1_Value_e;          // Expression: 1
                                          //  Referenced by: '<S3>/Constant1'

    real_T H_Value_o[21];              // Expression: pInitialization.H
                                          //  Referenced by: '<S113>/H'

    real_T G_Value_a[49];              // Expression: pInitialization.G
                                          //  Referenced by: '<S113>/G'

    real_T umin_scale1_Gain[3];        // Expression: MVscale
                                          //  Referenced by: '<S91>/umin_scale1'

    real_T Saturation_UpperSat_k;      // Expression: 1000
                                          //  Referenced by: '<S3>/Saturation'

    real_T Saturation_LowerSat_h;      // Expression: 0
                                          //  Referenced by: '<S3>/Saturation'

    real_T u_Y0_n;                     // Computed Parameter: u_Y0_n
                                          //  Referenced by: '<S4>/u'

    real_T yhat_Y0_k;                  // Computed Parameter: yhat_Y0_k
                                          //  Referenced by: '<S4>/yhat'

    real_T Constant_Value_p[2];        // Expression: [0;0]
                                          //  Referenced by: '<S4>/Constant'

    real_T DiscreteTimeIntegrator_gainva_b;
                          // Computed Parameter: DiscreteTimeIntegrator_gainva_b
                             //  Referenced by: '<S4>/Discrete-Time Integrator'

    real_T DiscreteTimeIntegrator_IC_n[2];// Expression: zeros(2, 1)
                                             //  Referenced by: '<S4>/Discrete-Time Integrator'

    real_T G_zero_Value_n;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S159>/G_zero'

    real_T ywt_zero_Value_n[5];        // Expression: zeros(5,1)
                                          //  Referenced by: '<S159>/y.wt_zero'

    real_T duwt_zero_Value_l[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S159>/du.wt_zero'

    real_T extmv_zero_Value_c[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S159>/ext.mv_zero'

    real_T extmv_scale_Gain_g[3];      // Expression: RMVscale
                                          //  Referenced by: '<S161>/ext.mv_scale'

    real_T mvtarget_zero_Value_e[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S159>/mv.target_zero'

    real_T extmv_scale1_Gain_b[3];     // Expression: RMVscale
                                          //  Referenced by: '<S161>/ext.mv_scale1'

    real_T last_mv_InitialCondition_b[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S161>/last_mv'

    real_T Constant3_Value_d[4];       // Expression: G1.A
                                          //  Referenced by: '<S4>/Constant3'

    real_T Constant4_Value_n[6];       // Expression: G1.B
                                          //  Referenced by: '<S4>/Constant4'

    real_T Constant12_Value_i[6];      // Expression: G1.C
                                          //  Referenced by: '<S4>/Constant12'

    real_T Constant13_Value_g[9];      // Expression: G1.D
                                          //  Referenced by: '<S4>/Constant13'

    real_T Constant2_Value_c[2];       // Expression: zeros(size(G1.A, 1), 1)
                                          //  Referenced by: '<S4>/Constant2'

    real_T Constant1_Value_p[6];       // Expression: zeros(size(God2.A,1),1)
                                          //  Referenced by: '<S160>/Constant1'

    real_T X0_Value_k[8];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S183>/X0'

    real_T ym_zero_Value_l[5];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S161>/ym_zero'

    real_T md_zero_Value_pu;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S159>/md_zero'

    real_T umin_zero_Value_e[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S159>/umin_zero'

    real_T ymin_zero_Value_g[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S159>/ymin_zero'

    real_T ymax_zero_Value_g[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S159>/ymax_zero'

    real_T E_zero_Value_b[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S159>/E_zero'

    real_T umin_scale4_Gain_g[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S161>/umin_scale4'

    real_T F_zero_Value_o[5];          // Expression: zeros(1,5)
                                          //  Referenced by: '<S159>/F_zero'

    real_T ymin_scale1_Gain_f[5];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S161>/ymin_scale1'

    real_T S_zero_Value_m;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S159>/S_zero'

    real_T ymin_scale2_Gain_g;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S161>/ymin_scale2'

    real_T switch_zero_Value_i;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S159>/switch_zero'

    real_T ecrwt_zero_Value_e;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S159>/ecr.wt_zero'

    real_T P0_Value_c[64];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S183>/P0'

    real_T Constant1_Value_pe;         // Expression: 1
                                          //  Referenced by: '<S4>/Constant1'

    real_T H_Value_k[24];              // Expression: pInitialization.H
                                          //  Referenced by: '<S183>/H'

    real_T G_Value_g[64];              // Expression: pInitialization.G
                                          //  Referenced by: '<S183>/G'

    real_T umin_scale1_Gain_p[3];      // Expression: MVscale
                                          //  Referenced by: '<S161>/umin_scale1'

    real_T Saturation_UpperSat_h;      // Expression: 1000
                                          //  Referenced by: '<S4>/Saturation'

    real_T Saturation_LowerSat_o;      // Expression: 0
                                          //  Referenced by: '<S4>/Saturation'

    real_T u_Y0_h;                     // Computed Parameter: u_Y0_h
                                          //  Referenced by: '<S5>/u'

    real_T yhat_Y0_f;                  // Computed Parameter: yhat_Y0_f
                                          //  Referenced by: '<S5>/yhat'

    real_T Constant_Value_e[2];        // Expression: [0;0]
                                          //  Referenced by: '<S5>/Constant'

    real_T DiscreteTimeIntegrator_gainva_k;
                          // Computed Parameter: DiscreteTimeIntegrator_gainva_k
                             //  Referenced by: '<S5>/Discrete-Time Integrator'

    real_T DiscreteTimeIntegrator_IC_c[2];// Expression: zeros(2, 1)
                                             //  Referenced by: '<S5>/Discrete-Time Integrator'

    real_T G_zero_Value_j;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S229>/G_zero'

    real_T ywt_zero_Value_g[5];        // Expression: zeros(5,1)
                                          //  Referenced by: '<S229>/y.wt_zero'

    real_T duwt_zero_Value_a[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S229>/du.wt_zero'

    real_T extmv_zero_Value_m[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S229>/ext.mv_zero'

    real_T extmv_scale_Gain_h[3];      // Expression: RMVscale
                                          //  Referenced by: '<S231>/ext.mv_scale'

    real_T mvtarget_zero_Value_d[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S229>/mv.target_zero'

    real_T extmv_scale1_Gain_e[3];     // Expression: RMVscale
                                          //  Referenced by: '<S231>/ext.mv_scale1'

    real_T last_mv_InitialCondition_i[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S231>/last_mv'

    real_T Constant3_Value_g[4];       // Expression: G2.A
                                          //  Referenced by: '<S5>/Constant3'

    real_T Constant4_Value_f[6];       // Expression: G2.B
                                          //  Referenced by: '<S5>/Constant4'

    real_T Constant12_Value_f[6];      // Expression: G2.C
                                          //  Referenced by: '<S5>/Constant12'

    real_T Constant13_Value_a[9];      // Expression: G2.D
                                          //  Referenced by: '<S5>/Constant13'

    real_T Constant2_Value_m[2];       // Expression: zeros(size(G2.A, 1), 1)
                                          //  Referenced by: '<S5>/Constant2'

    real_T Constant1_Value_h[6];       // Expression: zeros(size(God3.A,1),1)
                                          //  Referenced by: '<S230>/Constant1'

    real_T X0_Value_a[8];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S253>/X0'

    real_T ym_zero_Value_d[5];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S231>/ym_zero'

    real_T md_zero_Value_o;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S229>/md_zero'

    real_T umin_zero_Value_b[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S229>/umin_zero'

    real_T ymin_zero_Value_e[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S229>/ymin_zero'

    real_T ymax_zero_Value_d[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S229>/ymax_zero'

    real_T E_zero_Value_j[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S229>/E_zero'

    real_T umin_scale4_Gain_f[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S231>/umin_scale4'

    real_T F_zero_Value_n[5];          // Expression: zeros(1,5)
                                          //  Referenced by: '<S229>/F_zero'

    real_T ymin_scale1_Gain_e[5];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S231>/ymin_scale1'

    real_T S_zero_Value_i;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S229>/S_zero'

    real_T ymin_scale2_Gain_e;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S231>/ymin_scale2'

    real_T switch_zero_Value_k;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S229>/switch_zero'

    real_T ecrwt_zero_Value_j;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S229>/ecr.wt_zero'

    real_T P0_Value_m[64];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S253>/P0'

    real_T Constant1_Value_n;          // Expression: 1
                                          //  Referenced by: '<S5>/Constant1'

    real_T H_Value_oa[24];             // Expression: pInitialization.H
                                          //  Referenced by: '<S253>/H'

    real_T G_Value_h[64];              // Expression: pInitialization.G
                                          //  Referenced by: '<S253>/G'

    real_T umin_scale1_Gain_g[3];      // Expression: MVscale
                                          //  Referenced by: '<S231>/umin_scale1'

    real_T Saturation_UpperSat_c;      // Expression: 1000
                                          //  Referenced by: '<S5>/Saturation'

    real_T Saturation_LowerSat_b;      // Expression: 0
                                          //  Referenced by: '<S5>/Saturation'

    real_T ywt_Y0;                     // Computed Parameter: ywt_Y0
                                          //  Referenced by: '<S9>/ywt'

    real_T y_Y0;                       // Computed Parameter: y_Y0
                                          //  Referenced by: '<S9>/y_'

    real_T r_Y0;                       // Computed Parameter: r_Y0
                                          //  Referenced by: '<S9>/r_'

    real_T u_Y0_j;                     // Computed Parameter: u_Y0_j
                                          //  Referenced by: '<S6>/u'

    real_T Constant_Value_f[2];        // Expression: [0;0]
                                          //  Referenced by: '<S6>/Constant'

    real_T DiscreteTimeIntegrator_gainv_kx;
                          // Computed Parameter: DiscreteTimeIntegrator_gainv_kx
                             //  Referenced by: '<S6>/Discrete-Time Integrator'

    real_T DiscreteTimeIntegrator_IC_o[2];// Expression: zeros(2, 1)
                                             //  Referenced by: '<S6>/Discrete-Time Integrator'

    real_T G_zero_Value_a;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S299>/G_zero'

    real_T ywt_zero_Value_d[5];        // Expression: zeros(5,1)
                                          //  Referenced by: '<S299>/y.wt_zero'

    real_T duwt_zero_Value_lq[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S299>/du.wt_zero'

    real_T extmv_zero_Value_j[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S299>/ext.mv_zero'

    real_T extmv_scale_Gain_b[3];      // Expression: RMVscale
                                          //  Referenced by: '<S300>/ext.mv_scale'

    real_T mvtarget_zero_Value_n[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S299>/mv.target_zero'

    real_T extmv_scale1_Gain_j[3];     // Expression: RMVscale
                                          //  Referenced by: '<S300>/ext.mv_scale1'

    real_T last_mv_InitialCondition_g[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S300>/last_mv'

    real_T last_x_InitialCondition[11];// Expression: lastx+xoff
                                          //  Referenced by: '<S300>/last_x'

    real_T md_zero_Value_c;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S299>/md_zero'

    real_T umin_zero_Value_k[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S299>/umin_zero'

    real_T ymin_zero_Value_a[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S299>/ymin_zero'

    real_T ymax_zero_Value_k[5];       // Expression: zeros(5,1)
                                          //  Referenced by: '<S299>/ymax_zero'

    real_T E_zero_Value_l[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S299>/E_zero'

    real_T umin_scale4_Gain_m[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S300>/umin_scale4'

    real_T F_zero_Value_l[5];          // Expression: zeros(1,5)
                                          //  Referenced by: '<S299>/F_zero'

    real_T ymin_scale1_Gain_o[5];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S300>/ymin_scale1'

    real_T S_zero_Value_h;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S299>/S_zero'

    real_T ymin_scale2_Gain_n;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S300>/ymin_scale2'

    real_T switch_zero_Value_j;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S299>/switch_zero'

    real_T ecrwt_zero_Value_o0;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S299>/ecr.wt_zero'

    real_T umin_scale1_Gain_l[3];      // Expression: MVscale
                                          //  Referenced by: '<S300>/umin_scale1'

    real_T Saturation_UpperSat_cb;     // Expression: 1000
                                          //  Referenced by: '<S6>/Saturation'

    real_T Saturation_LowerSat_p;      // Expression: 0
                                          //  Referenced by: '<S6>/Saturation'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S41>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[246];// Expression: iA
                                              //  Referenced by: '<S13>/Memory'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S84>/isSqrtUsed'

    boolean_T Memory_InitialCondition_f[166];// Expression: iA
                                                //  Referenced by: '<S91>/Memory'

    boolean_T isSqrtUsed_Value_d;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S154>/isSqrtUsed'

    boolean_T Memory_InitialCondition_j[206];// Expression: iA
                                                //  Referenced by: '<S161>/Memory'

    boolean_T isSqrtUsed_Value_a;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S224>/isSqrtUsed'

    boolean_T Memory_InitialCondition_b[206];// Expression: iA
                                                //  Referenced by: '<S231>/Memory'

    boolean_T isSqrtUsed_Value_p;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S294>/isSqrtUsed'

    boolean_T Memory_InitialCondition_o[206];// Expression: iA
                                                //  Referenced by: '<S300>/Memory'

    P_MeasurementUpdate MeasurementUpdate_c;// '<S272>/MeasurementUpdate'
    P_MeasurementUpdate MeasurementUpdate_j;// '<S202>/MeasurementUpdate'
    P_paramEst1 paramEst2;             // '<S1>/paramEst2'
    P_paramEst1 paramEst1_o;           // '<S1>/paramEst1'
  };

  // Real-time Model Data Structure
  struct RT_MODEL {
    //
    //  DataMapInfo:
    //  The following substructure contains information regarding
    //  structures generated in the model's C API.

    struct {
      rtwCAPI_ModelMappingInfo mmi;
      void* dataAddress[238];
      int32_T* vardimsAddress[238];
      RTWLoggingFcnPtr loggingPtrs[238];
    } DataMapInfo;
  };

  // Copy Constructor
  SupervisoryController(SupervisoryController const&) = delete;

  // Assignment Operator
  SupervisoryController& operator= (SupervisoryController const&) & = delete;

  // Move Constructor
  SupervisoryController(SupervisoryController &&) = delete;

  // Move Assignment Operator
  SupervisoryController& operator= (SupervisoryController &&) = delete;

  // Real-Time Model get method
  SupervisoryController::RT_MODEL * getRTM();

  // External inputs
  ExtU rtU;

  // External outputs
  ExtY rtY;

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  SupervisoryController();

  // Destructor
  ~SupervisoryController();

  // private data and function members
 private:
  // Block states
  DW rtDW;

  // Tunable parameters
  static P rtP;

  // Triggered events
  PrevZCX rtPrevZCX;

  // private member function(s) for subsystem '<S1>/paramEst1'
  static void paramEst1_Init(real_T rty_theta[12], real_T rty_P[144], real_T
    rty_err[3], DW_paramEst1 *localDW, P_paramEst1 *localP);
  void paramEst1(const real_T rtu_y[3], const real_T rtu_y0[3], const real_T
                 rtu_u[3], const real_T rtu_u0[3], const boolean_T rtu_EN[3],
                 const real_T rtu_theta0[12], const real_T rtu_thetaSgn[12],
                 boolean_T rtu_rstTheta, const real_T rtu_P0[144], boolean_T
                 rtu_rstP, real_T rtu_p_, real_T rtu_lambda, real_T rty_theta[12],
                 real_T rty_P[144], real_T rty_err[3], DW_paramEst1 *localDW,
                 ZCE_paramEst1 *localZCE);
  void binary_expand_op(real_T in1[36], int32_T in2, int32_T in3, int32_T in4,
                        const real_T in5[12], int32_T in6, int32_T in7, const
                        real_T in8[12]);

  // private member function(s) for subsystem '<S113>/ScalarExpansionR'
  static void ScalarExpansionR(const real_T rtu_u[9], real_T rty_y[9]);

  // private member function(s) for subsystem '<S183>/CalculatePL'
  void CalculatePL(const real_T rtu_Ak[64], const real_T rtu_Ck[24], const
                   real_T rtu_Qbark[64], const real_T rtu_Rbark[9], const real_T
                   rtu_Nbark[24], boolean_T rtu_Enablek, const real_T rtu_Pk[64],
                   real_T rty_Mk[24], real_T rty_Lk[24], real_T rty_Zk[64],
                   real_T rty_Pk1[64]);
  void mrdiv(const real_T A[24], const real_T B_0[9], real_T Y[24]);

  // private member function(s) for subsystem '<S224>/SqrtUsedFcn'
  static void SqrtUsedFcn(const real_T rtu_u[64], boolean_T rtu_isSqrtUsed,
    real_T rty_P[64]);

  // private member function(s) for subsystem '<S202>/MeasurementUpdate'
  static void MeasurementUpdate_Init(real_T rty_Lykyhatkk1[8],
    P_MeasurementUpdate *localP);
  static void MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[8],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);
  void MeasurementUpdate(boolean_T rtu_Enable, const real_T rtu_Lk[24], const
    real_T rtu_yk[3], const real_T rtu_Ck[24], const real_T rtu_xhatkk1[8],
    const real_T rtu_Dk[9], const real_T rtu_uk[3], real_T rty_Lykyhatkk1[8],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);

  // private member function(s) for subsystem '<S183>/ReducedQRN'
  static void ReducedQRN(const real_T rtu_G[64], const real_T rtu_H[24], const
    real_T rtu_Q[64], const real_T rtu_R[9], const real_T rtu_N[24], real_T
    rty_Qbar[64], real_T rty_Rbar[9], real_T rty_Nbar[24]);

  // private member function(s) for subsystem '<S183>/ScalarExpansionQ'
  static void ScalarExpansionQ(const real_T rtu_u[64], real_T rty_y[64]);

  // private member function(s) for subsystem '<Root>'
  void do_vectors(const real_T b_data[], const int32_T *b_size, real_T c_data[],
                  int32_T c_size[2], int32_T ia_data[], int32_T *ia_size,
                  int32_T *ib_size);
  void ppval(const siswYcTR8LLamuD4YWmtXHC *pp, const real_T x_data[], const
             int32_T x_size[2], real_T v_data[], int32_T v_size[2]);
  void trajGen(const event_bus *event, const real_T y_[6], real_T trajectory
               [14400], uint16_T *trajectorySize);
  void handleEvent(real_T *holdT, boolean_T *eventDone, uint16_T *waypt_) const;
  boolean_T any(const real_T x[3]);
  real_T gainSchSig(const real_T ywt_[6]);
  int32_T xpotrf(real_T b_A[16]);
  real_T minimum(const real_T x[4]);
  void trisolve(const real_T b_A[16], real_T B_1[16]);
  real_T norm(const real_T x[4]);
  real_T maximum(const real_T x[4]);
  real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
             real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
             [4], real_T b_A[16], int32_T ia0);
  void KWIKfactor_o(const real_T b_Ac[824], const int32_T iC[206], int32_T nA,
                    const real_T b_Linv[16], real_T D[16], real_T b_H[16],
                    int32_T n, real_T RLinv[16], real_T *Status);
  void DropConstraint_f(int32_T kDrop, boolean_T iA[206], int32_T *nA, int32_T
                        iC[206]);
  void qpkwik_o(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T
                f[4], const real_T b_Ac[824], const real_T b[206], boolean_T iA
                [206], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
                lambda[206], int32_T *status);
  void KWIKfactor(const real_T b_Ac[664], const int32_T iC[166], int32_T nA,
                  const real_T b_Linv[16], real_T D[16], real_T b_H[16], int32_T
                  n, real_T RLinv[16], real_T *Status);
  void DropConstraint(int32_T kDrop, boolean_T iA[166], int32_T *nA, int32_T iC
                      [166]);
  void qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T f[4],
              const real_T b_Ac[664], const real_T b[166], boolean_T iA[166],
              int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[166],
              int32_T *status);
  void mrdiv_c(const real_T A[21], const real_T B_2[9], real_T Y[21]);

  // Real-Time Model
  RT_MODEL rtM;
};

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S13>/Floor' : Unused code path elimination
//  Block '<S13>/Floor1' : Unused code path elimination
//  Block '<S14>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S15>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S16>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S17>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S18>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S19>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S20>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S21>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S22>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S23>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S24>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S25>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S26>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S27>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S28>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S29>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S30>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S31>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S32>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S33>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S34>/Vector Dimension Check' : Unused code path elimination
//  Block '<S35>/Vector Dimension Check' : Unused code path elimination
//  Block '<S36>/Vector Dimension Check' : Unused code path elimination
//  Block '<S37>/Vector Dimension Check' : Unused code path elimination
//  Block '<S38>/Vector Dimension Check' : Unused code path elimination
//  Block '<S39>/Vector Dimension Check' : Unused code path elimination
//  Block '<S13>/last_x' : Unused code path elimination
//  Block '<S40>/Vector Dimension Check' : Unused code path elimination
//  Block '<S13>/useq_scale' : Unused code path elimination
//  Block '<S13>/useq_scale1' : Unused code path elimination
//  Block '<S10>/m_zero' : Unused code path elimination
//  Block '<S10>/p_zero' : Unused code path elimination
//  Block '<S12>/Display' : Unused code path elimination
//  Block '<S12>/Display1' : Unused code path elimination
//  Block '<S52>/Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Duplicate' : Unused code path elimination
//  Block '<S55>/Data Type Duplicate' : Unused code path elimination
//  Block '<S56>/Data Type Duplicate' : Unused code path elimination
//  Block '<S59>/Data Type Duplicate' : Unused code path elimination
//  Block '<S60>/Data Type Duplicate' : Unused code path elimination
//  Block '<S68>/CheckSignalProperties' : Unused code path elimination
//  Block '<S69>/CheckSignalProperties' : Unused code path elimination
//  Block '<S70>/CheckSignalProperties' : Unused code path elimination
//  Block '<S71>/CheckSignalProperties' : Unused code path elimination
//  Block '<S72>/CheckSignalProperties' : Unused code path elimination
//  Block '<S75>/CheckSignalProperties' : Unused code path elimination
//  Block '<S77>/CheckSignalProperties' : Unused code path elimination
//  Block '<S78>/CheckSignalProperties' : Unused code path elimination
//  Block '<S79>/CheckSignalProperties' : Unused code path elimination
//  Block '<S81>/CheckSignalProperties' : Unused code path elimination
//  Block '<S82>/CheckSignalProperties' : Unused code path elimination
//  Block '<S91>/Floor' : Unused code path elimination
//  Block '<S91>/Floor1' : Unused code path elimination
//  Block '<S92>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S93>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S94>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S95>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S96>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S97>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S98>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S99>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S100>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S101>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S102>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S103>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S104>/Vector Dimension Check' : Unused code path elimination
//  Block '<S105>/Vector Dimension Check' : Unused code path elimination
//  Block '<S106>/Vector Dimension Check' : Unused code path elimination
//  Block '<S107>/Vector Dimension Check' : Unused code path elimination
//  Block '<S108>/Vector Dimension Check' : Unused code path elimination
//  Block '<S109>/Vector Dimension Check' : Unused code path elimination
//  Block '<S91>/constant' : Unused code path elimination
//  Block '<S91>/last_x' : Unused code path elimination
//  Block '<S110>/Vector Dimension Check' : Unused code path elimination
//  Block '<S91>/umin_scale2' : Unused code path elimination
//  Block '<S91>/umin_scale3' : Unused code path elimination
//  Block '<S91>/umin_scale5' : Unused code path elimination
//  Block '<S89>/m_zero' : Unused code path elimination
//  Block '<S89>/p_zero' : Unused code path elimination
//  Block '<S122>/Data Type Duplicate' : Unused code path elimination
//  Block '<S123>/Data Type Duplicate' : Unused code path elimination
//  Block '<S125>/Data Type Duplicate' : Unused code path elimination
//  Block '<S126>/Data Type Duplicate' : Unused code path elimination
//  Block '<S129>/Data Type Duplicate' : Unused code path elimination
//  Block '<S130>/Data Type Duplicate' : Unused code path elimination
//  Block '<S138>/CheckSignalProperties' : Unused code path elimination
//  Block '<S139>/CheckSignalProperties' : Unused code path elimination
//  Block '<S140>/CheckSignalProperties' : Unused code path elimination
//  Block '<S141>/CheckSignalProperties' : Unused code path elimination
//  Block '<S142>/CheckSignalProperties' : Unused code path elimination
//  Block '<S145>/CheckSignalProperties' : Unused code path elimination
//  Block '<S147>/CheckSignalProperties' : Unused code path elimination
//  Block '<S148>/CheckSignalProperties' : Unused code path elimination
//  Block '<S149>/CheckSignalProperties' : Unused code path elimination
//  Block '<S151>/CheckSignalProperties' : Unused code path elimination
//  Block '<S152>/CheckSignalProperties' : Unused code path elimination
//  Block '<S161>/Floor' : Unused code path elimination
//  Block '<S161>/Floor1' : Unused code path elimination
//  Block '<S162>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S163>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S164>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S165>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S166>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S167>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S168>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S169>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S170>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S171>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S172>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S173>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S174>/Vector Dimension Check' : Unused code path elimination
//  Block '<S175>/Vector Dimension Check' : Unused code path elimination
//  Block '<S176>/Vector Dimension Check' : Unused code path elimination
//  Block '<S177>/Vector Dimension Check' : Unused code path elimination
//  Block '<S178>/Vector Dimension Check' : Unused code path elimination
//  Block '<S179>/Vector Dimension Check' : Unused code path elimination
//  Block '<S161>/constant' : Unused code path elimination
//  Block '<S161>/last_x' : Unused code path elimination
//  Block '<S180>/Vector Dimension Check' : Unused code path elimination
//  Block '<S161>/umin_scale2' : Unused code path elimination
//  Block '<S161>/umin_scale3' : Unused code path elimination
//  Block '<S161>/umin_scale5' : Unused code path elimination
//  Block '<S159>/m_zero' : Unused code path elimination
//  Block '<S159>/p_zero' : Unused code path elimination
//  Block '<S192>/Data Type Duplicate' : Unused code path elimination
//  Block '<S193>/Data Type Duplicate' : Unused code path elimination
//  Block '<S195>/Data Type Duplicate' : Unused code path elimination
//  Block '<S196>/Data Type Duplicate' : Unused code path elimination
//  Block '<S199>/Data Type Duplicate' : Unused code path elimination
//  Block '<S200>/Data Type Duplicate' : Unused code path elimination
//  Block '<S208>/CheckSignalProperties' : Unused code path elimination
//  Block '<S209>/CheckSignalProperties' : Unused code path elimination
//  Block '<S210>/CheckSignalProperties' : Unused code path elimination
//  Block '<S211>/CheckSignalProperties' : Unused code path elimination
//  Block '<S212>/CheckSignalProperties' : Unused code path elimination
//  Block '<S215>/CheckSignalProperties' : Unused code path elimination
//  Block '<S217>/CheckSignalProperties' : Unused code path elimination
//  Block '<S218>/CheckSignalProperties' : Unused code path elimination
//  Block '<S219>/CheckSignalProperties' : Unused code path elimination
//  Block '<S221>/CheckSignalProperties' : Unused code path elimination
//  Block '<S222>/CheckSignalProperties' : Unused code path elimination
//  Block '<S231>/Floor' : Unused code path elimination
//  Block '<S231>/Floor1' : Unused code path elimination
//  Block '<S232>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S233>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S234>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S235>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S236>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S237>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S238>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S239>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S240>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S241>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S242>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S243>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S244>/Vector Dimension Check' : Unused code path elimination
//  Block '<S245>/Vector Dimension Check' : Unused code path elimination
//  Block '<S246>/Vector Dimension Check' : Unused code path elimination
//  Block '<S247>/Vector Dimension Check' : Unused code path elimination
//  Block '<S248>/Vector Dimension Check' : Unused code path elimination
//  Block '<S249>/Vector Dimension Check' : Unused code path elimination
//  Block '<S231>/constant' : Unused code path elimination
//  Block '<S231>/last_x' : Unused code path elimination
//  Block '<S250>/Vector Dimension Check' : Unused code path elimination
//  Block '<S231>/umin_scale2' : Unused code path elimination
//  Block '<S231>/umin_scale3' : Unused code path elimination
//  Block '<S231>/umin_scale5' : Unused code path elimination
//  Block '<S229>/m_zero' : Unused code path elimination
//  Block '<S229>/p_zero' : Unused code path elimination
//  Block '<S262>/Data Type Duplicate' : Unused code path elimination
//  Block '<S263>/Data Type Duplicate' : Unused code path elimination
//  Block '<S265>/Data Type Duplicate' : Unused code path elimination
//  Block '<S266>/Data Type Duplicate' : Unused code path elimination
//  Block '<S269>/Data Type Duplicate' : Unused code path elimination
//  Block '<S270>/Data Type Duplicate' : Unused code path elimination
//  Block '<S278>/CheckSignalProperties' : Unused code path elimination
//  Block '<S279>/CheckSignalProperties' : Unused code path elimination
//  Block '<S280>/CheckSignalProperties' : Unused code path elimination
//  Block '<S281>/CheckSignalProperties' : Unused code path elimination
//  Block '<S282>/CheckSignalProperties' : Unused code path elimination
//  Block '<S285>/CheckSignalProperties' : Unused code path elimination
//  Block '<S287>/CheckSignalProperties' : Unused code path elimination
//  Block '<S288>/CheckSignalProperties' : Unused code path elimination
//  Block '<S289>/CheckSignalProperties' : Unused code path elimination
//  Block '<S291>/CheckSignalProperties' : Unused code path elimination
//  Block '<S292>/CheckSignalProperties' : Unused code path elimination
//  Block '<S300>/Floor' : Unused code path elimination
//  Block '<S300>/Floor1' : Unused code path elimination
//  Block '<S301>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S302>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S303>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S304>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S305>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S306>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S307>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S308>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S309>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S310>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S311>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S312>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S313>/Vector Dimension Check' : Unused code path elimination
//  Block '<S314>/Vector Dimension Check' : Unused code path elimination
//  Block '<S315>/Vector Dimension Check' : Unused code path elimination
//  Block '<S316>/Vector Dimension Check' : Unused code path elimination
//  Block '<S317>/Vector Dimension Check' : Unused code path elimination
//  Block '<S318>/Vector Dimension Check' : Unused code path elimination
//  Block '<S300>/constant' : Unused code path elimination
//  Block '<S319>/Vector Dimension Check' : Unused code path elimination
//  Block '<S300>/umin_scale2' : Unused code path elimination
//  Block '<S300>/umin_scale3' : Unused code path elimination
//  Block '<S300>/umin_scale5' : Unused code path elimination
//  Block '<S300>/ym_zero' : Unused code path elimination
//  Block '<S299>/m_zero' : Unused code path elimination
//  Block '<S299>/p_zero' : Unused code path elimination
//  Block '<S322>/Display' : Unused code path elimination
//  Block '<S322>/Display1' : Unused code path elimination
//  Block '<S322>/Display2' : Unused code path elimination
//  Block '<S324>/Display' : Unused code path elimination
//  Block '<S7>/Scope' : Unused code path elimination
//  Block '<S326>/Display' : Unused code path elimination
//  Block '<S326>/Display1' : Unused code path elimination
//  Block '<S326>/Display2' : Unused code path elimination
//  Block '<S328>/Display' : Unused code path elimination
//  Block '<S8>/Scope' : Unused code path elimination
//  Block '<S13>/Reshape' : Reshape block reduction
//  Block '<S13>/Reshape1' : Reshape block reduction
//  Block '<S13>/Reshape2' : Reshape block reduction
//  Block '<S13>/Reshape3' : Reshape block reduction
//  Block '<S13>/Reshape4' : Reshape block reduction
//  Block '<S13>/Reshape5' : Reshape block reduction
//  Block '<S55>/Conversion' : Eliminate redundant data type conversion
//  Block '<S59>/Conversion' : Eliminate redundant data type conversion
//  Block '<S62>/Reshape' : Reshape block reduction
//  Block '<S43>/ReshapeX0' : Reshape block reduction
//  Block '<S43>/Reshapeu' : Reshape block reduction
//  Block '<S43>/Reshapexhat' : Reshape block reduction
//  Block '<S43>/Reshapey' : Reshape block reduction
//  Block '<S43>/Reshapeyhat' : Reshape block reduction
//  Block '<S91>/Reshape' : Reshape block reduction
//  Block '<S91>/Reshape1' : Reshape block reduction
//  Block '<S91>/Reshape2' : Reshape block reduction
//  Block '<S91>/Reshape3' : Reshape block reduction
//  Block '<S91>/Reshape4' : Reshape block reduction
//  Block '<S91>/Reshape5' : Reshape block reduction
//  Block '<S125>/Conversion' : Eliminate redundant data type conversion
//  Block '<S129>/Conversion' : Eliminate redundant data type conversion
//  Block '<S132>/Reshape' : Reshape block reduction
//  Block '<S113>/ReshapeX0' : Reshape block reduction
//  Block '<S113>/Reshapeu' : Reshape block reduction
//  Block '<S113>/Reshapexhat' : Reshape block reduction
//  Block '<S113>/Reshapey' : Reshape block reduction
//  Block '<S113>/Reshapeyhat' : Reshape block reduction
//  Block '<S161>/Reshape' : Reshape block reduction
//  Block '<S161>/Reshape1' : Reshape block reduction
//  Block '<S161>/Reshape2' : Reshape block reduction
//  Block '<S161>/Reshape3' : Reshape block reduction
//  Block '<S161>/Reshape4' : Reshape block reduction
//  Block '<S161>/Reshape5' : Reshape block reduction
//  Block '<S195>/Conversion' : Eliminate redundant data type conversion
//  Block '<S199>/Conversion' : Eliminate redundant data type conversion
//  Block '<S202>/Reshape' : Reshape block reduction
//  Block '<S183>/ReshapeX0' : Reshape block reduction
//  Block '<S183>/Reshapeu' : Reshape block reduction
//  Block '<S183>/Reshapexhat' : Reshape block reduction
//  Block '<S183>/Reshapey' : Reshape block reduction
//  Block '<S183>/Reshapeyhat' : Reshape block reduction
//  Block '<S231>/Reshape' : Reshape block reduction
//  Block '<S231>/Reshape1' : Reshape block reduction
//  Block '<S231>/Reshape2' : Reshape block reduction
//  Block '<S231>/Reshape3' : Reshape block reduction
//  Block '<S231>/Reshape4' : Reshape block reduction
//  Block '<S231>/Reshape5' : Reshape block reduction
//  Block '<S265>/Conversion' : Eliminate redundant data type conversion
//  Block '<S269>/Conversion' : Eliminate redundant data type conversion
//  Block '<S272>/Reshape' : Reshape block reduction
//  Block '<S253>/ReshapeX0' : Reshape block reduction
//  Block '<S253>/Reshapeu' : Reshape block reduction
//  Block '<S253>/Reshapexhat' : Reshape block reduction
//  Block '<S253>/Reshapey' : Reshape block reduction
//  Block '<S253>/Reshapeyhat' : Reshape block reduction
//  Block '<S300>/Reshape' : Reshape block reduction
//  Block '<S300>/Reshape1' : Reshape block reduction
//  Block '<S300>/Reshape2' : Reshape block reduction
//  Block '<S300>/Reshape3' : Reshape block reduction
//  Block '<S300>/Reshape4' : Reshape block reduction
//  Block '<S300>/Reshape5' : Reshape block reduction


//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Note that this particular code originates from a subsystem build,
//  and has its own system numbers different from the parent model.
//  Refer to the system hierarchy for this subsystem below, and use the
//  MATLAB hilite_system command to trace the generated code back
//  to the parent model.  For example,
//
//  hilite_system('T_junction_mpc/SupervisoryController')    - opens subsystem T_junction_mpc/SupervisoryController
//  hilite_system('T_junction_mpc/SupervisoryController/Kp') - opens and selects block Kp
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'T_junction_mpc'
//  '<S1>'   : 'T_junction_mpc/SupervisoryController'
//  '<S2>'   : 'T_junction_mpc/SupervisoryController/ampc'
//  '<S3>'   : 'T_junction_mpc/SupervisoryController/mpc1'
//  '<S4>'   : 'T_junction_mpc/SupervisoryController/mpc2'
//  '<S5>'   : 'T_junction_mpc/SupervisoryController/mpc3'
//  '<S6>'   : 'T_junction_mpc/SupervisoryController/mpcg'
//  '<S7>'   : 'T_junction_mpc/SupervisoryController/paramEst1'
//  '<S8>'   : 'T_junction_mpc/SupervisoryController/paramEst2'
//  '<S9>'   : 'T_junction_mpc/SupervisoryController/wtMod'
//  '<S10>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller'
//  '<S11>'  : 'T_junction_mpc/SupervisoryController/ampc/MATLAB Function2'
//  '<S12>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)'
//  '<S13>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC'
//  '<S14>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S15>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S16>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S17>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S18>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S19>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S20>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S21>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S22>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S23>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S24>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S25>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S26>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S27>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S28>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S29>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S30>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S31>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S32>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S33>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S34>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S35>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S36>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S37>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S38>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S39>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S40>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/moorx'
//  '<S41>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer'
//  '<S42>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S43>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2'
//  '<S44>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/MATLAB Function'
//  '<S45>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculatePL'
//  '<S46>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculateYhat'
//  '<S47>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S48>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S49>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S50>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S51>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S52>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S53>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S54>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S55>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S56>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S57>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S58>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S59>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S60>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S61>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S62>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer'
//  '<S63>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ReducedQRN'
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkA'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkB'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkC'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkD'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkEnable'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkG'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkH'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkN'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkP0'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkQ'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkR'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkReset'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkX0'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checku'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checky'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Matrix Signal Check'
//  '<S93>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Matrix Signal Check1'
//  '<S94>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Matrix Signal Check2'
//  '<S95>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check'
//  '<S96>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check1'
//  '<S97>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check2'
//  '<S98>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check3'
//  '<S99>'  : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check4'
//  '<S100>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check5'
//  '<S101>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check6'
//  '<S102>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check7'
//  '<S103>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Preview Signal Check8'
//  '<S104>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Scalar Signal Check'
//  '<S105>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Scalar Signal Check1'
//  '<S106>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Scalar Signal Check2'
//  '<S107>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Vector Signal Check'
//  '<S108>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Vector Signal Check1'
//  '<S109>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/MPC Vector Signal Check6'
//  '<S110>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/moorx'
//  '<S111>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/optimizer'
//  '<S112>' : 'T_junction_mpc/SupervisoryController/mpc1/MPC Controller1/MPC/optimizer/optimizer'
//  '<S113>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2'
//  '<S114>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/MATLAB Function'
//  '<S115>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CalculatePL'
//  '<S116>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CalculateYhat'
//  '<S117>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S118>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S119>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S120>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S121>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S122>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S123>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S124>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S125>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S126>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S127>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S128>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S129>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S130>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S131>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S132>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/Observer'
//  '<S133>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ReducedQRN'
//  '<S134>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S135>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S136>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S137>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S138>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkA'
//  '<S139>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkB'
//  '<S140>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkC'
//  '<S141>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkD'
//  '<S142>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkEnable'
//  '<S143>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkG'
//  '<S144>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkH'
//  '<S145>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkN'
//  '<S146>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkP0'
//  '<S147>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkQ'
//  '<S148>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkR'
//  '<S149>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkReset'
//  '<S150>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checkX0'
//  '<S151>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checku'
//  '<S152>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/checky'
//  '<S153>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S154>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S155>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S156>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S157>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S158>' : 'T_junction_mpc/SupervisoryController/mpc1/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S159>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1'
//  '<S160>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)'
//  '<S161>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC'
//  '<S162>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Matrix Signal Check'
//  '<S163>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Matrix Signal Check1'
//  '<S164>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Matrix Signal Check2'
//  '<S165>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check'
//  '<S166>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check1'
//  '<S167>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check2'
//  '<S168>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check3'
//  '<S169>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check4'
//  '<S170>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check5'
//  '<S171>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check6'
//  '<S172>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check7'
//  '<S173>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Preview Signal Check8'
//  '<S174>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Scalar Signal Check'
//  '<S175>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Scalar Signal Check1'
//  '<S176>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Scalar Signal Check2'
//  '<S177>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Vector Signal Check'
//  '<S178>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Vector Signal Check1'
//  '<S179>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/MPC Vector Signal Check6'
//  '<S180>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/moorx'
//  '<S181>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/optimizer'
//  '<S182>' : 'T_junction_mpc/SupervisoryController/mpc2/MPC Controller1/MPC/optimizer/optimizer'
//  '<S183>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2'
//  '<S184>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/MATLAB Function'
//  '<S185>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CalculatePL'
//  '<S186>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CalculateYhat'
//  '<S187>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S188>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S189>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S190>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S191>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S192>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S193>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S194>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S195>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S196>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S197>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S198>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S199>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S200>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S201>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S202>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/Observer'
//  '<S203>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ReducedQRN'
//  '<S204>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S205>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S206>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S207>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S208>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkA'
//  '<S209>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkB'
//  '<S210>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkC'
//  '<S211>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkD'
//  '<S212>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkEnable'
//  '<S213>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkG'
//  '<S214>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkH'
//  '<S215>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkN'
//  '<S216>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkP0'
//  '<S217>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkQ'
//  '<S218>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkR'
//  '<S219>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkReset'
//  '<S220>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checkX0'
//  '<S221>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checku'
//  '<S222>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/checky'
//  '<S223>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S224>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S225>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S226>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S227>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S228>' : 'T_junction_mpc/SupervisoryController/mpc2/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S229>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1'
//  '<S230>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)'
//  '<S231>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC'
//  '<S232>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Matrix Signal Check'
//  '<S233>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Matrix Signal Check1'
//  '<S234>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Matrix Signal Check2'
//  '<S235>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check'
//  '<S236>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check1'
//  '<S237>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check2'
//  '<S238>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check3'
//  '<S239>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check4'
//  '<S240>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check5'
//  '<S241>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check6'
//  '<S242>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check7'
//  '<S243>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Preview Signal Check8'
//  '<S244>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Scalar Signal Check'
//  '<S245>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Scalar Signal Check1'
//  '<S246>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Scalar Signal Check2'
//  '<S247>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Vector Signal Check'
//  '<S248>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Vector Signal Check1'
//  '<S249>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/MPC Vector Signal Check6'
//  '<S250>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/moorx'
//  '<S251>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/optimizer'
//  '<S252>' : 'T_junction_mpc/SupervisoryController/mpc3/MPC Controller1/MPC/optimizer/optimizer'
//  '<S253>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2'
//  '<S254>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/MATLAB Function'
//  '<S255>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CalculatePL'
//  '<S256>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CalculateYhat'
//  '<S257>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S258>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S259>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S260>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S261>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S262>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S263>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S264>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S265>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S266>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S267>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S268>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S269>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S270>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S271>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S272>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/Observer'
//  '<S273>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ReducedQRN'
//  '<S274>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S275>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S276>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S277>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S278>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkA'
//  '<S279>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkB'
//  '<S280>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkC'
//  '<S281>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkD'
//  '<S282>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkEnable'
//  '<S283>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkG'
//  '<S284>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkH'
//  '<S285>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkN'
//  '<S286>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkP0'
//  '<S287>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkQ'
//  '<S288>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkR'
//  '<S289>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkReset'
//  '<S290>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checkX0'
//  '<S291>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checku'
//  '<S292>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/checky'
//  '<S293>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S294>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S295>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S296>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S297>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S298>' : 'T_junction_mpc/SupervisoryController/mpc3/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S299>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1'
//  '<S300>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC'
//  '<S301>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Matrix Signal Check'
//  '<S302>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Matrix Signal Check1'
//  '<S303>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Matrix Signal Check2'
//  '<S304>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check'
//  '<S305>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check1'
//  '<S306>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check2'
//  '<S307>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check3'
//  '<S308>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check4'
//  '<S309>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check5'
//  '<S310>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check6'
//  '<S311>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check7'
//  '<S312>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Preview Signal Check8'
//  '<S313>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Scalar Signal Check'
//  '<S314>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Scalar Signal Check1'
//  '<S315>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Scalar Signal Check2'
//  '<S316>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Vector Signal Check'
//  '<S317>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Vector Signal Check1'
//  '<S318>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/MPC Vector Signal Check6'
//  '<S319>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/moorx'
//  '<S320>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/optimizer'
//  '<S321>' : 'T_junction_mpc/SupervisoryController/mpcg/MPC Controller1/MPC/optimizer/optimizer'
//  '<S322>' : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)'
//  '<S323>' : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1'
//  '<S324>' : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/RLS'
//  '<S325>' : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function'
//  '<S326>' : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)'
//  '<S327>' : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/MATLAB Function1'
//  '<S328>' : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/RLS'
//  '<S329>' : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/RLS/MATLAB Function'
//  '<S330>' : 'T_junction_mpc/SupervisoryController/wtMod/MATLAB Function'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
