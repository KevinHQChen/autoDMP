//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.959
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May 15 05:05:46 2023
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
#include "rtwtypes.h"
#include <stddef.h>
#include "zero_crossing_types.h"
#ifndef DEFINED_TYPEDEF_FOR_event_bus_
#define DEFINED_TYPEDEF_FOR_event_bus_

struct event_bus
{
  real_T srcState;
  real_T destState;
  real_T destPos[3];
  real_T moveTime;
  real_T holdTime;
  boolean_T chs[3];
  boolean_T nextChs[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_mdl0_bus_
#define DEFINED_TYPEDEF_FOR_mdl0_bus_

struct mdl0_bus
{
  real_T A;
  real_T B[3];
  real_T C;
  real_T D[3];
  real_T U[3];
  real_T Y;
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
  real_T C[4];
  real_T D[6];
  real_T U[3];
  real_T Y[2];
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

#ifndef DEFINED_TYPEDEF_FOR_struct_k8WKy8tDBVlN8BL9RXVTiF_
#define DEFINED_TYPEDEF_FOR_struct_k8WKy8tDBVlN8BL9RXVTiF_

struct struct_k8WKy8tDBVlN8BL9RXVTiF
{
  uint8_T estimationMethod;
  int32_T nParameters;
  boolean_T isUsingFrames;
  int32_T windowLength;
};

#endif

// Custom Type definition for MATLAB Function: '<S38>/RLS'
#ifndef struct_c_controllib_internal_blocks_rl
#define struct_c_controllib_internal_blocks_rl

struct c_controllib_internal_blocks_rl
{
  int32_T IteratorPosition;
};

#endif                                // struct_c_controllib_internal_blocks_rl

#ifndef struct_d_controllib_internal_blocks_rl
#define struct_d_controllib_internal_blocks_rl

struct d_controllib_internal_blocks_rl
{
  c_controllib_internal_blocks_rl DataIterator;
};

#endif                                // struct_d_controllib_internal_blocks_rl

// Custom Type definition for Chart: '<Root>/SupervisoryController'
#ifndef struct_cell_wrap_5
#define struct_cell_wrap_5

struct cell_wrap_5
{
  real_T f1[4];
};

#endif                                 // struct_cell_wrap_5

#ifndef struct_emxArray_real_T_9x3
#define struct_emxArray_real_T_9x3

struct emxArray_real_T_9x3
{
  real_T data[27];
  int32_T size[2];
};

#endif                                 // struct_emxArray_real_T_9x3

#ifndef struct_cell_wrap_6
#define struct_cell_wrap_6

struct cell_wrap_6
{
  emxArray_real_T_9x3 f1;
};

#endif                                 // struct_cell_wrap_6

#ifndef struct_emxArray_real_T_3x5x3
#define struct_emxArray_real_T_3x5x3

struct emxArray_real_T_3x5x3
{
  real_T data[45];
  int32_T size[3];
};

#endif                                 // struct_emxArray_real_T_3x5x3

// Custom Type definition for Chart: '<Root>/SupervisoryController'
#ifndef struct_s_vjEZ2dxatR8VOmLd9oOqoD
#define struct_s_vjEZ2dxatR8VOmLd9oOqoD

struct s_vjEZ2dxatR8VOmLd9oOqoD
{
  real_T breaks[6];
  emxArray_real_T_3x5x3 coefs;
};

#endif                                 // struct_s_vjEZ2dxatR8VOmLd9oOqoD

#ifndef struct_emxArray_cell_wrap_6_3
#define struct_emxArray_cell_wrap_6_3

struct emxArray_cell_wrap_6_3
{
  cell_wrap_6 data[3];
  int32_T size;
};

#endif                                 // struct_emxArray_cell_wrap_6_3

#ifndef struct_emxArray_s_vjEZ2dxatR8VOmLd9oOq
#define struct_emxArray_s_vjEZ2dxatR8VOmLd9oOq

struct emxArray_s_vjEZ2dxatR8VOmLd9oOq
{
  s_vjEZ2dxatR8VOmLd9oOqoD data[3];
  int32_T size;
};

#endif                                // struct_emxArray_s_vjEZ2dxatR8VOmLd9oOq

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
  // Block signals and states (default storage) for system '<S1>/State0.controlLaw.AMPC0' 
  struct DW_State0controlLawAMPC0 {
    d_controllib_internal_blocks_rl rlsEstimator;// '<S38>/RLS'
    real_T Product3[2];                // '<S108>/Product3'
    real_T last_mv_DSTATE[3];          // '<S8>/last_mv'
    real_T UnitDelay2_DSTATE[3];       // '<S6>/Unit Delay2'
    real_T delayTheta_DSTATE[2];       // '<S38>/delayTheta'
    real_T delayL_DSTATE[4];           // '<S38>/delayL'
    real_T MemoryX_DSTATE[2];          // '<S65>/MemoryX'
    real_T MemoryP_DSTATE[4];          // '<S65>/MemoryP'
    real_T NextOutput[3];              // '<S2>/Measurement Noise'
    real_T dv[966];
    real_T Su[1200];
    real_T UnitDelay3_DSTATE;          // '<S6>/Unit Delay3'
    uint32_T RandSeed[3];              // '<S2>/Measurement Noise'
    boolean_T Memory_PreviousInput[46];// '<S8>/Memory'
    boolean_T icLoad;                  // '<S38>/delayBuffery'
    boolean_T icLoad_n;                // '<S38>/delayBufferH'
    boolean_T icLoad_m;                // '<S38>/delayTheta'
    boolean_T icLoad_e;                // '<S38>/delayL'
    boolean_T icLoad_k;                // '<S65>/MemoryX'
    boolean_T icLoad_f;                // '<S65>/MemoryP'
    boolean_T rlsEstimator_not_empty;  // '<S38>/RLS'
    boolean_T MeasurementUpdate_MODE;  // '<S84>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state for system '<S1>/State0.controlLaw.AMPC0'
  struct ZCE_State0controlLawAMPC0 {
    ZCSigState MemoryX_Reset_ZCE_j;    // '<S65>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_a;    // '<S65>/MemoryP'
  };

  // Block signals and states (default storage) for system '<S145>/RLS'
  struct DW_RLS {
    d_controllib_internal_blocks_rl rlsEstimator;// '<S145>/RLS'
    boolean_T rlsEstimator_not_empty;  // '<S145>/RLS'
  };

  // Block signals and states (default storage) for system '<S218>/MeasurementUpdate' 
  struct DW_MeasurementUpdate {
    boolean_T MeasurementUpdate_MODE;  // '<S218>/MeasurementUpdate'
  };

  // Block signals and states (default storage) for system '<S1>/State1.controlLaw.AMPC1' 
  struct DW_State1controlLawAMPC1 {
    DW_MeasurementUpdate MeasurementUpdate_n;// '<S218>/MeasurementUpdate'
    DW_RLS sf_RLS_g;                   // '<S146>/RLS'
    DW_RLS sf_RLS;                     // '<S145>/RLS'
    real_T Product3[4];                // '<S242>/Product3'
    real_T last_mv_DSTATE[3];          // '<S114>/last_mv'
    real_T UnitDelay3_DSTATE[2];       // '<S112>/Unit Delay3'
    real_T UnitDelay2_DSTATE[3];       // '<S112>/Unit Delay2'
    real_T delayTheta_DSTATE[5];       // '<S145>/delayTheta'
    real_T delayL_DSTATE[25];          // '<S145>/delayL'
    real_T delayTheta_DSTATE_i[5];     // '<S146>/delayTheta'
    real_T delayL_DSTATE_o[25];        // '<S146>/delayL'
    real_T MemoryX_DSTATE[4];          // '<S199>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S199>/MemoryP'
    real_T NextOutput[3];              // '<S3>/Measurement Noise'
    real_T dv[1806];
    real_T b_Hv[840];
    real_T Su[2400];
    real_T delayBuffery_DSTATE;        // '<S145>/delayBuffery'
    real_T delayBufferH_DSTATE;        // '<S145>/delayBufferH'
    real_T delayBuffery_DSTATE_i;      // '<S146>/delayBuffery'
    real_T delayBufferH_DSTATE_h;      // '<S146>/delayBufferH'
    uint32_T RandSeed[3];              // '<S3>/Measurement Noise'
    boolean_T Memory_PreviousInput[86];// '<S114>/Memory'
    boolean_T icLoad;                  // '<S145>/delayBuffery'
    boolean_T icLoad_p;                // '<S145>/delayBufferH'
    boolean_T icLoad_b;                // '<S145>/delayTheta'
    boolean_T icLoad_d;                // '<S145>/delayL'
    boolean_T icLoad_bx;               // '<S146>/delayBuffery'
    boolean_T icLoad_dn;               // '<S146>/delayBufferH'
    boolean_T icLoad_c;                // '<S146>/delayTheta'
    boolean_T icLoad_e;                // '<S146>/delayL'
    boolean_T icLoad_f;                // '<S199>/MemoryX'
    boolean_T icLoad_a;                // '<S199>/MemoryP'
  };

  // Zero-crossing (trigger) state for system '<S1>/State1.controlLaw.AMPC1'
  struct ZCE_State1controlLawAMPC1 {
    ZCSigState MemoryX_Reset_ZCE_b;    // '<S199>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_k;    // '<S199>/MemoryP'
  };

  // Block signals and states (default storage) for system '<S279>/RLS'
  struct DW_RLS_k {
    d_controllib_internal_blocks_rl rlsEstimator;// '<S279>/RLS'
    boolean_T rlsEstimator_not_empty;  // '<S279>/RLS'
  };

  // Block signals and states (default storage) for system '<S1>/State2.controlLaw.AMPC2' 
  struct DW_State2controlLawAMPC2 {
    DW_MeasurementUpdate MeasurementUpdate_a;// '<S352>/MeasurementUpdate'
    DW_RLS_k sf_RLS_i;                 // '<S280>/RLS'
    DW_RLS_k sf_RLS;                   // '<S279>/RLS'
    real_T Product3[4];                // '<S376>/Product3'
    real_T last_mv_DSTATE[3];          // '<S248>/last_mv'
    real_T UnitDelay1_DSTATE[3];       // '<S246>/Unit Delay1'
    real_T UnitDelay7_DSTATE[2];       // '<S246>/Unit Delay7'
    real_T delayTheta_DSTATE[3];       // '<S279>/delayTheta'
    real_T delayL_DSTATE[9];           // '<S279>/delayL'
    real_T delayTheta_DSTATE_m[3];     // '<S280>/delayTheta'
    real_T delayL_DSTATE_c[9];         // '<S280>/delayL'
    real_T MemoryX_DSTATE[4];          // '<S333>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S333>/MemoryP'
    real_T NextOutput[3];              // '<S4>/Measurement Noise'
    real_T dv[1806];
    real_T b_Hv[840];
    real_T Su[2400];
    real_T delayBuffery_DSTATE;        // '<S279>/delayBuffery'
    real_T delayBufferH_DSTATE;        // '<S279>/delayBufferH'
    real_T delayBuffery_DSTATE_l;      // '<S280>/delayBuffery'
    real_T delayBufferH_DSTATE_i;      // '<S280>/delayBufferH'
    uint32_T RandSeed[3];              // '<S4>/Measurement Noise'
    boolean_T Memory_PreviousInput[86];// '<S248>/Memory'
    boolean_T icLoad;                  // '<S279>/delayBuffery'
    boolean_T icLoad_a;                // '<S279>/delayBufferH'
    boolean_T icLoad_g;                // '<S279>/delayTheta'
    boolean_T icLoad_c;                // '<S279>/delayL'
    boolean_T icLoad_o;                // '<S280>/delayBuffery'
    boolean_T icLoad_n;                // '<S280>/delayBufferH'
    boolean_T icLoad_cf;               // '<S280>/delayTheta'
    boolean_T icLoad_l;                // '<S280>/delayL'
    boolean_T icLoad_l1;               // '<S333>/MemoryX'
    boolean_T icLoad_h;                // '<S333>/MemoryP'
  };

  // Zero-crossing (trigger) state for system '<S1>/State2.controlLaw.AMPC2'
  struct ZCE_State2controlLawAMPC2 {
    ZCSigState MemoryX_Reset_ZCE;      // '<S333>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S333>/MemoryP'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    DW_State2controlLawAMPC2 State2controlLawAMPC2_l;// '<S1>/State2.controlLaw.AMPC2' 
    DW_State1controlLawAMPC1 State1controlLawAMPC1_g;// '<S1>/State1.controlLaw.AMPC1' 
    DW_State0controlLawAMPC0 State0controlLawAMPC0_n;// '<S1>/State0.controlLaw.AMPC0' 
    real_T traj[7200];                 // '<Root>/SupervisoryController'
    real_T yhat1[2];                   // '<Root>/SupervisoryController'
    real_T yhat2[2];                   // '<Root>/SupervisoryController'
    real_T b_data[3597];
    real_T c_data[3597];
    real_T d_data[3597];
    real_T t_data[1199];
    real_T tmp_data[3597];
    real_T holdT;                      // '<Root>/SupervisoryController'
    real_T yhat0;                      // '<Root>/SupervisoryController'
    uint16_T chs1[2];                  // '<Root>/SupervisoryController'
    uint16_T chs2[2];                  // '<Root>/SupervisoryController'
    uint16_T waypt;                    // '<Root>/SupervisoryController'
    uint16_T trajSize;                 // '<Root>/SupervisoryController'
    uint16_T chs0;                     // '<Root>/SupervisoryController'
    uint8_T is_c6_SupervisoryController;// '<Root>/SupervisoryController'
    uint8_T is_EventHandler;           // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_n;         // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_k;         // '<Root>/SupervisoryController'
    uint8_T is_active_c6_SupervisoryControl;// '<Root>/SupervisoryController'
    boolean_T evDone;                  // '<Root>/SupervisoryController'
  };

  // Zero-crossing (trigger) state
  struct PrevZCX {
    ZCSigState SupervisoryController_Trig_ZCE[2];// '<Root>/SupervisoryController' 
    ZCE_State2controlLawAMPC2 State2controlLawAMPC2_l;// '<S1>/State2.controlLaw.AMPC2' 
    ZCE_State1controlLawAMPC1 State1controlLawAMPC1_g;// '<S1>/State1.controlLaw.AMPC1' 
    ZCE_State0controlLawAMPC0 State0controlLawAMPC0_n;// '<S1>/State0.controlLaw.AMPC0' 
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T y[3];                       // '<Root>/y'
    real_T y_max[3];                   // '<Root>/y_max'
    real_T y_o[3];                     // '<Root>/y_o'
    real_T u_o[3];                     // '<Root>/u_o'
    event_bus nextEv;                  // '<Root>/nextEv'
    real_T y_range[2];                 // '<Root>/y_range'
    boolean_T enAdapt[3];              // '<Root>/enAdapt'
    real_T excitation;                 // '<Root>/excitation'
    real_T inputevents[2];             // '<Root>/input events'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T u[3];                       // '<Root>/u'
    boolean_T inTransRegion;           // '<Root>/inTransRegion'
    event_bus currEv;                  // '<Root>/currEv'
    boolean_T requestEvent;            // '<Root>/requestEvent'
    real_T currTraj[3];                // '<Root>/currTraj'
    real_T yhat[3];                    // '<Root>/yhat'
    real_T B_o[9];                     // '<Root>/B'
  };

  // Parameters for system: '<S1>/State0.controlLaw.AMPC0'
  struct P_State0controlLawAMPC0 {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S108>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

    real_T params_Y0;                  // Computed Parameter: params_Y0
                                          //  Referenced by: '<S2>/params'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/G_zero'

    real_T LastPcov_InitialCondition[4];// Expression: lastPcov
                                           //  Referenced by: '<S8>/LastPcov'

    real_T ywt_zero_Value;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S8>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S8>/last_mv'

    real_T UnitDelay2_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S6>/Unit Delay2'

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S6>/Unit Delay3'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S38>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S38>/Normalization Bias'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S38>/InitialOutputs'

    real_T InitialRegressors_Value;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S38>/InitialRegressors'

    real_T Constant_Value;             // Expression: 1e4
                                          //  Referenced by: '<S6>/Constant'

    real_T Constant13_Value[3];        // Expression: G0.D
                                          //  Referenced by: '<S6>/Constant13'

    real_T Constant2_Value;            // Expression: G0.A
                                          //  Referenced by: '<S6>/Constant2'

    real_T Constant12_Value;           // Expression: G0.C
                                          //  Referenced by: '<S6>/Constant12'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S6>/Constant11'

    real_T Constant10_Value;           // Expression: 0
                                          //  Referenced by: '<S6>/Constant10'

    real_T Constant_Value_n;           // Expression: 0
                                          //  Referenced by: '<S7>/Constant'

    real_T X0_Value[2];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S65>/X0'

    real_T ym_zero_Value;              // Expression: zeros(nym,1)
                                          //  Referenced by: '<S8>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/umin_zero'

    real_T umax_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/umax_zero'

    real_T ymin_zero_Value;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/ymin_zero'

    real_T ymax_zero_Value;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/ymax_zero'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S5>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S8>/umin_scale4'

    real_T F_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/F_zero'

    real_T ymin_scale1_Gain;        // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S8>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S8>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S8>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/ecr.wt_zero'

    real_T H_Value[2];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S65>/H'

    real_T G_Value[4];                 // Expression: pInitialization.G
                                          //  Referenced by: '<S65>/G'

    real_T Constant_Value_b;           // Expression: 1
                                          //  Referenced by: '<S2>/Constant'

    real_T P0_Value[4];                // Expression: pInitialization.P0
                                          //  Referenced by: '<S65>/P0'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S8>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S2>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S2>/Measurement Noise'

    real_T MeasurementNoise_Seed;      // Expression: 12345
                                          //  Referenced by: '<S2>/Measurement Noise'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S36>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[46];// Expression: iA
                                             //  Referenced by: '<S8>/Memory'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S56>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S106>/isSqrtUsed'

    boolean_T Constant_Value_i;        // Expression: false()
                                          //  Referenced by: '<S56>/Constant'

  };

  // Parameters for system: '<S218>/MeasurementUpdate'
  struct P_MeasurementUpdate {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S242>/L*(y[k]-yhat[k|k-1])'

  };

  // Parameters for system: '<S1>/State1.controlLaw.AMPC1'
  struct P_State1controlLawAMPC1 {
    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S3>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S3>/yhat'

    real_T params_Y0;                  // Computed Parameter: params_Y0
                                          //  Referenced by: '<S3>/params'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S111>/G_zero'

    real_T LastPcov_InitialCondition[16];// Expression: lastPcov
                                            //  Referenced by: '<S114>/LastPcov'

    real_T ywt_zero_Value[2];          // Expression: zeros(2,1)
                                          //  Referenced by: '<S111>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S114>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S114>/last_mv'

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S112>/Unit Delay3'

    real_T UnitDelay2_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S112>/Unit Delay2'

    real_T Constant2_Value;            // Expression: 0
                                          //  Referenced by: '<S3>/Constant2'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S145>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S145>/Normalization Bias'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S145>/InitialOutputs'

    real_T InitialRegressors_Value;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S145>/InitialRegressors'

    real_T Constant3_Value[2];         // Expression: [1;1e-5]
                                          //  Referenced by: '<S112>/Constant3'

    real_T Constant_Value;             // Expression: 1e4
                                          //  Referenced by: '<S112>/Constant'

    real_T ForgettingFactor_Value_k;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S146>/Forgetting Factor'

    real_T NormalizationBias_Value_m;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S146>/Normalization Bias'

    real_T InitialOutputs_Value_h;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S146>/InitialOutputs'

    real_T InitialRegressors_Value_e;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S146>/InitialRegressors'

    real_T Constant12_Value[4];        // Expression: G1.C
                                          //  Referenced by: '<S112>/Constant12'

    real_T Constant13_Value[6];        // Expression: G1.D
                                          //  Referenced by: '<S112>/Constant13'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S112>/Constant11'

    real_T Constant2_Value_j[2];       // Expression: [0;0]
                                          //  Referenced by: '<S112>/Constant2'

    real_T Constant_Value_d[2];        // Expression: [0;0]
                                          //  Referenced by: '<S113>/Constant'

    real_T X0_Value[4];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S199>/X0'

    real_T ym_zero_Value[2];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S114>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S111>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/umin_zero'

    real_T umax_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/umax_zero'

    real_T ymin_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S111>/ymin_zero'

    real_T ymax_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S111>/ymax_zero'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S111>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S114>/umin_scale4'

    real_T F_zero_Value[2];            // Expression: zeros(1,2)
                                          //  Referenced by: '<S111>/F_zero'

    real_T ymin_scale1_Gain[2];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S114>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S111>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S114>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S111>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S111>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S114>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S111>/ecr.wt_zero'

    real_T P0_Value[16];               // Expression: pInitialization.P0
                                          //  Referenced by: '<S199>/P0'

    real_T Constant_Value_f;           // Expression: 1
                                          //  Referenced by: '<S3>/Constant'

    real_T H_Value[8];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S199>/H'

    real_T G_Value[16];                // Expression: pInitialization.G
                                          //  Referenced by: '<S199>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S114>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_Seed;      // Expression: 12345
                                          //  Referenced by: '<S3>/Measurement Noise'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S142>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[86];// Expression: iA
                                             //  Referenced by: '<S114>/Memory'

    boolean_T Constant1_Value;         // Expression: true
                                          //  Referenced by: '<S3>/Constant1'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S164>/Delay'

    boolean_T Delay_InitialCondition_l;// Expression: true()
                                          //  Referenced by: '<S190>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S240>/isSqrtUsed'

    boolean_T Constant_Value_n;        // Expression: false()
                                          //  Referenced by: '<S164>/Constant'

    boolean_T Constant_Value_a;        // Expression: false()
                                          //  Referenced by: '<S190>/Constant'

    P_MeasurementUpdate MeasurementUpdate_n;// '<S218>/MeasurementUpdate'
  };

  // Parameters for system: '<S1>/State2.controlLaw.AMPC2'
  struct P_State2controlLawAMPC2 {
    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S4>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S4>/yhat'

    real_T params_Y0;                  // Computed Parameter: params_Y0
                                          //  Referenced by: '<S4>/params'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S245>/G_zero'

    real_T LastPcov_InitialCondition[16];// Expression: lastPcov
                                            //  Referenced by: '<S248>/LastPcov'

    real_T ywt_zero_Value[2];          // Expression: zeros(2,1)
                                          //  Referenced by: '<S245>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S248>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S248>/last_mv'

    real_T Constant10_Value[4];        // Expression: G2.A
                                          //  Referenced by: '<S246>/Constant10'

    real_T UnitDelay1_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S246>/Unit Delay1'

    real_T UnitDelay7_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S246>/Unit Delay7'

    real_T Constant2_Value;            // Expression: 0
                                          //  Referenced by: '<S4>/Constant2'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S279>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S279>/Normalization Bias'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S279>/InitialOutputs'

    real_T InitialRegressors_Value;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S279>/InitialRegressors'

    real_T Constant_Value;             // Expression: 1e4
                                          //  Referenced by: '<S246>/Constant'

    real_T ForgettingFactor_Value_a;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S280>/Forgetting Factor'

    real_T NormalizationBias_Value_e;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S280>/Normalization Bias'

    real_T InitialOutputs_Value_a;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S280>/InitialOutputs'

    real_T InitialRegressors_Value_j;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S280>/InitialRegressors'

    real_T Constant4_Value[4];         // Expression: G2.C
                                          //  Referenced by: '<S246>/Constant4'

    real_T Constant5_Value[6];         // Expression: G2.D
                                          //  Referenced by: '<S246>/Constant5'

    real_T Constant3_Value;            // Expression: 1
                                          //  Referenced by: '<S246>/Constant3'

    real_T Constant6_Value[2];         // Expression: [0;0]
                                          //  Referenced by: '<S246>/Constant6'

    real_T Constant_Value_m[2];        // Expression: [0;0]
                                          //  Referenced by: '<S247>/Constant'

    real_T X0_Value[4];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S333>/X0'

    real_T ym_zero_Value[2];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S248>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S245>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/umin_zero'

    real_T umax_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/umax_zero'

    real_T ymin_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S245>/ymin_zero'

    real_T ymax_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S245>/ymax_zero'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S245>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S248>/umin_scale4'

    real_T F_zero_Value[2];            // Expression: zeros(1,2)
                                          //  Referenced by: '<S245>/F_zero'

    real_T ymin_scale1_Gain[2];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S248>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S245>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S248>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S245>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S245>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S248>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S245>/ecr.wt_zero'

    real_T P0_Value[16];               // Expression: pInitialization.P0
                                          //  Referenced by: '<S333>/P0'

    real_T Constant_Value_g;           // Expression: 1
                                          //  Referenced by: '<S4>/Constant'

    real_T H_Value[8];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S333>/H'

    real_T G_Value[16];                // Expression: pInitialization.G
                                          //  Referenced by: '<S333>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S248>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_Seed;      // Expression: 12345
                                          //  Referenced by: '<S4>/Measurement Noise'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S276>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[86];// Expression: iA
                                             //  Referenced by: '<S248>/Memory'

    boolean_T Constant1_Value;         // Expression: true
                                          //  Referenced by: '<S4>/Constant1'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S298>/Delay'

    boolean_T Delay_InitialCondition_c;// Expression: true()
                                          //  Referenced by: '<S324>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S374>/isSqrtUsed'

    boolean_T Constant_Value_l;        // Expression: false()
                                          //  Referenced by: '<S298>/Constant'

    boolean_T Constant_Value_e;        // Expression: false()
                                          //  Referenced by: '<S324>/Constant'

    P_MeasurementUpdate MeasurementUpdate_a;// '<S352>/MeasurementUpdate'
  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod0;                       // Variable: Aod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Aod1[4];                    // Variable: Aod1
                                          //  Referenced by: '<S113>/MATLAB Function'

    real_T Aod2[4];                    // Variable: Aod2
                                          //  Referenced by: '<S247>/MATLAB Function'

    real_T Bod0;                       // Variable: Bod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Bod1[4];                    // Variable: Bod1
                                          //  Referenced by: '<S113>/MATLAB Function'

    real_T Bod2[4];                    // Variable: Bod2
                                          //  Referenced by: '<S247>/MATLAB Function'

    real_T Cod0;                       // Variable: Cod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Cod1[4];                    // Variable: Cod1
                                          //  Referenced by: '<S113>/MATLAB Function'

    real_T Cod2[4];                    // Variable: Cod2
                                          //  Referenced by: '<S247>/MATLAB Function'

    real_T Dmn0;                       // Variable: Dmn0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Dmn1[4];                    // Variable: Dmn1
                                          //  Referenced by: '<S113>/MATLAB Function'

    real_T Dmn2[4];                    // Variable: Dmn2
                                          //  Referenced by: '<S247>/MATLAB Function'

    real_T Dod0;                       // Variable: Dod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Dod1[4];                    // Variable: Dod1
                                          //  Referenced by: '<S113>/MATLAB Function'

    real_T Dod2[4];                    // Variable: Dod2
                                          //  Referenced by: '<S247>/MATLAB Function'

    real_T dt;                         // Variable: dt
                                          //  Referenced by: '<Root>/SupervisoryController'

    P_State2controlLawAMPC2 State2controlLawAMPC2_l;// '<S1>/State2.controlLaw.AMPC2' 
    P_State1controlLawAMPC1 State1controlLawAMPC1_g;// '<S1>/State1.controlLaw.AMPC1' 
    P_State0controlLawAMPC0 State0controlLawAMPC0_n;// '<S1>/State0.controlLaw.AMPC0' 
  };

  // Copy Constructor
  SupervisoryController(SupervisoryController const&) = delete;

  // Assignment Operator
  SupervisoryController& operator= (SupervisoryController const&) & = delete;

  // Move Constructor
  SupervisoryController(SupervisoryController &&) = delete;

  // Move Assignment Operator
  SupervisoryController& operator= (SupervisoryController &&) = delete;

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

  // private member function(s) for subsystem '<S1>/State0.controlLaw.AMPC0'
  void State0controlLawAMPC0_Init(real_T rty_u[3], real_T *rty_yhat, real_T
    rty_params[3], DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0
    *localP);
  static void State0controlLawAMPC0_Disable(DW_State0controlLawAMPC0 *localDW,
    P_State0controlLawAMPC0 *localP);
  void State0controlLawAMPC0(real_T rtu_r, real_T rtu_y, real_T rtu_y0, const
    real_T rtu_u0[3], const real_T rtu_initParam[3], real_T rtu_paramEst, real_T
    rtu_excitationVal, real_T rty_u[3], real_T *rty_yhat, real_T rty_params[3],
    DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0 *localP, P *rtP,
    ZCE_State0controlLawAMPC0 *localZCE);
  real_T xnrm2_g(int32_T n, const real_T x[3], int32_T ix0);
  real_T qrFactor(const real_T A[2], const real_T S[4], real_T Ns);
  void trisolve_m(real_T A, real_T B_1[2]);
  real_T xnrm2_gg(int32_T n, const real_T x[6], int32_T ix0);
  void xgemv_m(int32_T m, int32_T n, const real_T A[6], int32_T ia0, const
               real_T x[6], int32_T ix0, real_T y[2]);
  void xgerc_n(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [2], real_T A[6], int32_T ia0);
  void sqrtMeasurementUpdate(real_T L[4], const real_T H[2], real_T a0, real_T
    K[2]);
  int32_T xpotrf(real_T b_A[16]);
  real_T minimum(const real_T x[4]);
  void trisolve(const real_T b_A[16], real_T b_B[16]);
  real_T norm(const real_T x[4]);
  real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
             real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
             [4], real_T b_A[16], int32_T ia0);
  void KWIKfactor(const real_T b_Ac[184], const int32_T iC[46], int32_T nA,
                  const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                  int32_T n, real_T RLinv[16], real_T *Status);
  real_T maximum(const real_T x[4]);
  void DropConstraint(int32_T kDrop, boolean_T iA[46], int32_T *nA, int32_T iC
                      [46]);
  void qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T f[4],
              const real_T b_Ac[184], const real_T b[46], boolean_T iA[46],
              int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[46],
              int32_T *status);
  void mpcblock_optimizer(const real_T rseq[20], const real_T vseq[21], const
    real_T x[2], const real_T old_u[3], const boolean_T iA[46], const real_T
    b_Mlim[46], real_T b_Mx[92], real_T b_Mu1[138], real_T b_Mv[966], const
    real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16], real_T b_Ac
    [184], const real_T b_Wdu[3], const real_T b_Jm[180], const real_T b_Wu[3],
    const real_T b_I1[180], const real_T b_A[4], const real_T Bu[126], const
    real_T Bv[42], const real_T b_C[2], const real_T Dv[21], const int32_T
    b_Mrows[46], real_T u[3], real_T useq[63], real_T *status, boolean_T iAout
    [46], DW_State0controlLawAMPC0 *localDW);

  // private member function(s) for subsystem '<S145>/ProcessInitialCovariance'
  static void ProcessInitialCovariance(real_T rtu_u, real_T rty_y[25]);

  // private member function(s) for subsystem '<S145>/RLS'
  void RLS(const real_T rtu_H[2], const real_T rtu_H_h[3], real_T rtu_y,
           boolean_T rtu_isEnabled, real_T rtu_adg1, real_T rtu_yBuffer, real_T
           rtu_HBuffer, const real_T rtu_x[5], const real_T rtu_L[25], real_T
           *rty_e, real_T *rty_yBuffer, real_T *rty_HBuffer, real_T rty_x[5],
           real_T rty_L[25], DW_RLS *localDW);
  real_T xnrm2_m(int32_T n, const real_T x[6], int32_T ix0);
  real_T qrFactor_h(const real_T A[5], const real_T S[25], real_T Ns);
  void trisolve_h(real_T A, real_T B_3[5]);
  real_T xnrm2_mw(int32_T n, const real_T x[30], int32_T ix0);
  void xgemv_p(int32_T m, int32_T n, const real_T A[30], int32_T ia0, const
               real_T x[30], int32_T ix0, real_T y[5]);
  void xgerc_d(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [5], real_T A[30], int32_T ia0);
  void sqrtMeasurementUpdate_c(real_T L[25], const real_T H[5], real_T a0,
    real_T K[5]);

  // private member function(s) for subsystem '<S199>/CalculatePL'
  void CalculatePL(const real_T rtu_Ak[16], const real_T rtu_Ck[8], const real_T
                   rtu_Qbark[16], const real_T rtu_Rbark[4], const real_T
                   rtu_Nbark[8], boolean_T rtu_Enablek, const real_T rtu_Pk[16],
                   real_T rty_Mk[8], real_T rty_Lk[8], real_T rty_Zk[16], real_T
                   rty_Pk1[16]);
  void mrdiv(const real_T A[8], const real_T B_4[4], real_T Y[8]);

  // private member function(s) for subsystem '<S240>/SqrtUsedFcn'
  static void SqrtUsedFcn(const real_T rtu_u[16], boolean_T rtu_isSqrtUsed,
    real_T rty_P[16]);

  // private member function(s) for subsystem '<S218>/MeasurementUpdate'
  static void MeasurementUpdate_Init(real_T rty_Lykyhatkk1[4],
    P_MeasurementUpdate *localP);
  static void MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);
  void MeasurementUpdate(boolean_T rtu_Enable, const real_T rtu_Lk[8], const
    real_T rtu_yk[2], const real_T rtu_Ck[8], const real_T rtu_xhatkk1[4], const
    real_T rtu_Dk[6], const real_T rtu_uk[3], real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);

  // private member function(s) for subsystem '<S199>/ReducedQRN'
  static void ReducedQRN(const real_T rtu_G[16], const real_T rtu_H[8], const
    real_T rtu_Q[16], const real_T rtu_R[4], const real_T rtu_N[8], real_T
    rty_Qbar[16], real_T rty_Rbar[4], real_T rty_Nbar[8]);

  // private member function(s) for subsystem '<S199>/ScalarExpansionQ'
  static void ScalarExpansionQ(const real_T rtu_u[16], real_T rty_y[16]);

  // private member function(s) for subsystem '<S199>/ScalarExpansionR'
  static void ScalarExpansionR(const real_T rtu_u[4], real_T rty_y[4]);

  // private member function(s) for subsystem '<S1>/State1.controlLaw.AMPC1'
  void State1controlLawAMPC1_Init(real_T rty_u[3], real_T rty_yhat[2], real_T
    rty_params[6], DW_State1controlLawAMPC1 *localDW, P_State1controlLawAMPC1
    *localP);
  static void State1controlLawAMPC1_Disable(DW_State1controlLawAMPC1 *localDW,
    P_State1controlLawAMPC1 *localP);
  void State1controlLawAMPC1(const real_T rtu_r[2], const real_T rtu_y[2], const
    real_T rtu_y0[2], const real_T rtu_u0[3], const real_T rtu_initParam[6],
    const real_T rtu_paramEst[2], real_T rtu_excitationVal, real_T rty_u[3],
    real_T rty_yhat[2], real_T rty_params[6], DW_State1controlLawAMPC1 *localDW,
    P_State1controlLawAMPC1 *localP, P *rtP, ZCE_State1controlLawAMPC1 *localZCE);
  int32_T xpotrf_n(real_T b_A[16]);
  real_T minimum_i(const real_T x[4]);
  void trisolve_i(const real_T b_A[16], real_T b_B[16]);
  real_T norm_h(const real_T x[4]);
  real_T xnrm2_n(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv_o(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
               real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc_j(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T
               y[4], real_T b_A[16], int32_T ia0);
  void KWIKfactor_h(const real_T b_Ac[344], const int32_T iC[86], int32_T nA,
                    const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                    int32_T n, real_T RLinv[16], real_T *Status);
  real_T maximum_l(const real_T x[4]);
  void DropConstraint_j(int32_T kDrop, boolean_T iA[86], int32_T *nA, int32_T
                        iC[86]);
  void qpkwik_i(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T
                f[4], const real_T b_Ac[344], const real_T b[86], boolean_T iA
                [86], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
                lambda[86], int32_T *status);
  void mpcblock_optimizer_d(const real_T rseq[40], const real_T vseq[21], const
    real_T x[4], const real_T old_u[3], const boolean_T iA[86], const real_T
    b_Mlim[86], real_T b_Mx[344], real_T b_Mu1[258], real_T b_Mv[1806], const
    real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16], real_T b_Ac
    [344], const real_T b_Wy[2], const real_T b_Wdu[3], const real_T b_Jm[180],
    const real_T b_Wu[3], const real_T b_I1[180], const real_T b_A[16], const
    real_T Bu[252], const real_T Bv[84], const real_T b_C[8], const real_T Dv[42],
    const int32_T b_Mrows[86], real_T u[3], real_T useq[63], real_T *status,
    boolean_T iAout[86], DW_State1controlLawAMPC1 *localDW);

  // private member function(s) for subsystem '<S279>/ProcessInitialCovariance'
  static void ProcessInitialCovariance_i(real_T rtu_u, real_T rty_y[9]);

  // private member function(s) for subsystem '<S279>/RLS'
  void RLS_h(const real_T rtu_H[3], real_T rtu_y, boolean_T rtu_isEnabled,
             real_T rtu_adg1, real_T *rty_e, real_T rty_x[3], real_T rty_L[9],
             DW_RLS_k *localDW);
  real_T xnrm2_p(int32_T n, const real_T x[4], int32_T ix0);
  real_T qrFactor_o(const real_T A[3], const real_T S[9], real_T Ns);
  void trisolve_o(real_T A, real_T B_6[3]);
  real_T xnrm2_pa(int32_T n, const real_T x[12], int32_T ix0);
  void xgemv_l(int32_T m, int32_T n, const real_T A[12], int32_T ia0, const
               real_T x[12], int32_T ix0, real_T y[3]);
  void xgerc_i(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [3], real_T A[12], int32_T ia0);
  void sqrtMeasurementUpdate_a(real_T L[9], const real_T H[3], real_T a0, real_T
    K[3]);

  // private member function(s) for subsystem '<S1>/State2.controlLaw.AMPC2'
  void State2controlLawAMPC2_Init(real_T rty_u[3], real_T rty_yhat[2], real_T
    rty_params[6], DW_State2controlLawAMPC2 *localDW, P_State2controlLawAMPC2
    *localP);
  static void State2controlLawAMPC2_Disable(DW_State2controlLawAMPC2 *localDW,
    P_State2controlLawAMPC2 *localP);
  void State2controlLawAMPC2(const real_T rtu_r[2], const real_T rtu_y[2], const
    real_T rtu_y0[2], const real_T rtu_u0[3], const real_T rtu_initParam[6],
    const real_T rtu_paramEst[2], real_T rtu_excitationVal, real_T rty_u[3],
    real_T rty_yhat[2], real_T rty_params[6], DW_State2controlLawAMPC2 *localDW,
    P_State2controlLawAMPC2 *localP, P *rtP, ZCE_State2controlLawAMPC2 *localZCE);
  int32_T xpotrf_b(real_T b_A[16]);
  real_T minimum_k(const real_T x[4]);
  void trisolve_a(const real_T b_A[16], real_T b_B[16]);
  real_T norm_i(const real_T x[4]);
  real_T xnrm2_i(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv_b(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
               real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc_h(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T
               y[4], real_T b_A[16], int32_T ia0);
  void KWIKfactor_g(const real_T b_Ac[344], const int32_T iC[86], int32_T nA,
                    const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                    int32_T n, real_T RLinv[16], real_T *Status);
  real_T maximum_m(const real_T x[4]);
  void DropConstraint_a(int32_T kDrop, boolean_T iA[86], int32_T *nA, int32_T
                        iC[86]);
  void qpkwik_l(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T
                f[4], const real_T b_Ac[344], const real_T b[86], boolean_T iA
                [86], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
                lambda[86], int32_T *status);
  void mpcblock_optimizer_l(const real_T rseq[40], const real_T vseq[21], const
    real_T x[4], const real_T old_u[3], const boolean_T iA[86], const real_T
    b_Mlim[86], real_T b_Mx[344], real_T b_Mu1[258], real_T b_Mv[1806], const
    real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16], real_T b_Ac
    [344], const real_T b_Wy[2], const real_T b_Wdu[3], const real_T b_Jm[180],
    const real_T b_Wu[3], const real_T b_I1[180], const real_T b_A[16], const
    real_T Bu[252], const real_T Bv[84], const real_T b_C[8], const real_T Dv[42],
    const int32_T b_Mrows[86], real_T u[3], real_T useq[63], real_T *status,
    boolean_T iAout[86], DW_State2controlLawAMPC2 *localDW);

  // private member function(s) for subsystem '<Root>'
  boolean_T isequal(const event_bus varargin_1, const event_bus varargin_2);
  void computeProfileParams(real_T i, const real_T wayPoints_data[], const
    int32_T wayPoints_size[2], const real_T Vel_data[], const int32_T *Vel_size,
    real_T *vParam, real_T *aParam, real_T *tAParam, real_T *tFParam);
  void processPolynomialResults(const real_T breakMat_data[], const int32_T
    breakMat_size[2], const real_T coeffMat_data[], const int32_T coeffMat_size
    [2], boolean_T hasMultipleBreaks, cell_wrap_5 breaksCell_data[], int32_T
    *breaksCell_size, cell_wrap_6 coeffCell_data[], int32_T *coeffCell_size);
  void linspace(real_T d2, uint16_T n, real_T y_data[], int32_T y_size[2]);
  void ppval(const s_vjEZ2dxatR8VOmLd9oOqoD *pp, const real_T x_data[], const
             int32_T x_size[2], real_T v_data[], int32_T v_size[2]);
  void generateTrajectoriesFromCoefs(const real_T breaks[4], const real_T
    coeffs_data[], const int32_T coeffs_size[2], real_T dim, const real_T
    t_data[], const int32_T t_size[2], real_T q_data[], int32_T q_size[2],
    real_T qd_data[], int32_T qd_size[2], real_T qdd_data[], int32_T qdd_size[2],
    real_T pp_breaks[6], real_T pp_coefs_data[], int32_T pp_coefs_size[3]);
  void trapveltraj(const real_T wayPoints_data[], const int32_T wayPoints_size[2],
                   uint16_T numSamples, real_T varargin_2, real_T q_data[],
                   int32_T q_size[2]);
  void trajGen(const event_bus event, const real_T y_i[3], real_T trajectory
               [7200], uint16_T *trajectorySize);
  void handleEvent(const event_bus event, boolean_T *inTransitionRegion,
                   boolean_T *eventDone, uint16_T *waypoint, real_T *holdTime)
    const;
  void enter_internal_State1(void);
  void enter_internal_State0(void);
  void State1(const int32_T *sfEvent);
  void chartstep_c6_SupervisoryControl(const int32_T *sfEvent);
};

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S8>/Floor' : Unused code path elimination
//  Block '<S8>/Floor1' : Unused code path elimination
//  Block '<S9>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S10>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S11>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S12>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S13>/Matrix Dimension Check' : Unused code path elimination
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
//  Block '<S29>/Vector Dimension Check' : Unused code path elimination
//  Block '<S30>/Vector Dimension Check' : Unused code path elimination
//  Block '<S31>/Vector Dimension Check' : Unused code path elimination
//  Block '<S32>/Vector Dimension Check' : Unused code path elimination
//  Block '<S33>/Vector Dimension Check' : Unused code path elimination
//  Block '<S34>/Vector Dimension Check' : Unused code path elimination
//  Block '<S8>/last_x' : Unused code path elimination
//  Block '<S35>/Vector Dimension Check' : Unused code path elimination
//  Block '<S8>/useq_scale' : Unused code path elimination
//  Block '<S8>/useq_scale1' : Unused code path elimination
//  Block '<S5>/m_zero' : Unused code path elimination
//  Block '<S5>/p_zero' : Unused code path elimination
//  Block '<S39>/S-Function' : Unused code path elimination
//  Block '<S40>/Dimension' : Unused code path elimination
//  Block '<S42>/Dimension' : Unused code path elimination
//  Block '<S38>/Check Same Ts' : Unused code path elimination
//  Block '<S44>/Output Dimension' : Unused code path elimination
//  Block '<S44>/Regressors Dimension' : Unused code path elimination
//  Block '<S44>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S45>/Data Type Duplicate' : Unused code path elimination
//  Block '<S46>/Data Type Duplicate' : Unused code path elimination
//  Block '<S47>/Data Type Duplicate' : Unused code path elimination
//  Block '<S48>/Data Type Duplicate' : Unused code path elimination
//  Block '<S49>/Data Type Duplicate' : Unused code path elimination
//  Block '<S50>/Data Type Duplicate' : Unused code path elimination
//  Block '<S51>/Data Type Duplicate' : Unused code path elimination
//  Block '<S52>/Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Duplicate' : Unused code path elimination
//  Block '<S54>/Data Type Duplicate' : Unused code path elimination
//  Block '<S63>/S-Function' : Unused code path elimination
//  Block '<S62>/Gain' : Unused code path elimination
//  Block '<S62>/Selector' : Unused code path elimination
//  Block '<S74>/Data Type Duplicate' : Unused code path elimination
//  Block '<S75>/Data Type Duplicate' : Unused code path elimination
//  Block '<S77>/Data Type Duplicate' : Unused code path elimination
//  Block '<S78>/Data Type Duplicate' : Unused code path elimination
//  Block '<S81>/Data Type Duplicate' : Unused code path elimination
//  Block '<S82>/Data Type Duplicate' : Unused code path elimination
//  Block '<S90>/CheckSignalProperties' : Unused code path elimination
//  Block '<S91>/CheckSignalProperties' : Unused code path elimination
//  Block '<S92>/CheckSignalProperties' : Unused code path elimination
//  Block '<S93>/CheckSignalProperties' : Unused code path elimination
//  Block '<S94>/CheckSignalProperties' : Unused code path elimination
//  Block '<S97>/CheckSignalProperties' : Unused code path elimination
//  Block '<S99>/CheckSignalProperties' : Unused code path elimination
//  Block '<S100>/CheckSignalProperties' : Unused code path elimination
//  Block '<S101>/CheckSignalProperties' : Unused code path elimination
//  Block '<S103>/CheckSignalProperties' : Unused code path elimination
//  Block '<S104>/CheckSignalProperties' : Unused code path elimination
//  Block '<S114>/Floor' : Unused code path elimination
//  Block '<S114>/Floor1' : Unused code path elimination
//  Block '<S115>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S116>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S117>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S118>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S119>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S120>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S121>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S122>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S123>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S124>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S125>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S126>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S127>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S128>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S129>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S130>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S131>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S132>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S133>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S134>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S135>/Vector Dimension Check' : Unused code path elimination
//  Block '<S136>/Vector Dimension Check' : Unused code path elimination
//  Block '<S137>/Vector Dimension Check' : Unused code path elimination
//  Block '<S138>/Vector Dimension Check' : Unused code path elimination
//  Block '<S139>/Vector Dimension Check' : Unused code path elimination
//  Block '<S140>/Vector Dimension Check' : Unused code path elimination
//  Block '<S114>/last_x' : Unused code path elimination
//  Block '<S141>/Vector Dimension Check' : Unused code path elimination
//  Block '<S114>/useq_scale' : Unused code path elimination
//  Block '<S114>/useq_scale1' : Unused code path elimination
//  Block '<S111>/m_zero' : Unused code path elimination
//  Block '<S111>/p_zero' : Unused code path elimination
//  Block '<S147>/S-Function' : Unused code path elimination
//  Block '<S148>/Dimension' : Unused code path elimination
//  Block '<S150>/Dimension' : Unused code path elimination
//  Block '<S145>/Check Same Ts' : Unused code path elimination
//  Block '<S152>/Output Dimension' : Unused code path elimination
//  Block '<S152>/Regressors Dimension' : Unused code path elimination
//  Block '<S152>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S153>/Data Type Duplicate' : Unused code path elimination
//  Block '<S154>/Data Type Duplicate' : Unused code path elimination
//  Block '<S155>/Data Type Duplicate' : Unused code path elimination
//  Block '<S156>/Data Type Duplicate' : Unused code path elimination
//  Block '<S157>/Data Type Duplicate' : Unused code path elimination
//  Block '<S158>/Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Duplicate' : Unused code path elimination
//  Block '<S161>/Data Type Duplicate' : Unused code path elimination
//  Block '<S162>/Data Type Duplicate' : Unused code path elimination
//  Block '<S171>/S-Function' : Unused code path elimination
//  Block '<S170>/Gain' : Unused code path elimination
//  Block '<S170>/Selector' : Unused code path elimination
//  Block '<S173>/S-Function' : Unused code path elimination
//  Block '<S174>/Dimension' : Unused code path elimination
//  Block '<S176>/Dimension' : Unused code path elimination
//  Block '<S146>/Check Same Ts' : Unused code path elimination
//  Block '<S178>/Output Dimension' : Unused code path elimination
//  Block '<S178>/Regressors Dimension' : Unused code path elimination
//  Block '<S178>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S179>/Data Type Duplicate' : Unused code path elimination
//  Block '<S180>/Data Type Duplicate' : Unused code path elimination
//  Block '<S181>/Data Type Duplicate' : Unused code path elimination
//  Block '<S182>/Data Type Duplicate' : Unused code path elimination
//  Block '<S183>/Data Type Duplicate' : Unused code path elimination
//  Block '<S184>/Data Type Duplicate' : Unused code path elimination
//  Block '<S185>/Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Duplicate' : Unused code path elimination
//  Block '<S187>/Data Type Duplicate' : Unused code path elimination
//  Block '<S188>/Data Type Duplicate' : Unused code path elimination
//  Block '<S197>/S-Function' : Unused code path elimination
//  Block '<S196>/Gain' : Unused code path elimination
//  Block '<S196>/Selector' : Unused code path elimination
//  Block '<S208>/Data Type Duplicate' : Unused code path elimination
//  Block '<S209>/Data Type Duplicate' : Unused code path elimination
//  Block '<S211>/Data Type Duplicate' : Unused code path elimination
//  Block '<S212>/Data Type Duplicate' : Unused code path elimination
//  Block '<S215>/Data Type Duplicate' : Unused code path elimination
//  Block '<S216>/Data Type Duplicate' : Unused code path elimination
//  Block '<S224>/CheckSignalProperties' : Unused code path elimination
//  Block '<S225>/CheckSignalProperties' : Unused code path elimination
//  Block '<S226>/CheckSignalProperties' : Unused code path elimination
//  Block '<S227>/CheckSignalProperties' : Unused code path elimination
//  Block '<S228>/CheckSignalProperties' : Unused code path elimination
//  Block '<S231>/CheckSignalProperties' : Unused code path elimination
//  Block '<S233>/CheckSignalProperties' : Unused code path elimination
//  Block '<S234>/CheckSignalProperties' : Unused code path elimination
//  Block '<S235>/CheckSignalProperties' : Unused code path elimination
//  Block '<S237>/CheckSignalProperties' : Unused code path elimination
//  Block '<S238>/CheckSignalProperties' : Unused code path elimination
//  Block '<S248>/Floor' : Unused code path elimination
//  Block '<S248>/Floor1' : Unused code path elimination
//  Block '<S249>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S250>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S251>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S252>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S253>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S254>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S255>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S256>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S257>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S258>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S259>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S260>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S261>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S262>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S263>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S264>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S265>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S266>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S267>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S268>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S269>/Vector Dimension Check' : Unused code path elimination
//  Block '<S270>/Vector Dimension Check' : Unused code path elimination
//  Block '<S271>/Vector Dimension Check' : Unused code path elimination
//  Block '<S272>/Vector Dimension Check' : Unused code path elimination
//  Block '<S273>/Vector Dimension Check' : Unused code path elimination
//  Block '<S274>/Vector Dimension Check' : Unused code path elimination
//  Block '<S248>/last_x' : Unused code path elimination
//  Block '<S275>/Vector Dimension Check' : Unused code path elimination
//  Block '<S248>/useq_scale' : Unused code path elimination
//  Block '<S248>/useq_scale1' : Unused code path elimination
//  Block '<S245>/m_zero' : Unused code path elimination
//  Block '<S245>/p_zero' : Unused code path elimination
//  Block '<S281>/S-Function' : Unused code path elimination
//  Block '<S282>/Dimension' : Unused code path elimination
//  Block '<S284>/Dimension' : Unused code path elimination
//  Block '<S279>/Check Same Ts' : Unused code path elimination
//  Block '<S286>/Output Dimension' : Unused code path elimination
//  Block '<S286>/Regressors Dimension' : Unused code path elimination
//  Block '<S286>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S287>/Data Type Duplicate' : Unused code path elimination
//  Block '<S288>/Data Type Duplicate' : Unused code path elimination
//  Block '<S289>/Data Type Duplicate' : Unused code path elimination
//  Block '<S290>/Data Type Duplicate' : Unused code path elimination
//  Block '<S291>/Data Type Duplicate' : Unused code path elimination
//  Block '<S292>/Data Type Duplicate' : Unused code path elimination
//  Block '<S293>/Data Type Duplicate' : Unused code path elimination
//  Block '<S294>/Data Type Duplicate' : Unused code path elimination
//  Block '<S295>/Data Type Duplicate' : Unused code path elimination
//  Block '<S296>/Data Type Duplicate' : Unused code path elimination
//  Block '<S305>/S-Function' : Unused code path elimination
//  Block '<S304>/Gain' : Unused code path elimination
//  Block '<S304>/Selector' : Unused code path elimination
//  Block '<S307>/S-Function' : Unused code path elimination
//  Block '<S308>/Dimension' : Unused code path elimination
//  Block '<S310>/Dimension' : Unused code path elimination
//  Block '<S280>/Check Same Ts' : Unused code path elimination
//  Block '<S312>/Output Dimension' : Unused code path elimination
//  Block '<S312>/Regressors Dimension' : Unused code path elimination
//  Block '<S312>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S313>/Data Type Duplicate' : Unused code path elimination
//  Block '<S314>/Data Type Duplicate' : Unused code path elimination
//  Block '<S315>/Data Type Duplicate' : Unused code path elimination
//  Block '<S316>/Data Type Duplicate' : Unused code path elimination
//  Block '<S317>/Data Type Duplicate' : Unused code path elimination
//  Block '<S318>/Data Type Duplicate' : Unused code path elimination
//  Block '<S319>/Data Type Duplicate' : Unused code path elimination
//  Block '<S320>/Data Type Duplicate' : Unused code path elimination
//  Block '<S321>/Data Type Duplicate' : Unused code path elimination
//  Block '<S322>/Data Type Duplicate' : Unused code path elimination
//  Block '<S331>/S-Function' : Unused code path elimination
//  Block '<S330>/Gain' : Unused code path elimination
//  Block '<S330>/Selector' : Unused code path elimination
//  Block '<S342>/Data Type Duplicate' : Unused code path elimination
//  Block '<S343>/Data Type Duplicate' : Unused code path elimination
//  Block '<S345>/Data Type Duplicate' : Unused code path elimination
//  Block '<S346>/Data Type Duplicate' : Unused code path elimination
//  Block '<S349>/Data Type Duplicate' : Unused code path elimination
//  Block '<S350>/Data Type Duplicate' : Unused code path elimination
//  Block '<S358>/CheckSignalProperties' : Unused code path elimination
//  Block '<S359>/CheckSignalProperties' : Unused code path elimination
//  Block '<S360>/CheckSignalProperties' : Unused code path elimination
//  Block '<S361>/CheckSignalProperties' : Unused code path elimination
//  Block '<S362>/CheckSignalProperties' : Unused code path elimination
//  Block '<S365>/CheckSignalProperties' : Unused code path elimination
//  Block '<S367>/CheckSignalProperties' : Unused code path elimination
//  Block '<S368>/CheckSignalProperties' : Unused code path elimination
//  Block '<S369>/CheckSignalProperties' : Unused code path elimination
//  Block '<S371>/CheckSignalProperties' : Unused code path elimination
//  Block '<S372>/CheckSignalProperties' : Unused code path elimination
//  Block '<S8>/Reshape' : Reshape block reduction
//  Block '<S8>/Reshape1' : Reshape block reduction
//  Block '<S8>/Reshape2' : Reshape block reduction
//  Block '<S8>/Reshape3' : Reshape block reduction
//  Block '<S8>/Reshape4' : Reshape block reduction
//  Block '<S8>/Reshape5' : Reshape block reduction
//  Block '<S47>/Conversion' : Eliminate redundant data type conversion
//  Block '<S49>/Conversion' : Eliminate redundant data type conversion
//  Block '<S51>/Conversion' : Eliminate redundant data type conversion
//  Block '<S52>/Conversion' : Eliminate redundant data type conversion
//  Block '<S53>/Conversion' : Eliminate redundant data type conversion
//  Block '<S54>/Conversion' : Eliminate redundant data type conversion
//  Block '<S77>/Conversion' : Eliminate redundant data type conversion
//  Block '<S81>/Conversion' : Eliminate redundant data type conversion
//  Block '<S84>/Reshape' : Reshape block reduction
//  Block '<S65>/ReshapeX0' : Reshape block reduction
//  Block '<S65>/Reshapeu' : Reshape block reduction
//  Block '<S65>/Reshapexhat' : Reshape block reduction
//  Block '<S65>/Reshapey' : Reshape block reduction
//  Block '<S65>/Reshapeyhat' : Reshape block reduction
//  Block '<S114>/Reshape' : Reshape block reduction
//  Block '<S114>/Reshape1' : Reshape block reduction
//  Block '<S114>/Reshape2' : Reshape block reduction
//  Block '<S114>/Reshape3' : Reshape block reduction
//  Block '<S114>/Reshape4' : Reshape block reduction
//  Block '<S114>/Reshape5' : Reshape block reduction
//  Block '<S155>/Conversion' : Eliminate redundant data type conversion
//  Block '<S157>/Conversion' : Eliminate redundant data type conversion
//  Block '<S159>/Conversion' : Eliminate redundant data type conversion
//  Block '<S160>/Conversion' : Eliminate redundant data type conversion
//  Block '<S161>/Conversion' : Eliminate redundant data type conversion
//  Block '<S162>/Conversion' : Eliminate redundant data type conversion
//  Block '<S145>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S181>/Conversion' : Eliminate redundant data type conversion
//  Block '<S183>/Conversion' : Eliminate redundant data type conversion
//  Block '<S185>/Conversion' : Eliminate redundant data type conversion
//  Block '<S186>/Conversion' : Eliminate redundant data type conversion
//  Block '<S187>/Conversion' : Eliminate redundant data type conversion
//  Block '<S188>/Conversion' : Eliminate redundant data type conversion
//  Block '<S146>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S112>/Reshape' : Reshape block reduction
//  Block '<S211>/Conversion' : Eliminate redundant data type conversion
//  Block '<S215>/Conversion' : Eliminate redundant data type conversion
//  Block '<S218>/Reshape' : Reshape block reduction
//  Block '<S199>/ReshapeX0' : Reshape block reduction
//  Block '<S199>/Reshapeu' : Reshape block reduction
//  Block '<S199>/Reshapexhat' : Reshape block reduction
//  Block '<S199>/Reshapey' : Reshape block reduction
//  Block '<S199>/Reshapeyhat' : Reshape block reduction
//  Block '<S248>/Reshape' : Reshape block reduction
//  Block '<S248>/Reshape1' : Reshape block reduction
//  Block '<S248>/Reshape2' : Reshape block reduction
//  Block '<S248>/Reshape3' : Reshape block reduction
//  Block '<S248>/Reshape4' : Reshape block reduction
//  Block '<S248>/Reshape5' : Reshape block reduction
//  Block '<S289>/Conversion' : Eliminate redundant data type conversion
//  Block '<S291>/Conversion' : Eliminate redundant data type conversion
//  Block '<S293>/Conversion' : Eliminate redundant data type conversion
//  Block '<S294>/Conversion' : Eliminate redundant data type conversion
//  Block '<S295>/Conversion' : Eliminate redundant data type conversion
//  Block '<S296>/Conversion' : Eliminate redundant data type conversion
//  Block '<S279>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S315>/Conversion' : Eliminate redundant data type conversion
//  Block '<S317>/Conversion' : Eliminate redundant data type conversion
//  Block '<S319>/Conversion' : Eliminate redundant data type conversion
//  Block '<S320>/Conversion' : Eliminate redundant data type conversion
//  Block '<S321>/Conversion' : Eliminate redundant data type conversion
//  Block '<S322>/Conversion' : Eliminate redundant data type conversion
//  Block '<S280>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S246>/Reshape' : Reshape block reduction
//  Block '<S345>/Conversion' : Eliminate redundant data type conversion
//  Block '<S349>/Conversion' : Eliminate redundant data type conversion
//  Block '<S352>/Reshape' : Reshape block reduction
//  Block '<S333>/ReshapeX0' : Reshape block reduction
//  Block '<S333>/Reshapeu' : Reshape block reduction
//  Block '<S333>/Reshapexhat' : Reshape block reduction
//  Block '<S333>/Reshapey' : Reshape block reduction
//  Block '<S333>/Reshapeyhat' : Reshape block reduction


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
//  '<S2>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0'
//  '<S3>'   : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1'
//  '<S4>'   : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2'
//  '<S5>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller'
//  '<S6>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator'
//  '<S7>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)'
//  '<S8>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC'
//  '<S9>'   : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S10>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S11>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S12>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S13>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S14>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S15>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S16>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S17>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S18>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S19>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S20>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S21>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S22>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S23>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S24>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S25>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S26>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S27>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S28>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S29>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S30>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S31>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S32>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S33>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S34>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S35>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/moorx'
//  '<S36>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/optimizer'
//  '<S37>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S38>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1'
//  '<S39>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Enable Signal'
//  '<S40>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Covariance'
//  '<S41>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Outputs'
//  '<S42>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Parameters'
//  '<S43>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Regressors'
//  '<S44>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Signals'
//  '<S45>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S46>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S47>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialCovariance'
//  '<S48>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialOutputs'
//  '<S49>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialParameters'
//  '<S50>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialRegressors'
//  '<S51>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - L'
//  '<S52>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - Theta'
//  '<S53>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - bufferH'
//  '<S54>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - buffery'
//  '<S55>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/MultiplyWithTranspose'
//  '<S56>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Process Reset'
//  '<S57>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance'
//  '<S58>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialOutputs'
//  '<S59>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialParameters'
//  '<S60>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialRegressors'
//  '<S61>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/RLS'
//  '<S62>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Reset'
//  '<S63>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Process Reset/Check Reset'
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance/ScalarExpansion'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/MATLAB Function'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S93>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S94>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S95>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S96>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S97>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S98>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S99>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S100>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S101>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S102>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S103>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checku'
//  '<S104>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checky'
//  '<S105>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S106>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S107>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S108>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S109>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S110>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S111>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller'
//  '<S112>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator'
//  '<S113>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)'
//  '<S114>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC'
//  '<S115>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S116>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S117>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S118>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S119>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S120>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S121>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S122>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S123>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S124>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S125>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S126>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S127>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S128>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S129>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S130>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S131>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S132>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S133>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S134>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S135>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S136>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S137>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S138>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S139>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S140>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S141>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/moorx'
//  '<S142>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer'
//  '<S143>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S144>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/MATLAB Function1'
//  '<S145>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1'
//  '<S146>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2'
//  '<S147>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Enable Signal'
//  '<S148>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Covariance'
//  '<S149>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Outputs'
//  '<S150>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Parameters'
//  '<S151>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Regressors'
//  '<S152>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Signals'
//  '<S153>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S154>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S155>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialCovariance'
//  '<S156>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialOutputs'
//  '<S157>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialParameters'
//  '<S158>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialRegressors'
//  '<S159>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - L'
//  '<S160>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - Theta'
//  '<S161>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - bufferH'
//  '<S162>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - buffery'
//  '<S163>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/MultiplyWithTranspose'
//  '<S164>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Process Reset'
//  '<S165>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance'
//  '<S166>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialOutputs'
//  '<S167>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialParameters'
//  '<S168>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialRegressors'
//  '<S169>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/RLS'
//  '<S170>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Reset'
//  '<S171>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Process Reset/Check Reset'
//  '<S172>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance/ScalarExpansion'
//  '<S173>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Enable Signal'
//  '<S174>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Covariance'
//  '<S175>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Outputs'
//  '<S176>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Parameters'
//  '<S177>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Regressors'
//  '<S178>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Signals'
//  '<S179>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - AdaptationParameter1'
//  '<S180>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - AdaptationParameter2'
//  '<S181>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialCovariance'
//  '<S182>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialOutputs'
//  '<S183>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialParameters'
//  '<S184>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialRegressors'
//  '<S185>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - L'
//  '<S186>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - Theta'
//  '<S187>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - bufferH'
//  '<S188>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - buffery'
//  '<S189>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/MultiplyWithTranspose'
//  '<S190>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Process Reset'
//  '<S191>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialCovariance'
//  '<S192>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialOutputs'
//  '<S193>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialParameters'
//  '<S194>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialRegressors'
//  '<S195>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/RLS'
//  '<S196>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Reset'
//  '<S197>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Process Reset/Check Reset'
//  '<S198>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialCovariance/ScalarExpansion'
//  '<S199>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2'
//  '<S200>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/MATLAB Function'
//  '<S201>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S202>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S203>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S204>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S205>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S206>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S207>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S208>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S209>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S210>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S211>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S212>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S213>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S214>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S215>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S216>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S217>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S218>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S219>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S220>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S221>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S222>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S223>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S224>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S225>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S226>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S227>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S228>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S229>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S230>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S231>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S232>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S233>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S234>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S235>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S236>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S237>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checku'
//  '<S238>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checky'
//  '<S239>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S240>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S241>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S242>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S243>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S244>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S245>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller'
//  '<S246>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator'
//  '<S247>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)'
//  '<S248>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC'
//  '<S249>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S250>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S251>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S252>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S253>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S254>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S255>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S256>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S257>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S258>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S259>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S260>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S261>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S262>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S263>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S264>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S265>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S266>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S267>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S268>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S269>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S270>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S271>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S272>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S273>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S274>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S275>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/moorx'
//  '<S276>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer'
//  '<S277>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S278>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/MATLAB Function2'
//  '<S279>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3'
//  '<S280>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4'
//  '<S281>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Enable Signal'
//  '<S282>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Covariance'
//  '<S283>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Outputs'
//  '<S284>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Parameters'
//  '<S285>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Regressors'
//  '<S286>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Signals'
//  '<S287>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - AdaptationParameter1'
//  '<S288>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - AdaptationParameter2'
//  '<S289>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialCovariance'
//  '<S290>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialOutputs'
//  '<S291>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialParameters'
//  '<S292>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialRegressors'
//  '<S293>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - L'
//  '<S294>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - Theta'
//  '<S295>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - bufferH'
//  '<S296>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - buffery'
//  '<S297>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/MultiplyWithTranspose'
//  '<S298>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Process Reset'
//  '<S299>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialCovariance'
//  '<S300>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialOutputs'
//  '<S301>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialParameters'
//  '<S302>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialRegressors'
//  '<S303>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/RLS'
//  '<S304>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Reset'
//  '<S305>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Process Reset/Check Reset'
//  '<S306>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialCovariance/ScalarExpansion'
//  '<S307>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Enable Signal'
//  '<S308>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Covariance'
//  '<S309>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Outputs'
//  '<S310>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Parameters'
//  '<S311>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Regressors'
//  '<S312>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Signals'
//  '<S313>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - AdaptationParameter1'
//  '<S314>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - AdaptationParameter2'
//  '<S315>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialCovariance'
//  '<S316>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialOutputs'
//  '<S317>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialParameters'
//  '<S318>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialRegressors'
//  '<S319>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - L'
//  '<S320>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - Theta'
//  '<S321>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - bufferH'
//  '<S322>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - buffery'
//  '<S323>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/MultiplyWithTranspose'
//  '<S324>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Process Reset'
//  '<S325>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialCovariance'
//  '<S326>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialOutputs'
//  '<S327>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialParameters'
//  '<S328>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialRegressors'
//  '<S329>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/RLS'
//  '<S330>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Reset'
//  '<S331>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Process Reset/Check Reset'
//  '<S332>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialCovariance/ScalarExpansion'
//  '<S333>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2'
//  '<S334>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/MATLAB Function'
//  '<S335>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S336>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S337>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S338>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S339>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S340>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S341>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S342>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S343>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S344>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S345>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S346>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S347>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S348>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S349>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S350>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S351>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S352>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S353>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S354>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S355>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S356>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S357>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S358>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S359>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S360>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S361>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S362>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S363>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S364>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S365>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S366>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S367>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S368>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S369>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S370>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S371>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checku'
//  '<S372>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checky'
//  '<S373>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S374>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S375>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S376>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S377>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S378>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
