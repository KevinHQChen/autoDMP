//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.738
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Thu May  4 15:04:01 2023
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

#ifndef DEFINED_TYPEDEF_FOR_mdl2_bus_
#define DEFINED_TYPEDEF_FOR_mdl2_bus_

struct mdl2_bus
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

#ifndef DEFINED_TYPEDEF_FOR_struct_KFAcLAwxCsqlX6o2IlMd8D_
#define DEFINED_TYPEDEF_FOR_struct_KFAcLAwxCsqlX6o2IlMd8D_

struct struct_KFAcLAwxCsqlX6o2IlMd8D
{
  uint8_T estimationMethod;
  uint16_T na;
  uint16_T nb[3];
  uint16_T nc;
  uint16_T nd;
  uint16_T nf[3];
  uint16_T nk[3];
  uint16_T nParameters;
  uint16_T necessaryDataForPhi[4];
  uint16_T b0PosInPhi[3];
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

// Custom Type definition for MATLAB Function: '<S169>/rarx'
#ifndef struct_c_recursiveEstimation_internal_
#define struct_c_recursiveEstimation_internal_

struct c_recursiveEstimation_internal_
{
  int32_T IteratorPosition;
};

#endif                                // struct_c_recursiveEstimation_internal_

#ifndef struct_d_recursiveEstimation_internal_
#define struct_d_recursiveEstimation_internal_

struct d_recursiveEstimation_internal_
{
  c_recursiveEstimation_internal_ OutputDataIterator;
  c_recursiveEstimation_internal_ InputDataIterator;
};

#endif                                // struct_d_recursiveEstimation_internal_

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
    real_T Product3[2];                // '<S107>/Product3'
    real_T last_mv_DSTATE[3];          // '<S8>/last_mv'
    real_T UnitDelay2_DSTATE[3];       // '<S6>/Unit Delay2'
    real_T delayTheta_DSTATE[4];       // '<S38>/delayTheta'
    real_T delayL_DSTATE[16];          // '<S38>/delayL'
    real_T MemoryX_DSTATE[2];          // '<S64>/MemoryX'
    real_T MemoryP_DSTATE[4];          // '<S64>/MemoryP'
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
    boolean_T icLoad_k;                // '<S64>/MemoryX'
    boolean_T icLoad_f;                // '<S64>/MemoryP'
    boolean_T rlsEstimator_not_empty;  // '<S38>/RLS'
    boolean_T MeasurementUpdate_MODE;  // '<S83>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state for system '<S1>/State0.controlLaw.AMPC0'
  struct ZCE_State0controlLawAMPC0 {
    ZCSigState MemoryX_Reset_ZCE_j;    // '<S64>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_a;    // '<S64>/MemoryP'
  };

  // Block signals and states (default storage) for system '<S169>/rarx'
  struct DW_rarx {
    d_recursiveEstimation_internal_ arxEstimator;// '<S169>/rarx'
    boolean_T arxEstimator_not_empty;  // '<S169>/rarx'
  };

  // Block signals and states (default storage) for system '<S227>/MeasurementUpdate' 
  struct DW_MeasurementUpdate {
    boolean_T MeasurementUpdate_MODE;  // '<S227>/MeasurementUpdate'
  };

  // Block signals and states (default storage) for system '<S1>/State1.controlLaw.AMPC1' 
  struct DW_State1controlLawAMPC1 {
    DW_MeasurementUpdate MeasurementUpdate_n;// '<S227>/MeasurementUpdate'
    DW_rarx sf_rarx_l;                 // '<S198>/rarx'
    DW_rarx sf_rarx;                   // '<S169>/rarx'
    real_T Product3[4];                // '<S251>/Product3'
    real_T last_mv_DSTATE[3];          // '<S113>/last_mv'
    real_T UnitDelay5_DSTATE[3];       // '<S111>/Unit Delay5'
    real_T UnitDelay4_DSTATE[2];       // '<S111>/Unit Delay4'
    real_T delayTheta_DSTATE[4];       // '<S146>/delayTheta'
    real_T delayPhi_DSTATE[7];         // '<S146>/delayPhi'
    real_T delayL_DSTATE[16];          // '<S146>/delayL'
    real_T delayTheta_DSTATE_l[4];     // '<S147>/delayTheta'
    real_T delayPhi_DSTATE_j[7];       // '<S147>/delayPhi'
    real_T delayL_DSTATE_j[16];        // '<S147>/delayL'
    real_T MemoryX_DSTATE[4];          // '<S208>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S208>/MemoryP'
    real_T NextOutput[3];              // '<S3>/Measurement Noise'
    real_T dv[1806];
    real_T b_Hv[840];
    real_T Su[2400];
    real_T delayBufferU_DSTATE;        // '<S146>/delayBufferU'
    real_T delayBufferY_DSTATE;        // '<S146>/delayBufferY'
    real_T delayBufferU_DSTATE_c;      // '<S147>/delayBufferU'
    real_T delayBufferY_DSTATE_j;      // '<S147>/delayBufferY'
    uint32_T RandSeed[3];              // '<S3>/Measurement Noise'
    boolean_T Memory_PreviousInput[86];// '<S113>/Memory'
    boolean_T icLoad;                  // '<S146>/delayTheta'
    boolean_T icLoad_o;                // '<S146>/delayPhi'
    boolean_T icLoad_a;                // '<S146>/delayBufferU'
    boolean_T icLoad_p;                // '<S146>/delayBufferY'
    boolean_T icLoad_n;                // '<S146>/delayL'
    boolean_T icLoad_pm;               // '<S147>/delayTheta'
    boolean_T icLoad_d;                // '<S147>/delayPhi'
    boolean_T icLoad_nq;               // '<S147>/delayBufferU'
    boolean_T icLoad_k;                // '<S147>/delayBufferY'
    boolean_T icLoad_g;                // '<S147>/delayL'
    boolean_T icLoad_f;                // '<S208>/MemoryX'
    boolean_T icLoad_a3;               // '<S208>/MemoryP'
  };

  // Zero-crossing (trigger) state for system '<S1>/State1.controlLaw.AMPC1'
  struct ZCE_State1controlLawAMPC1 {
    ZCSigState MemoryX_Reset_ZCE_b;    // '<S208>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_k;    // '<S208>/MemoryP'
  };

  // Block signals and states (default storage) for system '<S1>/State2.controlLaw.AMPC2' 
  struct DW_State2controlLawAMPC2 {
    DW_MeasurementUpdate MeasurementUpdate_a;// '<S371>/MeasurementUpdate'
    DW_rarx sf_rarx_o;                 // '<S342>/rarx'
    DW_rarx sf_rarx;                   // '<S313>/rarx'
    real_T Product3[4];                // '<S395>/Product3'
    real_T last_mv_DSTATE[3];          // '<S257>/last_mv'
    real_T UnitDelay5_DSTATE[3];       // '<S255>/Unit Delay5'
    real_T UnitDelay4_DSTATE[2];       // '<S255>/Unit Delay4'
    real_T delayTheta_DSTATE[4];       // '<S290>/delayTheta'
    real_T delayPhi_DSTATE[7];         // '<S290>/delayPhi'
    real_T delayL_DSTATE[16];          // '<S290>/delayL'
    real_T delayTheta_DSTATE_h[4];     // '<S291>/delayTheta'
    real_T delayPhi_DSTATE_c[7];       // '<S291>/delayPhi'
    real_T delayL_DSTATE_i[16];        // '<S291>/delayL'
    real_T MemoryX_DSTATE[4];          // '<S352>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S352>/MemoryP'
    real_T NextOutput[3];              // '<S4>/Measurement Noise'
    real_T dv[1806];
    real_T b_Hv[840];
    real_T Su[2400];
    real_T delayBufferU_DSTATE;        // '<S290>/delayBufferU'
    real_T delayBufferY_DSTATE;        // '<S290>/delayBufferY'
    real_T delayBufferU_DSTATE_l;      // '<S291>/delayBufferU'
    real_T delayBufferY_DSTATE_f;      // '<S291>/delayBufferY'
    uint32_T RandSeed[3];              // '<S4>/Measurement Noise'
    boolean_T Memory_PreviousInput[86];// '<S257>/Memory'
    boolean_T icLoad;                  // '<S290>/delayTheta'
    boolean_T icLoad_o;                // '<S290>/delayPhi'
    boolean_T icLoad_e;                // '<S290>/delayBufferU'
    boolean_T icLoad_d;                // '<S290>/delayBufferY'
    boolean_T icLoad_e2;               // '<S290>/delayL'
    boolean_T icLoad_f;                // '<S291>/delayTheta'
    boolean_T icLoad_ej;               // '<S291>/delayPhi'
    boolean_T icLoad_j;                // '<S291>/delayBufferU'
    boolean_T icLoad_i;                // '<S291>/delayBufferY'
    boolean_T icLoad_jb;               // '<S291>/delayL'
    boolean_T icLoad_l;                // '<S352>/MemoryX'
    boolean_T icLoad_h;                // '<S352>/MemoryP'
  };

  // Zero-crossing (trigger) state for system '<S1>/State2.controlLaw.AMPC2'
  struct ZCE_State2controlLawAMPC2 {
    ZCSigState MemoryX_Reset_ZCE;      // '<S352>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S352>/MemoryP'
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
  };

  // Parameters for system: '<S1>/State0.controlLaw.AMPC0'
  struct P_State0controlLawAMPC0 {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S107>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

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

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S6>/Unit Delay3'

    real_T UnitDelay2_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S6>/Unit Delay2'

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

    real_T InitialParameters_Value[4];// Expression: initializationParams.theta0
                                         //  Referenced by: '<S38>/InitialParameters'

    real_T InitialCovariance_Value[16];// Expression: initializationParams.L0
                                          //  Referenced by: '<S38>/InitialCovariance'

    real_T Constant13_Value[3];        // Expression: G0.D
                                          //  Referenced by: '<S6>/Constant13'

    real_T Constant12_Value;           // Expression: G0.C
                                          //  Referenced by: '<S6>/Constant12'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S6>/Constant11'

    real_T Constant10_Value;           // Expression: 0
                                          //  Referenced by: '<S6>/Constant10'

    real_T Constant_Value;             // Expression: 0
                                          //  Referenced by: '<S7>/Constant'

    real_T X0_Value[2];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S64>/X0'

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
                                          //  Referenced by: '<S64>/H'

    real_T G_Value[4];                 // Expression: pInitialization.G
                                          //  Referenced by: '<S64>/G'

    real_T Constant_Value_b;           // Expression: 1
                                          //  Referenced by: '<S2>/Constant'

    real_T P0_Value[4];                // Expression: pInitialization.P0
                                          //  Referenced by: '<S64>/P0'

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
                                          //  Referenced by: '<S105>/isSqrtUsed'

    boolean_T Constant_Value_i;        // Expression: false()
                                          //  Referenced by: '<S56>/Constant'

  };

  // Parameters for system: '<S227>/MeasurementUpdate'
  struct P_MeasurementUpdate {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S251>/L*(y[k]-yhat[k|k-1])'

  };

  // Parameters for system: '<S1>/State1.controlLaw.AMPC1'
  struct P_State1controlLawAMPC1 {
    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S3>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S3>/yhat'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S110>/G_zero'

    real_T LastPcov_InitialCondition[16];// Expression: lastPcov
                                            //  Referenced by: '<S113>/LastPcov'

    real_T ywt_zero_Value[2];          // Expression: zeros(2,1)
                                          //  Referenced by: '<S110>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S113>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S113>/last_mv'

    real_T UnitDelay5_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S111>/Unit Delay5'

    real_T UnitDelay4_InitialCondition[2];// Expression: [0; 0]
                                             //  Referenced by: '<S111>/Unit Delay4'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S146>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S146>/Normalization Bias'

    real_T InitialParameters_Value[4];// Expression: initializationParams.theta0
                                         //  Referenced by: '<S146>/InitialParameters'

    real_T InitialPhiMemory_Value[7];
                                  // Expression: initializationParams.phiMemory0
                                     //  Referenced by: '<S146>/InitialPhiMemory'

    real_T InitialInputs_Value;// Expression: initializationParams.initialInputs
                                  //  Referenced by: '<S146>/InitialInputs'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S146>/InitialOutputs'

    real_T InitialCovariance_Value[16];// Expression: initializationParams.L0
                                          //  Referenced by: '<S146>/InitialCovariance'

    real_T ForgettingFactor_Value_o;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S147>/Forgetting Factor'

    real_T NormalizationBias_Value_e;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S147>/Normalization Bias'

    real_T InitialParameters_Value_n[4];
                                      // Expression: initializationParams.theta0
                                         //  Referenced by: '<S147>/InitialParameters'

    real_T InitialPhiMemory_Value_o[7];
                                  // Expression: initializationParams.phiMemory0
                                     //  Referenced by: '<S147>/InitialPhiMemory'

    real_T InitialInputs_Value_n;
                               // Expression: initializationParams.initialInputs
                                  //  Referenced by: '<S147>/InitialInputs'

    real_T InitialOutputs_Value_n;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S147>/InitialOutputs'

    real_T InitialCovariance_Value_a[16];// Expression: initializationParams.L0
                                            //  Referenced by: '<S147>/InitialCovariance'

    real_T Constant1_Value;            // Expression: 1
                                          //  Referenced by: '<S111>/Constant1'

    real_T Constant9_Value[2];         // Expression: [0;0]
                                          //  Referenced by: '<S111>/Constant9'

    real_T Constant_Value[2];          // Expression: [0;0]
                                          //  Referenced by: '<S112>/Constant'

    real_T X0_Value[4];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S208>/X0'

    real_T ym_zero_Value[2];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S113>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S110>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/umin_zero'

    real_T umax_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/umax_zero'

    real_T ymin_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S110>/ymin_zero'

    real_T ymax_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S110>/ymax_zero'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S110>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S113>/umin_scale4'

    real_T F_zero_Value[2];            // Expression: zeros(1,2)
                                          //  Referenced by: '<S110>/F_zero'

    real_T ymin_scale1_Gain[2];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S113>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S110>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S113>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S110>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S110>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S113>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S110>/ecr.wt_zero'

    real_T P0_Value[16];               // Expression: pInitialization.P0
                                          //  Referenced by: '<S208>/P0'

    real_T Constant_Value_f;           // Expression: 1
                                          //  Referenced by: '<S3>/Constant'

    real_T H_Value[8];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S208>/H'

    real_T G_Value[16];                // Expression: pInitialization.G
                                          //  Referenced by: '<S208>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S113>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_Seed;      // Expression: 12345
                                          //  Referenced by: '<S3>/Measurement Noise'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S141>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[86];// Expression: iA
                                             //  Referenced by: '<S113>/Memory'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S171>/Delay'

    boolean_T Delay_InitialCondition_b;// Expression: true()
                                          //  Referenced by: '<S200>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S249>/isSqrtUsed'

    boolean_T Constant_Value_c;        // Expression: false()
                                          //  Referenced by: '<S171>/Constant'

    boolean_T Constant_Value_m;        // Expression: false()
                                          //  Referenced by: '<S200>/Constant'

    P_MeasurementUpdate MeasurementUpdate_n;// '<S227>/MeasurementUpdate'
  };

  // Parameters for system: '<S1>/State2.controlLaw.AMPC2'
  struct P_State2controlLawAMPC2 {
    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S4>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S4>/yhat'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S254>/G_zero'

    real_T LastPcov_InitialCondition[16];// Expression: lastPcov
                                            //  Referenced by: '<S257>/LastPcov'

    real_T ywt_zero_Value[2];          // Expression: zeros(2,1)
                                          //  Referenced by: '<S254>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S257>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S257>/last_mv'

    real_T UnitDelay5_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S255>/Unit Delay5'

    real_T UnitDelay4_InitialCondition[2];// Expression: [0; 0]
                                             //  Referenced by: '<S255>/Unit Delay4'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S290>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S290>/Normalization Bias'

    real_T InitialParameters_Value[4];// Expression: initializationParams.theta0
                                         //  Referenced by: '<S290>/InitialParameters'

    real_T InitialPhiMemory_Value[7];
                                  // Expression: initializationParams.phiMemory0
                                     //  Referenced by: '<S290>/InitialPhiMemory'

    real_T InitialInputs_Value;// Expression: initializationParams.initialInputs
                                  //  Referenced by: '<S290>/InitialInputs'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S290>/InitialOutputs'

    real_T InitialCovariance_Value[16];// Expression: initializationParams.L0
                                          //  Referenced by: '<S290>/InitialCovariance'

    real_T ForgettingFactor_Value_p;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S291>/Forgetting Factor'

    real_T NormalizationBias_Value_o;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S291>/Normalization Bias'

    real_T InitialParameters_Value_j[4];
                                      // Expression: initializationParams.theta0
                                         //  Referenced by: '<S291>/InitialParameters'

    real_T InitialPhiMemory_Value_n[7];
                                  // Expression: initializationParams.phiMemory0
                                     //  Referenced by: '<S291>/InitialPhiMemory'

    real_T InitialInputs_Value_f;
                               // Expression: initializationParams.initialInputs
                                  //  Referenced by: '<S291>/InitialInputs'

    real_T InitialOutputs_Value_j;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S291>/InitialOutputs'

    real_T InitialCovariance_Value_m[16];// Expression: initializationParams.L0
                                            //  Referenced by: '<S291>/InitialCovariance'

    real_T Constant1_Value;            // Expression: 1
                                          //  Referenced by: '<S255>/Constant1'

    real_T Constant9_Value[2];         // Expression: [0;0]
                                          //  Referenced by: '<S255>/Constant9'

    real_T Constant_Value[2];          // Expression: [0;0]
                                          //  Referenced by: '<S256>/Constant'

    real_T X0_Value[4];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S352>/X0'

    real_T ym_zero_Value[2];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S257>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S254>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/umin_zero'

    real_T umax_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/umax_zero'

    real_T ymin_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S254>/ymin_zero'

    real_T ymax_zero_Value[2];         // Expression: zeros(2,1)
                                          //  Referenced by: '<S254>/ymax_zero'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S254>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S257>/umin_scale4'

    real_T F_zero_Value[2];            // Expression: zeros(1,2)
                                          //  Referenced by: '<S254>/F_zero'

    real_T ymin_scale1_Gain[2];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S257>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S254>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S257>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S254>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S254>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S257>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S254>/ecr.wt_zero'

    real_T P0_Value[16];               // Expression: pInitialization.P0
                                          //  Referenced by: '<S352>/P0'

    real_T Constant_Value_g;           // Expression: 1
                                          //  Referenced by: '<S4>/Constant'

    real_T H_Value[8];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S352>/H'

    real_T G_Value[16];                // Expression: pInitialization.G
                                          //  Referenced by: '<S352>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S257>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_Seed;      // Expression: 12345
                                          //  Referenced by: '<S4>/Measurement Noise'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S285>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[86];// Expression: iA
                                             //  Referenced by: '<S257>/Memory'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S315>/Delay'

    boolean_T Delay_InitialCondition_c;// Expression: true()
                                          //  Referenced by: '<S344>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S393>/isSqrtUsed'

    boolean_T Constant_Value_c;        // Expression: false()
                                          //  Referenced by: '<S315>/Constant'

    boolean_T Constant_Value_d;        // Expression: false()
                                          //  Referenced by: '<S344>/Constant'

    P_MeasurementUpdate MeasurementUpdate_a;// '<S371>/MeasurementUpdate'
  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod0;                       // Variable: Aod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Aod1[4];                    // Variable: Aod1
                                          //  Referenced by: '<S112>/MATLAB Function'

    real_T Aod2[4];                    // Variable: Aod2
                                          //  Referenced by: '<S256>/MATLAB Function'

    real_T Bod0;                       // Variable: Bod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Bod1[4];                    // Variable: Bod1
                                          //  Referenced by: '<S112>/MATLAB Function'

    real_T Bod2[4];                    // Variable: Bod2
                                          //  Referenced by: '<S256>/MATLAB Function'

    real_T Cod0;                       // Variable: Cod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Cod1[4];                    // Variable: Cod1
                                          //  Referenced by: '<S112>/MATLAB Function'

    real_T Cod2[4];                    // Variable: Cod2
                                          //  Referenced by: '<S256>/MATLAB Function'

    real_T Dmn0;                       // Variable: Dmn0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Dmn1[4];                    // Variable: Dmn1
                                          //  Referenced by: '<S112>/MATLAB Function'

    real_T Dmn2[4];                    // Variable: Dmn2
                                          //  Referenced by: '<S256>/MATLAB Function'

    real_T Dod0;                       // Variable: Dod0
                                          //  Referenced by: '<S7>/MATLAB Function'

    real_T Dod1[4];                    // Variable: Dod1
                                          //  Referenced by: '<S112>/MATLAB Function'

    real_T Dod2[4];                    // Variable: Dod2
                                          //  Referenced by: '<S256>/MATLAB Function'

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
  void State0controlLawAMPC0_Init(real_T rty_u[3], real_T *rty_yhat,
    DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0 *localP);
  static void State0controlLawAMPC0_Disable(DW_State0controlLawAMPC0 *localDW,
    P_State0controlLawAMPC0 *localP);
  void State0controlLawAMPC0(real_T rtu_r, real_T rtu_y, real_T rtu_y0, const
    real_T rtu_u0[3], real_T rtu_paramEst, real_T rty_u[3], real_T *rty_yhat,
    DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0 *localP, P *rtP,
    ZCE_State0controlLawAMPC0 *localZCE);
  real_T xnrm2_g(int32_T n, const real_T x[5], int32_T ix0);
  real_T qrFactor(const real_T A[4], const real_T S[16], real_T Ns);
  void trisolve_m(real_T A, real_T B_1[4]);
  real_T xnrm2_gg(int32_T n, const real_T x[20], int32_T ix0);
  void xgemv_m(int32_T m, int32_T n, const real_T A[20], int32_T ia0, const
               real_T x[20], int32_T ix0, real_T y[4]);
  void xgerc_n(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [4], real_T A[20], int32_T ia0);
  void sqrtMeasurementUpdate(real_T L[16], const real_T H[4], real_T a0, real_T
    K[4]);
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

  // private member function(s) for subsystem '<S111>/MATLAB Function'
  static void MATLABFunction(real_T rtu_A1, const real_T rtu_B1[3], real_T
    rtu_C1, const real_T rtu_D1[3], real_T rtu_A2, const real_T rtu_B2[3],
    real_T rtu_C2, const real_T rtu_D2[3], real_T rty_A[4], real_T rty_B[6],
    real_T rty_C[4], real_T rty_D[6]);

  // private member function(s) for subsystem '<S144>/Model Type Converter'
  static void ModelTypeConverter(const real_T rtu_num[6], const real_T rtu_den[2],
    real_T *rty_A, real_T rty_B[3], real_T *rty_C, real_T rty_D[3]);

  // private member function(s) for subsystem '<S169>/rarx'
  void rarx(const real_T rtu_uMeas[3], real_T rtu_yMeas, boolean_T rtu_isEnabled,
            real_T rtu_adg1, real_T rty_A[2], real_T rty_B[6], real_T *rty_e,
            real_T rty_theta[4], real_T rty_L[16], real_T rty_phiMemory[7],
            DW_rarx *localDW);
  real_T xnrm2_j(int32_T n, const real_T x[5], int32_T ix0);
  real_T qrFactor_b(const real_T A[4], const real_T S[16], real_T Ns);
  void trisolve_b(real_T A, real_T B_3[4]);
  real_T xnrm2_ji(int32_T n, const real_T x[20], int32_T ix0);
  void xgemv_mj(int32_T m, int32_T n, const real_T A[20], int32_T ia0, const
                real_T x[20], int32_T ix0, real_T y[4]);
  void xgerc_c(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [4], real_T A[20], int32_T ia0);
  void sqrtMeasurementUpdate_g(real_T L[16], const real_T H[4], real_T a0,
    real_T K[4]);

  // private member function(s) for subsystem '<S208>/CalculatePL'
  void CalculatePL(const real_T rtu_Ak[16], const real_T rtu_Ck[8], const real_T
                   rtu_Qbark[16], const real_T rtu_Rbark[4], const real_T
                   rtu_Nbark[8], boolean_T rtu_Enablek, const real_T rtu_Pk[16],
                   real_T rty_Mk[8], real_T rty_Lk[8], real_T rty_Zk[16], real_T
                   rty_Pk1[16]);
  void mrdiv(const real_T A[8], const real_T B_4[4], real_T Y[8]);

  // private member function(s) for subsystem '<S249>/SqrtUsedFcn'
  static void SqrtUsedFcn(const real_T rtu_u[16], boolean_T rtu_isSqrtUsed,
    real_T rty_P[16]);

  // private member function(s) for subsystem '<S227>/MeasurementUpdate'
  static void MeasurementUpdate_Init(real_T rty_Lykyhatkk1[4],
    P_MeasurementUpdate *localP);
  static void MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);
  void MeasurementUpdate(boolean_T rtu_Enable, const real_T rtu_Lk[8], const
    real_T rtu_yk[2], const real_T rtu_Ck[8], const real_T rtu_xhatkk1[4], const
    real_T rtu_Dk[6], const real_T rtu_uk[3], real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);

  // private member function(s) for subsystem '<S208>/ReducedQRN'
  static void ReducedQRN(const real_T rtu_G[16], const real_T rtu_H[8], const
    real_T rtu_Q[16], const real_T rtu_R[4], const real_T rtu_N[8], real_T
    rty_Qbar[16], real_T rty_Rbar[4], real_T rty_Nbar[8]);

  // private member function(s) for subsystem '<S208>/ScalarExpansionQ'
  static void ScalarExpansionQ(const real_T rtu_u[16], real_T rty_y[16]);

  // private member function(s) for subsystem '<S208>/ScalarExpansionR'
  static void ScalarExpansionR(const real_T rtu_u[4], real_T rty_y[4]);

  // private member function(s) for subsystem '<S1>/State1.controlLaw.AMPC1'
  void State1controlLawAMPC1_Init(real_T rty_u[3], real_T rty_yhat[2],
    DW_State1controlLawAMPC1 *localDW, P_State1controlLawAMPC1 *localP);
  static void State1controlLawAMPC1_Disable(DW_State1controlLawAMPC1 *localDW,
    P_State1controlLawAMPC1 *localP);
  void State1controlLawAMPC1(const real_T rtu_r[2], const real_T rtu_y[2], const
    real_T rtu_y0[2], const real_T rtu_u0[3], real_T rtu_paramEst, real_T rty_u
    [3], real_T rty_yhat[2], DW_State1controlLawAMPC1 *localDW,
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

  // private member function(s) for subsystem '<S1>/State2.controlLaw.AMPC2'
  void State2controlLawAMPC2_Init(real_T rty_u[3], real_T rty_yhat[2],
    DW_State2controlLawAMPC2 *localDW, P_State2controlLawAMPC2 *localP);
  static void State2controlLawAMPC2_Disable(DW_State2controlLawAMPC2 *localDW,
    P_State2controlLawAMPC2 *localP);
  void State2controlLawAMPC2(const real_T rtu_r[2], const real_T rtu_y[2], const
    real_T rtu_y0[2], const real_T rtu_u0[3], real_T rtu_paramEst, real_T rty_u
    [3], real_T rty_yhat[2], DW_State2controlLawAMPC2 *localDW,
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
//  Block '<S73>/Data Type Duplicate' : Unused code path elimination
//  Block '<S74>/Data Type Duplicate' : Unused code path elimination
//  Block '<S76>/Data Type Duplicate' : Unused code path elimination
//  Block '<S77>/Data Type Duplicate' : Unused code path elimination
//  Block '<S80>/Data Type Duplicate' : Unused code path elimination
//  Block '<S81>/Data Type Duplicate' : Unused code path elimination
//  Block '<S89>/CheckSignalProperties' : Unused code path elimination
//  Block '<S90>/CheckSignalProperties' : Unused code path elimination
//  Block '<S91>/CheckSignalProperties' : Unused code path elimination
//  Block '<S92>/CheckSignalProperties' : Unused code path elimination
//  Block '<S93>/CheckSignalProperties' : Unused code path elimination
//  Block '<S96>/CheckSignalProperties' : Unused code path elimination
//  Block '<S98>/CheckSignalProperties' : Unused code path elimination
//  Block '<S99>/CheckSignalProperties' : Unused code path elimination
//  Block '<S100>/CheckSignalProperties' : Unused code path elimination
//  Block '<S102>/CheckSignalProperties' : Unused code path elimination
//  Block '<S103>/CheckSignalProperties' : Unused code path elimination
//  Block '<S113>/Floor' : Unused code path elimination
//  Block '<S113>/Floor1' : Unused code path elimination
//  Block '<S114>/Matrix Dimension Check' : Unused code path elimination
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
//  Block '<S134>/Vector Dimension Check' : Unused code path elimination
//  Block '<S135>/Vector Dimension Check' : Unused code path elimination
//  Block '<S136>/Vector Dimension Check' : Unused code path elimination
//  Block '<S137>/Vector Dimension Check' : Unused code path elimination
//  Block '<S138>/Vector Dimension Check' : Unused code path elimination
//  Block '<S139>/Vector Dimension Check' : Unused code path elimination
//  Block '<S113>/last_x' : Unused code path elimination
//  Block '<S140>/Vector Dimension Check' : Unused code path elimination
//  Block '<S113>/useq_scale' : Unused code path elimination
//  Block '<S113>/useq_scale1' : Unused code path elimination
//  Block '<S110>/m_zero' : Unused code path elimination
//  Block '<S110>/p_zero' : Unused code path elimination
//  Block '<S144>/Data Type Duplicate' : Unused code path elimination
//  Block '<S145>/Data Type Duplicate' : Unused code path elimination
//  Block '<S150>/S-Function' : Unused code path elimination
//  Block '<S146>/Check Same Ts' : Unused code path elimination
//  Block '<S155>/Inputs Dimension' : Unused code path elimination
//  Block '<S155>/Output Dimension' : Unused code path elimination
//  Block '<S155>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S156>/Data Type Duplicate' : Unused code path elimination
//  Block '<S157>/Data Type Duplicate' : Unused code path elimination
//  Block '<S158>/Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Duplicate' : Unused code path elimination
//  Block '<S161>/Data Type Duplicate' : Unused code path elimination
//  Block '<S162>/Data Type Duplicate' : Unused code path elimination
//  Block '<S163>/Data Type Duplicate' : Unused code path elimination
//  Block '<S164>/Data Type Duplicate' : Unused code path elimination
//  Block '<S165>/Data Type Duplicate' : Unused code path elimination
//  Block '<S166>/Data Type Duplicate' : Unused code path elimination
//  Block '<S167>/Data Type Duplicate' : Unused code path elimination
//  Block '<S168>/Data Type Duplicate' : Unused code path elimination
//  Block '<S169>/Constant' : Unused code path elimination
//  Block '<S169>/Probe' : Unused code path elimination
//  Block '<S169>/delayDummy' : Unused code path elimination
//  Block '<S178>/S-Function' : Unused code path elimination
//  Block '<S176>/Gain' : Unused code path elimination
//  Block '<S176>/Selector' : Unused code path elimination
//  Block '<S179>/S-Function' : Unused code path elimination
//  Block '<S147>/Check Same Ts' : Unused code path elimination
//  Block '<S184>/Inputs Dimension' : Unused code path elimination
//  Block '<S184>/Output Dimension' : Unused code path elimination
//  Block '<S184>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S185>/Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Duplicate' : Unused code path elimination
//  Block '<S187>/Data Type Duplicate' : Unused code path elimination
//  Block '<S188>/Data Type Duplicate' : Unused code path elimination
//  Block '<S189>/Data Type Duplicate' : Unused code path elimination
//  Block '<S190>/Data Type Duplicate' : Unused code path elimination
//  Block '<S191>/Data Type Duplicate' : Unused code path elimination
//  Block '<S192>/Data Type Duplicate' : Unused code path elimination
//  Block '<S193>/Data Type Duplicate' : Unused code path elimination
//  Block '<S194>/Data Type Duplicate' : Unused code path elimination
//  Block '<S195>/Data Type Duplicate' : Unused code path elimination
//  Block '<S196>/Data Type Duplicate' : Unused code path elimination
//  Block '<S197>/Data Type Duplicate' : Unused code path elimination
//  Block '<S198>/Constant' : Unused code path elimination
//  Block '<S198>/Probe' : Unused code path elimination
//  Block '<S198>/delayDummy' : Unused code path elimination
//  Block '<S207>/S-Function' : Unused code path elimination
//  Block '<S205>/Gain' : Unused code path elimination
//  Block '<S205>/Selector' : Unused code path elimination
//  Block '<S217>/Data Type Duplicate' : Unused code path elimination
//  Block '<S218>/Data Type Duplicate' : Unused code path elimination
//  Block '<S220>/Data Type Duplicate' : Unused code path elimination
//  Block '<S221>/Data Type Duplicate' : Unused code path elimination
//  Block '<S224>/Data Type Duplicate' : Unused code path elimination
//  Block '<S225>/Data Type Duplicate' : Unused code path elimination
//  Block '<S233>/CheckSignalProperties' : Unused code path elimination
//  Block '<S234>/CheckSignalProperties' : Unused code path elimination
//  Block '<S235>/CheckSignalProperties' : Unused code path elimination
//  Block '<S236>/CheckSignalProperties' : Unused code path elimination
//  Block '<S237>/CheckSignalProperties' : Unused code path elimination
//  Block '<S240>/CheckSignalProperties' : Unused code path elimination
//  Block '<S242>/CheckSignalProperties' : Unused code path elimination
//  Block '<S243>/CheckSignalProperties' : Unused code path elimination
//  Block '<S244>/CheckSignalProperties' : Unused code path elimination
//  Block '<S246>/CheckSignalProperties' : Unused code path elimination
//  Block '<S247>/CheckSignalProperties' : Unused code path elimination
//  Block '<S257>/Floor' : Unused code path elimination
//  Block '<S257>/Floor1' : Unused code path elimination
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
//  Block '<S269>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S270>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S271>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S272>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S273>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S274>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S275>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S276>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S277>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S278>/Vector Dimension Check' : Unused code path elimination
//  Block '<S279>/Vector Dimension Check' : Unused code path elimination
//  Block '<S280>/Vector Dimension Check' : Unused code path elimination
//  Block '<S281>/Vector Dimension Check' : Unused code path elimination
//  Block '<S282>/Vector Dimension Check' : Unused code path elimination
//  Block '<S283>/Vector Dimension Check' : Unused code path elimination
//  Block '<S257>/last_x' : Unused code path elimination
//  Block '<S284>/Vector Dimension Check' : Unused code path elimination
//  Block '<S257>/useq_scale' : Unused code path elimination
//  Block '<S257>/useq_scale1' : Unused code path elimination
//  Block '<S254>/m_zero' : Unused code path elimination
//  Block '<S254>/p_zero' : Unused code path elimination
//  Block '<S288>/Data Type Duplicate' : Unused code path elimination
//  Block '<S289>/Data Type Duplicate' : Unused code path elimination
//  Block '<S294>/S-Function' : Unused code path elimination
//  Block '<S290>/Check Same Ts' : Unused code path elimination
//  Block '<S299>/Inputs Dimension' : Unused code path elimination
//  Block '<S299>/Output Dimension' : Unused code path elimination
//  Block '<S299>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S300>/Data Type Duplicate' : Unused code path elimination
//  Block '<S301>/Data Type Duplicate' : Unused code path elimination
//  Block '<S302>/Data Type Duplicate' : Unused code path elimination
//  Block '<S303>/Data Type Duplicate' : Unused code path elimination
//  Block '<S304>/Data Type Duplicate' : Unused code path elimination
//  Block '<S305>/Data Type Duplicate' : Unused code path elimination
//  Block '<S306>/Data Type Duplicate' : Unused code path elimination
//  Block '<S307>/Data Type Duplicate' : Unused code path elimination
//  Block '<S308>/Data Type Duplicate' : Unused code path elimination
//  Block '<S309>/Data Type Duplicate' : Unused code path elimination
//  Block '<S310>/Data Type Duplicate' : Unused code path elimination
//  Block '<S311>/Data Type Duplicate' : Unused code path elimination
//  Block '<S312>/Data Type Duplicate' : Unused code path elimination
//  Block '<S313>/Constant' : Unused code path elimination
//  Block '<S313>/Probe' : Unused code path elimination
//  Block '<S313>/delayDummy' : Unused code path elimination
//  Block '<S322>/S-Function' : Unused code path elimination
//  Block '<S320>/Gain' : Unused code path elimination
//  Block '<S320>/Selector' : Unused code path elimination
//  Block '<S323>/S-Function' : Unused code path elimination
//  Block '<S291>/Check Same Ts' : Unused code path elimination
//  Block '<S328>/Inputs Dimension' : Unused code path elimination
//  Block '<S328>/Output Dimension' : Unused code path elimination
//  Block '<S328>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S329>/Data Type Duplicate' : Unused code path elimination
//  Block '<S330>/Data Type Duplicate' : Unused code path elimination
//  Block '<S331>/Data Type Duplicate' : Unused code path elimination
//  Block '<S332>/Data Type Duplicate' : Unused code path elimination
//  Block '<S333>/Data Type Duplicate' : Unused code path elimination
//  Block '<S334>/Data Type Duplicate' : Unused code path elimination
//  Block '<S335>/Data Type Duplicate' : Unused code path elimination
//  Block '<S336>/Data Type Duplicate' : Unused code path elimination
//  Block '<S337>/Data Type Duplicate' : Unused code path elimination
//  Block '<S338>/Data Type Duplicate' : Unused code path elimination
//  Block '<S339>/Data Type Duplicate' : Unused code path elimination
//  Block '<S340>/Data Type Duplicate' : Unused code path elimination
//  Block '<S341>/Data Type Duplicate' : Unused code path elimination
//  Block '<S342>/Constant' : Unused code path elimination
//  Block '<S342>/Probe' : Unused code path elimination
//  Block '<S342>/delayDummy' : Unused code path elimination
//  Block '<S351>/S-Function' : Unused code path elimination
//  Block '<S349>/Gain' : Unused code path elimination
//  Block '<S349>/Selector' : Unused code path elimination
//  Block '<S361>/Data Type Duplicate' : Unused code path elimination
//  Block '<S362>/Data Type Duplicate' : Unused code path elimination
//  Block '<S364>/Data Type Duplicate' : Unused code path elimination
//  Block '<S365>/Data Type Duplicate' : Unused code path elimination
//  Block '<S368>/Data Type Duplicate' : Unused code path elimination
//  Block '<S369>/Data Type Duplicate' : Unused code path elimination
//  Block '<S377>/CheckSignalProperties' : Unused code path elimination
//  Block '<S378>/CheckSignalProperties' : Unused code path elimination
//  Block '<S379>/CheckSignalProperties' : Unused code path elimination
//  Block '<S380>/CheckSignalProperties' : Unused code path elimination
//  Block '<S381>/CheckSignalProperties' : Unused code path elimination
//  Block '<S384>/CheckSignalProperties' : Unused code path elimination
//  Block '<S386>/CheckSignalProperties' : Unused code path elimination
//  Block '<S387>/CheckSignalProperties' : Unused code path elimination
//  Block '<S388>/CheckSignalProperties' : Unused code path elimination
//  Block '<S390>/CheckSignalProperties' : Unused code path elimination
//  Block '<S391>/CheckSignalProperties' : Unused code path elimination
//  Block '<S8>/Reshape' : Reshape block reduction
//  Block '<S8>/Reshape1' : Reshape block reduction
//  Block '<S8>/Reshape2' : Reshape block reduction
//  Block '<S8>/Reshape3' : Reshape block reduction
//  Block '<S8>/Reshape4' : Reshape block reduction
//  Block '<S8>/Reshape5' : Reshape block reduction
//  Block '<S51>/Conversion' : Eliminate redundant data type conversion
//  Block '<S52>/Conversion' : Eliminate redundant data type conversion
//  Block '<S53>/Conversion' : Eliminate redundant data type conversion
//  Block '<S54>/Conversion' : Eliminate redundant data type conversion
//  Block '<S76>/Conversion' : Eliminate redundant data type conversion
//  Block '<S80>/Conversion' : Eliminate redundant data type conversion
//  Block '<S83>/Reshape' : Reshape block reduction
//  Block '<S64>/ReshapeX0' : Reshape block reduction
//  Block '<S64>/Reshapeu' : Reshape block reduction
//  Block '<S64>/Reshapexhat' : Reshape block reduction
//  Block '<S64>/Reshapey' : Reshape block reduction
//  Block '<S64>/Reshapeyhat' : Reshape block reduction
//  Block '<S113>/Reshape' : Reshape block reduction
//  Block '<S113>/Reshape1' : Reshape block reduction
//  Block '<S113>/Reshape2' : Reshape block reduction
//  Block '<S113>/Reshape3' : Reshape block reduction
//  Block '<S113>/Reshape4' : Reshape block reduction
//  Block '<S113>/Reshape5' : Reshape block reduction
//  Block '<S163>/Conversion' : Eliminate redundant data type conversion
//  Block '<S164>/Conversion' : Eliminate redundant data type conversion
//  Block '<S165>/Conversion' : Eliminate redundant data type conversion
//  Block '<S166>/Conversion' : Eliminate redundant data type conversion
//  Block '<S167>/Conversion' : Eliminate redundant data type conversion
//  Block '<S168>/Conversion' : Eliminate redundant data type conversion
//  Block '<S192>/Conversion' : Eliminate redundant data type conversion
//  Block '<S193>/Conversion' : Eliminate redundant data type conversion
//  Block '<S194>/Conversion' : Eliminate redundant data type conversion
//  Block '<S195>/Conversion' : Eliminate redundant data type conversion
//  Block '<S196>/Conversion' : Eliminate redundant data type conversion
//  Block '<S197>/Conversion' : Eliminate redundant data type conversion
//  Block '<S220>/Conversion' : Eliminate redundant data type conversion
//  Block '<S224>/Conversion' : Eliminate redundant data type conversion
//  Block '<S227>/Reshape' : Reshape block reduction
//  Block '<S208>/ReshapeX0' : Reshape block reduction
//  Block '<S208>/Reshapeu' : Reshape block reduction
//  Block '<S208>/Reshapexhat' : Reshape block reduction
//  Block '<S208>/Reshapey' : Reshape block reduction
//  Block '<S208>/Reshapeyhat' : Reshape block reduction
//  Block '<S257>/Reshape' : Reshape block reduction
//  Block '<S257>/Reshape1' : Reshape block reduction
//  Block '<S257>/Reshape2' : Reshape block reduction
//  Block '<S257>/Reshape3' : Reshape block reduction
//  Block '<S257>/Reshape4' : Reshape block reduction
//  Block '<S257>/Reshape5' : Reshape block reduction
//  Block '<S307>/Conversion' : Eliminate redundant data type conversion
//  Block '<S308>/Conversion' : Eliminate redundant data type conversion
//  Block '<S309>/Conversion' : Eliminate redundant data type conversion
//  Block '<S310>/Conversion' : Eliminate redundant data type conversion
//  Block '<S311>/Conversion' : Eliminate redundant data type conversion
//  Block '<S312>/Conversion' : Eliminate redundant data type conversion
//  Block '<S336>/Conversion' : Eliminate redundant data type conversion
//  Block '<S337>/Conversion' : Eliminate redundant data type conversion
//  Block '<S338>/Conversion' : Eliminate redundant data type conversion
//  Block '<S339>/Conversion' : Eliminate redundant data type conversion
//  Block '<S340>/Conversion' : Eliminate redundant data type conversion
//  Block '<S341>/Conversion' : Eliminate redundant data type conversion
//  Block '<S364>/Conversion' : Eliminate redundant data type conversion
//  Block '<S368>/Conversion' : Eliminate redundant data type conversion
//  Block '<S371>/Reshape' : Reshape block reduction
//  Block '<S352>/ReshapeX0' : Reshape block reduction
//  Block '<S352>/Reshapeu' : Reshape block reduction
//  Block '<S352>/Reshapexhat' : Reshape block reduction
//  Block '<S352>/Reshapey' : Reshape block reduction
//  Block '<S352>/Reshapeyhat' : Reshape block reduction


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
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/MATLAB Function'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S93>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S94>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S95>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S96>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S97>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S98>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S99>'  : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S100>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S101>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S102>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checku'
//  '<S103>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checky'
//  '<S104>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S105>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S106>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S107>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S108>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S109>' : 'T_junction_mpc/SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S110>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller'
//  '<S111>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator'
//  '<S112>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)'
//  '<S113>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC'
//  '<S114>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S115>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S116>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S117>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S118>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S119>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S120>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S121>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S122>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S123>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S124>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S125>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S126>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S127>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S128>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S129>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S130>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S131>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S132>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S133>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S134>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S135>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S136>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S137>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S138>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S139>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S140>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/moorx'
//  '<S141>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer'
//  '<S142>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S143>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/MATLAB Function'
//  '<S144>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Model Type Converter1'
//  '<S145>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Model Type Converter2'
//  '<S146>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1'
//  '<S147>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2'
//  '<S148>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Model Type Converter1/Model Type Converter'
//  '<S149>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Model Type Converter2/Model Type Converter'
//  '<S150>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Enable Signal'
//  '<S151>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Covariance'
//  '<S152>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Inputs'
//  '<S153>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Outputs'
//  '<S154>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Parameters'
//  '<S155>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Check Signals'
//  '<S156>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S157>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S158>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialCovariance'
//  '<S159>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialInputs'
//  '<S160>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialOutputs'
//  '<S161>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialParameters'
//  '<S162>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialPhi'
//  '<S163>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - L'
//  '<S164>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - Phi'
//  '<S165>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - Theta'
//  '<S166>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - uMemory'
//  '<S167>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - yMemory'
//  '<S168>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion Inherited0'
//  '<S169>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Estimator'
//  '<S170>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/MultiplyWithTranspose'
//  '<S171>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Process Reset'
//  '<S172>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialCovariance'
//  '<S173>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialInputs'
//  '<S174>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialOutputs'
//  '<S175>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialParameters'
//  '<S176>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Reset'
//  '<S177>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Estimator/rarx'
//  '<S178>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator1/Process Reset/Check Reset'
//  '<S179>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Enable Signal'
//  '<S180>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Covariance'
//  '<S181>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Inputs'
//  '<S182>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Outputs'
//  '<S183>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Parameters'
//  '<S184>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Check Signals'
//  '<S185>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - AdaptationParameter1'
//  '<S186>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - AdaptationParameter2'
//  '<S187>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialCovariance'
//  '<S188>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialInputs'
//  '<S189>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialOutputs'
//  '<S190>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialParameters'
//  '<S191>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialPhi'
//  '<S192>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - L'
//  '<S193>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - Phi'
//  '<S194>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - Theta'
//  '<S195>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - uMemory'
//  '<S196>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - yMemory'
//  '<S197>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion Inherited0'
//  '<S198>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Estimator'
//  '<S199>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/MultiplyWithTranspose'
//  '<S200>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Process Reset'
//  '<S201>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialCovariance'
//  '<S202>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialInputs'
//  '<S203>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialOutputs'
//  '<S204>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialParameters'
//  '<S205>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Reset'
//  '<S206>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Estimator/rarx'
//  '<S207>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/Recursive Polynomial Model Estimator2/Process Reset/Check Reset'
//  '<S208>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2'
//  '<S209>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/MATLAB Function'
//  '<S210>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S211>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S212>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S213>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S214>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S215>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S216>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S217>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S218>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S219>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S220>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S221>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S222>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S223>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S224>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S225>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S226>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S227>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S228>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S229>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S230>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S231>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S232>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S233>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S234>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S235>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S236>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S237>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S238>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S239>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S240>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S241>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S242>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S243>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S244>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S245>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S246>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checku'
//  '<S247>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checky'
//  '<S248>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S249>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S250>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S251>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S252>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S253>' : 'T_junction_mpc/SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S254>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller'
//  '<S255>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator'
//  '<S256>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)'
//  '<S257>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC'
//  '<S258>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S259>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S260>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S261>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S262>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S263>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S264>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S265>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S266>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S267>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S268>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S269>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S270>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S271>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S272>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S273>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S274>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S275>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S276>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S277>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S278>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S279>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S280>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S281>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S282>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S283>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S284>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/moorx'
//  '<S285>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer'
//  '<S286>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S287>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/MATLAB Function'
//  '<S288>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Model Type Converter1'
//  '<S289>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Model Type Converter2'
//  '<S290>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1'
//  '<S291>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2'
//  '<S292>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Model Type Converter1/Model Type Converter'
//  '<S293>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Model Type Converter2/Model Type Converter'
//  '<S294>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Enable Signal'
//  '<S295>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Covariance'
//  '<S296>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Inputs'
//  '<S297>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Outputs'
//  '<S298>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Initial Parameters'
//  '<S299>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Check Signals'
//  '<S300>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S301>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S302>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialCovariance'
//  '<S303>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialInputs'
//  '<S304>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialOutputs'
//  '<S305>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialParameters'
//  '<S306>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - InitialPhi'
//  '<S307>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - L'
//  '<S308>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - Phi'
//  '<S309>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - Theta'
//  '<S310>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - uMemory'
//  '<S311>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion - yMemory'
//  '<S312>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Data Type Conversion Inherited0'
//  '<S313>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Estimator'
//  '<S314>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/MultiplyWithTranspose'
//  '<S315>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Process Reset'
//  '<S316>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialCovariance'
//  '<S317>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialInputs'
//  '<S318>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialOutputs'
//  '<S319>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/ProcessInitialParameters'
//  '<S320>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Reset'
//  '<S321>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Estimator/rarx'
//  '<S322>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator1/Process Reset/Check Reset'
//  '<S323>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Enable Signal'
//  '<S324>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Covariance'
//  '<S325>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Inputs'
//  '<S326>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Outputs'
//  '<S327>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Initial Parameters'
//  '<S328>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Check Signals'
//  '<S329>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - AdaptationParameter1'
//  '<S330>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - AdaptationParameter2'
//  '<S331>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialCovariance'
//  '<S332>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialInputs'
//  '<S333>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialOutputs'
//  '<S334>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialParameters'
//  '<S335>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - InitialPhi'
//  '<S336>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - L'
//  '<S337>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - Phi'
//  '<S338>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - Theta'
//  '<S339>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - uMemory'
//  '<S340>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion - yMemory'
//  '<S341>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Data Type Conversion Inherited0'
//  '<S342>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Estimator'
//  '<S343>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/MultiplyWithTranspose'
//  '<S344>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Process Reset'
//  '<S345>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialCovariance'
//  '<S346>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialInputs'
//  '<S347>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialOutputs'
//  '<S348>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/ProcessInitialParameters'
//  '<S349>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Reset'
//  '<S350>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Estimator/rarx'
//  '<S351>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/Model Estimator/Recursive Polynomial Model Estimator2/Process Reset/Check Reset'
//  '<S352>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2'
//  '<S353>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/MATLAB Function'
//  '<S354>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S355>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S356>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S357>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S358>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S359>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S360>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S361>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S362>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S363>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S364>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S365>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S366>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S367>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S368>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S369>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S370>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S371>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S372>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S373>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S374>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S375>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S376>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S377>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S378>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S379>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S380>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S381>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S382>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S383>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S384>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S385>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S386>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S387>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S388>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S389>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S390>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checku'
//  '<S391>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checky'
//  '<S392>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S393>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S394>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S395>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S396>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S397>' : 'T_junction_mpc/SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
