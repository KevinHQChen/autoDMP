//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.1222
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May 22 00:10:52 2023
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

// Custom Type definition for MATLAB Function: '<S147>/RLS'
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
  // Block signals and states (default storage) for system '<S147>/RLS'
  struct DW_RLS {
    d_controllib_internal_blocks_rl rlsEstimator;// '<S147>/RLS'
    boolean_T rlsEstimator_not_empty;  // '<S147>/RLS'
  };

  // Block signals and states (default storage) for system '<S220>/MeasurementUpdate' 
  struct DW_MeasurementUpdate {
    boolean_T MeasurementUpdate_MODE;  // '<S220>/MeasurementUpdate'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    DW_MeasurementUpdate MeasurementUpdate_a;// '<S355>/MeasurementUpdate'
    DW_RLS sf_RLS_h;                   // '<S283>/RLS'
    DW_RLS sf_RLS_er;                  // '<S282>/RLS'
    DW_MeasurementUpdate MeasurementUpdate_cg;// '<S220>/MeasurementUpdate'
    DW_RLS sf_RLS_f;                   // '<S148>/RLS'
    DW_RLS sf_RLS_e;                   // '<S147>/RLS'
    d_controllib_internal_blocks_rl rlsEstimator;// '<S39>/RLS'
    real_T Product3[4];                // '<S379>/Product3'
    real_T Product3_p[4];              // '<S244>/Product3'
    real_T Product3_l[2];              // '<S109>/Product3'
    real_T last_mv_DSTATE[3];          // '<S251>/last_mv'
    real_T UnitDelay1_DSTATE[3];       // '<S249>/Unit Delay1'
    real_T UnitDelay7_DSTATE[2];       // '<S249>/Unit Delay7'
    real_T delayTheta_DSTATE[3];       // '<S282>/delayTheta'
    real_T delayL_DSTATE[9];           // '<S282>/delayL'
    real_T delayTheta_DSTATE_k[3];     // '<S283>/delayTheta'
    real_T delayL_DSTATE_h[9];         // '<S283>/delayL'
    real_T MemoryX_DSTATE[4];          // '<S336>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S336>/MemoryP'
    real_T DiscreteFilter_states[118]; // '<S4>/Discrete Filter'
    real_T last_mv_DSTATE_j[3];        // '<S116>/last_mv'
    real_T UnitDelay2_DSTATE[3];       // '<S114>/Unit Delay2'
    real_T UnitDelay3_DSTATE[2];       // '<S114>/Unit Delay3'
    real_T delayTheta_DSTATE_d[3];     // '<S147>/delayTheta'
    real_T delayL_DSTATE_a[9];         // '<S147>/delayL'
    real_T delayTheta_DSTATE_kz[3];    // '<S148>/delayTheta'
    real_T delayL_DSTATE_l[9];         // '<S148>/delayL'
    real_T MemoryX_DSTATE_m[4];        // '<S201>/MemoryX'
    real_T MemoryP_DSTATE_a[16];       // '<S201>/MemoryP'
    real_T DiscreteFilter_states_b[118];// '<S3>/Discrete Filter'
    real_T last_mv_DSTATE_h[3];        // '<S9>/last_mv'
    real_T UnitDelay2_DSTATE_g[3];     // '<S7>/Unit Delay2'
    real_T delayTheta_DSTATE_p[2];     // '<S39>/delayTheta'
    real_T delayL_DSTATE_a2[4];        // '<S39>/delayL'
    real_T MemoryX_DSTATE_a[2];        // '<S66>/MemoryX'
    real_T MemoryP_DSTATE_d[4];        // '<S66>/MemoryP'
    real_T DiscreteFilter_states_d[59];// '<S2>/Discrete Filter'
    real_T traj[7200];                 // '<Root>/SupervisoryController'
    real_T ymax1[2];                   // '<Root>/SupervisoryController'
    real_T ymax2[2];                   // '<Root>/SupervisoryController'
    real_T uclean[3];                  // '<Root>/SupervisoryController'
    real_T B0[3];                      // '<Root>/SupervisoryController'
    real_T B1[6];                      // '<Root>/SupervisoryController'
    real_T B2[6];                      // '<Root>/SupervisoryController'
    real_T NextOutput[3];              // '<S4>/Measurement Noise'
    real_T NextOutput_i[3];            // '<S3>/Measurement Noise'
    real_T NextOutput_j[3];            // '<S2>/Measurement Noise'
    real_T dv[966];
    real_T Su[1200];
    real_T b_Hv[840];
    real_T Su_m[2400];
    real_T b_Hv_c[840];
    real_T Su_k[2400];
    real_T dv1[1806];
    real_T dv2[1806];
    real_T b_data[3597];
    real_T c_data[3597];
    real_T d_data[3597];
    real_T t_data[1199];
    real_T tmp_data[3597];
    real_T delayBuffery_DSTATE;        // '<S282>/delayBuffery'
    real_T delayBufferH_DSTATE;        // '<S282>/delayBufferH'
    real_T delayBuffery_DSTATE_e;      // '<S283>/delayBuffery'
    real_T delayBufferH_DSTATE_o;      // '<S283>/delayBufferH'
    real_T delayBuffery_DSTATE_l;      // '<S147>/delayBuffery'
    real_T delayBufferH_DSTATE_m;      // '<S147>/delayBufferH'
    real_T delayBuffery_DSTATE_m;      // '<S148>/delayBuffery'
    real_T delayBufferH_DSTATE_c;      // '<S148>/delayBufferH'
    real_T UnitDelay3_DSTATE_c;        // '<S7>/Unit Delay3'
    real_T holdT;                      // '<Root>/SupervisoryController'
    real_T ymax0;                      // '<Root>/SupervisoryController'
    uint32_T RandSeed[3];              // '<S4>/Measurement Noise'
    uint32_T RandSeed_o[3];            // '<S3>/Measurement Noise'
    uint32_T RandSeed_h[3];            // '<S2>/Measurement Noise'
    uint16_T waypt;                    // '<Root>/SupervisoryController'
    uint16_T trajSize;                 // '<Root>/SupervisoryController'
    uint8_T is_c6_SupervisoryController;// '<Root>/SupervisoryController'
    uint8_T is_EventHandler;           // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_n;         // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_k;         // '<Root>/SupervisoryController'
    uint8_T is_active_c6_SupervisoryControl;// '<Root>/SupervisoryController'
    boolean_T Memory_PreviousInput[86];// '<S251>/Memory'
    boolean_T Memory_PreviousInput_b[86];// '<S116>/Memory'
    boolean_T Memory_PreviousInput_c[46];// '<S9>/Memory'
    boolean_T evDone;                  // '<Root>/SupervisoryController'
    boolean_T icLoad;                  // '<S282>/delayBuffery'
    boolean_T icLoad_k;                // '<S282>/delayBufferH'
    boolean_T icLoad_n;                // '<S282>/delayTheta'
    boolean_T icLoad_kv;               // '<S282>/delayL'
    boolean_T icLoad_e;                // '<S283>/delayBuffery'
    boolean_T icLoad_k0;               // '<S283>/delayBufferH'
    boolean_T icLoad_p;                // '<S283>/delayTheta'
    boolean_T icLoad_nc;               // '<S283>/delayL'
    boolean_T icLoad_j;                // '<S336>/MemoryX'
    boolean_T icLoad_d;                // '<S336>/MemoryP'
    boolean_T icLoad_a;                // '<S147>/delayBuffery'
    boolean_T icLoad_c;                // '<S147>/delayBufferH'
    boolean_T icLoad_l;                // '<S147>/delayTheta'
    boolean_T icLoad_m;                // '<S147>/delayL'
    boolean_T icLoad_j2;               // '<S148>/delayBuffery'
    boolean_T icLoad_p5;               // '<S148>/delayBufferH'
    boolean_T icLoad_h;                // '<S148>/delayTheta'
    boolean_T icLoad_i;                // '<S148>/delayL'
    boolean_T icLoad_do;               // '<S201>/MemoryX'
    boolean_T icLoad_jj;               // '<S201>/MemoryP'
    boolean_T icLoad_ho;               // '<S39>/delayBuffery'
    boolean_T icLoad_nb;               // '<S39>/delayBufferH'
    boolean_T icLoad_mt;               // '<S39>/delayTheta'
    boolean_T icLoad_et;               // '<S39>/delayL'
    boolean_T icLoad_k0i;              // '<S66>/MemoryX'
    boolean_T icLoad_f;                // '<S66>/MemoryP'
    boolean_T rlsEstimator_not_empty;  // '<S39>/RLS'
    boolean_T MeasurementUpdate_MODE;  // '<S85>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state
  struct PrevZCX {
    ZCSigState SupervisoryController_Trig_ZCE;// '<Root>/SupervisoryController'
    ZCSigState MemoryX_Reset_ZCE;      // '<S336>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S336>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_h;    // '<S201>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_c;    // '<S201>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_j;    // '<S66>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_a;    // '<S66>/MemoryP'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T ymeas[3];                   // '<Root>/y'
    real_T y_max[3];                   // '<Root>/y_max'
    real_T y0[3];                      // '<Root>/y0'
    real_T u0[3];                      // '<Root>/u0'
    event_bus nextEv;                  // '<Root>/nextEv'
    boolean_T enAdapt[3];              // '<Root>/enAdapt'
    real_T excitation;                 // '<Root>/excitation'
    real_T measAvail;                  // '<Root>/measAvail'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T u[3];                       // '<Root>/u'
    event_bus currEv;                  // '<Root>/currEv'
    boolean_T requestEvent;            // '<Root>/requestEvent'
    real_T currTraj[3];                // '<Root>/currTraj'
    real_T yhat[3];                    // '<Root>/yhat'
    real_T B_b[9];                     // '<Root>/B'
    real_T uref[3];                    // '<Root>/uref'
    real_T uoffset[3];                 // '<Root>/uoffset'
    real_T paramEstErr[3];             // '<Root>/paramEstErr'
  };

  // Parameters for system: '<S220>/MeasurementUpdate'
  struct P_MeasurementUpdate {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S244>/L*(y[k]-yhat[k|k-1])'

  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod0;                       // Variable: Aod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Aod1[4];                    // Variable: Aod1
                                          //  Referenced by: '<S115>/MATLAB Function'

    real_T Aod2[4];                    // Variable: Aod2
                                          //  Referenced by: '<S250>/MATLAB Function'

    real_T Bod0;                       // Variable: Bod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Bod1[4];                    // Variable: Bod1
                                          //  Referenced by: '<S115>/MATLAB Function'

    real_T Bod2[4];                    // Variable: Bod2
                                          //  Referenced by: '<S250>/MATLAB Function'

    real_T Cod0;                       // Variable: Cod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Cod1[4];                    // Variable: Cod1
                                          //  Referenced by: '<S115>/MATLAB Function'

    real_T Cod2[4];                    // Variable: Cod2
                                          //  Referenced by: '<S250>/MATLAB Function'

    real_T Dmn0;                       // Variable: Dmn0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dmn1[4];                    // Variable: Dmn1
                                          //  Referenced by: '<S115>/MATLAB Function'

    real_T Dmn2[4];                    // Variable: Dmn2
                                          //  Referenced by: '<S250>/MATLAB Function'

    real_T Dod0;                       // Variable: Dod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dod1[4];                    // Variable: Dod1
                                          //  Referenced by: '<S115>/MATLAB Function'

    real_T Dod2[4];                    // Variable: Dod2
                                          //  Referenced by: '<S250>/MATLAB Function'

    real_T dt;                         // Variable: dt
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T lpfDen;                     // Variable: lpfDen
                                          //  Referenced by:
                                          //    '<S2>/Discrete Filter'
                                          //    '<S3>/Discrete Filter'
                                          //    '<S4>/Discrete Filter'

    real_T lpfNum[60];                 // Variable: lpfNum
                                          //  Referenced by:
                                          //    '<S2>/Discrete Filter'
                                          //    '<S3>/Discrete Filter'
                                          //    '<S4>/Discrete Filter'

    uint16_T chs0;                     // Variable: chs0
                                          //  Referenced by: '<Root>/SupervisoryController'

    uint16_T chs1[2];                  // Variable: chs1
                                          //  Referenced by: '<Root>/SupervisoryController'

    uint16_T chs2[2];                  // Variable: chs2
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S109>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

    real_T params_Y0;                  // Computed Parameter: params_Y0
                                          //  Referenced by: '<S2>/params'

    real_T uref_Y0;                    // Computed Parameter: uref_Y0
                                          //  Referenced by: '<S2>/uref'

    real_T paramErr_Y0;                // Computed Parameter: paramErr_Y0
                                          //  Referenced by: '<S2>/paramErr'

    real_T uclean_Y0;                  // Computed Parameter: uclean_Y0
                                          //  Referenced by: '<S2>/uclean'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/G_zero'

    real_T LastPcov_InitialCondition[4];// Expression: lastPcov
                                           //  Referenced by: '<S9>/LastPcov'

    real_T ywt_zero_Value;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/y.wt_zero'

    real_T uwt_zero_Value[3];          // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/u.wt_zero'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S9>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S9>/last_mv'

    real_T UnitDelay2_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S7>/Unit Delay2'

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S7>/Unit Delay3'

    real_T ForgettingFactor_Value;     // Expression: initializationParams.adg1
                                          //  Referenced by: '<S39>/Forgetting Factor'

    real_T NormalizationBias_Value;    // Expression: initializationParams.adg2
                                          //  Referenced by: '<S39>/Normalization Bias'

    real_T InitialOutputs_Value;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S39>/InitialOutputs'

    real_T InitialRegressors_Value;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S39>/InitialRegressors'

    real_T Constant_Value;             // Expression: 1e4
                                          //  Referenced by: '<S7>/Constant'

    real_T Constant13_Value[3];        // Expression: G0.D
                                          //  Referenced by: '<S7>/Constant13'

    real_T Constant2_Value;            // Expression: G0.A
                                          //  Referenced by: '<S7>/Constant2'

    real_T Constant12_Value;           // Expression: G0.C
                                          //  Referenced by: '<S7>/Constant12'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S7>/Constant11'

    real_T Constant10_Value;           // Expression: 0
                                          //  Referenced by: '<S7>/Constant10'

    real_T Constant_Value_n;           // Expression: 0
                                          //  Referenced by: '<S8>/Constant'

    real_T X0_Value[2];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S66>/X0'

    real_T ym_zero_Value;              // Expression: zeros(nym,1)
                                          //  Referenced by: '<S9>/ym_zero'

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
                                      //  Referenced by: '<S9>/umin_scale4'

    real_T F_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/F_zero'

    real_T ymin_scale1_Gain;        // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S9>/ymin_scale1'

    real_T S_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/S_zero'

    real_T ymin_scale2_Gain;       // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S9>/ymin_scale2'

    real_T switch_zero_Value;          // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/switch_zero'

    real_T mvtarget_zero_Value[3];     // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/mv.target_zero'

    real_T uref_scale_Gain[3];         // Expression: RMVscale
                                          //  Referenced by: '<S9>/uref_scale'

    real_T ecrwt_zero_Value;           // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/ecr.wt_zero'

    real_T H_Value[2];                 // Expression: pInitialization.H
                                          //  Referenced by: '<S66>/H'

    real_T G_Value[4];                 // Expression: pInitialization.G
                                          //  Referenced by: '<S66>/G'

    real_T Constant_Value_b;           // Expression: 1
                                          //  Referenced by: '<S2>/Constant'

    real_T P0_Value[4];                // Expression: pInitialization.P0
                                          //  Referenced by: '<S66>/P0'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S9>/u_scale'

    real_T MeasurementNoise_Mean[3];   // Expression: [0 0 0]
                                          //  Referenced by: '<S2>/Measurement Noise'

    real_T MeasurementNoise_StdDev[3];
                                  // Computed Parameter: MeasurementNoise_StdDev
                                     //  Referenced by: '<S2>/Measurement Noise'

    real_T MeasurementNoise_Seed[3];   // Expression: [12345 12345 12345]
                                          //  Referenced by: '<S2>/Measurement Noise'

    real_T DiscreteFilter_InitialStates;// Expression: 0
                                           //  Referenced by: '<S2>/Discrete Filter'

    real_T Constant1_Value;            // Expression: 0
                                          //  Referenced by: '<S2>/Constant1'

    real_T Switch_Threshold;           // Expression: 0
                                          //  Referenced by: '<S2>/Switch'

    real_T u_Y0_c;                     // Computed Parameter: u_Y0_c
                                          //  Referenced by: '<S3>/u'

    real_T yhat_Y0_f;                  // Computed Parameter: yhat_Y0_f
                                          //  Referenced by: '<S3>/yhat'

    real_T params_Y0_l;                // Computed Parameter: params_Y0_l
                                          //  Referenced by: '<S3>/params'

    real_T uref_Y0_k;                  // Computed Parameter: uref_Y0_k
                                          //  Referenced by: '<S3>/uref'

    real_T paramErr_Y0_k;              // Computed Parameter: paramErr_Y0_k
                                          //  Referenced by: '<S3>/paramErr'

    real_T uclean_Y0_e;                // Computed Parameter: uclean_Y0_e
                                          //  Referenced by: '<S3>/uclean'

    real_T G_zero_Value_m;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S112>/G_zero'

    real_T LastPcov_InitialCondition_a[16];// Expression: lastPcov
                                              //  Referenced by: '<S116>/LastPcov'

    real_T ywt_zero_Value_l[2];        // Expression: zeros(2,1)
                                          //  Referenced by: '<S112>/y.wt_zero'

    real_T uwt_zero_Value_o[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/u.wt_zero'

    real_T duwt_zero_Value_n[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/du.wt_zero'

    real_T extmv_zero_Value_k[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/ext.mv_zero'

    real_T extmv_scale_Gain_g[3];      // Expression: RMVscale
                                          //  Referenced by: '<S116>/ext.mv_scale'

    real_T last_mv_InitialCondition_o[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S116>/last_mv'

    real_T Constant3_Value[4];         // Expression: G1.A
                                          //  Referenced by: '<S114>/Constant3'

    real_T UnitDelay2_InitialCondition_f;// Expression: 0
                                            //  Referenced by: '<S114>/Unit Delay2'

    real_T UnitDelay3_InitialCondition_g;// Expression: 0
                                            //  Referenced by: '<S114>/Unit Delay3'

    real_T ForgettingFactor_Value_h;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S147>/Forgetting Factor'

    real_T NormalizationBias_Value_g;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S147>/Normalization Bias'

    real_T InitialOutputs_Value_j;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S147>/InitialOutputs'

    real_T InitialRegressors_Value_l;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S147>/InitialRegressors'

    real_T Constant_Value_bh;          // Expression: 1e4
                                          //  Referenced by: '<S114>/Constant'

    real_T ForgettingFactor_Value_a;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S148>/Forgetting Factor'

    real_T NormalizationBias_Value_m;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S148>/Normalization Bias'

    real_T InitialOutputs_Value_m;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S148>/InitialOutputs'

    real_T InitialRegressors_Value_o;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S148>/InitialRegressors'

    real_T Constant12_Value_a[4];      // Expression: G1.C
                                          //  Referenced by: '<S114>/Constant12'

    real_T Constant13_Value_m[6];      // Expression: G1.D
                                          //  Referenced by: '<S114>/Constant13'

    real_T Constant11_Value_c;         // Expression: 1
                                          //  Referenced by: '<S114>/Constant11'

    real_T Constant2_Value_c[2];       // Expression: [0;0]
                                          //  Referenced by: '<S114>/Constant2'

    real_T Constant_Value_l[2];        // Expression: [0;0]
                                          //  Referenced by: '<S115>/Constant'

    real_T X0_Value_m[4];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S201>/X0'

    real_T ym_zero_Value_a[2];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S116>/ym_zero'

    real_T md_zero_Value_d;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S112>/md_zero'

    real_T umin_zero_Value_i[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/umin_zero'

    real_T umax_zero_Value_b[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/umax_zero'

    real_T ymin_zero_Value_c[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S112>/ymin_zero'

    real_T ymax_zero_Value_g[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S112>/ymax_zero'

    real_T E_zero_Value_c[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S112>/E_zero'

    real_T umin_scale4_Gain_b[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S116>/umin_scale4'

    real_T F_zero_Value_k[2];          // Expression: zeros(1,2)
                                          //  Referenced by: '<S112>/F_zero'

    real_T ymin_scale1_Gain_h[2];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S116>/ymin_scale1'

    real_T S_zero_Value_b;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S112>/S_zero'

    real_T ymin_scale2_Gain_g;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S116>/ymin_scale2'

    real_T switch_zero_Value_g;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S112>/switch_zero'

    real_T mvtarget_zero_Value_c[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S112>/mv.target_zero'

    real_T uref_scale_Gain_m[3];       // Expression: RMVscale
                                          //  Referenced by: '<S116>/uref_scale'

    real_T ecrwt_zero_Value_m;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S112>/ecr.wt_zero'

    real_T P0_Value_a[16];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S201>/P0'

    real_T Constant_Value_nl;          // Expression: 1
                                          //  Referenced by: '<S3>/Constant'

    real_T H_Value_o[8];               // Expression: pInitialization.H
                                          //  Referenced by: '<S201>/H'

    real_T G_Value_a[16];              // Expression: pInitialization.G
                                          //  Referenced by: '<S201>/G'

    real_T u_scale_Gain_g[3];          // Expression: MVscale
                                          //  Referenced by: '<S116>/u_scale'

    real_T Constant2_Value_j;          // Expression: 0
                                          //  Referenced by: '<S3>/Constant2'

    real_T MeasurementNoise_Mean_o[3]; // Expression: [0 0 0]
                                          //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_StdDev_k[3];
                                // Computed Parameter: MeasurementNoise_StdDev_k
                                   //  Referenced by: '<S3>/Measurement Noise'

    real_T MeasurementNoise_Seed_b[3]; // Expression: [12345 12346 12347]
                                          //  Referenced by: '<S3>/Measurement Noise'

    real_T DiscreteFilter_InitialStates_d;// Expression: 0
                                             //  Referenced by: '<S3>/Discrete Filter'

    real_T Constant3_Value_k;          // Expression: 0
                                          //  Referenced by: '<S3>/Constant3'

    real_T Switch_Threshold_j;         // Expression: 0
                                          //  Referenced by: '<S3>/Switch'

    real_T u_Y0_l;                     // Computed Parameter: u_Y0_l
                                          //  Referenced by: '<S4>/u'

    real_T yhat_Y0_a;                  // Computed Parameter: yhat_Y0_a
                                          //  Referenced by: '<S4>/yhat'

    real_T params_Y0_b;                // Computed Parameter: params_Y0_b
                                          //  Referenced by: '<S4>/params'

    real_T uref_Y0_i;                  // Computed Parameter: uref_Y0_i
                                          //  Referenced by: '<S4>/uref'

    real_T paramErr_Y0_kg;             // Computed Parameter: paramErr_Y0_kg
                                          //  Referenced by: '<S4>/paramErr'

    real_T uclean_Y0_m;                // Computed Parameter: uclean_Y0_m
                                          //  Referenced by: '<S4>/uclean'

    real_T G_zero_Value_c;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S247>/G_zero'

    real_T LastPcov_InitialCondition_d[16];// Expression: lastPcov
                                              //  Referenced by: '<S251>/LastPcov'

    real_T ywt_zero_Value_l3[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S247>/y.wt_zero'

    real_T uwt_zero_Value_b[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/u.wt_zero'

    real_T duwt_zero_Value_p[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/du.wt_zero'

    real_T extmv_zero_Value_m[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/ext.mv_zero'

    real_T extmv_scale_Gain_n[3];      // Expression: RMVscale
                                          //  Referenced by: '<S251>/ext.mv_scale'

    real_T last_mv_InitialCondition_f[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S251>/last_mv'

    real_T Constant10_Value_h[4];      // Expression: G2.A
                                          //  Referenced by: '<S249>/Constant10'

    real_T UnitDelay1_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S249>/Unit Delay1'

    real_T UnitDelay7_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S249>/Unit Delay7'

    real_T Constant2_Value_m;          // Expression: 0
                                          //  Referenced by: '<S4>/Constant2'

    real_T ForgettingFactor_Value_o;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S282>/Forgetting Factor'

    real_T NormalizationBias_Value_b;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S282>/Normalization Bias'

    real_T InitialOutputs_Value_l;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S282>/InitialOutputs'

    real_T InitialRegressors_Value_d;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S282>/InitialRegressors'

    real_T Constant_Value_o;           // Expression: 1e4
                                          //  Referenced by: '<S249>/Constant'

    real_T ForgettingFactor_Value_i;   // Expression: initializationParams.adg1
                                          //  Referenced by: '<S283>/Forgetting Factor'

    real_T NormalizationBias_Value_h;  // Expression: initializationParams.adg2
                                          //  Referenced by: '<S283>/Normalization Bias'

    real_T InitialOutputs_Value_f;
                              // Expression: initializationParams.initialOutputs
                                 //  Referenced by: '<S283>/InitialOutputs'

    real_T InitialRegressors_Value_n;
                           // Expression: initializationParams.initialRegressors
                              //  Referenced by: '<S283>/InitialRegressors'

    real_T Constant4_Value[4];         // Expression: G2.C
                                          //  Referenced by: '<S249>/Constant4'

    real_T Constant5_Value[6];         // Expression: G2.D
                                          //  Referenced by: '<S249>/Constant5'

    real_T Constant3_Value_m;          // Expression: 1
                                          //  Referenced by: '<S249>/Constant3'

    real_T Constant6_Value[2];         // Expression: [0;0]
                                          //  Referenced by: '<S249>/Constant6'

    real_T Constant_Value_h[2];        // Expression: [0;0]
                                          //  Referenced by: '<S250>/Constant'

    real_T X0_Value_a[4];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S336>/X0'

    real_T ym_zero_Value_o[2];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S251>/ym_zero'

    real_T md_zero_Value_c;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S247>/md_zero'

    real_T umin_zero_Value_o[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/umin_zero'

    real_T umax_zero_Value_j[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/umax_zero'

    real_T ymin_zero_Value_f[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S247>/ymin_zero'

    real_T ymax_zero_Value_e[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S247>/ymax_zero'

    real_T E_zero_Value_b[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S247>/E_zero'

    real_T umin_scale4_Gain_p[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S251>/umin_scale4'

    real_T F_zero_Value_b[2];          // Expression: zeros(1,2)
                                          //  Referenced by: '<S247>/F_zero'

    real_T ymin_scale1_Gain_e[2];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S251>/ymin_scale1'

    real_T S_zero_Value_o;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S247>/S_zero'

    real_T ymin_scale2_Gain_gp;    // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S251>/ymin_scale2'

    real_T switch_zero_Value_f;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S247>/switch_zero'

    real_T mvtarget_zero_Value_e[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S247>/mv.target_zero'

    real_T uref_scale_Gain_g[3];       // Expression: RMVscale
                                          //  Referenced by: '<S251>/uref_scale'

    real_T ecrwt_zero_Value_p;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S247>/ecr.wt_zero'

    real_T P0_Value_n[16];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S336>/P0'

    real_T Constant_Value_j;           // Expression: 1
                                          //  Referenced by: '<S4>/Constant'

    real_T H_Value_f[8];               // Expression: pInitialization.H
                                          //  Referenced by: '<S336>/H'

    real_T G_Value_o[16];              // Expression: pInitialization.G
                                          //  Referenced by: '<S336>/G'

    real_T u_scale_Gain_c[3];          // Expression: MVscale
                                          //  Referenced by: '<S251>/u_scale'

    real_T MeasurementNoise_Mean_j[3]; // Expression: [0 0 0]
                                          //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_StdDev_l[3];
                                // Computed Parameter: MeasurementNoise_StdDev_l
                                   //  Referenced by: '<S4>/Measurement Noise'

    real_T MeasurementNoise_Seed_n;    // Expression: 12345
                                          //  Referenced by: '<S4>/Measurement Noise'

    real_T DiscreteFilter_InitialStates_c;// Expression: 0
                                             //  Referenced by: '<S4>/Discrete Filter'

    real_T Constant3_Value_c;          // Expression: 0
                                          //  Referenced by: '<S4>/Constant3'

    real_T Switch_Threshold_i;         // Expression: 0
                                          //  Referenced by: '<S4>/Switch'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S37>/FixedHorizonOptimizer'

    int32_T FixedHorizonOptimizer_Ndis_b;// Expression: Ndis
                                            //  Referenced by: '<S144>/FixedHorizonOptimizer'

    int32_T FixedHorizonOptimizer_Ndis_c;// Expression: Ndis
                                            //  Referenced by: '<S279>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[46];// Expression: iA
                                             //  Referenced by: '<S9>/Memory'

    boolean_T Delay_InitialCondition;  // Expression: true()
                                          //  Referenced by: '<S57>/Delay'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S107>/isSqrtUsed'

    boolean_T Constant_Value_i;        // Expression: false()
                                          //  Referenced by: '<S57>/Constant'

    boolean_T Memory_InitialCondition_p[86];// Expression: iA
                                               //  Referenced by: '<S116>/Memory'

    boolean_T Delay_InitialCondition_f;// Expression: true()
                                          //  Referenced by: '<S166>/Delay'

    boolean_T Delay_InitialCondition_b;// Expression: true()
                                          //  Referenced by: '<S192>/Delay'

    boolean_T isSqrtUsed_Value_l;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S242>/isSqrtUsed'

    boolean_T Constant1_Value_h;       // Expression: true
                                          //  Referenced by: '<S3>/Constant1'

    boolean_T Constant_Value_nq;       // Expression: false()
                                          //  Referenced by: '<S166>/Constant'

    boolean_T Constant_Value_f;        // Expression: false()
                                          //  Referenced by: '<S192>/Constant'

    boolean_T Memory_InitialCondition_l[86];// Expression: iA
                                               //  Referenced by: '<S251>/Memory'

    boolean_T Constant1_Value_k;       // Expression: true
                                          //  Referenced by: '<S4>/Constant1'

    boolean_T Delay_InitialCondition_f1;// Expression: true()
                                           //  Referenced by: '<S301>/Delay'

    boolean_T Delay_InitialCondition_bj;// Expression: true()
                                           //  Referenced by: '<S327>/Delay'

    boolean_T isSqrtUsed_Value_d;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S377>/isSqrtUsed'

    boolean_T Constant_Value_p;        // Expression: false()
                                          //  Referenced by: '<S301>/Constant'

    boolean_T Constant_Value_e;        // Expression: false()
                                          //  Referenced by: '<S327>/Constant'

    P_MeasurementUpdate MeasurementUpdate_a;// '<S355>/MeasurementUpdate'
    P_MeasurementUpdate MeasurementUpdate_cg;// '<S220>/MeasurementUpdate'
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

  // private member function(s) for subsystem '<S3>/MATLAB Function'
  static void MATLABFunction(const real_T rtu_e[2], real_T *rty_decay);

  // private member function(s) for subsystem '<S147>/ProcessInitialCovariance'
  static void ProcessInitialCovariance(real_T rtu_u, real_T rty_y[9]);

  // private member function(s) for subsystem '<S147>/RLS'
  void RLS(const real_T rtu_H[3], real_T rtu_y, boolean_T rtu_isEnabled, real_T
           rtu_adg1, real_T *rty_e, real_T rty_x[3], real_T rty_L[9], DW_RLS
           *localDW);
  real_T xnrm2(int32_T n, const real_T x[4], int32_T ix0);
  real_T qrFactor(const real_T A[3], const real_T S[9], real_T Ns);
  void trisolve(real_T A, real_T B_1[3]);
  real_T xnrm2_c(int32_T n, const real_T x[12], int32_T ix0);
  void xgemv(int32_T m, int32_T n, const real_T A[12], int32_T ia0, const real_T
             x[12], int32_T ix0, real_T y[3]);
  void xgerc(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y[3],
             real_T A[12], int32_T ia0);
  void sqrtMeasurementUpdate(real_T L[9], const real_T H[3], real_T a0, real_T
    K[3]);

  // private member function(s) for subsystem '<S201>/CalculatePL'
  void CalculatePL(const real_T rtu_Ak[16], const real_T rtu_Ck[8], const real_T
                   rtu_Qbark[16], const real_T rtu_Rbark[4], const real_T
                   rtu_Nbark[8], boolean_T rtu_Enablek, const real_T rtu_Pk[16],
                   real_T rty_Mk[8], real_T rty_Lk[8], real_T rty_Zk[16], real_T
                   rty_Pk1[16]);
  void mrdiv(const real_T A[8], const real_T B_2[4], real_T Y[8]);

  // private member function(s) for subsystem '<S242>/SqrtUsedFcn'
  static void SqrtUsedFcn(const real_T rtu_u[16], boolean_T rtu_isSqrtUsed,
    real_T rty_P[16]);

  // private member function(s) for subsystem '<S220>/MeasurementUpdate'
  static void MeasurementUpdate_Init(real_T rty_Lykyhatkk1[4],
    P_MeasurementUpdate *localP);
  static void MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);
  void MeasurementUpdate(boolean_T rtu_Enable, const real_T rtu_Lk[8], const
    real_T rtu_yk[2], const real_T rtu_Ck[8], const real_T rtu_xhatkk1[4], const
    real_T rtu_Dk[6], const real_T rtu_uk[3], real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);

  // private member function(s) for subsystem '<S201>/ReducedQRN'
  static void ReducedQRN(const real_T rtu_G[16], const real_T rtu_H[8], const
    real_T rtu_Q[16], const real_T rtu_R[4], const real_T rtu_N[8], real_T
    rty_Qbar[16], real_T rty_Rbar[4], real_T rty_Nbar[8]);

  // private member function(s) for subsystem '<S201>/ScalarExpansionQ'
  static void ScalarExpansionQ(const real_T rtu_u[16], real_T rty_y[16]);

  // private member function(s) for subsystem '<S201>/ScalarExpansionR'
  static void ScalarExpansionR(const real_T rtu_u[4], real_T rty_y[4]);

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
  void handleEvent(const event_bus event, boolean_T *eventDone, uint16_T
                   *waypoint, real_T *holdTime) const;
  int32_T xpotrf(real_T b_A[16]);
  real_T minimum(const real_T x[4]);
  void trisolve_d(const real_T b_A[16], real_T b_B[16]);
  real_T norm(const real_T x[4]);
  real_T maximum(const real_T x[4]);
  real_T xnrm2_o(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv_g(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
               real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc_h(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T
               y[4], real_T b_A[16], int32_T ia0);
  void KWIKfactor_d(const real_T b_Ac[344], const int32_T iC[86], int32_T nA,
                    const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                    int32_T n, real_T RLinv[16], real_T *Status);
  void DropConstraint_l(int32_T kDrop, boolean_T iA[86], int32_T *nA, int32_T
                        iC[86]);
  void qpkwik_m(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T
                f[4], const real_T b_Ac[344], const real_T b[86], boolean_T iA
                [86], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
                lambda[86], int32_T *status);
  void mpcblock_optimizer_p(const real_T rseq[40], const real_T vseq[21], const
    real_T x[4], const real_T old_u[3], const boolean_T iA[86], const real_T
    b_Mlim[86], real_T b_Mx[344], real_T b_Mu1[258], real_T b_Mv[1806], const
    real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16], real_T b_Ac
    [344], const real_T b_Wy[2], const real_T b_Wdu[3], const real_T b_Jm[180],
    const real_T b_Wu[3], const real_T b_I1[180], const real_T b_A[16], const
    real_T Bu[252], const real_T Bv[84], const real_T b_C[8], const real_T Dv[42],
    const int32_T b_Mrows[86], real_T u[3], real_T useq[63], real_T *status,
    boolean_T iAout[86]);
  void State2(void);
  real_T xnrm2_g(int32_T n, const real_T x[3], int32_T ix0);
  real_T qrFactor_f(const real_T A[2], const real_T S[4], real_T Ns);
  void trisolve_m(real_T A, real_T B_4[2]);
  real_T xnrm2_gg(int32_T n, const real_T x[6], int32_T ix0);
  void xgemv_m(int32_T m, int32_T n, const real_T A[6], int32_T ia0, const
               real_T x[6], int32_T ix0, real_T y[2]);
  void xgerc_n(int32_T m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
               [2], real_T A[6], int32_T ia0);
  void sqrtMeasurementUpdate_i(real_T L[4], const real_T H[2], real_T a0, real_T
    K[2]);
  void KWIKfactor(const real_T b_Ac[184], const int32_T iC[46], int32_T nA,
                  const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                  int32_T n, real_T RLinv[16], real_T *Status);
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
    [46]);
  void State1(void);
};

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S9>/Floor' : Unused code path elimination
//  Block '<S9>/Floor1' : Unused code path elimination
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
//  Block '<S29>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S30>/Vector Dimension Check' : Unused code path elimination
//  Block '<S31>/Vector Dimension Check' : Unused code path elimination
//  Block '<S32>/Vector Dimension Check' : Unused code path elimination
//  Block '<S33>/Vector Dimension Check' : Unused code path elimination
//  Block '<S34>/Vector Dimension Check' : Unused code path elimination
//  Block '<S35>/Vector Dimension Check' : Unused code path elimination
//  Block '<S9>/last_x' : Unused code path elimination
//  Block '<S36>/Vector Dimension Check' : Unused code path elimination
//  Block '<S9>/useq_scale' : Unused code path elimination
//  Block '<S9>/useq_scale1' : Unused code path elimination
//  Block '<S5>/m_zero' : Unused code path elimination
//  Block '<S5>/p_zero' : Unused code path elimination
//  Block '<S40>/S-Function' : Unused code path elimination
//  Block '<S41>/Dimension' : Unused code path elimination
//  Block '<S43>/Dimension' : Unused code path elimination
//  Block '<S39>/Check Same Ts' : Unused code path elimination
//  Block '<S45>/Output Dimension' : Unused code path elimination
//  Block '<S45>/Regressors Dimension' : Unused code path elimination
//  Block '<S45>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S46>/Data Type Duplicate' : Unused code path elimination
//  Block '<S47>/Data Type Duplicate' : Unused code path elimination
//  Block '<S48>/Data Type Duplicate' : Unused code path elimination
//  Block '<S49>/Data Type Duplicate' : Unused code path elimination
//  Block '<S50>/Data Type Duplicate' : Unused code path elimination
//  Block '<S51>/Data Type Duplicate' : Unused code path elimination
//  Block '<S52>/Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Duplicate' : Unused code path elimination
//  Block '<S54>/Data Type Duplicate' : Unused code path elimination
//  Block '<S55>/Data Type Duplicate' : Unused code path elimination
//  Block '<S64>/S-Function' : Unused code path elimination
//  Block '<S63>/Gain' : Unused code path elimination
//  Block '<S63>/Selector' : Unused code path elimination
//  Block '<S75>/Data Type Duplicate' : Unused code path elimination
//  Block '<S76>/Data Type Duplicate' : Unused code path elimination
//  Block '<S78>/Data Type Duplicate' : Unused code path elimination
//  Block '<S79>/Data Type Duplicate' : Unused code path elimination
//  Block '<S82>/Data Type Duplicate' : Unused code path elimination
//  Block '<S83>/Data Type Duplicate' : Unused code path elimination
//  Block '<S91>/CheckSignalProperties' : Unused code path elimination
//  Block '<S92>/CheckSignalProperties' : Unused code path elimination
//  Block '<S93>/CheckSignalProperties' : Unused code path elimination
//  Block '<S94>/CheckSignalProperties' : Unused code path elimination
//  Block '<S95>/CheckSignalProperties' : Unused code path elimination
//  Block '<S98>/CheckSignalProperties' : Unused code path elimination
//  Block '<S100>/CheckSignalProperties' : Unused code path elimination
//  Block '<S101>/CheckSignalProperties' : Unused code path elimination
//  Block '<S102>/CheckSignalProperties' : Unused code path elimination
//  Block '<S104>/CheckSignalProperties' : Unused code path elimination
//  Block '<S105>/CheckSignalProperties' : Unused code path elimination
//  Block '<S116>/Floor' : Unused code path elimination
//  Block '<S116>/Floor1' : Unused code path elimination
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
//  Block '<S135>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S136>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S137>/Vector Dimension Check' : Unused code path elimination
//  Block '<S138>/Vector Dimension Check' : Unused code path elimination
//  Block '<S139>/Vector Dimension Check' : Unused code path elimination
//  Block '<S140>/Vector Dimension Check' : Unused code path elimination
//  Block '<S141>/Vector Dimension Check' : Unused code path elimination
//  Block '<S142>/Vector Dimension Check' : Unused code path elimination
//  Block '<S116>/last_x' : Unused code path elimination
//  Block '<S143>/Vector Dimension Check' : Unused code path elimination
//  Block '<S116>/useq_scale' : Unused code path elimination
//  Block '<S116>/useq_scale1' : Unused code path elimination
//  Block '<S112>/m_zero' : Unused code path elimination
//  Block '<S112>/p_zero' : Unused code path elimination
//  Block '<S149>/S-Function' : Unused code path elimination
//  Block '<S150>/Dimension' : Unused code path elimination
//  Block '<S152>/Dimension' : Unused code path elimination
//  Block '<S147>/Check Same Ts' : Unused code path elimination
//  Block '<S154>/Output Dimension' : Unused code path elimination
//  Block '<S154>/Regressors Dimension' : Unused code path elimination
//  Block '<S154>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S155>/Data Type Duplicate' : Unused code path elimination
//  Block '<S156>/Data Type Duplicate' : Unused code path elimination
//  Block '<S157>/Data Type Duplicate' : Unused code path elimination
//  Block '<S158>/Data Type Duplicate' : Unused code path elimination
//  Block '<S159>/Data Type Duplicate' : Unused code path elimination
//  Block '<S160>/Data Type Duplicate' : Unused code path elimination
//  Block '<S161>/Data Type Duplicate' : Unused code path elimination
//  Block '<S162>/Data Type Duplicate' : Unused code path elimination
//  Block '<S163>/Data Type Duplicate' : Unused code path elimination
//  Block '<S164>/Data Type Duplicate' : Unused code path elimination
//  Block '<S173>/S-Function' : Unused code path elimination
//  Block '<S172>/Gain' : Unused code path elimination
//  Block '<S172>/Selector' : Unused code path elimination
//  Block '<S175>/S-Function' : Unused code path elimination
//  Block '<S176>/Dimension' : Unused code path elimination
//  Block '<S178>/Dimension' : Unused code path elimination
//  Block '<S148>/Check Same Ts' : Unused code path elimination
//  Block '<S180>/Output Dimension' : Unused code path elimination
//  Block '<S180>/Regressors Dimension' : Unused code path elimination
//  Block '<S180>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S181>/Data Type Duplicate' : Unused code path elimination
//  Block '<S182>/Data Type Duplicate' : Unused code path elimination
//  Block '<S183>/Data Type Duplicate' : Unused code path elimination
//  Block '<S184>/Data Type Duplicate' : Unused code path elimination
//  Block '<S185>/Data Type Duplicate' : Unused code path elimination
//  Block '<S186>/Data Type Duplicate' : Unused code path elimination
//  Block '<S187>/Data Type Duplicate' : Unused code path elimination
//  Block '<S188>/Data Type Duplicate' : Unused code path elimination
//  Block '<S189>/Data Type Duplicate' : Unused code path elimination
//  Block '<S190>/Data Type Duplicate' : Unused code path elimination
//  Block '<S199>/S-Function' : Unused code path elimination
//  Block '<S198>/Gain' : Unused code path elimination
//  Block '<S198>/Selector' : Unused code path elimination
//  Block '<S210>/Data Type Duplicate' : Unused code path elimination
//  Block '<S211>/Data Type Duplicate' : Unused code path elimination
//  Block '<S213>/Data Type Duplicate' : Unused code path elimination
//  Block '<S214>/Data Type Duplicate' : Unused code path elimination
//  Block '<S217>/Data Type Duplicate' : Unused code path elimination
//  Block '<S218>/Data Type Duplicate' : Unused code path elimination
//  Block '<S226>/CheckSignalProperties' : Unused code path elimination
//  Block '<S227>/CheckSignalProperties' : Unused code path elimination
//  Block '<S228>/CheckSignalProperties' : Unused code path elimination
//  Block '<S229>/CheckSignalProperties' : Unused code path elimination
//  Block '<S230>/CheckSignalProperties' : Unused code path elimination
//  Block '<S233>/CheckSignalProperties' : Unused code path elimination
//  Block '<S235>/CheckSignalProperties' : Unused code path elimination
//  Block '<S236>/CheckSignalProperties' : Unused code path elimination
//  Block '<S237>/CheckSignalProperties' : Unused code path elimination
//  Block '<S239>/CheckSignalProperties' : Unused code path elimination
//  Block '<S240>/CheckSignalProperties' : Unused code path elimination
//  Block '<S251>/Floor' : Unused code path elimination
//  Block '<S251>/Floor1' : Unused code path elimination
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
//  Block '<S269>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S270>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S271>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S272>/Vector Dimension Check' : Unused code path elimination
//  Block '<S273>/Vector Dimension Check' : Unused code path elimination
//  Block '<S274>/Vector Dimension Check' : Unused code path elimination
//  Block '<S275>/Vector Dimension Check' : Unused code path elimination
//  Block '<S276>/Vector Dimension Check' : Unused code path elimination
//  Block '<S277>/Vector Dimension Check' : Unused code path elimination
//  Block '<S251>/last_x' : Unused code path elimination
//  Block '<S278>/Vector Dimension Check' : Unused code path elimination
//  Block '<S251>/useq_scale' : Unused code path elimination
//  Block '<S251>/useq_scale1' : Unused code path elimination
//  Block '<S247>/m_zero' : Unused code path elimination
//  Block '<S247>/p_zero' : Unused code path elimination
//  Block '<S284>/S-Function' : Unused code path elimination
//  Block '<S285>/Dimension' : Unused code path elimination
//  Block '<S287>/Dimension' : Unused code path elimination
//  Block '<S282>/Check Same Ts' : Unused code path elimination
//  Block '<S289>/Output Dimension' : Unused code path elimination
//  Block '<S289>/Regressors Dimension' : Unused code path elimination
//  Block '<S289>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S290>/Data Type Duplicate' : Unused code path elimination
//  Block '<S291>/Data Type Duplicate' : Unused code path elimination
//  Block '<S292>/Data Type Duplicate' : Unused code path elimination
//  Block '<S293>/Data Type Duplicate' : Unused code path elimination
//  Block '<S294>/Data Type Duplicate' : Unused code path elimination
//  Block '<S295>/Data Type Duplicate' : Unused code path elimination
//  Block '<S296>/Data Type Duplicate' : Unused code path elimination
//  Block '<S297>/Data Type Duplicate' : Unused code path elimination
//  Block '<S298>/Data Type Duplicate' : Unused code path elimination
//  Block '<S299>/Data Type Duplicate' : Unused code path elimination
//  Block '<S308>/S-Function' : Unused code path elimination
//  Block '<S307>/Gain' : Unused code path elimination
//  Block '<S307>/Selector' : Unused code path elimination
//  Block '<S310>/S-Function' : Unused code path elimination
//  Block '<S311>/Dimension' : Unused code path elimination
//  Block '<S313>/Dimension' : Unused code path elimination
//  Block '<S283>/Check Same Ts' : Unused code path elimination
//  Block '<S315>/Output Dimension' : Unused code path elimination
//  Block '<S315>/Regressors Dimension' : Unused code path elimination
//  Block '<S315>/Sample Times and Data Type' : Unused code path elimination
//  Block '<S316>/Data Type Duplicate' : Unused code path elimination
//  Block '<S317>/Data Type Duplicate' : Unused code path elimination
//  Block '<S318>/Data Type Duplicate' : Unused code path elimination
//  Block '<S319>/Data Type Duplicate' : Unused code path elimination
//  Block '<S320>/Data Type Duplicate' : Unused code path elimination
//  Block '<S321>/Data Type Duplicate' : Unused code path elimination
//  Block '<S322>/Data Type Duplicate' : Unused code path elimination
//  Block '<S323>/Data Type Duplicate' : Unused code path elimination
//  Block '<S324>/Data Type Duplicate' : Unused code path elimination
//  Block '<S325>/Data Type Duplicate' : Unused code path elimination
//  Block '<S334>/S-Function' : Unused code path elimination
//  Block '<S333>/Gain' : Unused code path elimination
//  Block '<S333>/Selector' : Unused code path elimination
//  Block '<S345>/Data Type Duplicate' : Unused code path elimination
//  Block '<S346>/Data Type Duplicate' : Unused code path elimination
//  Block '<S348>/Data Type Duplicate' : Unused code path elimination
//  Block '<S349>/Data Type Duplicate' : Unused code path elimination
//  Block '<S352>/Data Type Duplicate' : Unused code path elimination
//  Block '<S353>/Data Type Duplicate' : Unused code path elimination
//  Block '<S361>/CheckSignalProperties' : Unused code path elimination
//  Block '<S362>/CheckSignalProperties' : Unused code path elimination
//  Block '<S363>/CheckSignalProperties' : Unused code path elimination
//  Block '<S364>/CheckSignalProperties' : Unused code path elimination
//  Block '<S365>/CheckSignalProperties' : Unused code path elimination
//  Block '<S368>/CheckSignalProperties' : Unused code path elimination
//  Block '<S370>/CheckSignalProperties' : Unused code path elimination
//  Block '<S371>/CheckSignalProperties' : Unused code path elimination
//  Block '<S372>/CheckSignalProperties' : Unused code path elimination
//  Block '<S374>/CheckSignalProperties' : Unused code path elimination
//  Block '<S375>/CheckSignalProperties' : Unused code path elimination
//  Block '<S9>/Reshape' : Reshape block reduction
//  Block '<S9>/Reshape1' : Reshape block reduction
//  Block '<S9>/Reshape2' : Reshape block reduction
//  Block '<S9>/Reshape3' : Reshape block reduction
//  Block '<S9>/Reshape4' : Reshape block reduction
//  Block '<S9>/Reshape5' : Reshape block reduction
//  Block '<S48>/Conversion' : Eliminate redundant data type conversion
//  Block '<S50>/Conversion' : Eliminate redundant data type conversion
//  Block '<S52>/Conversion' : Eliminate redundant data type conversion
//  Block '<S53>/Conversion' : Eliminate redundant data type conversion
//  Block '<S54>/Conversion' : Eliminate redundant data type conversion
//  Block '<S55>/Conversion' : Eliminate redundant data type conversion
//  Block '<S78>/Conversion' : Eliminate redundant data type conversion
//  Block '<S82>/Conversion' : Eliminate redundant data type conversion
//  Block '<S85>/Reshape' : Reshape block reduction
//  Block '<S66>/ReshapeX0' : Reshape block reduction
//  Block '<S66>/Reshapeu' : Reshape block reduction
//  Block '<S66>/Reshapexhat' : Reshape block reduction
//  Block '<S66>/Reshapey' : Reshape block reduction
//  Block '<S66>/Reshapeyhat' : Reshape block reduction
//  Block '<S116>/Reshape' : Reshape block reduction
//  Block '<S116>/Reshape1' : Reshape block reduction
//  Block '<S116>/Reshape2' : Reshape block reduction
//  Block '<S116>/Reshape3' : Reshape block reduction
//  Block '<S116>/Reshape4' : Reshape block reduction
//  Block '<S116>/Reshape5' : Reshape block reduction
//  Block '<S157>/Conversion' : Eliminate redundant data type conversion
//  Block '<S159>/Conversion' : Eliminate redundant data type conversion
//  Block '<S161>/Conversion' : Eliminate redundant data type conversion
//  Block '<S162>/Conversion' : Eliminate redundant data type conversion
//  Block '<S163>/Conversion' : Eliminate redundant data type conversion
//  Block '<S164>/Conversion' : Eliminate redundant data type conversion
//  Block '<S183>/Conversion' : Eliminate redundant data type conversion
//  Block '<S185>/Conversion' : Eliminate redundant data type conversion
//  Block '<S187>/Conversion' : Eliminate redundant data type conversion
//  Block '<S188>/Conversion' : Eliminate redundant data type conversion
//  Block '<S189>/Conversion' : Eliminate redundant data type conversion
//  Block '<S190>/Conversion' : Eliminate redundant data type conversion
//  Block '<S114>/Reshape' : Reshape block reduction
//  Block '<S213>/Conversion' : Eliminate redundant data type conversion
//  Block '<S217>/Conversion' : Eliminate redundant data type conversion
//  Block '<S220>/Reshape' : Reshape block reduction
//  Block '<S201>/ReshapeX0' : Reshape block reduction
//  Block '<S201>/Reshapeu' : Reshape block reduction
//  Block '<S201>/Reshapexhat' : Reshape block reduction
//  Block '<S201>/Reshapey' : Reshape block reduction
//  Block '<S201>/Reshapeyhat' : Reshape block reduction
//  Block '<S251>/Reshape' : Reshape block reduction
//  Block '<S251>/Reshape1' : Reshape block reduction
//  Block '<S251>/Reshape2' : Reshape block reduction
//  Block '<S251>/Reshape3' : Reshape block reduction
//  Block '<S251>/Reshape4' : Reshape block reduction
//  Block '<S251>/Reshape5' : Reshape block reduction
//  Block '<S292>/Conversion' : Eliminate redundant data type conversion
//  Block '<S294>/Conversion' : Eliminate redundant data type conversion
//  Block '<S296>/Conversion' : Eliminate redundant data type conversion
//  Block '<S297>/Conversion' : Eliminate redundant data type conversion
//  Block '<S298>/Conversion' : Eliminate redundant data type conversion
//  Block '<S299>/Conversion' : Eliminate redundant data type conversion
//  Block '<S282>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S318>/Conversion' : Eliminate redundant data type conversion
//  Block '<S320>/Conversion' : Eliminate redundant data type conversion
//  Block '<S322>/Conversion' : Eliminate redundant data type conversion
//  Block '<S323>/Conversion' : Eliminate redundant data type conversion
//  Block '<S324>/Conversion' : Eliminate redundant data type conversion
//  Block '<S325>/Conversion' : Eliminate redundant data type conversion
//  Block '<S283>/DataTypeConversionEnable' : Eliminate redundant data type conversion
//  Block '<S249>/Reshape' : Reshape block reduction
//  Block '<S348>/Conversion' : Eliminate redundant data type conversion
//  Block '<S352>/Conversion' : Eliminate redundant data type conversion
//  Block '<S355>/Reshape' : Reshape block reduction
//  Block '<S336>/ReshapeX0' : Reshape block reduction
//  Block '<S336>/Reshapeu' : Reshape block reduction
//  Block '<S336>/Reshapexhat' : Reshape block reduction
//  Block '<S336>/Reshapey' : Reshape block reduction
//  Block '<S336>/Reshapeyhat' : Reshape block reduction


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
//  '<S2>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0'
//  '<S3>'   : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1'
//  '<S4>'   : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2'
//  '<S5>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller'
//  '<S6>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/MATLAB Function'
//  '<S7>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator'
//  '<S8>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)'
//  '<S9>'   : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC'
//  '<S10>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S11>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S12>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S13>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S14>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S15>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S16>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S17>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S18>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S19>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S20>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S21>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S22>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S23>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S24>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S25>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S26>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S27>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S28>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S29>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S30>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S31>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S32>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S33>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S34>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S35>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S36>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/moorx'
//  '<S37>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/optimizer'
//  '<S38>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S39>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1'
//  '<S40>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Enable Signal'
//  '<S41>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Covariance'
//  '<S42>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Outputs'
//  '<S43>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Parameters'
//  '<S44>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Initial Regressors'
//  '<S45>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Check Signals'
//  '<S46>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S47>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S48>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialCovariance'
//  '<S49>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialOutputs'
//  '<S50>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialParameters'
//  '<S51>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialRegressors'
//  '<S52>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - L'
//  '<S53>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - Theta'
//  '<S54>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - bufferH'
//  '<S55>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - buffery'
//  '<S56>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/MultiplyWithTranspose'
//  '<S57>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Process Reset'
//  '<S58>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance'
//  '<S59>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialOutputs'
//  '<S60>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialParameters'
//  '<S61>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialRegressors'
//  '<S62>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/RLS'
//  '<S63>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Reset'
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/Process Reset/Check Reset'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance/ScalarExpansion'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/MATLAB Function'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S93>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S94>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S95>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S96>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S97>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S98>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S99>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S100>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S101>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S102>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S103>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S104>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checku'
//  '<S105>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checky'
//  '<S106>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S107>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S108>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S109>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S110>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S111>' : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S112>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller'
//  '<S113>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/MATLAB Function'
//  '<S114>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator'
//  '<S115>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)'
//  '<S116>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC'
//  '<S117>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S118>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S119>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S120>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S121>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S122>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S123>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S124>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S125>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S126>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S127>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S128>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S129>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S130>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S131>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S132>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S133>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S134>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S135>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S136>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S137>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S138>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S139>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S140>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S141>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S142>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S143>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/moorx'
//  '<S144>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer'
//  '<S145>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S146>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/MATLAB Function1'
//  '<S147>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1'
//  '<S148>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2'
//  '<S149>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Enable Signal'
//  '<S150>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Covariance'
//  '<S151>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Outputs'
//  '<S152>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Parameters'
//  '<S153>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Initial Regressors'
//  '<S154>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Check Signals'
//  '<S155>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter1'
//  '<S156>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - AdaptationParameter2'
//  '<S157>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialCovariance'
//  '<S158>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialOutputs'
//  '<S159>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialParameters'
//  '<S160>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - InitialRegressors'
//  '<S161>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - L'
//  '<S162>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - Theta'
//  '<S163>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - bufferH'
//  '<S164>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Data Type Conversion - buffery'
//  '<S165>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/MultiplyWithTranspose'
//  '<S166>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Process Reset'
//  '<S167>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance'
//  '<S168>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialOutputs'
//  '<S169>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialParameters'
//  '<S170>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialRegressors'
//  '<S171>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/RLS'
//  '<S172>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Reset'
//  '<S173>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/Process Reset/Check Reset'
//  '<S174>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator1/ProcessInitialCovariance/ScalarExpansion'
//  '<S175>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Enable Signal'
//  '<S176>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Covariance'
//  '<S177>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Outputs'
//  '<S178>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Parameters'
//  '<S179>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Initial Regressors'
//  '<S180>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Check Signals'
//  '<S181>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - AdaptationParameter1'
//  '<S182>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - AdaptationParameter2'
//  '<S183>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialCovariance'
//  '<S184>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialOutputs'
//  '<S185>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialParameters'
//  '<S186>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - InitialRegressors'
//  '<S187>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - L'
//  '<S188>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - Theta'
//  '<S189>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - bufferH'
//  '<S190>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Data Type Conversion - buffery'
//  '<S191>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/MultiplyWithTranspose'
//  '<S192>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Process Reset'
//  '<S193>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialCovariance'
//  '<S194>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialOutputs'
//  '<S195>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialParameters'
//  '<S196>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialRegressors'
//  '<S197>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/RLS'
//  '<S198>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Reset'
//  '<S199>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/Process Reset/Check Reset'
//  '<S200>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/Recursive Least Squares Estimator2/ProcessInitialCovariance/ScalarExpansion'
//  '<S201>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2'
//  '<S202>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/MATLAB Function'
//  '<S203>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S204>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S205>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S206>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S207>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S208>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S209>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S210>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S211>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S212>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S213>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S214>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S215>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S216>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S217>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S218>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S219>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S220>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S221>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S222>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S223>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S224>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S225>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S226>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S227>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S228>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S229>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S230>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S231>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S232>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S233>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S234>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S235>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S236>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S237>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S238>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S239>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checku'
//  '<S240>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checky'
//  '<S241>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S242>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S243>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S244>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S245>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S246>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S247>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller'
//  '<S248>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/MATLAB Function'
//  '<S249>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator'
//  '<S250>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)'
//  '<S251>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC'
//  '<S252>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S253>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S254>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S255>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S256>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S257>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S258>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S259>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S260>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S261>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S262>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S263>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S264>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S265>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S266>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S267>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S268>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S269>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S270>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S271>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S272>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S273>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S274>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S275>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S276>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S277>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S278>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/moorx'
//  '<S279>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer'
//  '<S280>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S281>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/MATLAB Function2'
//  '<S282>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3'
//  '<S283>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4'
//  '<S284>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Enable Signal'
//  '<S285>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Covariance'
//  '<S286>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Outputs'
//  '<S287>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Parameters'
//  '<S288>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Initial Regressors'
//  '<S289>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Check Signals'
//  '<S290>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - AdaptationParameter1'
//  '<S291>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - AdaptationParameter2'
//  '<S292>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialCovariance'
//  '<S293>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialOutputs'
//  '<S294>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialParameters'
//  '<S295>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - InitialRegressors'
//  '<S296>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - L'
//  '<S297>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - Theta'
//  '<S298>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - bufferH'
//  '<S299>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Data Type Conversion - buffery'
//  '<S300>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/MultiplyWithTranspose'
//  '<S301>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Process Reset'
//  '<S302>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialCovariance'
//  '<S303>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialOutputs'
//  '<S304>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialParameters'
//  '<S305>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialRegressors'
//  '<S306>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/RLS'
//  '<S307>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Reset'
//  '<S308>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/Process Reset/Check Reset'
//  '<S309>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator3/ProcessInitialCovariance/ScalarExpansion'
//  '<S310>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Enable Signal'
//  '<S311>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Covariance'
//  '<S312>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Outputs'
//  '<S313>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Parameters'
//  '<S314>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Initial Regressors'
//  '<S315>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Check Signals'
//  '<S316>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - AdaptationParameter1'
//  '<S317>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - AdaptationParameter2'
//  '<S318>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialCovariance'
//  '<S319>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialOutputs'
//  '<S320>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialParameters'
//  '<S321>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - InitialRegressors'
//  '<S322>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - L'
//  '<S323>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - Theta'
//  '<S324>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - bufferH'
//  '<S325>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Data Type Conversion - buffery'
//  '<S326>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/MultiplyWithTranspose'
//  '<S327>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Process Reset'
//  '<S328>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialCovariance'
//  '<S329>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialOutputs'
//  '<S330>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialParameters'
//  '<S331>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialRegressors'
//  '<S332>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/RLS'
//  '<S333>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Reset'
//  '<S334>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/Process Reset/Check Reset'
//  '<S335>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/Recursive Least Squares Estimator4/ProcessInitialCovariance/ScalarExpansion'
//  '<S336>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2'
//  '<S337>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/MATLAB Function'
//  '<S338>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S339>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S340>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S341>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S342>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S343>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S344>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S345>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S346>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S347>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S348>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S349>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S350>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S351>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S352>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S353>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S354>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S355>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S356>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S357>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S358>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S359>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S360>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S361>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S362>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S363>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S364>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S365>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S366>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S367>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S368>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S369>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S370>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S371>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S372>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S373>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S374>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checku'
//  '<S375>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checky'
//  '<S376>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S377>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S378>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S379>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S380>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S381>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
