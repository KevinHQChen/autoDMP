//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.1475
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May 29 01:50:12 2023
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
  // Block signals and states (default storage) for system '<S145>/MeasurementUpdate' 
  struct DW_MeasurementUpdate {
    boolean_T MeasurementUpdate_MODE;  // '<S145>/MeasurementUpdate'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    DW_MeasurementUpdate MeasurementUpdate_h;// '<S230>/MeasurementUpdate'
    DW_MeasurementUpdate MeasurementUpdate_cg;// '<S145>/MeasurementUpdate'
    real_T Product3[4];                // '<S254>/Product3'
    real_T Product3_p[4];              // '<S169>/Product3'
    real_T Product3_l[2];              // '<S84>/Product3'
    real_T last_mv_DSTATE[3];          // '<S176>/last_mv'
    real_T Delay_DSTATE[3];            // '<S207>/Delay'
    real_T UnitDelay2_DSTATE[3];       // '<S174>/Unit Delay2'
    real_T UnitDelay3_DSTATE[2];       // '<S174>/Unit Delay3'
    real_T Delay1_DSTATE[9];           // '<S207>/Delay1'
    real_T Delay_DSTATE_c[3];          // '<S208>/Delay'
    real_T Delay1_DSTATE_c[9];         // '<S208>/Delay1'
    real_T MemoryX_DSTATE[4];          // '<S211>/MemoryX'
    real_T MemoryP_DSTATE[16];         // '<S211>/MemoryP'
    real_T DiscreteFilter_states[118]; // '<S4>/Discrete Filter'
    real_T last_mv_DSTATE_j[3];        // '<S91>/last_mv'
    real_T Delay_DSTATE_l[3];          // '<S122>/Delay'
    real_T UnitDelay2_DSTATE_m[3];     // '<S89>/Unit Delay2'
    real_T UnitDelay3_DSTATE_g[2];     // '<S89>/Unit Delay3'
    real_T Delay1_DSTATE_cq[9];        // '<S122>/Delay1'
    real_T Delay_DSTATE_m[3];          // '<S123>/Delay'
    real_T Delay1_DSTATE_g[9];         // '<S123>/Delay1'
    real_T MemoryX_DSTATE_m[4];        // '<S126>/MemoryX'
    real_T MemoryP_DSTATE_a[16];       // '<S126>/MemoryP'
    real_T DiscreteFilter_states_b[118];// '<S3>/Discrete Filter'
    real_T last_mv_DSTATE_h[3];        // '<S9>/last_mv'
    real_T Delay_DSTATE_g[2];          // '<S39>/Delay'
    real_T UnitDelay2_DSTATE_j[3];     // '<S7>/Unit Delay2'
    real_T Delay1_DSTATE_h[4];         // '<S39>/Delay1'
    real_T MemoryX_DSTATE_a[2];        // '<S41>/MemoryX'
    real_T MemoryP_DSTATE_d[4];        // '<S41>/MemoryP'
    real_T DiscreteFilter_states_d[59];// '<S2>/Discrete Filter'
    real_T traj[7200];                 // '<Root>/SupervisoryController'
    real_T ymax1[2];                   // '<Root>/SupervisoryController'
    real_T ymax2[2];                   // '<Root>/SupervisoryController'
    real_T uclean[3];                  // '<Root>/SupervisoryController'
    real_T B_0[3];                     // '<Root>/SupervisoryController'
    real_T B_1[6];                     // '<Root>/SupervisoryController'
    real_T B_2[6];                     // '<Root>/SupervisoryController'
    real_T NextOutput[3];              // '<S4>/excitation'
    real_T NextOutput_i[3];            // '<S3>/excitation'
    real_T NextOutput_j[3];            // '<S2>/excitation'
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
    real_T UnitDelay3_DSTATE_m;        // '<S7>/Unit Delay3'
    real_T holdT;                      // '<Root>/SupervisoryController'
    real_T ymax0;                      // '<Root>/SupervisoryController'
    uint32_T RandSeed[3];              // '<S4>/excitation'
    uint32_T RandSeed_o[3];            // '<S3>/excitation'
    uint32_T RandSeed_h[3];            // '<S2>/excitation'
    uint16_T waypt;                    // '<Root>/SupervisoryController'
    uint16_T trajSize;                 // '<Root>/SupervisoryController'
    uint8_T is_c6_SupervisoryController;// '<Root>/SupervisoryController'
    uint8_T is_EventHandler;           // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_n;         // '<Root>/SupervisoryController'
    uint8_T is_EventHandler_k;         // '<Root>/SupervisoryController'
    uint8_T is_active_c6_SupervisoryControl;// '<Root>/SupervisoryController'
    boolean_T Memory_PreviousInput[86];// '<S176>/Memory'
    boolean_T Memory_PreviousInput_b[86];// '<S91>/Memory'
    boolean_T Memory_PreviousInput_c[46];// '<S9>/Memory'
    boolean_T evDone;                  // '<Root>/SupervisoryController'
    boolean_T rstP;                    // '<Root>/SupervisoryController'
    boolean_T icLoad;                  // '<S207>/Delay'
    boolean_T icLoad_m;                // '<S207>/Delay1'
    boolean_T icLoad_a;                // '<S208>/Delay'
    boolean_T icLoad_p;                // '<S208>/Delay1'
    boolean_T icLoad_i;                // '<S211>/MemoryX'
    boolean_T icLoad_j;                // '<S211>/MemoryP'
    boolean_T icLoad_i5;               // '<S122>/Delay'
    boolean_T icLoad_f;                // '<S122>/Delay1'
    boolean_T icLoad_c;                // '<S123>/Delay'
    boolean_T icLoad_pp;               // '<S123>/Delay1'
    boolean_T icLoad_d;                // '<S126>/MemoryX'
    boolean_T icLoad_jj;               // '<S126>/MemoryP'
    boolean_T icLoad_b;                // '<S39>/Delay'
    boolean_T icLoad_fo;               // '<S39>/Delay1'
    boolean_T icLoad_k;                // '<S41>/MemoryX'
    boolean_T icLoad_fq;               // '<S41>/MemoryP'
    boolean_T MeasurementUpdate_MODE;  // '<S60>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state
  struct PrevZCX {
    ZCSigState SupervisoryController_Trig_ZCE;// '<Root>/SupervisoryController'
    ZCSigState Delay1_Reset_ZCE;       // '<S207>/Delay1'
    ZCSigState Delay1_Reset_ZCE_g;     // '<S208>/Delay1'
    ZCSigState MemoryX_Reset_ZCE;      // '<S211>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S211>/MemoryP'
    ZCSigState Delay1_Reset_ZCE_d;     // '<S122>/Delay1'
    ZCSigState Delay1_Reset_ZCE_o;     // '<S123>/Delay1'
    ZCSigState MemoryX_Reset_ZCE_h;    // '<S126>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_c;    // '<S126>/MemoryP'
    ZCSigState MemoryX_Reset_ZCE_j;    // '<S41>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE_a;    // '<S41>/MemoryP'
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
    real_T dPmod_;                     // '<Root>/dPmod_'
    real_T p_;                         // '<Root>/p_'
    real_T measAvail;                  // '<Root>/measAvail'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY {
    real_T u[3];                       // '<Root>/u'
    event_bus currEv;                  // '<Root>/currEv'
    boolean_T requestEvent;            // '<Root>/requestEvent'
    real_T currTraj[3];                // '<Root>/currTraj'
    real_T yhat[3];                    // '<Root>/yhat'
    real_T B_a[9];                     // '<Root>/B'
    real_T uref[3];                    // '<Root>/uref'
    real_T uoffset[3];                 // '<Root>/uoffset'
    real_T paramEstErr[3];             // '<Root>/paramEstErr'
  };

  // Parameters for system: '<S145>/MeasurementUpdate'
  struct P_MeasurementUpdate {
    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S169>/L*(y[k]-yhat[k|k-1])'

  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod0;                       // Variable: Aod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Aod1[4];                    // Variable: Aod1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S175>/MATLAB Function'

    real_T Bod0;                       // Variable: Bod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Bod1[4];                    // Variable: Bod1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S175>/MATLAB Function'

    real_T Cod0;                       // Variable: Cod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Cod1[4];                    // Variable: Cod1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S175>/MATLAB Function'

    real_T Dmn0;                       // Variable: Dmn0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dmn1[4];                    // Variable: Dmn1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S175>/MATLAB Function'

    real_T Dod0;                       // Variable: Dod0
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dod1[4];                    // Variable: Dod1
                                          //  Referenced by:
                                          //    '<S90>/MATLAB Function'
                                          //    '<S175>/MATLAB Function'

    real_T dt;                         // Variable: dt
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T forgettingFactor;           // Variable: forgettingFactor
                                          //  Referenced by:
                                          //    '<S39>/addLambda'
                                          //    '<S39>/divByLambda'
                                          //    '<S122>/addLambda'
                                          //    '<S122>/divByLambda'
                                          //    '<S123>/addLambda'
                                          //    '<S123>/divByLambda'
                                          //    '<S207>/addLambda'
                                          //    '<S207>/divByLambda'
                                          //    '<S208>/addLambda'
                                          //    '<S208>/divByLambda'

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
                                          //  Referenced by: '<S84>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

    real_T prms_Y0;                    // Computed Parameter: prms_Y0
                                          //  Referenced by: '<S2>/prms'

    real_T uref_Y0;                    // Computed Parameter: uref_Y0
                                          //  Referenced by: '<S2>/uref'

    real_T prmErr_Y0;                  // Computed Parameter: prmErr_Y0
                                          //  Referenced by: '<S2>/prmErr'

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

    real_T Constant4_Value[4];         // Expression: 1e4*eye(2)
                                          //  Referenced by: '<S7>/Constant4'

    real_T Constant13_Value[3];        // Expression: G0.D
                                          //  Referenced by: '<S7>/Constant13'

    real_T Constant3_Value;            // Expression: G0.A
                                          //  Referenced by: '<S7>/Constant3'

    real_T Constant12_Value;           // Expression: G0.C
                                          //  Referenced by: '<S7>/Constant12'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S7>/Constant11'

    real_T Constant2_Value;            // Expression: 0
                                          //  Referenced by: '<S7>/Constant2'

    real_T Constant_Value;             // Expression: 0
                                          //  Referenced by: '<S8>/Constant'

    real_T X0_Value[2];                // Expression: pInitialization.X0
                                          //  Referenced by: '<S41>/X0'

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
                                          //  Referenced by: '<S41>/H'

    real_T G_Value[4];                 // Expression: pInitialization.G
                                          //  Referenced by: '<S41>/G'

    real_T Constant_Value_b;           // Expression: 1
                                          //  Referenced by: '<S2>/Constant'

    real_T P0_Value[4];                // Expression: pInitialization.P0
                                          //  Referenced by: '<S41>/P0'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S9>/u_scale'

    real_T excitation_Mean[3];         // Expression: [0 0 0]
                                          //  Referenced by: '<S2>/excitation'

    real_T excitation_StdDev[3];       // Computed Parameter: excitation_StdDev
                                          //  Referenced by: '<S2>/excitation'

    real_T excitation_Seed[3];         // Expression: [12345 12345 12345]
                                          //  Referenced by: '<S2>/excitation'

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

    real_T prms_Y0_l;                  // Computed Parameter: prms_Y0_l
                                          //  Referenced by: '<S3>/prms'

    real_T uref_Y0_k;                  // Computed Parameter: uref_Y0_k
                                          //  Referenced by: '<S3>/uref'

    real_T prmErr_Y0_k;                // Computed Parameter: prmErr_Y0_k
                                          //  Referenced by: '<S3>/prmErr'

    real_T uclean_Y0_e;                // Computed Parameter: uclean_Y0_e
                                          //  Referenced by: '<S3>/uclean'

    real_T G_zero_Value_m;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S87>/G_zero'

    real_T LastPcov_InitialCondition_a[16];// Expression: lastPcov
                                              //  Referenced by: '<S91>/LastPcov'

    real_T ywt_zero_Value_l[2];        // Expression: zeros(2,1)
                                          //  Referenced by: '<S87>/y.wt_zero'

    real_T uwt_zero_Value_o[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/u.wt_zero'

    real_T duwt_zero_Value_n[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/du.wt_zero'

    real_T extmv_zero_Value_k[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/ext.mv_zero'

    real_T extmv_scale_Gain_g[3];      // Expression: RMVscale
                                          //  Referenced by: '<S91>/ext.mv_scale'

    real_T last_mv_InitialCondition_o[3];// Expression: lastu+uoff
                                            //  Referenced by: '<S91>/last_mv'

    real_T Constant3_Value_j[4];       // Expression: G1.A
                                          //  Referenced by: '<S89>/Constant3'

    real_T UnitDelay2_InitialCondition_f;// Expression: 0
                                            //  Referenced by: '<S89>/Unit Delay2'

    real_T UnitDelay3_InitialCondition_g;// Expression: 0
                                            //  Referenced by: '<S89>/Unit Delay3'

    real_T Constant4_Value_k[9];       // Expression: 1e4*eye(3)
                                          //  Referenced by: '<S89>/Constant4'

    real_T Constant5_Value[9];         // Expression: 1e4*eye(3)
                                          //  Referenced by: '<S89>/Constant5'

    real_T Constant12_Value_a[4];      // Expression: G1.C
                                          //  Referenced by: '<S89>/Constant12'

    real_T Constant13_Value_m[6];      // Expression: G1.D
                                          //  Referenced by: '<S89>/Constant13'

    real_T Constant11_Value_c;         // Expression: 1
                                          //  Referenced by: '<S89>/Constant11'

    real_T Constant2_Value_c[2];       // Expression: [0;0]
                                          //  Referenced by: '<S89>/Constant2'

    real_T Constant_Value_l[2];        // Expression: [0;0]
                                          //  Referenced by: '<S90>/Constant'

    real_T X0_Value_m[4];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S126>/X0'

    real_T ym_zero_Value_a[2];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S91>/ym_zero'

    real_T md_zero_Value_d;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S87>/md_zero'

    real_T umin_zero_Value_i[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/umin_zero'

    real_T umax_zero_Value_b[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/umax_zero'

    real_T ymin_zero_Value_c[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S87>/ymin_zero'

    real_T ymax_zero_Value_g[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S87>/ymax_zero'

    real_T E_zero_Value_c[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S87>/E_zero'

    real_T umin_scale4_Gain_b[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S91>/umin_scale4'

    real_T F_zero_Value_k[2];          // Expression: zeros(1,2)
                                          //  Referenced by: '<S87>/F_zero'

    real_T ymin_scale1_Gain_h[2];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S91>/ymin_scale1'

    real_T S_zero_Value_b;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S87>/S_zero'

    real_T ymin_scale2_Gain_g;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S91>/ymin_scale2'

    real_T switch_zero_Value_g;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S87>/switch_zero'

    real_T mvtarget_zero_Value_c[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S87>/mv.target_zero'

    real_T uref_scale_Gain_m[3];       // Expression: RMVscale
                                          //  Referenced by: '<S91>/uref_scale'

    real_T ecrwt_zero_Value_m;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S87>/ecr.wt_zero'

    real_T P0_Value_a[16];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S126>/P0'

    real_T Constant_Value_n;           // Expression: 1
                                          //  Referenced by: '<S3>/Constant'

    real_T H_Value_o[8];               // Expression: pInitialization.H
                                          //  Referenced by: '<S126>/H'

    real_T G_Value_a[16];              // Expression: pInitialization.G
                                          //  Referenced by: '<S126>/G'

    real_T u_scale_Gain_g[3];          // Expression: MVscale
                                          //  Referenced by: '<S91>/u_scale'

    real_T Constant2_Value_j;          // Expression: 0
                                          //  Referenced by: '<S3>/Constant2'

    real_T excitation_Mean_o[3];       // Expression: [0 0 0]
                                          //  Referenced by: '<S3>/excitation'

    real_T excitation_StdDev_k[3];    // Computed Parameter: excitation_StdDev_k
                                         //  Referenced by: '<S3>/excitation'

    real_T excitation_Seed_b[3];       // Expression: [12345 12346 12347]
                                          //  Referenced by: '<S3>/excitation'

    real_T DiscreteFilter_InitialStates_d;// Expression: 0
                                             //  Referenced by: '<S3>/Discrete Filter'

    real_T Constant3_Value_k;          // Expression: 0
                                          //  Referenced by: '<S3>/Constant3'

    real_T Switch_Threshold_j;         // Expression: 0
                                          //  Referenced by: '<S3>/Switch'

    real_T u_Y0_h;                     // Computed Parameter: u_Y0_h
                                          //  Referenced by: '<S4>/u'

    real_T yhat_Y0_i;                  // Computed Parameter: yhat_Y0_i
                                          //  Referenced by: '<S4>/yhat'

    real_T prms_Y0_e;                  // Computed Parameter: prms_Y0_e
                                          //  Referenced by: '<S4>/prms'

    real_T uref_Y0_km;                 // Computed Parameter: uref_Y0_km
                                          //  Referenced by: '<S4>/uref'

    real_T prmErr_Y0_c;                // Computed Parameter: prmErr_Y0_c
                                          //  Referenced by: '<S4>/prmErr'

    real_T uclean_Y0_k;                // Computed Parameter: uclean_Y0_k
                                          //  Referenced by: '<S4>/uclean'

    real_T G_zero_Value_j;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S172>/G_zero'

    real_T LastPcov_InitialCondition_k[16];// Expression: lastPcov
                                              //  Referenced by: '<S176>/LastPcov'

    real_T ywt_zero_Value_h[2];        // Expression: zeros(2,1)
                                          //  Referenced by: '<S172>/y.wt_zero'

    real_T uwt_zero_Value_d[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/u.wt_zero'

    real_T duwt_zero_Value_l[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/du.wt_zero'

    real_T extmv_zero_Value_f[3];      // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/ext.mv_zero'

    real_T extmv_scale_Gain_gs[3];     // Expression: RMVscale
                                          //  Referenced by: '<S176>/ext.mv_scale'

    real_T last_mv_InitialCondition_oo[3];// Expression: lastu+uoff
                                             //  Referenced by: '<S176>/last_mv'

    real_T Constant3_Value_h[4];       // Expression: G2.A
                                          //  Referenced by: '<S174>/Constant3'

    real_T UnitDelay2_InitialCondition_o;// Expression: 0
                                            //  Referenced by: '<S174>/Unit Delay2'

    real_T UnitDelay3_InitialCondition_m;// Expression: 0
                                            //  Referenced by: '<S174>/Unit Delay3'

    real_T Constant4_Value_j[9];       // Expression: 1e4*eye(3)
                                          //  Referenced by: '<S174>/Constant4'

    real_T Constant5_Value_l[9];       // Expression: 1e4*eye(3)
                                          //  Referenced by: '<S174>/Constant5'

    real_T Constant12_Value_l[4];      // Expression: G2.C
                                          //  Referenced by: '<S174>/Constant12'

    real_T Constant13_Value_k[6];      // Expression: G2.D
                                          //  Referenced by: '<S174>/Constant13'

    real_T Constant11_Value_o;         // Expression: 1
                                          //  Referenced by: '<S174>/Constant11'

    real_T Constant2_Value_d[2];       // Expression: [0;0]
                                          //  Referenced by: '<S174>/Constant2'

    real_T Constant_Value_c[2];        // Expression: [0;0]
                                          //  Referenced by: '<S175>/Constant'

    real_T X0_Value_g[4];              // Expression: pInitialization.X0
                                          //  Referenced by: '<S211>/X0'

    real_T ym_zero_Value_m[2];         // Expression: zeros(nym,1)
                                          //  Referenced by: '<S176>/ym_zero'

    real_T md_zero_Value_f;            // Expression: zeros(1,1)
                                          //  Referenced by: '<S172>/md_zero'

    real_T umin_zero_Value_b[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/umin_zero'

    real_T umax_zero_Value_a[3];       // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/umax_zero'

    real_T ymin_zero_Value_p[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S172>/ymin_zero'

    real_T ymax_zero_Value_j[2];       // Expression: zeros(2,1)
                                          //  Referenced by: '<S172>/ymax_zero'

    real_T E_zero_Value_j[3];          // Expression: zeros(1,3)
                                          //  Referenced by: '<S172>/E_zero'

    real_T umin_scale4_Gain_a[3];  // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S176>/umin_scale4'

    real_T F_zero_Value_kj[2];         // Expression: zeros(1,2)
                                          //  Referenced by: '<S172>/F_zero'

    real_T ymin_scale1_Gain_c[2];   // Expression: Yscale(:,ones(1,max(nCC,1)))'
                                       //  Referenced by: '<S176>/ymin_scale1'

    real_T S_zero_Value_h;             // Expression: zeros(1,1)
                                          //  Referenced by: '<S172>/S_zero'

    real_T ymin_scale2_Gain_e;     // Expression: MDscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S176>/ymin_scale2'

    real_T switch_zero_Value_b;        // Expression: zeros(1,1)
                                          //  Referenced by: '<S172>/switch_zero'

    real_T mvtarget_zero_Value_g[3];   // Expression: zeros(3,1)
                                          //  Referenced by: '<S172>/mv.target_zero'

    real_T uref_scale_Gain_p[3];       // Expression: RMVscale
                                          //  Referenced by: '<S176>/uref_scale'

    real_T ecrwt_zero_Value_j;         // Expression: zeros(1,1)
                                          //  Referenced by: '<S172>/ecr.wt_zero'

    real_T P0_Value_m[16];             // Expression: pInitialization.P0
                                          //  Referenced by: '<S211>/P0'

    real_T Constant_Value_i;           // Expression: 1
                                          //  Referenced by: '<S4>/Constant'

    real_T H_Value_c[8];               // Expression: pInitialization.H
                                          //  Referenced by: '<S211>/H'

    real_T G_Value_d[16];              // Expression: pInitialization.G
                                          //  Referenced by: '<S211>/G'

    real_T u_scale_Gain_p[3];          // Expression: MVscale
                                          //  Referenced by: '<S176>/u_scale'

    real_T Constant2_Value_l;          // Expression: 0
                                          //  Referenced by: '<S4>/Constant2'

    real_T excitation_Mean_k[3];       // Expression: [0 0 0]
                                          //  Referenced by: '<S4>/excitation'

    real_T excitation_StdDev_h[3];    // Computed Parameter: excitation_StdDev_h
                                         //  Referenced by: '<S4>/excitation'

    real_T excitation_Seed_n[3];       // Expression: [12345 12346 12347]
                                          //  Referenced by: '<S4>/excitation'

    real_T DiscreteFilter_InitialStates_h;// Expression: 0
                                             //  Referenced by: '<S4>/Discrete Filter'

    real_T Constant3_Value_i;          // Expression: 0
                                          //  Referenced by: '<S4>/Constant3'

    real_T Switch_Threshold_e;         // Expression: 0
                                          //  Referenced by: '<S4>/Switch'

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S37>/FixedHorizonOptimizer'

    int32_T FixedHorizonOptimizer_Ndis_b;// Expression: Ndis
                                            //  Referenced by: '<S119>/FixedHorizonOptimizer'

    int32_T FixedHorizonOptimizer_Ndis_h;// Expression: Ndis
                                            //  Referenced by: '<S204>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[46];// Expression: iA
                                             //  Referenced by: '<S9>/Memory'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S82>/isSqrtUsed'

    boolean_T Memory_InitialCondition_p[86];// Expression: iA
                                               //  Referenced by: '<S91>/Memory'

    boolean_T isSqrtUsed_Value_l;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S167>/isSqrtUsed'

    boolean_T Constant1_Value_h;       // Expression: true
                                          //  Referenced by: '<S3>/Constant1'

    boolean_T Memory_InitialCondition_j[86];// Expression: iA
                                               //  Referenced by: '<S176>/Memory'

    boolean_T isSqrtUsed_Value_d;      // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S252>/isSqrtUsed'

    boolean_T Constant1_Value_m;       // Expression: true
                                          //  Referenced by: '<S4>/Constant1'

    P_MeasurementUpdate MeasurementUpdate_h;// '<S230>/MeasurementUpdate'
    P_MeasurementUpdate MeasurementUpdate_cg;// '<S145>/MeasurementUpdate'
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

  // private member function(s) for subsystem '<S89>/MATLAB Function1'
  static void MATLABFunction1(const real_T rtu_B1[3], const real_T rtu_B2[3],
    real_T rty_B[6]);

  // private member function(s) for subsystem '<S122>/MATLAB Function'
  static void MATLABFunction_n(const real_T rtu_theta[3], const real_T rtu_P[9],
    real_T rtu_epsil, const real_T rtu_phi[3], real_T rtu_ms_, boolean_T rtu_EN,
    real_T rtu_p_, real_T rtu_dPmod_, real_T rty_dtheta[3], real_T rty_dP[9]);

  // private member function(s) for subsystem '<S126>/CalculatePL'
  void CalculatePL(const real_T rtu_Ak[16], const real_T rtu_Ck[8], const real_T
                   rtu_Qbark[16], const real_T rtu_Rbark[4], const real_T
                   rtu_Nbark[8], boolean_T rtu_Enablek, const real_T rtu_Pk[16],
                   real_T rty_Mk[8], real_T rty_Lk[8], real_T rty_Zk[16], real_T
                   rty_Pk1[16]);
  void mrdiv(const real_T A[8], const real_T B_0[4], real_T Y[8]);

  // private member function(s) for subsystem '<S167>/SqrtUsedFcn'
  static void SqrtUsedFcn(const real_T rtu_u[16], boolean_T rtu_isSqrtUsed,
    real_T rty_P[16]);

  // private member function(s) for subsystem '<S145>/MeasurementUpdate'
  static void MeasurementUpdate_Init(real_T rty_Lykyhatkk1[4],
    P_MeasurementUpdate *localP);
  static void MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);
  void MeasurementUpdate(boolean_T rtu_Enable, const real_T rtu_Lk[8], const
    real_T rtu_yk[2], const real_T rtu_Ck[8], const real_T rtu_xhatkk1[4], const
    real_T rtu_Dk[6], const real_T rtu_uk[3], real_T rty_Lykyhatkk1[4],
    DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP);

  // private member function(s) for subsystem '<S126>/ReducedQRN'
  static void ReducedQRN(const real_T rtu_G[16], const real_T rtu_H[8], const
    real_T rtu_Q[16], const real_T rtu_R[4], const real_T rtu_N[8], real_T
    rty_Qbar[16], real_T rty_Rbar[4], real_T rty_Nbar[8]);

  // private member function(s) for subsystem '<S126>/ScalarExpansionQ'
  static void ScalarExpansionQ(const real_T rtu_u[16], real_T rty_y[16]);

  // private member function(s) for subsystem '<S126>/ScalarExpansionR'
  static void ScalarExpansionR(const real_T rtu_u[4], real_T rty_y[4]);

  // private member function(s) for subsystem '<S90>/MATLAB Function'
  static void MATLABFunction_c(const real_T rtu_Ap[4], const real_T rtu_Bp[6],
    const real_T rtu_Cp[4], const real_T rtu_Dp[6], real_T rty_A[16], real_T
    rty_B[12], real_T rty_C[8], real_T rty_D[6], real_T rty_Q[16], real_T rty_R
    [4], real_T rty_N[8], P *rtP);

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
  void trisolve(const real_T b_A[16], real_T b_B[16]);
  real_T norm(const real_T x[4]);
  real_T maximum(const real_T x[4]);
  real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
             real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
             [4], real_T b_A[16], int32_T ia0);
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
//  Block '<S7>/P' : Unused code path elimination
//  Block '<S50>/Data Type Duplicate' : Unused code path elimination
//  Block '<S51>/Data Type Duplicate' : Unused code path elimination
//  Block '<S53>/Data Type Duplicate' : Unused code path elimination
//  Block '<S54>/Data Type Duplicate' : Unused code path elimination
//  Block '<S57>/Data Type Duplicate' : Unused code path elimination
//  Block '<S58>/Data Type Duplicate' : Unused code path elimination
//  Block '<S66>/CheckSignalProperties' : Unused code path elimination
//  Block '<S67>/CheckSignalProperties' : Unused code path elimination
//  Block '<S68>/CheckSignalProperties' : Unused code path elimination
//  Block '<S69>/CheckSignalProperties' : Unused code path elimination
//  Block '<S70>/CheckSignalProperties' : Unused code path elimination
//  Block '<S73>/CheckSignalProperties' : Unused code path elimination
//  Block '<S75>/CheckSignalProperties' : Unused code path elimination
//  Block '<S76>/CheckSignalProperties' : Unused code path elimination
//  Block '<S77>/CheckSignalProperties' : Unused code path elimination
//  Block '<S79>/CheckSignalProperties' : Unused code path elimination
//  Block '<S80>/CheckSignalProperties' : Unused code path elimination
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
//  Block '<S104>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S105>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S106>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S107>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S108>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S109>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S110>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S111>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S112>/Vector Dimension Check' : Unused code path elimination
//  Block '<S113>/Vector Dimension Check' : Unused code path elimination
//  Block '<S114>/Vector Dimension Check' : Unused code path elimination
//  Block '<S115>/Vector Dimension Check' : Unused code path elimination
//  Block '<S116>/Vector Dimension Check' : Unused code path elimination
//  Block '<S117>/Vector Dimension Check' : Unused code path elimination
//  Block '<S91>/last_x' : Unused code path elimination
//  Block '<S118>/Vector Dimension Check' : Unused code path elimination
//  Block '<S91>/useq_scale' : Unused code path elimination
//  Block '<S91>/useq_scale1' : Unused code path elimination
//  Block '<S87>/m_zero' : Unused code path elimination
//  Block '<S87>/p_zero' : Unused code path elimination
//  Block '<S89>/P' : Unused code path elimination
//  Block '<S135>/Data Type Duplicate' : Unused code path elimination
//  Block '<S136>/Data Type Duplicate' : Unused code path elimination
//  Block '<S138>/Data Type Duplicate' : Unused code path elimination
//  Block '<S139>/Data Type Duplicate' : Unused code path elimination
//  Block '<S142>/Data Type Duplicate' : Unused code path elimination
//  Block '<S143>/Data Type Duplicate' : Unused code path elimination
//  Block '<S151>/CheckSignalProperties' : Unused code path elimination
//  Block '<S152>/CheckSignalProperties' : Unused code path elimination
//  Block '<S153>/CheckSignalProperties' : Unused code path elimination
//  Block '<S154>/CheckSignalProperties' : Unused code path elimination
//  Block '<S155>/CheckSignalProperties' : Unused code path elimination
//  Block '<S158>/CheckSignalProperties' : Unused code path elimination
//  Block '<S160>/CheckSignalProperties' : Unused code path elimination
//  Block '<S161>/CheckSignalProperties' : Unused code path elimination
//  Block '<S162>/CheckSignalProperties' : Unused code path elimination
//  Block '<S164>/CheckSignalProperties' : Unused code path elimination
//  Block '<S165>/CheckSignalProperties' : Unused code path elimination
//  Block '<S176>/Floor' : Unused code path elimination
//  Block '<S176>/Floor1' : Unused code path elimination
//  Block '<S177>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S178>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S179>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S180>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S181>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S182>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S183>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S184>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S185>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S186>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S187>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S188>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S189>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S190>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S191>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S192>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S193>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S194>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S195>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S196>/Matrix Dimension Check' : Unused code path elimination
//  Block '<S197>/Vector Dimension Check' : Unused code path elimination
//  Block '<S198>/Vector Dimension Check' : Unused code path elimination
//  Block '<S199>/Vector Dimension Check' : Unused code path elimination
//  Block '<S200>/Vector Dimension Check' : Unused code path elimination
//  Block '<S201>/Vector Dimension Check' : Unused code path elimination
//  Block '<S202>/Vector Dimension Check' : Unused code path elimination
//  Block '<S176>/last_x' : Unused code path elimination
//  Block '<S203>/Vector Dimension Check' : Unused code path elimination
//  Block '<S176>/useq_scale' : Unused code path elimination
//  Block '<S176>/useq_scale1' : Unused code path elimination
//  Block '<S172>/m_zero' : Unused code path elimination
//  Block '<S172>/p_zero' : Unused code path elimination
//  Block '<S174>/P' : Unused code path elimination
//  Block '<S220>/Data Type Duplicate' : Unused code path elimination
//  Block '<S221>/Data Type Duplicate' : Unused code path elimination
//  Block '<S223>/Data Type Duplicate' : Unused code path elimination
//  Block '<S224>/Data Type Duplicate' : Unused code path elimination
//  Block '<S227>/Data Type Duplicate' : Unused code path elimination
//  Block '<S228>/Data Type Duplicate' : Unused code path elimination
//  Block '<S236>/CheckSignalProperties' : Unused code path elimination
//  Block '<S237>/CheckSignalProperties' : Unused code path elimination
//  Block '<S238>/CheckSignalProperties' : Unused code path elimination
//  Block '<S239>/CheckSignalProperties' : Unused code path elimination
//  Block '<S240>/CheckSignalProperties' : Unused code path elimination
//  Block '<S243>/CheckSignalProperties' : Unused code path elimination
//  Block '<S245>/CheckSignalProperties' : Unused code path elimination
//  Block '<S246>/CheckSignalProperties' : Unused code path elimination
//  Block '<S247>/CheckSignalProperties' : Unused code path elimination
//  Block '<S249>/CheckSignalProperties' : Unused code path elimination
//  Block '<S250>/CheckSignalProperties' : Unused code path elimination
//  Block '<S9>/Reshape' : Reshape block reduction
//  Block '<S9>/Reshape1' : Reshape block reduction
//  Block '<S9>/Reshape2' : Reshape block reduction
//  Block '<S9>/Reshape3' : Reshape block reduction
//  Block '<S9>/Reshape4' : Reshape block reduction
//  Block '<S9>/Reshape5' : Reshape block reduction
//  Block '<S39>/Reshape1' : Reshape block reduction
//  Block '<S53>/Conversion' : Eliminate redundant data type conversion
//  Block '<S57>/Conversion' : Eliminate redundant data type conversion
//  Block '<S60>/Reshape' : Reshape block reduction
//  Block '<S41>/ReshapeX0' : Reshape block reduction
//  Block '<S41>/Reshapeu' : Reshape block reduction
//  Block '<S41>/Reshapexhat' : Reshape block reduction
//  Block '<S41>/Reshapey' : Reshape block reduction
//  Block '<S41>/Reshapeyhat' : Reshape block reduction
//  Block '<S91>/Reshape' : Reshape block reduction
//  Block '<S91>/Reshape1' : Reshape block reduction
//  Block '<S91>/Reshape2' : Reshape block reduction
//  Block '<S91>/Reshape3' : Reshape block reduction
//  Block '<S91>/Reshape4' : Reshape block reduction
//  Block '<S91>/Reshape5' : Reshape block reduction
//  Block '<S122>/Reshape1' : Reshape block reduction
//  Block '<S123>/Reshape1' : Reshape block reduction
//  Block '<S89>/Reshape' : Reshape block reduction
//  Block '<S138>/Conversion' : Eliminate redundant data type conversion
//  Block '<S142>/Conversion' : Eliminate redundant data type conversion
//  Block '<S145>/Reshape' : Reshape block reduction
//  Block '<S126>/ReshapeX0' : Reshape block reduction
//  Block '<S126>/Reshapeu' : Reshape block reduction
//  Block '<S126>/Reshapexhat' : Reshape block reduction
//  Block '<S126>/Reshapey' : Reshape block reduction
//  Block '<S126>/Reshapeyhat' : Reshape block reduction
//  Block '<S176>/Reshape' : Reshape block reduction
//  Block '<S176>/Reshape1' : Reshape block reduction
//  Block '<S176>/Reshape2' : Reshape block reduction
//  Block '<S176>/Reshape3' : Reshape block reduction
//  Block '<S176>/Reshape4' : Reshape block reduction
//  Block '<S176>/Reshape5' : Reshape block reduction
//  Block '<S207>/Reshape1' : Reshape block reduction
//  Block '<S208>/Reshape1' : Reshape block reduction
//  Block '<S174>/Reshape' : Reshape block reduction
//  Block '<S223>/Conversion' : Eliminate redundant data type conversion
//  Block '<S227>/Conversion' : Eliminate redundant data type conversion
//  Block '<S230>/Reshape' : Reshape block reduction
//  Block '<S211>/ReshapeX0' : Reshape block reduction
//  Block '<S211>/Reshapeu' : Reshape block reduction
//  Block '<S211>/Reshapexhat' : Reshape block reduction
//  Block '<S211>/Reshapey' : Reshape block reduction
//  Block '<S211>/Reshapeyhat' : Reshape block reduction


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
//  '<S39>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/RLS'
//  '<S40>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/Model Estimator/RLS/MATLAB Function'
//  '<S41>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2'
//  '<S42>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/MATLAB Function'
//  '<S43>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S44>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S45>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S46>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S47>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S48>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S49>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S50>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S51>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S52>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S53>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S54>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S55>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S56>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S57>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S58>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S59>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S60>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S61>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S62>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S63>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checku'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/checky'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/State0.ControlLaw.AMPC0/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/MATLAB Function'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S93>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S94>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S95>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S96>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S97>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S98>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S99>'  : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S100>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S101>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S102>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S103>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S104>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S105>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S106>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S107>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S108>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S109>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S110>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S111>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S112>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S113>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S114>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S115>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S116>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S117>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S118>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/moorx'
//  '<S119>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer'
//  '<S120>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S121>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/MATLAB Function1'
//  '<S122>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/RLS1'
//  '<S123>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/RLS2'
//  '<S124>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/RLS1/MATLAB Function'
//  '<S125>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/Model Estimator/RLS2/MATLAB Function'
//  '<S126>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2'
//  '<S127>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/MATLAB Function'
//  '<S128>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S129>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S130>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S131>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S132>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S133>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S134>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S135>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S136>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S137>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S138>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S139>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S140>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S141>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S142>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S143>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S144>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S145>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S146>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S147>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S148>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S149>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S150>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S151>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S152>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S153>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S154>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S155>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S156>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S157>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S158>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S159>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S160>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S161>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S162>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S163>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S164>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checku'
//  '<S165>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/checky'
//  '<S166>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S167>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S168>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S169>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S170>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S171>' : 'T_junction_mpc/SupervisoryController/State1.ControlLaw.AMPC1/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S172>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller'
//  '<S173>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/MATLAB Function'
//  '<S174>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator'
//  '<S175>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)'
//  '<S176>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC'
//  '<S177>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S178>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S179>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S180>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S181>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S182>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S183>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S184>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S185>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S186>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S187>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S188>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S189>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S190>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S191>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S192>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S193>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S194>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S195>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S196>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S197>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S198>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S199>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S200>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S201>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S202>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S203>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/moorx'
//  '<S204>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer'
//  '<S205>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S206>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/MATLAB Function1'
//  '<S207>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/RLS1'
//  '<S208>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/RLS2'
//  '<S209>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/RLS1/MATLAB Function'
//  '<S210>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/Model Estimator/RLS2/MATLAB Function'
//  '<S211>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2'
//  '<S212>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/MATLAB Function'
//  '<S213>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL'
//  '<S214>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculateYhat'
//  '<S215>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S216>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S217>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S218>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S219>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S220>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S221>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S222>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S223>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S224>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S225>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S226>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S227>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S228>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S229>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S230>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer'
//  '<S231>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ReducedQRN'
//  '<S232>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S233>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S234>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S235>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S236>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkA'
//  '<S237>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkB'
//  '<S238>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkC'
//  '<S239>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkD'
//  '<S240>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkEnable'
//  '<S241>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkG'
//  '<S242>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkH'
//  '<S243>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkN'
//  '<S244>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkP0'
//  '<S245>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkQ'
//  '<S246>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkR'
//  '<S247>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkReset'
//  '<S248>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checkX0'
//  '<S249>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checku'
//  '<S250>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/checky'
//  '<S251>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S252>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S253>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S254>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S255>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S256>' : 'T_junction_mpc/SupervisoryController/State2.ControlLaw.AMPC2/State Estimator (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
