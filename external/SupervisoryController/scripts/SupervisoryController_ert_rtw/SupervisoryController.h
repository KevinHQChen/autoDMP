//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.h
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2224
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Thu Jul  6 16:16:39 2023
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
#define SS_UINT64                      23
#endif

#ifndef SS_INT64
#define SS_INT64                       24
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
    real_T Delay1_DSTATE[144];         // '<S87>/Delay1'
    real_T UnitDelay3_DSTATE[3];       // '<S85>/Unit Delay3'
    real_T Delay_DSTATE[12];           // '<S87>/Delay'
    boolean_T icLoad;                  // '<S87>/Delay1'
    boolean_T icLoad_n;                // '<S87>/Delay'
  };

  // Zero-crossing (trigger) state for system '<S1>/paramEst1'
  struct ZCE_paramEst1 {
    ZCSigState Delay1_Reset_ZCE;       // '<S87>/Delay1'
    ZCSigState Delay_Reset_ZCE;        // '<S87>/Delay'
  };

  // Block signals and states (default storage) for system '<Root>'
  struct DW {
    DW_paramEst1 paramEst2;            // '<S1>/paramEst2'
    DW_paramEst1 paramEst1_o;          // '<S1>/paramEst1'
    real_T Product3[18];               // '<S82>/Product3'
    real_T Delay_DSTATE[6];            // '<S2>/Delay'
    real_T last_mv_DSTATE[3];          // '<S9>/last_mv'
    real_T MemoryX_DSTATE[18];         // '<S39>/MemoryX'
    real_T MemoryP_DSTATE[324];        // '<S39>/MemoryP'
    real_T DiscreteFilter1_states[177];// '<S2>/Discrete Filter1'
    real_T traj[14400];                // '<Root>/SupervisoryController'
    real_T P0_1[144];                  // '<Root>/SupervisoryController'
    real_T P0_2[144];                  // '<Root>/SupervisoryController'
    real_T theta0_1[12];               // '<Root>/SupervisoryController'
    real_T theta0_2[12];               // '<Root>/SupervisoryController'
    real_T thetaSgn[24];               // '<Root>/SupervisoryController'
    real_T NextOutput[3];              // '<S2>/excitation'
    real_T Bu[1134];
    real_T dv[5166];
    real_T f[4428];
    real_T g[738];
    real_T k[984];
    real_T b_Sx[2160];
    real_T b_Hv[2520];
    real_T Su[7200];
    real_T e_data[1199];
    real_T t_data[1199];
    real_T tmp_data[1199];
    real_T SFunction_o11;              // '<Root>/SupervisoryController'
    real_T SFunction_o12;              // '<Root>/SupervisoryController'
    real_T SFunction_o13;              // '<Root>/SupervisoryController'
    real_T SFunction_o14;              // '<Root>/SupervisoryController'
    real_T SFunction_o15;              // '<Root>/SupervisoryController'
    real_T SFunction_o16;              // '<Root>/SupervisoryController'
    real_T SFunction_o17;              // '<Root>/SupervisoryController'
    real_T SFunction_o18;              // '<Root>/SupervisoryController'
    real_T SFunction_o19;              // '<Root>/SupervisoryController'
    real_T SFunction_o20;              // '<Root>/SupervisoryController'
    real_T SFunction_o21;              // '<Root>/SupervisoryController'
    real_T SFunction_o22;              // '<Root>/SupervisoryController'
    uint32_T RandSeed[3];              // '<S2>/excitation'
    uint16_T waypt;                    // '<Root>/SupervisoryController'
    uint16_T trajSize;                 // '<Root>/SupervisoryController'
    uint8_T is_EventHandler;           // '<Root>/SupervisoryController'
    uint8_T is_active_c6_SupervisoryControl;// '<Root>/SupervisoryController'
    boolean_T enAdapt_[6];             // '<Root>/SupervisoryController'
    boolean_T Memory_PreviousInput[246];// '<S9>/Memory'
    boolean_T evDone;                  // '<Root>/SupervisoryController'
    boolean_T icLoad;                  // '<S39>/MemoryX'
    boolean_T icLoad_j;                // '<S39>/MemoryP'
    boolean_T MeasurementUpdate_MODE;  // '<S58>/MeasurementUpdate'
  };

  // Zero-crossing (trigger) state
  struct PrevZCX {
    ZCSigState SupervisoryController_Trig_ZCE;// '<Root>/SupervisoryController'
    ZCE_paramEst1 paramEst2;           // '<S1>/paramEst2'
    ZCSigState MemoryX_Reset_ZCE;      // '<S39>/MemoryX'
    ZCSigState MemoryP_Reset_ZCE;      // '<S39>/MemoryP'
    ZCE_paramEst1 paramEst1_o;         // '<S1>/paramEst1'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU {
    real_T y[6];                       // '<Root>/y'
    real_T ymax[6];                    // '<Root>/ymax'
    real_T y0[6];                      // '<Root>/y0'
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
    real_T P_j[576];                   // '<Root>/P'
    boolean_T requestEvent;            // '<Root>/requestEvent'
  };

  // Parameters for system: '<S1>/paramEst1'
  struct P_paramEst1 {
    real_T theta_Y0;                   // Computed Parameter: theta_Y0
                                          //  Referenced by: '<S3>/theta'

    real_T P_Y0;                       // Computed Parameter: P_Y0
                                          //  Referenced by: '<S3>/P'

    real_T err_Y0;                     // Computed Parameter: err_Y0
                                          //  Referenced by: '<S3>/err'

    real_T UnitDelay3_InitialCondition;// Expression: 0
                                          //  Referenced by: '<S85>/Unit Delay3'

  };

  // Parameters (default storage)
  struct P {
    event_bus nullEv;                  // Variable: nullEv
                                          //  Referenced by: '<Root>/SupervisoryController'

    real_T Aod[144];                   // Variable: Aod
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Bod[72];                    // Variable: Bod
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Cod[72];                    // Variable: Cod
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dmn[36];                    // Variable: Dmn
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T Dod[36];                    // Variable: Dod
                                          //  Referenced by: '<S8>/MATLAB Function'

    real_T beta;                       // Variable: beta
                                          //  Referenced by:
                                          //    '<S2>/Gain'
                                          //    '<S2>/Gain1'

    real_T dt;                         // Variable: dt
                                          //  Referenced by:
                                          //    '<Root>/SupervisoryController'
                                          //    '<S2>/MATLAB Function'

    real_T lpfDen;                     // Variable: lpfDen
                                          //  Referenced by: '<S2>/Discrete Filter1'

    real_T lpfNum[60];                 // Variable: lpfNum
                                          //  Referenced by: '<S2>/Discrete Filter1'

    real_T uwt0[3];                    // Variable: uwt0
                                          //  Referenced by: '<S2>/Delay1'

    real_T ywt0[6];                    // Variable: ywt0
                                          //  Referenced by: '<S2>/Delay'

    real_T Lykyhatkk1_Y0;              // Expression: 0
                                          //  Referenced by: '<S82>/L*(y[k]-yhat[k|k-1])'

    real_T u_Y0;                       // Computed Parameter: u_Y0
                                          //  Referenced by: '<S2>/u'

    real_T ywt_Y0;                     // Computed Parameter: ywt_Y0
                                          //  Referenced by: '<S2>/ywt'

    real_T yhat_Y0;                    // Computed Parameter: yhat_Y0
                                          //  Referenced by: '<S2>/yhat'

    real_T r_Y0;                       // Computed Parameter: r_Y0
                                          //  Referenced by: '<S2>/r_'

    real_T G_zero_Value;               // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/G_zero'

    real_T LastPcov_InitialCondition[324];// Expression: lastPcov
                                             //  Referenced by: '<S9>/LastPcov'

    real_T duwt_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/du.wt_zero'

    real_T extmv_zero_Value[3];        // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/ext.mv_zero'

    real_T extmv_scale_Gain[3];        // Expression: RMVscale
                                          //  Referenced by: '<S9>/ext.mv_scale'

    real_T last_mv_InitialCondition[3];// Expression: lastu+uoff
                                          //  Referenced by: '<S9>/last_mv'

    real_T Constant12_Value[36];       // Expression: G.C
                                          //  Referenced by: '<S2>/Constant12'

    real_T Constant13_Value[18];       // Expression: G.D
                                          //  Referenced by: '<S2>/Constant13'

    real_T Constant11_Value;           // Expression: 1
                                          //  Referenced by: '<S2>/Constant11'

    real_T Constant2_Value[6];         // Expression: zeros(2*no, 1)
                                          //  Referenced by: '<S2>/Constant2'

    real_T Constant1_Value[12];        // Expression: zeros(size(God.A,1),1)
                                          //  Referenced by: '<S8>/Constant1'

    real_T X0_Value[18];               // Expression: pInitialization.X0
                                          //  Referenced by: '<S39>/X0'

    real_T ym_zero_Value[6];           // Expression: zeros(nym,1)
                                          //  Referenced by: '<S9>/ym_zero'

    real_T md_zero_Value;              // Expression: zeros(1,1)
                                          //  Referenced by: '<S5>/md_zero'

    real_T umin_zero_Value[3];         // Expression: zeros(3,1)
                                          //  Referenced by: '<S5>/umin_zero'

    real_T Gain2_Gain;                 // Expression: -2
                                          //  Referenced by: '<S2>/Gain2'

    real_T Gain3_Gain;                 // Expression: 2
                                          //  Referenced by: '<S2>/Gain3'

    real_T E_zero_Value[3];            // Expression: zeros(1,3)
                                          //  Referenced by: '<S5>/E_zero'

    real_T umin_scale4_Gain[3];    // Expression: MVscale(:,ones(1,max(nCC,1)))'
                                      //  Referenced by: '<S9>/umin_scale4'

    real_T F_zero_Value[6];            // Expression: zeros(1,6)
                                          //  Referenced by: '<S5>/F_zero'

    real_T ymin_scale1_Gain[6];     // Expression: Yscale(:,ones(1,max(nCC,1)))'
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

    real_T P0_Value[324];              // Expression: pInitialization.P0
                                          //  Referenced by: '<S39>/P0'

    real_T Constant1_Value_c;          // Expression: 1
                                          //  Referenced by: '<S2>/Constant1'

    real_T H_Value[108];               // Expression: pInitialization.H
                                          //  Referenced by: '<S39>/H'

    real_T G_Value[324];               // Expression: pInitialization.G
                                          //  Referenced by: '<S39>/G'

    real_T u_scale_Gain[3];            // Expression: MVscale
                                          //  Referenced by: '<S9>/u_scale'

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

    int32_T FixedHorizonOptimizer_Ndis;// Expression: Ndis
                                          //  Referenced by: '<S37>/FixedHorizonOptimizer'

    boolean_T Memory_InitialCondition[246];// Expression: iA
                                              //  Referenced by: '<S9>/Memory'

    boolean_T isSqrtUsed_Value;        // Expression: pInitialization.isSqrtUsed
                                          //  Referenced by: '<S80>/isSqrtUsed'

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
      void* dataAddress[66];
      int32_T* vardimsAddress[66];
      RTWLoggingFcnPtr loggingPtrs[66];
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
  void binary_expand_op_n(real_T in1[24], int32_T in2, const real_T in3[24],
    int32_T in4, int32_T in5, const real_T in6[24], int32_T in7, int32_T in8);
  int32_T xpotrf(real_T b_A[16]);
  real_T minimum(const real_T x[4]);
  void mpc_checkhessian(real_T b_H[16], real_T L[16], real_T *BadH);
  void trisolve(const real_T b_A[16], real_T b_B[16]);
  real_T norm(const real_T x[4]);
  real_T maximum(const real_T x[4]);
  real_T xnrm2(int32_T n, const real_T x[16], int32_T ix0);
  void xgemv(int32_T b_m, int32_T n, const real_T b_A[16], int32_T ia0, const
             real_T x[16], int32_T ix0, real_T y[4]);
  void xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const real_T y
             [4], real_T b_A[16], int32_T ia0);
  void KWIKfactor(const real_T b_Ac[984], const int32_T iC[246], int32_T nA,
                  const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
                  int32_T n, real_T RLinv[16], real_T *Status);
  void DropConstraint(int32_T kDrop, boolean_T iA[246], int32_T *nA, int32_T iC
                      [246]);
  void qpkwik(const real_T b_Linv[16], const real_T b_Hinv[16], const real_T f[4],
              const real_T b_Ac[984], const real_T b[246], boolean_T iA[246],
              int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[246],
              int32_T *status);
  void mpcblock_optimizer(const real_T rseq[120], const real_T vseq[21], const
    real_T umax[3], const real_T ymin[6], const real_T ymax[6], const real_T x
    [18], const real_T old_u[3], const boolean_T iA[246], const real_T b_Mlim
    [246], real_T b_Mx[4428], real_T b_Mu1[738], real_T b_Mv[5166], const real_T
    b_utarget[60], const real_T b_uoff[3], const real_T b_yoff[6], real_T b_H[16],
    real_T b_Ac[984], const real_T ywt[6], const real_T uwt[3], const real_T
    b_Wdu[3], const real_T b_Jm[180], const real_T b_I1[180], const real_T b_A
    [324], const real_T Bu[1134], const real_T Bv[378], const real_T b_C[108],
    const real_T Dv[126], const int32_T b_Mrows[246], const real_T b_RYscale[6],
    const real_T b_RMVscale[3], real_T u[3], real_T useq[63], real_T *status,
    boolean_T iAout[246]);
  void mrdiv(const real_T A[108], const real_T B_0[36], real_T Y[108]);

  // Real-Time Model
  RT_MODEL rtM;
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
//  Block '<S8>/Display' : Unused code path elimination
//  Block '<S8>/Display1' : Unused code path elimination
//  Block '<S48>/Data Type Duplicate' : Unused code path elimination
//  Block '<S49>/Data Type Duplicate' : Unused code path elimination
//  Block '<S51>/Data Type Duplicate' : Unused code path elimination
//  Block '<S52>/Data Type Duplicate' : Unused code path elimination
//  Block '<S55>/Data Type Duplicate' : Unused code path elimination
//  Block '<S56>/Data Type Duplicate' : Unused code path elimination
//  Block '<S64>/CheckSignalProperties' : Unused code path elimination
//  Block '<S65>/CheckSignalProperties' : Unused code path elimination
//  Block '<S66>/CheckSignalProperties' : Unused code path elimination
//  Block '<S67>/CheckSignalProperties' : Unused code path elimination
//  Block '<S68>/CheckSignalProperties' : Unused code path elimination
//  Block '<S71>/CheckSignalProperties' : Unused code path elimination
//  Block '<S73>/CheckSignalProperties' : Unused code path elimination
//  Block '<S74>/CheckSignalProperties' : Unused code path elimination
//  Block '<S75>/CheckSignalProperties' : Unused code path elimination
//  Block '<S77>/CheckSignalProperties' : Unused code path elimination
//  Block '<S78>/CheckSignalProperties' : Unused code path elimination
//  Block '<S85>/Display' : Unused code path elimination
//  Block '<S85>/Display1' : Unused code path elimination
//  Block '<S85>/Display2' : Unused code path elimination
//  Block '<S87>/Display' : Unused code path elimination
//  Block '<S3>/Scope' : Unused code path elimination
//  Block '<S89>/Display' : Unused code path elimination
//  Block '<S89>/Display1' : Unused code path elimination
//  Block '<S89>/Display2' : Unused code path elimination
//  Block '<S91>/Display' : Unused code path elimination
//  Block '<S4>/Scope' : Unused code path elimination
//  Block '<S9>/Reshape' : Reshape block reduction
//  Block '<S9>/Reshape1' : Reshape block reduction
//  Block '<S9>/Reshape2' : Reshape block reduction
//  Block '<S9>/Reshape3' : Reshape block reduction
//  Block '<S9>/Reshape4' : Reshape block reduction
//  Block '<S9>/Reshape5' : Reshape block reduction
//  Block '<S51>/Conversion' : Eliminate redundant data type conversion
//  Block '<S55>/Conversion' : Eliminate redundant data type conversion
//  Block '<S58>/Reshape' : Reshape block reduction
//  Block '<S39>/ReshapeX0' : Reshape block reduction
//  Block '<S39>/Reshapeu' : Reshape block reduction
//  Block '<S39>/Reshapexhat' : Reshape block reduction
//  Block '<S39>/Reshapey' : Reshape block reduction
//  Block '<S39>/Reshapeyhat' : Reshape block reduction


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
//  '<S3>'   : 'T_junction_mpc/SupervisoryController/paramEst1'
//  '<S4>'   : 'T_junction_mpc/SupervisoryController/paramEst2'
//  '<S5>'   : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller'
//  '<S6>'   : 'T_junction_mpc/SupervisoryController/ampc/MATLAB Function'
//  '<S7>'   : 'T_junction_mpc/SupervisoryController/ampc/MATLAB Function2'
//  '<S8>'   : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)'
//  '<S9>'   : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC'
//  '<S10>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check'
//  '<S11>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check A'
//  '<S12>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check B'
//  '<S13>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check C'
//  '<S14>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check D'
//  '<S15>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check DX'
//  '<S16>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check U'
//  '<S17>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check X'
//  '<S18>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check Y'
//  '<S19>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check1'
//  '<S20>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Matrix Signal Check2'
//  '<S21>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check'
//  '<S22>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check1'
//  '<S23>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check2'
//  '<S24>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check3'
//  '<S25>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check4'
//  '<S26>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check5'
//  '<S27>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check6'
//  '<S28>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check7'
//  '<S29>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Preview Signal Check8'
//  '<S30>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check'
//  '<S31>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check1'
//  '<S32>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Scalar Signal Check2'
//  '<S33>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check'
//  '<S34>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check1'
//  '<S35>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/MPC Vector Signal Check6'
//  '<S36>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/moorx'
//  '<S37>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer'
//  '<S38>'  : 'T_junction_mpc/SupervisoryController/ampc/Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer'
//  '<S39>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2'
//  '<S40>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/MATLAB Function'
//  '<S41>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculatePL'
//  '<S42>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculateYhat'
//  '<S43>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator'
//  '<S44>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionA'
//  '<S45>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionB'
//  '<S46>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionC'
//  '<S47>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionD'
//  '<S48>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionG'
//  '<S49>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionH'
//  '<S50>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionN'
//  '<S51>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP'
//  '<S52>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionP0'
//  '<S53>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionQ'
//  '<S54>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionR'
//  '<S55>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX'
//  '<S56>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionX0'
//  '<S57>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/DataTypeConversionu'
//  '<S58>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer'
//  '<S59>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ReducedQRN'
//  '<S60>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionP0'
//  '<S61>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ'
//  '<S62>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR'
//  '<S63>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/UseCurrentEstimator'
//  '<S64>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkA'
//  '<S65>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkB'
//  '<S66>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkC'
//  '<S67>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkD'
//  '<S68>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkEnable'
//  '<S69>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkG'
//  '<S70>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkH'
//  '<S71>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkN'
//  '<S72>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkP0'
//  '<S73>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkQ'
//  '<S74>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkR'
//  '<S75>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkReset'
//  '<S76>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checkX0'
//  '<S77>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checku'
//  '<S78>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/checky'
//  '<S79>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CalculatePL/Discrete-Time KF - Calculate PLMZ'
//  '<S80>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput'
//  '<S81>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn'
//  '<S82>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/Observer/MeasurementUpdate'
//  '<S83>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionQ/ScalarExpansion'
//  '<S84>'  : 'T_junction_mpc/SupervisoryController/ampc/State Estimator OD (KF)/Kalman Filter2/ScalarExpansionR/ScalarExpansion'
//  '<S85>'  : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)'
//  '<S86>'  : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1'
//  '<S87>'  : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/RLS'
//  '<S88>'  : 'T_junction_mpc/SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function'
//  '<S89>'  : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)'
//  '<S90>'  : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/MATLAB Function1'
//  '<S91>'  : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/RLS'
//  '<S92>'  : 'T_junction_mpc/SupervisoryController/paramEst2/Param Estimator (RLS)/RLS/MATLAB Function'


//-
//  Requirements for '<Root>': SupervisoryController


#endif                                 // RTW_HEADER_SupervisoryController_h_

//
// File trailer for generated code.
//
// [EOF]
//
