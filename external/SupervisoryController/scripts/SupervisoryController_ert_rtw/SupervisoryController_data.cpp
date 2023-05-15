//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.854
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sun May 14 23:58:51 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#include "SupervisoryController.h"

// Block parameters (default storage)
SupervisoryController::P SupervisoryController::rtP{
  // Variable: nullEv
  //  Referenced by: '<Root>/SupervisoryController'

  {
    0.0,
    0.0,

    { 0.0, 0.0, 0.0 },
    0.0,
    0.0,

    { false, false, false },

    { false, false, false }
  },

  // Variable: Aod0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Aod1
  //  Referenced by: '<S113>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Aod2
  //  Referenced by: '<S247>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Bod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.025,

  // Variable: Bod1
  //  Referenced by: '<S113>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Bod2
  //  Referenced by: '<S247>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Cod0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Cod1
  //  Referenced by: '<S113>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Cod2
  //  Referenced by: '<S247>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Dmn1
  //  Referenced by: '<S113>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn2
  //  Referenced by: '<S247>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.0,

  // Variable: Dod1
  //  Referenced by: '<S113>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod2
  //  Referenced by: '<S247>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0 },

  // Variable: dt
  //  Referenced by: '<Root>/SupervisoryController'

  0.025,

  // Start of '<S1>/State2.controlLaw.AMPC2'
  {
    // Computed Parameter: u_Y0
    //  Referenced by: '<S4>/u'

    0.0,

    // Computed Parameter: yhat_Y0
    //  Referenced by: '<S4>/yhat'

    0.0,

    // Computed Parameter: params_Y0
    //  Referenced by: '<S4>/params'

    0.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S245>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S248>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S245>/y.wt_zero'

    { 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S248>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S248>/last_mv'

    { 0.0, 0.0, 0.0 },

    // Expression: G2.A
    //  Referenced by: '<S246>/Constant10'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: 0
    //  Referenced by: '<S246>/Unit Delay7'

    0.0,

    // Expression: 0
    //  Referenced by: '<S246>/Unit Delay1'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S279>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S279>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S279>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S279>/InitialRegressors'

    0.0,

    // Expression: 1e4
    //  Referenced by: '<S246>/Constant'

    10000.0,

    // Expression: G2.B(1,:)
    //  Referenced by: '<S246>/Constant7'

    { 1.0, -0.5, -0.5 },

    // Expression: initializationParams.adg1
    //  Referenced by: '<S280>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S280>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S280>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S280>/InitialRegressors'

    0.0,

    // Expression: G2.B(2,:)
    //  Referenced by: '<S246>/Constant8'

    { -0.5, -0.5, 1.0 },

    // Expression: G2.C
    //  Referenced by: '<S246>/Constant4'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: G2.D
    //  Referenced by: '<S246>/Constant5'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: 1
    //  Referenced by: '<S246>/Constant3'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S246>/Constant6'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S247>/Constant'

    { 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S333>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S248>/ym_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S245>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S245>/ymin_zero'

    { 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S245>/ymax_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,3)
    //  Referenced by: '<S245>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S248>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S245>/F_zero'

    { 0.0, 0.0 },

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S248>/ymin_scale1'

    { 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S245>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S248>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S245>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S245>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S248>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S245>/ecr.wt_zero'

    0.0,

    // Expression: pInitialization.P0
    //  Referenced by: '<S333>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: 1
    //  Referenced by: '<S4>/Constant'

    1.0,

    // Expression: pInitialization.H
    //  Referenced by: '<S333>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S333>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: MVscale
    //  Referenced by: '<S248>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: 0
    //  Referenced by: '<S4>/Constant2'

    0.0,

    // Expression: [0 0 0]
    //  Referenced by: '<S4>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S4>/Measurement Noise'

    { 1.0, 1.0, 1.0 },

    // Expression: 12345
    //  Referenced by: '<S4>/Measurement Noise'

    12345.0,

    // Expression: Ndis
    //  Referenced by: '<S276>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S248>/Memory'

    { false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false },

    // Expression: true()
    //  Referenced by: '<S298>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S324>/Delay'

    true,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S374>/isSqrtUsed'

    false,

    // Expression: true
    //  Referenced by: '<S4>/Constant1'

    true,

    // Expression: false()
    //  Referenced by: '<S298>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S324>/Constant'

    false,

    // Start of '<S352>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S376>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S352>/MeasurementUpdate'
  }
  ,

  // End of '<S1>/State2.controlLaw.AMPC2'

  // Start of '<S1>/State1.controlLaw.AMPC1'
  {
    // Computed Parameter: u_Y0
    //  Referenced by: '<S3>/u'

    0.0,

    // Computed Parameter: yhat_Y0
    //  Referenced by: '<S3>/yhat'

    0.0,

    // Computed Parameter: params_Y0
    //  Referenced by: '<S3>/params'

    0.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S111>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S114>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S111>/y.wt_zero'

    { 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S114>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S114>/last_mv'

    { 0.0, 0.0, 0.0 },

    // Expression: G1.A
    //  Referenced by: '<S112>/Constant5'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: 0
    //  Referenced by: '<S112>/Unit Delay3'

    0.0,

    // Expression: 0
    //  Referenced by: '<S112>/Unit Delay2'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S145>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S145>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S145>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S145>/InitialRegressors'

    0.0,

    // Expression: 1e4
    //  Referenced by: '<S112>/Constant'

    10000.0,

    // Expression: G1.B(1,:)
    //  Referenced by: '<S112>/Constant3'

    { -0.5, 1.0, -0.5 },

    // Expression: initializationParams.adg1
    //  Referenced by: '<S146>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S146>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S146>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S146>/InitialRegressors'

    0.0,

    // Expression: G1.B(2,:)
    //  Referenced by: '<S112>/Constant4'

    { -0.5, -0.5, 1.0 },

    // Expression: G1.C
    //  Referenced by: '<S112>/Constant12'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: G1.D
    //  Referenced by: '<S112>/Constant13'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: 1
    //  Referenced by: '<S112>/Constant11'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S112>/Constant2'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S113>/Constant'

    { 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S199>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S114>/ym_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S111>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S111>/ymin_zero'

    { 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S111>/ymax_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,3)
    //  Referenced by: '<S111>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S114>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S111>/F_zero'

    { 0.0, 0.0 },

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S114>/ymin_scale1'

    { 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S111>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S114>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S111>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S111>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S114>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S111>/ecr.wt_zero'

    0.0,

    // Expression: pInitialization.P0
    //  Referenced by: '<S199>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: 1
    //  Referenced by: '<S3>/Constant'

    1.0,

    // Expression: pInitialization.H
    //  Referenced by: '<S199>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S199>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: MVscale
    //  Referenced by: '<S114>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: 0
    //  Referenced by: '<S3>/Constant2'

    0.0,

    // Expression: [0 0 0]
    //  Referenced by: '<S3>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S3>/Measurement Noise'

    { 1.0, 1.0, 1.0 },

    // Expression: 12345
    //  Referenced by: '<S3>/Measurement Noise'

    12345.0,

    // Expression: Ndis
    //  Referenced by: '<S142>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S114>/Memory'

    { false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false },

    // Expression: true()
    //  Referenced by: '<S164>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S190>/Delay'

    true,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S240>/isSqrtUsed'

    false,

    // Expression: true
    //  Referenced by: '<S3>/Constant1'

    true,

    // Expression: false()
    //  Referenced by: '<S164>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S190>/Constant'

    false,

    // Start of '<S218>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S242>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S218>/MeasurementUpdate'
  }
  ,

  // End of '<S1>/State1.controlLaw.AMPC1'

  // Start of '<S1>/State0.controlLaw.AMPC0'
  {
    // Expression: 0
    //  Referenced by: '<S108>/L*(y[k]-yhat[k|k-1])'

    0.0,

    // Computed Parameter: u_Y0
    //  Referenced by: '<S2>/u'

    0.0,

    // Computed Parameter: yhat_Y0
    //  Referenced by: '<S2>/yhat'

    0.0,

    // Computed Parameter: params_Y0
    //  Referenced by: '<S2>/params'

    0.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S8>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/y.wt_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S8>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S8>/last_mv'

    { 0.0, 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S6>/Unit Delay2'

    0.0,

    // Expression: 0
    //  Referenced by: '<S6>/Unit Delay3'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S38>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S38>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S38>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S38>/InitialRegressors'

    0.0,

    // Expression: 1e4
    //  Referenced by: '<S6>/Constant'

    10000.0,

    // Expression: G0.B
    //  Referenced by: '<S6>/Constant1'

    { 1.0, -0.5, -0.5 },

    // Expression: G0.D
    //  Referenced by: '<S6>/Constant13'

    { 0.0, 0.0, 0.0 },

    // Expression: G0.A
    //  Referenced by: '<S6>/Constant2'

    1.0,

    // Expression: G0.C
    //  Referenced by: '<S6>/Constant12'

    1.0,

    // Expression: 1
    //  Referenced by: '<S6>/Constant11'

    1.0,

    // Expression: 0
    //  Referenced by: '<S6>/Constant10'

    0.0,

    // Expression: 0
    //  Referenced by: '<S7>/Constant'

    0.0,

    // Expression: pInitialization.X0
    //  Referenced by: '<S65>/X0'

    { 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S8>/ym_zero'

    0.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/ymin_zero'

    0.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/ymax_zero'

    0.0,

    // Expression: zeros(1,3)
    //  Referenced by: '<S5>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S8>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/F_zero'

    0.0,

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S8>/ymin_scale1'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S8>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S5>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S8>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/ecr.wt_zero'

    0.0,

    // Expression: pInitialization.H
    //  Referenced by: '<S65>/H'

    { 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S65>/G'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: 1
    //  Referenced by: '<S2>/Constant'

    1.0,

    // Expression: pInitialization.P0
    //  Referenced by: '<S65>/P0'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: MVscale
    //  Referenced by: '<S8>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: [0 0 0]
    //  Referenced by: '<S2>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S2>/Measurement Noise'

    { 1.0, 1.0, 1.0 },

    // Expression: 12345
    //  Referenced by: '<S2>/Measurement Noise'

    12345.0,

    // Expression: Ndis
    //  Referenced by: '<S36>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S8>/Memory'

    { false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false, false, false, false, false,
      false, false, false, false, false, false },

    // Expression: true()
    //  Referenced by: '<S56>/Delay'

    true,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S106>/isSqrtUsed'

    false,

    // Expression: false()
    //  Referenced by: '<S56>/Constant'

    false
  }
  // End of '<S1>/State0.controlLaw.AMPC0'
};

//
// File trailer for generated code.
//
// [EOF]
//
