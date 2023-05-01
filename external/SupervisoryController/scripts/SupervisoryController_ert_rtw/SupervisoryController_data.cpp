//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.714
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 13:52:35 2023
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
  // Variable: Aod0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Aod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Aod2
  //  Referenced by: '<S256>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Bod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.025,

  // Variable: Bod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Bod2
  //  Referenced by: '<S256>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Cod0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Cod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Cod2
  //  Referenced by: '<S256>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Dmn1
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn2
  //  Referenced by: '<S256>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.0,

  // Variable: Dod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod2
  //  Referenced by: '<S256>/MATLAB Function'

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

    // Expression: zeros(1,3)
    //  Referenced by: '<S254>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S254>/F_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S254>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S257>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S254>/y.wt_zero'

    { 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S257>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S257>/last_mv'

    { 40.0, 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S352>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S255>/Unit Delay5'

    0.0,

    // Expression: [0; 0]
    //  Referenced by: '<S255>/Unit Delay4'

    { 0.0, 0.0 },

    // Expression: initializationParams.adg1
    //  Referenced by: '<S290>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S290>/Normalization Bias'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S290>/InitialParameters'

    { -0.99996947351391963, 0.70746101068200651, -0.34920656360027097,
      -0.34258808392522866 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S290>/InitialCovariance'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: initializationParams.phiMemory0
    //  Referenced by: '<S290>/InitialPhiMemory'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: initializationParams.initialInputs
    //  Referenced by: '<S290>/InitialInputs'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S290>/InitialOutputs'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S291>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S291>/Normalization Bias'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S291>/InitialParameters'

    { -0.99996947351391963, -0.35521300513038595, -0.074930066385699018,
      0.43035445788055782 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S291>/InitialCovariance'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: initializationParams.phiMemory0
    //  Referenced by: '<S291>/InitialPhiMemory'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: initializationParams.initialInputs
    //  Referenced by: '<S291>/InitialInputs'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S291>/InitialOutputs'

    0.0,

    // Expression: 1
    //  Referenced by: '<S255>/Constant1'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S255>/Constant9'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S256>/Constant'

    { 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S257>/ym_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S254>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S254>/ymin_zero'

    { 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S254>/ymax_zero'

    { 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S257>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S257>/ymin_scale1'

    { 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S254>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S257>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S254>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S254>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S257>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S254>/ecr.wt_zero'

    0.0,

    // Expression: MVscale
    //  Referenced by: '<S257>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: 1
    //  Referenced by: '<S4>/Constant'

    1.0,

    // Expression: [0 0 0]
    //  Referenced by: '<S4>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S4>/Measurement Noise'

    { 0.70710678118654757, 0.70710678118654757, 0.70710678118654757 },

    // Expression: 12345
    //  Referenced by: '<S4>/Measurement Noise'

    12345.0,

    // Expression: pInitialization.G
    //  Referenced by: '<S352>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: pInitialization.H
    //  Referenced by: '<S352>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.P0
    //  Referenced by: '<S352>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: Ndis
    //  Referenced by: '<S285>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S257>/Memory'

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
    //  Referenced by: '<S315>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S344>/Delay'

    true,

    // Expression: false()
    //  Referenced by: '<S315>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S344>/Constant'

    false,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S393>/isSqrtUsed'

    false,

    // Start of '<S371>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S395>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S371>/MeasurementUpdate'
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

    // Expression: zeros(1,3)
    //  Referenced by: '<S110>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S110>/F_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S110>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S113>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S110>/y.wt_zero'

    { 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S113>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S113>/last_mv'

    { 40.0, 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S208>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S111>/Unit Delay5'

    0.0,

    // Expression: [0; 0]
    //  Referenced by: '<S111>/Unit Delay4'

    { 0.0, 0.0 },

    // Expression: initializationParams.adg1
    //  Referenced by: '<S146>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S146>/Normalization Bias'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S146>/InitialParameters'

    { -1.0000306206339091, -0.34290108520429963, 0.43337106536580761,
      -0.079433657162736954 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S146>/InitialCovariance'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: initializationParams.phiMemory0
    //  Referenced by: '<S146>/InitialPhiMemory'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: initializationParams.initialInputs
    //  Referenced by: '<S146>/InitialInputs'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S146>/InitialOutputs'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S147>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S147>/Normalization Bias'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S147>/InitialParameters'

    { -1.0000306206339091, -0.3552101510579076, -0.0749598690059731,
      0.43031643660853247 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S147>/InitialCovariance'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: initializationParams.phiMemory0
    //  Referenced by: '<S147>/InitialPhiMemory'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: initializationParams.initialInputs
    //  Referenced by: '<S147>/InitialInputs'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S147>/InitialOutputs'

    0.0,

    // Expression: 1
    //  Referenced by: '<S111>/Constant1'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S111>/Constant9'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S112>/Constant'

    { 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S113>/ym_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S110>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S110>/ymin_zero'

    { 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S110>/ymax_zero'

    { 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S113>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S113>/ymin_scale1'

    { 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S110>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S113>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S110>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S110>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S113>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S110>/ecr.wt_zero'

    0.0,

    // Expression: MVscale
    //  Referenced by: '<S113>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: 1
    //  Referenced by: '<S3>/Constant'

    1.0,

    // Expression: [0 0 0]
    //  Referenced by: '<S3>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S3>/Measurement Noise'

    { 0.70710678118654757, 0.70710678118654757, 0.70710678118654757 },

    // Expression: 12345
    //  Referenced by: '<S3>/Measurement Noise'

    12345.0,

    // Expression: pInitialization.G
    //  Referenced by: '<S208>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: pInitialization.H
    //  Referenced by: '<S208>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.P0
    //  Referenced by: '<S208>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: Ndis
    //  Referenced by: '<S141>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S113>/Memory'

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
    //  Referenced by: '<S171>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S200>/Delay'

    true,

    // Expression: false()
    //  Referenced by: '<S171>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S200>/Constant'

    false,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S249>/isSqrtUsed'

    false,

    // Start of '<S227>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S251>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S227>/MeasurementUpdate'
  }
  ,

  // End of '<S1>/State1.controlLaw.AMPC1'

  // Start of '<S1>/State0.controlLaw.AMPC0'
  {
    // Expression: 0
    //  Referenced by: '<S107>/L*(y[k]-yhat[k|k-1])'

    0.0,

    // Computed Parameter: u_Y0
    //  Referenced by: '<S2>/u'

    0.0,

    // Computed Parameter: yhat_Y0
    //  Referenced by: '<S2>/yhat'

    0.0,

    // Expression: zeros(1,3)
    //  Referenced by: '<S5>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S5>/F_zero'

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

    { 40.0, 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S64>/X0'

    { 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S6>/Unit Delay3'

    0.0,

    // Expression: 0
    //  Referenced by: '<S6>/Unit Delay2'

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

    // Expression: initializationParams.theta0
    //  Referenced by: '<S38>/InitialParameters'

    { 1.000032235800572, 0.7084943675087555, -0.34693844974269455,
      -0.3458694000362138 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S38>/InitialCovariance'

    { 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
      0.0, 100.0 },

    // Expression: G0.C
    //  Referenced by: '<S6>/Constant12'

    1.0,

    // Expression: G0.D
    //  Referenced by: '<S6>/Constant13'

    { 0.0, 0.0, 0.0 },

    // Expression: 1
    //  Referenced by: '<S6>/Constant11'

    1.0,

    // Expression: 0
    //  Referenced by: '<S6>/Constant10'

    0.0,

    // Expression: 0
    //  Referenced by: '<S7>/Constant'

    0.0,

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

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S8>/umin_scale4'

    { 1.0, 1.0, 1.0 },

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

    // Expression: MVscale
    //  Referenced by: '<S8>/u_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: 1
    //  Referenced by: '<S2>/Constant'

    1.0,

    // Expression: [0 0 0]
    //  Referenced by: '<S2>/Measurement Noise'

    { 0.0, 0.0, 0.0 },

    // Computed Parameter: MeasurementNoise_StdDev
    //  Referenced by: '<S2>/Measurement Noise'

    { 0.70710678118654757, 0.70710678118654757, 0.70710678118654757 },

    // Expression: 12345
    //  Referenced by: '<S2>/Measurement Noise'

    12345.0,

    // Expression: pInitialization.G
    //  Referenced by: '<S64>/G'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: pInitialization.H
    //  Referenced by: '<S64>/H'

    { 0.0, 0.0 },

    // Expression: pInitialization.P0
    //  Referenced by: '<S64>/P0'

    { 1.0, 0.0, 0.0, 1.0 },

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

    // Expression: false()
    //  Referenced by: '<S56>/Constant'

    false,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S105>/isSqrtUsed'

    false
  }
  // End of '<S1>/State0.controlLaw.AMPC0'
};

//
// File trailer for generated code.
//
// [EOF]
//
