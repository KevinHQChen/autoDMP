//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.783
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Fri May 12 13:35:40 2023
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
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Aod2
  //  Referenced by: '<S244>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Bod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.025,

  // Variable: Bod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Bod2
  //  Referenced by: '<S244>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Cod0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Cod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Cod2
  //  Referenced by: '<S244>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn0
  //  Referenced by: '<S7>/MATLAB Function'

  1.0,

  // Variable: Dmn1
  //  Referenced by: '<S112>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn2
  //  Referenced by: '<S244>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dod0
  //  Referenced by: '<S7>/MATLAB Function'

  0.0,

  // Variable: Dod1
  //  Referenced by: '<S112>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod2
  //  Referenced by: '<S244>/MATLAB Function'

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

    // Expression: zeros(1,1)
    //  Referenced by: '<S242>/G_zero'

    0.0,

    // Expression: lastPcov
    //  Referenced by: '<S245>/LastPcov'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
      0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S242>/y.wt_zero'

    { 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/u.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/du.wt_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/ext.mv_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S245>/ext.mv_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: lastu+uoff
    //  Referenced by: '<S245>/last_mv'

    { 0.0, 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S243>/Unit Delay3'

    0.0,

    // Expression: 0
    //  Referenced by: '<S243>/Unit Delay2'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S276>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S276>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S276>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S276>/InitialRegressors'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S276>/InitialParameters'

    { 0.99996947351391963, 0.70746101068200651, -0.34920656360027097,
      -0.34258808392522866 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S276>/InitialCovariance'

    { 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
      0.0, 100.0 },

    // Expression: 0
    //  Referenced by: '<S243>/Unit Delay6'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S277>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S277>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S277>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S277>/InitialRegressors'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S277>/InitialParameters'

    { 0.99995521109485208, -0.35521300513038595, -0.074930066385699018,
      0.43035445788055782 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S277>/InitialCovariance'

    { 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
      0.0, 100.0 },

    // Expression: G2.C
    //  Referenced by: '<S243>/Constant12'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: G2.D
    //  Referenced by: '<S243>/Constant13'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: 1
    //  Referenced by: '<S243>/Constant11'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S243>/Constant2'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S244>/Constant'

    { 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S328>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

    // Expression: zeros(nym,1)
    //  Referenced by: '<S245>/ym_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S242>/md_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/umin_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/umax_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S242>/ymin_zero'

    { 0.0, 0.0 },

    // Expression: zeros(2,1)
    //  Referenced by: '<S242>/ymax_zero'

    { 0.0, 0.0 },

    // Expression: zeros(1,3)
    //  Referenced by: '<S242>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S245>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S242>/F_zero'

    { 0.0, 0.0 },

    // Expression: Yscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S245>/ymin_scale1'

    { 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S242>/S_zero'

    0.0,

    // Expression: MDscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S245>/ymin_scale2'

    1.0,

    // Expression: zeros(1,1)
    //  Referenced by: '<S242>/switch_zero'

    0.0,

    // Expression: zeros(3,1)
    //  Referenced by: '<S242>/mv.target_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: RMVscale
    //  Referenced by: '<S245>/uref_scale'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,1)
    //  Referenced by: '<S242>/ecr.wt_zero'

    0.0,

    // Expression: pInitialization.P0
    //  Referenced by: '<S328>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: 1
    //  Referenced by: '<S4>/Constant'

    1.0,

    // Expression: pInitialization.H
    //  Referenced by: '<S328>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S328>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: MVscale
    //  Referenced by: '<S245>/u_scale'

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
    //  Referenced by: '<S273>/FixedHorizonOptimizer'

    0,

    // Expression: iA
    //  Referenced by: '<S245>/Memory'

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
    //  Referenced by: '<S295>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S320>/Delay'

    true,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S369>/isSqrtUsed'

    false,

    // Expression: true
    //  Referenced by: '<S4>/Constant1'

    true,

    // Expression: false()
    //  Referenced by: '<S295>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S320>/Constant'

    false,

    // Start of '<S347>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S371>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S347>/MeasurementUpdate'
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

    { 0.0, 0.0, 0.0 },

    // Expression: 0
    //  Referenced by: '<S111>/Unit Delay3'

    0.0,

    // Expression: 0
    //  Referenced by: '<S111>/Unit Delay2'

    0.0,

    // Expression: initializationParams.adg1
    //  Referenced by: '<S144>/Forgetting Factor'

    1.0,

    // Expression: initializationParams.adg2
    //  Referenced by: '<S144>/Normalization Bias'

    0.0,

    // Expression: initializationParams.initialOutputs
    //  Referenced by: '<S144>/InitialOutputs'

    0.0,

    // Expression: initializationParams.initialRegressors
    //  Referenced by: '<S144>/InitialRegressors'

    0.0,

    // Expression: initializationParams.theta0
    //  Referenced by: '<S144>/InitialParameters'

    { 1.0000306206339091, -0.34290108520429963, 0.43337106536580761,
      -0.079433657162736954 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S144>/InitialCovariance'

    { 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
      0.0, 100.0 },

    // Expression: 0
    //  Referenced by: '<S111>/Unit Delay6'

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

    // Expression: initializationParams.theta0
    //  Referenced by: '<S145>/InitialParameters'

    { 0.99997520047355148, -0.3552101510579076, -0.0749598690059731,
      0.43031643660853247 },

    // Expression: initializationParams.L0
    //  Referenced by: '<S145>/InitialCovariance'

    { 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0,
      0.0, 100.0 },

    // Expression: 0
    //  Referenced by: '<S111>/Switch1'

    0.0,

    // Expression: 0
    //  Referenced by: '<S111>/Switch'

    0.0,

    // Expression: G1.C
    //  Referenced by: '<S111>/Constant12'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: G1.D
    //  Referenced by: '<S111>/Constant13'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: 1
    //  Referenced by: '<S111>/Constant11'

    1.0,

    // Expression: [0;0]
    //  Referenced by: '<S111>/Constant2'

    { 0.0, 0.0 },

    // Expression: [0;0]
    //  Referenced by: '<S112>/Constant'

    { 0.0, 0.0 },

    // Expression: pInitialization.X0
    //  Referenced by: '<S196>/X0'

    { 0.0, 0.0, 0.0, 0.0 },

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

    // Expression: zeros(1,3)
    //  Referenced by: '<S110>/E_zero'

    { 0.0, 0.0, 0.0 },

    // Expression: MVscale(:,ones(1,max(nCC,1)))'
    //  Referenced by: '<S113>/umin_scale4'

    { 1.0, 1.0, 1.0 },

    // Expression: zeros(1,2)
    //  Referenced by: '<S110>/F_zero'

    { 0.0, 0.0 },

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

    // Expression: pInitialization.P0
    //  Referenced by: '<S196>/P0'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: 1
    //  Referenced by: '<S3>/Constant'

    1.0,

    // Expression: pInitialization.H
    //  Referenced by: '<S196>/H'

    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S196>/G'

    { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
      1.0 },

    // Expression: MVscale
    //  Referenced by: '<S113>/u_scale'

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
    //  Referenced by: '<S163>/Delay'

    true,

    // Expression: true()
    //  Referenced by: '<S188>/Delay'

    true,

    // Expression: pInitialization.isSqrtUsed
    //  Referenced by: '<S237>/isSqrtUsed'

    false,

    // Expression: true
    //  Referenced by: '<S3>/Constant1'

    true,

    // Expression: false()
    //  Referenced by: '<S163>/Constant'

    false,

    // Expression: false()
    //  Referenced by: '<S188>/Constant'

    false,

    // Start of '<S215>/MeasurementUpdate'
    {
      // Expression: 0
      //  Referenced by: '<S239>/L*(y[k]-yhat[k|k-1])'

      0.0
    }
    // End of '<S215>/MeasurementUpdate'
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

    // Expression: G0.D
    //  Referenced by: '<S6>/Constant13'

    { 0.0, 0.0, 0.0 },

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
    //  Referenced by: '<S64>/X0'

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
    //  Referenced by: '<S64>/H'

    { 0.0, 0.0 },

    // Expression: pInitialization.G
    //  Referenced by: '<S64>/G'

    { 1.0, 0.0, 0.0, 1.0 },

    // Expression: 1
    //  Referenced by: '<S2>/Constant'

    1.0,

    // Expression: pInitialization.P0
    //  Referenced by: '<S64>/P0'

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
    //  Referenced by: '<S105>/isSqrtUsed'

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
