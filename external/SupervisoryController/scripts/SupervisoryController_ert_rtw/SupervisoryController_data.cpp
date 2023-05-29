//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.1475
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May 29 02:28:03 2023
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

    { true, false, false },

    { true, false, false }
  },

  // Variable: Aod0
  //  Referenced by: '<S8>/MATLAB Function'

  1.0,

  // Variable: Aod1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S175>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Bod0
  //  Referenced by: '<S8>/MATLAB Function'

  0.025,

  // Variable: Bod1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S175>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.025 },

  // Variable: Cod0
  //  Referenced by: '<S8>/MATLAB Function'

  1.0,

  // Variable: Cod1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S175>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn0
  //  Referenced by: '<S8>/MATLAB Function'

  1.0,

  // Variable: Dmn1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S175>/MATLAB Function'

  { 1.0, 0.0, 0.0, 1.0 },

  // Variable: Dod0
  //  Referenced by: '<S8>/MATLAB Function'

  0.0,

  // Variable: Dod1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S175>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0 },

  // Variable: dt
  //  Referenced by: '<Root>/SupervisoryController'

  0.025,

  // Variable: forgettingFactor
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

  0.9975,

  // Variable: lpfDen
  //  Referenced by:
  //    '<S2>/Discrete Filter'
  //    '<S3>/Discrete Filter'
  //    '<S4>/Discrete Filter'

  1.0,

  // Variable: lpfNum
  //  Referenced by:
  //    '<S2>/Discrete Filter'
  //    '<S3>/Discrete Filter'
  //    '<S4>/Discrete Filter'

  { -0.0001395194037777206, 0.0001062405464448335, 0.00053156288262108274,
    0.00061498512134384043, -0.0001072340498580264, -0.0012869873067591802,
    -0.0016761062869331484, -0.00023795039586547393, 0.002358892157775442,
    0.0036213661265000964, 0.0013467247781293011, -0.0035727765134688375,
    -0.0067466463590714672, -0.0038194621120323639, 0.0045377009594818243,
    0.011326022139335446, 0.0085042852608373139, -0.0045719847433271926,
    -0.017685942779495263, -0.01676882841526425, 0.0025008456098981418,
    0.026580380967637394, 0.031629264420877129, 0.0042313858248733812,
    -0.040884247859049747, -0.063732792845230007, -0.025454456431879591,
    0.079274985582173219, 0.20972814939304377, 0.29979214373103996,
    0.29979214373103996, 0.20972814939304377, 0.079274985582173219,
    -0.025454456431879591, -0.063732792845230007, -0.040884247859049747,
    0.0042313858248733812, 0.031629264420877129, 0.026580380967637394,
    0.0025008456098981418, -0.01676882841526425, -0.017685942779495263,
    -0.0045719847433271926, 0.0085042852608373139, 0.011326022139335446,
    0.0045377009594818243, -0.0038194621120323639, -0.0067466463590714672,
    -0.0035727765134688375, 0.0013467247781293011, 0.0036213661265000964,
    0.002358892157775442, -0.00023795039586547393, -0.0016761062869331484,
    -0.0012869873067591802, -0.0001072340498580264, 0.00061498512134384043,
    0.00053156288262108274, 0.0001062405464448335, -0.0001395194037777206 },

  // Variable: chs0
  //  Referenced by: '<Root>/SupervisoryController'

  1U,

  // Variable: chs1
  //  Referenced by: '<Root>/SupervisoryController'

  { 2U, 3U },

  // Variable: chs2
  //  Referenced by: '<Root>/SupervisoryController'

  { 1U, 3U },

  // Expression: 0
  //  Referenced by: '<S84>/L*(y[k]-yhat[k|k-1])'

  0.0,

  // Computed Parameter: u_Y0
  //  Referenced by: '<S2>/u'

  0.0,

  // Computed Parameter: yhat_Y0
  //  Referenced by: '<S2>/yhat'

  0.0,

  // Computed Parameter: prms_Y0
  //  Referenced by: '<S2>/prms'

  0.0,

  // Computed Parameter: uref_Y0
  //  Referenced by: '<S2>/uref'

  0.0,

  // Computed Parameter: prmErr_Y0
  //  Referenced by: '<S2>/prmErr'

  0.0,

  // Computed Parameter: uclean_Y0
  //  Referenced by: '<S2>/uclean'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S5>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S9>/LastPcov'

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
  //  Referenced by: '<S9>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S9>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: 0
  //  Referenced by: '<S7>/Unit Delay2'

  0.0,

  // Expression: 0
  //  Referenced by: '<S7>/Unit Delay3'

  0.0,

  // Expression: 1e4*eye(2)
  //  Referenced by: '<S7>/Constant4'

  { 10000.0, 0.0, 0.0, 10000.0 },

  // Expression: G0.D
  //  Referenced by: '<S7>/Constant13'

  { 0.0, 0.0, 0.0 },

  // Expression: G0.A
  //  Referenced by: '<S7>/Constant3'

  1.0,

  // Expression: G0.C
  //  Referenced by: '<S7>/Constant12'

  1.0,

  // Expression: 1
  //  Referenced by: '<S7>/Constant11'

  1.0,

  // Expression: 0
  //  Referenced by: '<S7>/Constant2'

  0.0,

  // Expression: 0
  //  Referenced by: '<S8>/Constant'

  0.0,

  // Expression: pInitialization.X0
  //  Referenced by: '<S41>/X0'

  { 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S9>/ym_zero'

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
  //  Referenced by: '<S9>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S5>/F_zero'

  0.0,

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S9>/ymin_scale1'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S5>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S9>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S5>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S5>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S9>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S5>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S41>/H'

  { 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S41>/G'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S2>/Constant'

  1.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S41>/P0'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S9>/u_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: [0 0 0]
  //  Referenced by: '<S2>/excitation'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: excitation_StdDev
  //  Referenced by: '<S2>/excitation'

  { 1.0, 1.0, 1.0 },

  // Expression: [12345 12345 12345]
  //  Referenced by: '<S2>/excitation'

  { 12345.0, 12345.0, 12345.0 },

  // Expression: 0
  //  Referenced by: '<S2>/Discrete Filter'

  0.0,

  // Expression: 0
  //  Referenced by: '<S2>/Constant1'

  0.0,

  // Expression: 0
  //  Referenced by: '<S2>/Switch'

  0.0,

  // Computed Parameter: u_Y0_c
  //  Referenced by: '<S3>/u'

  0.0,

  // Computed Parameter: yhat_Y0_f
  //  Referenced by: '<S3>/yhat'

  0.0,

  // Computed Parameter: prms_Y0_l
  //  Referenced by: '<S3>/prms'

  0.0,

  // Computed Parameter: uref_Y0_k
  //  Referenced by: '<S3>/uref'

  0.0,

  // Computed Parameter: prmErr_Y0_k
  //  Referenced by: '<S3>/prmErr'

  0.0,

  // Computed Parameter: uclean_Y0_e
  //  Referenced by: '<S3>/uclean'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S87>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S91>/LastPcov'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S87>/y.wt_zero'

  { 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/u.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S91>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S91>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G1.A
  //  Referenced by: '<S89>/Constant3'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: 0
  //  Referenced by: '<S89>/Unit Delay2'

  0.0,

  // Expression: 0
  //  Referenced by: '<S89>/Unit Delay3'

  0.0,

  // Expression: 1e4*eye(3)
  //  Referenced by: '<S89>/Constant4'

  { 10000.0, 0.0, 0.0, 0.0, 10000.0, 0.0, 0.0, 0.0, 10000.0 },

  // Expression: 1e4*eye(3)
  //  Referenced by: '<S89>/Constant5'

  { 10000.0, 0.0, 0.0, 0.0, 10000.0, 0.0, 0.0, 0.0, 10000.0 },

  // Expression: G1.C
  //  Referenced by: '<S89>/Constant12'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: G1.D
  //  Referenced by: '<S89>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: 1
  //  Referenced by: '<S89>/Constant11'

  1.0,

  // Expression: [0;0]
  //  Referenced by: '<S89>/Constant2'

  { 0.0, 0.0 },

  // Expression: [0;0]
  //  Referenced by: '<S90>/Constant'

  { 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S126>/X0'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S91>/ym_zero'

  { 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S87>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/umax_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S87>/ymin_zero'

  { 0.0, 0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S87>/ymax_zero'

  { 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S87>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,2)
  //  Referenced by: '<S87>/F_zero'

  { 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/ymin_scale1'

  { 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S87>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S87>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S87>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S91>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S87>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S126>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: 1
  //  Referenced by: '<S3>/Constant'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S126>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S126>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: MVscale
  //  Referenced by: '<S91>/u_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: 0
  //  Referenced by: '<S3>/Constant2'

  0.0,

  // Expression: [0 0 0]
  //  Referenced by: '<S3>/excitation'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: excitation_StdDev_k
  //  Referenced by: '<S3>/excitation'

  { 1.0, 1.0, 1.0 },

  // Expression: [12345 12346 12347]
  //  Referenced by: '<S3>/excitation'

  { 12345.0, 12346.0, 12347.0 },

  // Expression: 0
  //  Referenced by: '<S3>/Discrete Filter'

  0.0,

  // Expression: 0
  //  Referenced by: '<S3>/Constant3'

  0.0,

  // Expression: 0
  //  Referenced by: '<S3>/Switch'

  0.0,

  // Computed Parameter: u_Y0_h
  //  Referenced by: '<S4>/u'

  0.0,

  // Computed Parameter: yhat_Y0_i
  //  Referenced by: '<S4>/yhat'

  0.0,

  // Computed Parameter: prms_Y0_e
  //  Referenced by: '<S4>/prms'

  0.0,

  // Computed Parameter: uref_Y0_km
  //  Referenced by: '<S4>/uref'

  0.0,

  // Computed Parameter: prmErr_Y0_c
  //  Referenced by: '<S4>/prmErr'

  0.0,

  // Computed Parameter: uclean_Y0_k
  //  Referenced by: '<S4>/uclean'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S172>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S176>/LastPcov'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S172>/y.wt_zero'

  { 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/u.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S176>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S176>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G2.A
  //  Referenced by: '<S174>/Constant3'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: 0
  //  Referenced by: '<S174>/Unit Delay2'

  0.0,

  // Expression: 0
  //  Referenced by: '<S174>/Unit Delay3'

  0.0,

  // Expression: 1e4*eye(3)
  //  Referenced by: '<S174>/Constant4'

  { 10000.0, 0.0, 0.0, 0.0, 10000.0, 0.0, 0.0, 0.0, 10000.0 },

  // Expression: 1e4*eye(3)
  //  Referenced by: '<S174>/Constant5'

  { 10000.0, 0.0, 0.0, 0.0, 10000.0, 0.0, 0.0, 0.0, 10000.0 },

  // Expression: G2.C
  //  Referenced by: '<S174>/Constant12'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: G2.D
  //  Referenced by: '<S174>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: 1
  //  Referenced by: '<S174>/Constant11'

  1.0,

  // Expression: [0;0]
  //  Referenced by: '<S174>/Constant2'

  { 0.0, 0.0 },

  // Expression: [0;0]
  //  Referenced by: '<S175>/Constant'

  { 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S211>/X0'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S176>/ym_zero'

  { 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S172>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/umax_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S172>/ymin_zero'

  { 0.0, 0.0 },

  // Expression: zeros(2,1)
  //  Referenced by: '<S172>/ymax_zero'

  { 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S172>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S176>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,2)
  //  Referenced by: '<S172>/F_zero'

  { 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S176>/ymin_scale1'

  { 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S172>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S176>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S172>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S172>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S176>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S172>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S211>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: 1
  //  Referenced by: '<S4>/Constant'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S211>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S211>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: MVscale
  //  Referenced by: '<S176>/u_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: 0
  //  Referenced by: '<S4>/Constant2'

  0.0,

  // Expression: [0 0 0]
  //  Referenced by: '<S4>/excitation'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: excitation_StdDev_h
  //  Referenced by: '<S4>/excitation'

  { 1.0, 1.0, 1.0 },

  // Expression: [12345 12346 12347]
  //  Referenced by: '<S4>/excitation'

  { 12345.0, 12346.0, 12347.0 },

  // Expression: 0
  //  Referenced by: '<S4>/Discrete Filter'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Constant3'

  0.0,

  // Expression: 0
  //  Referenced by: '<S4>/Switch'

  0.0,

  // Expression: Ndis
  //  Referenced by: '<S37>/FixedHorizonOptimizer'

  0,

  // Expression: Ndis
  //  Referenced by: '<S119>/FixedHorizonOptimizer'

  0,

  // Expression: Ndis
  //  Referenced by: '<S204>/FixedHorizonOptimizer'

  0,

  // Expression: iA
  //  Referenced by: '<S9>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S82>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S91>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S167>/isSqrtUsed'

  false,

  // Expression: true
  //  Referenced by: '<S3>/Constant1'

  true,

  // Expression: iA
  //  Referenced by: '<S176>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S252>/isSqrtUsed'

  false,

  // Expression: true
  //  Referenced by: '<S4>/Constant1'

  true,

  // Start of '<S230>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S254>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  ,

  // End of '<S230>/MeasurementUpdate'

  // Start of '<S145>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S169>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  // End of '<S145>/MeasurementUpdate'
};

//
// File trailer for generated code.
//
// [EOF]
//
