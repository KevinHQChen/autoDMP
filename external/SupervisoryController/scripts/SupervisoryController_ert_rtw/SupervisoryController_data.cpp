//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2337
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Fri Aug  4 06:09:14 2023
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
    { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    0.0,
    0.0,
    0.0
  },

  // Variable: Aod
  //  Referenced by: '<S9>/MATLAB Function'

  { 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Bod
  //  Referenced by: '<S9>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.05, 0.00062500000000000012 },

  // Variable: Cod
  //  Referenced by: '<S9>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5 },

  // Variable: Dmn
  //  Referenced by: '<S9>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dod
  //  Referenced by: '<S9>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: beta
  //  Referenced by:
  //    '<S2>/Gain'
  //    '<S2>/Gain1'
  //    '<S3>/Gain'
  //    '<S3>/Gain1'

  0.13534,

  // Variable: dt
  //  Referenced by:
  //    '<Root>/SupervisoryController'
  //    '<S2>/MATLAB Function'
  //    '<S2>/MATLAB Function2'
  //    '<S3>/MATLAB Function'
  //    '<S158>/MATLAB Function1'
  //    '<S162>/MATLAB Function1'

  0.025,

  // Variable: lpfDen
  //  Referenced by: '<S2>/Discrete Filter1'

  1.0,

  // Variable: lpfNum
  //  Referenced by: '<S2>/Discrete Filter1'

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

  // Variable: mdlNum
  //  Referenced by:
  //    '<S2>/MATLAB Function2'
  //    '<S158>/MATLAB Function1'
  //    '<S162>/MATLAB Function1'

  1.0,

  // Variable: uwt0
  //  Referenced by:
  //    '<S2>/Delay1'
  //    '<S3>/Delay1'

  { 0.0, 0.0, 0.0 },

  // Variable: ywt0
  //  Referenced by:
  //    '<S2>/Delay'
  //    '<S3>/Delay'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: 0
  //  Referenced by: '<S83>/L*(y[k]-yhat[k|k-1])'

  0.0,

  // Computed Parameter: u_Y0
  //  Referenced by: '<S2>/u'

  0.0,

  // Computed Parameter: ywt_Y0
  //  Referenced by: '<S2>/ywt'

  0.0,

  // Computed Parameter: yhat_Y0
  //  Referenced by: '<S2>/yhat'

  0.0,

  // Computed Parameter: r_Y0
  //  Referenced by: '<S2>/r_'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S6>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S10>/LastPcov'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S6>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S6>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S10>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S10>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G.C
  //  Referenced by: '<S2>/Constant12'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: G.D
  //  Referenced by: '<S2>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: zeros(2*ns, 1)
  //  Referenced by: '<S2>/Constant2'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(size(God.A,1),1)
  //  Referenced by: '<S9>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S40>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S10>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S6>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S6>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: -2
  //  Referenced by: '<S2>/Gain2'

  -2.0,

  // Expression: 2
  //  Referenced by: '<S2>/Gain3'

  2.0,

  // Expression: zeros(1,3)
  //  Referenced by: '<S6>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S10>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S6>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S10>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S6>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S10>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S6>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S6>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S10>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S6>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S40>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S2>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S40>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S40>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S10>/u_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: [0 0 0]
  //  Referenced by: '<S2>/excitation'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: excitation_StdDev
  //  Referenced by: '<S2>/excitation'

  { 1.0, 1.0, 1.0 },

  // Expression: [12345 12346 12347]
  //  Referenced by: '<S2>/excitation'

  { 12345.0, 12346.0, 12347.0 },

  // Expression: 0
  //  Referenced by: '<S2>/Discrete Filter1'

  0.0,

  // Expression: 1000
  //  Referenced by: '<S2>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S2>/Saturation'

  0.0,

  // Expression: MVscale
  //  Referenced by: '<S92>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S114>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S136>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S92>/constant'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale
  //  Referenced by: '<S92>/umin_scale2'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S114>/constant'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale
  //  Referenced by: '<S114>/umin_scale2'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S136>/constant'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale
  //  Referenced by: '<S136>/umin_scale2'

  { 1.0, 1.0, 1.0 },

  // Expression: mmpc_u0
  //  Referenced by: '<S88>/u0_zero'

  { 0.0, 0.0, 0.0 },

  // Computed Parameter: u_Y0_j
  //  Referenced by: '<S3>/u'

  0.0,

  // Computed Parameter: ywt_Y0_j
  //  Referenced by: '<S3>/ywt'

  0.0,

  // Computed Parameter: r_Y0_g
  //  Referenced by: '<S3>/r_'

  0.0,

  // Expression: -2
  //  Referenced by: '<S3>/Gain2'

  -2.0,

  // Expression: 2
  //  Referenced by: '<S3>/Gain3'

  2.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/G_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: 0.0
  //  Referenced by: '<S3>/Delay2'

  0.0,

  // Expression: RMVscale
  //  Referenced by: '<S92>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S92>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S92>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: lastx+xoff
  //  Referenced by: '<S92>/last_x'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S89>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S92>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S89>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S92>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S92>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/ecr.wt_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S90>/G_zero'

  0.0,

  // Expression: RMVscale
  //  Referenced by: '<S114>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: RMVscale
  //  Referenced by: '<S114>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S114>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: lastx+xoff
  //  Referenced by: '<S114>/last_x'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S90>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S114>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S90>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S114>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S90>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S114>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S91>/G_zero'

  0.0,

  // Expression: RMVscale
  //  Referenced by: '<S136>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: RMVscale
  //  Referenced by: '<S136>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S136>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: lastx+xoff
  //  Referenced by: '<S136>/last_x'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S91>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S136>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S91>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S136>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S91>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S136>/ymin_scale2'

  1.0,

  // Expression: 1000
  //  Referenced by: '<S3>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S3>/Saturation'

  0.0,

  // Expression: Ndis
  //  Referenced by: '<S38>/FixedHorizonOptimizer'

  0,

  // Expression: iA
  //  Referenced by: '<S10>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S81>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S92>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false },

  // Expression: iA
  //  Referenced by: '<S114>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false },

  // Expression: iA
  //  Referenced by: '<S136>/Memory'

  { false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false },

  // Start of '<S1>/paramEst2'
  {
    // Computed Parameter: theta_Y0
    //  Referenced by: '<S5>/theta'

    0.0,

    // Computed Parameter: P_Y0
    //  Referenced by: '<S5>/P'

    0.0,

    // Computed Parameter: err_Y0
    //  Referenced by: '<S5>/err'

    0.0,

    // Expression: 0
    //  Referenced by: '<S162>/Unit Delay3'

    0.0
  }
  ,

  // End of '<S1>/paramEst2'

  // Start of '<S1>/paramEst1'
  {
    // Computed Parameter: theta_Y0
    //  Referenced by: '<S4>/theta'

    0.0,

    // Computed Parameter: P_Y0
    //  Referenced by: '<S4>/P'

    0.0,

    // Computed Parameter: err_Y0
    //  Referenced by: '<S4>/err'

    0.0,

    // Expression: 0
    //  Referenced by: '<S158>/Unit Delay3'

    0.0
  }
  // End of '<S1>/paramEst1'
};

//
// File trailer for generated code.
//
// [EOF]
//
