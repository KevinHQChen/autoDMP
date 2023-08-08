//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2467
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon Aug  7 22:26:20 2023
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
  //  Referenced by: '<S11>/MATLAB Function'

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

  // Variable: Aod1
  //  Referenced by: '<S89>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Aod2
  //  Referenced by: '<S159>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Aod3
  //  Referenced by: '<S229>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Bod
  //  Referenced by: '<S11>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.05, 0.00062500000000000012 },

  // Variable: Bod1
  //  Referenced by: '<S89>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.025 },

  // Variable: Bod2
  //  Referenced by: '<S159>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.025 },

  // Variable: Bod3
  //  Referenced by: '<S229>/MATLAB Function'

  { 0.025, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.025 },

  // Variable: Cod
  //  Referenced by: '<S11>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5 },

  // Variable: Cod1
  //  Referenced by: '<S89>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Cod2
  //  Referenced by: '<S159>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Cod3
  //  Referenced by: '<S229>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn
  //  Referenced by: '<S11>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn1
  //  Referenced by:
  //    '<S89>/MATLAB Function'
  //    '<S159>/MATLAB Function'
  //    '<S229>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dod
  //  Referenced by: '<S11>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod1
  //  Referenced by: '<S89>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod2
  //  Referenced by: '<S159>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod3
  //  Referenced by: '<S229>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: beta
  //  Referenced by:
  //    '<Root>/SupervisoryController'
  //    '<S2>/Gain1'
  //    '<S8>/Gain'

  0.13534,

  // Variable: dt
  //  Referenced by:
  //    '<Root>/SupervisoryController'
  //    '<S2>/MATLAB Function2'
  //    '<S3>/Gain2'
  //    '<S4>/Gain2'
  //    '<S5>/Gain2'
  //    '<S8>/MATLAB Function'
  //    '<S298>/MATLAB Function1'
  //    '<S302>/MATLAB Function1'

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
  //    '<S298>/MATLAB Function1'
  //    '<S302>/MATLAB Function1'

  1.0,

  // Variable: uwt0
  //  Referenced by: '<S8>/Delay1'

  { 0.0, 0.0, 0.0 },

  // Variable: ywt0
  //  Referenced by: '<S8>/Delay'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: 0
  //  Referenced by: '<S85>/L*(y[k]-yhat[k|k-1])'

  0.0,

  // Computed Parameter: u_Y0
  //  Referenced by: '<S2>/u'

  0.0,

  // Computed Parameter: yhat_Y0
  //  Referenced by: '<S2>/yhat'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S9>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S12>/LastPcov'

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
  //  Referenced by: '<S9>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S9>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S12>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S12>/last_mv'

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
  //  Referenced by: '<S11>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S42>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S12>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S9>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S9>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: -2
  //  Referenced by: '<S2>/Gain2'

  -2.0,

  // Expression: 2
  //  Referenced by: '<S2>/Gain3'

  2.0,

  // Expression: zeros(1,3)
  //  Referenced by: '<S9>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S12>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S9>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S12>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S9>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S12>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S9>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S9>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S12>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S9>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S42>/P0'

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
  //  Referenced by: '<S42>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S42>/G'

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
  //  Referenced by: '<S12>/u_scale'

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

  // Expression: 0
  //  Referenced by: '<S155>/L*(y[k]-yhat[k|k-1])'

  0.0,

  // Computed Parameter: u_Y0_b
  //  Referenced by: '<S3>/u'

  0.0,

  // Computed Parameter: yhat_Y0_d
  //  Referenced by: '<S3>/yhat'

  0.0,

  // Expression: 0
  //  Referenced by: '<S3>/Constant'

  0.0,

  // Computed Parameter: DiscreteTimeIntegrator_gainval
  //  Referenced by: '<S3>/Discrete-Time Integrator'

  1.0,

  // Expression: 0
  //  Referenced by: '<S3>/Discrete-Time Integrator'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/G_zero'

  0.0,

  // Expression: zeros(4,1)
  //  Referenced by: '<S88>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/u.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S90>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S90>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S90>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G0.B
  //  Referenced by: '<S3>/Constant4'

  { 0.22599754912864506, -0.087367622049086546, -0.087367622049086546 },

  // Expression: G0.C
  //  Referenced by: '<S3>/Constant12'

  { 1.0, 0.0, 0.0 },

  // Expression: G0.D
  //  Referenced by: '<S3>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: G0.A
  //  Referenced by: '<S3>/Constant3'

  0.9977021834978721,

  // Expression: zeros(size(G0.A, 1), 1)
  //  Referenced by: '<S3>/Constant2'

  0.0,

  // Expression: zeros(size(God1.A,1),1)
  //  Referenced by: '<S89>/Constant1'

  { 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S112>/X0'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S90>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S88>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(4,1)
  //  Referenced by: '<S88>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(4,1)
  //  Referenced by: '<S88>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S88>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S90>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,4)
  //  Referenced by: '<S88>/F_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S90>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S90>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S88>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S112>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: 1
  //  Referenced by: '<S3>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S112>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S112>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    1.0 },

  // Expression: MVscale
  //  Referenced by: '<S90>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: 1000
  //  Referenced by: '<S3>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S3>/Saturation'

  0.0,

  // Computed Parameter: u_Y0_n
  //  Referenced by: '<S4>/u'

  0.0,

  // Computed Parameter: yhat_Y0_k
  //  Referenced by: '<S4>/yhat'

  0.0,

  // Expression: [0;0]
  //  Referenced by: '<S4>/Constant'

  { 0.0, 0.0 },

  // Computed Parameter: DiscreteTimeIntegrator_gainva_b
  //  Referenced by: '<S4>/Discrete-Time Integrator'

  1.0,

  // Expression: zeros(2, 1)
  //  Referenced by: '<S4>/Discrete-Time Integrator'

  { 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S158>/G_zero'

  0.0,

  // Expression: zeros(5,1)
  //  Referenced by: '<S158>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S158>/u.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S158>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S158>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S160>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S158>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S160>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S160>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G1.A
  //  Referenced by: '<S4>/Constant3'

  { 1.0019018232567545, -0.0016184163979325317, 0.0010492809235588263,
    0.99949680855411027 },

  // Expression: G1.B
  //  Referenced by: '<S4>/Constant4'

  { 0.0402193096070776, 0.010618210823094141, 0.1047915945578364,
    -0.11013071843533731, -0.029336890171075523, 0.09778708813999859 },

  // Expression: G1.C
  //  Referenced by: '<S4>/Constant12'

  { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: G1.D
  //  Referenced by: '<S4>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(size(G1.A, 1), 1)
  //  Referenced by: '<S4>/Constant2'

  { 0.0, 0.0 },

  // Expression: zeros(size(God2.A,1),1)
  //  Referenced by: '<S159>/Constant1'

  { 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S182>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S160>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S158>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S158>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S158>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S158>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S158>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S160>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,5)
  //  Referenced by: '<S158>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S160>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S158>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S160>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S158>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S158>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S182>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S4>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S182>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S182>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S160>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: 1000
  //  Referenced by: '<S4>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S4>/Saturation'

  0.0,

  // Computed Parameter: u_Y0_h
  //  Referenced by: '<S5>/u'

  0.0,

  // Computed Parameter: yhat_Y0_f
  //  Referenced by: '<S5>/yhat'

  0.0,

  // Expression: [0;0]
  //  Referenced by: '<S5>/Constant'

  { 0.0, 0.0 },

  // Computed Parameter: DiscreteTimeIntegrator_gainva_k
  //  Referenced by: '<S5>/Discrete-Time Integrator'

  1.0,

  // Expression: zeros(2, 1)
  //  Referenced by: '<S5>/Discrete-Time Integrator'

  { 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S228>/G_zero'

  0.0,

  // Expression: zeros(5,1)
  //  Referenced by: '<S228>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S228>/u.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S228>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S228>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S230>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S228>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S230>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S230>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G2.A
  //  Referenced by: '<S5>/Constant3'

  { 0.98979340248168846, -0.036210018833518431, 0.00030104618687214896,
    1.0007224384072255 },

  // Expression: G2.B
  //  Referenced by: '<S5>/Constant4'

  { 0.025118718490858984, 0.016819265413734148, -0.0189087824026894,
    0.066906613140151408, 0.032190816623613963, -0.026334884546125613 },

  // Expression: G2.C
  //  Referenced by: '<S5>/Constant12'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 },

  // Expression: G2.D
  //  Referenced by: '<S5>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(size(G2.A, 1), 1)
  //  Referenced by: '<S5>/Constant2'

  { 0.0, 0.0 },

  // Expression: zeros(size(God3.A,1),1)
  //  Referenced by: '<S229>/Constant1'

  { 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S252>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S230>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S228>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S228>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S228>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S228>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S228>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S230>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,5)
  //  Referenced by: '<S228>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S230>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S228>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S230>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S228>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S228>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S252>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S5>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S252>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S252>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S230>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: 1000
  //  Referenced by: '<S5>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S5>/Saturation'

  0.0,

  // Computed Parameter: ywt_Y0
  //  Referenced by: '<S8>/ywt'

  0.0,

  // Computed Parameter: y_Y0
  //  Referenced by: '<S8>/y_'

  0.0,

  // Computed Parameter: r_Y0
  //  Referenced by: '<S8>/r_'

  0.0,

  // Expression: Ndis
  //  Referenced by: '<S40>/FixedHorizonOptimizer'

  0,

  // Expression: iA
  //  Referenced by: '<S12>/Memory'

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
  //  Referenced by: '<S83>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S90>/Memory'

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
    false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S153>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S160>/Memory'

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
    false, false, false, false, false, false, false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S223>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S230>/Memory'

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
    false, false, false, false, false },

  // Expression: pInitialization.isSqrtUsed
  //  Referenced by: '<S293>/isSqrtUsed'

  false,

  // Start of '<S271>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S295>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  ,

  // End of '<S271>/MeasurementUpdate'

  // Start of '<S201>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S225>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  ,

  // End of '<S201>/MeasurementUpdate'

  // Start of '<S1>/paramEst2'
  {
    // Computed Parameter: theta_Y0
    //  Referenced by: '<S7>/theta'

    0.0,

    // Computed Parameter: P_Y0
    //  Referenced by: '<S7>/P'

    0.0,

    // Computed Parameter: err_Y0
    //  Referenced by: '<S7>/err'

    0.0,

    // Expression: 0
    //  Referenced by: '<S302>/Unit Delay3'

    0.0
  }
  ,

  // End of '<S1>/paramEst2'

  // Start of '<S1>/paramEst1'
  {
    // Computed Parameter: theta_Y0
    //  Referenced by: '<S6>/theta'

    0.0,

    // Computed Parameter: P_Y0
    //  Referenced by: '<S6>/P'

    0.0,

    // Computed Parameter: err_Y0
    //  Referenced by: '<S6>/err'

    0.0,

    // Expression: 0
    //  Referenced by: '<S298>/Unit Delay3'

    0.0
  }
  // End of '<S1>/paramEst1'
};

//
// File trailer for generated code.
//
// [EOF]
//
