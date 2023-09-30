//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController_data.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2514
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sat Sep 30 09:52:54 2023
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
  //  Referenced by: '<S12>/MATLAB Function'

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
  //  Referenced by: '<S90>/MATLAB Function'

  { 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Aod2
  //  Referenced by: '<S160>/MATLAB Function'

  { 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Aod3
  //  Referenced by: '<S230>/MATLAB Function'

  { 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Bod
  //  Referenced by: '<S12>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.05, 0.00062500000000000012 },

  // Variable: Bod1
  //  Referenced by: '<S90>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012 },

  // Variable: Bod2
  //  Referenced by: '<S160>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012 },

  // Variable: Bod3
  //  Referenced by: '<S230>/MATLAB Function'

  { 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012 },

  // Variable: Cod
  //  Referenced by: '<S12>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5 },

  // Variable: Cod1
  //  Referenced by: '<S90>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0 },

  // Variable: Cod2
  //  Referenced by: '<S160>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0 },

  // Variable: Cod3
  //  Referenced by: '<S230>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0 },

  // Variable: Dmn
  //  Referenced by: '<S12>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dmn1
  //  Referenced by:
  //    '<S90>/MATLAB Function'
  //    '<S160>/MATLAB Function'
  //    '<S230>/MATLAB Function'

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },

  // Variable: Dod
  //  Referenced by: '<S12>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod1
  //  Referenced by: '<S90>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod2
  //  Referenced by: '<S160>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: Dod3
  //  Referenced by: '<S230>/MATLAB Function'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Variable: beta
  //  Referenced by:
  //    '<Root>/SupervisoryController'
  //    '<S2>/Gain1'
  //    '<S3>/Gain1'
  //    '<S4>/Gain1'
  //    '<S5>/Gain1'
  //    '<S6>/Gain1'
  //    '<S9>/Gain'

  0.13534,

  // Variable: dt
  //  Referenced by:
  //    '<Root>/SupervisoryController'
  //    '<S2>/MATLAB Function2'
  //    '<S3>/Gain2'
  //    '<S4>/Gain2'
  //    '<S5>/Gain2'
  //    '<S6>/Gain2'
  //    '<S9>/MATLAB Function'
  //    '<S322>/MATLAB Function1'
  //    '<S326>/MATLAB Function1'

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
  //    '<S322>/MATLAB Function1'
  //    '<S326>/MATLAB Function1'

  1.0,

  // Variable: uwt0
  //  Referenced by: '<S9>/Delay1'

  { 0.0, 0.0, 0.0 },

  // Variable: ywt0
  //  Referenced by: '<S9>/Delay'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: 0
  //  Referenced by: '<S86>/L*(y[k]-yhat[k|k-1])'

  0.0,

  // Computed Parameter: u_Y0
  //  Referenced by: '<S2>/u'

  0.0,

  // Computed Parameter: yhat_Y0
  //  Referenced by: '<S2>/yhat'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S10>/G_zero'

  0.0,

  // Expression: lastPcov
  //  Referenced by: '<S13>/LastPcov'

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
  //  Referenced by: '<S10>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S10>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S13>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S13>/last_mv'

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
  //  Referenced by: '<S12>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S43>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S13>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S10>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S10>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: -2
  //  Referenced by: '<S2>/Gain2'

  -2.0,

  // Expression: 2
  //  Referenced by: '<S2>/Gain3'

  2.0,

  // Expression: zeros(1,3)
  //  Referenced by: '<S10>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S13>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,6)
  //  Referenced by: '<S10>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S13>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S10>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S13>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S10>/switch_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S10>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S13>/uref_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S10>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S43>/P0'

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
  //  Referenced by: '<S43>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S43>/G'

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
  //  Referenced by: '<S13>/u_scale'

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
  //  Referenced by: '<S156>/L*(y[k]-yhat[k|k-1])'

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
  //  Referenced by: '<S89>/G_zero'

  0.0,

  // Expression: zeros(4,1)
  //  Referenced by: '<S89>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S89>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S89>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S91>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S89>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S91>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S91>/last_mv'

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
  //  Referenced by: '<S90>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S113>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S91>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S89>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(4,1)
  //  Referenced by: '<S89>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(4,1)
  //  Referenced by: '<S89>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S89>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,4)
  //  Referenced by: '<S89>/F_zero'

  { 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S91>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S89>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S113>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S3>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S113>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S113>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S91>/umin_scale1'

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
  //  Referenced by: '<S159>/G_zero'

  0.0,

  // Expression: zeros(5,1)
  //  Referenced by: '<S159>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S159>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S159>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S161>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S159>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S161>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S161>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G1.A
  //  Referenced by: '<S4>/Constant3'

  { 1.0001226651251065, -0.0016804698402870472, -6.4394278387945068E-5,
    0.999403766481181 },

  // Expression: G1.B
  //  Referenced by: '<S4>/Constant4'

  { -0.0071318388640683669, -0.012368697539837475, 0.13659152296751884,
    -0.063538603002101665, -0.04588890657825119, 0.10453377294976968 },

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
  //  Referenced by: '<S160>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S183>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S161>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S159>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S159>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S159>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S159>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S159>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S161>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,5)
  //  Referenced by: '<S159>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S161>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S159>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S161>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S159>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S159>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S183>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S4>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S183>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S183>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S161>/umin_scale1'

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
  //  Referenced by: '<S229>/G_zero'

  0.0,

  // Expression: zeros(5,1)
  //  Referenced by: '<S229>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S229>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S229>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S231>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S229>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S231>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S231>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: G2.A
  //  Referenced by: '<S5>/Constant3'

  { 1.0, 0.0, 0.0, 1.0 },

  // Expression: G2.B
  //  Referenced by: '<S5>/Constant4'

  { 0.00025, -0.099108107281090041, -0.00025, 0.020123801605080698,
    -0.019490220104579719, -0.000125 },

  // Expression: G2.C
  //  Referenced by: '<S5>/Constant12'

  { 1.0, 0.0, 0.0, 0.0, 2.0, 0.0 },

  // Expression: G2.D
  //  Referenced by: '<S5>/Constant13'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(size(G2.A, 1), 1)
  //  Referenced by: '<S5>/Constant2'

  { 0.0, 0.0 },

  // Expression: zeros(size(God3.A,1),1)
  //  Referenced by: '<S230>/Constant1'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.X0
  //  Referenced by: '<S253>/X0'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(nym,1)
  //  Referenced by: '<S231>/ym_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S229>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S229>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S229>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S229>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S229>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S231>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,5)
  //  Referenced by: '<S229>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S231>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S229>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S231>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S229>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S229>/ecr.wt_zero'

  0.0,

  // Expression: pInitialization.P0
  //  Referenced by: '<S253>/P0'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: 1
  //  Referenced by: '<S5>/Constant1'

  1.0,

  // Expression: pInitialization.H
  //  Referenced by: '<S253>/H'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: pInitialization.G
  //  Referenced by: '<S253>/G'

  { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 },

  // Expression: MVscale
  //  Referenced by: '<S231>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: 1000
  //  Referenced by: '<S5>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S5>/Saturation'

  0.0,

  // Computed Parameter: ywt_Y0
  //  Referenced by: '<S9>/ywt'

  0.0,

  // Computed Parameter: y_Y0
  //  Referenced by: '<S9>/y_'

  0.0,

  // Computed Parameter: r_Y0
  //  Referenced by: '<S9>/r_'

  0.0,

  // Computed Parameter: u_Y0_j
  //  Referenced by: '<S6>/u'

  0.0,

  // Computed Parameter: yhat_Y0_b
  //  Referenced by: '<S6>/yhat'

  0.0,

  // Expression: [0;0]
  //  Referenced by: '<S6>/Constant'

  { 0.0, 0.0 },

  // Computed Parameter: DiscreteTimeIntegrator_gainv_kx
  //  Referenced by: '<S6>/Discrete-Time Integrator'

  1.0,

  // Expression: zeros(2, 1)
  //  Referenced by: '<S6>/Discrete-Time Integrator'

  { 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S299>/G_zero'

  0.0,

  // Expression: zeros(5,1)
  //  Referenced by: '<S299>/y.wt_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S299>/du.wt_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S299>/ext.mv_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S300>/ext.mv_scale'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(3,1)
  //  Referenced by: '<S299>/mv.target_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: RMVscale
  //  Referenced by: '<S300>/ext.mv_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: lastu+uoff
  //  Referenced by: '<S300>/last_mv'

  { 0.0, 0.0, 0.0 },

  // Expression: lastx+xoff
  //  Referenced by: '<S300>/last_x'

  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S299>/md_zero'

  0.0,

  // Expression: zeros(3,1)
  //  Referenced by: '<S299>/umin_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S299>/ymin_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(5,1)
  //  Referenced by: '<S299>/ymax_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: zeros(1,3)
  //  Referenced by: '<S299>/E_zero'

  { 0.0, 0.0, 0.0 },

  // Expression: MVscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S300>/umin_scale4'

  { 1.0, 1.0, 1.0 },

  // Expression: zeros(1,5)
  //  Referenced by: '<S299>/F_zero'

  { 0.0, 0.0, 0.0, 0.0, 0.0 },

  // Expression: Yscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S300>/ymin_scale1'

  { 1.0, 1.0, 1.0, 1.0, 1.0 },

  // Expression: zeros(1,1)
  //  Referenced by: '<S299>/S_zero'

  0.0,

  // Expression: MDscale(:,ones(1,max(nCC,1)))'
  //  Referenced by: '<S300>/ymin_scale2'

  1.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S299>/switch_zero'

  0.0,

  // Expression: zeros(1,1)
  //  Referenced by: '<S299>/ecr.wt_zero'

  0.0,

  // Expression: MVscale
  //  Referenced by: '<S300>/umin_scale1'

  { 1.0, 1.0, 1.0 },

  // Expression: 1000
  //  Referenced by: '<S6>/Saturation'

  1000.0,

  // Expression: 0
  //  Referenced by: '<S6>/Saturation'

  0.0,

  // Expression: Ndis
  //  Referenced by: '<S41>/FixedHorizonOptimizer'

  0,

  // Expression: iA
  //  Referenced by: '<S13>/Memory'

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
  //  Referenced by: '<S84>/isSqrtUsed'

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
  //  Referenced by: '<S154>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S161>/Memory'

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
  //  Referenced by: '<S224>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S231>/Memory'

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
  //  Referenced by: '<S294>/isSqrtUsed'

  false,

  // Expression: iA
  //  Referenced by: '<S300>/Memory'

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

  // Start of '<S272>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S296>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  ,

  // End of '<S272>/MeasurementUpdate'

  // Start of '<S202>/MeasurementUpdate'
  {
    // Expression: 0
    //  Referenced by: '<S226>/L*(y[k]-yhat[k|k-1])'

    0.0
  }
  ,

  // End of '<S202>/MeasurementUpdate'

  // Start of '<S1>/paramEst2'
  {
    // Computed Parameter: theta_Y0
    //  Referenced by: '<S8>/theta'

    0.0,

    // Computed Parameter: P_Y0
    //  Referenced by: '<S8>/P'

    0.0,

    // Computed Parameter: err_Y0
    //  Referenced by: '<S8>/err'

    0.0,

    // Expression: 0
    //  Referenced by: '<S326>/Unit Delay3'

    0.0
  }
  ,

  // End of '<S1>/paramEst2'

  // Start of '<S1>/paramEst1'
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
    //  Referenced by: '<S322>/Unit Delay3'

    0.0
  }
  // End of '<S1>/paramEst1'
};

//
// File trailer for generated code.
//
// [EOF]
//
