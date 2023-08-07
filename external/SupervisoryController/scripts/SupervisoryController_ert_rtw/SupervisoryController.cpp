//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2463
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon Aug  7 01:52:35 2023
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives:
//    1. Execution efficiency
//    2. RAM efficiency
// Validation result: Not run
//
#include "SupervisoryController.h"
#include "omp.h"
#include "rtwtypes.h"
#include <cstring>
#include <emmintrin.h>
#include <cmath>
#include "SupervisoryController_capi.h"
#include "zero_crossing_types.h"
#include <stddef.h>
#include "solver_zc.h"

// Named constants for MATLAB Function: '<S40>/FixedHorizonOptimizer'
const int32_T degrees{ 4 };

const int32_T nu{ 3 };

// Named constants for Chart: '<Root>/SupervisoryController'
const uint8_T IN_HandleEvent{ 1U };

const uint8_T IN_RequestEvent{ 2U };

#define NumBitsPerChar                 8U
#include "solver_zc.h"
#ifndef slZcHadEvent
#define slZcHadEvent(ev, zcsDir)       (((ev) & (zcsDir)) != 0x00 )
#endif

#ifndef slZcUnAliasEvents
#define slZcUnAliasEvents(evL, evR)    ((((slZcHadEvent((evL), (SL_ZCS_EVENT_N2Z)) && slZcHadEvent((evR), (SL_ZCS_EVENT_Z2P))) || (slZcHadEvent((evL), (SL_ZCS_EVENT_P2Z)) && slZcHadEvent((evR), (SL_ZCS_EVENT_Z2N)))) ? (SL_ZCS_EVENT_NUL) : (evR)))
#endif

extern real_T rt_powd_snf(real_T u0, real_T u1);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u);
static int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator);
extern "C"
{
  static ZCEventType rt_ZCFcn(ZCDirection zcDir, ZCSigState *prevZc, real_T
    currValue);
}                                      // extern "C"

extern "C"
{
  real_T rtInf;
  real_T rtMinusInf;
  real_T rtNaN;
  real32_T rtInfF;
  real32_T rtMinusInfF;
  real32_T rtNaNF;
}

extern "C"
{
  //
  // Initialize rtInf needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetInf(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T inf{ 0.0 };

    if (bitsPerReal == 32U) {
      inf = rtGetInfF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0x7FF00000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      inf = tmpVal.fltVal;
    }

    return inf;
  }

  //
  // Initialize rtInfF needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetInfF(void)
  {
    IEEESingle infF;
    infF.wordL.wordLuint = 0x7F800000U;
    return infF.wordL.wordLreal;
  }

  //
  // Initialize rtMinusInf needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetMinusInf(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T minf{ 0.0 };

    if (bitsPerReal == 32U) {
      minf = rtGetMinusInfF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0xFFF00000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      minf = tmpVal.fltVal;
    }

    return minf;
  }

  //
  // Initialize rtMinusInfF needed by the generated code.
  // Inf is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetMinusInfF(void)
  {
    IEEESingle minfF;
    minfF.wordL.wordLuint = 0xFF800000U;
    return minfF.wordL.wordLreal;
  }
}

extern "C"
{
  //
  // Initialize rtNaN needed by the generated code.
  // NaN is initialized as non-signaling. Assumes IEEE.
  //
  static real_T rtGetNaN(void)
  {
    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    real_T nan{ 0.0 };

    if (bitsPerReal == 32U) {
      nan = rtGetNaNF();
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.bitVal.words.wordH = 0xFFF80000U;
      tmpVal.bitVal.words.wordL = 0x00000000U;
      nan = tmpVal.fltVal;
    }

    return nan;
  }

  //
  // Initialize rtNaNF needed by the generated code.
  // NaN is initialized as non-signaling. Assumes IEEE.
  //
  static real32_T rtGetNaNF(void)
  {
    IEEESingle nanF{ { 0.0F } };

    nanF.wordL.wordLuint = 0xFFC00000U;
    return nanF.wordL.wordLreal;
  }
}

extern "C"
{
  // Detect zero crossings events.
  static ZCEventType rt_ZCFcn(ZCDirection zcDir, ZCSigState *prevZc, real_T
    currValue)
  {
    slZcEventType zcsDir;
    slZcEventType tempEv;
    ZCEventType zcEvent{ NO_ZCEVENT }; // assume

    // zcEvent matrix
    static const slZcEventType eventMatrix[4][4]{
      //          ZER              POS              NEG              UNK
      { SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_Z2P, SL_ZCS_EVENT_Z2N, SL_ZCS_EVENT_NUL },// ZER 

      { SL_ZCS_EVENT_P2Z, SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_P2N, SL_ZCS_EVENT_NUL },// POS 

      { SL_ZCS_EVENT_N2Z, SL_ZCS_EVENT_N2P, SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_NUL },// NEG 

      { SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_NUL, SL_ZCS_EVENT_NUL }// UNK 
    };

    // get prevZcEvent and prevZcSign from prevZc
    const slZcEventType prevEv{ (slZcEventType)(((uint8_T)(*prevZc)) >> 2) };

    const slZcSignalSignType prevSign{ (slZcSignalSignType)(((uint8_T)(*prevZc))
      & (uint8_T)0x03) };

    // get current zcSignal sign from current zcSignal value
    const slZcSignalSignType currSign{ (slZcSignalSignType)((currValue) > 0.0 ?
      SL_ZCS_SIGN_POS :
      ((currValue) < 0.0 ? SL_ZCS_SIGN_NEG : SL_ZCS_SIGN_ZERO)) };

    // get current zcEvent based on prev and current zcSignal value
    slZcEventType currEv { eventMatrix[prevSign][currSign] };

    // get slZcEventType from ZCDirection
    switch (zcDir) {
     case ANY_ZERO_CROSSING:
      zcsDir = SL_ZCS_EVENT_ALL;
      break;

     case FALLING_ZERO_CROSSING:
      zcsDir = SL_ZCS_EVENT_ALL_DN;
      break;

     case RISING_ZERO_CROSSING:
      zcsDir = SL_ZCS_EVENT_ALL_UP;
      break;

     default:
      zcsDir = SL_ZCS_EVENT_NUL;
      break;
    }

    //had event, check if double zc happend remove double detection.
    if (slZcHadEvent(currEv, zcsDir)) {
      currEv = (slZcEventType)(slZcUnAliasEvents(prevEv, currEv));
    } else {
      currEv = SL_ZCS_EVENT_NUL;
    }

    // Update prevZc
    tempEv = (slZcEventType)(currEv << 2);// shift left by 2 bits
    *prevZc = (ZCSigState)((currSign) | (tempEv));
    if ((currEv & SL_ZCS_EVENT_ALL_DN) != 0) {
      zcEvent = FALLING_ZCEVENT;
    } else if ((currEv & SL_ZCS_EVENT_ALL_UP) != 0) {
      zcEvent = RISING_ZCEVENT;
    } else {
      zcEvent = NO_ZCEVENT;
    }

    return zcEvent;
  }                                    // rt_ZCFcn
}

extern "C"
{
  //
  // Initialize the rtInf, rtMinusInf, and rtNaN needed by the
  // generated code. NaN is initialized as non-signaling. Assumes IEEE.
  //
  static void rt_InitInfAndNaN(size_t realSize)
  {
    (void) (realSize);
    rtNaN = rtGetNaN();
    rtNaNF = rtGetNaNF();
    rtInf = rtGetInf();
    rtInfF = rtGetInfF();
    rtMinusInf = rtGetMinusInf();
    rtMinusInfF = rtGetMinusInfF();
  }

  // Test if value is infinite
  static boolean_T rtIsInf(real_T value)
  {
    return (boolean_T)((value==rtInf || value==rtMinusInf) ? 1U : 0U);
  }

  // Test if single-precision value is infinite
  static boolean_T rtIsInfF(real32_T value)
  {
    return (boolean_T)(((value)==rtInfF || (value)==rtMinusInfF) ? 1U : 0U);
  }

  // Test if value is not a number
  static boolean_T rtIsNaN(real_T value)
  {
    boolean_T result{ (boolean_T) 0 };

    size_t bitsPerReal{ sizeof(real_T) * (NumBitsPerChar) };

    if (bitsPerReal == 32U) {
      result = rtIsNaNF((real32_T)value);
    } else {
      union {
        LittleEndianIEEEDouble bitVal;
        real_T fltVal;
      } tmpVal;

      tmpVal.fltVal = value;
      result = (boolean_T)((tmpVal.bitVal.words.wordH & 0x7FF00000) ==
                           0x7FF00000 &&
                           ( (tmpVal.bitVal.words.wordH & 0x000FFFFF) != 0 ||
                            (tmpVal.bitVal.words.wordL != 0) ));
    }

    return result;
  }

  // Test if single-precision value is not a number
  static boolean_T rtIsNaNF(real32_T value)
  {
    IEEESingle tmp;
    tmp.wordL.wordLreal = value;
    return (boolean_T)( (tmp.wordL.wordLuint & 0x7F800000) == 0x7F800000 &&
                       (tmp.wordL.wordLuint & 0x007FFFFF) != 0 );
  }
}

static int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator)
{
  return static_cast<int32_T>(((numerator < 0) != (denominator < 0)) &&
    (numerator % denominator != 0) ? -1 : 0) + numerator / denominator;
}

void SupervisoryController::binary_expand_op(real_T in1[36], int32_T in2,
  int32_T in3, int32_T in4, const real_T in5[12], int32_T in6, int32_T in7,
  const real_T in8[12])
{
  int32_T in3_idx_0;
  int32_T stride_0_0;

  // MATLAB Function: '<S298>/MATLAB Function1'
  in3_idx_0 = in3 - in2;
  stride_0_0 = (in7 - in6) + 1 != 1 ? static_cast<int32_T>(1) : static_cast<
    int32_T>(0);
  for (int32_T i{0}; i < in3_idx_0; i++) {
    in1[(in2 + i) + 12 * in4] = in5[i * stride_0_0 + in6] * in8[(in4 << 2UL) + i];
  }

  // End of MATLAB Function: '<S298>/MATLAB Function1'
}

//
// System initialize for function-call system:
//    '<S1>/paramEst1'
//    '<S1>/paramEst2'
//
void SupervisoryController::paramEst1_Init(real_T rty_theta[12], real_T rty_P
  [144], real_T rty_err[3], DW_paramEst1 *localDW, P_paramEst1 *localP)
{
  // InitializeConditions for Delay: '<S300>/Delay1'
  localDW->icLoad = true;

  // InitializeConditions for UnitDelay: '<S298>/Unit Delay3'
  localDW->UnitDelay3_DSTATE[0] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[1] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[2] = localP->UnitDelay3_InitialCondition;

  // InitializeConditions for Delay: '<S300>/Delay'
  localDW->icLoad_n = true;

  // SystemInitialize for Outport: '<S6>/theta'
  for (int32_T i{0}; i < 12; i++) {
    rty_theta[i] = localP->theta_Y0;
  }

  // End of SystemInitialize for Outport: '<S6>/theta'

  // SystemInitialize for Outport: '<S6>/P'
  for (int32_T i{0}; i < 144; i++) {
    rty_P[i] = localP->P_Y0;
  }

  // End of SystemInitialize for Outport: '<S6>/P'

  // SystemInitialize for Outport: '<S6>/err'
  rty_err[0] = localP->err_Y0;
  rty_err[1] = localP->err_Y0;
  rty_err[2] = localP->err_Y0;
}

//
// Output and update for function-call system:
//    '<S1>/paramEst1'
//    '<S1>/paramEst2'
//
void SupervisoryController::paramEst1(const real_T rtu_y[3], const real_T
  rtu_y0[3], const real_T rtu_u[3], const real_T rtu_u0[3], const boolean_T
  rtu_EN[3], const real_T rtu_theta0[12], const real_T rtu_thetaSgn[12],
  boolean_T rtu_rstTheta, const real_T rtu_P0[144], boolean_T rtu_rstP, real_T
  rtu_p_, real_T rtu_lambda, real_T rty_theta[12], real_T rty_P[144], real_T
  rty_err[3], DW_paramEst1 *localDW, ZCE_paramEst1 *localZCE)
{
  static const int8_T b_b_0[9]{ 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  __m128d tmp_1;
  __m128d tmp_2;
  real_T rtb_L_k[144];
  real_T rtb_L_o[36];
  real_T rtb_phi[36];
  real_T tmp[36];
  real_T tmp_0[36];
  real_T regs[12];
  real_T b_b[9];
  real_T b_x[9];
  real_T rtb_Add1_pv[3];
  real_T rtb_Add3_idx_0;
  real_T rtb_Add3_idx_1;
  real_T rtb_Add3_idx_2;
  int32_T i;
  int32_T ibtile;
  int32_T ibtile_tmp;
  int32_T p1;
  int32_T p2;
  int32_T p2_tmp;
  int32_T scalarLB;

  // Delay: '<S300>/Delay1'
  localDW->icLoad = ((rtu_rstP && (static_cast<uint32_T>
    (localZCE->Delay1_Reset_ZCE) != POS_ZCSIG)) || localDW->icLoad);
  localZCE->Delay1_Reset_ZCE = rtu_rstP ? static_cast<ZCSigState>(1) :
    static_cast<ZCSigState>(0);
  if (localDW->icLoad) {
    (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rtu_P0[0], 144U * sizeof
                      (real_T));
  }

  // Sum: '<S298>/Add1'
  rtb_Add1_pv[0] = rtu_y[0] - rtu_y0[0];

  // Sum: '<S298>/Add3'
  rtb_Add3_idx_0 = rtu_u[0] - rtu_u0[0];

  // Sum: '<S298>/Add1'
  rtb_Add1_pv[1] = rtu_y[1] - rtu_y0[1];

  // Sum: '<S298>/Add3'
  rtb_Add3_idx_1 = rtu_u[1] - rtu_u0[1];

  // Sum: '<S298>/Add1'
  rtb_Add1_pv[2] = rtu_y[2] - rtu_y0[2];

  // Sum: '<S298>/Add3'
  rtb_Add3_idx_2 = rtu_u[2] - rtu_u0[2];

  // MATLAB Function: '<S298>/MATLAB Function1' incorporates:
  //   Sum: '<S298>/Add1'
  //   Sum: '<S298>/Add3'

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1': '<S299>:1' 
  // '<S299>:1:2' [z, phi] = getRegressors_(y, yPrev, u, sign_, no, ni, np, dt, mdlNum); 
  // 'getRegressors_:3' z = zeros(no, 1);
  // 'getRegressors_:4' phi = zeros(no*np, no);
  (void)std::memset(&rtb_phi[0], 0, 36U * sizeof(real_T));

  // 'getRegressors_:6' regs = [y'; repmat(u, 1, no)];
  for (p1 = 0; p1 < 3; p1++) {
    ibtile = p1 * 3;
    b_b[ibtile] = rtb_Add3_idx_0;
    b_b[ibtile + 1] = rtb_Add3_idx_1;
    b_b[ibtile + 2] = rtb_Add3_idx_2;
    i = p1 << 2UL;
    regs[i] = rtb_Add1_pv[p1];
    regs[i + 1] = b_b[3 * p1];
    regs[i + 2] = b_b[3 * p1 + 1];
    regs[i + 3] = b_b[3 * p1 + 2];
  }

  //  normalized regressor matrix
  // 'getRegressors_:7' for i=1:no
  for (p1 = 0; p1 < 3; p1++) {
    //  concatenate each col vector of regs into blkdiag form and apply sign
    // 'getRegressors_:8' phi( ((i-1)*np + 1):i*np, i ) = sign_( ((i-1)*np + 1):i*np ).*regs(:,i); 
    p2_tmp = p1 << 2UL;
    p2 = p2_tmp;
    ibtile_tmp = (p1 + 1) << 2UL;
    if (p2_tmp + 1 > ibtile_tmp) {
      p2 = 0;
      i = 0;
    } else {
      i = ibtile_tmp;
    }

    ibtile = p2_tmp;
    if (p2_tmp + 1 > ibtile_tmp) {
      ibtile = 0;
      ibtile_tmp = 0;
    }

    if (i - p2 == 4) {
      int32_T vectorUB;
      ibtile_tmp -= ibtile;
      scalarLB = (ibtile_tmp / 2) << 1UL;
      vectorUB = scalarLB - 2;
      for (i = 0; i <= vectorUB; i += 2) {
        tmp_1 = _mm_loadu_pd(&regs[p2_tmp + i]);
        tmp_2 = _mm_loadu_pd(&rtu_thetaSgn[p2 + i]);
        (void)_mm_storeu_pd(&rtb_phi[(ibtile + i) + 12 * p1], _mm_mul_pd(tmp_1,
          tmp_2));
      }

      for (i = scalarLB; i < ibtile_tmp; i++) {
        rtb_phi[(ibtile + i) + 12 * p1] = regs[(p1 << 2UL) + i] *
          rtu_thetaSgn[p2 + i];
      }
    } else {
      binary_expand_op(rtb_phi, ibtile, ibtile_tmp, p1, rtu_thetaSgn, p2, i - 1,
                       regs);
    }
  }

  // Delay: '<S300>/Delay'
  // 'getRegressors_:11' z = y - yPrev;
  localDW->icLoad_n = ((rtu_rstTheta && (static_cast<uint32_T>
    (localZCE->Delay_Reset_ZCE) != POS_ZCSIG)) || localDW->icLoad_n);
  localZCE->Delay_Reset_ZCE = rtu_rstTheta ? static_cast<ZCSigState>(1) :
    static_cast<ZCSigState>(0);
  if (localDW->icLoad_n) {
    (void)std::memcpy(&localDW->Delay_DSTATE[0], &rtu_theta0[0], 12U * sizeof
                      (real_T));
  }

  for (i = 0; i < 12; i++) {
    // Math: '<S300>/Transpose' incorporates:
    //   MATLAB Function: '<S300>/MATLAB Function'

    tmp[3 * i] = rtb_phi[i];
    tmp[3 * i + 1] = rtb_phi[i + 12];
    tmp[3 * i + 2] = rtb_phi[i + 24];
  }

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function': '<S301>:1' 
  // '<S301>:1:2' [dtheta, dP, L] = rls_(theta, phi, epsil, EN, p_, dPmod_, lambda, P, no, ni, np); 
  // 'rls_:3' dtheta = zeros(no*np, 1);
  // 'rls_:4' dP = zeros(no*np,no*np);
  // 'rls_:5' L = zeros(no*np, no);
  // 'rls_:7' L = P*phi*inv(lambda*eye(no) + phi'*P*phi);
  for (i = 0; i < 3; i++) {
    // Sum: '<S300>/Sum2' incorporates:
    //   Delay: '<S300>/Delay'
    //   MATLAB Function: '<S298>/MATLAB Function1'
    //   Math: '<S300>/Transpose'
    //   Sum: '<S298>/Add1'
    //   UnitDelay: '<S298>/Unit Delay3'

    rtb_Add3_idx_0 = 0.0;
    for (p2 = 0; p2 < 12; p2++) {
      // MATLAB Function: '<S300>/MATLAB Function'
      p1 = 3 * p2 + i;
      rtb_Add3_idx_0 += tmp[p1] * localDW->Delay_DSTATE[p2];

      // MATLAB Function: '<S300>/MATLAB Function' incorporates:
      //   Delay: '<S300>/Delay'
      //   Delay: '<S300>/Delay1'

      tmp_0[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        tmp_0[p1] += tmp[3 * ibtile + i] * localDW->Delay1_DSTATE[12 * p2 +
          ibtile];
      }
    }

    rty_err[i] = (rtb_Add1_pv[i] - localDW->UnitDelay3_DSTATE[i]) -
      rtb_Add3_idx_0;

    // End of Sum: '<S300>/Sum2'

    // MATLAB Function: '<S300>/MATLAB Function'
    for (p2 = 0; p2 < 3; p2++) {
      rtb_Add3_idx_0 = 0.0;
      for (p1 = 0; p1 < 12; p1++) {
        rtb_Add3_idx_0 += tmp_0[3 * p1 + i] * rtb_phi[12 * p2 + p1];
      }

      p1 = 3 * p2 + i;
      b_b[p1] = static_cast<real_T>(b_b_0[p1]) * rtu_lambda + rtb_Add3_idx_0;
    }
  }

  // MATLAB Function: '<S300>/MATLAB Function' incorporates:
  //   Delay: '<S300>/Delay'
  //   Delay: '<S300>/Delay1'

  (void)std::memcpy(&b_x[0], &b_b[0], 9U * sizeof(real_T));
  p1 = 0;
  p2 = 3;
  ibtile = 6;
  rtb_Add3_idx_0 = std::abs(b_b[0]);
  rtb_Add3_idx_1 = std::abs(b_b[1]);
  rtb_Add3_idx_2 = std::abs(b_b[2]);
  if ((rtb_Add3_idx_1 > rtb_Add3_idx_0) && (rtb_Add3_idx_1 > rtb_Add3_idx_2)) {
    p1 = 3;
    p2 = 0;
    b_x[0] = b_b[1];
    b_x[1] = b_b[0];
    b_x[3] = b_b[4];
    b_x[4] = b_b[3];
    b_x[6] = b_b[7];
    b_x[7] = b_b[6];
  } else if (rtb_Add3_idx_2 > rtb_Add3_idx_0) {
    p1 = 6;
    ibtile = 0;
    b_x[0] = b_b[2];
    b_x[2] = b_b[0];
    b_x[3] = b_b[5];
    b_x[5] = b_b[3];
    b_x[6] = b_b[8];
    b_x[8] = b_b[6];
  } else {
    // no actions
  }

  b_x[1] /= b_x[0];
  b_x[2] /= b_x[0];
  b_x[4] -= b_x[1] * b_x[3];
  b_x[5] -= b_x[2] * b_x[3];
  b_x[7] -= b_x[1] * b_x[6];
  b_x[8] -= b_x[2] * b_x[6];
  if (std::abs(b_x[5]) > std::abs(b_x[4])) {
    i = p2;
    p2 = ibtile;
    ibtile = i;
    rtb_Add3_idx_0 = b_x[1];
    b_x[1] = b_x[2];
    b_x[2] = rtb_Add3_idx_0;
    rtb_Add3_idx_0 = b_x[4];
    b_x[4] = b_x[5];
    b_x[5] = rtb_Add3_idx_0;
    rtb_Add3_idx_0 = b_x[7];
    b_x[7] = b_x[8];
    b_x[8] = rtb_Add3_idx_0;
  }

  b_x[5] /= b_x[4];
  b_x[8] -= b_x[5] * b_x[7];
  rtb_Add3_idx_0 = (b_x[1] * b_x[5] - b_x[2]) / b_x[8];
  rtb_Add3_idx_1 = -(b_x[7] * rtb_Add3_idx_0 + b_x[1]) / b_x[4];
  b_b[p1] = ((1.0 - b_x[3] * rtb_Add3_idx_1) - b_x[6] * rtb_Add3_idx_0) / b_x[0];
  b_b[p1 + 1] = rtb_Add3_idx_1;
  b_b[p1 + 2] = rtb_Add3_idx_0;
  rtb_Add3_idx_0 = -b_x[5] / b_x[8];
  rtb_Add3_idx_1 = (1.0 - b_x[7] * rtb_Add3_idx_0) / b_x[4];
  b_b[p2] = -(b_x[3] * rtb_Add3_idx_1 + b_x[6] * rtb_Add3_idx_0) / b_x[0];
  b_b[p2 + 1] = rtb_Add3_idx_1;
  b_b[p2 + 2] = rtb_Add3_idx_0;
  rtb_Add3_idx_0 = 1.0 / b_x[8];
  rtb_Add3_idx_1 = -b_x[7] * rtb_Add3_idx_0 / b_x[4];
  b_b[ibtile] = -(b_x[3] * rtb_Add3_idx_1 + b_x[6] * rtb_Add3_idx_0) / b_x[0];
  b_b[ibtile + 1] = rtb_Add3_idx_1;
  b_b[ibtile + 2] = rtb_Add3_idx_0;

  // 'rls_:9' dtheta(1:no*np, 1) = L*epsil;
  // 'rls_:10' dP(1:no*np, 1:no*np) = L*phi'*P;
  for (i = 0; i < 12; i++) {
    for (p2 = 0; p2 < 3; p2++) {
      p1 = 12 * p2 + i;
      tmp_0[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        tmp_0[p1] += localDW->Delay1_DSTATE[12 * ibtile + i] * rtb_phi[12 * p2 +
          ibtile];
      }
    }

    regs[i] = 0.0;
    for (p2 = 0; p2 < 3; p2++) {
      p1 = 12 * p2 + i;
      rtb_L_o[p1] = 0.0;
      rtb_L_o[p1] += b_b[3 * p2] * tmp_0[i];
      rtb_L_o[p1] += b_b[3 * p2 + 1] * tmp_0[i + 12];
      rtb_L_o[p1] += b_b[3 * p2 + 2] * tmp_0[i + 24];
      regs[i] += rtb_L_o[p1] * rty_err[p2];
    }

    rty_theta[i] = regs[i];
    for (p2 = 0; p2 < 12; p2++) {
      p1 = 12 * p2 + i;
      rtb_L_k[p1] = 0.0;
      rtb_L_k[p1] += tmp[3 * p2] * rtb_L_o[i];
      rtb_L_k[p1] += tmp[3 * p2 + 1] * rtb_L_o[i + 12];
      rtb_L_k[p1] += tmp[3 * p2 + 2] * rtb_L_o[i + 24];
    }

    for (p2 = 0; p2 < 12; p2++) {
      p1 = 12 * p2 + i;
      rty_P[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        rty_P[p1] += rtb_L_k[12 * ibtile + i] * localDW->Delay1_DSTATE[12 * p2 +
          ibtile];
      }
    }
  }

  //  parameter projection
  // 'rls_:13' for i = 1:no*np
  for (p1 = 0; p1 < 12; p1++) {
    // 'rls_:14' if mod(i-1, np) == 0
    if (p1 == 0) {
      i = 0;
    } else {
      i = static_cast<int32_T>(std::fmod((static_cast<real_T>(p1) + 1.0) - 1.0,
        4.0));
    }

    if (i == 0) {
      //  a_i requires lower and upper bounds
      // 'rls_:15' if ~( (theta(i) + dtheta(i) > p_) || ... % lower bound
      // 'rls_:16'               ( (theta(i) + dtheta(i) == p_) && (dtheta(i) >= 0) ) ) ... 
      // 'rls_:17'          && ~( (theta(i) + dtheta(i) < p_+0.5) || ... % upper bound 
      // 'rls_:18'               ( (theta(i) + dtheta(i) == p_+0.5) && (dtheta(i) <= 0) ) ) 
      rtb_Add3_idx_0 = localDW->Delay_DSTATE[p1] + rty_theta[p1];
      if ((!(rtb_Add3_idx_0 > rtu_p_)) && ((!(rtb_Add3_idx_0 == rtu_p_)) ||
           (!(rty_theta[p1] >= 0.0)))) {
        rtb_Add3_idx_0 = localDW->Delay_DSTATE[p1] + rty_theta[p1];
        if ((!(rtb_Add3_idx_0 < rtu_p_ + 0.5)) && ((!(rtu_p_ + 0.5 ==
               rtb_Add3_idx_0)) || (!(rty_theta[p1] <= 0.0)))) {
          //  % lower bound
          //  % upper bound
          // 'rls_:19' dtheta(i) = 0;
          rty_theta[p1] = 0.0;

          // 'rls_:20' dP(i,i) = 0;
          // 'rls_:21' dP(i,i) = 0;
          rty_P[p1 + 12 * p1] = 0.0;

          //  dP(i, i:i+np-1) = 0;
          //  dP(i:i+np-1, i) = 0;
        }
      }
    } else {
      // 'rls_:25' else
      //  b_ij only requires lower bound
      // 'rls_:26' if ~( (theta(i) + dtheta(i) > p_) || ...
      // 'rls_:27'               ( (theta(i) + dtheta(i) == p_) && (dtheta(i) >= 0) ) ) 
      rtb_Add3_idx_0 = localDW->Delay_DSTATE[p1] + rty_theta[p1];
      if ((!(rtb_Add3_idx_0 > rtu_p_)) && ((!(rtb_Add3_idx_0 == rtu_p_)) ||
           (!(rty_theta[p1] >= 0.0)))) {
        // 'rls_:28' dtheta(i) = 0;
        rty_theta[p1] = 0.0;

        // 'rls_:29' dP(i,i) = 0;
        // 'rls_:30' dP(i,i) = 0;
        rty_P[p1 + 12 * p1] = 0.0;
      }
    }
  }

  // 'rls_:35' for i = 1:no
  for (p1 = 0; p1 < 3; p1++) {
    // 'rls_:36' if ~EN(i)
    if (!rtu_EN[p1]) {
      //  set all elements of dtheta and dP to 0
      // 'rls_:37' dtheta( ((i-1)*np + 1):i*np ) = 0;
      p2_tmp = p1 << 2UL;
      p2 = p2_tmp;
      ibtile_tmp = (p1 + 1) << 2UL;
      ibtile = ibtile_tmp;
      if (p2_tmp + 1 > ibtile_tmp) {
        p2 = 0;
        ibtile = 0;
      }

      scalarLB = ibtile - p2;
      if (scalarLB - 1 >= 0) {
        (void)std::memset(&rty_theta[p2], 0, static_cast<uint32_T>
                          (static_cast<int32_T>((scalarLB + p2) - p2)) * sizeof
                          (real_T));
      }

      // 'rls_:38' dP( ((i-1)*np + 1):i*np,: ) = 0;
      ibtile = p2_tmp;
      p2 = ibtile_tmp;
      if (p2_tmp + 1 > ibtile_tmp) {
        ibtile = 0;
        p2 = 0;
      }

      scalarLB = p2 - ibtile;
      for (i = 0; i < 12; i++) {
        if (scalarLB - 1 >= 0) {
          (void)std::memset(&rty_P[i * 12 + ibtile], 0, static_cast<uint32_T>(
            static_cast<int32_T>((scalarLB + ibtile) - ibtile)) * sizeof(real_T));
        }
      }

      // 'rls_:39' dP( :,((i-1)*np + 1):i*np ) = 0;
      ibtile = p2_tmp;
      p2 = ibtile_tmp;
      if (p2_tmp + 1 > ibtile_tmp) {
        ibtile = 0;
        p2 = 0;
      }

      scalarLB = p2 - ibtile;
      for (i = 0; i < scalarLB; i++) {
        (void)std::memset(&rty_P[i * 12 + ibtile * 12], 0, 12U * sizeof(real_T));
      }
    }
  }

  //  if ~( all(theta+dtheta > p_) || ...
  //        (all(theta+dtheta >= p_) && all(b_(find(a_)))) )
  //      dtheta = zeros(numChs,numChs);
  //      dP = zeros(numChs,numChs);
  //      % dP = - dPmod_*eye(numChs); % covariance "kick"
  //  end
  //  covariance resetting (doesn't work well)
  //  if norm(P - dP) < 5e-2
  //      dP = - ( P + 5e-1*eye(3) );
  //  end
  //  if min(svd(P - dP)) < 1e-3
  //      dP = - ( P + 1*eye(3) );
  //  end
  for (i = 0; i <= 142; i += 2) {
    // Delay: '<S300>/Delay1'
    tmp_1 = _mm_loadu_pd(&localDW->Delay1_DSTATE[i]);

    // Product: '<S300>/Product1' incorporates:
    //   Delay: '<S300>/Delay1'

    tmp_2 = _mm_loadu_pd(&rty_P[i]);
    (void)_mm_storeu_pd(&rty_P[i], _mm_div_pd(_mm_sub_pd(tmp_1, tmp_2),
      _mm_set1_pd(rtu_lambda)));
  }

  for (i = 0; i <= 10; i += 2) {
    // Delay: '<S300>/Delay'
    tmp_1 = _mm_loadu_pd(&localDW->Delay_DSTATE[i]);

    // Sum: '<S300>/Sum' incorporates:
    //   Delay: '<S300>/Delay'

    tmp_2 = _mm_loadu_pd(&rty_theta[i]);
    (void)_mm_storeu_pd(&rty_theta[i], _mm_add_pd(tmp_1, tmp_2));
  }

  // Update for Delay: '<S300>/Delay1'
  localDW->icLoad = false;
  (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rty_P[0], 144U * sizeof(real_T));

  // Update for UnitDelay: '<S298>/Unit Delay3' incorporates:
  //   Sum: '<S298>/Add1'

  localDW->UnitDelay3_DSTATE[0] = rtb_Add1_pv[0];
  localDW->UnitDelay3_DSTATE[1] = rtb_Add1_pv[1];
  localDW->UnitDelay3_DSTATE[2] = rtb_Add1_pv[2];

  // Update for Delay: '<S300>/Delay'
  localDW->icLoad_n = false;
  (void)std::memcpy(&localDW->Delay_DSTATE[0], &rty_theta[0], 12U * sizeof
                    (real_T));
}

//
// Output and update for atomic system:
//    '<S112>/ScalarExpansionR'
//    '<S182>/ScalarExpansionR'
//    '<S252>/ScalarExpansionR'
//
void SupervisoryController::ScalarExpansionR(const real_T rtu_u[9], real_T
  rty_y[9])
{
  int32_T tmp;

  // MATLAB Function: '<S135>/ScalarExpansion'
  //  ctrlScalarExpansion Helper function for scalar expansion.
  //
  //  	y  = ctrlScalarExpansion(u,n)
  //
  //    An n-ny-n matrix y is created based on u. If u is a scalar, y has u
  //    on its diagonals. If u is a vector, y has the elements of u on its
  //    diagonals. If u is a matrix y = (u+u.')/2.
  //
  //    When u is scalar or vector, we enforce symmetric positive-definiteness.
  //    When u is a matrix, we enly enforce symmetry.
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S157>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S157>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  tmp = 0;
  for (int32_T i{0}; i < 3; i++) {
    rty_y[tmp] = (rtu_u[tmp] + rtu_u[i]) / 2.0;
    rty_y[tmp + 1] = (rtu_u[tmp + 1] + rtu_u[i + 3]) / 2.0;
    rty_y[tmp + 2] = (rtu_u[tmp + 2] + rtu_u[i + 6]) / 2.0;
    tmp += 3;
  }

  // End of MATLAB Function: '<S135>/ScalarExpansion'
}

// Function for MATLAB Function: '<S184>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv(const real_T A[15], const real_T B_0[9],
  real_T Y[15])
{
  real_T b_A[9];
  real_T a21;
  real_T maxval;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  int32_T rtemp;
  (void)std::memcpy(&b_A[0], &B_0[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = std::abs(B_0[0]);
  a21 = std::abs(B_0[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (std::abs(B_0[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = B_0[r2] / B_0[r1];
  b_A[r3] /= b_A[r1];
  b_A[r2 + 3] -= b_A[r1 + 3] * b_A[r2];
  b_A[r3 + 3] -= b_A[r1 + 3] * b_A[r3];
  b_A[r2 + 6] -= b_A[r1 + 6] * b_A[r2];
  b_A[r3 + 6] -= b_A[r1 + 6] * b_A[r3];
  if (std::abs(b_A[r3 + 3]) > std::abs(b_A[r2 + 3])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[r3 + 3] /= b_A[r2 + 3];
  b_A[r3 + 6] -= b_A[r3 + 3] * b_A[r2 + 6];
  for (rtemp = 0; rtemp < 5; rtemp++) {
    int32_T Y_tmp;
    int32_T Y_tmp_0;
    int32_T Y_tmp_1;
    Y_tmp = 5 * r1 + rtemp;
    Y[Y_tmp] = A[rtemp] / b_A[r1];
    Y_tmp_0 = 5 * r2 + rtemp;
    Y[Y_tmp_0] = A[rtemp + 5] - b_A[r1 + 3] * Y[Y_tmp];
    Y_tmp_1 = 5 * r3 + rtemp;
    Y[Y_tmp_1] = A[rtemp + 10] - b_A[r1 + 6] * Y[Y_tmp];
    Y[Y_tmp_0] /= b_A[r2 + 3];
    Y[Y_tmp_1] -= b_A[r2 + 6] * Y[Y_tmp_0];
    Y[Y_tmp_1] /= b_A[r3 + 6];
    Y[Y_tmp_0] -= b_A[r3 + 3] * Y[Y_tmp_1];
    Y[Y_tmp] -= Y[Y_tmp_1] * b_A[r3];
    Y[Y_tmp] -= Y[Y_tmp_0] * b_A[r2];
  }
}

//
// Output and update for atomic system:
//    '<S182>/CalculatePL'
//    '<S252>/CalculatePL'
//
void SupervisoryController::CalculatePL(const real_T rtu_Ak[25], const real_T
  rtu_Ck[15], const real_T rtu_Qbark[25], const real_T rtu_Rbark[9], const
  real_T rtu_Nbark[15], boolean_T rtu_Enablek, const real_T rtu_Pk[25], real_T
  rty_Mk[15], real_T rty_Lk[15], real_T rty_Zk[25], real_T rty_Pk1[25])
{
  real_T Abar[25];
  real_T Abar_0[25];
  real_T Abar_1[25];
  real_T rty_Mk_0[25];
  real_T NRinv[15];
  real_T rtu_Ck_0[15];
  real_T yCov[9];
  int8_T b_I[25];

  // MATLAB Function: '<S184>/Discrete-Time KF - Calculate PLMZ'
  //  See help of ctrlKalmanFilterDTCalculatePL.m
  // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S222>:1' 
  //    Copyright 2014 The MathWorks, Inc.
  // '<S222>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
  if (rtu_Enablek) {
    __m128d tmp_0;
    __m128d tmp_1;
    real_T tmp;
    int32_T i;
    int32_T k;
    int32_T rtu_Ak_tmp;
    int32_T rtu_Pk_tmp;
    int32_T rty_Mk_tmp;
    int32_T yCov_tmp;
    i = 0;
    for (k = 0; k < 3; k++) {
      rty_Mk_tmp = 0;
      rtu_Ak_tmp = 0;
      for (yCov_tmp = 0; yCov_tmp < 5; yCov_tmp++) {
        int32_T NRinv_tmp;
        NRinv_tmp = rty_Mk_tmp + k;
        NRinv[yCov_tmp + i] = rtu_Ck[NRinv_tmp];
        rtu_Ck_0[NRinv_tmp] = 0.0;
        rtu_Pk_tmp = 0;
        for (int32_T i_0{0}; i_0 < 5; i_0++) {
          rtu_Ck_0[NRinv_tmp] += rtu_Ck[rtu_Pk_tmp + k] * rtu_Pk[i_0 +
            rtu_Ak_tmp];
          rtu_Pk_tmp += 3;
        }

        rty_Mk_tmp += 3;
        rtu_Ak_tmp += 5;
      }

      i += 5;
    }

    for (i = 0; i < 3; i++) {
      k = 0;
      rty_Mk_tmp = 0;
      for (rtu_Ak_tmp = 0; rtu_Ak_tmp < 3; rtu_Ak_tmp++) {
        tmp = 0.0;
        yCov_tmp = 0;
        for (rtu_Pk_tmp = 0; rtu_Pk_tmp < 5; rtu_Pk_tmp++) {
          tmp += rtu_Ck_0[yCov_tmp + i] * NRinv[rtu_Pk_tmp + rty_Mk_tmp];
          yCov_tmp += 3;
        }

        yCov_tmp = k + i;
        yCov[yCov_tmp] = rtu_Rbark[yCov_tmp] + tmp;
        k += 3;
        rty_Mk_tmp += 5;
      }
    }

    for (i = 0; i < 5; i++) {
      for (k = 0; k < 5; k++) {
        rtu_Ak_tmp = 5 * k + i;
        Abar[rtu_Ak_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          Abar[rtu_Ak_tmp] += rtu_Ak[5 * rty_Mk_tmp + i] * rtu_Pk[5 * k +
            rty_Mk_tmp];
        }
      }

      for (k = 0; k < 3; k++) {
        tmp = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          tmp += Abar[5 * rty_Mk_tmp + i] * NRinv[5 * k + rty_Mk_tmp];
        }

        rtu_Ak_tmp = 5 * k + i;
        rtu_Ck_0[rtu_Ak_tmp] = rtu_Nbark[rtu_Ak_tmp] + tmp;
      }
    }

    mrdiv(rtu_Ck_0, yCov, rty_Lk);
    i = 0;
    for (k = 0; k < 3; k++) {
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
        rtu_Pk_tmp = rty_Mk_tmp + i;
        rtu_Ck_0[rtu_Pk_tmp] = 0.0;
        rtu_Ak_tmp = 0;
        for (yCov_tmp = 0; yCov_tmp < 5; yCov_tmp++) {
          rtu_Ck_0[rtu_Pk_tmp] += rtu_Pk[rtu_Ak_tmp + rty_Mk_tmp] *
            NRinv[yCov_tmp + i];
          rtu_Ak_tmp += 5;
        }
      }

      i += 5;
    }

    mrdiv(rtu_Ck_0, yCov, rty_Mk);
    for (i = 0; i < 25; i++) {
      b_I[i] = 0;
    }

    k = 0;
    for (i = 0; i < 5; i++) {
      b_I[k] = 1;
      k += 6;
    }

    for (i = 0; i < 5; i++) {
      for (k = 0; k < 5; k++) {
        rtu_Pk_tmp = 5 * k + i;
        Abar[rtu_Pk_tmp] = static_cast<real_T>(b_I[rtu_Pk_tmp]) - ((rtu_Ck[3 * k
          + 1] * rty_Mk[i + 5] + rtu_Ck[3 * k] * rty_Mk[i]) + rtu_Ck[3 * k + 2] *
          rty_Mk[i + 10]);
      }

      for (k = 0; k < 5; k++) {
        rtu_Pk_tmp = 5 * k + i;
        Abar_0[rtu_Pk_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          Abar_0[rtu_Pk_tmp] += Abar[5 * rty_Mk_tmp + i] * rtu_Pk[5 * k +
            rty_Mk_tmp];
        }
      }

      for (k = 0; k < 3; k++) {
        rty_Mk_tmp = 5 * k + i;
        NRinv[rty_Mk_tmp] = 0.0;
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k] * rty_Mk[i];
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k + 1] * rty_Mk[i + 5];
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k + 2] * rty_Mk[i + 10];
      }
    }

    for (i = 0; i < 5; i++) {
      k = 0;
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
        rtu_Pk_tmp = k + i;
        Abar_1[rtu_Pk_tmp] = 0.0;
        rtu_Ak_tmp = 0;
        for (yCov_tmp = 0; yCov_tmp < 5; yCov_tmp++) {
          Abar_1[rtu_Pk_tmp] += Abar_0[rtu_Ak_tmp + i] * Abar[rtu_Ak_tmp +
            rty_Mk_tmp];
          rtu_Ak_tmp += 5;
        }

        rty_Mk_0[rtu_Pk_tmp] = 0.0;
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i] * rty_Mk[rty_Mk_tmp];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 5] * rty_Mk[rty_Mk_tmp + 5];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 10] * rty_Mk[rty_Mk_tmp + 10];
        k += 5;
      }
    }

    for (i = 0; i <= 22; i += 2) {
      tmp_0 = _mm_loadu_pd(&Abar_1[i]);
      tmp_1 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Zk[i], _mm_add_pd(tmp_0, tmp_1));
    }

    for (i = 24; i < 25; i++) {
      rty_Zk[i] = Abar_1[i] + rty_Mk_0[i];
    }

    mrdiv(rtu_Nbark, rtu_Rbark, NRinv);
    for (i = 0; i < 5; i++) {
      for (k = 0; k < 5; k++) {
        rtu_Pk_tmp = 5 * k + i;
        Abar[rtu_Pk_tmp] = rtu_Ak[rtu_Pk_tmp] - ((rtu_Ck[3 * k + 1] * NRinv[i +
          5] + rtu_Ck[3 * k] * NRinv[i]) + rtu_Ck[3 * k + 2] * NRinv[i + 10]);
      }

      for (k = 0; k < 5; k++) {
        rtu_Pk_tmp = 5 * k + i;
        Abar_0[rtu_Pk_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          Abar_0[rtu_Pk_tmp] += Abar[5 * rty_Mk_tmp + i] * rty_Zk[5 * k +
            rty_Mk_tmp];
        }
      }
    }

    for (i = 0; i < 5; i++) {
      k = 0;
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
        tmp = 0.0;
        rtu_Ak_tmp = 0;
        for (yCov_tmp = 0; yCov_tmp < 5; yCov_tmp++) {
          tmp += Abar_0[rtu_Ak_tmp + i] * Abar[rtu_Ak_tmp + rty_Mk_tmp];
          rtu_Ak_tmp += 5;
        }

        rtu_Pk_tmp = k + i;
        Abar_1[rtu_Pk_tmp] = rtu_Qbark[rtu_Pk_tmp] + tmp;
        rty_Mk_0[rtu_Pk_tmp] = 0.0;
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i] * rtu_Nbark[rty_Mk_tmp];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 5] * rtu_Nbark[rty_Mk_tmp + 5];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 10] * rtu_Nbark[rty_Mk_tmp + 10];
        k += 5;
      }
    }

    for (i = 0; i <= 22; i += 2) {
      tmp_0 = _mm_loadu_pd(&Abar_1[i]);
      tmp_1 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Pk1[i], _mm_sub_pd(tmp_0, tmp_1));
    }

    for (i = 24; i < 25; i++) {
      rty_Pk1[i] = Abar_1[i] - rty_Mk_0[i];
    }
  } else {
    (void)std::memset(&rty_Lk[0], 0, 15U * sizeof(real_T));
    (void)std::memset(&rty_Mk[0], 0, 15U * sizeof(real_T));
    (void)std::memcpy(&rty_Zk[0], &rtu_Pk[0], 25U * sizeof(real_T));
    for (int32_T i{0}; i < 5; i++) {
      int32_T rty_Mk_tmp;
      for (int32_T k{0}; k < 5; k++) {
        int32_T rtu_Ak_tmp;
        rtu_Ak_tmp = 5 * k + i;
        Abar[rtu_Ak_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          Abar[rtu_Ak_tmp] += rtu_Ak[5 * rty_Mk_tmp + i] * rtu_Pk[5 * k +
            rty_Mk_tmp];
        }
      }

      for (int32_T k{0}; k < 5; k++) {
        real_T tmp;
        tmp = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 5; rty_Mk_tmp++) {
          tmp += Abar[5 * rty_Mk_tmp + i] * rtu_Ak[5 * rty_Mk_tmp + k];
        }

        rty_Mk_tmp = 5 * k + i;
        rty_Pk1[rty_Mk_tmp] = rtu_Qbark[rty_Mk_tmp] + tmp;
      }
    }
  }

  // End of MATLAB Function: '<S184>/Discrete-Time KF - Calculate PLMZ'
}

//
// Output and update for atomic system:
//    '<S223>/SqrtUsedFcn'
//    '<S293>/SqrtUsedFcn'
//
void SupervisoryController::SqrtUsedFcn(const real_T rtu_u[25], boolean_T
  rtu_isSqrtUsed, real_T rty_P[25])
{
  //  Determine if the Square-Root algorithm was used
  // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S224>:1' 
  // '<S224>:1:4' if isSqrtUsed
  if (rtu_isSqrtUsed) {
    // '<S224>:1:5' P = u*u.';
    for (int32_T i{0}; i < 5; i++) {
      int32_T tmp;
      tmp = 0;
      for (int32_T i_0{0}; i_0 < 5; i_0++) {
        int32_T tmp_0;
        int32_T tmp_1;
        tmp_1 = tmp + i;
        rty_P[tmp_1] = 0.0;
        tmp_0 = 0;
        for (int32_T i_1{0}; i_1 < 5; i_1++) {
          rty_P[tmp_1] += rtu_u[tmp_0 + i] * rtu_u[tmp_0 + i_0];
          tmp_0 += 5;
        }

        tmp += 5;
      }
    }
  } else {
    // '<S224>:1:6' else
    // '<S224>:1:7' P = u;
    (void)std::memcpy(&rty_P[0], &rtu_u[0], 25U * sizeof(real_T));
  }
}

//
// System initialize for enable system:
//    '<S201>/MeasurementUpdate'
//    '<S271>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Init(real_T rty_Lykyhatkk1[5],
  P_MeasurementUpdate *localP)
{
  // SystemInitialize for Outport: '<S225>/L*(y[k]-yhat[k|k-1])'
  for (int32_T i{0}; i < 5; i++) {
    rty_Lykyhatkk1[i] = localP->Lykyhatkk1_Y0;
  }

  // End of SystemInitialize for Outport: '<S225>/L*(y[k]-yhat[k|k-1])'
}

//
// Disable for enable system:
//    '<S201>/MeasurementUpdate'
//    '<S271>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[5],
  DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  // Disable for Outport: '<S225>/L*(y[k]-yhat[k|k-1])'
  for (int32_T i{0}; i < 5; i++) {
    // Outputs for Enabled SubSystem: '<S201>/MeasurementUpdate' incorporates:
    //   EnablePort: '<S225>/Enable'

    rty_Lykyhatkk1[i] = localP->Lykyhatkk1_Y0;

    // End of Outputs for SubSystem: '<S201>/MeasurementUpdate'
  }

  // End of Disable for Outport: '<S225>/L*(y[k]-yhat[k|k-1])'
  localDW->MeasurementUpdate_MODE = false;
}

//
// Output and update for enable system:
//    '<S201>/MeasurementUpdate'
//    '<S271>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate(boolean_T rtu_Enable, const real_T
  rtu_Lk[15], const real_T rtu_yk[3], const real_T rtu_Ck[15], const real_T
  rtu_xhatkk1[5], const real_T rtu_Dk[9], const real_T rtu_uk[3], real_T
  rty_Lykyhatkk1[5], DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  real_T rtu_Ck_0[3];
  real_T rtu_Dk_0[3];
  real_T rtu_yk_0[3];

  // Outputs for Enabled SubSystem: '<S201>/MeasurementUpdate' incorporates:
  //   EnablePort: '<S225>/Enable'

  if (rtu_Enable) {
    localDW->MeasurementUpdate_MODE = true;
    for (int32_T i{0}; i < 3; i++) {
      int32_T tmp;

      // Product: '<S225>/C[k]*xhat[k|k-1]'
      rtu_Ck_0[i] = 0.0;
      tmp = 0;
      for (int32_T i_0{0}; i_0 < 5; i_0++) {
        rtu_Ck_0[i] += rtu_Ck[tmp + i] * rtu_xhatkk1[i_0];
        tmp += 3;
      }

      // Product: '<S225>/D[k]*u[k]' incorporates:
      //   Product: '<S225>/C[k]*xhat[k|k-1]'

      rtu_Dk_0[i] = 0.0;
      rtu_Dk_0[i] += rtu_Dk[i] * rtu_uk[0];
      rtu_Dk_0[i] += rtu_Dk[i + 3] * rtu_uk[1];
      rtu_Dk_0[i] += rtu_Dk[i + 6] * rtu_uk[2];

      // Sum: '<S225>/Sum' incorporates:
      //   Product: '<S225>/C[k]*xhat[k|k-1]'
      //   Sum: '<S225>/Add1'

      rtu_yk_0[i] = rtu_yk[i] - (rtu_Ck_0[i] + rtu_Dk_0[i]);
    }

    for (int32_T i{0}; i <= 2; i += 2) {
      __m128d tmp_0;
      __m128d tmp_1;

      // Product: '<S225>/Product3'
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[0])), tmp_1));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i + 5]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[1])), tmp_1));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i + 10]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[2])), tmp_1));
    }

    // Product: '<S225>/Product3'
    for (int32_T i{4}; i < 5; i++) {
      rty_Lykyhatkk1[i] = 0.0;
      rty_Lykyhatkk1[i] += rtu_Lk[i] * rtu_yk_0[0];
      rty_Lykyhatkk1[i] += rtu_Lk[i + 5] * rtu_yk_0[1];
      rty_Lykyhatkk1[i] += rtu_Lk[i + 10] * rtu_yk_0[2];
    }
  } else if (localDW->MeasurementUpdate_MODE) {
    MeasurementUpdate_Disable(rty_Lykyhatkk1, localDW, localP);
  } else {
    // no actions
  }

  // End of Outputs for SubSystem: '<S201>/MeasurementUpdate'
}

//
// Output and update for atomic system:
//    '<S182>/ReducedQRN'
//    '<S252>/ReducedQRN'
//
void SupervisoryController::ReducedQRN(const real_T rtu_G[25], const real_T
  rtu_H[15], const real_T rtu_Q[25], const real_T rtu_R[9], const real_T rtu_N
  [15], real_T rty_Qbar[25], real_T rty_Rbar[9], real_T rty_Nbar[15])
{
  real_T rtu_Q_0[25];
  real_T rtb_Add_a[15];
  real_T rtb_Transpose2_k[15];
  real_T rtu_H_0[9];
  real_T rtu_N_0[9];
  int32_T i;
  int32_T i_0;
  int32_T i_1;
  int32_T rtb_Add_ih_tmp;
  int32_T rtu_H_tmp;
  int32_T rtu_Q_tmp;

  // Product: '<S202>/Product' incorporates:
  //   Math: '<S202>/Transpose1'

  for (i = 0; i < 5; i++) {
    i_1 = 0;
    for (i_0 = 0; i_0 < 5; i_0++) {
      rtu_Q_tmp = i_1 + i;
      rtu_Q_0[rtu_Q_tmp] = 0.0;
      rtb_Add_ih_tmp = 0;
      for (rtu_H_tmp = 0; rtu_H_tmp < 5; rtu_H_tmp++) {
        rtu_Q_0[rtu_Q_tmp] += rtu_Q[rtb_Add_ih_tmp + i] * rtu_G[rtb_Add_ih_tmp +
          i_0];
        rtb_Add_ih_tmp += 5;
      }

      i_1 += 5;
    }
  }

  i = 0;
  for (i_1 = 0; i_1 < 5; i_1++) {
    for (i_0 = 0; i_0 < 5; i_0++) {
      rtb_Add_ih_tmp = i_0 + i;
      rty_Qbar[rtb_Add_ih_tmp] = 0.0;
      rtu_H_tmp = 0;
      for (rtu_Q_tmp = 0; rtu_Q_tmp < 5; rtu_Q_tmp++) {
        rty_Qbar[rtb_Add_ih_tmp] += rtu_G[rtu_H_tmp + i_0] * rtu_Q_0[rtu_Q_tmp +
          i];
        rtu_H_tmp += 5;
      }
    }

    i += 5;
  }

  // End of Product: '<S202>/Product'

  // Math: '<S202>/Transpose2'
  i = 0;
  for (i_1 = 0; i_1 < 3; i_1++) {
    i_0 = 0;
    for (rtb_Add_ih_tmp = 0; rtb_Add_ih_tmp < 5; rtb_Add_ih_tmp++) {
      rtb_Transpose2_k[rtb_Add_ih_tmp + i] = rtu_H[i_0 + i_1];
      i_0 += 3;
    }

    i += 5;
  }

  // End of Math: '<S202>/Transpose2'

  // Sum: '<S202>/Add' incorporates:
  //   Math: '<S202>/Transpose2'
  //   Product: '<S202>/Product1'

  for (i = 0; i < 5; i++) {
    i_1 = 0;
    for (i_0 = 0; i_0 < 3; i_0++) {
      real_T tmp;
      tmp = 0.0;
      rtb_Add_ih_tmp = 0;
      for (rtu_H_tmp = 0; rtu_H_tmp < 5; rtu_H_tmp++) {
        tmp += rtu_Q[rtb_Add_ih_tmp + i] * rtb_Transpose2_k[rtu_H_tmp + i_1];
        rtb_Add_ih_tmp += 5;
      }

      rtb_Add_ih_tmp = i_1 + i;
      rtb_Add_a[rtb_Add_ih_tmp] = rtu_N[rtb_Add_ih_tmp] + tmp;
      i_1 += 5;
    }
  }

  // End of Sum: '<S202>/Add'
  for (i = 0; i < 3; i++) {
    // Product: '<S202>/Product2' incorporates:
    //   Sum: '<S202>/Add'

    for (i_1 = 0; i_1 < 5; i_1++) {
      i_0 = 5 * i + i_1;
      rty_Nbar[i_0] = 0.0;
      for (rtb_Add_ih_tmp = 0; rtb_Add_ih_tmp < 5; rtb_Add_ih_tmp++) {
        rty_Nbar[i_0] += rtu_G[5 * rtb_Add_ih_tmp + i_1] * rtb_Add_a[5 * i +
          rtb_Add_ih_tmp];
      }
    }

    // End of Product: '<S202>/Product2'
    for (i_1 = 0; i_1 < 3; i_1++) {
      // Product: '<S202>/Product3' incorporates:
      //   Product: '<S202>/Product4'

      rtb_Add_ih_tmp = 3 * i_1 + i;
      rtu_H_0[rtb_Add_ih_tmp] = 0.0;

      // Product: '<S202>/Product4'
      rtu_N_0[rtb_Add_ih_tmp] = 0.0;
      for (i_0 = 0; i_0 < 5; i_0++) {
        // Product: '<S202>/Product3' incorporates:
        //   Product: '<S202>/Product4'
        //   Sum: '<S202>/Add'

        rtu_H_tmp = 5 * i_1 + i_0;
        rtu_H_0[rtb_Add_ih_tmp] += rtu_H[3 * i_0 + i] * rtb_Add_a[rtu_H_tmp];

        // Product: '<S202>/Product4' incorporates:
        //   Math: '<S202>/Transpose'
        //   Math: '<S202>/Transpose2'

        rtu_N_0[rtb_Add_ih_tmp] += rtu_N[5 * i + i_0] *
          rtb_Transpose2_k[rtu_H_tmp];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
    __m128d tmp_0;
    __m128d tmp_1;
    __m128d tmp_2;

    // Sum: '<S202>/Add1'
    tmp_0 = _mm_loadu_pd(&rtu_H_0[i]);
    tmp_1 = _mm_loadu_pd(&rtu_N_0[i]);
    tmp_2 = _mm_loadu_pd(&rtu_R[i]);
    (void)_mm_storeu_pd(&rty_Rbar[i], _mm_add_pd(_mm_add_pd(tmp_0, tmp_1), tmp_2));
  }

  // Sum: '<S202>/Add1'
  for (i = 8; i < 9; i++) {
    rty_Rbar[i] = (rtu_H_0[i] + rtu_N_0[i]) + rtu_R[i];
  }
}

//
// Output and update for atomic system:
//    '<S182>/ScalarExpansionQ'
//    '<S252>/ScalarExpansionQ'
//
void SupervisoryController::ScalarExpansionQ(const real_T rtu_u[25], real_T
  rty_y[25])
{
  int32_T tmp;

  // MATLAB Function: '<S204>/ScalarExpansion'
  //  ctrlScalarExpansion Helper function for scalar expansion.
  //
  //  	y  = ctrlScalarExpansion(u,n)
  //
  //    An n-ny-n matrix y is created based on u. If u is a scalar, y has u
  //    on its diagonals. If u is a vector, y has the elements of u on its
  //    diagonals. If u is a matrix y = (u+u.')/2.
  //
  //    When u is scalar or vector, we enforce symmetric positive-definiteness.
  //    When u is a matrix, we enly enforce symmetry.
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S226>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S226>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  tmp = 0;
  for (int32_T i_0{0}; i_0 < 5; i_0++) {
    int32_T tmp_0;
    tmp_0 = 0;
    for (int32_T i{0}; i < 5; i++) {
      int32_T tmp_1;
      tmp_1 = i + tmp;
      rty_y[tmp_1] = (rtu_u[tmp_0 + i_0] + rtu_u[tmp_1]) / 2.0;
      tmp_0 += 5;
    }

    tmp += 5;
  }

  // End of MATLAB Function: '<S204>/ScalarExpansion'
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::do_vectors(const real_T b_data[], const int32_T
  *b_size, real_T c_data[], int32_T c_size[2], int32_T ia_data[], int32_T
  *ia_size, int32_T *ib_size)
{
  int32_T iafirst;
  int32_T ialast;
  int32_T iblast;
  int32_T nc;
  int32_T nia;
  c_size[0] = 1;
  *ib_size = 0;
  nc = 0;
  nia = -1;
  iafirst = 0;
  ialast = 1;
  iblast = 1;
  while ((ialast <= 6) && (iblast <= *b_size)) {
    real_T bk;
    int32_T ak;
    int32_T b_ialast;
    b_ialast = ialast;
    ak = ialast;
    while ((b_ialast < 6) && (b_ialast + 1 == ialast)) {
      b_ialast++;
    }

    ialast = b_ialast;
    bk = b_data[iblast - 1];
    while ((iblast < *b_size) && (b_data[iblast] == bk)) {
      iblast++;
    }

    if (static_cast<real_T>(ak) == bk) {
      ialast = b_ialast + 1;
      iafirst = b_ialast;
      iblast++;
    } else if (static_cast<real_T>(ak) < bk) {
      nc++;
      nia++;
      c_data[nc - 1] = static_cast<real_T>(ak);
      ia_data[nia] = iafirst + 1;
      ialast = b_ialast + 1;
      iafirst = b_ialast;
    } else {
      iblast++;
    }
  }

  while (ialast <= 6) {
    iblast = ialast;
    while ((iblast < 6) && (iblast + 1 == ialast)) {
      iblast++;
    }

    nc++;
    nia++;
    c_data[nc - 1] = (static_cast<real_T>(ialast) - 1.0) + 1.0;
    ia_data[nia] = iafirst + 1;
    ialast = iblast + 1;
    iafirst = iblast;
  }

  if (nia + 1 < 1) {
    nia = -1;
  }

  *ia_size = nia + 1;
  if (nc < 1) {
    c_size[1] = 0;
  } else {
    c_size[1] = nc;
  }
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (std::isinf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
    }
  }

  return y;
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::ppval(const siswYcTR8LLamuD4YWmtXHC *pp, const
  real_T x_data[], const int32_T x_size[2], real_T v_data[], int32_T v_size[2])
{
  real_T xloc;
  int32_T b;
  int32_T high_i;
  int32_T low_i;
  int32_T low_ip1;
  int32_T mid_i;
  v_size[0] = 1;
  v_size[1] = x_size[1];
  b = x_size[1];

#pragma omp parallel for num_threads(omp_get_max_threads()) private(xloc,low_i,low_ip1,high_i,mid_i)

  for (int32_T b_ix = 0; b_ix < b; b_ix++) {
    xloc = x_data[b_ix];
    if (std::isnan(xloc)) {
      v_data[b_ix] = xloc;
    } else {
      low_i = 0;
      low_ip1 = 2;
      high_i = 6;
      while (high_i > low_ip1) {
        mid_i = ((low_i + high_i) + 1) >> 1UL;
        if (xloc >= pp->breaks[mid_i - 1]) {
          low_i = mid_i - 1;
          low_ip1 = mid_i + 1;
        } else {
          high_i = mid_i;
        }
      }

      xloc -= pp->breaks[low_i];
      v_data[b_ix] = (xloc * pp->coefs[low_i] + pp->coefs[low_i + 5]) * xloc +
        pp->coefs[low_i + 10];
    }
  }
}

//
// Function for Chart: '<Root>/SupervisoryController'
// function [trajectory, trajectorySize] = trajGen(event, y_)
//
void SupervisoryController::trajGen(const event_bus *event, const real_T y_[6],
  real_T trajectory[14400], uint16_T *trajectorySize)
{
  cell_wrap_9 coeffsCell;
  siswYcTR8LLamuD4YWmtXHC expl_temp;
  real_T b_coefs[12];
  real_T coefs[9];
  real_T directChs_data[6];
  real_T inferredChs_data[6];
  real_T yDest[6];
  real_T breaks[5];
  real_T newSegmentCoeffs[3];
  real_T holdPoint;
  int32_T ii_data[6];
  int32_T inferredChs_size[2];
  int32_T t_size[2];
  int32_T deltaSign;
  int32_T directChs_size;
  int32_T i;
  int32_T loop_ub;
  int8_T j_data[3];
  boolean_T exitg1;

  // MATLAB Function 'trajGen': '<S1>:527'
  // '<S1>:527:2' [trajectory, trajectorySize] = trajGen_(event, y_, ymax, dt, no); 
  // 'trajGen_:3' traj = zeros(2*no, 2400);
  (void)std::memset(&trajectory[0], 0, 14400U * sizeof(real_T));

  // 'trajGen_:4' yDest = zeros(2*no, 1);
  for (i = 0; i < 6; i++) {
    yDest[i] = 0.0;
  }

  // 'trajGen_:5' directChs = find(event.r);
  i = 0;
  deltaSign = 0;
  exitg1 = false;
  while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
          false) && (deltaSign < 6)) {
    if (event->r[deltaSign] != 0.0) {
      i++;
      ii_data[i - 1] = deltaSign + 1;
      if (i >= 6) {
        exitg1 = true;
      } else {
        deltaSign++;
      }
    } else {
      deltaSign++;
    }
  }

  if (i < 1) {
    i = 0;
  }

  directChs_size = i;
  for (deltaSign = 0; deltaSign < i; deltaSign++) {
    directChs_data[deltaSign] = static_cast<real_T>(ii_data[deltaSign]);
  }

  //  direct channels
  // 'trajGen_:6' inferredChs = setdiff(1:2*no, directChs);
  do_vectors(directChs_data, &directChs_size, inferredChs_data, inferredChs_size,
             ii_data, &i, &deltaSign);

  //  inferred channels
  // 'trajGen_:8' numWaypts = cast(event.moveT/dt, "uint16");
  holdPoint = std::round(event->moveT / rtP.dt);
  if (holdPoint < 65536.0) {
    if (holdPoint >= 0.0) {
      *trajectorySize = static_cast<uint16_T>(holdPoint);
    } else {
      *trajectorySize = 0U;
    }
  } else {
    *trajectorySize = MAX_uint16_T;
  }

  // 'trajGen_:9' assert(numWaypts < 1200);
  // 'trajGen_:11' for i=1:length(directChs)
  for (deltaSign = 0; deltaSign < directChs_size; deltaSign++) {
    loop_ub = static_cast<int32_T>(directChs_data[deltaSign]);

    // Inport: '<Root>/ymax'
    // 'trajGen_:12' yDest(directChs(i)) = ymax(directChs(i)).*event.r(directChs(i)); 
    yDest[loop_ub - 1] = rtU.ymax[loop_ub - 1] * event->r[loop_ub - 1];
  }

  // 'trajGen_:14' for i=1:length(inferredChs)
  i = inferredChs_size[1];
  for (deltaSign = 0; deltaSign < i; deltaSign++) {
    //  apply KCL/conservation of charge at junction
    // 'trajGen_:15' yDest(inferredChs(i)) = -sum(yDest(directChs)) / length(inferredChs); 
    if (directChs_size == 0) {
      holdPoint = 0.0;
    } else {
      holdPoint = yDest[static_cast<int32_T>(directChs_data[0]) - 1];
      for (loop_ub = 2; loop_ub <= directChs_size; loop_ub++) {
        holdPoint += yDest[static_cast<int32_T>(directChs_data[loop_ub - 1]) - 1];
      }
    }

    yDest[static_cast<int32_T>(inferredChs_data[deltaSign]) - 1] = -holdPoint /
      static_cast<real_T>(inferredChs_size[1]);
  }

  // 'trajGen_:18' for i=1:size(y,1)
  loop_ub = static_cast<int32_T>(*trajectorySize);
  j_data[0] = 1;
  j_data[1] = 2;
  j_data[2] = 3;
  breaks[0] = -1.0;
  breaks[1] = 0.0;
  for (i = 0; i < 6; i++) {
    real_T delta1;
    real_T s0;
    real_T sF;
    real_T segATime;
    real_T segAcc;
    real_T yDest_0;
    real_T y__0;
    yDest_0 = yDest[i];
    y__0 = y_[i];

    // 'trajGen_:19' traj(i, 1:numWaypts) = trapveltraj(...
    // 'trajGen_:20'         [y(i), yDest(i)],...
    // 'trajGen_:21'         numWaypts,...
    // 'trajGen_:22'         PeakVelocity=abs( y(i) - yDest(i) ) / event.moveT); 
    delta1 = std::abs(y__0 - yDest_0) / event->moveT;
    if (static_cast<int32_T>(*trajectorySize) - 1 >= 0) {
      (void)std::memset(&rtDW.e_data[0], 0, static_cast<uint32_T>
                        (*trajectorySize) * sizeof(real_T));
    }

    (void)std::memset(&coeffsCell.f1[0], 0, 9U * sizeof(real_T));
    s0 = y__0;
    sF = yDest_0;
    deltaSign = 1;
    if (yDest_0 < y__0) {
      s0 = yDest_0;
      sF = y__0;
      deltaSign = -1;
    }

    holdPoint = (sF - s0) * 1.5 / delta1;
    segATime = ((s0 - sF) + delta1 * holdPoint) / delta1;
    segAcc = delta1 / segATime;
    if (s0 == sF) {
      segAcc = 0.0;
      delta1 = 0.0;
      if (std::isnan(holdPoint)) {
        holdPoint = 1.0;
      } else if (holdPoint == 0.0) {
        holdPoint = 1.0;
      } else {
        // no actions
      }

      segATime = holdPoint / 3.0;
    }

    delta1 *= static_cast<real_T>(deltaSign);
    segAcc *= static_cast<real_T>(deltaSign);
    (void)std::memset(&coefs[0], 0, 9U * sizeof(real_T));
    if (delta1 == 0.0) {
      coefs[6] = y__0;
      coefs[7] = y__0;
      coefs[8] = y__0;
    } else {
      coefs[0] = segAcc / 2.0;
      coefs[3] = 0.0;
      coefs[6] = y__0;
      coefs[1] = 0.0;
      coefs[4] = delta1;
      s0 = segAcc / 2.0 * (segATime * segATime);
      coefs[7] = s0 + y__0;
      coefs[2] = -segAcc / 2.0;
      coefs[5] = delta1;
      coefs[8] = (s0 + yDest_0) - delta1 * segATime;
    }

    for (deltaSign = 0; deltaSign < 3; deltaSign++) {
      for (directChs_size = 0; directChs_size < 3; directChs_size++) {
        coeffsCell.f1[(static_cast<int32_T>(j_data[directChs_size]) + 3 *
                       deltaSign) - 1] = coefs[3 * deltaSign + directChs_size];
      }
    }

    t_size[0] = 1;
    t_size[1] = static_cast<int32_T>(*trajectorySize);
    if (*trajectorySize >= 1UL) {
      rtDW.t_data[static_cast<int32_T>(*trajectorySize) - 1] = holdPoint;
      if (*trajectorySize >= 2UL) {
        rtDW.t_data[0] = 0.0;
        if (*trajectorySize >= 3UL) {
          if (-holdPoint == 0.0) {
            delta1 = holdPoint / (static_cast<real_T>(*trajectorySize) - 1.0);
            deltaSign = static_cast<int32_T>(*trajectorySize) - 1;
            for (directChs_size = 2; directChs_size <= deltaSign; directChs_size
                 ++) {
              rtDW.t_data[directChs_size - 1] = (static_cast<real_T>(
                static_cast<int32_T>((directChs_size << 1UL) -
                static_cast<int32_T>(*trajectorySize))) - 1.0) * delta1;
            }

            if ((static_cast<int32_T>(*trajectorySize) & 1) == 1) {
              rtDW.t_data[*trajectorySize >> 1UL] = 0.0;
            }
          } else {
            boolean_T guard1{ false };

            guard1 = false;
            if (holdPoint < 0.0) {
              if (std::abs(holdPoint) > 8.9884656743115785E+307) {
                delta1 = holdPoint / (static_cast<real_T>(*trajectorySize) - 1.0);
                deltaSign = static_cast<int32_T>(*trajectorySize);
                for (directChs_size = 0; directChs_size <= deltaSign - 3;
                     directChs_size++) {
                  rtDW.t_data[directChs_size + 1] = (static_cast<real_T>
                    (directChs_size) + 1.0) * delta1;
                }
              } else {
                guard1 = true;
              }
            } else {
              guard1 = true;
            }

            if (guard1) {
              delta1 = holdPoint / (static_cast<real_T>(*trajectorySize) - 1.0);
              deltaSign = static_cast<int32_T>(*trajectorySize);
              for (directChs_size = 0; directChs_size <= deltaSign - 3;
                   directChs_size++) {
                rtDW.t_data[directChs_size + 1] = (static_cast<real_T>
                  (directChs_size) + 1.0) * delta1;
              }
            }
          }
        }
      }
    }

    newSegmentCoeffs[0] = 0.0;
    newSegmentCoeffs[1] = 0.0;
    newSegmentCoeffs[2] = (coeffsCell.f1[0] * 0.0 + coeffsCell.f1[3] * 0.0) +
      coeffsCell.f1[6];
    breaks[2] = segATime;
    breaks[3] = holdPoint - segATime;
    breaks[4] = holdPoint;
    holdPoint -= breaks[3];
    yDest_0 = 0.0;
    deltaSign = 0;
    directChs_size = 0;
    for (int32_T deltaSign_0{0}; deltaSign_0 < 3; deltaSign_0++) {
      b_coefs[deltaSign] = newSegmentCoeffs[deltaSign_0];
      b_coefs[deltaSign + 1] = coeffsCell.f1[directChs_size];
      b_coefs[deltaSign + 2] = coeffsCell.f1[directChs_size + 1];
      b_coefs[deltaSign + 3] = coeffsCell.f1[directChs_size + 2];
      newSegmentCoeffs[deltaSign_0] = 0.0;
      yDest_0 += rt_powd_snf(holdPoint, 3.0 - (static_cast<real_T>(deltaSign_0)
        + 1.0)) * b_coefs[deltaSign + 3];
      deltaSign += 4;
      directChs_size += 3;
    }

    newSegmentCoeffs[2] = yDest_0;
    (void)std::memset(&expl_temp.coefs[0], 0, 15U * sizeof(real_T));
    deltaSign = 0;
    directChs_size = 0;
    for (int32_T deltaSign_0{0}; deltaSign_0 < 3; deltaSign_0++) {
      expl_temp.coefs[deltaSign] = b_coefs[directChs_size];
      expl_temp.coefs[deltaSign + 1] = b_coefs[directChs_size + 1];
      expl_temp.coefs[deltaSign + 2] = b_coefs[directChs_size + 2];
      expl_temp.coefs[deltaSign + 3] = b_coefs[directChs_size + 3];
      expl_temp.coefs[deltaSign + 4] = newSegmentCoeffs[deltaSign_0];
      deltaSign += 5;
      directChs_size += 4;
    }

    for (deltaSign = 0; deltaSign < 5; deltaSign++) {
      expl_temp.breaks[deltaSign] = breaks[deltaSign];
    }

    expl_temp.breaks[5] = breaks[4] + 1.0;
    ppval(&expl_temp, rtDW.t_data, t_size, rtDW.tmp_data, inferredChs_size);
    directChs_size = inferredChs_size[1];
    if (directChs_size - 1 >= 0) {
      (void)std::memcpy(&rtDW.e_data[0], &rtDW.tmp_data[0], static_cast<uint32_T>
                        (directChs_size) * sizeof(real_T));
    }

    for (deltaSign = 0; deltaSign < loop_ub; deltaSign++) {
      trajectory[i + 6 * deltaSign] = rtDW.e_data[deltaSign];
    }
  }
}

//
// Function for Chart: '<Root>/SupervisoryController'
// function [eventDone, waypt_, holdT] = handleEvent(event, holdT)
//
void SupervisoryController::handleEvent(real_T *holdT, boolean_T *eventDone,
  uint16_T *waypt_) const
{
  real_T b_holdT;

  // MATLAB Function 'handleEvent': '<S1>:532'
  // '<S1>:532:2' [eventDone, waypt_, holdT] = handleEvent_(event, waypt, holdT, y, ymax, trajSize, dt); 
  *waypt_ = rtDW.waypt;
  b_holdT = *holdT;

  //  initialize flags
  // 'handleEvent_:4' evDone = false;
  *eventDone = false;

  //  increment waypoint
  // 'handleEvent_:7' if waypt < numWaypts
  if (rtDW.waypt < rtDW.trajSize) {
    // 'handleEvent_:8' waypt = waypt + 1;
    *waypt_ = static_cast<uint16_T>(static_cast<uint32_T>(rtDW.waypt) + 1U);
  } else {
    // 'handleEvent_:9' else
    // 'handleEvent_:10' postT = postT - dt;
    b_holdT = *holdT - rtP.dt;
  }

  // 'handleEvent_:13' if postT < 0
  if (b_holdT < 0.0) {
    // 'handleEvent_:14' evDone = true;
    *eventDone = true;
  }

  *holdT = b_holdT;
}

// Function for MATLAB Function: '<S8>/MATLAB Function'
boolean_T SupervisoryController::any(const real_T x[3])
{
  int32_T k;
  boolean_T exitg1;
  boolean_T y;
  y = false;
  k = 0;
  exitg1 = false;
  while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
          false) && (k < 3)) {
    if ((!(x[k] == 0.0)) && (!std::isnan(x[k]))) {
      y = true;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

//
// Function for Chart: '<Root>/SupervisoryController'
// function sig_ = gainSchSig(ywt_)
//
real_T SupervisoryController::gainSchSig(const real_T ywt_[6])
{
  static const int8_T b[5]{ 0, 1, 2, 4, 5 };

  real_T ywt[6];
  real_T sig_;
  int32_T k;
  boolean_T x[5];
  boolean_T b_x[4];
  boolean_T exitg1;
  boolean_T guard1{ false };

  boolean_T guard2{ false };

  boolean_T guard3{ false };

  boolean_T guard4{ false };

  boolean_T guard5{ false };

  boolean_T y;

  // MATLAB Function 'gainSchSig': '<S1>:907'
  // '<S1>:907:2' sig_ = gainSchSig_(ywt_/beta);
  for (k = 0; k <= 4; k += 2) {
    (void)_mm_storeu_pd(&ywt[k], _mm_div_pd(_mm_loadu_pd(&ywt_[k]), _mm_set1_pd
      (rtP.beta)));
  }

  // 'gainSchSig_:3' if isempty(sigPrev)
  // 'gainSchSig_:7' if ( ywt(1) > 0.5 && all(ywt(2:end) < 0.5) )...
  // 'gainSchSig_:8'         || ( ywt(4) > 0.5 && all(ywt([1:3,5:end]) < 0.5) )
  guard1 = false;
  guard2 = false;
  guard3 = false;
  guard4 = false;
  guard5 = false;
  if (ywt[0] > 0.5) {
    for (k = 0; k < 5; k++) {
      x[k] = (ywt[k + 1] < 0.5);
    }

    y = true;
    k = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (k <= 4)) {
      if (!x[k]) {
        y = false;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (y) {
      // 'gainSchSig_:9' sig = 1;
      sig_ = 1.0;
    } else {
      guard5 = true;
    }
  } else {
    guard5 = true;
  }

  if (guard5) {
    if (ywt[3] > 0.5) {
      for (k = 0; k < 5; k++) {
        x[k] = (ywt[b[k]] < 0.5);
      }

      y = true;
      k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= 4)) {
        if (!x[k]) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }

      if (y) {
        // 'gainSchSig_:9' sig = 1;
        sig_ = 1.0;
      } else {
        guard4 = true;
      }
    } else {
      guard4 = true;
    }
  }

  if (guard4) {
    if ((ywt[1] > 0.5) && (ywt[2] > 0.5)) {
      b_x[0] = (ywt[0] < 0.5);
      b_x[1] = (ywt[3] < 0.5);
      b_x[2] = (ywt[4] < 0.5);
      b_x[3] = (ywt[5] < 0.5);
      y = true;
      k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= 3)) {
        if (!b_x[k]) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }

      if (y) {
        // 'gainSchSig_:10' elseif ( ywt(2) > 0.5 && ywt(3) > 0.5 && all(ywt([1,4:end]) < 0.5) )... 
        // 'gainSchSig_:11'         || ( ywt(5) > 0.5 && ywt(6) > 0.5 && all(ywt([1:4]) < 0.5) ) 
        // 'gainSchSig_:12' sig = 2;
        sig_ = 2.0;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }
  }

  if (guard3) {
    if ((ywt[4] > 0.5) && (ywt[5] > 0.5)) {
      b_x[0] = (ywt[0] < 0.5);
      b_x[1] = (ywt[1] < 0.5);
      b_x[2] = (ywt[2] < 0.5);
      b_x[3] = (ywt[3] < 0.5);
      y = true;
      k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= 3)) {
        if (!b_x[k]) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }

      if (y) {
        // 'gainSchSig_:10' elseif ( ywt(2) > 0.5 && ywt(3) > 0.5 && all(ywt([1,4:end]) < 0.5) )... 
        // 'gainSchSig_:11'         || ( ywt(5) > 0.5 && ywt(6) > 0.5 && all(ywt([1:4]) < 0.5) ) 
        // 'gainSchSig_:12' sig = 2;
        sig_ = 2.0;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }
  }

  if (guard2) {
    if ((ywt[0] > 0.5) && (ywt[1] > 0.5)) {
      b_x[0] = (ywt[2] < 0.5);
      b_x[1] = (ywt[3] < 0.5);
      b_x[2] = (ywt[4] < 0.5);
      b_x[3] = (ywt[5] < 0.5);
      y = true;
      k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= 3)) {
        if (!b_x[k]) {
          y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }

      if (y) {
        // 'gainSchSig_:13' elseif ( ywt(1) > 0.5 && ywt(2) > 0.5 && all(ywt([3:end]) < 0.5) )... 
        // 'gainSchSig_:14'         || ( ywt(4) > 0.5 && ywt(5) > 0.5 && all(ywt([1:3,6]) < 0.5) ) 
        // 'gainSchSig_:15' sig = 3;
        sig_ = 3.0;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    if (ywt[3] > 0.5) {
      if (ywt[4] > 0.5) {
        b_x[0] = (ywt[0] < 0.5);
        b_x[1] = (ywt[1] < 0.5);
        b_x[2] = (ywt[2] < 0.5);
        b_x[3] = (ywt[5] < 0.5);
        y = true;
        k = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (k <= 3)) {
          if (!b_x[k]) {
            y = false;
            exitg1 = true;
          } else {
            k++;
          }
        }

        if (y) {
          // 'gainSchSig_:13' elseif ( ywt(1) > 0.5 && ywt(2) > 0.5 && all(ywt([3:end]) < 0.5) )... 
          // 'gainSchSig_:14'         || ( ywt(4) > 0.5 && ywt(5) > 0.5 && all(ywt([1:3,6]) < 0.5) ) 
          // 'gainSchSig_:15' sig = 3;
          sig_ = 3.0;
        } else {
          // 'gainSchSig_:16' else
          // 'gainSchSig_:17' sig = sigPrev;
          sig_ = rtDW.sigPrev;
        }
      } else {
        // 'gainSchSig_:16' else
        // 'gainSchSig_:17' sig = sigPrev;
        sig_ = rtDW.sigPrev;
      }
    } else {
      // 'gainSchSig_:16' else
      // 'gainSchSig_:17' sig = sigPrev;
      sig_ = rtDW.sigPrev;
    }
  }

  // 'gainSchSig_:19' sigPrev = sig;
  rtDW.sigPrev = sig_;
  return sig_;
}

// Function for MATLAB Function: '<S110>/optimizer'
real_T SupervisoryController::norm(const real_T x[4])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[3]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

// Function for MATLAB Function: '<S110>/optimizer'
real_T SupervisoryController::maximum(const real_T x[4])
{
  real_T ex;
  int32_T idx;
  int32_T k;
  if (!std::isnan(x[0])) {
    idx = 1;
  } else {
    boolean_T exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (k < 5)) {
      if (!std::isnan(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    for (k = idx + 1; k < 5; k++) {
      real_T tmp;
      tmp = x[k - 1];
      if (ex < tmp) {
        ex = tmp;
      }
    }
  }

  return ex;
}

// Function for MATLAB Function: '<S110>/optimizer'
real_T SupervisoryController::xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (int32_T k{ix0}; k <= kend; k++) {
        real_T absxk;
        absxk = std::abs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * std::sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = std::sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = std::sqrt(b * b + 1.0) * a;
  } else if (std::isnan(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S110>/optimizer'
void SupervisoryController::xgemv(int32_T b_m, int32_T n, const real_T b_A[16],
  int32_T ia0, const real_T x[16], int32_T ix0, real_T y[4])
{
  if ((b_m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      (void)std::memset(&y[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    b = ((n - 1) << 2UL) + ia0;
    for (int32_T b_iy{ia0}; b_iy <= b; b_iy += 4) {
      real_T c;
      int32_T d;
      int32_T iyend;
      c = 0.0;
      d = b_iy + b_m;
      for (iyend = b_iy; iyend < d; iyend++) {
        c += x[((ix0 + iyend) - b_iy) - 1] * b_A[iyend - 1];
      }

      iyend = (b_iy - ia0) >> 2UL;
      y[iyend] += c;
    }
  }
}

// Function for MATLAB Function: '<S110>/optimizer'
void SupervisoryController::xgerc(int32_T b_m, int32_T n, real_T alpha1, int32_T
  ix0, const real_T y[4], real_T b_A[16], int32_T ia0)
{
  if (!(alpha1 == 0.0)) {
    int32_T jA;
    jA = ia0;
    for (int32_T j{0}; j < n; j++) {
      real_T temp;
      temp = y[j];
      if (temp != 0.0) {
        int32_T b;
        temp *= alpha1;
        b = b_m + jA;
        for (int32_T ijA{jA}; ijA < b; ijA++) {
          b_A[ijA - 1] += b_A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 4;
    }
  }
}

// Function for MATLAB Function: '<S250>/optimizer'
void SupervisoryController::KWIKfactor_ow(const real_T b_Ac[504], const int32_T
  iC[126], int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16],
  int32_T n, real_T RLinv[16], real_T *Status)
{
  __m128d tmp;
  real_T Q[16];
  real_T R[16];
  real_T TL[16];
  real_T b_A[16];
  real_T tau[4];
  real_T work[4];
  int32_T b_coltop;
  int32_T b_lastv;
  int32_T coltop;
  int32_T exitg1;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  int32_T scalarLB;
  int32_T vectorUB;
  boolean_T exitg2;
  *Status = 1.0;
  (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
  for (k_i = 0; k_i < nA; k_i++) {
    b_lastv = iC[k_i];
    for (b_coltop = 0; b_coltop < 4; b_coltop++) {
      knt = (k_i << 2UL) + b_coltop;
      RLinv[knt] = 0.0;
      RLinv[knt] += b_Ac[b_lastv - 1] * b_Linv[b_coltop];
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 125];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 251];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 377];
    }
  }

  (void)std::memcpy(&b_A[0], &RLinv[0], sizeof(real_T) << 4UL);
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  tau[3] = 0.0;
  work[3] = 0.0;
  for (k_i = 0; k_i < 4; k_i++) {
    ii = (k_i << 2UL) + k_i;
    if (k_i + 1 < 4) {
      real_T atmp;
      real_T beta1;
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = xnrm2(3 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        beta1 = rt_hypotd_snf(b_A[ii], beta1);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          coltop = (ii - k_i) + 4;
          do {
            knt++;
            scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
            vectorUB = scalarLB - 2;
            for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
              tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
              (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp,
                _mm_set1_pd(9.9792015476736E+291)));
            }

            for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
              b_A[b_coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((std::abs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, xnrm2(3 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          coltop = (ii - k_i) + 4;
          scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
            tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
            (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (atmp)));
          }

          for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
            b_A[b_coltop - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A[ii]) / beta1;
          atmp = 1.0 / (b_A[ii] - beta1);
          b_coltop = (ii - k_i) + 4;
          scalarLB = (((((b_coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (knt = b_lastv; knt <= vectorUB; knt += 2) {
            tmp = _mm_loadu_pd(&b_A[knt - 1]);
            (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
          }

          for (knt = scalarLB; knt <= b_coltop; knt++) {
            b_A[knt - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
      beta1 = b_A[ii];
      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        knt = (ii - k_i) + 3;
        while ((b_lastv > 0) && (b_A[knt] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (knt > 0)) {
          b_coltop = (((knt - 1) << 2UL) + ii) + 4;
          coltop = b_coltop;
          do {
            exitg1 = 0;
            if (coltop + 1 <= b_coltop + b_lastv) {
              if (b_A[coltop] != 0.0) {
                exitg1 = 1;
              } else {
                coltop++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        xgemv(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        xgerc(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
      }

      b_A[ii] = beta1;
    } else {
      tau[3] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      R[ii + (k_i << 2UL)] = b_A[(k_i << 2UL) + ii];
    }

    for (ii = k_i + 2; ii < 5; ii++) {
      R[(ii + (k_i << 2UL)) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 3; k_i >= 0; k_i--) {
    b_lastv = ((k_i << 2UL) + k_i) + 5;
    if (k_i + 1 < 4) {
      b_A[b_lastv - 5] = 1.0;
      if (tau[k_i] != 0.0) {
        knt = 4 - k_i;
        b_coltop = b_lastv - k_i;
        while ((knt > 0) && (b_A[b_coltop - 2] == 0.0)) {
          knt--;
          b_coltop--;
        }

        b_coltop = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (b_coltop > 0)) {
          coltop = ((b_coltop - 1) << 2UL) + b_lastv;
          ii = coltop;
          do {
            exitg1 = 0;
            if (ii <= (coltop + knt) - 1) {
              if (b_A[ii - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ii++;
              }
            } else {
              b_coltop--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        knt = 0;
        b_coltop = 0;
      }

      if (knt > 0) {
        xgemv(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 4, work);
        xgerc(knt, b_coltop, -tau[k_i], b_lastv - 4, work, b_A, b_lastv);
      }

      b_coltop = (b_lastv - k_i) - 1;
      scalarLB = (((((b_coltop - b_lastv) + 4) / 2) << 1UL) + b_lastv) - 3;
      vectorUB = scalarLB - 2;
      for (knt = b_lastv - 3; knt <= vectorUB; knt += 2) {
        tmp = _mm_loadu_pd(&b_A[knt - 1]);
        (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(-tau[k_i])));
      }

      for (knt = scalarLB; knt <= b_coltop; knt++) {
        b_A[knt - 1] *= -tau[k_i];
      }
    }

    b_A[b_lastv - 5] = 1.0 - tau[k_i];
    for (knt = 0; knt < k_i; knt++) {
      b_A[(b_lastv - knt) - 6] = 0.0;
    }
  }

  knt = 0;
  for (k_i = 0; k_i < 4; k_i++) {
    Q[knt] = b_A[knt];
    Q[knt + 1] = b_A[knt + 1];
    Q[knt + 2] = b_A[knt + 2];
    Q[knt + 3] = b_A[knt + 3];
    knt += 4;
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (std::abs(R[(k_i << 2UL) + k_i]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = 0;
        for (ii = 0; ii < n; ii++) {
          TL[coltop + k_i] = ((b_Linv[knt + 1] * Q[coltop + 1] + b_Linv[knt] *
                               Q[coltop]) + b_Linv[knt + 2] * Q[coltop + 2]) +
            b_Linv[knt + 3] * Q[coltop + 3];
          coltop += 4;
        }

        knt += 4;
      }

      (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
      for (k_i = nA; k_i >= 1; k_i--) {
        b_coltop = (k_i - 1) << 2UL;
        knt = (k_i + b_coltop) - 1;
        RLinv[knt] = 1.0;
        for (ii = k_i; ii <= nA; ii++) {
          coltop = (((ii - 1) << 2UL) + k_i) - 1;
          RLinv[coltop] /= R[knt];
        }

        if (k_i > 1) {
          for (ii = 0; ii <= k_i - 2; ii++) {
            for (b_lastv = k_i; b_lastv <= nA; b_lastv++) {
              knt = (b_lastv - 1) << 2UL;
              coltop = knt + ii;
              RLinv[coltop] -= RLinv[(knt + k_i) - 1] * R[b_coltop + ii];
            }
          }
        }
      }

      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = (k_i + 1) << 2UL;
        for (ii = k_i + 1; ii <= n; ii++) {
          b_coltop = (coltop + k_i) - 4;
          b_H[b_coltop] = 0.0;
          scalarLB = (nA + 1) << 2UL;
          for (b_lastv = nA + 1; b_lastv <= n; b_lastv++) {
            b_H[b_coltop] -= TL[(scalarLB + ii) - 5] * TL[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }

          b_H[(ii + knt) - 1] = b_H[b_coltop];
          coltop += 4;
        }

        knt += 4;
      }

      knt = 0;
      for (k_i = 0; k_i < nA; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_coltop = ii + knt;
          D[b_coltop] = 0.0;
          scalarLB = (k_i + 1) << 2UL;
          for (b_lastv = k_i + 1; b_lastv <= nA; b_lastv++) {
            D[b_coltop] += TL[(scalarLB + ii) - 4] * RLinv[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }
        }

        knt += 4;
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S250>/optimizer'
void SupervisoryController::DropConstraint_m(int32_T kDrop, boolean_T iA[126],
  int32_T *nA, int32_T iC[126])
{
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      for (int32_T i{kDrop}; i < *nA; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

// Function for MATLAB Function: '<S250>/optimizer'
void SupervisoryController::qpkwik_f(const real_T b_Linv[16], const real_T
  b_Hinv[16], const real_T f[4], const real_T b_Ac[504], const real_T b[126],
  boolean_T iA[126], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
  lambda[126], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[126];
  real_T D[16];
  real_T RLinv[16];
  real_T U[16];
  real_T b_H[16];
  real_T Opt[8];
  real_T Rhs[8];
  real_T r[4];
  real_T z[4];
  real_T Xnorm0;
  real_T cMin;
  real_T rMin;
  int32_T iC[126];
  int32_T b_exponent;
  int32_T exponent;
  int32_T i;
  int32_T iC_0;
  int32_T iSave;
  int32_T nA;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T guard1{ false };

  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 126; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 126; tmp++) {
    if (iA[tmp]) {
      nA++;
      iC[nA - 1] = tmp + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    int32_T exitg3;
    (void)std::memset(&Opt[0], 0, sizeof(real_T) << 3UL);
    Rhs[0] = f[0];
    Rhs[4] = 0.0;
    Rhs[1] = f[1];
    Rhs[5] = 0.0;
    Rhs[2] = f[2];
    Rhs[6] = 0.0;
    Rhs[3] = f[3];
    Rhs[7] = 0.0;
    DualFeasible = false;
    tmp = static_cast<int32_T>(std::round(0.3 * static_cast<real_T>(nA)));
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        KWIKfactor_ow(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            (void)std::memset(&iA[0], 0, 126U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 126U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          int32_T U_tmp;
          for (i = 0; i < nA; i++) {
            Rhs[i + 4] = b[iC[i] - 1];
            for (iSave = i + 1; iSave <= nA; iSave++) {
              U[(iSave + (i << 2UL)) - 1] = 0.0;
              for (iC_0 = 0; iC_0 < nA; iC_0++) {
                int32_T U_tmp_0;
                U_tmp = iC_0 << 2UL;
                U_tmp_0 = ((i << 2UL) + iSave) - 1;
                U[U_tmp_0] += RLinv[(U_tmp + iSave) - 1] * RLinv[U_tmp + i];
              }

              U[i + ((iSave - 1) << 2UL)] = U[((i << 2UL) + iSave) - 1];
            }
          }

          for (i = 0; i < 4; i++) {
            Opt[i] = ((b_H[i + 4] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 8] *
                      Rhs[2]) + b_H[i + 12] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i] += D[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }
          }

          U_tmp = 0;
          for (i = 0; i < nA; i++) {
            Opt[i + 4] = ((D[U_tmp + 1] * Rhs[1] + D[U_tmp] * Rhs[0]) + D[U_tmp
                          + 2] * Rhs[2]) + D[U_tmp + 3] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i + 4] += U[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }

            U_tmp += 4;
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (iSave = 0; iSave < nA; iSave++) {
            lambda[iC[iSave] - 1] = Opt[iSave + 4];
            cMin = Opt[iSave + 4];
            if ((cMin < Xnorm0) && (iSave + 1 <= nA)) {
              i = iSave;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
            x[3] = Opt[3];
          } else {
            (*status)++;
            if (tmp <= 5) {
              iC_0 = 5;
            } else {
              iC_0 = tmp;
            }

            if (*status > iC_0) {
              nA = 0;
              (void)std::memset(&iA[0], 0, 126U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 126U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint_m(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 126U * sizeof(real_T));
          for (tmp = 0; tmp <= 2; tmp += 2) {
            tmp_3 = _mm_set1_pd(-1.0);
            (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd
              (_mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3),
                          _mm_set1_pd(f[1])), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
              (&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f[0]))), _mm_mul_pd(_mm_mul_pd
              (_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3), _mm_set1_pd(f[2]))),
              _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3),
                         _mm_set1_pd(f[3]))));
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    for (tmp = 0; tmp <= 2; tmp += 2) {
      tmp_3 = _mm_set1_pd(-1.0);
      (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
        (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3), _mm_set1_pd(f[1])),
        _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f
        [0]))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3),
                           _mm_set1_pd(f[2]))), _mm_mul_pd(_mm_mul_pd
        (_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3), _mm_set1_pd(f[3]))));
    }

    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    Xnorm0 = norm(x);
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (*status <= maxiter)) {
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 126; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 126] * x[1]);
          z[2] = std::abs(b_Ac[i + 252] * x[2]);
          z[3] = std::abs(b_Ac[i + 378] * x[3]);
          t = std::fmax(t, maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 126] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 252] * x
                    [2]) + b_Ac[i + 378] * x[3]) - b[i]) / t;
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }

        cTol[i] = t;
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        int32_T exitg1;
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            boolean_T guard2{ false };

            guard2 = false;
            if (nA == 0) {
              for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 378])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 252])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 126])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0]), _mm_set1_pd
                  (b_Ac[tmp])), _mm_set1_pd(0.0))))));
              }

              guard2 = true;
            } else {
              KWIKfactor_ow(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (iC_0 = 0; iC_0 <= 14; iC_0 += 2) {
                  tmp_3 = _mm_loadu_pd(&b_H[iC_0]);
                  (void)_mm_storeu_pd(&U[iC_0], _mm_mul_pd(tmp_3, _mm_set1_pd
                    (-1.0)));
                }

                for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                  __m128d tmp_0;
                  __m128d tmp_1;
                  __m128d tmp_2;
                  tmp_3 = _mm_loadu_pd(&U[iC_0]);
                  tmp_0 = _mm_loadu_pd(&U[iC_0 + 4]);
                  tmp_1 = _mm_loadu_pd(&U[iC_0 + 8]);
                  tmp_2 = _mm_loadu_pd(&U[iC_0 + 12]);
                  (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(tmp_2,
                    _mm_set1_pd(b_Ac[tmp + 378])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 252])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 126])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((D[iSave + 1] * b_Ac[tmp + 126] + D[iSave] * b_Ac[tmp])
                          + D[iSave + 2] * b_Ac[tmp + 252]) + D[iSave + 3] *
                    b_Ac[tmp + 378];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              real_T cVal_tmp;
              real_T cVal_tmp_0;
              boolean_T exitg4;
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                iSave = 0;
                exitg4 = false;
                while (((exitg4 ? static_cast<uint32_T>(1U) :
                         static_cast<uint32_T>(0U)) == false) && (iSave <= nA -
                        1)) {
                  if (r[iSave] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    iSave++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (iSave = 0; iSave < nA; iSave++) {
                  cVal = r[iSave];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[iSave] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = iSave + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              t = b_Ac[tmp + 126];
              cVal_tmp = b_Ac[tmp + 252];
              cVal_tmp_0 = b_Ac[tmp + 378];
              cVal = ((t * z[1] + z[0] * b_Ac[tmp]) + cVal_tmp * z[2]) +
                cVal_tmp_0 * z[3];
              if (cVal <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - (((t * x[1] + b_Ac[tmp] * x[0]) + cVal_tmp * x
                                   [2]) + cVal_tmp_0 * x[3])) / cVal;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (iSave = 0; iSave < nA; iSave++) {
                  iC_0 = iC[iSave];
                  lambda[iC_0 - 1] -= t * r[iSave];
                  if ((iC_0 <= 126) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                (void)std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint_m(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  (void)std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == static_cast<int32_T>(degrees)) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while (((exitg4 ? static_cast<uint32_T>(1U) : static_cast<
                               uint32_T>(0U)) == false) && (i + 1 > 1)) {
                        iC_0 = iC[i - 1];
                        if (iC[i] > iC_0) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = iC_0;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 126; tmp++) {
                cTol[tmp] = std::fmax(std::abs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S180>/optimizer'
void SupervisoryController::KWIKfactor_o(const real_T b_Ac[824], const int32_T
  iC[206], int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16],
  int32_T n, real_T RLinv[16], real_T *Status)
{
  __m128d tmp;
  real_T Q[16];
  real_T R[16];
  real_T TL[16];
  real_T b_A[16];
  real_T tau[4];
  real_T work[4];
  int32_T b_coltop;
  int32_T b_lastv;
  int32_T coltop;
  int32_T exitg1;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  int32_T scalarLB;
  int32_T vectorUB;
  boolean_T exitg2;
  *Status = 1.0;
  (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
  for (k_i = 0; k_i < nA; k_i++) {
    b_lastv = iC[k_i];
    for (b_coltop = 0; b_coltop < 4; b_coltop++) {
      knt = (k_i << 2UL) + b_coltop;
      RLinv[knt] = 0.0;
      RLinv[knt] += b_Ac[b_lastv - 1] * b_Linv[b_coltop];
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 205];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 411];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 617];
    }
  }

  (void)std::memcpy(&b_A[0], &RLinv[0], sizeof(real_T) << 4UL);
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  tau[3] = 0.0;
  work[3] = 0.0;
  for (k_i = 0; k_i < 4; k_i++) {
    ii = (k_i << 2UL) + k_i;
    if (k_i + 1 < 4) {
      real_T atmp;
      real_T beta1;
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = xnrm2(3 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        beta1 = rt_hypotd_snf(b_A[ii], beta1);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          coltop = (ii - k_i) + 4;
          do {
            knt++;
            scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
            vectorUB = scalarLB - 2;
            for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
              tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
              (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp,
                _mm_set1_pd(9.9792015476736E+291)));
            }

            for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
              b_A[b_coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((std::abs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, xnrm2(3 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          coltop = (ii - k_i) + 4;
          scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
            tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
            (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (atmp)));
          }

          for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
            b_A[b_coltop - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A[ii]) / beta1;
          atmp = 1.0 / (b_A[ii] - beta1);
          b_coltop = (ii - k_i) + 4;
          scalarLB = (((((b_coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (knt = b_lastv; knt <= vectorUB; knt += 2) {
            tmp = _mm_loadu_pd(&b_A[knt - 1]);
            (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
          }

          for (knt = scalarLB; knt <= b_coltop; knt++) {
            b_A[knt - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
      beta1 = b_A[ii];
      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        knt = (ii - k_i) + 3;
        while ((b_lastv > 0) && (b_A[knt] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (knt > 0)) {
          b_coltop = (((knt - 1) << 2UL) + ii) + 4;
          coltop = b_coltop;
          do {
            exitg1 = 0;
            if (coltop + 1 <= b_coltop + b_lastv) {
              if (b_A[coltop] != 0.0) {
                exitg1 = 1;
              } else {
                coltop++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        xgemv(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        xgerc(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
      }

      b_A[ii] = beta1;
    } else {
      tau[3] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      R[ii + (k_i << 2UL)] = b_A[(k_i << 2UL) + ii];
    }

    for (ii = k_i + 2; ii < 5; ii++) {
      R[(ii + (k_i << 2UL)) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 3; k_i >= 0; k_i--) {
    b_lastv = ((k_i << 2UL) + k_i) + 5;
    if (k_i + 1 < 4) {
      b_A[b_lastv - 5] = 1.0;
      if (tau[k_i] != 0.0) {
        knt = 4 - k_i;
        b_coltop = b_lastv - k_i;
        while ((knt > 0) && (b_A[b_coltop - 2] == 0.0)) {
          knt--;
          b_coltop--;
        }

        b_coltop = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (b_coltop > 0)) {
          coltop = ((b_coltop - 1) << 2UL) + b_lastv;
          ii = coltop;
          do {
            exitg1 = 0;
            if (ii <= (coltop + knt) - 1) {
              if (b_A[ii - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ii++;
              }
            } else {
              b_coltop--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        knt = 0;
        b_coltop = 0;
      }

      if (knt > 0) {
        xgemv(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 4, work);
        xgerc(knt, b_coltop, -tau[k_i], b_lastv - 4, work, b_A, b_lastv);
      }

      b_coltop = (b_lastv - k_i) - 1;
      scalarLB = (((((b_coltop - b_lastv) + 4) / 2) << 1UL) + b_lastv) - 3;
      vectorUB = scalarLB - 2;
      for (knt = b_lastv - 3; knt <= vectorUB; knt += 2) {
        tmp = _mm_loadu_pd(&b_A[knt - 1]);
        (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(-tau[k_i])));
      }

      for (knt = scalarLB; knt <= b_coltop; knt++) {
        b_A[knt - 1] *= -tau[k_i];
      }
    }

    b_A[b_lastv - 5] = 1.0 - tau[k_i];
    for (knt = 0; knt < k_i; knt++) {
      b_A[(b_lastv - knt) - 6] = 0.0;
    }
  }

  knt = 0;
  for (k_i = 0; k_i < 4; k_i++) {
    Q[knt] = b_A[knt];
    Q[knt + 1] = b_A[knt + 1];
    Q[knt + 2] = b_A[knt + 2];
    Q[knt + 3] = b_A[knt + 3];
    knt += 4;
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (std::abs(R[(k_i << 2UL) + k_i]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = 0;
        for (ii = 0; ii < n; ii++) {
          TL[coltop + k_i] = ((b_Linv[knt + 1] * Q[coltop + 1] + b_Linv[knt] *
                               Q[coltop]) + b_Linv[knt + 2] * Q[coltop + 2]) +
            b_Linv[knt + 3] * Q[coltop + 3];
          coltop += 4;
        }

        knt += 4;
      }

      (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
      for (k_i = nA; k_i >= 1; k_i--) {
        b_coltop = (k_i - 1) << 2UL;
        knt = (k_i + b_coltop) - 1;
        RLinv[knt] = 1.0;
        for (ii = k_i; ii <= nA; ii++) {
          coltop = (((ii - 1) << 2UL) + k_i) - 1;
          RLinv[coltop] /= R[knt];
        }

        if (k_i > 1) {
          for (ii = 0; ii <= k_i - 2; ii++) {
            for (b_lastv = k_i; b_lastv <= nA; b_lastv++) {
              knt = (b_lastv - 1) << 2UL;
              coltop = knt + ii;
              RLinv[coltop] -= RLinv[(knt + k_i) - 1] * R[b_coltop + ii];
            }
          }
        }
      }

      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = (k_i + 1) << 2UL;
        for (ii = k_i + 1; ii <= n; ii++) {
          b_coltop = (coltop + k_i) - 4;
          b_H[b_coltop] = 0.0;
          scalarLB = (nA + 1) << 2UL;
          for (b_lastv = nA + 1; b_lastv <= n; b_lastv++) {
            b_H[b_coltop] -= TL[(scalarLB + ii) - 5] * TL[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }

          b_H[(ii + knt) - 1] = b_H[b_coltop];
          coltop += 4;
        }

        knt += 4;
      }

      knt = 0;
      for (k_i = 0; k_i < nA; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_coltop = ii + knt;
          D[b_coltop] = 0.0;
          scalarLB = (k_i + 1) << 2UL;
          for (b_lastv = k_i + 1; b_lastv <= nA; b_lastv++) {
            D[b_coltop] += TL[(scalarLB + ii) - 4] * RLinv[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }
        }

        knt += 4;
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S180>/optimizer'
void SupervisoryController::DropConstraint_f(int32_T kDrop, boolean_T iA[206],
  int32_T *nA, int32_T iC[206])
{
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      for (int32_T i{kDrop}; i < *nA; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

// Function for MATLAB Function: '<S180>/optimizer'
void SupervisoryController::qpkwik_o(const real_T b_Linv[16], const real_T
  b_Hinv[16], const real_T f[4], const real_T b_Ac[824], const real_T b[206],
  boolean_T iA[206], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
  lambda[206], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[206];
  real_T D[16];
  real_T RLinv[16];
  real_T U[16];
  real_T b_H[16];
  real_T Opt[8];
  real_T Rhs[8];
  real_T r[4];
  real_T z[4];
  real_T Xnorm0;
  real_T cMin;
  real_T rMin;
  int32_T iC[206];
  int32_T b_exponent;
  int32_T exponent;
  int32_T i;
  int32_T iC_0;
  int32_T iSave;
  int32_T nA;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T guard1{ false };

  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 206; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 206; tmp++) {
    if (iA[tmp]) {
      nA++;
      iC[nA - 1] = tmp + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    int32_T exitg3;
    (void)std::memset(&Opt[0], 0, sizeof(real_T) << 3UL);
    Rhs[0] = f[0];
    Rhs[4] = 0.0;
    Rhs[1] = f[1];
    Rhs[5] = 0.0;
    Rhs[2] = f[2];
    Rhs[6] = 0.0;
    Rhs[3] = f[3];
    Rhs[7] = 0.0;
    DualFeasible = false;
    tmp = static_cast<int32_T>(std::round(0.3 * static_cast<real_T>(nA)));
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        KWIKfactor_o(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            (void)std::memset(&iA[0], 0, 206U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 206U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          int32_T U_tmp;
          for (i = 0; i < nA; i++) {
            Rhs[i + 4] = b[iC[i] - 1];
            for (iSave = i + 1; iSave <= nA; iSave++) {
              U[(iSave + (i << 2UL)) - 1] = 0.0;
              for (iC_0 = 0; iC_0 < nA; iC_0++) {
                int32_T U_tmp_0;
                U_tmp = iC_0 << 2UL;
                U_tmp_0 = ((i << 2UL) + iSave) - 1;
                U[U_tmp_0] += RLinv[(U_tmp + iSave) - 1] * RLinv[U_tmp + i];
              }

              U[i + ((iSave - 1) << 2UL)] = U[((i << 2UL) + iSave) - 1];
            }
          }

          for (i = 0; i < 4; i++) {
            Opt[i] = ((b_H[i + 4] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 8] *
                      Rhs[2]) + b_H[i + 12] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i] += D[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }
          }

          U_tmp = 0;
          for (i = 0; i < nA; i++) {
            Opt[i + 4] = ((D[U_tmp + 1] * Rhs[1] + D[U_tmp] * Rhs[0]) + D[U_tmp
                          + 2] * Rhs[2]) + D[U_tmp + 3] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i + 4] += U[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }

            U_tmp += 4;
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (iSave = 0; iSave < nA; iSave++) {
            lambda[iC[iSave] - 1] = Opt[iSave + 4];
            cMin = Opt[iSave + 4];
            if ((cMin < Xnorm0) && (iSave + 1 <= nA)) {
              i = iSave;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
            x[3] = Opt[3];
          } else {
            (*status)++;
            if (tmp <= 5) {
              iC_0 = 5;
            } else {
              iC_0 = tmp;
            }

            if (*status > iC_0) {
              nA = 0;
              (void)std::memset(&iA[0], 0, 206U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 206U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint_f(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 206U * sizeof(real_T));
          for (tmp = 0; tmp <= 2; tmp += 2) {
            tmp_3 = _mm_set1_pd(-1.0);
            (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd
              (_mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3),
                          _mm_set1_pd(f[1])), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
              (&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f[0]))), _mm_mul_pd(_mm_mul_pd
              (_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3), _mm_set1_pd(f[2]))),
              _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3),
                         _mm_set1_pd(f[3]))));
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    for (tmp = 0; tmp <= 2; tmp += 2) {
      tmp_3 = _mm_set1_pd(-1.0);
      (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
        (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3), _mm_set1_pd(f[1])),
        _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f
        [0]))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3),
                           _mm_set1_pd(f[2]))), _mm_mul_pd(_mm_mul_pd
        (_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3), _mm_set1_pd(f[3]))));
    }

    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    Xnorm0 = norm(x);
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (*status <= maxiter)) {
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 206; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 206] * x[1]);
          z[2] = std::abs(b_Ac[i + 412] * x[2]);
          z[3] = std::abs(b_Ac[i + 618] * x[3]);
          t = std::fmax(t, maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 206] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 412] * x
                    [2]) + b_Ac[i + 618] * x[3]) - b[i]) / t;
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }

        cTol[i] = t;
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        int32_T exitg1;
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            boolean_T guard2{ false };

            guard2 = false;
            if (nA == 0) {
              for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 618])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 412])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 206])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0]), _mm_set1_pd
                  (b_Ac[tmp])), _mm_set1_pd(0.0))))));
              }

              guard2 = true;
            } else {
              KWIKfactor_o(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (iC_0 = 0; iC_0 <= 14; iC_0 += 2) {
                  tmp_3 = _mm_loadu_pd(&b_H[iC_0]);
                  (void)_mm_storeu_pd(&U[iC_0], _mm_mul_pd(tmp_3, _mm_set1_pd
                    (-1.0)));
                }

                for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                  __m128d tmp_0;
                  __m128d tmp_1;
                  __m128d tmp_2;
                  tmp_3 = _mm_loadu_pd(&U[iC_0]);
                  tmp_0 = _mm_loadu_pd(&U[iC_0 + 4]);
                  tmp_1 = _mm_loadu_pd(&U[iC_0 + 8]);
                  tmp_2 = _mm_loadu_pd(&U[iC_0 + 12]);
                  (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(tmp_2,
                    _mm_set1_pd(b_Ac[tmp + 618])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 412])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 206])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((D[iSave + 1] * b_Ac[tmp + 206] + D[iSave] * b_Ac[tmp])
                          + D[iSave + 2] * b_Ac[tmp + 412]) + D[iSave + 3] *
                    b_Ac[tmp + 618];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              real_T cVal_tmp;
              real_T cVal_tmp_0;
              boolean_T exitg4;
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                iSave = 0;
                exitg4 = false;
                while (((exitg4 ? static_cast<uint32_T>(1U) :
                         static_cast<uint32_T>(0U)) == false) && (iSave <= nA -
                        1)) {
                  if (r[iSave] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    iSave++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (iSave = 0; iSave < nA; iSave++) {
                  cVal = r[iSave];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[iSave] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = iSave + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              t = b_Ac[tmp + 206];
              cVal_tmp = b_Ac[tmp + 412];
              cVal_tmp_0 = b_Ac[tmp + 618];
              cVal = ((t * z[1] + z[0] * b_Ac[tmp]) + cVal_tmp * z[2]) +
                cVal_tmp_0 * z[3];
              if (cVal <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - (((t * x[1] + b_Ac[tmp] * x[0]) + cVal_tmp * x
                                   [2]) + cVal_tmp_0 * x[3])) / cVal;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (iSave = 0; iSave < nA; iSave++) {
                  iC_0 = iC[iSave];
                  lambda[iC_0 - 1] -= t * r[iSave];
                  if ((iC_0 <= 206) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                (void)std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint_f(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  (void)std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == static_cast<int32_T>(degrees)) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while (((exitg4 ? static_cast<uint32_T>(1U) : static_cast<
                               uint32_T>(0U)) == false) && (i + 1 > 1)) {
                        iC_0 = iC[i - 1];
                        if (iC[i] > iC_0) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = iC_0;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 206; tmp++) {
                cTol[tmp] = std::fmax(std::abs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S110>/optimizer'
void SupervisoryController::KWIKfactor(const real_T b_Ac[664], const int32_T iC
  [166], int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16],
  int32_T n, real_T RLinv[16], real_T *Status)
{
  __m128d tmp;
  real_T Q[16];
  real_T R[16];
  real_T TL[16];
  real_T b_A[16];
  real_T tau[4];
  real_T work[4];
  int32_T b_coltop;
  int32_T b_lastv;
  int32_T coltop;
  int32_T exitg1;
  int32_T ii;
  int32_T k_i;
  int32_T knt;
  int32_T scalarLB;
  int32_T vectorUB;
  boolean_T exitg2;
  *Status = 1.0;
  (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
  for (k_i = 0; k_i < nA; k_i++) {
    b_lastv = iC[k_i];
    for (b_coltop = 0; b_coltop < 4; b_coltop++) {
      knt = (k_i << 2UL) + b_coltop;
      RLinv[knt] = 0.0;
      RLinv[knt] += b_Ac[b_lastv - 1] * b_Linv[b_coltop];
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 165];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 331];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 497];
    }
  }

  (void)std::memcpy(&b_A[0], &RLinv[0], sizeof(real_T) << 4UL);
  tau[0] = 0.0;
  work[0] = 0.0;
  tau[1] = 0.0;
  work[1] = 0.0;
  tau[2] = 0.0;
  work[2] = 0.0;
  tau[3] = 0.0;
  work[3] = 0.0;
  for (k_i = 0; k_i < 4; k_i++) {
    ii = (k_i << 2UL) + k_i;
    if (k_i + 1 < 4) {
      real_T atmp;
      real_T beta1;
      atmp = b_A[ii];
      b_lastv = ii + 2;
      tau[k_i] = 0.0;
      beta1 = xnrm2(3 - k_i, b_A, ii + 2);
      if (beta1 != 0.0) {
        beta1 = rt_hypotd_snf(b_A[ii], beta1);
        if (b_A[ii] >= 0.0) {
          beta1 = -beta1;
        }

        if (std::abs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          coltop = (ii - k_i) + 4;
          do {
            knt++;
            scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
            vectorUB = scalarLB - 2;
            for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
              tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
              (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp,
                _mm_set1_pd(9.9792015476736E+291)));
            }

            for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
              b_A[b_coltop - 1] *= 9.9792015476736E+291;
            }

            beta1 *= 9.9792015476736E+291;
            atmp *= 9.9792015476736E+291;
          } while ((std::abs(beta1) < 1.0020841800044864E-292) && (knt < 20));

          beta1 = rt_hypotd_snf(atmp, xnrm2(3 - k_i, b_A, ii + 2));
          if (atmp >= 0.0) {
            beta1 = -beta1;
          }

          tau[k_i] = (beta1 - atmp) / beta1;
          atmp = 1.0 / (atmp - beta1);
          coltop = (ii - k_i) + 4;
          scalarLB = (((((coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (b_coltop = b_lastv; b_coltop <= vectorUB; b_coltop += 2) {
            tmp = _mm_loadu_pd(&b_A[b_coltop - 1]);
            (void)_mm_storeu_pd(&b_A[b_coltop - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (atmp)));
          }

          for (b_coltop = scalarLB; b_coltop <= coltop; b_coltop++) {
            b_A[b_coltop - 1] *= atmp;
          }

          for (b_lastv = 0; b_lastv < knt; b_lastv++) {
            beta1 *= 1.0020841800044864E-292;
          }

          atmp = beta1;
        } else {
          tau[k_i] = (beta1 - b_A[ii]) / beta1;
          atmp = 1.0 / (b_A[ii] - beta1);
          b_coltop = (ii - k_i) + 4;
          scalarLB = (((((b_coltop - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB - 2;
          for (knt = b_lastv; knt <= vectorUB; knt += 2) {
            tmp = _mm_loadu_pd(&b_A[knt - 1]);
            (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
          }

          for (knt = scalarLB; knt <= b_coltop; knt++) {
            b_A[knt - 1] *= atmp;
          }

          atmp = beta1;
        }
      }

      b_A[ii] = atmp;
      beta1 = b_A[ii];
      b_A[ii] = 1.0;
      if (tau[k_i] != 0.0) {
        b_lastv = 4 - k_i;
        knt = (ii - k_i) + 3;
        while ((b_lastv > 0) && (b_A[knt] == 0.0)) {
          b_lastv--;
          knt--;
        }

        knt = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (knt > 0)) {
          b_coltop = (((knt - 1) << 2UL) + ii) + 4;
          coltop = b_coltop;
          do {
            exitg1 = 0;
            if (coltop + 1 <= b_coltop + b_lastv) {
              if (b_A[coltop] != 0.0) {
                exitg1 = 1;
              } else {
                coltop++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        b_lastv = 0;
        knt = 0;
      }

      if (b_lastv > 0) {
        xgemv(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        xgerc(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
      }

      b_A[ii] = beta1;
    } else {
      tau[3] = 0.0;
    }
  }

  for (k_i = 0; k_i < 4; k_i++) {
    for (ii = 0; ii <= k_i; ii++) {
      R[ii + (k_i << 2UL)] = b_A[(k_i << 2UL) + ii];
    }

    for (ii = k_i + 2; ii < 5; ii++) {
      R[(ii + (k_i << 2UL)) - 1] = 0.0;
    }

    work[k_i] = 0.0;
  }

  for (k_i = 3; k_i >= 0; k_i--) {
    b_lastv = ((k_i << 2UL) + k_i) + 5;
    if (k_i + 1 < 4) {
      b_A[b_lastv - 5] = 1.0;
      if (tau[k_i] != 0.0) {
        knt = 4 - k_i;
        b_coltop = b_lastv - k_i;
        while ((knt > 0) && (b_A[b_coltop - 2] == 0.0)) {
          knt--;
          b_coltop--;
        }

        b_coltop = 3 - k_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (b_coltop > 0)) {
          coltop = ((b_coltop - 1) << 2UL) + b_lastv;
          ii = coltop;
          do {
            exitg1 = 0;
            if (ii <= (coltop + knt) - 1) {
              if (b_A[ii - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ii++;
              }
            } else {
              b_coltop--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        knt = 0;
        b_coltop = 0;
      }

      if (knt > 0) {
        xgemv(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 4, work);
        xgerc(knt, b_coltop, -tau[k_i], b_lastv - 4, work, b_A, b_lastv);
      }

      b_coltop = (b_lastv - k_i) - 1;
      scalarLB = (((((b_coltop - b_lastv) + 4) / 2) << 1UL) + b_lastv) - 3;
      vectorUB = scalarLB - 2;
      for (knt = b_lastv - 3; knt <= vectorUB; knt += 2) {
        tmp = _mm_loadu_pd(&b_A[knt - 1]);
        (void)_mm_storeu_pd(&b_A[knt - 1], _mm_mul_pd(tmp, _mm_set1_pd(-tau[k_i])));
      }

      for (knt = scalarLB; knt <= b_coltop; knt++) {
        b_A[knt - 1] *= -tau[k_i];
      }
    }

    b_A[b_lastv - 5] = 1.0 - tau[k_i];
    for (knt = 0; knt < k_i; knt++) {
      b_A[(b_lastv - knt) - 6] = 0.0;
    }
  }

  knt = 0;
  for (k_i = 0; k_i < 4; k_i++) {
    Q[knt] = b_A[knt];
    Q[knt + 1] = b_A[knt + 1];
    Q[knt + 2] = b_A[knt + 2];
    Q[knt + 3] = b_A[knt + 3];
    knt += 4;
  }

  k_i = 0;
  do {
    exitg1 = 0;
    if (k_i <= nA - 1) {
      if (std::abs(R[(k_i << 2UL) + k_i]) < 1.0E-12) {
        *Status = -2.0;
        exitg1 = 1;
      } else {
        k_i++;
      }
    } else {
      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = 0;
        for (ii = 0; ii < n; ii++) {
          TL[coltop + k_i] = ((b_Linv[knt + 1] * Q[coltop + 1] + b_Linv[knt] *
                               Q[coltop]) + b_Linv[knt + 2] * Q[coltop + 2]) +
            b_Linv[knt + 3] * Q[coltop + 3];
          coltop += 4;
        }

        knt += 4;
      }

      (void)std::memset(&RLinv[0], 0, sizeof(real_T) << 4UL);
      for (k_i = nA; k_i >= 1; k_i--) {
        b_coltop = (k_i - 1) << 2UL;
        knt = (k_i + b_coltop) - 1;
        RLinv[knt] = 1.0;
        for (ii = k_i; ii <= nA; ii++) {
          coltop = (((ii - 1) << 2UL) + k_i) - 1;
          RLinv[coltop] /= R[knt];
        }

        if (k_i > 1) {
          for (ii = 0; ii <= k_i - 2; ii++) {
            for (b_lastv = k_i; b_lastv <= nA; b_lastv++) {
              knt = (b_lastv - 1) << 2UL;
              coltop = knt + ii;
              RLinv[coltop] -= RLinv[(knt + k_i) - 1] * R[b_coltop + ii];
            }
          }
        }
      }

      knt = 0;
      for (k_i = 0; k_i < n; k_i++) {
        coltop = (k_i + 1) << 2UL;
        for (ii = k_i + 1; ii <= n; ii++) {
          b_coltop = (coltop + k_i) - 4;
          b_H[b_coltop] = 0.0;
          scalarLB = (nA + 1) << 2UL;
          for (b_lastv = nA + 1; b_lastv <= n; b_lastv++) {
            b_H[b_coltop] -= TL[(scalarLB + ii) - 5] * TL[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }

          b_H[(ii + knt) - 1] = b_H[b_coltop];
          coltop += 4;
        }

        knt += 4;
      }

      knt = 0;
      for (k_i = 0; k_i < nA; k_i++) {
        for (ii = 0; ii < n; ii++) {
          b_coltop = ii + knt;
          D[b_coltop] = 0.0;
          scalarLB = (k_i + 1) << 2UL;
          for (b_lastv = k_i + 1; b_lastv <= nA; b_lastv++) {
            D[b_coltop] += TL[(scalarLB + ii) - 4] * RLinv[(scalarLB + k_i) - 4];
            scalarLB += 4;
          }
        }

        knt += 4;
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S110>/optimizer'
void SupervisoryController::DropConstraint(int32_T kDrop, boolean_T iA[166],
  int32_T *nA, int32_T iC[166])
{
  if (kDrop > 0) {
    iA[iC[kDrop - 1] - 1] = false;
    if (kDrop < *nA) {
      for (int32_T i{kDrop}; i < *nA; i++) {
        iC[i - 1] = iC[i];
      }
    }

    iC[*nA - 1] = 0;
    (*nA)--;
  }
}

// Function for MATLAB Function: '<S110>/optimizer'
void SupervisoryController::qpkwik(const real_T b_Linv[16], const real_T b_Hinv
  [16], const real_T f[4], const real_T b_Ac[664], const real_T b[166],
  boolean_T iA[166], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
  lambda[166], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[166];
  real_T D[16];
  real_T RLinv[16];
  real_T U[16];
  real_T b_H[16];
  real_T Opt[8];
  real_T Rhs[8];
  real_T r[4];
  real_T z[4];
  real_T Xnorm0;
  real_T cMin;
  real_T rMin;
  int32_T iC[166];
  int32_T b_exponent;
  int32_T exponent;
  int32_T i;
  int32_T iC_0;
  int32_T iSave;
  int32_T nA;
  int32_T tmp;
  boolean_T ColdReset;
  boolean_T DualFeasible;
  boolean_T cTolComputed;
  boolean_T guard1{ false };

  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  x[3] = 0.0;
  *status = 1;
  r[0] = 0.0;
  r[1] = 0.0;
  r[2] = 0.0;
  r[3] = 0.0;
  rMin = 0.0;
  cTolComputed = false;
  for (i = 0; i < 166; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 166; tmp++) {
    if (iA[tmp]) {
      nA++;
      iC[nA - 1] = tmp + 1;
    }
  }

  guard1 = false;
  if (nA > 0) {
    int32_T exitg3;
    (void)std::memset(&Opt[0], 0, sizeof(real_T) << 3UL);
    Rhs[0] = f[0];
    Rhs[4] = 0.0;
    Rhs[1] = f[1];
    Rhs[5] = 0.0;
    Rhs[2] = f[2];
    Rhs[6] = 0.0;
    Rhs[3] = f[3];
    Rhs[7] = 0.0;
    DualFeasible = false;
    tmp = static_cast<int32_T>(std::round(0.3 * static_cast<real_T>(nA)));
    ColdReset = false;
    do {
      exitg3 = 0;
      if ((!DualFeasible) && (nA > 0) && (*status <= maxiter)) {
        KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            (void)std::memset(&iA[0], 0, 166U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 166U * sizeof(int32_T));
            ColdReset = true;
          }
        } else {
          int32_T U_tmp;
          for (i = 0; i < nA; i++) {
            Rhs[i + 4] = b[iC[i] - 1];
            for (iSave = i + 1; iSave <= nA; iSave++) {
              U[(iSave + (i << 2UL)) - 1] = 0.0;
              for (iC_0 = 0; iC_0 < nA; iC_0++) {
                int32_T U_tmp_0;
                U_tmp = iC_0 << 2UL;
                U_tmp_0 = ((i << 2UL) + iSave) - 1;
                U[U_tmp_0] += RLinv[(U_tmp + iSave) - 1] * RLinv[U_tmp + i];
              }

              U[i + ((iSave - 1) << 2UL)] = U[((i << 2UL) + iSave) - 1];
            }
          }

          for (i = 0; i < 4; i++) {
            Opt[i] = ((b_H[i + 4] * Rhs[1] + b_H[i] * Rhs[0]) + b_H[i + 8] *
                      Rhs[2]) + b_H[i + 12] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i] += D[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }
          }

          U_tmp = 0;
          for (i = 0; i < nA; i++) {
            Opt[i + 4] = ((D[U_tmp + 1] * Rhs[1] + D[U_tmp] * Rhs[0]) + D[U_tmp
                          + 2] * Rhs[2]) + D[U_tmp + 3] * Rhs[3];
            iC_0 = 0;
            for (iSave = 0; iSave < nA; iSave++) {
              Opt[i + 4] += U[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }

            U_tmp += 4;
          }

          Xnorm0 = -1.0E-12;
          i = -1;
          for (iSave = 0; iSave < nA; iSave++) {
            lambda[iC[iSave] - 1] = Opt[iSave + 4];
            cMin = Opt[iSave + 4];
            if ((cMin < Xnorm0) && (iSave + 1 <= nA)) {
              i = iSave;
              Xnorm0 = cMin;
            }
          }

          if (i + 1 <= 0) {
            DualFeasible = true;
            x[0] = Opt[0];
            x[1] = Opt[1];
            x[2] = Opt[2];
            x[3] = Opt[3];
          } else {
            (*status)++;
            if (tmp <= 5) {
              iC_0 = 5;
            } else {
              iC_0 = tmp;
            }

            if (*status > iC_0) {
              nA = 0;
              (void)std::memset(&iA[0], 0, 166U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 166U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 166U * sizeof(real_T));
          for (tmp = 0; tmp <= 2; tmp += 2) {
            tmp_3 = _mm_set1_pd(-1.0);
            (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd
              (_mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3),
                          _mm_set1_pd(f[1])), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
              (&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f[0]))), _mm_mul_pd(_mm_mul_pd
              (_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3), _mm_set1_pd(f[2]))),
              _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3),
                         _mm_set1_pd(f[3]))));
          }
        }

        exitg3 = 1;
      }
    } while (exitg3 == 0);

    if (exitg3 == 1) {
      guard1 = true;
    }
  } else {
    for (tmp = 0; tmp <= 2; tmp += 2) {
      tmp_3 = _mm_set1_pd(-1.0);
      (void)_mm_storeu_pd(&x[tmp], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
        (_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 4]), tmp_3), _mm_set1_pd(f[1])),
        _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp]), tmp_3), _mm_set1_pd(f
        [0]))), _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[tmp + 8]), tmp_3),
                           _mm_set1_pd(f[2]))), _mm_mul_pd(_mm_mul_pd
        (_mm_loadu_pd(&b_Hinv[tmp + 12]), tmp_3), _mm_set1_pd(f[3]))));
    }

    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    Xnorm0 = norm(x);
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (*status <= maxiter)) {
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 166; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 166] * x[1]);
          z[2] = std::abs(b_Ac[i + 332] * x[2]);
          z[3] = std::abs(b_Ac[i + 498] * x[3]);
          t = std::fmax(t, maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 166] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 332] * x
                    [2]) + b_Ac[i + 498] * x[3]) - b[i]) / t;
          if (cVal < cMin) {
            cMin = cVal;
            tmp = i;
          }
        }

        cTol[i] = t;
      }

      cTolComputed = true;
      if (tmp + 1 <= 0) {
        exitg2 = true;
      } else if (*status == maxiter) {
        *status = 0;
        exitg2 = true;
      } else {
        int32_T exitg1;
        do {
          exitg1 = 0;
          if ((tmp + 1 > 0) && (*status <= maxiter)) {
            boolean_T guard2{ false };

            guard2 = false;
            if (nA == 0) {
              for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 498])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 332])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 166])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0]), _mm_set1_pd
                  (b_Ac[tmp])), _mm_set1_pd(0.0))))));
              }

              guard2 = true;
            } else {
              KWIKfactor(b_Ac, iC, nA, b_Linv, D, b_H, degrees, RLinv, &cMin);
              if (cMin <= 0.0) {
                *status = -2;
                exitg1 = 1;
              } else {
                for (iC_0 = 0; iC_0 <= 14; iC_0 += 2) {
                  tmp_3 = _mm_loadu_pd(&b_H[iC_0]);
                  (void)_mm_storeu_pd(&U[iC_0], _mm_mul_pd(tmp_3, _mm_set1_pd
                    (-1.0)));
                }

                for (iC_0 = 0; iC_0 <= 2; iC_0 += 2) {
                  __m128d tmp_0;
                  __m128d tmp_1;
                  __m128d tmp_2;
                  tmp_3 = _mm_loadu_pd(&U[iC_0]);
                  tmp_0 = _mm_loadu_pd(&U[iC_0 + 4]);
                  tmp_1 = _mm_loadu_pd(&U[iC_0 + 8]);
                  tmp_2 = _mm_loadu_pd(&U[iC_0 + 12]);
                  (void)_mm_storeu_pd(&z[iC_0], _mm_add_pd(_mm_mul_pd(tmp_2,
                    _mm_set1_pd(b_Ac[tmp + 498])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 332])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 166])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((D[iSave + 1] * b_Ac[tmp + 166] + D[iSave] * b_Ac[tmp])
                          + D[iSave + 2] * b_Ac[tmp + 332]) + D[iSave + 3] *
                    b_Ac[tmp + 498];
                }

                guard2 = true;
              }
            }

            if (guard2) {
              real_T cVal_tmp;
              real_T cVal_tmp_0;
              boolean_T exitg4;
              i = 0;
              cMin = 0.0;
              DualFeasible = true;
              ColdReset = true;
              if (nA > 0) {
                iSave = 0;
                exitg4 = false;
                while (((exitg4 ? static_cast<uint32_T>(1U) :
                         static_cast<uint32_T>(0U)) == false) && (iSave <= nA -
                        1)) {
                  if (r[iSave] >= 1.0E-12) {
                    ColdReset = false;
                    exitg4 = true;
                  } else {
                    iSave++;
                  }
                }
              }

              if ((nA != 0) && (!ColdReset)) {
                for (iSave = 0; iSave < nA; iSave++) {
                  cVal = r[iSave];
                  if (cVal > 1.0E-12) {
                    cVal = lambda[iC[iSave] - 1] / cVal;
                    if ((i == 0) || (cVal < rMin)) {
                      rMin = cVal;
                      i = iSave + 1;
                    }
                  }
                }

                if (i > 0) {
                  cMin = rMin;
                  DualFeasible = false;
                }
              }

              t = b_Ac[tmp + 166];
              cVal_tmp = b_Ac[tmp + 332];
              cVal_tmp_0 = b_Ac[tmp + 498];
              cVal = ((t * z[1] + z[0] * b_Ac[tmp]) + cVal_tmp * z[2]) +
                cVal_tmp_0 * z[3];
              if (cVal <= 0.0) {
                cVal = 0.0;
                ColdReset = true;
              } else {
                cVal = (b[tmp] - (((t * x[1] + b_Ac[tmp] * x[0]) + cVal_tmp * x
                                   [2]) + cVal_tmp_0 * x[3])) / cVal;
                ColdReset = false;
              }

              if (DualFeasible && ColdReset) {
                *status = -1;
                exitg1 = 1;
              } else {
                if (ColdReset) {
                  t = cMin;
                } else if (DualFeasible) {
                  t = cVal;
                } else if (cMin < cVal) {
                  t = cMin;
                } else {
                  t = cVal;
                }

                for (iSave = 0; iSave < nA; iSave++) {
                  iC_0 = iC[iSave];
                  lambda[iC_0 - 1] -= t * r[iSave];
                  if ((iC_0 <= 166) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                (void)std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  (void)std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == static_cast<int32_T>(degrees)) {
                      *status = -1;
                      exitg1 = 1;
                    } else {
                      nA++;
                      iC[nA - 1] = tmp + 1;
                      i = nA - 1;
                      exitg4 = false;
                      while (((exitg4 ? static_cast<uint32_T>(1U) : static_cast<
                               uint32_T>(0U)) == false) && (i + 1 > 1)) {
                        iC_0 = iC[i - 1];
                        if (iC[i] > iC_0) {
                          exitg4 = true;
                        } else {
                          iSave = iC[i];
                          iC[i] = iC_0;
                          iC[i - 1] = iSave;
                          i--;
                        }
                      }

                      iA[tmp] = true;
                      tmp = -1;
                      (*status)++;
                    }
                  } else {
                    (*status)++;
                  }
                } else {
                  (*status)++;
                }
              }
            }
          } else {
            cMin = norm(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 166; tmp++) {
                cTol[tmp] = std::fmax(std::abs(b[tmp]), 1.0);
              }

              cTolComputed = false;
            }

            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S114>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv_c(const real_T A[12], const real_T B_1[9],
  real_T Y[12])
{
  real_T b_A[9];
  real_T Y_tmp_1;
  real_T Y_tmp_2;
  real_T Y_tmp_3;
  real_T Y_tmp_4;
  real_T a21;
  real_T maxval;
  int32_T Y_tmp;
  int32_T Y_tmp_0;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  int32_T rtemp;
  (void)std::memcpy(&b_A[0], &B_1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = std::abs(B_1[0]);
  a21 = std::abs(B_1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (std::abs(B_1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = B_1[r2] / B_1[r1];
  b_A[r3] /= b_A[r1];
  b_A[r2 + 3] -= b_A[r1 + 3] * b_A[r2];
  b_A[r3 + 3] -= b_A[r1 + 3] * b_A[r3];
  b_A[r2 + 6] -= b_A[r1 + 6] * b_A[r2];
  b_A[r3 + 6] -= b_A[r1 + 6] * b_A[r3];
  if (std::abs(b_A[r3 + 3]) > std::abs(b_A[r2 + 3])) {
    rtemp = r2;
    r2 = r3;
    r3 = rtemp;
  }

  b_A[r3 + 3] /= b_A[r2 + 3];
  b_A[r3 + 6] -= b_A[r3 + 3] * b_A[r2 + 6];
  rtemp = r1 << 2UL;
  Y[rtemp] = A[0] / b_A[r1];
  Y_tmp = r2 << 2UL;
  maxval = b_A[r1 + 3];
  Y[Y_tmp] = A[4] - Y[rtemp] * maxval;
  Y_tmp_0 = r3 << 2UL;
  a21 = b_A[r1 + 6];
  Y[Y_tmp_0] = A[8] - Y[rtemp] * a21;
  Y_tmp_1 = b_A[r2 + 3];
  Y[Y_tmp] /= Y_tmp_1;
  Y_tmp_2 = b_A[r2 + 6];
  Y[Y_tmp_0] -= Y[Y_tmp] * Y_tmp_2;
  Y_tmp_3 = b_A[r3 + 6];
  Y[Y_tmp_0] /= Y_tmp_3;
  Y_tmp_4 = b_A[r3 + 3];
  Y[Y_tmp] -= Y[Y_tmp_0] * Y_tmp_4;
  Y[rtemp] -= Y[Y_tmp_0] * b_A[r3];
  Y[rtemp] -= Y[Y_tmp] * b_A[r2];
  Y[rtemp + 1] = A[1] / b_A[r1];
  Y[Y_tmp + 1] = A[5] - Y[rtemp + 1] * maxval;
  Y[Y_tmp_0 + 1] = A[9] - Y[rtemp + 1] * a21;
  Y[Y_tmp + 1] /= Y_tmp_1;
  Y[Y_tmp_0 + 1] -= Y[Y_tmp + 1] * Y_tmp_2;
  Y[Y_tmp_0 + 1] /= Y_tmp_3;
  Y[Y_tmp + 1] -= Y[Y_tmp_0 + 1] * Y_tmp_4;
  Y[rtemp + 1] -= Y[Y_tmp_0 + 1] * b_A[r3];
  Y[rtemp + 1] -= Y[Y_tmp + 1] * b_A[r2];
  Y[rtemp + 2] = A[2] / b_A[r1];
  Y[Y_tmp + 2] = A[6] - Y[rtemp + 2] * maxval;
  Y[Y_tmp_0 + 2] = A[10] - Y[rtemp + 2] * a21;
  Y[Y_tmp + 2] /= Y_tmp_1;
  Y[Y_tmp_0 + 2] -= Y[Y_tmp + 2] * Y_tmp_2;
  Y[Y_tmp_0 + 2] /= Y_tmp_3;
  Y[Y_tmp + 2] -= Y[Y_tmp_0 + 2] * Y_tmp_4;
  Y[rtemp + 2] -= Y[Y_tmp_0 + 2] * b_A[r3];
  Y[rtemp + 2] -= Y[Y_tmp + 2] * b_A[r2];
  Y[rtemp + 3] = A[3] / b_A[r1];
  Y[Y_tmp + 3] = A[7] - Y[rtemp + 3] * maxval;
  Y[Y_tmp_0 + 3] = A[11] - Y[rtemp + 3] * a21;
  Y[Y_tmp + 3] /= Y_tmp_1;
  Y[Y_tmp_0 + 3] -= Y[Y_tmp + 3] * Y_tmp_2;
  Y[Y_tmp_0 + 3] /= Y_tmp_3;
  Y[Y_tmp + 3] -= Y[Y_tmp_0 + 3] * Y_tmp_4;
  Y[rtemp + 3] -= Y[Y_tmp_0 + 3] * b_A[r3];
  Y[rtemp + 3] -= Y[Y_tmp + 3] * b_A[r2];
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T hi;
  uint32_T lo;

  // Uniform random number generator (random number between 0 and 1)

  // #define IA      16807                      magic multiplier = 7^5
  // #define IM      2147483647                 modulus = 2^31-1
  // #define IQ      127773                     IM div IA
  // #define IR      2836                       IM modulo IA
  // #define S       4.656612875245797e-10      reciprocal of 2^31-1
  // test = IA * (seed % IQ) - IR * (seed/IQ)
  // seed = test < 0 ? (test + IM) : test
  // return (seed*S)

  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return static_cast<real_T>(*u) * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T si;
  real_T sr;
  real_T y;

  // Normal (Gaussian) random number generator
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = std::sqrt(-2.0 * std::log(si) / si) * sr;
  return y;
}

// Model step function
void SupervisoryController::step()
{
  static const real_T b_a_0[1442]{ -0.0, -0.99801496201563555,
    -0.00012982817451900192, -0.025, -0.0, -0.0, -0.99603367687094413,
    -0.00025924037690453866, -0.049950374050390892, -3.2457043629750481E-6, -0.0,
    -0.99405613771726609, -0.00038823760144369459, -0.074851215972164487,
    -9.726713785588515E-6, -0.0, -0.99208233771810028, -0.00051682084032238652,
    -0.099702619415096144, -1.9432653821680881E-5, -0.0, -0.99011227004908309,
    -0.00064499108362950342, -0.12450467785804867, -3.2353174829740542E-5, -0.0,
    -0.9881459278979674, -0.00077274931936103868, -0.14925748460927574,
    -4.8477951920478122E-5, -0.0, -0.98618330446460145, -0.00090009653342421417,
    -0.17396113280672493, -6.7796684904504082E-5, -0.0, -0.984224392960908,
    -0.0010270337096415967, -0.19861571541833997, -9.0299098240109431E-5, -0.0,
    -0.98226918661086315, -0.001153561829755207, -0.22322132524236266,
    -0.00011597494098114935, -0.0, -0.98031767865047559, -0.0012796818734306208,
    -0.24777805490763424, -0.00014481398672502952, -0.0, -0.97836986232776568,
    -0.0014053948182610615, -0.27228599687389615, -0.00017680603356079504, -0.0,
    -0.97642573090274454, -0.0015307016397714866, -0.2967452434320903,
    -0.00021194090401732159, -0.0, -0.97448527764739334, -0.0016556033114226655,
    -0.32115588670465894, -0.00025020844501160875, -0.0, -0.97254849584564251,
    -0.0017801008046152487, -0.34551801864584375, -0.00029159852779717539, -0.0,
    -0.97061537879335091, -0.0019041950886938322, -0.3698317310419848,
    -0.00033610104791255663, -0.0, -0.968685919798285, -0.0020278871309510108,
    -0.39409711551181864, -0.00038370592512990245, -0.0, -0.96676011218009861,
    -0.0021511778966314264, -0.41831426350677575, -0.00043440310340367769, -0.0,
    -0.96483794927031175, -0.0022740683489358075, -0.44248326631127827,
    -0.00048818255081946337, -0.0, -0.9629194244122905, -0.0023965594490250011,
    -0.46660421504303612, -0.00054503425954285857, -0.0, -0.96100453096122618,
    -0.0025186521560239968, -0.49067720065334341, -0.00060494824576848359, 0.0,
    0.99801496201563555, 0.00012982817451900192, 0.025, 0.0, 0.0,
    0.99603367687094413, 0.00025924037690453866, 0.049950374050390892,
    3.2457043629750481E-6, 0.0, 0.99405613771726609, 0.00038823760144369459,
    0.074851215972164487, 9.726713785588515E-6, 0.0, 0.99208233771810028,
    0.00051682084032238652, 0.099702619415096144, 1.9432653821680881E-5, 0.0,
    0.99011227004908309, 0.00064499108362950342, 0.12450467785804867,
    3.2353174829740542E-5, 0.0, 0.9881459278979674, 0.00077274931936103868,
    0.14925748460927574, 4.8477951920478122E-5, 0.0, 0.98618330446460145,
    0.00090009653342421417, 0.17396113280672493, 6.7796684904504082E-5, 0.0,
    0.984224392960908, 0.0010270337096415967, 0.19861571541833997,
    9.0299098240109431E-5, 0.0, 0.98226918661086315, 0.001153561829755207,
    0.22322132524236266, 0.00011597494098114935, 0.0, 0.98031767865047559,
    0.0012796818734306208, 0.24777805490763424, 0.00014481398672502952, 0.0,
    0.97836986232776568, 0.0014053948182610615, 0.27228599687389615,
    0.00017680603356079504, 0.0, 0.97642573090274454, 0.0015307016397714866,
    0.2967452434320903, 0.00021194090401732159, 0.0, 0.97448527764739334,
    0.0016556033114226655, 0.32115588670465894, 0.00025020844501160875, 0.0,
    0.97254849584564251, 0.0017801008046152487, 0.34551801864584375,
    0.00029159852779717539, 0.0, 0.97061537879335091, 0.0019041950886938322,
    0.3698317310419848, 0.00033610104791255663, 0.0, 0.968685919798285,
    0.0020278871309510108, 0.39409711551181864, 0.00038370592512990245, 0.0,
    0.96676011218009861, 0.0021511778966314264, 0.41831426350677575,
    0.00043440310340367769, 0.0, 0.96483794927031175, 0.0022740683489358075,
    0.44248326631127827, 0.00048818255081946337, 0.0, 0.9629194244122905,
    0.0023965594490250011, 0.46660421504303612, 0.00054503425954285857, 0.0,
    0.96100453096122618, 0.0025186521560239968, 0.49067720065334341,
    0.00060494824576848359, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    0.001444494825713849, -0.99878101747798131, -0.0, -0.025, -0.0,
    0.0028843614603847466, -0.99756333333822533, 3.6112370642846224E-5,
    -0.049969525436949533, -0.0, 0.0043196109666541723, -0.99634694659889633,
    0.00010822140715246489, -0.074908608770405172, -0.0, 0.0057502543837855847,
    -0.99513185627791911, 0.00021621168131881919, -0.099817282435377575, -0.0,
    0.0071763027277104837, -0.99391806139298267, 0.00035996804091345877,
    -0.12469557884232554, -0.0, 0.0085977669910743838, -0.99270556096154294,
    0.0005393756091062209, -0.1495435303771501, -0.0, 0.0100146581432827,
    -0.99149435400082664, 0.00075431978388308045, -0.17436116940118868, -0.0,
    0.011426987130546547, -0.99028443952783429, 0.0010046862374651479,
    -0.19914852825120935, -0.0, 0.012834764875928464, -0.9890758165593434,
    0.0012903609157288115, -0.22390563923940521, -0.0, 0.014238002279388031,
    -0.987868484111912, 0.0016112300376270232, -0.24863253465338878, -0.0,
    0.015636710217827417, -0.98666244120188173, 0.0019671800946117241,
    -0.27332924675618658, -0.0, 0.017030899545136844, -0.98545768684538126,
    0.0023580978500574094, -0.29799580778623364, -0.0, 0.018420581092239943,
    -0.98425422005832941, 0.0027838703386858306, -0.32263224995736817, -0.0,
    0.019805765667139059, -0.98305203985643841, 0.0032443848659918294,
    -0.3472386054588264, -0.0, 0.021186464054960434, -0.98185114525521711,
    0.0037395290076703059, -0.3718149064552374, -0.0, 0.022562687017999343,
    -0.98065153526997417, 0.0042691906090443167, -0.39636118508661783, -0.0,
    0.023934445295765105, -0.97945320891582144, 0.0048332577844943008,
    -0.42087747346836718, -0.0, 0.025301749605026048, -0.97825616520767711,
    0.0054316189168884289, -0.44536380369126272, -0.0, 0.026664610639854355,
    -0.97706040316026854, 0.0060641626570140808, -0.46982020782145467, -0.0,
    0.028023039071670849, -0.975865921788136, 0.00673077792301044,
    -0.49424671790046143, 0.0, -0.001444494825713849, 0.99878101747798131, 0.0,
    0.025, 0.0, -0.0028843614603847466, 0.99756333333822533,
    -3.6112370642846224E-5, 0.049969525436949533, 0.0, -0.0043196109666541723,
    0.99634694659889633, -0.00010822140715246489, 0.074908608770405172, 0.0,
    -0.0057502543837855847, 0.99513185627791911, -0.00021621168131881919,
    0.099817282435377575, 0.0, -0.0071763027277104837, 0.99391806139298267,
    -0.00035996804091345877, 0.12469557884232554, 0.0, -0.0085977669910743838,
    0.99270556096154294, -0.0005393756091062209, 0.1495435303771501, 0.0,
    -0.0100146581432827, 0.99149435400082664, -0.00075431978388308045,
    0.17436116940118868, 0.0, -0.011426987130546547, 0.99028443952783429,
    -0.0010046862374651479, 0.19914852825120935, 0.0, -0.012834764875928464,
    0.9890758165593434, -0.0012903609157288115, 0.22390563923940521, 0.0,
    -0.014238002279388031, 0.987868484111912, -0.0016112300376270232,
    0.24863253465338878, 0.0, -0.015636710217827417, 0.98666244120188173,
    -0.0019671800946117241, 0.27332924675618658, 0.0, -0.017030899545136844,
    0.98545768684538126, -0.0023580978500574094, 0.29799580778623364, 0.0,
    -0.018420581092239943, 0.98425422005832941, -0.0027838703386858306,
    0.32263224995736817, 0.0, -0.019805765667139059, 0.98305203985643841,
    -0.0032443848659918294, 0.3472386054588264, 0.0, -0.021186464054960434,
    0.98185114525521711, -0.0037395290076703059, 0.3718149064552374, 0.0,
    -0.022562687017999343, 0.98065153526997417, -0.0042691906090443167,
    0.39636118508661783, 0.0, -0.023934445295765105, 0.97945320891582144,
    -0.0048332577844943008, 0.42087747346836718, 0.0, -0.025301749605026048,
    0.97825616520767711, -0.0054316189168884289, 0.44536380369126272, 0.0,
    -0.026664610639854355, 0.97706040316026854, -0.0060641626570140808,
    0.46982020782145467, 0.0, -0.028023039071670849, 0.975865921788136,
    -0.00673077792301044, 0.49424671790046143, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_a_1[882]{ -0.99747143658641013, -0.0011849245467832059,
    -0.0, -0.9949489107331112, -0.0023674262193844361, -0.0,
    -0.99243240790182863, -0.0035475109996374703, -0.0, -0.9899219135892513,
    -0.0047251848550434227, -0.0, -0.98741741332694732, -0.0059004537388052374,
    -0.0, -0.98491889268128, -0.0070733235898620978, -0.0, -0.98242633725332451,
    -0.0082438003329237571, -0.0, -0.97993973267878431, -0.009411889878504786,
    -0.0, -0.97745906462790788, -0.010577598122958734, -0.0,
    -0.97498431880540515, -0.011740930948512213, -0.0, -0.97251548095036522,
    -0.012901894223298897, -0.0, -0.97005253683617321, -0.014060493801393438,
    -0.0, -0.9675954722704283, -0.015216735522845306, -0.0, -0.96514427309486084,
    -0.01637062521371254, -0.0, -0.96269892518525058, -0.017522168686095428,
    -0.0, -0.96025941445134477, -0.018671371738170094, -0.0,
    -0.95782572683677636, -0.019818240154222011, -0.0, -0.95539784831898256,
    -0.020962779704679434, -0.0, -0.95297576490912339, -0.022104996146146753,
    -0.0, -0.95055946265200075, -0.023244895221437765, -0.0, 0.99747143658641013,
    0.0011849245467832059, 0.0, 0.9949489107331112, 0.0023674262193844361, 0.0,
    0.99243240790182863, 0.0035475109996374703, 0.0, 0.9899219135892513,
    0.0047251848550434227, 0.0, 0.98741741332694732, 0.0059004537388052374, 0.0,
    0.98491889268128, 0.0070733235898620978, 0.0, 0.98242633725332451,
    0.0082438003329237571, 0.0, 0.97993973267878431, 0.009411889878504786, 0.0,
    0.97745906462790788, 0.010577598122958734, 0.0, 0.97498431880540515,
    0.011740930948512213, 0.0, 0.97251548095036522, 0.012901894223298897, 0.0,
    0.97005253683617321, 0.014060493801393438, 0.0, 0.9675954722704283,
    0.015216735522845306, 0.0, 0.96514427309486084, 0.01637062521371254, 0.0,
    0.96269892518525058, 0.017522168686095428, 0.0, 0.96025941445134477,
    0.018671371738170094, 0.0, 0.95782572683677636, 0.019818240154222011, 0.0,
    0.95539784831898256, 0.020962779704679434, 0.0, 0.95297576490912339,
    0.022104996146146753, 0.0, 0.95055946265200075, 0.023244895221437765, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00030050237926615431, -1.0004838136541614,
    -0.0, 0.00060039030636459647, -1.0009675053113292, -0.0,
    0.00089966529831621029, -1.001451075640559, -0.0, 0.0011983288685070489,
    -1.0019345253094327, -0.0, 0.0014963825266970827, -1.0024178549840617, -0.0,
    0.0017938277790289267, -1.0029010653290908, -0.0, 0.0020906661280365455,
    -1.003384157007702, -0.0, 0.0023868990726539398, -1.0038671306816174, -0.0,
    0.0026825281082238088, -1.004349987011103, -0.0, 0.0029775547265061954,
    -1.004832726654973, -0.0, 0.003271980415687107, -1.005315350270592, -0.0,
    0.0035658066603871182, -1.0057978585138798, -0.0, 0.003859034941669952,
    -1.0062802520393141, -0.0, 0.00415166673705104, -1.0067625314999344, -0.0,
    0.0044437035205060612, -1.0072446975473455, -0.0, 0.0047351467624794641,
    -1.0077267508317205, -0.0, 0.0050259979298929629, -1.008208692001805, -0.0,
    0.0053162584861540171, -1.0086905217049205, -0.0, 0.0056059298911642881,
    -1.0091722405869672, -0.0, 0.0058950136013280795, -1.0096538492924283, -0.0,
    -0.00030050237926615431, 1.0004838136541614, 0.0, -0.00060039030636459647,
    1.0009675053113292, 0.0, -0.00089966529831621029, 1.001451075640559, 0.0,
    -0.0011983288685070489, 1.0019345253094327, 0.0, -0.0014963825266970827,
    1.0024178549840617, 0.0, -0.0017938277790289267, 1.0029010653290908, 0.0,
    -0.0020906661280365455, 1.003384157007702, 0.0, -0.0023868990726539398,
    1.0038671306816174, 0.0, -0.0026825281082238088, 1.004349987011103, 0.0,
    -0.0029775547265061954, 1.004832726654973, 0.0, -0.003271980415687107,
    1.005315350270592, 0.0, -0.0035658066603871182, 1.0057978585138798, 0.0,
    -0.003859034941669952, 1.0062802520393141, 0.0, -0.00415166673705104,
    1.0067625314999344, 0.0, -0.0044437035205060612, 1.0072446975473455, 0.0,
    -0.0047351467624794641, 1.0077267508317205, 0.0, -0.0050259979298929629,
    1.008208692001805, 0.0, -0.0053162584861540171, 1.0086905217049205, 0.0,
    -0.0056059298911642881, 1.0091722405869672, 0.0, -0.0058950136013280795,
    1.0096538492924283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_a[830]{ -0.99636289448383575, -0.0, -0.0, -0.025,
    -0.99273901750420723, -0.0, -0.0, -0.049909072362095894,
    -0.98912832094753123, -0.0, -0.0, -0.074727547799701075,
    -0.98553075687521863, -0.0, -0.0, -0.099455755823389363, -0.9819462775230382,
    -0.0, -0.0, -0.12409402474526984, -0.97837483530048219, -0.0, -0.0,
    -0.1486426816833458, -0.97481638279013449, -0.0, -0.0, -0.17310205256585784,
    -0.97127087274704116, -0.0, -0.0, -0.1974724621356112, -0.96773825809808323,
    -0.0, -0.0, -0.22175423395428723, -0.96421849194135145, -0.0, -0.0,
    -0.24594769040673931, -0.960711527545524, -0.0, -0.0, -0.27005315270527308,
    -0.95721731834924562, -0.0, -0.0, -0.29407094089391123, -0.95373581796050966,
    -0.0, -0.0, -0.31800137385264238, -0.950266980156042, -0.0, -0.0,
    -0.34184476930165514, -0.94681075888068778, -0.0, -0.0, -0.3656014438055562,
    -0.94336710824679915, -0.0, -0.0, -0.38927171277757344, -0.93993598253362676,
    -0.0, -0.0, -0.41285589048374344, -0.9365173361867124, -0.0, -0.0,
    -0.43635429004708415, -0.93311112381728423, -0.0, -0.0, -0.45976722345175197,
    -0.92971730020165411, -0.0, -0.0, -0.4830950015471841, 0.99636289448383575,
    0.0, 0.0, 0.025, 0.99273901750420723, 0.0, 0.0, 0.049909072362095894,
    0.98912832094753123, 0.0, 0.0, 0.074727547799701075, 0.98553075687521863,
    0.0, 0.0, 0.099455755823389363, 0.9819462775230382, 0.0, 0.0,
    0.12409402474526984, 0.97837483530048219, 0.0, 0.0, 0.1486426816833458,
    0.97481638279013449, 0.0, 0.0, 0.17310205256585784, 0.97127087274704116, 0.0,
    0.0, 0.1974724621356112, 0.96773825809808323, 0.0, 0.0, 0.22175423395428723,
    0.96421849194135145, 0.0, 0.0, 0.24594769040673931, 0.960711527545524, 0.0,
    0.0, 0.27005315270527308, 0.95721731834924562, 0.0, 0.0, 0.29407094089391123,
    0.95373581796050966, 0.0, 0.0, 0.31800137385264238, 0.950266980156042, 0.0,
    0.0, 0.34184476930165514, 0.94681075888068778, 0.0, 0.0, 0.3656014438055562,
    0.94336710824679915, 0.0, 0.0, 0.38927171277757344, 0.93993598253362676, 0.0,
    0.0, 0.41285589048374344, 0.9365173361867124, 0.0, 0.0, 0.43635429004708415,
    0.93311112381728423, 0.0, 0.0, 0.45976722345175197, 0.92971730020165411, 0.0,
    0.0, 0.4830950015471841, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Ac_0[824]{ -0.0, 0.1888968412469153,
    0.060492001339316988, -0.0, -0.0, -0.0, 0.3773313347058973,
    0.1209347881183509, 0.0047224210311728829, 0.0015123000334829248, -0.0,
    0.56530446924504474, 0.18132836030302127, 0.014155704398820317,
    0.0045356697364416975, -0.0, 0.75281723176956472, 0.24167271798767209,
    0.028288316129946436, 0.0090688787440172287, -0.0, 0.939870607225484,
    0.30196786139466048, 0.047108746924185554, 0.015110696693709031, -0.0,
    1.1264655786033526, 0.36221379087394617, 0.07060551210482266,
    0.022659893228575541, -0.0, 1.3126031269419416, 0.42241050690268223,
    0.098767151569906481, 0.031715238000424195, -0.0, 1.4982842313319336,
    0.48255801008480659, 0.13158222974345501, 0.042275500672991254, -0.0,
    1.6835098689196064, 0.54265630115063446, 0.16903933552675335,
    0.054339450925111416, -0.0, 1.8682810149105107, 0.6027053809564521,
    0.2111270822497435, 0.067905858453877282, -0.0, 2.05259864257314,
    0.66270525048411089, 0.25783410762250625, 0.082973492977788582, -0.0,
    2.2364637232425952, 0.72265591084062319, 0.30914907368683475,
    0.099541124239891351, -0.0, 2.41987722632424, 0.78255736325775849,
    0.36506066676789961, 0.11760752201090693, -0.0, 2.6028401192973529,
    0.842409609091641, 0.42555759742600563, 0.1371714560923509, -0.0,
    2.7853533677187707, 0.90221264982234783, 0.49062860040843947,
    0.15823169631964193, -0.0, 2.9674179352265249, 0.96196648705350851,
    0.56026243460140879, 0.18078701256520063, -0.0, 3.1490347835434722,
    1.0216711225119051, 0.63444788298207189, 0.20483617474153834, -0.0,
    3.3302048724809192, 1.081326558047073, 0.71317375257065874,
    0.23037795280433598, -0.0, 3.5109291599422385, 1.1409327956309046,
    0.79642887438268173, 0.25741111675551281, -0.0, 3.69120860192648,
    1.2004898373572503, 0.88420210338123773, 0.28593443664628543, 0.0,
    -0.1888968412469153, -0.060492001339316988, 0.0, 0.0, 0.0,
    -0.3773313347058973, -0.1209347881183509, -0.0047224210311728829,
    -0.0015123000334829248, 0.0, -0.56530446924504474, -0.18132836030302127,
    -0.014155704398820317, -0.0045356697364416975, 0.0, -0.75281723176956472,
    -0.24167271798767209, -0.028288316129946436, -0.0090688787440172287, 0.0,
    -0.939870607225484, -0.30196786139466048, -0.047108746924185554,
    -0.015110696693709031, 0.0, -1.1264655786033526, -0.36221379087394617,
    -0.07060551210482266, -0.022659893228575541, 0.0, -1.3126031269419416,
    -0.42241050690268223, -0.098767151569906481, -0.031715238000424195, 0.0,
    -1.4982842313319336, -0.48255801008480659, -0.13158222974345501,
    -0.042275500672991254, 0.0, -1.6835098689196064, -0.54265630115063446,
    -0.16903933552675335, -0.054339450925111416, 0.0, -1.8682810149105107,
    -0.6027053809564521, -0.2111270822497435, -0.067905858453877282, 0.0,
    -2.05259864257314, -0.66270525048411089, -0.25783410762250625,
    -0.082973492977788582, 0.0, -2.2364637232425952, -0.72265591084062319,
    -0.30914907368683475, -0.099541124239891351, 0.0, -2.41987722632424,
    -0.78255736325775849, -0.36506066676789961, -0.11760752201090693, 0.0,
    -2.6028401192973529, -0.842409609091641, -0.42555759742600563,
    -0.1371714560923509, 0.0, -2.7853533677187707, -0.90221264982234783,
    -0.49062860040843947, -0.15823169631964193, 0.0, -2.9674179352265249,
    -0.96196648705350851, -0.56026243460140879, -0.18078701256520063, 0.0,
    -3.1490347835434722, -1.0216711225119051, -0.63444788298207189,
    -0.20483617474153834, 0.0, -3.3302048724809192, -1.081326558047073,
    -0.71317375257065874, -0.23037795280433598, 0.0, -3.5109291599422385,
    -1.1409327956309046, -0.79642887438268173, -0.25741111675551281, 0.0,
    -3.69120860192648, -1.2004898373572503, -0.88420210338123773,
    -0.28593443664628543, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, -0.0,
    -0.1599184353223182, 0.089801925435736424, -0.0, -0.0, -0.0,
    -0.31964914489275253, 0.17947362197539324, -0.0039979608830579549,
    0.0022450481358934106, -0.0, -0.47919231323908462, 0.26901527273781028,
    -0.011989189505376769, 0.0067318886852782414, -0.0, -0.63854812478731549,
    0.35842706059465185, -0.023968997336353881, 0.013457270503723497, -0.0,
    -0.797716763861511, 0.44770916817072121, -0.039932700456036765,
    0.022417947018589793, -0.0, -0.956698414683647, 0.53686177784427469,
    -0.059875619552574547, 0.033610676222857822, -0.0, -1.1154932613734565,
    0.62588507174733543, -0.083793079919665725, 0.047032220668964689, -0.0,
    -1.2741014879482764, 0.71477923176600711, -0.11168041145400215,
    0.062679347462648069, -0.0, -1.4325232783228956, 0.80354443954078658,
    -0.14353294865270905, 0.080548828256798238, -0.0, -1.5907588163094044,
    0.89218087646687683, -0.17934603061078144, 0.1006374392453179, -0.0,
    -1.748808285617043, 0.98068872369449933, -0.21911500101851655,
    0.12294196115698981, -0.0, -1.9066718698520528, 1.069068162129206,
    -0.26283520815894262, 0.14745917924935228, -0.0, -2.0643497525175265,
    1.1573193724321906, -0.31050200490524393, 0.17418588330258244, -0.0,
    -2.221842117013261, 1.2454425350206004, -0.36211074871818211,
    0.20311886761338721, -0.0, -2.3791491466356085, 1.3334378300678469,
    -0.41765680164351365, 0.23425493098890221, -0.0, -2.536271024577331,
    1.421305437503916, -0.47713553030940387, 0.26759087674059839, -0.0,
    -2.6932079339274537, 1.5090455370156788, -0.54054230592383712,
    0.30312351267819632, -0.0, -2.8499600576711197, 1.596658308047201,
    -0.60787250427202344, 0.34084965110358828, -0.0, -3.0065275786894463,
    1.6841439298000525, -0.67912150571380148, 0.3807661088047683, -0.0,
    -3.16291067975938, 1.7715025812336165, -0.75428469518103769,
    0.42286970704976962, 0.0, 0.1599184353223182, -0.089801925435736424, 0.0,
    0.0, 0.0, 0.31964914489275253, -0.17947362197539324, 0.0039979608830579549,
    -0.0022450481358934106, 0.0, 0.47919231323908462, -0.26901527273781028,
    0.011989189505376769, -0.0067318886852782414, 0.0, 0.63854812478731549,
    -0.35842706059465185, 0.023968997336353881, -0.013457270503723497, 0.0,
    0.797716763861511, -0.44770916817072121, 0.039932700456036765,
    -0.022417947018589793, 0.0, 0.956698414683647, -0.53686177784427469,
    0.059875619552574547, -0.033610676222857822, 0.0, 1.1154932613734565,
    -0.62588507174733543, 0.083793079919665725, -0.047032220668964689, 0.0,
    1.2741014879482764, -0.71477923176600711, 0.11168041145400215,
    -0.062679347462648069, 0.0, 1.4325232783228956, -0.80354443954078658,
    0.14353294865270905, -0.080548828256798238, 0.0, 1.5907588163094044,
    -0.89218087646687683, 0.17934603061078144, -0.1006374392453179, 0.0,
    1.748808285617043, -0.98068872369449933, 0.21911500101851655,
    -0.12294196115698981, 0.0, 1.9066718698520528, -1.069068162129206,
    0.26283520815894262, -0.14745917924935228, 0.0, 2.0643497525175265,
    -1.1573193724321906, 0.31050200490524393, -0.17418588330258244, 0.0,
    2.221842117013261, -1.2454425350206004, 0.36211074871818211,
    -0.20311886761338721, 0.0, 2.3791491466356085, -1.3334378300678469,
    0.41765680164351365, -0.23425493098890221, 0.0, 2.536271024577331,
    -1.421305437503916, 0.47713553030940387, -0.26759087674059839, 0.0,
    2.6932079339274537, -1.5090455370156788, 0.54054230592383712,
    -0.30312351267819632, 0.0, 2.8499600576711197, -1.596658308047201,
    0.60787250427202344, -0.34084965110358828, 0.0, 3.0065275786894463,
    -1.6841439298000525, 0.67912150571380148, -0.3807661088047683, 0.0,
    3.16291067975938, -1.7715025812336165, 0.75428469518103769,
    -0.42286970704976962, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, -0.0,
    0.10536597514810161, -0.1119360731420971, -0.0, -0.0, -0.0,
    0.21068448591173913, -0.22372201869524067, 0.0026341493787025404,
    -0.0027984018285524275, -0.0, 0.31595540965099317, -0.3353580258245521,
    0.0079012615264960184, -0.008391452295933445, -0.0, 0.4211786242426373,
    -0.446844283480486, 0.015800146767770848, -0.016775402941547252, -0.0,
    0.52635400807880228, -0.55818098039902453, 0.026329612373836783,
    -0.027946510028559402, -0.0, 0.63148144006564322, -0.6693683051018724,
    0.039488462575806842, -0.041901034538535019, -0.0, 0.73656079962200982,
    -0.78040644589665065, 0.055275498577447929, -0.058635242166081826, -0.0,
    0.84159196667811931, -0.89129559087709187, 0.073689518567998186,
    -0.078145403313498091, -0.0, 0.94657482167423235, -1.0020359279232338,
    0.094729317734951177, -0.10042779308542539, -0.0, 1.051509245559332,
    -1.1126276447016141, 0.11839368827680699, -0.12547869128350625, -0.0,
    1.1563951197898057, -1.2230709286654644, 0.1446814194157903,
    -0.15329438240104659, -0.0, 1.2612323263281291, -1.3333659670549041,
    0.17359129741053544, -0.18387115561768319, -0.0, 1.3660207476415547,
    -1.4435129468971351, 0.20512210556873869, -0.21720530479405578, -0.0,
    1.4707602667008017, -1.5535120550066355, 0.23927262425977758,
    -0.25329312846648416, -0.0, 1.5754507669787494, -1.6633634779853532,
    0.27604163092729761, -0.29213092984165007, -0.0, 1.6800921324491334,
    -1.7730674022229005, 0.31542790010176636, -0.3337150167912839, -0.0,
    1.7846842475852445, -1.8826240138967472, 0.3574302034129947,
    -0.37804170184685643, -0.0, 1.8892269973586309, -1.9920334989724147,
    0.4020473096026258, -0.42510730219427512, -0.0, 1.9937202672378032,
    -2.1012960432036696, 0.44927798453659162, -0.47490813966858547, -0.0,
    2.098163943186941, -2.2104118321327175, 0.49912099121753672,
    -0.52744054074867719, 0.0, -0.10536597514810161, 0.1119360731420971, 0.0,
    0.0, 0.0, -0.21068448591173913, 0.22372201869524067, -0.0026341493787025404,
    0.0027984018285524275, 0.0, -0.31595540965099317, 0.3353580258245521,
    -0.0079012615264960184, 0.008391452295933445, 0.0, -0.4211786242426373,
    0.446844283480486, -0.015800146767770848, 0.016775402941547252, 0.0,
    -0.52635400807880228, 0.55818098039902453, -0.026329612373836783,
    0.027946510028559402, 0.0, -0.63148144006564322, 0.6693683051018724,
    -0.039488462575806842, 0.041901034538535019, 0.0, -0.73656079962200982,
    0.78040644589665065, -0.055275498577447929, 0.058635242166081826, 0.0,
    -0.84159196667811931, 0.89129559087709187, -0.073689518567998186,
    0.078145403313498091, 0.0, -0.94657482167423235, 1.0020359279232338,
    -0.094729317734951177, 0.10042779308542539, 0.0, -1.051509245559332,
    1.1126276447016141, -0.11839368827680699, 0.12547869128350625, 0.0,
    -1.1563951197898057, 1.2230709286654644, -0.1446814194157903,
    0.15329438240104659, 0.0, -1.2612323263281291, 1.3333659670549041,
    -0.17359129741053544, 0.18387115561768319, 0.0, -1.3660207476415547,
    1.4435129468971351, -0.20512210556873869, 0.21720530479405578, 0.0,
    -1.4707602667008017, 1.5535120550066355, -0.23927262425977758,
    0.25329312846648416, 0.0, -1.5754507669787494, 1.6633634779853532,
    -0.27604163092729761, 0.29213092984165007, 0.0, -1.6800921324491334,
    1.7730674022229005, -0.31542790010176636, 0.3337150167912839, 0.0,
    -1.7846842475852445, 1.8826240138967472, -0.3574302034129947,
    0.37804170184685643, 0.0, -1.8892269973586309, 1.9920334989724147,
    -0.4020473096026258, 0.42510730219427512, 0.0, -1.9937202672378032,
    2.1012960432036696, -0.44927798453659162, 0.47490813966858547, 0.0,
    -2.098163943186941, 2.2104118321327175, -0.49912099121753672,
    0.52744054074867719, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Ac[664]{ -0.20979775256106059, -0.0, -0.0, -0.0,
    -0.41883244855900248, -0.0, -0.0, -0.0052449438140265148,
    -0.62710686331106058, -0.0, -0.0, -0.015715755028001578, -0.834623762040348,
    -0.0, -0.0, -0.0313934266107781, -1.04138589991257, -0.0, -0.0,
    -0.052259020661786795, -1.247396022072603, -0.0, -0.0, -0.078293668159601051,
    -1.452656863680942, -0.0, -0.0, -0.10947856871141612, -1.6571711499500148,
    -0.0, -0.0, -0.14579499030343968, -1.8609415961803639, -0.0, -0.0,
    -0.18722426905219003, -2.0639709077966972, -0.0, -0.0, -0.23374780895669911,
    -2.2662617803838079, -0.0, -0.0, -0.28534708165161654, -2.4678168997223624,
    -0.0, -0.0, -0.34200362616121172, -2.6686389418245593, -0.0, -0.0,
    -0.40369904865427075, -2.8687305729696591, -0.0, -0.0, -0.47041502219988474,
    -3.0680944497393825, -0.0, -0.0, -0.54213328652412618, -3.266733219053183,
    -0.0, -0.0, -0.6188356477676108, -3.4646495182033883, -0.0, -0.0,
    -0.7005039782439404, -3.6618459748902157, -0.0, -0.0, -0.78712021619902517,
    -3.8583252072566592, -0.0, -0.0, -0.87866636557128064, -4.0540898239232508,
    -0.0, -0.0, -0.97512449575269711, 0.20979775256106059, 0.0, 0.0, 0.0,
    0.41883244855900248, 0.0, 0.0, 0.0052449438140265148, 0.62710686331106058,
    0.0, 0.0, 0.015715755028001578, 0.834623762040348, 0.0, 0.0,
    0.0313934266107781, 1.04138589991257, 0.0, 0.0, 0.052259020661786795,
    1.247396022072603, 0.0, 0.0, 0.078293668159601051, 1.452656863680942, 0.0,
    0.0, 0.10947856871141612, 1.6571711499500148, 0.0, 0.0, 0.14579499030343968,
    1.8609415961803639, 0.0, 0.0, 0.18722426905219003, 2.0639709077966972, 0.0,
    0.0, 0.23374780895669911, 2.2662617803838079, 0.0, 0.0, 0.28534708165161654,
    2.4678168997223624, 0.0, 0.0, 0.34200362616121172, 2.6686389418245593, 0.0,
    0.0, 0.40369904865427075, 2.8687305729696591, 0.0, 0.0, 0.47041502219988474,
    3.0680944497393825, 0.0, 0.0, 0.54213328652412618, 3.266733219053183, 0.0,
    0.0, 0.6188356477676108, 3.4646495182033883, 0.0, 0.0, 0.7005039782439404,
    3.6618459748902157, 0.0, 0.0, 0.78712021619902517, 3.8583252072566592, 0.0,
    0.0, 0.87866636557128064, 4.0540898239232508, 0.0, 0.0, 0.97512449575269711,
    -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.22599738357714264, -0.0, -0.0, -0.0,
    0.45117279082383815, -0.0, -0.0, 0.0056499345894285664, 0.67552921135473221,
    -0.0, -0.0, 0.016929254360024522, 0.89906962391092649, -0.0, -0.0,
    0.033817484643892823, 1.1217969963995269, -0.0, -0.0, 0.056294225241665982,
    1.3437142859330482, -0.0, -0.0, 0.084339150151654163, 1.5648244388686749,
    -0.0, -0.0, 0.11793200729998038, 1.7851303908473797, -0.0, -0.0,
    0.15705261827169725, 2.0046350668328987, -0.0, -0.0, 0.20168087804288176,
    2.2233413811505671, -0.0, -0.0, 0.25179675471370422, 2.4412522375260108,
    -0.0, -0.0, 0.30738028924246841, 2.6583705291236992, -0.0, -0.0,
    0.36841159518061867, 2.8746991385853575, -0.0, -0.0, 0.43487085840871115,
    3.0902409380682387, -0.0, -0.0, 0.50673833687334513, 3.3049987892832569,
    -0.0, -0.0, 0.58399436032505114, 3.518975543532981, -0.0, -0.0,
    0.66661933005713259, 3.7321740417494929, -0.0, -0.0, 0.7545937186454571,
    3.9445971145321033, -0.0, -0.0, 0.84789806968919446, 4.1562475821849354,
    -0.0, -0.0, 0.94651299755249707, 4.3671282547543688, -0.0, -0.0,
    1.0504191871071205, -0.22599738357714264, 0.0, 0.0, 0.0,
    -0.45117279082383815, 0.0, 0.0, -0.0056499345894285664, -0.67552921135473221,
    0.0, 0.0, -0.016929254360024522, -0.89906962391092649, 0.0, 0.0,
    -0.033817484643892823, -1.1217969963995269, 0.0, 0.0, -0.056294225241665982,
    -1.3437142859330482, 0.0, 0.0, -0.084339150151654163, -1.5648244388686749,
    0.0, 0.0, -0.11793200729998038, -1.7851303908473797, 0.0, 0.0,
    -0.15705261827169725, -2.0046350668328987, 0.0, 0.0, -0.20168087804288176,
    -2.2233413811505671, 0.0, 0.0, -0.25179675471370422, -2.4412522375260108,
    0.0, 0.0, -0.30738028924246841, -2.6583705291236992, 0.0, 0.0,
    -0.36841159518061867, -2.8746991385853575, 0.0, 0.0, -0.43487085840871115,
    -3.0902409380682387, 0.0, 0.0, -0.50673833687334513, -3.3049987892832569,
    0.0, 0.0, -0.58399436032505114, -3.518975543532981, 0.0, 0.0,
    -0.66661933005713259, -3.7321740417494929, 0.0, 0.0, -0.7545937186454571,
    -3.9445971145321033, 0.0, 0.0, -0.84789806968919446, -4.1562475821849354,
    0.0, 0.0, -0.94651299755249707, -4.3671282547543688, 0.0, 0.0,
    -1.0504191871071205, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.030657969889276472,
    -0.0, -0.0, -0.0, 0.061204433507154259, -0.0, -0.0, 0.0007664492472319118,
    0.091639796413708147, -0.0, -0.0, 0.0022965600849107681, 0.12196446269394814,
    -0.0, -0.0, 0.0045875549952534723, 0.15217883496318446, -0.0, -0.0,
    0.0076366665626021769, 0.18228331437237288, -0.0, -0.0, 0.011441137436681788,
    0.21227830061344088, -0.0, -0.0, 0.015998220295991111, 0.24216419192459424,
    -0.0, -0.0, 0.021305177811327132, 0.27194138509560428, -0.0, -0.0,
    0.027359282609441989, 0.30161027547307617, -0.0, -0.0, 0.034157817236832096,
    0.33117125696569771, -0.0, -0.0, 0.041698074123659, 0.36062472204946922,
    -0.0, -0.0, 0.049977355547801444, 0.38997106177291435, -0.0, -0.0,
    0.058992973599038177, 0.41921066576227212, -0.0, -0.0, 0.068742250143361036,
    0.44834392222666974, -0.0, -0.0, 0.079222516787417846, 0.47737121796327686,
    -0.0, -0.0, 0.0904311148430846, 0.506292938362441, -0.0, -0.0,
    0.10236539529216652, 0.53510946741280452, -0.0, -0.0, 0.11502271875122755,
    0.5638211877064021, -0.0, -0.0, 0.12840045543654766, 0.5924284804437413,
    -0.0, -0.0, 0.14249598512920772, -0.030657969889276472, 0.0, 0.0, 0.0,
    -0.061204433507154259, 0.0, 0.0, -0.0007664492472319118,
    -0.091639796413708147, 0.0, 0.0, -0.0022965600849107681,
    -0.12196446269394814, 0.0, 0.0, -0.0045875549952534723, -0.15217883496318446,
    0.0, 0.0, -0.0076366665626021769, -0.18228331437237288, 0.0, 0.0,
    -0.011441137436681788, -0.21227830061344088, 0.0, 0.0, -0.015998220295991111,
    -0.24216419192459424, 0.0, 0.0, -0.021305177811327132, -0.27194138509560428,
    0.0, 0.0, -0.027359282609441989, -0.30161027547307617, 0.0, 0.0,
    -0.034157817236832096, -0.33117125696569771, 0.0, 0.0, -0.041698074123659,
    -0.36062472204946922, 0.0, 0.0, -0.049977355547801444, -0.38997106177291435,
    0.0, 0.0, -0.058992973599038177, -0.41921066576227212, 0.0, 0.0,
    -0.068742250143361036, -0.44834392222666974, 0.0, 0.0, -0.079222516787417846,
    -0.47737121796327686, 0.0, 0.0, -0.0904311148430846, -0.506292938362441, 0.0,
    0.0, -0.10236539529216652, -0.53510946741280452, 0.0, 0.0,
    -0.11502271875122755, -0.5638211877064021, 0.0, 0.0, -0.12840045543654766,
    -0.5924284804437413, 0.0, 0.0, -0.14249598512920772, -0.0, -0.0, -1.0, 0.0,
    0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T a_0[618]{ -0.0, 0.1888968412469153, 0.060492001339316988,
    -0.0, -0.0, -0.0, 0.3773313347058973, 0.1209347881183509,
    0.0047224210311728829, 0.0015123000334829248, -0.0, 0.56530446924504474,
    0.18132836030302127, 0.014155704398820317, 0.0045356697364416975, -0.0,
    0.75281723176956472, 0.24167271798767209, 0.028288316129946436,
    0.0090688787440172287, -0.0, 0.939870607225484, 0.30196786139466048,
    0.047108746924185554, 0.015110696693709031, -0.0, 1.1264655786033526,
    0.36221379087394617, 0.07060551210482266, 0.022659893228575541, -0.0,
    1.3126031269419416, 0.42241050690268223, 0.098767151569906481,
    0.031715238000424195, -0.0, 1.4982842313319336, 0.48255801008480659,
    0.13158222974345501, 0.042275500672991254, -0.0, 1.6835098689196064,
    0.54265630115063446, 0.16903933552675335, 0.054339450925111416, -0.0,
    1.8682810149105107, 0.6027053809564521, 0.2111270822497435,
    0.067905858453877282, -0.0, 2.05259864257314, 0.66270525048411089,
    0.25783410762250625, 0.082973492977788582, -0.0, 2.2364637232425952,
    0.72265591084062319, 0.30914907368683475, 0.099541124239891351, -0.0,
    2.41987722632424, 0.78255736325775849, 0.36506066676789961,
    0.11760752201090693, -0.0, 2.6028401192973529, 0.842409609091641,
    0.42555759742600563, 0.1371714560923509, -0.0, 2.7853533677187707,
    0.90221264982234783, 0.49062860040843947, 0.15823169631964193, -0.0,
    2.9674179352265249, 0.96196648705350851, 0.56026243460140879,
    0.18078701256520063, -0.0, 3.1490347835434722, 1.0216711225119051,
    0.63444788298207189, 0.20483617474153834, -0.0, 3.3302048724809192,
    1.081326558047073, 0.71317375257065874, 0.23037795280433598, -0.0,
    3.5109291599422385, 1.1409327956309046, 0.79642887438268173,
    0.25741111675551281, -0.0, 3.69120860192648, 1.2004898373572503,
    0.88420210338123773, 0.28593443664628543, 0.0, -0.1888968412469153,
    -0.060492001339316988, 0.0, 0.0, 0.0, -0.3773313347058973,
    -0.1209347881183509, -0.0047224210311728829, -0.0015123000334829248, 0.0,
    -0.56530446924504474, -0.18132836030302127, -0.014155704398820317,
    -0.0045356697364416975, 0.0, -0.75281723176956472, -0.24167271798767209,
    -0.028288316129946436, -0.0090688787440172287, 0.0, -0.939870607225484,
    -0.30196786139466048, -0.047108746924185554, -0.015110696693709031, 0.0,
    -1.1264655786033526, -0.36221379087394617, -0.07060551210482266,
    -0.022659893228575541, 0.0, -1.3126031269419416, -0.42241050690268223,
    -0.098767151569906481, -0.031715238000424195, 0.0, -1.4982842313319336,
    -0.48255801008480659, -0.13158222974345501, -0.042275500672991254, 0.0,
    -1.6835098689196064, -0.54265630115063446, -0.16903933552675335,
    -0.054339450925111416, 0.0, -1.8682810149105107, -0.6027053809564521,
    -0.2111270822497435, -0.067905858453877282, 0.0, -2.05259864257314,
    -0.66270525048411089, -0.25783410762250625, -0.082973492977788582, 0.0,
    -2.2364637232425952, -0.72265591084062319, -0.30914907368683475,
    -0.099541124239891351, 0.0, -2.41987722632424, -0.78255736325775849,
    -0.36506066676789961, -0.11760752201090693, 0.0, -2.6028401192973529,
    -0.842409609091641, -0.42555759742600563, -0.1371714560923509, 0.0,
    -2.7853533677187707, -0.90221264982234783, -0.49062860040843947,
    -0.15823169631964193, 0.0, -2.9674179352265249, -0.96196648705350851,
    -0.56026243460140879, -0.18078701256520063, 0.0, -3.1490347835434722,
    -1.0216711225119051, -0.63444788298207189, -0.20483617474153834, 0.0,
    -3.3302048724809192, -1.081326558047073, -0.71317375257065874,
    -0.23037795280433598, 0.0, -3.5109291599422385, -1.1409327956309046,
    -0.79642887438268173, -0.25741111675551281, 0.0, -3.69120860192648,
    -1.2004898373572503, -0.88420210338123773, -0.28593443664628543, -1.0, -0.0,
    -0.0, 1.0, 0.0, 0.0, -0.0, -0.1599184353223182, 0.089801925435736424, -0.0,
    -0.0, -0.0, -0.31964914489275253, 0.17947362197539324,
    -0.0039979608830579549, 0.0022450481358934106, -0.0, -0.47919231323908462,
    0.26901527273781028, -0.011989189505376769, 0.0067318886852782414, -0.0,
    -0.63854812478731549, 0.35842706059465185, -0.023968997336353881,
    0.013457270503723497, -0.0, -0.797716763861511, 0.44770916817072121,
    -0.039932700456036765, 0.022417947018589793, -0.0, -0.956698414683647,
    0.53686177784427469, -0.059875619552574547, 0.033610676222857822, -0.0,
    -1.1154932613734565, 0.62588507174733543, -0.083793079919665725,
    0.047032220668964689, -0.0, -1.2741014879482764, 0.71477923176600711,
    -0.11168041145400215, 0.062679347462648069, -0.0, -1.4325232783228956,
    0.80354443954078658, -0.14353294865270905, 0.080548828256798238, -0.0,
    -1.5907588163094044, 0.89218087646687683, -0.17934603061078144,
    0.1006374392453179, -0.0, -1.748808285617043, 0.98068872369449933,
    -0.21911500101851655, 0.12294196115698981, -0.0, -1.9066718698520528,
    1.069068162129206, -0.26283520815894262, 0.14745917924935228, -0.0,
    -2.0643497525175265, 1.1573193724321906, -0.31050200490524393,
    0.17418588330258244, -0.0, -2.221842117013261, 1.2454425350206004,
    -0.36211074871818211, 0.20311886761338721, -0.0, -2.3791491466356085,
    1.3334378300678469, -0.41765680164351365, 0.23425493098890221, -0.0,
    -2.536271024577331, 1.421305437503916, -0.47713553030940387,
    0.26759087674059839, -0.0, -2.6932079339274537, 1.5090455370156788,
    -0.54054230592383712, 0.30312351267819632, -0.0, -2.8499600576711197,
    1.596658308047201, -0.60787250427202344, 0.34084965110358828, -0.0,
    -3.0065275786894463, 1.6841439298000525, -0.67912150571380148,
    0.3807661088047683, -0.0, -3.16291067975938, 1.7715025812336165,
    -0.75428469518103769, 0.42286970704976962, 0.0, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.0, 0.0, 0.31964914489275253,
    -0.17947362197539324, 0.0039979608830579549, -0.0022450481358934106, 0.0,
    0.47919231323908462, -0.26901527273781028, 0.011989189505376769,
    -0.0067318886852782414, 0.0, 0.63854812478731549, -0.35842706059465185,
    0.023968997336353881, -0.013457270503723497, 0.0, 0.797716763861511,
    -0.44770916817072121, 0.039932700456036765, -0.022417947018589793, 0.0,
    0.956698414683647, -0.53686177784427469, 0.059875619552574547,
    -0.033610676222857822, 0.0, 1.1154932613734565, -0.62588507174733543,
    0.083793079919665725, -0.047032220668964689, 0.0, 1.2741014879482764,
    -0.71477923176600711, 0.11168041145400215, -0.062679347462648069, 0.0,
    1.4325232783228956, -0.80354443954078658, 0.14353294865270905,
    -0.080548828256798238, 0.0, 1.5907588163094044, -0.89218087646687683,
    0.17934603061078144, -0.1006374392453179, 0.0, 1.748808285617043,
    -0.98068872369449933, 0.21911500101851655, -0.12294196115698981, 0.0,
    1.9066718698520528, -1.069068162129206, 0.26283520815894262,
    -0.14745917924935228, 0.0, 2.0643497525175265, -1.1573193724321906,
    0.31050200490524393, -0.17418588330258244, 0.0, 2.221842117013261,
    -1.2454425350206004, 0.36211074871818211, -0.20311886761338721, 0.0,
    2.3791491466356085, -1.3334378300678469, 0.41765680164351365,
    -0.23425493098890221, 0.0, 2.536271024577331, -1.421305437503916,
    0.47713553030940387, -0.26759087674059839, 0.0, 2.6932079339274537,
    -1.5090455370156788, 0.54054230592383712, -0.30312351267819632, 0.0,
    2.8499600576711197, -1.596658308047201, 0.60787250427202344,
    -0.34084965110358828, 0.0, 3.0065275786894463, -1.6841439298000525,
    0.67912150571380148, -0.3807661088047683, 0.0, 3.16291067975938,
    -1.7715025812336165, 0.75428469518103769, -0.42286970704976962, -0.0, -1.0,
    -0.0, 0.0, 1.0, 0.0, -0.0, 0.10536597514810161, -0.1119360731420971, -0.0,
    -0.0, -0.0, 0.21068448591173913, -0.22372201869524067, 0.0026341493787025404,
    -0.0027984018285524275, -0.0, 0.31595540965099317, -0.3353580258245521,
    0.0079012615264960184, -0.008391452295933445, -0.0, 0.4211786242426373,
    -0.446844283480486, 0.015800146767770848, -0.016775402941547252, -0.0,
    0.52635400807880228, -0.55818098039902453, 0.026329612373836783,
    -0.027946510028559402, -0.0, 0.63148144006564322, -0.6693683051018724,
    0.039488462575806842, -0.041901034538535019, -0.0, 0.73656079962200982,
    -0.78040644589665065, 0.055275498577447929, -0.058635242166081826, -0.0,
    0.84159196667811931, -0.89129559087709187, 0.073689518567998186,
    -0.078145403313498091, -0.0, 0.94657482167423235, -1.0020359279232338,
    0.094729317734951177, -0.10042779308542539, -0.0, 1.051509245559332,
    -1.1126276447016141, 0.11839368827680699, -0.12547869128350625, -0.0,
    1.1563951197898057, -1.2230709286654644, 0.1446814194157903,
    -0.15329438240104659, -0.0, 1.2612323263281291, -1.3333659670549041,
    0.17359129741053544, -0.18387115561768319, -0.0, 1.3660207476415547,
    -1.4435129468971351, 0.20512210556873869, -0.21720530479405578, -0.0,
    1.4707602667008017, -1.5535120550066355, 0.23927262425977758,
    -0.25329312846648416, -0.0, 1.5754507669787494, -1.6633634779853532,
    0.27604163092729761, -0.29213092984165007, -0.0, 1.6800921324491334,
    -1.7730674022229005, 0.31542790010176636, -0.3337150167912839, -0.0,
    1.7846842475852445, -1.8826240138967472, 0.3574302034129947,
    -0.37804170184685643, -0.0, 1.8892269973586309, -1.9920334989724147,
    0.4020473096026258, -0.42510730219427512, -0.0, 1.9937202672378032,
    -2.1012960432036696, 0.44927798453659162, -0.47490813966858547, -0.0,
    2.098163943186941, -2.2104118321327175, 0.49912099121753672,
    -0.52744054074867719, 0.0, -0.10536597514810161, 0.1119360731420971, 0.0,
    0.0, 0.0, -0.21068448591173913, 0.22372201869524067, -0.0026341493787025404,
    0.0027984018285524275, 0.0, -0.31595540965099317, 0.3353580258245521,
    -0.0079012615264960184, 0.008391452295933445, 0.0, -0.4211786242426373,
    0.446844283480486, -0.015800146767770848, 0.016775402941547252, 0.0,
    -0.52635400807880228, 0.55818098039902453, -0.026329612373836783,
    0.027946510028559402, 0.0, -0.63148144006564322, 0.6693683051018724,
    -0.039488462575806842, 0.041901034538535019, 0.0, -0.73656079962200982,
    0.78040644589665065, -0.055275498577447929, 0.058635242166081826, 0.0,
    -0.84159196667811931, 0.89129559087709187, -0.073689518567998186,
    0.078145403313498091, 0.0, -0.94657482167423235, 1.0020359279232338,
    -0.094729317734951177, 0.10042779308542539, 0.0, -1.051509245559332,
    1.1126276447016141, -0.11839368827680699, 0.12547869128350625, 0.0,
    -1.1563951197898057, 1.2230709286654644, -0.1446814194157903,
    0.15329438240104659, 0.0, -1.2612323263281291, 1.3333659670549041,
    -0.17359129741053544, 0.18387115561768319, 0.0, -1.3660207476415547,
    1.4435129468971351, -0.20512210556873869, 0.21720530479405578, 0.0,
    -1.4707602667008017, 1.5535120550066355, -0.23927262425977758,
    0.25329312846648416, 0.0, -1.5754507669787494, 1.6633634779853532,
    -0.27604163092729761, 0.29213092984165007, 0.0, -1.6800921324491334,
    1.7730674022229005, -0.31542790010176636, 0.3337150167912839, 0.0,
    -1.7846842475852445, 1.8826240138967472, -0.3574302034129947,
    0.37804170184685643, 0.0, -1.8892269973586309, 1.9920334989724147,
    -0.4020473096026258, 0.42510730219427512, 0.0, -1.9937202672378032,
    2.1012960432036696, -0.44927798453659162, 0.47490813966858547, 0.0,
    -2.098163943186941, 2.2104118321327175, -0.49912099121753672,
    0.52744054074867719, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T b_Ac_1[504]{ -0.12647183082737606, 0.14396327290074137,
    -0.0, -0.25266713091650589, 0.28784633762178591, -0.0, -0.37858657539014817,
    0.43164948302538075, -0.0, -0.50423083775077449, 0.57537299731355918, -0.0,
    -0.62960058988446832, 0.71901716802974069, -0.0, -0.75469650206481409,
    0.862582282060328, -0.0, -0.87951924295677675, 1.0060686256362998, -0.0,
    -1.0040694796205727, 1.1494764843347998, -0.0, -1.1283478775155311,
    1.2928061430807216, -0.0, -1.2523551005039453, 1.4360578861482902, -0.0,
    -1.3760918108549158, 1.57923199716264, -0.0, -1.4995586692481842,
    1.7223287591013881, -0.0, -1.6227563347779563, 1.8653484542962038, -0.0,
    -1.7456854649567182, 2.0082913644343763, -0.0, -1.8683467157190412,
    2.1511577705603755, -0.0, -1.9907407414253795, 2.2939479530774123, -0.0,
    -2.112868194865857, 2.4366621917489919, -0.0, -2.2347297272640461,
    2.5793007657004661, -0.0, -2.3563259882807368, 2.7218639534205811, -0.0,
    -2.477657626017697, 2.86435203276302, -0.0, 0.12647183082737606,
    -0.14396327290074137, 0.0, 0.25266713091650589, -0.28784633762178591, 0.0,
    0.37858657539014817, -0.43164948302538075, 0.0, 0.50423083775077449,
    -0.57537299731355918, 0.0, 0.62960058988446832, -0.71901716802974069, 0.0,
    0.75469650206481409, -0.862582282060328, 0.0, 0.87951924295677675,
    -1.0060686256362998, 0.0, 1.0040694796205727, -1.1494764843347998, 0.0,
    1.1283478775155311, -1.2928061430807216, 0.0, 1.2523551005039453,
    -1.4360578861482902, 0.0, 1.3760918108549158, -1.57923199716264, 0.0,
    1.4995586692481842, -1.7223287591013881, 0.0, 1.6227563347779563,
    -1.8653484542962038, 0.0, 1.7456854649567182, -2.0082913644343763, 0.0,
    1.8683467157190412, -2.1511577705603755, 0.0, 1.9907407414253795,
    -2.2939479530774123, 0.0, 2.112868194865857, -2.4366621917489919, 0.0,
    2.2347297272640461, -2.5793007657004661, 0.0, 2.3563259882807368,
    -2.7218639534205811, 0.0, 2.477657626017697, -2.86435203276302, 0.0, -1.0,
    -0.0, -0.0, 1.0, 0.0, 0.0, 0.18891678043810189, -0.2603991729373682, -0.0,
    0.37743412338800142, -0.52070047841969547, -0.0, 0.56555300944331421,
    -0.78090434240065543, -0.0, 0.7532744168461627, -1.0410111898780745, -0.0,
    0.94059932149283521, -1.3010214448962558, -0.0, 1.1275286969394305,
    -1.560935530548297, -0.0, 1.3140635144074895, -1.8207538689784042, -0.0,
    1.5002047427896132, -2.0804768813841976, -0.0, 1.6859533486550675,
    -2.3401049880190148, -0.0, 1.8713102962553732, -2.5996386081942049, -0.0,
    2.0562765475298841, -2.8590781602814213, -0.0, 2.2408530621113516,
    -3.1184240617149053, -0.0, 2.4250407973314743, -3.3776767289937664, -0.0,
    2.6088407082264355, -3.6368365776842562, -0.0, 2.7922537475424285,
    -3.8959040224220369, -0.0, 2.975280865741166, -4.154879476914445, -0.0,
    3.1579230110053778, -4.4137633539427483, -0.0, 3.3401811292442956,
    -4.6725560653643994, -0.0, 3.5220561640991237, -4.9312580221152817, -0.0,
    3.7035490569484968, -5.1898696342119521, -0.0, -0.18891678043810189,
    0.2603991729373682, 0.0, -0.37743412338800142, 0.52070047841969547, 0.0,
    -0.56555300944331421, 0.78090434240065543, 0.0, -0.7532744168461627,
    1.0410111898780745, 0.0, -0.94059932149283521, 1.3010214448962558, 0.0,
    -1.1275286969394305, 1.560935530548297, 0.0, -1.3140635144074895,
    1.8207538689784042, 0.0, -1.5002047427896132, 2.0804768813841976, 0.0,
    -1.6859533486550675, 2.3401049880190148, 0.0, -1.8713102962553732,
    2.5996386081942049, 0.0, -2.0562765475298841, 2.8590781602814213, 0.0,
    -2.2408530621113516, 3.1184240617149053, 0.0, -2.4250407973314743,
    3.3776767289937664, 0.0, -2.6088407082264355, 3.6368365776842562, 0.0,
    -2.7922537475424285, 3.8959040224220369, 0.0, -2.975280865741166,
    4.154879476914445, 0.0, -3.1579230110053778, 4.4137633539427483, 0.0,
    -3.3401811292442956, 4.6725560653643994, 0.0, -3.5220561640991237,
    4.9312580221152817, 0.0, -3.7035490569484968, 5.1898696342119521, 0.0, -0.0,
    -1.0, -0.0, 0.0, 1.0, 0.0, 0.13919179616088076, 0.07490340326547, -0.0,
    0.27800912838762087, 0.15000797759616485, -0.0, 0.41645288308364015,
    0.22531337660988532, -0.0, 0.55452394451511977, 0.30081925480716887, -0.0,
    0.69222319481614147, 0.37652526756918425, -0.0, 0.82955151399381344,
    0.4524310711556313, -0.0, 0.96650977993338516, 0.52853632270264572, -0.0,
    1.1030988684033489, 0.60484068022070858, -0.0, 1.2393196530605293,
    0.68134380259256155, -0.0, 1.375173005455161, 0.75804534957112635, -0.0,
    1.5106597950359535, 0.83494498177742982, -0.0, 1.6457808891551446,
    0.91204236069853384, -0.0, 1.7805371530735408, 0.98933714868547, -0.0,
    1.914929449965546, 1.0668290089511794, -0.0, 2.048958640924178,
    1.1445176055684576, -0.0, 2.1826255849660736, 1.2224026034679043, -0.0,
    2.315931139036481, 1.3004836684358778, -0.0, 2.4488761580142389,
    1.3787604671124549, -0.0, 2.5814614947167467, 1.4572326669893951, -0.0,
    2.7136879999049208, 1.53589993640811, -0.0, -0.13919179616088076,
    -0.07490340326547, 0.0, -0.27800912838762087, -0.15000797759616485, 0.0,
    -0.41645288308364015, -0.22531337660988532, 0.0, -0.55452394451511977,
    -0.30081925480716887, 0.0, -0.69222319481614147, -0.37652526756918425, 0.0,
    -0.82955151399381344, -0.4524310711556313, 0.0, -0.96650977993338516,
    -0.52853632270264572, 0.0, -1.1030988684033489, -0.60484068022070858, 0.0,
    -1.2393196530605293, -0.68134380259256155, 0.0, -1.375173005455161,
    -0.75804534957112635, 0.0, -1.5106597950359535, -0.83494498177742982, 0.0,
    -1.6457808891551446, -0.91204236069853384, 0.0, -1.7805371530735408,
    -0.98933714868547, 0.0, -1.914929449965546, -1.0668290089511794, 0.0,
    -2.048958640924178, -1.1445176055684576, 0.0, -2.1826255849660736,
    -1.2224026034679043, 0.0, -2.315931139036481, -1.3004836684358778, 0.0,
    -2.4488761580142389, -1.3787604671124549, 0.0, -2.5814614947167467,
    -1.4572326669893951, 0.0, -2.7136879999049208, -1.53589993640811, 0.0, -0.0,
    -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T a[498]{ -0.20979775256106059, -0.0, -0.0, -0.0,
    -0.41883244855900248, -0.0, -0.0, -0.0052449438140265148,
    -0.62710686331106058, -0.0, -0.0, -0.015715755028001578, -0.834623762040348,
    -0.0, -0.0, -0.0313934266107781, -1.04138589991257, -0.0, -0.0,
    -0.052259020661786795, -1.247396022072603, -0.0, -0.0, -0.078293668159601051,
    -1.452656863680942, -0.0, -0.0, -0.10947856871141612, -1.6571711499500148,
    -0.0, -0.0, -0.14579499030343968, -1.8609415961803639, -0.0, -0.0,
    -0.18722426905219003, -2.0639709077966972, -0.0, -0.0, -0.23374780895669911,
    -2.2662617803838079, -0.0, -0.0, -0.28534708165161654, -2.4678168997223624,
    -0.0, -0.0, -0.34200362616121172, -2.6686389418245593, -0.0, -0.0,
    -0.40369904865427075, -2.8687305729696591, -0.0, -0.0, -0.47041502219988474,
    -3.0680944497393825, -0.0, -0.0, -0.54213328652412618, -3.266733219053183,
    -0.0, -0.0, -0.6188356477676108, -3.4646495182033883, -0.0, -0.0,
    -0.7005039782439404, -3.6618459748902157, -0.0, -0.0, -0.78712021619902517,
    -3.8583252072566592, -0.0, -0.0, -0.87866636557128064, -4.0540898239232508,
    -0.0, -0.0, -0.97512449575269711, 0.20979775256106059, 0.0, 0.0, 0.0,
    0.41883244855900248, 0.0, 0.0, 0.0052449438140265148, 0.62710686331106058,
    0.0, 0.0, 0.015715755028001578, 0.834623762040348, 0.0, 0.0,
    0.0313934266107781, 1.04138589991257, 0.0, 0.0, 0.052259020661786795,
    1.247396022072603, 0.0, 0.0, 0.078293668159601051, 1.452656863680942, 0.0,
    0.0, 0.10947856871141612, 1.6571711499500148, 0.0, 0.0, 0.14579499030343968,
    1.8609415961803639, 0.0, 0.0, 0.18722426905219003, 2.0639709077966972, 0.0,
    0.0, 0.23374780895669911, 2.2662617803838079, 0.0, 0.0, 0.28534708165161654,
    2.4678168997223624, 0.0, 0.0, 0.34200362616121172, 2.6686389418245593, 0.0,
    0.0, 0.40369904865427075, 2.8687305729696591, 0.0, 0.0, 0.47041502219988474,
    3.0680944497393825, 0.0, 0.0, 0.54213328652412618, 3.266733219053183, 0.0,
    0.0, 0.6188356477676108, 3.4646495182033883, 0.0, 0.0, 0.7005039782439404,
    3.6618459748902157, 0.0, 0.0, 0.78712021619902517, 3.8583252072566592, 0.0,
    0.0, 0.87866636557128064, 4.0540898239232508, 0.0, 0.0, 0.97512449575269711,
    -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.22599738357714264, -0.0, -0.0, -0.0,
    0.45117279082383815, -0.0, -0.0, 0.0056499345894285664, 0.67552921135473221,
    -0.0, -0.0, 0.016929254360024522, 0.89906962391092649, -0.0, -0.0,
    0.033817484643892823, 1.1217969963995269, -0.0, -0.0, 0.056294225241665982,
    1.3437142859330482, -0.0, -0.0, 0.084339150151654163, 1.5648244388686749,
    -0.0, -0.0, 0.11793200729998038, 1.7851303908473797, -0.0, -0.0,
    0.15705261827169725, 2.0046350668328987, -0.0, -0.0, 0.20168087804288176,
    2.2233413811505671, -0.0, -0.0, 0.25179675471370422, 2.4412522375260108,
    -0.0, -0.0, 0.30738028924246841, 2.6583705291236992, -0.0, -0.0,
    0.36841159518061867, 2.8746991385853575, -0.0, -0.0, 0.43487085840871115,
    3.0902409380682387, -0.0, -0.0, 0.50673833687334513, 3.3049987892832569,
    -0.0, -0.0, 0.58399436032505114, 3.518975543532981, -0.0, -0.0,
    0.66661933005713259, 3.7321740417494929, -0.0, -0.0, 0.7545937186454571,
    3.9445971145321033, -0.0, -0.0, 0.84789806968919446, 4.1562475821849354,
    -0.0, -0.0, 0.94651299755249707, 4.3671282547543688, -0.0, -0.0,
    1.0504191871071205, -0.22599738357714264, 0.0, 0.0, 0.0,
    -0.45117279082383815, 0.0, 0.0, -0.0056499345894285664, -0.67552921135473221,
    0.0, 0.0, -0.016929254360024522, -0.89906962391092649, 0.0, 0.0,
    -0.033817484643892823, -1.1217969963995269, 0.0, 0.0, -0.056294225241665982,
    -1.3437142859330482, 0.0, 0.0, -0.084339150151654163, -1.5648244388686749,
    0.0, 0.0, -0.11793200729998038, -1.7851303908473797, 0.0, 0.0,
    -0.15705261827169725, -2.0046350668328987, 0.0, 0.0, -0.20168087804288176,
    -2.2233413811505671, 0.0, 0.0, -0.25179675471370422, -2.4412522375260108,
    0.0, 0.0, -0.30738028924246841, -2.6583705291236992, 0.0, 0.0,
    -0.36841159518061867, -2.8746991385853575, 0.0, 0.0, -0.43487085840871115,
    -3.0902409380682387, 0.0, 0.0, -0.50673833687334513, -3.3049987892832569,
    0.0, 0.0, -0.58399436032505114, -3.518975543532981, 0.0, 0.0,
    -0.66661933005713259, -3.7321740417494929, 0.0, 0.0, -0.7545937186454571,
    -3.9445971145321033, 0.0, 0.0, -0.84789806968919446, -4.1562475821849354,
    0.0, 0.0, -0.94651299755249707, -4.3671282547543688, 0.0, 0.0,
    -1.0504191871071205, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.030657969889276472,
    -0.0, -0.0, -0.0, 0.061204433507154259, -0.0, -0.0, 0.0007664492472319118,
    0.091639796413708147, -0.0, -0.0, 0.0022965600849107681, 0.12196446269394814,
    -0.0, -0.0, 0.0045875549952534723, 0.15217883496318446, -0.0, -0.0,
    0.0076366665626021769, 0.18228331437237288, -0.0, -0.0, 0.011441137436681788,
    0.21227830061344088, -0.0, -0.0, 0.015998220295991111, 0.24216419192459424,
    -0.0, -0.0, 0.021305177811327132, 0.27194138509560428, -0.0, -0.0,
    0.027359282609441989, 0.30161027547307617, -0.0, -0.0, 0.034157817236832096,
    0.33117125696569771, -0.0, -0.0, 0.041698074123659, 0.36062472204946922,
    -0.0, -0.0, 0.049977355547801444, 0.38997106177291435, -0.0, -0.0,
    0.058992973599038177, 0.41921066576227212, -0.0, -0.0, 0.068742250143361036,
    0.44834392222666974, -0.0, -0.0, 0.079222516787417846, 0.47737121796327686,
    -0.0, -0.0, 0.0904311148430846, 0.506292938362441, -0.0, -0.0,
    0.10236539529216652, 0.53510946741280452, -0.0, -0.0, 0.11502271875122755,
    0.5638211877064021, -0.0, -0.0, 0.12840045543654766, 0.5924284804437413,
    -0.0, -0.0, 0.14249598512920772, -0.030657969889276472, 0.0, 0.0, 0.0,
    -0.061204433507154259, 0.0, 0.0, -0.0007664492472319118,
    -0.091639796413708147, 0.0, 0.0, -0.0022965600849107681,
    -0.12196446269394814, 0.0, 0.0, -0.0045875549952534723, -0.15217883496318446,
    0.0, 0.0, -0.0076366665626021769, -0.18228331437237288, 0.0, 0.0,
    -0.011441137436681788, -0.21227830061344088, 0.0, 0.0, -0.015998220295991111,
    -0.24216419192459424, 0.0, 0.0, -0.021305177811327132, -0.27194138509560428,
    0.0, 0.0, -0.027359282609441989, -0.30161027547307617, 0.0, 0.0,
    -0.034157817236832096, -0.33117125696569771, 0.0, 0.0, -0.041698074123659,
    -0.36062472204946922, 0.0, 0.0, -0.049977355547801444, -0.38997106177291435,
    0.0, 0.0, -0.058992973599038177, -0.41921066576227212, 0.0, 0.0,
    -0.068742250143361036, -0.44834392222666974, 0.0, 0.0, -0.079222516787417846,
    -0.47737121796327686, 0.0, 0.0, -0.0904311148430846, -0.506292938362441, 0.0,
    0.0, -0.10236539529216652, -0.53510946741280452, 0.0, 0.0,
    -0.11502271875122755, -0.5638211877064021, 0.0, 0.0, -0.12840045543654766,
    -0.5924284804437413, 0.0, 0.0, -0.14249598512920772, -0.0, -0.0, -1.0, 0.0,
    0.0, 1.0 };

  static const real_T a_1[378]{ -0.12647183082737606, 0.14396327290074137, -0.0,
    -0.25266713091650589, 0.28784633762178591, -0.0, -0.37858657539014817,
    0.43164948302538075, -0.0, -0.50423083775077449, 0.57537299731355918, -0.0,
    -0.62960058988446832, 0.71901716802974069, -0.0, -0.75469650206481409,
    0.862582282060328, -0.0, -0.87951924295677675, 1.0060686256362998, -0.0,
    -1.0040694796205727, 1.1494764843347998, -0.0, -1.1283478775155311,
    1.2928061430807216, -0.0, -1.2523551005039453, 1.4360578861482902, -0.0,
    -1.3760918108549158, 1.57923199716264, -0.0, -1.4995586692481842,
    1.7223287591013881, -0.0, -1.6227563347779563, 1.8653484542962038, -0.0,
    -1.7456854649567182, 2.0082913644343763, -0.0, -1.8683467157190412,
    2.1511577705603755, -0.0, -1.9907407414253795, 2.2939479530774123, -0.0,
    -2.112868194865857, 2.4366621917489919, -0.0, -2.2347297272640461,
    2.5793007657004661, -0.0, -2.3563259882807368, 2.7218639534205811, -0.0,
    -2.477657626017697, 2.86435203276302, -0.0, 0.12647183082737606,
    -0.14396327290074137, 0.0, 0.25266713091650589, -0.28784633762178591, 0.0,
    0.37858657539014817, -0.43164948302538075, 0.0, 0.50423083775077449,
    -0.57537299731355918, 0.0, 0.62960058988446832, -0.71901716802974069, 0.0,
    0.75469650206481409, -0.862582282060328, 0.0, 0.87951924295677675,
    -1.0060686256362998, 0.0, 1.0040694796205727, -1.1494764843347998, 0.0,
    1.1283478775155311, -1.2928061430807216, 0.0, 1.2523551005039453,
    -1.4360578861482902, 0.0, 1.3760918108549158, -1.57923199716264, 0.0,
    1.4995586692481842, -1.7223287591013881, 0.0, 1.6227563347779563,
    -1.8653484542962038, 0.0, 1.7456854649567182, -2.0082913644343763, 0.0,
    1.8683467157190412, -2.1511577705603755, 0.0, 1.9907407414253795,
    -2.2939479530774123, 0.0, 2.112868194865857, -2.4366621917489919, 0.0,
    2.2347297272640461, -2.5793007657004661, 0.0, 2.3563259882807368,
    -2.7218639534205811, 0.0, 2.477657626017697, -2.86435203276302, 0.0, -1.0,
    -0.0, -0.0, 1.0, 0.0, 0.0, 0.18891678043810189, -0.2603991729373682, -0.0,
    0.37743412338800142, -0.52070047841969547, -0.0, 0.56555300944331421,
    -0.78090434240065543, -0.0, 0.7532744168461627, -1.0410111898780745, -0.0,
    0.94059932149283521, -1.3010214448962558, -0.0, 1.1275286969394305,
    -1.560935530548297, -0.0, 1.3140635144074895, -1.8207538689784042, -0.0,
    1.5002047427896132, -2.0804768813841976, -0.0, 1.6859533486550675,
    -2.3401049880190148, -0.0, 1.8713102962553732, -2.5996386081942049, -0.0,
    2.0562765475298841, -2.8590781602814213, -0.0, 2.2408530621113516,
    -3.1184240617149053, -0.0, 2.4250407973314743, -3.3776767289937664, -0.0,
    2.6088407082264355, -3.6368365776842562, -0.0, 2.7922537475424285,
    -3.8959040224220369, -0.0, 2.975280865741166, -4.154879476914445, -0.0,
    3.1579230110053778, -4.4137633539427483, -0.0, 3.3401811292442956,
    -4.6725560653643994, -0.0, 3.5220561640991237, -4.9312580221152817, -0.0,
    3.7035490569484968, -5.1898696342119521, -0.0, -0.18891678043810189,
    0.2603991729373682, 0.0, -0.37743412338800142, 0.52070047841969547, 0.0,
    -0.56555300944331421, 0.78090434240065543, 0.0, -0.7532744168461627,
    1.0410111898780745, 0.0, -0.94059932149283521, 1.3010214448962558, 0.0,
    -1.1275286969394305, 1.560935530548297, 0.0, -1.3140635144074895,
    1.8207538689784042, 0.0, -1.5002047427896132, 2.0804768813841976, 0.0,
    -1.6859533486550675, 2.3401049880190148, 0.0, -1.8713102962553732,
    2.5996386081942049, 0.0, -2.0562765475298841, 2.8590781602814213, 0.0,
    -2.2408530621113516, 3.1184240617149053, 0.0, -2.4250407973314743,
    3.3776767289937664, 0.0, -2.6088407082264355, 3.6368365776842562, 0.0,
    -2.7922537475424285, 3.8959040224220369, 0.0, -2.975280865741166,
    4.154879476914445, 0.0, -3.1579230110053778, 4.4137633539427483, 0.0,
    -3.3401811292442956, 4.6725560653643994, 0.0, -3.5220561640991237,
    4.9312580221152817, 0.0, -3.7035490569484968, 5.1898696342119521, 0.0, -0.0,
    -1.0, -0.0, 0.0, 1.0, 0.0, 0.13919179616088076, 0.07490340326547, -0.0,
    0.27800912838762087, 0.15000797759616485, -0.0, 0.41645288308364015,
    0.22531337660988532, -0.0, 0.55452394451511977, 0.30081925480716887, -0.0,
    0.69222319481614147, 0.37652526756918425, -0.0, 0.82955151399381344,
    0.4524310711556313, -0.0, 0.96650977993338516, 0.52853632270264572, -0.0,
    1.1030988684033489, 0.60484068022070858, -0.0, 1.2393196530605293,
    0.68134380259256155, -0.0, 1.375173005455161, 0.75804534957112635, -0.0,
    1.5106597950359535, 0.83494498177742982, -0.0, 1.6457808891551446,
    0.91204236069853384, -0.0, 1.7805371530735408, 0.98933714868547, -0.0,
    1.914929449965546, 1.0668290089511794, -0.0, 2.048958640924178,
    1.1445176055684576, -0.0, 2.1826255849660736, 1.2224026034679043, -0.0,
    2.315931139036481, 1.3004836684358778, -0.0, 2.4488761580142389,
    1.3787604671124549, -0.0, 2.5814614947167467, 1.4572326669893951, -0.0,
    2.7136879999049208, 1.53589993640811, -0.0, -0.13919179616088076,
    -0.07490340326547, 0.0, -0.27800912838762087, -0.15000797759616485, 0.0,
    -0.41645288308364015, -0.22531337660988532, 0.0, -0.55452394451511977,
    -0.30081925480716887, 0.0, -0.69222319481614147, -0.37652526756918425, 0.0,
    -0.82955151399381344, -0.4524310711556313, 0.0, -0.96650977993338516,
    -0.52853632270264572, 0.0, -1.1030988684033489, -0.60484068022070858, 0.0,
    -1.2393196530605293, -0.68134380259256155, 0.0, -1.375173005455161,
    -0.75804534957112635, 0.0, -1.5106597950359535, -0.83494498177742982, 0.0,
    -1.6457808891551446, -0.91204236069853384, 0.0, -1.7805371530735408,
    -0.98933714868547, 0.0, -1.914929449965546, -1.0668290089511794, 0.0,
    -2.048958640924178, -1.1445176055684576, 0.0, -2.1826255849660736,
    -1.2224026034679043, 0.0, -2.315931139036481, -1.3004836684358778, 0.0,
    -2.4488761580142389, -1.3787604671124549, 0.0, -2.5814614947167467,
    -1.4572326669893951, 0.0, -2.7136879999049208, -1.53589993640811, 0.0, -0.0,
    -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T b_Kr_0[300]{ -0.0, 0.0034600074982263458,
    0.0011080268830073561, -0.0, -0.0, -0.0, 0.0069115462110432707,
    0.0022151523070677161, 0.0, 0.0, -0.0, 0.010354634251464278,
    0.0033213762715568304, 0.0, 0.0, -0.0, 0.013789289696548754,
    0.004426698778202791, 0.0, 0.0, -0.0, 0.017215530587469937,
    0.0055311198310784935, 0.0, 0.0, -0.0, 0.020633374929582771,
    0.0066346394365941214, 0.0, 0.0, -0.0, 0.024042840692491626,
    0.0077372576034896463, 0.0, 0.0, -0.0, 0.027443945810117897,
    0.00883897434282735, 0.0, 0.0, -0.0, 0.030836708180767489,
    0.0099397896679843528, 0.0, 0.0, -0.0, 0.034221145667198159,
    0.011039703594645178, 0.0, 0.0, -0.0, 0.037597276096686769,
    0.012138716140794316, 0.0, 0.0, -0.0, 0.040965117261096369,
    0.013236827326708818, 0.0, 0.0, -0.0, 0.04432468691694319,
    0.014334037174950902, 0.0, 0.0, -0.0, 0.047676002785463537,
    0.015430345710360577, 0.0, 0.0, -0.0, 0.051019082552680478,
    0.0165257529600483, 0.0, 0.0, -0.0, 0.054353943869470511,
    0.017620258953387605, 0.0, 0.0, -0.0, 0.057680604351630041,
    0.018713863722007822, 0.0, 0.0, -0.0, 0.060999081579941751,
    0.019806567299786734, 0.0, 0.0, -0.0, 0.064309393100240878,
    0.020898369722843323, 0.0, 0.0, -0.0, 0.067611556423481323,
    0.02198927102953048, 0.0, 0.0, -0.0, -0.0029292124826829609,
    0.0016448942889238771, -0.0, -0.0, -0.0, -0.0058549864086127187,
    0.0032874031861495827, -0.0, 0.0, -0.0, -0.0087773251577690746,
    0.0049275300458494512, -0.0, 0.0, -0.0, -0.011696232108267524,
    0.0065652782176683225, -0.0, 0.0, -0.0, -0.014611710636356426,
    0.008200651046729306, -0.0, 0.0, -0.0, -0.01752376411641416,
    0.0098336518736395276, -0.0, 0.0, -0.0, -0.020432395920946339,
    0.011464284034495885, -0.0, 0.0, -0.0, -0.023337609420582991,
    0.013092550860890789, -0.0, 0.0, -0.0, -0.026239407984075785,
    0.014718455679917888, -0.0, 0.0, -0.0, -0.029137794978295259,
    0.016342001814177807, -0.0, 0.0, -0.0, -0.032032773768228065,
    0.01796319258178386, -0.0, 0.0, -0.0, -0.034924347716974231,
    0.019582031296367777, -0.0, 0.0, -0.0, -0.037812520185744415,
    0.0211985212670854, -0.0, 0.0, -0.0, -0.040697294533857217,
    0.022812665798622379, -0.0, 0.0, -0.0, -0.043578674118736456,
    0.024424468191199891, -0.0, 0.0, -0.0, -0.046456662295908488,
    0.0260339317405803, -0.0, 0.0, -0.0, -0.049331262418999537,
    0.027641059738072859, -0.0, 0.0, -0.0, -0.052202477839733023,
    0.029245855470539376, -0.0, 0.0, -0.0, -0.055070311907926939,
    0.03084832222039988, -0.0, 0.0, -0.0, -0.057934767971491184,
    0.032448463265638293, -0.0, 0.0, -0.0, 0.0019299796738994745,
    -0.002050323604339219, -0.0, -0.0, -0.0, 0.0038590899466747143,
    -0.004097897334302345, 0.0, -0.0, -0.0, 0.0057873285719406663,
    -0.00614272465481094, 0.0, -0.0, -0.0, 0.0077146933127765,
    -0.0081848090268545354, 0.0, -0.0, -0.0, 0.0096411819417011376,
    -0.010224153907494185, 0.0, -0.0, -0.0, 0.011566792240648844,
    -0.012260762749866044, 0.0, -0.0, -0.0, 0.013491522000944865,
    -0.014294639003184914, 0.0, -0.0, -0.0, 0.015415369023281122,
    -0.01632578611274782, 0.0, -0.0, -0.0, 0.017338331117691962,
    -0.018354207519937556, 0.0, -0.0, -0.0, 0.019260406103529956,
    -0.020379906662226249, 0.0, -0.0, -0.0, 0.021181591809441756,
    -0.02240288697317893, 0.0, -0.0, -0.0, 0.023101886073343996,
    -0.024423151882457056, 0.0, -0.0, -0.0, 0.025021286742399253,
    -0.0264407048158221, 0.0, -0.0, -0.0, 0.026939791672992073,
    -0.028455549195139096, 0.0, -0.0, -0.0, 0.028857398730705017,
    -0.030467688438380169, 0.0, -0.0, -0.0, 0.030774105790294791,
    -0.032477125959628116, 0.0, -0.0, -0.0, 0.032689910735668425,
    -0.03448386516907994, 0.0, -0.0, -0.0, 0.034604811459859462,
    -0.0364879094730504, 0.0, -0.0, -0.0, 0.036518805865004278,
    -0.038489262273975561, 0.0, -0.0, -0.0, 0.038431891862318386,
    -0.040487926970416352, 0.0, -0.0 };

  static const real_T b_Kr_1[300]{ -0.0023165738510425252, 0.0026369631192226465,
    -0.0, -0.0, -0.0, -0.004628082511891788, 0.0052724570719873565, -0.0,
    -5.7914346276063136E-5, 6.5924077980566165E-5, -0.00693453834871438,
    0.00790648714935953, -0.0, -0.00017361640907335787, 0.00019773550478025009,
    -0.0092359536979982285, 0.010539058630311489, -0.0, -0.00034697986779121738,
    0.0003953976835142384, -0.011532340866624018, 0.013170176781751777, -0.0,
    -0.00057787871024117326, 0.00065887414927202567, -0.013823712131936424,
    0.015799846858554421, -0.0, -0.00086618723190677387, 0.00098812856881582,
    -0.016110079741815171, 0.0184280741035881, -0.0, -0.0012117800352051847,
    0.0013831247402796806, -0.018391455914745947, 0.021054863747745248, -0.0,
    -0.0016145320287505641, 0.0018438265928693827, -0.020667852839891119,
    0.023680221009971098, -0.0, -0.002074318426619213, 0.0023701981865630137,
    -0.022939282677160281, 0.026304151097292637, -0.0, -0.0025910147476164908,
    0.002962203711812291, -0.025205757557280654, 0.028926659204847514, -0.0,
    -0.003164496814545498, 0.0036198074892446067, -0.0274672895818673,
    0.031547750515912854, -0.0, -0.0037946407534775141, 0.0043429739693657942,
    -0.029723890823493165, 0.034167430201934017, -0.0, -0.004481322993024197,
    0.0051316677322636155, -0.031975573325758962, 0.036785703422553308, -0.0,
    -0.0052244202636115263, 0.0059858534873119661, -0.034222349103362866,
    0.03940257532563856, -0.0, -0.0060238095967555007, 0.006905496072875798,
    -0.036464230142170095, 0.042018051047311718, -0.0, -0.0068793683243395728,
    0.0078905604560167621, -0.03870122839928225, 0.044632135711977293, -0.0,
    -0.0077909740778938259, 0.0089410117321995566, -0.040933355803106544,
    0.047244834432350807, -0.0, -0.0087585047878758829, 0.010056815124998989,
    -0.043160624253424838, 0.04985615230948711, -0.0, -0.0097818386829535472,
    0.011237935985807759, -0.04538304562146251, 0.052466094432808669, -0.0,
    -0.010860854289289169, 0.012484339793544939, 0.0034603727227084428,
    -0.0047697096730035769, -0.0, -0.0, -0.0, 0.0069134289826580072,
    -0.0095376267160931817, -0.0, 8.650931806771108E-5, -0.00011924274182508943,
    0.010359186741299188, -0.014303758931426304, -0.0, 0.00025934504263416125,
    -0.00035768340972741896, 0.013797663917010379, -0.019068114103652261, -0.0,
    0.0005183247111666409, -0.0007152773830130765, 0.017228878385201526,
    -0.023830699999954765, -0.0, 0.00086326630909190052, -0.001191980235604383,
    0.020652847978417523, -0.028591524370094375, -0.0, 0.0012939882687219388,
    -0.0017877477356032523, 0.024069590486441365, -0.03335059494645088, -0.0,
    0.0018103094681823772, -0.0025025358448556114, 0.027479123656397049,
    -0.038107919444065552, -0.0, 0.0024120492303434113, -0.0033363007185168835,
    0.030881465192852241, -0.0428635055606833, -0.0, 0.0030990273217533378,
    -0.004288998704618522, 0.034276632757920658, -0.047617360976794712, -0.0,
    0.0038710639515746434, -0.005360586343635604, 0.037664643971364271,
    -0.05236949335567806, -0.0, 0.00472797977052266, -0.0065510203680554709,
    0.041045516410695176, -0.057119910343441105, -0.0, 0.0056695958698067669,
    -0.0078602577019474228, 0.044419267611277312, -0.06186861956906288, -0.0,
    0.0066957337800741465, -0.00928825546053345, 0.047785915066427841,
    -0.066615628644435348, -0.0, 0.00780621547035608, -0.010834970949760022,
    0.051145476227518361, -0.071360945164404943, -0.0, 0.009000863347016776,
    -0.012500361665870907, 0.054497968504075858, -0.076104576706814031, -0.0,
    0.010279500252704736, -0.014284385294981031, 0.057843409263883365,
    -0.080846530832542229, -0.0, 0.011641949465306632, -0.016186999712651383,
    0.061181815833080445, -0.085586815085547774, -0.0, 0.013088034696903717,
    -0.018208162983464938, 0.0645132054962634, -0.090325436992908542, -0.0,
    0.014617580092730729, -0.020347833360603633, 0.0678375954965852,
    -0.095062404064863187, -0.0, 0.016230410230137316, -0.022605969285426346,
    0.0025495643824912564, 0.0013719993157663782, -0.0, -0.0, -0.0,
    0.0050922697407056149, 0.002747683464955642, -0.0, 6.3739109562281415E-5,
    3.4299982894159452E-5, 0.0076281323108197035, 0.0041270461029142827, -0.0,
    0.00019104585307992176, 0.00010299206951805052, 0.010157168289862531,
    0.0055100809011578052, -0.0, 0.00038174916085041433, 0.00020616822209090755,
    0.012679393835809618, 0.0068967815473321642, -0.0, 0.00063567836809697768,
    0.00034392024461985264, 0.015194825067676897, 0.0082871417451752926, -0.0,
    0.00095266321399221823, 0.0005163397833031567, 0.017703478065614387,
    0.0096811552144787233, -0.0, 0.0013325338406841407, 0.000723518326932539,
    0.020205368870999645, 0.011078815691049307, -0.0, 0.0017751207923245005,
    0.00096554720729450706, 0.022700513486530991, 0.012480116926671008, -0.0,
    0.0022802550140994917, 0.0012425175995707397, 0.025188927876320519,
    0.013885052689066816, -0.0, 0.0028477678512627666, 0.0015545205227375149,
    0.027670627965986856, 0.015293616761860718, -0.0, 0.0034774910481707791,
    0.0019016468399641851, 0.030145629642747734, 0.0167058029445398, -0.0,
    0.0041692567473204508, 0.002283987259010703, 0.032613948755512323,
    0.0181216050524164, -0.0, 0.004922897488389145, 0.002701632332624198,
    0.035075601114973327, 0.019541016916590395, -0.0, 0.0057382462072769523,
    0.0031546724589346076, 0.03753060249369887, 0.020964032383911525, -0.0,
    0.0066151362351512854, 0.0036431978818493678, 0.039978968626224194,
    0.022390645316941866, -0.0, 0.0075534012974937574, 0.0041672986914471563,
    0.042420715209143083, 0.023820849593918354, -0.0, 0.008552875513149363,
    0.0047270648243707032, 0.044855857901199073, 0.025254639108715407, -0.0,
    0.00961339339337794, 0.0053225860642186623, 0.04728441232337649,
    0.026692007770807654, -0.0, 0.010734789840907918, 0.0059539520419365482,
    0.049706394058991236, 0.028132949505232714, -0.0, 0.011916900148992328,
    0.00662125223620674 };

  static const real_T b_Kr[240]{ -0.00384284772673063, -0.0, -0.0, -0.0,
    -0.0076717186107965888, -0.0, -0.0, -9.6071193168265748E-5,
    -0.011486663487449431, -0.0, -0.0, -0.00028786415843818052,
    -0.015287733007047537, -0.0, -0.0, -0.00057503074562441641,
    -0.019074977635728591, -0.0, -0.0, -0.00095722407080060475,
    -0.022848447656079605, -0.0, -0.0, -0.0014340985116938197,
    -0.026608193167804516, -0.0, -0.0, -0.0020053097030958097,
    -0.030354264088389361, -0.0, -0.0, -0.0026705145322909224,
    -0.034086710153765, -0.0, -0.0, -0.0034293711345006563,
    -0.037805580918967477, -0.0, -0.0, -0.0042815388883447811,
    -0.04151092575879594, -0.0, -0.0, -0.0052266784113189678,
    -0.04520279386846817, -0.0, -0.0, -0.006264451555288866,
    -0.048881234264273758, -0.0, -0.0, -0.00739452140200057,
    -0.052546295784224879, -0.0, -0.0, -0.0086165522586074145,
    -0.0561980270887047, -0.0, -0.0, -0.0099302096532130346,
    -0.059836476661113455, -0.0, -0.0, -0.011335160330430654,
    -0.063461692808512118, -0.0, -0.0, -0.012831072246958491,
    -0.067073723662263787, -0.0, -0.0, -0.014417614567171295,
    -0.070672617178672723, -0.0, -0.0, -0.016094457658727891,
    -0.074258421139621028, -0.0, -0.0, -0.017861273088194709,
    0.0041395750008033473, -0.0, -0.0, -0.0, 0.0082640939305366964, -0.0, -0.0,
    0.00010348937502008368, 0.01237361154971919, -0.0, -0.0,
    0.00031009172328350115, 0.016468182419700179, -0.0, -0.0,
    0.00061943201202648078, 0.020547860903383636, -0.0, -0.0,
    0.0010311365725189853, 0.024612701165949907, -0.0, -0.0,
    0.0015448330951035763, 0.028662757175574875, -0.0, -0.0,
    0.0021601506242523243, 0.03269808270414646, -0.0, -0.0,
    0.0028767195536416961, 0.036718731327978558, -0.0, -0.0,
    0.0036941716212453577, 0.040724756428522364, -0.0, -0.0,
    0.0046121399044448215, 0.044716211193075089, -0.0, -0.0,
    0.0056302588151578813, 0.048693148615486133, -0.0, -0.0,
    0.006748164094984758, 0.052655621496860691, -0.0, -0.0,
    0.0079654928103719113, 0.056603682446260746, -0.0, -0.0,
    0.0092818833477934289, 0.060537383881403589, -0.0, -0.0,
    0.010696975408949949, 0.064456778029357731, -0.0, -0.0, 0.012210410005985039,
    0.068361916927236324, -0.0, -0.0, 0.013821829456718981, 0.072252852422888064,
    -0.0, -0.0, 0.015530877379899891, 0.076129636175585519, -0.0, -0.0,
    0.017337198690472094, 0.079992319656711067, -0.0, -0.0, 0.019240439594861731,
    0.00056155944692921835, -0.0, -0.0, -0.0, 0.0011210764428963565, -0.0, -0.0,
    1.4038986173230459E-5, 0.0016785584165110745, -0.0, -0.0,
    4.2065897245639363E-5, 0.0022340127693643965, -0.0, -0.0,
    8.4029857658416242E-5, 0.0027874468761269785, -0.0, -0.0,
    0.00013988017689252618, 0.0033388680846470208, -0.0, -0.0,
    0.00020956634879570063, 0.0038882837160478245, -0.0, -0.0,
    0.00029303805091187614, 0.0044357010648249937, -0.0, -0.0,
    0.00039024514381307175, 0.0049811273989432804, -0.0, -0.0,
    0.00050113767043369654, 0.0055245699599330853, -0.0, -0.0,
    0.00062566585540727861, 0.0060660359629865959, -0.0, -0.0,
    0.00076378010440560572, 0.0066055325970535853, -0.0, -0.0,
    0.00091543100348027073, 0.0071430670249368576, -0.0, -0.0,
    0.0010805693184066103, 0.007678646383387347, -0.0, -0.0,
    0.0012591459940300318, 0.0082122777831988725, -0.0, -0.0,
    0.0014511121536147156, 0.0087439683093025449, -0.0, -0.0,
    0.0016564190981946875, 0.0092737250208608339, -0.0, -0.0,
    0.0018750183059272512, 0.0098015549513612889, -0.0, -0.0,
    0.0021068614314487719, 0.010327465108709924, -0.0, -0.0,
    0.0023519003052328044, 0.010851462475324258, -0.0, -0.0,
    0.0026100869329505525 };

  static const real_T b_Mlim_1[206]{ 612.0, 612.0, 612.0, 15.3, 15.3, 612.0,
    612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0,
    612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0,
    15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3,
    612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0,
    612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0,
    612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0,
    15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3,
    612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0,
    612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0,
    612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0,
    15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3,
    612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0,
    612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0,
    612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0,
    15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3,
    612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0,
    612.0, 612.0, 15.3, 15.3, 612.0, 612.0, 612.0, 15.3, 15.3, 612.0, 612.0,
    612.0, 15.3, 15.3, 80.0, 80.0, 80.0, -0.0, -0.0, -0.0 };

  static const real_T b_Mlim_0[166]{ 612.0, 612.0, 612.0, 15.3, 612.0, 612.0,
    612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0,
    612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3,
    612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0,
    15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0,
    612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0,
    612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3,
    612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0,
    15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0,
    612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0,
    612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3,
    612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0,
    15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0,
    612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0,
    612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3, 612.0, 612.0, 612.0, 15.3,
    80.0, 80.0, 80.0, -0.0, -0.0, -0.0 };

  static const real_T b_Kx_0[21]{ -0.69673840168599055, -0.21388804697315189,
    0.0, 0.0, -0.71544576846254548, -0.23148674875687272, 0.0,
    0.59366352289972235, -0.34839143910761772, 0.0, 0.0, 0.61058153197160281,
    -0.34227521261873245, 0.0, -0.3925750187968412, 0.42767458393192759, 0.0,
    0.0, -0.40412617467511669, 0.42693128172689154, 0.0 };

  static const real_T b_Kx_1[21]{ 0.48340300034693129, -0.59149873648563711,
    0.076274452387246269, -0.087461875059517052, 0.47981721719302906,
    -0.55183968617460621, 0.0, -0.71952158789002318, 1.0704675607877174,
    -0.11397174659709566, 0.15832856990069444, -0.7170540047020777,
    0.999301175481917, 0.0, -0.5451907758831589, -0.31069274462787111,
    -0.083754941124182636, -0.0461671145895243, -0.5263824000186843,
    -0.29298303895350225, 0.0 };

  static const real_T b_Hinv[16]{ 15.965974948812415, 14.371228479761776,
    1.9495477471050162, 0.0, 14.371228479761776, 13.826155390817547,
    -2.1000829829014531, 0.0, 1.9495477471050162, -2.1000829829014531,
    29.02217546333554, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Hinv_0[16]{ 2.0107517234261243, 4.5458833394575526,
    3.9511670626748625, 0.0, 4.5458833394575526, 14.618862953624468,
    13.733921672533578, 0.0, 3.9511670626748625, 13.733921672533578,
    13.867140525790775, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Hinv_1[16]{ 21.650975649402479, 12.603745696635549,
    2.4194444663932737, 0.0, 12.603745696635549, 7.5182342833276721,
    1.3680667678292235, 0.0, 2.4194444663932737, 1.3680667678292235,
    1.0117729338234265, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Linv[16]{ 0.65820686398155936, 3.9245100251700515,
    0.36188359365925382, 0.0, 0.0, 3.6978629787937276, -0.38982660361284555, 0.0,
    0.0, 0.0, 5.3872233537635639, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Linv_0[16]{ 0.70093703398182694, 0.62740114259684043,
    1.0610404347088109, 0.0, 0.0, 1.0084048430045782, 3.6880865806308534, 0.0,
    0.0, 0.0, 3.7238609702553043, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Linv_1[16]{ 0.70777919342616624, 3.9197482317002441,
    2.4053270178581321, 0.0, 0.0, 2.3808413457804276, 1.3600841038517122, 0.0,
    0.0, 0.0, 1.0058692429055709, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Kx[15]{ 0.79794110112653016, 0.12570901410986934,
    0.78870934466740916, 0.0, 0.0, -0.85955449427788844, -0.13541569408673249,
    -0.84960990345118026, 0.0, 0.0, -0.11660398623494962, -0.018369992635022737,
    -0.11525493979334633, 0.0, 0.0 };

  static const real_T b_Ku1[9]{ 2.2740878326680209, -2.4496825820763193,
    -0.33231488635329792, -2.4496825820763197, 2.6388359617084949,
    0.35797473482335179, -0.332314886353298, 0.35797473482335174,
    0.048561529640850407 };

  static const real_T b_Ku1_0[9]{ 2.0012420575262064, -1.2663459608291823,
    0.67424390104201049, -1.2663459608291823, 1.7371632768756062,
    -1.3934481605951545, 0.67424390104201049, -1.3934481605951548,
    1.2259405904556138 };

  static const real_T b_Ku1_1[9]{ 1.9620802154550925, -3.286488628112092,
    -0.32968191199989894, -3.286488628112092, 5.5530746742684762,
    0.30423759222659658, -0.329681911999899, 0.30423759222659663,
    1.3312340329844452 };

  static const int16_T b_Mrows_2[206]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
    90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106,
    107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
    122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
    137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151,
    152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166,
    167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
    182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
    197, 198, 199, 200, 201, 202, 203, 261, 262, 263 };

  static const int16_T b_Mlim_2[126]{ 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 80, 80, 80, 0, 0, 0 };

  static const int16_T b_Mrows_3[126]{ 1, 2, 3, 6, 7, 8, 11, 12, 13, 16, 17, 18,
    21, 22, 23, 26, 27, 28, 31, 32, 33, 36, 37, 38, 41, 42, 43, 46, 47, 48, 51,
    52, 53, 56, 57, 58, 61, 62, 63, 66, 67, 68, 71, 72, 73, 76, 77, 78, 81, 82,
    83, 86, 87, 88, 91, 92, 93, 96, 97, 98, 101, 102, 103, 106, 107, 108, 111,
    112, 113, 116, 117, 118, 121, 122, 123, 126, 127, 128, 131, 132, 133, 136,
    137, 138, 141, 142, 143, 146, 147, 148, 151, 152, 153, 156, 157, 158, 161,
    162, 163, 166, 167, 168, 171, 172, 173, 176, 177, 178, 181, 182, 183, 186,
    187, 188, 191, 192, 193, 196, 197, 198, 201, 202, 203, 261, 262, 263 };

  static const uint8_T b_Mrows_1[166]{ 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U, 10U,
    11U, 12U, 13U, 14U, 15U, 16U, 17U, 18U, 19U, 20U, 21U, 22U, 23U, 24U, 25U,
    26U, 27U, 28U, 29U, 30U, 31U, 32U, 33U, 34U, 35U, 36U, 37U, 38U, 39U, 40U,
    41U, 42U, 43U, 44U, 45U, 46U, 47U, 48U, 49U, 50U, 51U, 52U, 53U, 54U, 55U,
    56U, 57U, 58U, 59U, 60U, 61U, 62U, 63U, 64U, 65U, 66U, 67U, 68U, 69U, 70U,
    71U, 72U, 73U, 74U, 75U, 76U, 77U, 78U, 79U, 80U, 81U, 82U, 83U, 84U, 85U,
    86U, 87U, 88U, 89U, 90U, 91U, 92U, 93U, 94U, 95U, 96U, 97U, 98U, 99U, 100U,
    101U, 102U, 103U, 104U, 105U, 106U, 107U, 108U, 109U, 110U, 111U, 112U, 113U,
    114U, 115U, 116U, 117U, 118U, 119U, 120U, 121U, 122U, 123U, 124U, 125U, 126U,
    127U, 128U, 129U, 130U, 131U, 132U, 133U, 134U, 135U, 136U, 137U, 138U, 139U,
    140U, 141U, 142U, 143U, 144U, 145U, 146U, 147U, 148U, 149U, 150U, 151U, 152U,
    153U, 154U, 155U, 156U, 157U, 158U, 159U, 160U, 161U, 162U, 163U, 221U, 222U,
    223U };

  real_T Bc_0[206];
  real_T a__1_0[206];
  real_T Bc[166];
  real_T a__1[166];
  real_T Product1_j[144];
  real_T Bc_1[126];
  real_T a__1_1[126];
  real_T rseq_0[100];
  real_T rseq[80];
  real_T B_est_0[45];
  real_T B_est[36];
  real_T y_0[30];
  real_T D_est[27];
  real_T rtb_R_tmp[27];
  real_T rtb_A[25];
  real_T rtb_Product[25];
  real_T rtb_Q[25];
  real_T rtb_Z_e[25];
  real_T y[24];
  real_T Abar[16];
  real_T rtb_A_e[16];
  real_T rtb_Q_j[16];
  real_T rtb_Transpose2_0[16];
  real_T rtb_Z[16];
  real_T rtb_y_g[16];
  real_T rtb_y_m[16];
  real_T rtb_B[15];
  real_T rtb_C[15];
  real_T rtb_L[15];
  real_T rtb_N[15];
  real_T rtb_Product2[15];
  real_T Sum_h[12];
  real_T rtb_Add_k[12];
  real_T rtb_C_c[12];
  real_T rtb_N_f[12];
  real_T rtb_Product2_bg[12];
  real_T rtb_Transpose2[12];
  real_T rtb_N_h[9];
  real_T rtb_R[9];
  real_T rtb_y[9];
  real_T tmp_0[9];
  real_T rtb_xest[7];
  real_T rtb_ywt[6];
  real_T rtb_ywtT[6];
  real_T rtb_A_0[5];
  real_T rtb_TmpSignalConversionAtSFu_o4[5];
  real_T rtb_Sum2_f[4];
  real_T rtb_TmpSignalConversionAtSFu_ia[4];
  real_T Sum2_c[3];
  real_T rtb_C_0[3];
  real_T rtb_Sum1[3];
  real_T rtb_Sum6[3];
  real_T tmp[3];
  real_T umax_incr[3];
  real_T y__m[3];
  int32_T k;
  uint16_T waypt;
  int8_T b_I[16];
  int8_T P0_2_tmp[12];
  boolean_T umax_incr_flag[3];
  boolean_T rstP1;
  boolean_T rstP2;
  boolean_T rstTheta1;
  boolean_T rstTheta2;
  ZCEventType zcEvent;

  // Chart: '<Root>/SupervisoryController' incorporates:
  //   TriggerPort: '<S1>/measAvail'

  // Inport: '<Root>/measAvail'
  zcEvent = rt_ZCFcn(ANY_ZERO_CROSSING,&rtPrevZCX.SupervisoryController_Trig_ZCE,
                     (rtU.measAvail));
  if (zcEvent != NO_ZCEVENT) {
    // Gateway: SupervisoryController
    // Event: '<S1>:16'
    // During: SupervisoryController
    if (static_cast<uint32_T>(rtDW.is_active_c6_SupervisoryControl) == 0U) {
      // Entry: SupervisoryController
      rtDW.is_active_c6_SupervisoryControl = 1U;

      // Outport: '<Root>/currEv'
      // Entry Internal: SupervisoryController
      // Transition: '<S1>:2'
      // '<S1>:2:1' currEv = nullEv;
      rtY.currEv = rtP.nullEv;

      // Entry Internal 'Supervisor': '<S1>:1'
      // Entry Internal 'EventHandler': '<S1>:519'
      // Transition: '<S1>:533'
      rtDW.is_EventHandler = IN_RequestEvent;

      // Entry 'RequestEvent': '<S1>:520'
      // '<S1>:520:3' evDone = false;
      rtDW.evDone = false;

      // '<S1>:520:4' if waypt == 1
      if (rtDW.waypt == 1UL) {
        //  hold curr pos
        // '<S1>:520:5' traj(:, waypt) = y;
        for (int32_T k_0{0}; k_0 < 6; k_0++) {
          // Inport: '<Root>/y'
          rtDW.traj[k_0] = rtU.y[k_0];
        }
      } else {
        // '<S1>:520:6' else
        //  hold last waypoint pos
        // '<S1>:520:7' traj(:,1) = traj(:, waypt);
        for (int32_T k_0{0}; k_0 < 6; k_0++) {
          rtDW.traj[k_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6
            + k_0];
        }

        // '<S1>:520:8' waypt = 1;
        rtDW.waypt = 1U;
      }

      // Outport: '<Root>/requestEvent'
      // '<S1>:520:10' requestEvent = true;
      rtY.requestEvent = true;

      //  request new event
      // Entry 'ControlLaw': '<S1>:59'
    } else {
      __m128d tmp_3;
      real_T dwt;
      real_T umin_scale1_idx_0;
      real_T umin_scale1_idx_1;
      real_T umin_scale1_idx_2;
      int32_T i;
      int32_T ii;
      int32_T k_0;
      boolean_T c_y;
      boolean_T exitg1;
      boolean_T guard1{ false };

      // During 'Supervisor': '<S1>:1'
      // During 'EventHandler': '<S1>:519'
      if (static_cast<uint32_T>(rtDW.is_EventHandler) == IN_HandleEvent) {
        // During 'HandleEvent': '<S1>:521'
        // '<S1>:534:1' sf_internal_predicateOutput = evDone;
        if (rtDW.evDone) {
          // Transition: '<S1>:534'
          rtDW.is_EventHandler = IN_RequestEvent;

          // Entry 'RequestEvent': '<S1>:520'
          // '<S1>:520:3' evDone = false;
          rtDW.evDone = false;

          // '<S1>:520:4' if waypt == 1
          if (rtDW.waypt == 1UL) {
            //  hold curr pos
            // '<S1>:520:5' traj(:, waypt) = y;
            for (k_0 = 0; k_0 < 6; k_0++) {
              // Inport: '<Root>/y'
              rtDW.traj[k_0] = rtU.y[k_0];
            }
          } else {
            // '<S1>:520:6' else
            //  hold last waypoint pos
            // '<S1>:520:7' traj(:,1) = traj(:, waypt);
            for (k_0 = 0; k_0 < 6; k_0++) {
              rtDW.traj[k_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                6 + k_0];
            }

            // '<S1>:520:8' waypt = 1;
            rtDW.waypt = 1U;
          }

          // Outport: '<Root>/requestEvent'
          // '<S1>:520:10' requestEvent = true;
          rtY.requestEvent = true;

          //  request new event
        } else {
          // Outport: '<Root>/currEv'
          // '<S1>:521:7' [evDone, waypt, currEv.postT] = handleEvent(currEv, currEv.postT); 
          handleEvent(&rtY.currEv.postT, &rtDW.evDone, &waypt);
          rtDW.waypt = waypt;
        }
      } else {
        // During 'RequestEvent': '<S1>:520'
        // '<S1>:522:1' sf_internal_predicateOutput = ~isequal(nextEv, nullEv);
        rstP2 = false;

        // Inport: '<Root>/nextEv'
        if (rtU.nextEv.postT == rtP.nullEv.postT) {
          if (rtU.nextEv.moveT == rtP.nullEv.moveT) {
            if (rtU.nextEv.preT == rtP.nullEv.preT) {
              rstTheta2 = true;
              k = 0;
              exitg1 = false;
              while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                       (0U)) == false) && (k < 6)) {
                if (!(rtU.nextEv.r[k] == rtP.nullEv.r[k])) {
                  rstTheta2 = false;
                  exitg1 = true;
                } else {
                  k++;
                }
              }
            } else {
              rstTheta2 = false;
            }
          } else {
            rstTheta2 = false;
          }
        } else {
          rstTheta2 = false;
        }

        if (rstTheta2) {
          rstP2 = true;
        }

        if (!rstP2) {
          // Transition: '<S1>:522'
          // '<S1>:522:1' evDone = false;
          rtDW.evDone = false;

          // Outport: '<Root>/requestEvent'
          // Exit 'RequestEvent': '<S1>:520'
          // '<S1>:520:12' requestEvent = false;
          rtY.requestEvent = false;
          rtDW.is_EventHandler = IN_HandleEvent;

          // Outport: '<Root>/currEv' incorporates:
          //   Inport: '<Root>/nextEv'
          //   Inport: '<Root>/y'

          // Entry 'HandleEvent': '<S1>:521'
          // '<S1>:521:3' currEv = nextEv;
          rtY.currEv = rtU.nextEv;

          // '<S1>:521:4' [traj, trajSize] = trajGen(currEv, y);
          trajGen(&rtY.currEv, rtU.y, rtDW.traj, &rtDW.trajSize);

          // '<S1>:521:5' disp("received event");
          (void)printf("%s\n", "received event");
          (void)fflush(stdout);
        }
      }

      // Outport: '<Root>/currEv' incorporates:
      //   Outport: '<Root>/P'
      //   Outport: '<Root>/theta'

      // During 'ControlLaw': '<S1>:59'
      // '<S1>:59:3' if any(currEv.r(1:no)) && any(zeroCross(no+1:2*no))
      guard1 = false;
      if (any(&rtY.currEv.r[0])) {
        rstP2 = false;
        k = 3;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (k - 3 < 3)) {
          if (rtU.yo[k]) {
            rstP2 = true;
            exitg1 = true;
          } else {
            k++;
          }
        }

        if (rstP2) {
          // '<S1>:59:4' rstP2 = true;
          // '<S1>:59:4' rstTheta2 = true;
          rstTheta2 = true;

          // '<S1>:59:4' enAdapt_(no+1:2*no) = false;
          rtDW.enAdapt_[3] = false;
          rtDW.enAdapt_[4] = false;
          rtDW.enAdapt_[5] = false;

          // '<S1>:59:5' P0_2 = P(1:np*no, 1:np*no);
          for (k_0 = 0; k_0 < 12; k_0++) {
            P0_2_tmp[k_0] = static_cast<int8_T>(k_0 + 1);
          }

          // '<S1>:59:5' theta0_2 = theta(1:np*no);
          for (k = 0; k < 12; k++) {
            for (k_0 = 0; k_0 < 12; k_0++) {
              rtDW.P0_2[k_0 + 12 * k] = rtY.P_c[((P0_2_tmp[k] - 1) * 24 +
                P0_2_tmp[k_0]) - 1];
            }

            rtDW.theta0_2[k] = rtY.theta[k];
          }
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        // '<S1>:59:6' else
        // '<S1>:59:7' rstP2 = false;
        rstP2 = false;

        // '<S1>:59:7' rstTheta2 = false;
        rstTheta2 = false;

        // '<S1>:59:7' enAdapt_(no+1:2*no) = true;
        rtDW.enAdapt_[3] = true;
        rtDW.enAdapt_[4] = true;
        rtDW.enAdapt_[5] = true;
      }

      // '<S1>:59:9' if any(currEv.r(no+1:2*no)) && any(zeroCross(1:no))
      guard1 = false;
      if (any(&rtY.currEv.r[3])) {
        rstP1 = false;
        k = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (k < 3)) {
          if (rtU.yo[k]) {
            rstP1 = true;
            exitg1 = true;
          } else {
            k++;
          }
        }

        if (rstP1) {
          // '<S1>:59:10' rstP1 = true;
          // '<S1>:59:10' rstTheta1 = true;
          rstTheta1 = true;

          // '<S1>:59:10' enAdapt_(1:no) = false;
          rtDW.enAdapt_[0] = false;
          rtDW.enAdapt_[1] = false;
          rtDW.enAdapt_[2] = false;

          // '<S1>:59:11' P0_1 = P(np*no+1:2*np*no, np*no+1:2*np*no);
          for (k_0 = 0; k_0 < 12; k_0++) {
            P0_2_tmp[k_0] = static_cast<int8_T>(k_0 + 13);
          }

          // '<S1>:59:11' theta0_1 = theta(np*no+1:2*np*no);
          for (k = 0; k < 12; k++) {
            for (k_0 = 0; k_0 < 12; k_0++) {
              rtDW.P0_1[k_0 + 12 * k] = rtY.P_c[((P0_2_tmp[k] - 1) * 24 +
                P0_2_tmp[k_0]) - 1];
            }

            rtDW.theta0_1[k] = rtY.theta[k + 12];
          }
        } else {
          guard1 = true;
        }
      } else {
        guard1 = true;
      }

      if (guard1) {
        // '<S1>:59:12' else
        // '<S1>:59:13' rstP1 = false;
        rstP1 = false;

        // '<S1>:59:13' rstTheta1 = false;
        rstTheta1 = false;

        // '<S1>:59:13' enAdapt_(1:no) = true;
        rtDW.enAdapt_[0] = true;
        rtDW.enAdapt_[1] = true;
        rtDW.enAdapt_[2] = true;
      }

      // '<S1>:59:15' if ~all(enAdapt)
      c_y = true;
      k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k < 6)) {
        if (!rtU.enAdapt[k]) {
          c_y = false;
          exitg1 = true;
        } else {
          k++;
        }
      }

      if (!c_y) {
        // '<S1>:59:16' enAdapt_(:) = false;
        for (k = 0; k < 6; k++) {
          rtDW.enAdapt_[k] = false;
        }
      }

      // Outputs for Function Call SubSystem: '<S1>/paramEst1'
      // Inport: '<Root>/y' incorporates:
      //   Inport: '<Root>/lambda'
      //   Inport: '<Root>/p_'
      //   Inport: '<Root>/u0'
      //   Inport: '<Root>/y0'
      //   Outport: '<Root>/u'

      // '<S1>:59:18' [theta(1:np*no), P(1:np*no, 1:np*no), prmErr(1:no)] = ...
      // '<S1>:59:19'     paramEst1(y(1:no), y0(1:no), u, u0,...
      // '<S1>:59:20'     enAdapt_(1:no),...
      // '<S1>:59:21'     theta0_1, thetaSgn(1:np*no), rstTheta1, P0_1, rstP1,... 
      // '<S1>:59:22'     p_, dPmod_, lambda);
      // Simulink Function 'paramEst1': '<S1>:574'
      paramEst1(&rtU.y[0], &rtU.y0[0], rtY.u, rtU.u0, &rtDW.enAdapt_[0],
                rtDW.theta0_1, &rtDW.thetaSgn[0], rstTheta1, rtDW.P0_1, rstP1,
                rtU.p_, rtU.lambda, Sum_h, Product1_j, Sum2_c, &rtDW.paramEst1_o,
                &rtPrevZCX.paramEst1_o);

      // End of Outputs for SubSystem: '<S1>/paramEst1'
      k = 0;
      ii = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S300>/Product1'

        (void)std::memcpy(&rtY.P_c[k], &Product1_j[ii], 12U * sizeof(real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i] = Sum_h[i];
        k += 24;
        ii += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S300>/Sum2'

      rtY.prmErr[0] = Sum2_c[0];
      rtY.prmErr[1] = Sum2_c[1];
      rtY.prmErr[2] = Sum2_c[2];

      // Outputs for Function Call SubSystem: '<S1>/paramEst2'
      // Inport: '<Root>/y' incorporates:
      //   Inport: '<Root>/lambda'
      //   Inport: '<Root>/p_'
      //   Inport: '<Root>/u0'
      //   Inport: '<Root>/y0'
      //   Outport: '<Root>/u'

      // '<S1>:59:23' [theta(np*no+1:2*np*no), P(np*no+1:2*np*no, np*no+1:2*np*no), prmErr(no+1:2*no)] = ... 
      // '<S1>:59:24'     paramEst2(y(no+1:2*no), y0(no+1:2*no), u, u0,...
      // '<S1>:59:25'     enAdapt_(no+1:2*no),...
      // '<S1>:59:26'     theta0_2, thetaSgn(np*no+1:2*np*no), rstTheta2, P0_2, rstP2,... 
      // '<S1>:59:27'     p_, dPmod_, lambda);
      // Simulink Function 'paramEst2': '<S1>:844'
      paramEst1(&rtU.y[3], &rtU.y0[3], rtY.u, rtU.u0, &rtDW.enAdapt_[3],
                rtDW.theta0_2, &rtDW.thetaSgn[12], rstTheta2, rtDW.P0_2, rstP2,
                rtU.p_, rtU.lambda, Sum_h, Product1_j, Sum2_c, &rtDW.paramEst2,
                &rtPrevZCX.paramEst2);

      // End of Outputs for SubSystem: '<S1>/paramEst2'
      k = 0;
      ii = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S304>/Product1'

        (void)std::memcpy(&rtY.P_c[k + 300], &Product1_j[ii], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i + 12] = Sum_h[i];
        k += 24;
        ii += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S304>/Sum2'

      rtY.prmErr[3] = Sum2_c[0];
      rtY.prmErr[4] = Sum2_c[1];
      rtY.prmErr[5] = Sum2_c[2];

      // '<S1>:59:28' currTraj = traj(:,waypt);
      for (k_0 = 0; k_0 < 6; k_0++) {
        // Outport: '<Root>/currTraj'
        rtY.currTraj[k_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6
          + k_0];
      }

      // Outport: '<Root>/sig' incorporates:
      //   Outport: '<Root>/ywt'

      // '<S1>:59:29' sig = gainSchSig(ywt);
      rtY.sig = gainSchSig(rtY.ywt);

      // Outputs for Function Call SubSystem: '<S1>/wtMod'
      // MATLAB Function: '<S8>/MATLAB Function' incorporates:
      //   Inport: '<Root>/k_2'

      // '<S1>:59:30' [ywt, y__, r_] = wtMod(currEv.r, y, k_2, currTraj);
      // Simulink Function 'wtMod': '<S1>:911'
      // MATLAB Function 'SupervisoryController/wtMod/MATLAB Function': '<S306>:1' 
      // '<S306>:1:2' [ywt, ywtT, uwt, uwtT] = wtMod_(y, yDest, ywtT, uwtT, dt, no, ni, k_2); 
      // 'wtMod_:3' ywt = zeros(1,2*no);
      // 'wtMod_:4' uwt = dt*ones(1,ni);
      //  time-scaled sigmoid (around 0->1 in 0->k_2 seconds, k_2 being a time scale factor) 
      // 'wtMod_:7' r = 0.2;
      //  10% to 90% rise time
      // 'wtMod_:8' dwt = dt/k_2;
      dwt = rtP.dt / rtU.k_2;

      // 'wtMod_:9' k_1 = 2.197/(r*k_2);
      umin_scale1_idx_0 = 2.197 / (0.2 * rtU.k_2);

      // 'wtMod_:10' x0 = 0.5*k_2;
      umin_scale1_idx_1 = 0.5 * rtU.k_2;

      //  midpoint
      // 'wtMod_:11' sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));
      // 'wtMod_:13' for i = 1:2*no
      for (k = 0; k < 6; k++) {
        // Delay: '<S8>/Delay'
        umin_scale1_idx_2 = rtDW.Delay_DSTATE[k];

        // MATLAB Function: '<S8>/MATLAB Function' incorporates:
        //   Inport: '<Root>/k_2'
        //   Outport: '<Root>/currEv'

        // 'wtMod_:14' if yDest(i) ~= 0
        if (rtY.currEv.r[k] != 0.0) {
          //  drive ywt to 1
          // 'wtMod_:16' if (ywtT(i) <= 1)
          if (umin_scale1_idx_2 <= 1.0) {
            // 'wtMod_:17' ywtT(i) = ywtT(i) + dwt;
            umin_scale1_idx_2 += dwt;
          }

          // 'wtMod_:19' else
          //  drive ywt to 0
          // 'wtMod_:21' if (ywtT(i) > 0)
        } else if (umin_scale1_idx_2 > 0.0) {
          // 'wtMod_:22' ywtT(i) = ywtT(i) - dwt;
          umin_scale1_idx_2 -= dwt;
        } else {
          // no actions
        }

        // 'wtMod_:25' if ywtT(i) <= 0
        if (umin_scale1_idx_2 <= 0.0) {
          // 'wtMod_:26' ywt(i) = 0;
          rtb_ywt[k] = 0.0;
        } else {
          // 'wtMod_:27' else
          // 'wtMod_:28' ywt(i) = sigmoid(ywtT(i)*k_2);
          // 'wtMod_:11' @(x) 1/(1 + exp(-k_1*(x-x0)))
          rtb_ywt[k] = 1.0 / (std::exp((umin_scale1_idx_2 * rtU.k_2 -
            umin_scale1_idx_1) * -umin_scale1_idx_0) + 1.0);
        }

        // Delay: '<S8>/Delay'
        rtb_ywtT[k] = umin_scale1_idx_2;
      }

      // MATLAB Function: '<S8>/MATLAB Function' incorporates:
      //   Inport: '<Root>/y'
      //   Outport: '<Root>/currEv'
      //   Outport: '<Root>/currTraj'

      //  for i = 1:ni
      //      if yDest(i) ~= 0
      //          % drive uwt to 0
      //          if (uwtT(i) > 0)
      //              uwtT(i) = uwtT(i) - dwt;
      //          end
      //      else
      //          % drive uwt to 1
      //          if (uwtT(i) <= 1)
      //              uwtT(i) = uwtT(i) + dwt;
      //          end
      //      end
      //      uwt(i) = sigmoid(uwtT(i)*k_2);
      //  end
      // '<S306>:1:3' r_ = zeros(no, 1);
      // '<S306>:1:4' if any(yDest(1:no))
      if (any(&rtY.currEv.r[0])) {
        // '<S306>:1:5' r_ = r(1:no);
        // '<S306>:1:6' y_ = y(1:no);
        Sum2_c[0] = rtY.currTraj[0];
        y__m[0] = rtU.y[0];
        Sum2_c[1] = rtY.currTraj[1];
        y__m[1] = rtU.y[1];
        Sum2_c[2] = rtY.currTraj[2];
        y__m[2] = rtU.y[2];
      } else if (any(&rtY.currEv.r[3])) {
        // '<S306>:1:7' elseif any(yDest(no+1:2*no))
        // '<S306>:1:8' r_ = r(no+1:2*no);
        // '<S306>:1:9' y_ = y(no+1:2*no);
        Sum2_c[0] = rtY.currTraj[3];
        y__m[0] = rtU.y[3];
        Sum2_c[1] = rtY.currTraj[4];
        y__m[1] = rtU.y[4];
        Sum2_c[2] = rtY.currTraj[5];
        y__m[2] = rtU.y[5];
      } else {
        // '<S306>:1:10' else
        // '<S306>:1:11' r_ = zeros(no, 1);
        // '<S306>:1:12' y_ = y(1:no);
        Sum2_c[0] = 0.0;
        y__m[0] = rtU.y[0];
        Sum2_c[1] = 0.0;
        y__m[1] = rtU.y[1];
        Sum2_c[2] = 0.0;
        y__m[2] = rtU.y[2];
      }

      // End of Outputs for SubSystem: '<S1>/wtMod'
      for (k = 0; k <= 4; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/wtMod'
        // Update for Delay: '<S8>/Delay'
        tmp_3 = _mm_loadu_pd(&rtb_ywtT[k]);
        (void)_mm_storeu_pd(&rtDW.Delay_DSTATE[k], tmp_3);

        // Gain: '<S8>/Gain' incorporates:
        //   Delay: '<S8>/Delay'

        tmp_3 = _mm_loadu_pd(&rtb_ywt[k]);

        // Outport: '<Root>/ywt' incorporates:
        //   Delay: '<S8>/Delay'
        //   Gain: '<S8>/Gain'

        (void)_mm_storeu_pd(&rtY.ywt[k], _mm_mul_pd(_mm_set1_pd(rtP.beta), tmp_3));

        // End of Outputs for SubSystem: '<S1>/wtMod'
      }

      // Outport: '<Root>/sig' incorporates:
      //   Constant: '<S182>/G'
      //   Constant: '<S182>/H'
      //   Constant: '<S252>/G'
      //   Constant: '<S252>/H'
      //   Constant: '<S4>/Constant1'
      //   Constant: '<S4>/Constant13'
      //   Constant: '<S5>/Constant1'
      //   Constant: '<S5>/Constant13'
      //   DataTypeConversion: '<S182>/DataTypeConversionEnable'
      //   DataTypeConversion: '<S252>/DataTypeConversionEnable'
      //   Delay: '<S182>/MemoryP'
      //   Delay: '<S182>/MemoryX'
      //   Delay: '<S252>/MemoryP'
      //   Delay: '<S252>/MemoryX'
      //   Outport: '<Root>/yhat'
      //   Product: '<S185>/Product'
      //   Product: '<S185>/Product1'
      //   Product: '<S255>/Product'
      //   Product: '<S255>/Product1'
      //   Sum: '<S159>/Sum3'
      //   Sum: '<S229>/Sum3'
      //   Sum: '<S89>/Sum3'

      // [u, yhat] = ampc(currTraj, y, ymax, ywt, y0, x0, u0, umax, uwt, excitation, theta, thetaSgn); 
      // [u, ywt, currTraj] = gmpc(traj(:,waypt), currEv.r, y, ymax, umax, uwt, k_2); 
      // '<S1>:59:33' if sig == 1
      if (rtY.sig == 1.0) {
        __m128d tmp_1;
        __m128d tmp_2;
        real_T Bc_2;
        int32_T rtb_Q_tmp;
        int32_T rtb_y_tmp;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'

        // '<S1>:59:34' [u, yhat(1:no)] = mpc1(r_, y__, [0;0;0], 0, u0, umax, iRST); 
        // Simulink Function 'mpc1': '<S1>:882'
        if (rtU.iRST && (rtDW.DiscreteTimeIntegrator_PrevRe_b <= 0)) {
          rtDW.DiscreteTimeIntegrator_DSTATE_j = rtP.DiscreteTimeIntegrator_IC;
        }

        // Delay: '<S112>/MemoryX' incorporates:
        //   Constant: '<S112>/X0'
        //   DataTypeConversion: '<S112>/DataTypeConversionReset'

        rtDW.icLoad_n = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE_l) ==
                          POS_ZCSIG) || rtDW.icLoad_n);
        rtPrevZCX.MemoryX_Reset_ZCE_l = 0U;
        if (rtDW.icLoad_n) {
          rtDW.MemoryX_DSTATE_l[0] = rtP.X0_Value_f[0];
          rtDW.MemoryX_DSTATE_l[1] = rtP.X0_Value_f[1];
          rtDW.MemoryX_DSTATE_l[2] = rtP.X0_Value_f[2];
          rtDW.MemoryX_DSTATE_l[3] = rtP.X0_Value_f[3];
        }

        // SignalConversion generated from: '<S111>/ SFunction ' incorporates:
        //   MATLAB Function: '<S110>/optimizer'

        rtb_TmpSignalConversionAtSFu_ia[0] = Sum2_c[0];
        rtb_TmpSignalConversionAtSFu_ia[1] = Sum2_c[1];
        rtb_TmpSignalConversionAtSFu_ia[2] = Sum2_c[2];

        // MATLAB Function: '<S110>/optimizer' incorporates:
        //   Constant: '<S3>/Constant'
        //   Constant: '<S89>/Constant1'
        //   Delay: '<S112>/MemoryX'
        //   DiscreteIntegrator: '<S3>/Discrete-Time Integrator'
        //   Gain: '<S3>/Gain2'
        //   Inport: '<Root>/umax'
        //   SignalConversion generated from: '<S111>/ SFunction '
        //   Sum: '<S89>/Sum2'
        //   UnitDelay: '<S90>/last_mv'

        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S111>:1'
        // '<S111>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S111>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S111>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S111>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S111>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S111>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S111>:1:252' if isa(u,'double')
        // '<S111>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S111>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S111>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S111>:1:252' if isa(u,'double')
        // '<S111>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S111>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S111>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S111>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S111>:1:252' if isa(u,'double')
        // '<S111>:1:253' y = u;
        //    md:         measured disturbance
        // '<S111>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S111>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S111>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S111>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S111>:1:252' if isa(u,'double')
        // '<S111>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S111>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S111>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S111>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S111>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S111>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S111>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S111>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S111>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S111>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S111>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S111>:1:58' uwt = convertDataType(uwt0,isDouble);
        //    duwt:       run-time DMV weights
        // '<S111>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S111>:1:62' ewt = convertDataType(ewt0,isDouble);
        //    iA:         logical previous active set (for warm start)
        //  Outputs (in BlockDataType except iAout)
        //    xk1:        x[k+1|k] from built-in KF
        //    u:          optimal MV
        //    cost:       optimal cost
        //    useq:       optimal MV sequence
        //    xseq:       optimal state sequence
        //    yseq:       optimal OV sequence
        //    status:     QP exit flag
        //    xest:       x[k|k] from built-in KF
        //    iAout:      logical current active set
        //  Parameters (constant)
        //    dimensions (int32):
        //        nx, nxp, nup, nu, ny, degrees, p, nxQP, enable_value, Mrows, nCC, nv 
        //        myindex, mvindex, mdindex, unindex, nxid, m, Ndis, numdis, maxdis 
        //    MPC constants (BlockDataType):
        //        Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, Mx, Mu1, Mv, utarget
        //        H, Linv, Ac, Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1
        //        A, C, B, D, Cid, Did, Ecc, Fcc, Scc, Gcc
        //        RYscale, RMDscale, xoff, Uscale, Yscale
        //        uoff, voff, yoff, myoff, RMVscale, Mdis, Vdis
        //    configurations (logical):
        //        isQP, CustomSolver, CustomSolverCodeGen, UseSuboptimalSolution, UseActiveSetSolver 
        //        openloopflag, no_umin, no_umax, no_ymin, no_ymax, switch_inport, no_switch 
        //        return_cost, return_mvseq, return_xseq, return_ovseq, isLTV
        //        no_ywt, no_uwt, no_duwt, no_rhoeps, no_md, no_ref, no_uref, no_mv 
        //        CustomEstimation, no_cc, isHyb, isDouble
        //    ASOptions
        //    IPOptions
        //    MIQPOptions
        //  Parameters
        // '<S111>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S111>:1:96' isAdaptive = false;
        // '<S111>:1:97' isLTV = false;
        // '<S111>:1:98' ZERO = zeros('like',ref);
        // '<S111>:1:99' ONE = ones('like',ref);
        // '<S111>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S111>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S111>:1:119' if isSimulation
        // '<S111>:1:126' else
        //  When doing code generation, use M code directly
        // '<S111>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        k_0 = 0;
        for (k = 0; k < 20; k++) {
          rseq[k_0] = rtb_TmpSignalConversionAtSFu_ia[0];
          rseq[k_0 + 1] = rtb_TmpSignalConversionAtSFu_ia[1];
          rseq[k_0 + 2] = rtb_TmpSignalConversionAtSFu_ia[2];
          rseq[k_0 + 3] = rtP.Constant_Value;
          k_0 += 4;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S111>:1:133' old_u = old_u - uoff;
        // '<S111>:1:134' if no_mv
        // '<S111>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S111>:1:143' xk = xk - xoff;
        rtb_TmpSignalConversionAtSFu_o4[0] = rtDW.MemoryX_DSTATE_l[0];
        rtb_TmpSignalConversionAtSFu_o4[1] = rtP.dt *
          rtDW.DiscreteTimeIntegrator_DSTATE_j;
        umin_scale1_idx_0 = rtDW.last_mv_DSTATE_n[0];
        rtb_TmpSignalConversionAtSFu_o4[2] = rtP.Constant1_Value_j[0] +
          rtDW.MemoryX_DSTATE_l[1];
        umin_scale1_idx_1 = rtDW.last_mv_DSTATE_n[1];
        rtb_TmpSignalConversionAtSFu_o4[3] = rtP.Constant1_Value_j[1] +
          rtDW.MemoryX_DSTATE_l[2];
        umin_scale1_idx_2 = rtDW.last_mv_DSTATE_n[2];
        rtb_TmpSignalConversionAtSFu_o4[4] = rtP.Constant1_Value_j[2] +
          rtDW.MemoryX_DSTATE_l[3];

        //  Remove offset
        // '<S111>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S111>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S111>:1:162' if no_uref
        //  no external utarget
        // '<S111>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S111>:1:173' if ~no_cc
        // '<S111>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S111>:1:183' if isSimulation
        // '<S111>:1:214' else
        //  When doing code generation, use M code directly
        // '<S111>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S111>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S111>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S111>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S111>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S111>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S111>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S111>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S111>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S111>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S111>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S111>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S111>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        umax_incr_flag[0] = false;
        umax_incr[0] = 0.0;
        umax_incr_flag[1] = false;
        umax_incr[1] = 0.0;
        umax_incr_flag[2] = false;
        umax_incr[2] = 0.0;
        for (k = 0; k < 166; k++) {
          uint8_T b_Mrows;
          dwt = b_Mlim_0[k];
          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 5; k_0++) {
            Bc_2 += b_a[166 * k_0 + k] * rtb_TmpSignalConversionAtSFu_o4[k_0];
          }

          Bc_2 = -(((a[k + 166] * umin_scale1_idx_1 + a[k] * umin_scale1_idx_0)
                    + a[k + 332] * umin_scale1_idx_2) + (dwt + Bc_2));
          b_Mrows = b_Mrows_1[k];
          if ((b_Mrows > 80UL) && (b_Mrows > 160UL) && (b_Mrows <= 220UL)) {
            ii = (static_cast<int32_T>(b_Mrows) - div_nde_s32_floor(static_cast<
                   int32_T>(b_Mrows) - 161, static_cast<int32_T>(nu)) *
                  static_cast<int32_T>(nu)) - 161;
            rstP2 = umax_incr_flag[ii];
            if (!umax_incr_flag[ii]) {
              dwt = -rtU.umax[ii] - (-dwt);
              rstP2 = true;
            } else {
              dwt = umax_incr[ii];
            }

            umax_incr[ii] = dwt;
            umax_incr_flag[ii] = rstP2;
            Bc_2 += dwt;
          }

          Bc[k] = Bc_2;
        }

        rtb_Sum2_f[0] = 0.0;
        rtb_Sum2_f[1] = 0.0;
        rtb_Sum2_f[2] = 0.0;
        rtb_Sum2_f[3] = 0.0;
        for (k = 0; k < 3; k++) {
          dwt = 0.0;
          for (k_0 = 0; k_0 < 5; k_0++) {
            dwt += b_Kx[5 * k + k_0] * rtb_TmpSignalConversionAtSFu_o4[k_0];
          }

          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 80; k_0++) {
            Bc_2 += b_Kr[80 * k + k_0] * rseq[k_0];
          }

          rtb_Sum2_f[k] = ((b_Ku1[3 * k + 1] * umin_scale1_idx_1 + b_Ku1[3 * k] *
                            umin_scale1_idx_0) + b_Ku1[3 * k + 2] *
                           umin_scale1_idx_2) + (dwt + Bc_2);
        }

        // Update for Memory: '<S90>/Memory' incorporates:
        //   MATLAB Function: '<S110>/optimizer'

        qpkwik(b_Linv, b_Hinv, rtb_Sum2_f, b_Ac, Bc, rtDW.Memory_PreviousInput_d,
               680, 1.0E-6, rtb_TmpSignalConversionAtSFu_ia, a__1, &k);

        // MATLAB Function: '<S110>/optimizer' incorporates:
        //   UnitDelay: '<S90>/last_mv'

        if ((k < 0) || (k == 0)) {
          rtb_TmpSignalConversionAtSFu_ia[0] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[1] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[2] = 0.0;
        }

        umax_incr[0] = rtDW.last_mv_DSTATE_n[0] +
          rtb_TmpSignalConversionAtSFu_ia[0];
        umax_incr[1] = rtDW.last_mv_DSTATE_n[1] +
          rtb_TmpSignalConversionAtSFu_ia[1];
        umax_incr[2] = rtDW.last_mv_DSTATE_n[2] +
          rtb_TmpSignalConversionAtSFu_ia[2];

        // Delay: '<S112>/MemoryP' incorporates:
        //   Constant: '<S112>/P0'
        //   DataTypeConversion: '<S112>/DataTypeConversionReset'

        // '<S111>:1:231' if return_xseq || return_ovseq
        // '<S111>:1:233' else
        // '<S111>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S111>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S111>:1:238' if CustomEstimation
        // '<S111>:1:239' xk1 = xk;
        // '<S111>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S111>:1:247' xest = xest + xoff;
        rtDW.icLoad_h = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE_n) ==
                          POS_ZCSIG) || rtDW.icLoad_h);
        rtPrevZCX.MemoryP_Reset_ZCE_n = 0U;
        if (rtDW.icLoad_h) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE_e[0], &rtP.P0_Value_a[0],
                            sizeof(real_T) << 4UL);
        }

        // MATLAB Function: '<S89>/MATLAB Function' incorporates:
        //   BusCreator: '<S3>/Bus Creator1'
        //   Constant: '<S3>/Constant12'
        //   Constant: '<S3>/Constant13'
        //   Constant: '<S3>/Constant3'
        //   Constant: '<S3>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc1/State Estimator OD (KF)/MATLAB Function': '<S113>:1' 
        // '<S113>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod1, Bod1, Cod1, Dod1, Dmn1, 3, 3, 1); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&Sum_h[0], 0, 12U * sizeof(real_T));

        //  n_states  x n_inputs
        // 'stateEst_:9' C = zeros(no,ns);
        //  n_outputs x n_states
        // 'stateEst_:10' D = zeros(no,ni);
        //  n_outputs x n_inputs
        // 'stateEst_:11' Q = zeros(ns,ns);
        //  n_states  x n_states
        // 'stateEst_:12' G = eye(ns);
        //  n_states  x n_states
        // 'stateEst_:13' R = zeros(no,no);
        //  n_outputs x n_outputs
        // 'stateEst_:14' N = zeros(ns,no);
        //  n_states  x n_outputs
        // 'stateEst_:15' H = zeros(no,ns);
        //  n_outputs x n_states
        //  combine plant and output disturbance model
        //  (force the outputs to fit in preallocated memory)
        // 'stateEst_:19' A(1:ns, 1:ns) = blkdiag(Ap, Aod);
        (void)std::memset(&rtb_A_e[0], 0, sizeof(real_T) << 4UL);
        rtb_A_e[0] = rtP.Constant3_Value;

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        k = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          rtb_A_e[k + 5] = rtP.Aod1[ii];
          rtb_A_e[k + 6] = rtP.Aod1[ii + 1];
          rtb_A_e[k + 7] = rtP.Aod1[ii + 2];
          Sum_h[k] = rtP.Constant4_Value[i];
          rtb_C_c[i] = rtP.Constant12_Value_e[i];
          k += 4;
          ii += 3;
        }

        (void)std::memcpy(&rtb_C_c[3], &rtP.Cod1[0], 9U * sizeof(real_T));

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est[0], 0, 36U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y[0], 0, 24U * sizeof(real_T));
        k_0 = 0;
        k = 0;
        for (ii = 0; ii < 3; ii++) {
          y[k_0] = rtP.Constant4_Value[ii];
          y[k_0 + 13] = rtP.Bod1[k];
          y[k_0 + 14] = rtP.Bod1[k + 1];
          y[k_0 + 15] = rtP.Bod1[k + 2];
          k_0 += 4;
          k += 3;
        }

        k_0 = 0;
        for (k = 0; k < 6; k++) {
          B_est[k_0] = y[k_0];
          B_est[k_0 + 1] = y[k_0 + 1];
          B_est[k_0 + 2] = y[k_0 + 2];
          B_est[k_0 + 3] = y[k_0 + 3];
          k_0 += 4;
        }

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        for (k_0 = 0; k_0 < 9; k_0++) {
          D_est[k_0] = rtP.Constant13_Value_c[k_0];
          D_est[k_0 + 9] = rtP.Dod1[k_0];
          D_est[k_0 + 18] = rtP.Dmn1[k_0];
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (k_0 = 0; k_0 < 4; k_0++) {
          k = 0;
          for (ii = 0; ii < 4; ii++) {
            rtb_Q_tmp = k + k_0;
            rtb_Q_j[rtb_Q_tmp] = 0.0;
            i = 0;
            for (rtb_y_tmp = 0; rtb_y_tmp < 9; rtb_y_tmp++) {
              rtb_Q_j[rtb_Q_tmp] += B_est[i + k_0] * B_est[i + ii];
              i += 4;
            }

            k += 4;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        k_0 = 0;
        for (k = 0; k < 9; k++) {
          rtb_R_tmp[k] = D_est[k_0];
          rtb_R_tmp[k + 9] = D_est[k_0 + 1];
          rtb_R_tmp[k + 18] = D_est[k_0 + 2];
          k_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (k_0 = 0; k_0 < 3; k_0++) {
          for (k = 0; k < 3; k++) {
            i = 3 * k_0 + k;
            rtb_R[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_R[i] += D_est[3 * ii + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }

          for (k = 0; k < 4; k++) {
            i = (k_0 << 2UL) + k;
            rtb_N_f[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_N_f[i] += B_est[(ii << 2UL) + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }
        }

        // End of MATLAB Function: '<S89>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S112>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S112>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S112>/ScalarExpansionQ'
        // MATLAB Function: '<S134>/ScalarExpansion'
        //  ctrlScalarExpansion Helper function for scalar expansion.
        //
        //  	y  = ctrlScalarExpansion(u,n)
        //
        //    An n-ny-n matrix y is created based on u. If u is a scalar, y has u  
        //    on its diagonals. If u is a vector, y has the elements of u on its 
        //    diagonals. If u is a matrix y = (u+u.')/2.
        //
        //    When u is scalar or vector, we enforce symmetric positive-definiteness. 
        //    When u is a matrix, we enly enforce symmetry.
        // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S156>:1' 
        //    Copyright 2014-2015 The MathWorks, Inc.
        // '<S156>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
        k_0 = 0;
        for (k = 0; k < 4; k++) {
          rtb_y_g[k_0] = (rtb_Q_j[k_0] + rtb_Q_j[k]) / 2.0;
          rtb_y_g[k_0 + 1] = (rtb_Q_j[k_0 + 1] + rtb_Q_j[k + 4]) / 2.0;
          rtb_y_g[k_0 + 2] = (rtb_Q_j[k_0 + 2] + rtb_Q_j[k + 8]) / 2.0;
          rtb_y_g[k_0 + 3] = (rtb_Q_j[k_0 + 3] + rtb_Q_j[k + 12]) / 2.0;
          k_0 += 4;
        }

        // End of MATLAB Function: '<S134>/ScalarExpansion'
        // End of Outputs for SubSystem: '<S112>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S112>/ReducedQRN'
        for (k_0 = 0; k_0 < 4; k_0++) {
          // Product: '<S132>/Product' incorporates:
          //   Constant: '<S112>/G'
          //   Math: '<S132>/Transpose1'

          k = 0;
          for (ii = 0; ii < 4; ii++) {
            i = k + k_0;
            rtb_y_m[i] = 0.0;
            rtb_y_m[i] += rtb_y_g[k_0] * rtP.G_Value_a[ii];
            rtb_y_m[i] += rtb_y_g[k_0 + 4] * rtP.G_Value_a[ii + 4];
            rtb_y_m[i] += rtb_y_g[k_0 + 8] * rtP.G_Value_a[ii + 8];
            rtb_y_m[i] += rtb_y_g[k_0 + 12] * rtP.G_Value_a[ii + 12];
            k += 4;
          }
        }

        // Product: '<S132>/Product' incorporates:
        //   Constant: '<S112>/G'

        k_0 = 0;
        for (k = 0; k < 4; k++) {
          for (ii = 0; ii <= 2; ii += 2) {
            i = ii + k_0;
            (void)_mm_storeu_pd(&rtb_Q_j[i], _mm_set1_pd(0.0));
            tmp_3 = _mm_loadu_pd(&rtb_Q_j[i]);
            (void)_mm_storeu_pd(&rtb_Q_j[i], _mm_add_pd(tmp_3, _mm_mul_pd
              (_mm_set1_pd(rtb_y_m[k_0]), _mm_loadu_pd(&rtP.G_Value_a[ii]))));
            tmp_3 = _mm_loadu_pd(&rtb_Q_j[i]);
            (void)_mm_storeu_pd(&rtb_Q_j[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_y_m[k_0 + 1]), _mm_loadu_pd(&rtP.G_Value_a[ii + 4])), tmp_3));
            tmp_3 = _mm_loadu_pd(&rtb_Q_j[i]);
            (void)_mm_storeu_pd(&rtb_Q_j[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_y_m[k_0 + 2]), _mm_loadu_pd(&rtP.G_Value_a[ii + 8])), tmp_3));
            tmp_3 = _mm_loadu_pd(&rtb_Q_j[i]);
            (void)_mm_storeu_pd(&rtb_Q_j[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_y_m[k_0 + 3]), _mm_loadu_pd(&rtP.G_Value_a[ii + 12])), tmp_3));
          }

          k_0 += 4;
        }

        // Math: '<S132>/Transpose2' incorporates:
        //   Constant: '<S112>/H'

        k_0 = 0;
        for (k = 0; k < 3; k++) {
          rtb_Transpose2[k_0] = rtP.H_Value_o[k];
          rtb_Transpose2[k_0 + 1] = rtP.H_Value_o[k + 3];
          rtb_Transpose2[k_0 + 2] = rtP.H_Value_o[k + 6];
          rtb_Transpose2[k_0 + 3] = rtP.H_Value_o[k + 9];
          k_0 += 4;
        }

        // End of Math: '<S132>/Transpose2'
        for (k_0 = 0; k_0 < 4; k_0++) {
          // Sum: '<S132>/Add' incorporates:
          //   Math: '<S132>/Transpose2'
          //   Product: '<S132>/Product1'

          k = 0;
          for (ii = 0; ii < 3; ii++) {
            i = k + k_0;
            rtb_Add_k[i] = (((rtb_Transpose2[k + 1] * rtb_y_g[k_0 + 4] +
                              rtb_Transpose2[k] * rtb_y_g[k_0]) +
                             rtb_Transpose2[k + 2] * rtb_y_g[k_0 + 8]) +
                            rtb_Transpose2[k + 3] * rtb_y_g[k_0 + 12]) +
              rtb_N_f[i];
            k += 4;
          }

          // End of Sum: '<S132>/Add'
        }

        // Product: '<S132>/Product2' incorporates:
        //   Constant: '<S112>/G'
        //   Product: '<S132>/Product3'
        //   Product: '<S132>/Product4'
        //   Sum: '<S132>/Add'

        k_0 = 0;
        for (k = 0; k < 3; k++) {
          for (ii = 0; ii <= 2; ii += 2) {
            i = ii + k_0;
            (void)_mm_storeu_pd(&rtb_Product2_bg[i], _mm_set1_pd(0.0));
            tmp_3 = _mm_loadu_pd(&rtb_Product2_bg[i]);
            (void)_mm_storeu_pd(&rtb_Product2_bg[i], _mm_add_pd(tmp_3,
              _mm_mul_pd(_mm_set1_pd(rtb_Add_k[k_0]), _mm_loadu_pd
                         (&rtP.G_Value_a[ii]))));
            tmp_3 = _mm_loadu_pd(&rtb_Product2_bg[i]);
            (void)_mm_storeu_pd(&rtb_Product2_bg[i], _mm_add_pd(_mm_mul_pd
              (_mm_set1_pd(rtb_Add_k[k_0 + 1]), _mm_loadu_pd(&rtP.G_Value_a[ii +
              4])), tmp_3));
            tmp_3 = _mm_loadu_pd(&rtb_Product2_bg[i]);
            (void)_mm_storeu_pd(&rtb_Product2_bg[i], _mm_add_pd(_mm_mul_pd
              (_mm_set1_pd(rtb_Add_k[k_0 + 2]), _mm_loadu_pd(&rtP.G_Value_a[ii +
              8])), tmp_3));
            tmp_3 = _mm_loadu_pd(&rtb_Product2_bg[i]);
            (void)_mm_storeu_pd(&rtb_Product2_bg[i], _mm_add_pd(_mm_mul_pd
              (_mm_set1_pd(rtb_Add_k[k_0 + 3]), _mm_loadu_pd(&rtP.G_Value_a[ii +
              12])), tmp_3));
          }

          k_0 += 4;
        }

        // Product: '<S132>/Product3' incorporates:
        //   Product: '<S132>/Product4'

        k_0 = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          // Product: '<S132>/Product2' incorporates:
          //   Product: '<S132>/Product3'
          //   Product: '<S132>/Product4'

          rtb_y_tmp = 0;
          for (k = 0; k < 3; k++) {
            // Product: '<S132>/Product3' incorporates:
            //   Product: '<S132>/Product2'

            rtb_Q_tmp = k + k_0;
            tmp_0[rtb_Q_tmp] = 0.0;

            // Product: '<S132>/Product4'
            rtb_N_h[rtb_Q_tmp] = 0.0;
            tmp_0[rtb_Q_tmp] += rtb_Add_k[ii] * rtP.H_Value_o[k];
            rtb_N_h[rtb_Q_tmp] += rtb_N_f[rtb_y_tmp] * rtb_Transpose2[ii];
            tmp_0[rtb_Q_tmp] += rtb_Add_k[ii + 1] * rtP.H_Value_o[k + 3];
            rtb_N_h[rtb_Q_tmp] += rtb_N_f[rtb_y_tmp + 1] * rtb_Transpose2[ii + 1];
            tmp_0[rtb_Q_tmp] += rtb_Add_k[ii + 2] * rtP.H_Value_o[k + 6];
            rtb_N_h[rtb_Q_tmp] += rtb_N_f[rtb_y_tmp + 2] * rtb_Transpose2[ii + 2];
            tmp_0[rtb_Q_tmp] += rtb_Add_k[ii + 3] * rtP.H_Value_o[k + 9];
            rtb_N_h[rtb_Q_tmp] += rtb_N_f[rtb_y_tmp + 3] * rtb_Transpose2[ii + 3];

            // Product: '<S132>/Product3' incorporates:
            //   Constant: '<S112>/H'
            //   Math: '<S132>/Transpose'
            //   Math: '<S132>/Transpose2'
            //   Product: '<S132>/Product2'
            //   Product: '<S132>/Product4'
            //   Sum: '<S132>/Add'

            rtb_y_tmp += 4;
          }

          // Product: '<S132>/Product2' incorporates:
          //   Constant: '<S112>/H'
          //   Math: '<S132>/Transpose'
          //   Math: '<S132>/Transpose2'
          //   Product: '<S132>/Product3'
          //   Product: '<S132>/Product4'
          //   Sum: '<S132>/Add'

          k_0 += 3;
          ii += 4;
        }

        // End of Outputs for SubSystem: '<S112>/ReducedQRN'
        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (k_0 = 0; k_0 <= 6; k_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          // Outputs for Atomic SubSystem: '<S112>/ReducedQRN'
          tmp_3 = _mm_loadu_pd(&tmp_0[k_0]);
          tmp_1 = _mm_loadu_pd(&rtb_N_h[k_0]);
          tmp_2 = _mm_loadu_pd(&rtb_y[k_0]);
          (void)_mm_storeu_pd(&rtb_R[k_0], _mm_add_pd(_mm_add_pd(tmp_3, tmp_1),
            tmp_2));

          // End of Outputs for SubSystem: '<S112>/ReducedQRN'
          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Outputs for Atomic SubSystem: '<S112>/ReducedQRN'
        for (k_0 = 8; k_0 < 9; k_0++) {
          // Sum: '<S132>/Add1'
          rtb_R[k_0] = (tmp_0[k_0] + rtb_N_h[k_0]) + rtb_y[k_0];
        }

        // End of Outputs for SubSystem: '<S112>/ReducedQRN'

        // Gain: '<S90>/umin_scale1' incorporates:
        //   Sum: '<S132>/Add1'

        //  See help of ctrlKalmanFilterDTCalculatePL.m
        // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S152>:1' 
        //    Copyright 2014 The MathWorks, Inc.
        // '<S152>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
        //  Determine if the Square-Root algorithm was used
        // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S154>:1' 
        // '<S154>:1:4' if isSqrtUsed
        umin_scale1_idx_0 = rtP.umin_scale1_Gain[0] * umax_incr[0];

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = umin_scale1_idx_0 - rtU.u0[0];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[0] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Gain: '<S90>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain[1] * umax_incr[1];

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = umin_scale1_idx_0 - rtU.u0[1];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[1] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Gain: '<S90>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain[2] * umax_incr[2];

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = umin_scale1_idx_0 - rtU.u0[2];

        // Outputs for Enabled SubSystem: '<S131>/MeasurementUpdate' incorporates:
        //   EnablePort: '<S155>/Enable'

        // Outputs for Atomic SubSystem: '<S112>/CalculatePL'
        // MATLAB Function: '<S114>/Discrete-Time KF - Calculate PLMZ' incorporates:
        //   Constant: '<S3>/Constant1'
        //   Constant: '<S3>/Constant13'
        //   DataTypeConversion: '<S112>/DataTypeConversionEnable'
        //   Delay: '<S112>/MemoryP'
        //   Delay: '<S112>/MemoryX'
        //   Product: '<S132>/Product'
        //   Product: '<S132>/Product2'
        //   Product: '<S155>/C[k]*xhat[k|k-1]'
        //   Product: '<S155>/D[k]*u[k]'
        //   Product: '<S155>/Product3'
        //   Sum: '<S132>/Add1'
        //   Sum: '<S155>/Add1'
        //   Sum: '<S155>/Sum'

        if (rtP.Constant1_Value_e != 0.0) {
          k_0 = 0;
          for (k = 0; k < 3; k++) {
            ii = 0;
            i = 0;
            for (rtb_y_tmp = 0; rtb_y_tmp < 4; rtb_y_tmp++) {
              rtb_Q_tmp = ii + k;
              rtb_Transpose2[rtb_y_tmp + k_0] = rtb_C_c[rtb_Q_tmp];
              rtb_N_f[rtb_Q_tmp] = 0.0;
              rtb_N_f[rtb_Q_tmp] += rtDW.MemoryP_DSTATE_e[i] * rtb_C_c[k];
              rtb_N_f[rtb_Q_tmp] += rtDW.MemoryP_DSTATE_e[i + 1] * rtb_C_c[k + 3];
              rtb_N_f[rtb_Q_tmp] += rtDW.MemoryP_DSTATE_e[i + 2] * rtb_C_c[k + 6];
              rtb_N_f[rtb_Q_tmp] += rtDW.MemoryP_DSTATE_e[i + 3] * rtb_C_c[k + 9];
              ii += 3;
              i += 4;
            }

            k_0 += 4;
          }

          for (k_0 = 0; k_0 < 3; k_0++) {
            k = 0;
            ii = 0;
            for (i = 0; i < 3; i++) {
              rtb_y_tmp = k + k_0;
              rtb_y[rtb_y_tmp] = (((rtb_Transpose2[ii + 1] * rtb_N_f[k_0 + 3] +
                                    rtb_Transpose2[ii] * rtb_N_f[k_0]) +
                                   rtb_Transpose2[ii + 2] * rtb_N_f[k_0 + 6]) +
                                  rtb_Transpose2[ii + 3] * rtb_N_f[k_0 + 9]) +
                rtb_R[rtb_y_tmp];
              k += 3;
              ii += 4;
            }
          }

          for (k_0 = 0; k_0 < 4; k_0++) {
            for (k = 0; k < 4; k++) {
              i = k << 2UL;
              ii = k_0 + i;
              rtb_y_m[ii] = 0.0;
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i] * rtb_A_e[k_0];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 1] * rtb_A_e[k_0 + 4];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 2] * rtb_A_e[k_0 + 8];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 3] * rtb_A_e[k_0 + 12];
            }

            for (k = 0; k < 3; k++) {
              ii = k << 2UL;
              i = ii + k_0;
              rtb_Add_k[i] = (((rtb_Transpose2[ii + 1] * rtb_y_m[k_0 + 4] +
                                rtb_Transpose2[ii] * rtb_y_m[k_0]) +
                               rtb_Transpose2[ii + 2] * rtb_y_m[k_0 + 8]) +
                              rtb_Transpose2[ii + 3] * rtb_y_m[k_0 + 12]) +
                rtb_Product2_bg[i];
            }
          }

          mrdiv_c(rtb_Add_k, rtb_y, rtb_N_f);
          k_0 = 0;
          for (k = 0; k < 3; k++) {
            for (ii = 0; ii <= 2; ii += 2) {
              i = ii + k_0;
              (void)_mm_storeu_pd(&rtb_Add_k[i], _mm_set1_pd(0.0));
              tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE_e[ii]);
              tmp_1 = _mm_loadu_pd(&rtb_Add_k[i]);
              (void)_mm_storeu_pd(&rtb_Add_k[i], _mm_add_pd(tmp_1, _mm_mul_pd
                (_mm_set1_pd(rtb_Transpose2[k_0]), tmp_3)));
              tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE_e[ii + 4]);
              tmp_1 = _mm_loadu_pd(&rtb_Add_k[i]);
              (void)_mm_storeu_pd(&rtb_Add_k[i], _mm_add_pd(_mm_mul_pd
                (_mm_set1_pd(rtb_Transpose2[k_0 + 1]), tmp_3), tmp_1));
              tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE_e[ii + 8]);
              tmp_1 = _mm_loadu_pd(&rtb_Add_k[i]);
              (void)_mm_storeu_pd(&rtb_Add_k[i], _mm_add_pd(_mm_mul_pd
                (_mm_set1_pd(rtb_Transpose2[k_0 + 2]), tmp_3), tmp_1));
              tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE_e[ii + 12]);
              tmp_1 = _mm_loadu_pd(&rtb_Add_k[i]);
              (void)_mm_storeu_pd(&rtb_Add_k[i], _mm_add_pd(_mm_mul_pd
                (_mm_set1_pd(rtb_Transpose2[k_0 + 3]), tmp_3), tmp_1));
            }

            k_0 += 4;
          }

          mrdiv_c(rtb_Add_k, rtb_y, rtb_Transpose2);
          for (k_0 = 0; k_0 < 16; k_0++) {
            b_I[k_0] = 0;
          }

          b_I[0] = 1;
          b_I[5] = 1;
          b_I[10] = 1;
          b_I[15] = 1;
          for (k_0 = 0; k_0 < 4; k_0++) {
            for (k = 0; k < 4; k++) {
              i = (k << 2UL) + k_0;
              rtb_y_g[i] = static_cast<real_T>(b_I[i]) - ((rtb_C_c[3 * k + 1] *
                rtb_Transpose2[k_0 + 4] + rtb_C_c[3 * k] * rtb_Transpose2[k_0])
                + rtb_C_c[3 * k + 2] * rtb_Transpose2[k_0 + 8]);
            }

            for (k = 0; k < 4; k++) {
              ii = k << 2UL;
              i = k_0 + ii;
              Abar[i] = 0.0;
              Abar[i] += rtDW.MemoryP_DSTATE_e[ii] * rtb_y_g[k_0];
              Abar[i] += rtDW.MemoryP_DSTATE_e[ii + 1] * rtb_y_g[k_0 + 4];
              Abar[i] += rtDW.MemoryP_DSTATE_e[ii + 2] * rtb_y_g[k_0 + 8];
              Abar[i] += rtDW.MemoryP_DSTATE_e[ii + 3] * rtb_y_g[k_0 + 12];
            }

            for (k = 0; k < 3; k++) {
              rtb_Q_tmp = (k << 2UL) + k_0;
              rtb_Add_k[rtb_Q_tmp] = 0.0;
              rtb_Add_k[rtb_Q_tmp] += rtb_R[3 * k] * rtb_Transpose2[k_0];
              rtb_Add_k[rtb_Q_tmp] += rtb_R[3 * k + 1] * rtb_Transpose2[k_0 + 4];
              rtb_Add_k[rtb_Q_tmp] += rtb_R[3 * k + 2] * rtb_Transpose2[k_0 + 8];
            }
          }

          for (k_0 = 0; k_0 < 4; k_0++) {
            k = 0;
            for (ii = 0; ii < 4; ii++) {
              i = k + k_0;
              rtb_y_m[i] = 0.0;
              rtb_y_m[i] += Abar[k_0] * rtb_y_g[ii];
              rtb_y_m[i] += Abar[k_0 + 4] * rtb_y_g[ii + 4];
              rtb_y_m[i] += Abar[k_0 + 8] * rtb_y_g[ii + 8];
              rtb_y_m[i] += Abar[k_0 + 12] * rtb_y_g[ii + 12];
              rtb_Transpose2_0[i] = 0.0;
              rtb_Transpose2_0[i] += rtb_Add_k[k_0] * rtb_Transpose2[ii];
              rtb_Transpose2_0[i] += rtb_Add_k[k_0 + 4] * rtb_Transpose2[ii + 4];
              rtb_Transpose2_0[i] += rtb_Add_k[k_0 + 8] * rtb_Transpose2[ii + 8];
              k += 4;
            }
          }

          for (k_0 = 0; k_0 <= 14; k_0 += 2) {
            tmp_3 = _mm_loadu_pd(&rtb_y_m[k_0]);
            tmp_1 = _mm_loadu_pd(&rtb_Transpose2_0[k_0]);
            (void)_mm_storeu_pd(&rtb_Z[k_0], _mm_add_pd(tmp_3, tmp_1));
          }

          mrdiv_c(rtb_Product2_bg, rtb_R, rtb_Transpose2);
          for (k_0 = 0; k_0 < 4; k_0++) {
            for (k = 0; k < 4; k++) {
              i = (k << 2UL) + k_0;
              rtb_y_g[i] = rtb_A_e[i] - ((rtb_C_c[3 * k + 1] *
                rtb_Transpose2[k_0 + 4] + rtb_C_c[3 * k] * rtb_Transpose2[k_0])
                + rtb_C_c[3 * k + 2] * rtb_Transpose2[k_0 + 8]);
            }

            for (k = 0; k < 4; k++) {
              ii = k << 2UL;
              i = k_0 + ii;
              Abar[i] = 0.0;
              Abar[i] += rtb_Z[ii] * rtb_y_g[k_0];
              Abar[i] += rtb_Z[ii + 1] * rtb_y_g[k_0 + 4];
              Abar[i] += rtb_Z[ii + 2] * rtb_y_g[k_0 + 8];
              Abar[i] += rtb_Z[ii + 3] * rtb_y_g[k_0 + 12];
            }
          }

          for (k_0 = 0; k_0 < 4; k_0++) {
            k = 0;
            for (ii = 0; ii < 4; ii++) {
              i = k + k_0;
              rtb_y_m[i] = (((Abar[k_0 + 4] * rtb_y_g[ii + 4] + Abar[k_0] *
                              rtb_y_g[ii]) + Abar[k_0 + 8] * rtb_y_g[ii + 8]) +
                            Abar[k_0 + 12] * rtb_y_g[ii + 12]) + rtb_Q_j[i];
              rtb_Transpose2_0[i] = 0.0;
              rtb_Transpose2_0[i] += rtb_Transpose2[k_0] * rtb_Product2_bg[ii];
              rtb_Transpose2_0[i] += rtb_Transpose2[k_0 + 4] *
                rtb_Product2_bg[ii + 4];
              rtb_Transpose2_0[i] += rtb_Transpose2[k_0 + 8] *
                rtb_Product2_bg[ii + 8];
              k += 4;
            }
          }

          for (k_0 = 0; k_0 <= 14; k_0 += 2) {
            tmp_3 = _mm_loadu_pd(&rtb_y_m[k_0]);
            tmp_1 = _mm_loadu_pd(&rtb_Transpose2_0[k_0]);
            (void)_mm_storeu_pd(&rtb_y_g[k_0], _mm_sub_pd(tmp_3, tmp_1));
          }

          rtDW.MeasurementUpdate_MODE = true;
          for (k_0 = 0; k_0 <= 0; k_0 += 2) {
            tmp_3 = _mm_set1_pd(0.0);
            (void)_mm_storeu_pd(&rtb_C_0[k_0], tmp_3);
            tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0]);
            tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
            (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
              _mm_set1_pd(rtDW.MemoryX_DSTATE_l[0])), tmp_2));
            tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 3]);
            tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
            (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
              _mm_set1_pd(rtDW.MemoryX_DSTATE_l[1])), tmp_2));
            tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 6]);
            tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
            (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
              _mm_set1_pd(rtDW.MemoryX_DSTATE_l[2])), tmp_2));
            tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 9]);
            tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
            (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
              _mm_set1_pd(rtDW.MemoryX_DSTATE_l[3])), tmp_2));
            (void)_mm_storeu_pd(&tmp[k_0], tmp_3);
            tmp_3 = _mm_loadu_pd(&tmp[k_0]);
            (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
              (&rtP.Constant13_Value_c[k_0]), _mm_set1_pd(rtb_Sum1[0])), tmp_3));
            tmp_3 = _mm_loadu_pd(&tmp[k_0]);
            (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
              (&rtP.Constant13_Value_c[k_0 + 3]), _mm_set1_pd(rtb_Sum1[1])),
              tmp_3));
            tmp_3 = _mm_loadu_pd(&tmp[k_0]);
            (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
              (&rtP.Constant13_Value_c[k_0 + 6]), _mm_set1_pd(rtb_Sum1[2])),
              tmp_3));
            tmp_3 = _mm_loadu_pd(&rtb_C_0[k_0]);
            tmp_1 = _mm_loadu_pd(&tmp[k_0]);
            tmp_2 = _mm_loadu_pd(&y__m[k_0]);
            (void)_mm_storeu_pd(&rtb_Sum6[k_0], _mm_sub_pd(tmp_2, _mm_add_pd
              (tmp_3, tmp_1)));
          }

          for (k_0 = 2; k_0 < 3; k_0++) {
            rtb_C_0[k_0] = 0.0;
            rtb_C_0[k_0] += rtb_C_c[k_0] * rtDW.MemoryX_DSTATE_l[0];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 3] * rtDW.MemoryX_DSTATE_l[1];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 6] * rtDW.MemoryX_DSTATE_l[2];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 9] * rtDW.MemoryX_DSTATE_l[3];
            tmp[k_0] = 0.0;
            tmp[k_0] += rtP.Constant13_Value_c[k_0] * rtb_Sum1[0];
            tmp[k_0] += rtP.Constant13_Value_c[k_0 + 3] * rtb_Sum1[1];
            tmp[k_0] += rtP.Constant13_Value_c[k_0 + 6] * rtb_Sum1[2];
            rtb_Sum6[k_0] = y__m[k_0] - (rtb_C_0[k_0] + tmp[k_0]);
          }

          for (k_0 = 0; k_0 <= 2; k_0 += 2) {
            (void)_mm_storeu_pd(&rtDW.Product3_c[k_0], _mm_set1_pd(0.0));
            tmp_3 = _mm_loadu_pd(&rtb_N_f[k_0]);
            tmp_1 = _mm_loadu_pd(&rtDW.Product3_c[k_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[k_0], _mm_add_pd(_mm_mul_pd
              (tmp_3, _mm_set1_pd(rtb_Sum6[0])), tmp_1));
            tmp_3 = _mm_loadu_pd(&rtb_N_f[k_0 + 4]);
            tmp_1 = _mm_loadu_pd(&rtDW.Product3_c[k_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[k_0], _mm_add_pd(_mm_mul_pd
              (tmp_3, _mm_set1_pd(rtb_Sum6[1])), tmp_1));
            tmp_3 = _mm_loadu_pd(&rtb_N_f[k_0 + 8]);
            tmp_1 = _mm_loadu_pd(&rtDW.Product3_c[k_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[k_0], _mm_add_pd(_mm_mul_pd
              (tmp_3, _mm_set1_pd(rtb_Sum6[2])), tmp_1));
          }
        } else {
          for (k_0 = 0; k_0 < 4; k_0++) {
            for (k = 0; k < 4; k++) {
              i = k << 2UL;
              ii = k_0 + i;
              rtb_y_m[ii] = 0.0;
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i] * rtb_A_e[k_0];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 1] * rtb_A_e[k_0 + 4];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 2] * rtb_A_e[k_0 + 8];
              rtb_y_m[ii] += rtDW.MemoryP_DSTATE_e[i + 3] * rtb_A_e[k_0 + 12];
            }

            for (k = 0; k < 4; k++) {
              i = (k << 2UL) + k_0;
              rtb_y_g[i] = (((rtb_y_m[k_0 + 4] * rtb_A_e[k + 4] + rtb_y_m[k_0] *
                              rtb_A_e[k]) + rtb_y_m[k_0 + 8] * rtb_A_e[k + 8]) +
                            rtb_y_m[k_0 + 12] * rtb_A_e[k + 12]) + rtb_Q_j[i];
            }
          }

          if (rtDW.MeasurementUpdate_MODE) {
            // Disable for Product: '<S155>/Product3' incorporates:
            //   Outport: '<S155>/L*(y[k]-yhat[k|k-1])'
            //
            rtDW.Product3_c[0] = rtP.Lykyhatkk1_Y0_c;
            rtDW.Product3_c[1] = rtP.Lykyhatkk1_Y0_c;
            rtDW.Product3_c[2] = rtP.Lykyhatkk1_Y0_c;
            rtDW.Product3_c[3] = rtP.Lykyhatkk1_Y0_c;
            rtDW.MeasurementUpdate_MODE = false;
          }
        }

        // End of MATLAB Function: '<S114>/Discrete-Time KF - Calculate PLMZ'
        // End of Outputs for SubSystem: '<S112>/CalculatePL'
        // End of Outputs for SubSystem: '<S131>/MeasurementUpdate'

        // Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'
        //   Sum: '<S3>/Sum'

        rtDW.DiscreteTimeIntegrator_DSTATE_j += (y__m[0] - Sum2_c[0]) *
          rtP.DiscreteTimeIntegrator_gainval;
        rtDW.DiscreteTimeIntegrator_PrevRe_b = static_cast<int8_T>(rtU.iRST ? 1 :
          0);

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (k_0 = 0; k_0 <= 0; k_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_3 = _mm_set1_pd(0.0);
          (void)_mm_storeu_pd(&rtb_C_0[k_0], tmp_3);
          tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0]);
          tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
          (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtDW.MemoryX_DSTATE_l[0])), tmp_2));
          tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 3]);
          tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
          (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtDW.MemoryX_DSTATE_l[1])), tmp_2));
          tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 6]);
          tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
          (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtDW.MemoryX_DSTATE_l[2])), tmp_2));
          tmp_1 = _mm_loadu_pd(&rtb_C_c[k_0 + 9]);
          tmp_2 = _mm_loadu_pd(&rtb_C_0[k_0]);
          (void)_mm_storeu_pd(&rtb_C_0[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtDW.MemoryX_DSTATE_l[3])), tmp_2));
          (void)_mm_storeu_pd(&tmp[k_0], tmp_3);
          tmp_3 = _mm_loadu_pd(&tmp[k_0]);
          (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[k_0]), _mm_set1_pd(rtb_Sum1[0])), tmp_3));
          tmp_3 = _mm_loadu_pd(&tmp[k_0]);
          (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[k_0 + 3]), _mm_set1_pd(rtb_Sum1[1])), tmp_3));
          tmp_3 = _mm_loadu_pd(&tmp[k_0]);
          (void)_mm_storeu_pd(&tmp[k_0], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[k_0 + 6]), _mm_set1_pd(rtb_Sum1[2])), tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_C_0[k_0]);
          tmp_1 = _mm_loadu_pd(&tmp[k_0]);
          (void)_mm_storeu_pd(&rtb_Sum6[k_0], _mm_add_pd(tmp_3, tmp_1));
          tmp_3 = _mm_loadu_pd(&umax_incr[k_0]);
          (void)_mm_storeu_pd(&rtDW.last_mv_DSTATE_n[k_0], tmp_3);

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        for (k_0 = 2; k_0 < 3; k_0++) {
          // Product: '<S115>/Product'
          rtb_C_0[k_0] = 0.0;
          rtb_C_0[k_0] += rtb_C_c[k_0] * rtDW.MemoryX_DSTATE_l[0];
          rtb_C_0[k_0] += rtb_C_c[k_0 + 3] * rtDW.MemoryX_DSTATE_l[1];
          rtb_C_0[k_0] += rtb_C_c[k_0 + 6] * rtDW.MemoryX_DSTATE_l[2];
          rtb_C_0[k_0] += rtb_C_c[k_0 + 9] * rtDW.MemoryX_DSTATE_l[3];

          // Product: '<S115>/Product1' incorporates:
          //   Delay: '<S112>/MemoryX'
          //   Product: '<S115>/Product'

          tmp[k_0] = 0.0;
          tmp[k_0] += rtP.Constant13_Value_c[k_0] * rtb_Sum1[0];
          tmp[k_0] += rtP.Constant13_Value_c[k_0 + 3] * rtb_Sum1[1];
          tmp[k_0] += rtP.Constant13_Value_c[k_0 + 6] * rtb_Sum1[2];

          // Sum: '<S115>/Add1' incorporates:
          //   Constant: '<S3>/Constant13'
          //   Product: '<S115>/Product'
          //   Product: '<S115>/Product1'

          rtb_Sum6[k_0] = rtb_C_0[k_0] + tmp[k_0];

          // Update for UnitDelay: '<S90>/last_mv' incorporates:
          //   Product: '<S115>/Product'

          rtDW.last_mv_DSTATE_n[k_0] = umax_incr[k_0];
        }

        // Update for Delay: '<S112>/MemoryX' incorporates:
        //   Constant: '<S3>/Constant13'
        //   Product: '<S115>/Product'
        //   Product: '<S115>/Product1'
        //   Sum: '<S115>/Add1'
        //   UnitDelay: '<S90>/last_mv'

        rtDW.icLoad_n = false;

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (k_0 = 0; k_0 <= 2; k_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_3 = _mm_set1_pd(0.0);
          (void)_mm_storeu_pd(&rtb_Sum2_f[k_0], tmp_3);
          tmp_1 = _mm_loadu_pd(&Sum_h[k_0]);
          tmp_2 = _mm_loadu_pd(&rtb_Sum2_f[k_0]);
          (void)_mm_storeu_pd(&rtb_Sum2_f[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtb_Sum1[0])), tmp_2));
          tmp_1 = _mm_loadu_pd(&Sum_h[k_0 + 4]);
          tmp_2 = _mm_loadu_pd(&rtb_Sum2_f[k_0]);
          (void)_mm_storeu_pd(&rtb_Sum2_f[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtb_Sum1[1])), tmp_2));
          tmp_1 = _mm_loadu_pd(&Sum_h[k_0 + 8]);
          tmp_2 = _mm_loadu_pd(&rtb_Sum2_f[k_0]);
          (void)_mm_storeu_pd(&rtb_Sum2_f[k_0], _mm_add_pd(_mm_mul_pd(tmp_1,
            _mm_set1_pd(rtb_Sum1[2])), tmp_2));
          (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0], tmp_3);
          tmp_3 = _mm_loadu_pd(&rtb_A_e[k_0]);
          tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0]);
          (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0], _mm_add_pd
                              (_mm_mul_pd(tmp_3, _mm_set1_pd
            (rtDW.MemoryX_DSTATE_l[0])), tmp_1));
          tmp_3 = _mm_loadu_pd(&rtb_A_e[k_0 + 4]);
          tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0]);
          (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0], _mm_add_pd
                              (_mm_mul_pd(tmp_3, _mm_set1_pd
            (rtDW.MemoryX_DSTATE_l[1])), tmp_1));
          tmp_3 = _mm_loadu_pd(&rtb_A_e[k_0 + 8]);
          tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0]);
          (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0], _mm_add_pd
                              (_mm_mul_pd(tmp_3, _mm_set1_pd
            (rtDW.MemoryX_DSTATE_l[2])), tmp_1));
          tmp_3 = _mm_loadu_pd(&rtb_A_e[k_0 + 12]);
          tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0]);
          (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFu_ia[k_0], _mm_add_pd
                              (_mm_mul_pd(tmp_3, _mm_set1_pd
            (rtDW.MemoryX_DSTATE_l[3])), tmp_1));

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Update for Delay: '<S112>/MemoryX' incorporates:
        //   Product: '<S131>/A[k]*xhat[k|k-1]'
        //   Product: '<S131>/B[k]*u[k]'
        //   Sum: '<S131>/Add'

        rtDW.MemoryX_DSTATE_l[0] = (rtb_Sum2_f[0] +
          rtb_TmpSignalConversionAtSFu_ia[0]) + rtDW.Product3_c[0];
        rtDW.MemoryX_DSTATE_l[1] = (rtb_Sum2_f[1] +
          rtb_TmpSignalConversionAtSFu_ia[1]) + rtDW.Product3_c[1];
        rtDW.MemoryX_DSTATE_l[2] = (rtb_Sum2_f[2] +
          rtb_TmpSignalConversionAtSFu_ia[2]) + rtDW.Product3_c[2];
        rtDW.MemoryX_DSTATE_l[3] = (rtb_Sum2_f[3] +
          rtb_TmpSignalConversionAtSFu_ia[3]) + rtDW.Product3_c[3];

        // Update for Delay: '<S112>/MemoryP'
        rtDW.icLoad_h = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE_e[0], &rtb_y_g[0], sizeof(real_T)
                          << 4UL);
        rtY.yhat[0] = rtb_Sum6[0];
        rtY.yhat[1] = rtb_Sum6[1];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S89>/Sum3'

        rtY.u[2] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        rtY.yhat[2] = rtb_Sum6[2];

        // End of Outputs for SubSystem: '<S1>/mpc1'
      } else if (rtY.sig == 2.0) {
        real_T Bc_2;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // DiscreteIntegrator: '<S4>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'

        // '<S1>:59:35' elseif sig == 2
        // '<S1>:59:36' [u, yhat(1:no)] = mpc2(r_, y__, [0;0;0], [0;0], u0, umax, iRST); 
        // Simulink Function 'mpc2': '<S1>:902'
        if (rtU.iRST && (rtDW.DiscreteTimeIntegrator_PrevRe_f <= 0)) {
          rtDW.DiscreteTimeIntegrator_DSTATE_m[0] =
            rtP.DiscreteTimeIntegrator_IC_n[0];
          rtDW.DiscreteTimeIntegrator_DSTATE_m[1] =
            rtP.DiscreteTimeIntegrator_IC_n[1];
        }

        // Delay: '<S182>/MemoryX' incorporates:
        //   Constant: '<S182>/X0'
        //   DataTypeConversion: '<S182>/DataTypeConversionReset'

        rtDW.icLoad_a = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE_g) ==
                          POS_ZCSIG) || rtDW.icLoad_a);
        rtPrevZCX.MemoryX_Reset_ZCE_g = 0U;
        if (rtDW.icLoad_a) {
          for (k = 0; k < 5; k++) {
            rtDW.MemoryX_DSTATE_c[k] = rtP.X0_Value_k[k];
          }
        }

        // SignalConversion generated from: '<S181>/ SFunction ' incorporates:
        //   Constant: '<S4>/Constant'
        //   MATLAB Function: '<S180>/optimizer'

        rtb_TmpSignalConversionAtSFu_o4[0] = Sum2_c[0];
        rtb_TmpSignalConversionAtSFu_o4[1] = Sum2_c[1];
        rtb_TmpSignalConversionAtSFu_o4[2] = Sum2_c[2];
        rtb_TmpSignalConversionAtSFu_o4[3] = rtP.Constant_Value_p[0];
        rtb_TmpSignalConversionAtSFu_o4[4] = rtP.Constant_Value_p[1];

        // MATLAB Function: '<S180>/optimizer' incorporates:
        //   Constant: '<S159>/Constant1'
        //   Delay: '<S182>/MemoryX'
        //   DiscreteIntegrator: '<S4>/Discrete-Time Integrator'
        //   Gain: '<S4>/Gain2'
        //   Inport: '<Root>/umax'
        //   SignalConversion generated from: '<S181>/ SFunction '
        //   Sum: '<S159>/Sum2'
        //   UnitDelay: '<S160>/last_mv'

        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S181>:1'
        // '<S181>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S181>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S181>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S181>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S181>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S181>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S181>:1:252' if isa(u,'double')
        // '<S181>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S181>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S181>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S181>:1:252' if isa(u,'double')
        // '<S181>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S181>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S181>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S181>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S181>:1:252' if isa(u,'double')
        // '<S181>:1:253' y = u;
        //    md:         measured disturbance
        // '<S181>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S181>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S181>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S181>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S181>:1:252' if isa(u,'double')
        // '<S181>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S181>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S181>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S181>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S181>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S181>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S181>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S181>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S181>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S181>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S181>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S181>:1:58' uwt = convertDataType(uwt0,isDouble);
        //    duwt:       run-time DMV weights
        // '<S181>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S181>:1:62' ewt = convertDataType(ewt0,isDouble);
        //    iA:         logical previous active set (for warm start)
        //  Outputs (in BlockDataType except iAout)
        //    xk1:        x[k+1|k] from built-in KF
        //    u:          optimal MV
        //    cost:       optimal cost
        //    useq:       optimal MV sequence
        //    xseq:       optimal state sequence
        //    yseq:       optimal OV sequence
        //    status:     QP exit flag
        //    xest:       x[k|k] from built-in KF
        //    iAout:      logical current active set
        //  Parameters (constant)
        //    dimensions (int32):
        //        nx, nxp, nup, nu, ny, degrees, p, nxQP, enable_value, Mrows, nCC, nv 
        //        myindex, mvindex, mdindex, unindex, nxid, m, Ndis, numdis, maxdis 
        //    MPC constants (BlockDataType):
        //        Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, Mx, Mu1, Mv, utarget
        //        H, Linv, Ac, Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1
        //        A, C, B, D, Cid, Did, Ecc, Fcc, Scc, Gcc
        //        RYscale, RMDscale, xoff, Uscale, Yscale
        //        uoff, voff, yoff, myoff, RMVscale, Mdis, Vdis
        //    configurations (logical):
        //        isQP, CustomSolver, CustomSolverCodeGen, UseSuboptimalSolution, UseActiveSetSolver 
        //        openloopflag, no_umin, no_umax, no_ymin, no_ymax, switch_inport, no_switch 
        //        return_cost, return_mvseq, return_xseq, return_ovseq, isLTV
        //        no_ywt, no_uwt, no_duwt, no_rhoeps, no_md, no_ref, no_uref, no_mv 
        //        CustomEstimation, no_cc, isHyb, isDouble
        //    ASOptions
        //    IPOptions
        //    MIQPOptions
        //  Parameters
        // '<S181>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S181>:1:96' isAdaptive = false;
        // '<S181>:1:97' isLTV = false;
        // '<S181>:1:98' ZERO = zeros('like',ref);
        // '<S181>:1:99' ONE = ones('like',ref);
        // '<S181>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S181>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S181>:1:119' if isSimulation
        // '<S181>:1:126' else
        //  When doing code generation, use M code directly
        // '<S181>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        k_0 = 0;
        for (k = 0; k < 20; k++) {
          for (ii = 0; ii < 5; ii++) {
            rseq_0[ii + k_0] = rtb_TmpSignalConversionAtSFu_o4[ii];
          }

          k_0 += 5;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S181>:1:133' old_u = old_u - uoff;
        umin_scale1_idx_0 = rtDW.last_mv_DSTATE_i[0];
        umin_scale1_idx_1 = rtDW.last_mv_DSTATE_i[1];
        umin_scale1_idx_2 = rtDW.last_mv_DSTATE_i[2];

        // '<S181>:1:134' if no_mv
        // '<S181>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S181>:1:143' xk = xk - xoff;
        rtb_xest[0] = rtDW.MemoryX_DSTATE_c[0];
        rtb_xest[2] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE_m[0];
        rtb_xest[1] = rtDW.MemoryX_DSTATE_c[1];
        rtb_xest[3] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE_m[1];
        rtb_xest[4] = rtP.Constant1_Value_p[0] + rtDW.MemoryX_DSTATE_c[2];
        rtb_xest[5] = rtP.Constant1_Value_p[1] + rtDW.MemoryX_DSTATE_c[3];
        rtb_xest[6] = rtP.Constant1_Value_p[2] + rtDW.MemoryX_DSTATE_c[4];

        //  Remove offset
        // '<S181>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S181>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S181>:1:162' if no_uref
        //  no external utarget
        // '<S181>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S181>:1:173' if ~no_cc
        // '<S181>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S181>:1:183' if isSimulation
        // '<S181>:1:214' else
        //  When doing code generation, use M code directly
        // '<S181>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S181>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S181>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S181>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S181>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S181>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S181>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S181>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S181>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S181>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S181>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S181>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S181>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        umax_incr_flag[0] = false;
        umax_incr[0] = 0.0;
        umax_incr_flag[1] = false;
        umax_incr[1] = 0.0;
        umax_incr_flag[2] = false;
        umax_incr[2] = 0.0;
        for (k = 0; k < 206; k++) {
          int16_T b_Mrows_0;
          dwt = b_Mlim_1[k];
          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 7; k_0++) {
            Bc_2 += b_a_0[206 * k_0 + k] * rtb_xest[k_0];
          }

          Bc_2 = -(((a_0[k + 206] * umin_scale1_idx_1 + a_0[k] *
                     umin_scale1_idx_0) + a_0[k + 412] * umin_scale1_idx_2) +
                   (dwt + Bc_2));
          b_Mrows_0 = b_Mrows_2[k];
          if ((b_Mrows_0 > 100) && (b_Mrows_0 > 200) && (b_Mrows_0 <= 260)) {
            ii = (static_cast<int32_T>(b_Mrows_0) - div_nde_s32_floor(
                   static_cast<int32_T>(b_Mrows_0) - 201, static_cast<int32_T>
                   (nu)) * static_cast<int32_T>(nu)) - 201;
            rstP2 = umax_incr_flag[ii];
            if (!umax_incr_flag[ii]) {
              dwt = -rtU.umax[ii] - (-dwt);
              rstP2 = true;
            } else {
              dwt = umax_incr[ii];
            }

            umax_incr[ii] = dwt;
            umax_incr_flag[ii] = rstP2;
            Bc_2 += dwt;
          }

          Bc_0[k] = Bc_2;
        }

        rtb_Sum2_f[0] = 0.0;
        rtb_Sum2_f[1] = 0.0;
        rtb_Sum2_f[2] = 0.0;
        rtb_Sum2_f[3] = 0.0;
        for (k = 0; k < 3; k++) {
          dwt = 0.0;
          for (k_0 = 0; k_0 < 7; k_0++) {
            dwt += b_Kx_0[7 * k + k_0] * rtb_xest[k_0];
          }

          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 100; k_0++) {
            Bc_2 += b_Kr_0[100 * k + k_0] * rseq_0[k_0];
          }

          rtb_Sum2_f[k] = ((b_Ku1_0[3 * k + 1] * umin_scale1_idx_1 + b_Ku1_0[3 *
                            k] * umin_scale1_idx_0) + b_Ku1_0[3 * k + 2] *
                           umin_scale1_idx_2) + (dwt + Bc_2);
        }

        // Update for Memory: '<S160>/Memory' incorporates:
        //   MATLAB Function: '<S180>/optimizer'

        qpkwik_o(b_Linv_0, b_Hinv_0, rtb_Sum2_f, b_Ac_0, Bc_0,
                 rtDW.Memory_PreviousInput_c, 840, 1.0E-6,
                 rtb_TmpSignalConversionAtSFu_ia, a__1_0, &k);

        // MATLAB Function: '<S180>/optimizer' incorporates:
        //   UnitDelay: '<S160>/last_mv'

        if ((k < 0) || (k == 0)) {
          rtb_TmpSignalConversionAtSFu_ia[0] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[1] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[2] = 0.0;
        }

        umax_incr[0] = rtDW.last_mv_DSTATE_i[0] +
          rtb_TmpSignalConversionAtSFu_ia[0];
        umax_incr[1] = rtDW.last_mv_DSTATE_i[1] +
          rtb_TmpSignalConversionAtSFu_ia[1];
        umax_incr[2] = rtDW.last_mv_DSTATE_i[2] +
          rtb_TmpSignalConversionAtSFu_ia[2];

        // Delay: '<S182>/MemoryP' incorporates:
        //   Constant: '<S182>/P0'
        //   DataTypeConversion: '<S182>/DataTypeConversionReset'

        // '<S181>:1:231' if return_xseq || return_ovseq
        // '<S181>:1:233' else
        // '<S181>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S181>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S181>:1:238' if CustomEstimation
        // '<S181>:1:239' xk1 = xk;
        // '<S181>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S181>:1:247' xest = xest + xoff;
        rtDW.icLoad_p = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE_i) ==
                          POS_ZCSIG) || rtDW.icLoad_p);
        rtPrevZCX.MemoryP_Reset_ZCE_i = 0U;
        if (rtDW.icLoad_p) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE_h[0], &rtP.P0_Value_c[0], 25U *
                            sizeof(real_T));
        }

        // MATLAB Function: '<S159>/MATLAB Function' incorporates:
        //   BusCreator: '<S4>/Bus Creator1'
        //   Constant: '<S4>/Constant12'
        //   Constant: '<S4>/Constant13'
        //   Constant: '<S4>/Constant3'
        //   Constant: '<S4>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc2/State Estimator OD (KF)/MATLAB Function': '<S183>:1' 
        // '<S183>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod2, Bod2, Cod2, Dod2, Dmn1, 3, 3, 2); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&rtb_B[0], 0, 15U * sizeof(real_T));

        //  n_states  x n_inputs
        // 'stateEst_:9' C = zeros(no,ns);
        //  n_outputs x n_states
        // 'stateEst_:10' D = zeros(no,ni);
        //  n_outputs x n_inputs
        // 'stateEst_:11' Q = zeros(ns,ns);
        //  n_states  x n_states
        // 'stateEst_:12' G = eye(ns);
        //  n_states  x n_states
        // 'stateEst_:13' R = zeros(no,no);
        //  n_outputs x n_outputs
        // 'stateEst_:14' N = zeros(ns,no);
        //  n_states  x n_outputs
        // 'stateEst_:15' H = zeros(no,ns);
        //  n_outputs x n_states
        //  combine plant and output disturbance model
        //  (force the outputs to fit in preallocated memory)
        // 'stateEst_:19' A(1:ns, 1:ns) = blkdiag(Ap, Aod);
        (void)std::memset(&rtb_A[0], 0, 25U * sizeof(real_T));
        rtb_A[0] = rtP.Constant3_Value_d[0];
        rtb_A[1] = rtP.Constant3_Value_d[1];
        rtb_A[5] = rtP.Constant3_Value_d[2];
        rtb_A[6] = rtP.Constant3_Value_d[3];

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        k_0 = 0;
        k = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          rtb_A[k_0 + 12] = rtP.Aod2[k];
          rtb_A[k_0 + 13] = rtP.Aod2[k + 1];
          rtb_A[k_0 + 14] = rtP.Aod2[k + 2];
          rtb_B[k_0] = rtP.Constant4_Value_n[ii];
          rtb_B[k_0 + 1] = rtP.Constant4_Value_n[ii + 1];
          k_0 += 5;
          k += 3;
          ii += 2;
        }

        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        for (k = 0; k < 6; k++) {
          rtb_C[k] = rtP.Constant12_Value_i[k];
        }

        (void)std::memcpy(&rtb_C[6], &rtP.Cod2[0], 9U * sizeof(real_T));

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est_0[0], 0, 45U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y_0[0], 0, 30U * sizeof(real_T));
        k_0 = 0;
        k = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          y_0[k_0] = rtP.Constant4_Value_n[k];
          y_0[k_0 + 1] = rtP.Constant4_Value_n[k + 1];
          y_0[k_0 + 17] = rtP.Bod2[ii];
          y_0[k_0 + 18] = rtP.Bod2[ii + 1];
          y_0[k_0 + 19] = rtP.Bod2[ii + 2];
          k_0 += 5;
          k += 2;
          ii += 3;
        }

        (void)std::memcpy(&B_est_0[0], &y_0[0], 30U * sizeof(real_T));

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        for (k_0 = 0; k_0 < 9; k_0++) {
          D_est[k_0] = rtP.Constant13_Value_g[k_0];
          D_est[k_0 + 9] = rtP.Dod2[k_0];
          D_est[k_0 + 18] = rtP.Dmn1[k_0];
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (k_0 = 0; k_0 < 5; k_0++) {
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            int32_T rtb_Q_tmp;
            rtb_Q_tmp = k + k_0;
            rtb_Q[rtb_Q_tmp] = 0.0;
            i = 0;
            for (int32_T rtb_y_tmp{0}; rtb_y_tmp < 9; rtb_y_tmp++) {
              rtb_Q[rtb_Q_tmp] += B_est_0[i + k_0] * B_est_0[i + ii];
              i += 5;
            }

            k += 5;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        k_0 = 0;
        for (k = 0; k < 9; k++) {
          rtb_R_tmp[k] = D_est[k_0];
          rtb_R_tmp[k + 9] = D_est[k_0 + 1];
          rtb_R_tmp[k + 18] = D_est[k_0 + 2];
          k_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (k_0 = 0; k_0 < 3; k_0++) {
          for (k = 0; k < 3; k++) {
            i = 3 * k_0 + k;
            rtb_R[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_R[i] += D_est[3 * ii + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }

          for (k = 0; k < 5; k++) {
            i = 5 * k_0 + k;
            rtb_N[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_N[i] += B_est_0[5 * ii + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }
        }

        // End of MATLAB Function: '<S159>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S182>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S182>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S182>/ScalarExpansionQ'
        ScalarExpansionQ(rtb_Q, rtb_Z_e);

        // End of Outputs for SubSystem: '<S182>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S182>/ReducedQRN'
        ReducedQRN(rtP.G_Value_g, rtP.H_Value_k, rtb_Z_e, rtb_y, rtb_N,
                   rtb_Product, rtb_R, rtb_Product2);

        // End of Outputs for SubSystem: '<S182>/ReducedQRN'

        // Outputs for Atomic SubSystem: '<S182>/CalculatePL'
        CalculatePL(rtb_A, rtb_C, rtb_Product, rtb_R, rtb_Product2,
                    rtP.Constant1_Value_pe != 0.0, rtDW.MemoryP_DSTATE_h, rtb_N,
                    rtb_L, rtb_Z_e, rtb_Q);

        // End of Outputs for SubSystem: '<S182>/CalculatePL'

        // MATLAB Function: '<S223>/SqrtUsedFcn' incorporates:
        //   Constant: '<S182>/G'
        //   Constant: '<S182>/H'
        //   Constant: '<S223>/isSqrtUsed'
        //   Constant: '<S4>/Constant1'
        //   DataTypeConversion: '<S182>/DataTypeConversionEnable'
        //   Delay: '<S182>/MemoryP'

        SqrtUsedFcn(rtb_Z_e, rtP.isSqrtUsed_Value_a, rtb_Product);

        // Gain: '<S160>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_p[0] * umax_incr[0];

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = umin_scale1_idx_0 - rtU.u0[0];

        // Sum: '<S159>/Sum6'
        rtb_Sum6[0] = y__m[0];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[0] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Gain: '<S160>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_p[1] * umax_incr[1];

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = umin_scale1_idx_0 - rtU.u0[1];

        // Sum: '<S159>/Sum6'
        rtb_Sum6[1] = y__m[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[1] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Gain: '<S160>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_p[2] * umax_incr[2];

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = umin_scale1_idx_0 - rtU.u0[2];

        // Sum: '<S159>/Sum6'
        rtb_Sum6[2] = y__m[2];

        // Outputs for Enabled SubSystem: '<S201>/MeasurementUpdate'
        MeasurementUpdate(rtP.Constant1_Value_pe != 0.0, rtb_L, rtb_Sum6, rtb_C,
                          rtDW.MemoryX_DSTATE_c, rtP.Constant13_Value_g,
                          rtb_Sum1, rtDW.Product3_a, &rtDW.MeasurementUpdate_j,
                          &rtP.MeasurementUpdate_j);

        // End of Outputs for SubSystem: '<S201>/MeasurementUpdate'
        for (k_0 = 0; k_0 < 3; k_0++) {
          // Product: '<S185>/Product'
          rtb_C_0[k_0] = 0.0;
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            rtb_C_0[k_0] += rtb_C[k + k_0] * rtDW.MemoryX_DSTATE_c[ii];
            k += 3;
          }

          // Product: '<S185>/Product1' incorporates:
          //   Product: '<S185>/Product'

          tmp[k_0] = 0.0;
          tmp[k_0] += rtP.Constant13_Value_g[k_0] * rtb_Sum1[0];
          tmp[k_0] += rtP.Constant13_Value_g[k_0 + 3] * rtb_Sum1[1];
          tmp[k_0] += rtP.Constant13_Value_g[k_0 + 6] * rtb_Sum1[2];

          // Sum: '<S185>/Add1' incorporates:
          //   Product: '<S185>/Product'
          //   Product: '<S185>/Product1'

          rtb_Sum6[k_0] = rtb_C_0[k_0] + tmp[k_0];
        }

        // Update for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' incorporates:
        //   Constant: '<S4>/Constant1'
        //   Constant: '<S4>/Constant13'
        //   DataTypeConversion: '<S182>/DataTypeConversionEnable'
        //   Delay: '<S182>/MemoryX'
        //   Inport: '<Root>/iRST'
        //   Product: '<S185>/Product'
        //   Product: '<S185>/Product1'
        //   Sum: '<S4>/Sum'

        rtDW.DiscreteTimeIntegrator_DSTATE_m[0] += (y__m[1] - Sum2_c[1]) *
          rtP.DiscreteTimeIntegrator_gainva_b;
        rtDW.DiscreteTimeIntegrator_DSTATE_m[1] += (y__m[2] - Sum2_c[2]) *
          rtP.DiscreteTimeIntegrator_gainva_b;
        rtDW.DiscreteTimeIntegrator_PrevRe_f = static_cast<int8_T>(rtU.iRST ? 1 :
          0);

        // Update for UnitDelay: '<S160>/last_mv'
        rtDW.last_mv_DSTATE_i[0] = umax_incr[0];
        rtDW.last_mv_DSTATE_i[1] = umax_incr[1];
        rtDW.last_mv_DSTATE_i[2] = umax_incr[2];

        // Update for Delay: '<S182>/MemoryX'
        rtDW.icLoad_a = false;
        for (k_0 = 0; k_0 < 5; k_0++) {
          // Product: '<S201>/B[k]*u[k]'
          rtb_TmpSignalConversionAtSFu_o4[k_0] = 0.0;
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0] * rtb_Sum1[0];
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0 + 5] * rtb_Sum1[1];
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0 + 10] * rtb_Sum1[2];

          // Product: '<S201>/A[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S182>/MemoryX'
          //   Product: '<S201>/B[k]*u[k]'

          rtb_A_0[k_0] = 0.0;
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            rtb_A_0[k_0] += rtb_A[k + k_0] * rtDW.MemoryX_DSTATE_c[ii];
            k += 5;
          }

          // End of Product: '<S201>/A[k]*xhat[k|k-1]'
        }

        // End of Outputs for SubSystem: '<S1>/mpc2'
        for (k_0 = 0; k_0 <= 2; k_0 += 2) {
          __m128d tmp_1;
          __m128d tmp_2;

          // Outputs for Function Call SubSystem: '<S1>/mpc2'
          tmp_3 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_o4[k_0]);
          tmp_1 = _mm_loadu_pd(&rtb_A_0[k_0]);
          tmp_2 = _mm_loadu_pd(&rtDW.Product3_a[k_0]);
          (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE_c[k_0], _mm_add_pd(_mm_add_pd
            (tmp_3, tmp_1), tmp_2));

          // End of Outputs for SubSystem: '<S1>/mpc2'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        for (k_0 = 4; k_0 < 5; k_0++) {
          // Update for Delay: '<S182>/MemoryX' incorporates:
          //   Sum: '<S201>/Add'

          rtDW.MemoryX_DSTATE_c[k_0] = (rtb_TmpSignalConversionAtSFu_o4[k_0] +
            rtb_A_0[k_0]) + rtDW.Product3_a[k_0];
        }

        // Update for Delay: '<S182>/MemoryP' incorporates:
        //   Delay: '<S182>/MemoryX'
        //   Product: '<S201>/B[k]*u[k]'
        //   Sum: '<S201>/Add'

        rtDW.icLoad_p = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE_h[0], &rtb_Q[0], 25U * sizeof
                          (real_T));
        rtY.yhat[0] = rtb_Sum6[0];
        rtY.yhat[1] = rtb_Sum6[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S159>/Sum3'

        rtY.u[2] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        rtY.yhat[2] = rtb_Sum6[2];

        // End of Outputs for SubSystem: '<S1>/mpc2'
      } else if (rtY.sig == 3.0) {
        real_T Bc_2;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // DiscreteIntegrator: '<S5>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'

        // '<S1>:59:37' elseif sig == 3
        // '<S1>:59:38' [u, yhat(1:no)] = mpc3(r_, y__, [0;0;0], [0;0], u0, umax, iRST); 
        // Simulink Function 'mpc3': '<S1>:936'
        if (rtU.iRST && (rtDW.DiscreteTimeIntegrator_PrevRese <= 0)) {
          rtDW.DiscreteTimeIntegrator_DSTATE[0] =
            rtP.DiscreteTimeIntegrator_IC_c[0];
          rtDW.DiscreteTimeIntegrator_DSTATE[1] =
            rtP.DiscreteTimeIntegrator_IC_c[1];
        }

        // Delay: '<S252>/MemoryX' incorporates:
        //   Constant: '<S252>/X0'
        //   DataTypeConversion: '<S252>/DataTypeConversionReset'

        rtDW.icLoad = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE) ==
                        POS_ZCSIG) || rtDW.icLoad);
        rtPrevZCX.MemoryX_Reset_ZCE = 0U;
        if (rtDW.icLoad) {
          for (k = 0; k < 5; k++) {
            rtDW.MemoryX_DSTATE[k] = rtP.X0_Value_a[k];
          }
        }

        // SignalConversion generated from: '<S251>/ SFunction ' incorporates:
        //   Constant: '<S5>/Constant'
        //   MATLAB Function: '<S250>/optimizer'

        rtb_TmpSignalConversionAtSFu_o4[0] = Sum2_c[0];
        rtb_TmpSignalConversionAtSFu_o4[1] = Sum2_c[1];
        rtb_TmpSignalConversionAtSFu_o4[2] = Sum2_c[2];
        rtb_TmpSignalConversionAtSFu_o4[3] = rtP.Constant_Value_e[0];
        rtb_TmpSignalConversionAtSFu_o4[4] = rtP.Constant_Value_e[1];

        // MATLAB Function: '<S250>/optimizer' incorporates:
        //   Constant: '<S229>/Constant1'
        //   Delay: '<S252>/MemoryX'
        //   DiscreteIntegrator: '<S5>/Discrete-Time Integrator'
        //   Gain: '<S5>/Gain2'
        //   Inport: '<Root>/umax'
        //   SignalConversion generated from: '<S251>/ SFunction '
        //   Sum: '<S229>/Sum2'
        //   UnitDelay: '<S230>/last_mv'

        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S251>:1'
        // '<S251>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S251>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S251>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S251>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S251>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S251>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S251>:1:252' if isa(u,'double')
        // '<S251>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S251>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S251>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S251>:1:252' if isa(u,'double')
        // '<S251>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S251>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S251>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S251>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S251>:1:252' if isa(u,'double')
        // '<S251>:1:253' y = u;
        //    md:         measured disturbance
        // '<S251>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S251>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S251>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S251>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S251>:1:252' if isa(u,'double')
        // '<S251>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S251>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S251>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S251>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S251>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S251>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S251>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S251>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S251>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S251>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S251>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S251>:1:58' uwt = convertDataType(uwt0,isDouble);
        //    duwt:       run-time DMV weights
        // '<S251>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S251>:1:62' ewt = convertDataType(ewt0,isDouble);
        //    iA:         logical previous active set (for warm start)
        //  Outputs (in BlockDataType except iAout)
        //    xk1:        x[k+1|k] from built-in KF
        //    u:          optimal MV
        //    cost:       optimal cost
        //    useq:       optimal MV sequence
        //    xseq:       optimal state sequence
        //    yseq:       optimal OV sequence
        //    status:     QP exit flag
        //    xest:       x[k|k] from built-in KF
        //    iAout:      logical current active set
        //  Parameters (constant)
        //    dimensions (int32):
        //        nx, nxp, nup, nu, ny, degrees, p, nxQP, enable_value, Mrows, nCC, nv 
        //        myindex, mvindex, mdindex, unindex, nxid, m, Ndis, numdis, maxdis 
        //    MPC constants (BlockDataType):
        //        Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, Mx, Mu1, Mv, utarget
        //        H, Linv, Ac, Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1
        //        A, C, B, D, Cid, Did, Ecc, Fcc, Scc, Gcc
        //        RYscale, RMDscale, xoff, Uscale, Yscale
        //        uoff, voff, yoff, myoff, RMVscale, Mdis, Vdis
        //    configurations (logical):
        //        isQP, CustomSolver, CustomSolverCodeGen, UseSuboptimalSolution, UseActiveSetSolver 
        //        openloopflag, no_umin, no_umax, no_ymin, no_ymax, switch_inport, no_switch 
        //        return_cost, return_mvseq, return_xseq, return_ovseq, isLTV
        //        no_ywt, no_uwt, no_duwt, no_rhoeps, no_md, no_ref, no_uref, no_mv 
        //        CustomEstimation, no_cc, isHyb, isDouble
        //    ASOptions
        //    IPOptions
        //    MIQPOptions
        //  Parameters
        // '<S251>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S251>:1:96' isAdaptive = false;
        // '<S251>:1:97' isLTV = false;
        // '<S251>:1:98' ZERO = zeros('like',ref);
        // '<S251>:1:99' ONE = ones('like',ref);
        // '<S251>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S251>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S251>:1:119' if isSimulation
        // '<S251>:1:126' else
        //  When doing code generation, use M code directly
        // '<S251>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        k_0 = 0;
        for (k = 0; k < 20; k++) {
          for (ii = 0; ii < 5; ii++) {
            rseq_0[ii + k_0] = rtb_TmpSignalConversionAtSFu_o4[ii];
          }

          k_0 += 5;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S251>:1:133' old_u = old_u - uoff;
        umin_scale1_idx_0 = rtDW.last_mv_DSTATE[0];
        umin_scale1_idx_1 = rtDW.last_mv_DSTATE[1];
        umin_scale1_idx_2 = rtDW.last_mv_DSTATE[2];

        // '<S251>:1:134' if no_mv
        // '<S251>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S251>:1:143' xk = xk - xoff;
        rtb_xest[0] = rtDW.MemoryX_DSTATE[0];
        rtb_xest[2] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE[0];
        rtb_xest[1] = rtDW.MemoryX_DSTATE[1];
        rtb_xest[3] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE[1];
        rtb_xest[4] = rtP.Constant1_Value_h[0] + rtDW.MemoryX_DSTATE[2];
        rtb_xest[5] = rtP.Constant1_Value_h[1] + rtDW.MemoryX_DSTATE[3];
        rtb_xest[6] = rtP.Constant1_Value_h[2] + rtDW.MemoryX_DSTATE[4];

        //  Remove offset
        // '<S251>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S251>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S251>:1:162' if no_uref
        //  no external utarget
        // '<S251>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S251>:1:173' if ~no_cc
        // '<S251>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S251>:1:183' if isSimulation
        // '<S251>:1:214' else
        //  When doing code generation, use M code directly
        // '<S251>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S251>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S251>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S251>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S251>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S251>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S251>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S251>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S251>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S251>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S251>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S251>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S251>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        umax_incr_flag[0] = false;
        umax_incr[0] = 0.0;
        umax_incr_flag[1] = false;
        umax_incr[1] = 0.0;
        umax_incr_flag[2] = false;
        umax_incr[2] = 0.0;
        for (k = 0; k < 126; k++) {
          int16_T b_Mlim;
          int16_T b_Mrows_0;
          b_Mlim = b_Mlim_2[k];
          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 7; k_0++) {
            Bc_2 += b_a_1[126 * k_0 + k] * rtb_xest[k_0];
          }

          Bc_2 = -(((a_1[k + 126] * umin_scale1_idx_1 + a_1[k] *
                     umin_scale1_idx_0) + a_1[k + 252] * umin_scale1_idx_2) + (
                    static_cast<real_T>(b_Mlim) + Bc_2));
          b_Mrows_0 = b_Mrows_3[k];
          if ((b_Mrows_0 > 100) && (b_Mrows_0 > 200) && (b_Mrows_0 <= 260)) {
            ii = (static_cast<int32_T>(b_Mrows_0) - div_nde_s32_floor(
                   static_cast<int32_T>(b_Mrows_0) - 201, static_cast<int32_T>
                   (nu)) * static_cast<int32_T>(nu)) - 201;
            rstP2 = umax_incr_flag[ii];
            if (!umax_incr_flag[ii]) {
              dwt = -rtU.umax[ii] - (-static_cast<real_T>(b_Mlim));
              rstP2 = true;
            } else {
              dwt = umax_incr[ii];
            }

            umax_incr[ii] = dwt;
            umax_incr_flag[ii] = rstP2;
            Bc_2 += dwt;
          }

          Bc_1[k] = Bc_2;
        }

        rtb_Sum2_f[0] = 0.0;
        rtb_Sum2_f[1] = 0.0;
        rtb_Sum2_f[2] = 0.0;
        rtb_Sum2_f[3] = 0.0;
        for (k = 0; k < 3; k++) {
          dwt = 0.0;
          for (k_0 = 0; k_0 < 7; k_0++) {
            dwt += b_Kx_1[7 * k + k_0] * rtb_xest[k_0];
          }

          Bc_2 = 0.0;
          for (k_0 = 0; k_0 < 100; k_0++) {
            Bc_2 += b_Kr_1[100 * k + k_0] * rseq_0[k_0];
          }

          rtb_Sum2_f[k] = ((b_Ku1_1[3 * k + 1] * umin_scale1_idx_1 + b_Ku1_1[3 *
                            k] * umin_scale1_idx_0) + b_Ku1_1[3 * k + 2] *
                           umin_scale1_idx_2) + (dwt + Bc_2);
        }

        // Update for Memory: '<S230>/Memory' incorporates:
        //   MATLAB Function: '<S250>/optimizer'

        qpkwik_f(b_Linv_1, b_Hinv_1, rtb_Sum2_f, b_Ac_1, Bc_1,
                 rtDW.Memory_PreviousInput, 520, 1.0E-6,
                 rtb_TmpSignalConversionAtSFu_ia, a__1_1, &k);

        // MATLAB Function: '<S250>/optimizer' incorporates:
        //   UnitDelay: '<S230>/last_mv'

        if ((k < 0) || (k == 0)) {
          rtb_TmpSignalConversionAtSFu_ia[0] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[1] = 0.0;
          rtb_TmpSignalConversionAtSFu_ia[2] = 0.0;
        }

        umax_incr[0] = rtDW.last_mv_DSTATE[0] + rtb_TmpSignalConversionAtSFu_ia
          [0];
        umax_incr[1] = rtDW.last_mv_DSTATE[1] + rtb_TmpSignalConversionAtSFu_ia
          [1];
        umax_incr[2] = rtDW.last_mv_DSTATE[2] + rtb_TmpSignalConversionAtSFu_ia
          [2];

        // Delay: '<S252>/MemoryP' incorporates:
        //   Constant: '<S252>/P0'
        //   DataTypeConversion: '<S252>/DataTypeConversionReset'

        // '<S251>:1:231' if return_xseq || return_ovseq
        // '<S251>:1:233' else
        // '<S251>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S251>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S251>:1:238' if CustomEstimation
        // '<S251>:1:239' xk1 = xk;
        // '<S251>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S251>:1:247' xest = xest + xoff;
        rtDW.icLoad_e = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE) ==
                          POS_ZCSIG) || rtDW.icLoad_e);
        rtPrevZCX.MemoryP_Reset_ZCE = 0U;
        if (rtDW.icLoad_e) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtP.P0_Value_m[0], 25U *
                            sizeof(real_T));
        }

        // MATLAB Function: '<S229>/MATLAB Function' incorporates:
        //   BusCreator: '<S5>/Bus Creator1'
        //   Constant: '<S5>/Constant12'
        //   Constant: '<S5>/Constant13'
        //   Constant: '<S5>/Constant3'
        //   Constant: '<S5>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc3/State Estimator OD (KF)/MATLAB Function': '<S253>:1' 
        // '<S253>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod3, Bod3, Cod3, Dod3, Dmn1, 3, 3, 2); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&rtb_B[0], 0, 15U * sizeof(real_T));

        //  n_states  x n_inputs
        // 'stateEst_:9' C = zeros(no,ns);
        //  n_outputs x n_states
        // 'stateEst_:10' D = zeros(no,ni);
        //  n_outputs x n_inputs
        // 'stateEst_:11' Q = zeros(ns,ns);
        //  n_states  x n_states
        // 'stateEst_:12' G = eye(ns);
        //  n_states  x n_states
        // 'stateEst_:13' R = zeros(no,no);
        //  n_outputs x n_outputs
        // 'stateEst_:14' N = zeros(ns,no);
        //  n_states  x n_outputs
        // 'stateEst_:15' H = zeros(no,ns);
        //  n_outputs x n_states
        //  combine plant and output disturbance model
        //  (force the outputs to fit in preallocated memory)
        // 'stateEst_:19' A(1:ns, 1:ns) = blkdiag(Ap, Aod);
        (void)std::memset(&rtb_A[0], 0, 25U * sizeof(real_T));
        rtb_A[0] = rtP.Constant3_Value_g[0];
        rtb_A[1] = rtP.Constant3_Value_g[1];
        rtb_A[5] = rtP.Constant3_Value_g[2];
        rtb_A[6] = rtP.Constant3_Value_g[3];

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        k_0 = 0;
        k = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          rtb_A[k_0 + 12] = rtP.Aod3[k];
          rtb_A[k_0 + 13] = rtP.Aod3[k + 1];
          rtb_A[k_0 + 14] = rtP.Aod3[k + 2];
          rtb_B[k_0] = rtP.Constant4_Value_f[ii];
          rtb_B[k_0 + 1] = rtP.Constant4_Value_f[ii + 1];
          k_0 += 5;
          k += 3;
          ii += 2;
        }

        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        for (k = 0; k < 6; k++) {
          rtb_C[k] = rtP.Constant12_Value_f[k];
        }

        (void)std::memcpy(&rtb_C[6], &rtP.Cod3[0], 9U * sizeof(real_T));

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est_0[0], 0, 45U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y_0[0], 0, 30U * sizeof(real_T));
        k_0 = 0;
        k = 0;
        ii = 0;
        for (i = 0; i < 3; i++) {
          y_0[k_0] = rtP.Constant4_Value_f[k];
          y_0[k_0 + 1] = rtP.Constant4_Value_f[k + 1];
          y_0[k_0 + 17] = rtP.Bod3[ii];
          y_0[k_0 + 18] = rtP.Bod3[ii + 1];
          y_0[k_0 + 19] = rtP.Bod3[ii + 2];
          k_0 += 5;
          k += 2;
          ii += 3;
        }

        (void)std::memcpy(&B_est_0[0], &y_0[0], 30U * sizeof(real_T));

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        for (k_0 = 0; k_0 < 9; k_0++) {
          D_est[k_0] = rtP.Constant13_Value_a[k_0];
          D_est[k_0 + 9] = rtP.Dod3[k_0];
          D_est[k_0 + 18] = rtP.Dmn1[k_0];
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (k_0 = 0; k_0 < 5; k_0++) {
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            int32_T rtb_Q_tmp;
            rtb_Q_tmp = k + k_0;
            rtb_Q[rtb_Q_tmp] = 0.0;
            i = 0;
            for (int32_T rtb_y_tmp{0}; rtb_y_tmp < 9; rtb_y_tmp++) {
              rtb_Q[rtb_Q_tmp] += B_est_0[i + k_0] * B_est_0[i + ii];
              i += 5;
            }

            k += 5;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        k_0 = 0;
        for (k = 0; k < 9; k++) {
          rtb_R_tmp[k] = D_est[k_0];
          rtb_R_tmp[k + 9] = D_est[k_0 + 1];
          rtb_R_tmp[k + 18] = D_est[k_0 + 2];
          k_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (k_0 = 0; k_0 < 3; k_0++) {
          for (k = 0; k < 3; k++) {
            i = 3 * k_0 + k;
            rtb_R[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_R[i] += D_est[3 * ii + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }

          for (k = 0; k < 5; k++) {
            i = 5 * k_0 + k;
            rtb_N[i] = 0.0;
            for (ii = 0; ii < 9; ii++) {
              rtb_N[i] += B_est_0[5 * ii + k] * rtb_R_tmp[9 * k_0 + ii];
            }
          }
        }

        // End of MATLAB Function: '<S229>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S252>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S252>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S252>/ScalarExpansionQ'
        ScalarExpansionQ(rtb_Q, rtb_Z_e);

        // End of Outputs for SubSystem: '<S252>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S252>/ReducedQRN'
        ReducedQRN(rtP.G_Value_h, rtP.H_Value_oa, rtb_Z_e, rtb_y, rtb_N,
                   rtb_Product, rtb_R, rtb_Product2);

        // End of Outputs for SubSystem: '<S252>/ReducedQRN'

        // Outputs for Atomic SubSystem: '<S252>/CalculatePL'
        CalculatePL(rtb_A, rtb_C, rtb_Product, rtb_R, rtb_Product2,
                    rtP.Constant1_Value_n != 0.0, rtDW.MemoryP_DSTATE, rtb_N,
                    rtb_L, rtb_Z_e, rtb_Q);

        // End of Outputs for SubSystem: '<S252>/CalculatePL'

        // MATLAB Function: '<S293>/SqrtUsedFcn' incorporates:
        //   Constant: '<S252>/G'
        //   Constant: '<S252>/H'
        //   Constant: '<S293>/isSqrtUsed'
        //   Constant: '<S5>/Constant1'
        //   DataTypeConversion: '<S252>/DataTypeConversionEnable'
        //   Delay: '<S252>/MemoryP'

        SqrtUsedFcn(rtb_Z_e, rtP.isSqrtUsed_Value_p, rtb_Product);

        // Gain: '<S230>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_g[0] * umax_incr[0];

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = umin_scale1_idx_0 - rtU.u0[0];

        // Sum: '<S229>/Sum6'
        rtb_Sum6[0] = y__m[0];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[0] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Gain: '<S230>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_g[1] * umax_incr[1];

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = umin_scale1_idx_0 - rtU.u0[1];

        // Sum: '<S229>/Sum6'
        rtb_Sum6[1] = y__m[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[1] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Gain: '<S230>/umin_scale1'
        umin_scale1_idx_0 = rtP.umin_scale1_Gain_g[2] * umax_incr[2];

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = umin_scale1_idx_0 - rtU.u0[2];

        // Sum: '<S229>/Sum6'
        rtb_Sum6[2] = y__m[2];

        // Outputs for Enabled SubSystem: '<S271>/MeasurementUpdate'
        MeasurementUpdate(rtP.Constant1_Value_n != 0.0, rtb_L, rtb_Sum6, rtb_C,
                          rtDW.MemoryX_DSTATE, rtP.Constant13_Value_a, rtb_Sum1,
                          rtDW.Product3, &rtDW.MeasurementUpdate_c,
                          &rtP.MeasurementUpdate_c);

        // End of Outputs for SubSystem: '<S271>/MeasurementUpdate'
        for (k_0 = 0; k_0 < 3; k_0++) {
          // Product: '<S255>/Product'
          rtb_C_0[k_0] = 0.0;
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            rtb_C_0[k_0] += rtb_C[k + k_0] * rtDW.MemoryX_DSTATE[ii];
            k += 3;
          }

          // Product: '<S255>/Product1' incorporates:
          //   Product: '<S255>/Product'

          tmp[k_0] = 0.0;
          tmp[k_0] += rtP.Constant13_Value_a[k_0] * rtb_Sum1[0];
          tmp[k_0] += rtP.Constant13_Value_a[k_0 + 3] * rtb_Sum1[1];
          tmp[k_0] += rtP.Constant13_Value_a[k_0 + 6] * rtb_Sum1[2];

          // Sum: '<S255>/Add1' incorporates:
          //   Product: '<S255>/Product'
          //   Product: '<S255>/Product1'

          rtb_Sum6[k_0] = rtb_C_0[k_0] + tmp[k_0];
        }

        // Update for DiscreteIntegrator: '<S5>/Discrete-Time Integrator' incorporates:
        //   Constant: '<S5>/Constant1'
        //   Constant: '<S5>/Constant13'
        //   DataTypeConversion: '<S252>/DataTypeConversionEnable'
        //   Delay: '<S252>/MemoryX'
        //   Inport: '<Root>/iRST'
        //   Product: '<S255>/Product'
        //   Product: '<S255>/Product1'
        //   Sum: '<S5>/Sum'

        rtDW.DiscreteTimeIntegrator_DSTATE[0] += (y__m[0] - Sum2_c[0]) *
          rtP.DiscreteTimeIntegrator_gainva_k;
        rtDW.DiscreteTimeIntegrator_DSTATE[1] += (y__m[1] - Sum2_c[1]) *
          rtP.DiscreteTimeIntegrator_gainva_k;
        rtDW.DiscreteTimeIntegrator_PrevRese = static_cast<int8_T>(rtU.iRST ? 1 :
          0);

        // Update for UnitDelay: '<S230>/last_mv'
        rtDW.last_mv_DSTATE[0] = umax_incr[0];
        rtDW.last_mv_DSTATE[1] = umax_incr[1];
        rtDW.last_mv_DSTATE[2] = umax_incr[2];

        // Update for Delay: '<S252>/MemoryX'
        rtDW.icLoad = false;
        for (k_0 = 0; k_0 < 5; k_0++) {
          // Product: '<S271>/B[k]*u[k]'
          rtb_TmpSignalConversionAtSFu_o4[k_0] = 0.0;
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0] * rtb_Sum1[0];
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0 + 5] * rtb_Sum1[1];
          rtb_TmpSignalConversionAtSFu_o4[k_0] += rtb_B[k_0 + 10] * rtb_Sum1[2];

          // Product: '<S271>/A[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S252>/MemoryX'
          //   Product: '<S271>/B[k]*u[k]'

          rtb_A_0[k_0] = 0.0;
          k = 0;
          for (ii = 0; ii < 5; ii++) {
            rtb_A_0[k_0] += rtb_A[k + k_0] * rtDW.MemoryX_DSTATE[ii];
            k += 5;
          }

          // End of Product: '<S271>/A[k]*xhat[k|k-1]'
        }

        // End of Outputs for SubSystem: '<S1>/mpc3'
        for (k_0 = 0; k_0 <= 2; k_0 += 2) {
          __m128d tmp_1;
          __m128d tmp_2;

          // Outputs for Function Call SubSystem: '<S1>/mpc3'
          tmp_3 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFu_o4[k_0]);
          tmp_1 = _mm_loadu_pd(&rtb_A_0[k_0]);
          tmp_2 = _mm_loadu_pd(&rtDW.Product3[k_0]);
          (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE[k_0], _mm_add_pd(_mm_add_pd
            (tmp_3, tmp_1), tmp_2));

          // End of Outputs for SubSystem: '<S1>/mpc3'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        for (k_0 = 4; k_0 < 5; k_0++) {
          // Update for Delay: '<S252>/MemoryX' incorporates:
          //   Sum: '<S271>/Add'

          rtDW.MemoryX_DSTATE[k_0] = (rtb_TmpSignalConversionAtSFu_o4[k_0] +
            rtb_A_0[k_0]) + rtDW.Product3[k_0];
        }

        // Update for Delay: '<S252>/MemoryP' incorporates:
        //   Delay: '<S252>/MemoryX'
        //   Product: '<S271>/B[k]*u[k]'
        //   Sum: '<S271>/Add'

        rtDW.icLoad_e = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtb_Q[0], 25U * sizeof
                          (real_T));
        rtY.yhat[0] = rtb_Sum6[0];
        rtY.yhat[1] = rtb_Sum6[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S229>/Sum3'

        rtY.u[2] = umin_scale1_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        rtY.yhat[2] = rtb_Sum6[2];

        // End of Outputs for SubSystem: '<S1>/mpc3'
      } else {
        // no actions
      }
    }
  }
}

// Model initialize function
void SupervisoryController::initialize()
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  // Initialize DataMapInfo substructure containing ModelMap for C API
  SupervisoryController_InitializeDataMapInfo((&rtM), &rtP);

  {
    static const real_T tmp_2[24]{ 0.025, 0.028867513459481294,
      -0.014433756729740647, -0.014433756729740647, 0.025, -0.014433756729740647,
      0.028867513459481294, -0.014433756729740647, 0.025, -0.014433756729740647,
      -0.014433756729740647, 0.028867513459481294, 0.025, 0.028867513459481294,
      -0.014433756729740647, -0.014433756729740647, 0.025, -0.014433756729740647,
      0.028867513459481294, -0.014433756729740647, 0.025, -0.014433756729740647,
      -0.014433756729740647, 0.028867513459481294 };

    static const real_T tmp_0[12]{ 0.025, 0.028867513459481294,
      -0.014433756729740647, -0.014433756729740647, 0.025, -0.014433756729740647,
      0.028867513459481294, -0.014433756729740647, 0.025, -0.014433756729740647,
      -0.014433756729740647, 0.028867513459481294 };

    static const int8_T tmp_3[576]{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

    static const int8_T tmp[144]{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1 };

    static const int8_T tmp_1[24]{ 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, 1,
      1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1 };

    real_T Product1_j[144];
    real_T Sum_h[12];
    real_T Sum2_c[3];
    int32_T i;
    int32_T t;
    uint32_T tseed;
    rtPrevZCX.MemoryX_Reset_ZCE_o = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryP_Reset_ZCE_b = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryX_Reset_ZCE_l = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryP_Reset_ZCE_n = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryX_Reset_ZCE_g = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryP_Reset_ZCE_i = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryX_Reset_ZCE = UNINITIALIZED_ZCSIG;
    rtPrevZCX.MemoryP_Reset_ZCE = UNINITIALIZED_ZCSIG;
    rtPrevZCX.SupervisoryController_Trig_ZCE = UNINITIALIZED_ZCSIG;
    rtPrevZCX.paramEst2.Delay1_Reset_ZCE = POS_ZCSIG;
    rtPrevZCX.paramEst2.Delay_Reset_ZCE = POS_ZCSIG;
    rtPrevZCX.paramEst1_o.Delay1_Reset_ZCE = POS_ZCSIG;
    rtPrevZCX.paramEst1_o.Delay_Reset_ZCE = POS_ZCSIG;
    for (i = 0; i < 144; i++) {
      rtDW.P0_1[i] = static_cast<real_T>(tmp[i]);
    }

    (void)std::memcpy(&rtDW.theta0_1[0], &tmp_0[0], 12U * sizeof(real_T));
    for (i = 0; i < 24; i++) {
      rtDW.thetaSgn[i] = static_cast<real_T>(tmp_1[i]);

      // SystemInitialize for Outport: '<Root>/theta'
      rtY.theta[i] = tmp_2[i];
    }

    for (i = 0; i < 576; i++) {
      // SystemInitialize for Outport: '<Root>/P'
      rtY.P_c[i] = static_cast<real_T>(tmp_3[i]);
    }

    rtDW.waypt = 1U;
    rtDW.trajSize = 2400U;
    (void)std::memcpy(&rtDW.P0_2[0], &rtDW.P0_1[0], 144U * sizeof(real_T));
    (void)std::memcpy(&rtDW.theta0_2[0], &rtDW.theta0_1[0], 12U * sizeof(real_T));

    // SystemInitialize for Outport: '<Root>/currEv'
    rtY.currEv.preT = 0.0;
    rtY.currEv.moveT = 0.0;
    rtY.currEv.postT = 0.0;
    for (i = 0; i < 6; i++) {
      // SystemInitialize for Outport: '<Root>/currEv' incorporates:
      //   Outport: '<Root>/ywt'

      rtY.currEv.r[i] = 0.0;
    }

    // SystemInitialize for Outport: '<Root>/sig'
    rtY.sig = 1.0;

    // 'gainSchSig_:4' sigPrev = 1;
    rtDW.sigPrev = 1.0;

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/paramEst1'

    paramEst1_Init(Sum_h, Product1_j, Sum2_c, &rtDW.paramEst1_o,
                   &rtP.paramEst1_o);

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/paramEst2'

    paramEst1_Init(Sum_h, Product1_j, Sum2_c, &rtDW.paramEst2, &rtP.paramEst2);

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/ampc'

    // InitializeConditions for RandomNumber: '<S2>/excitation'
    i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[0]) >> 16UL));
    t = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[0]) & 32768U));
    tseed = ((((static_cast<uint32_T>(rtP.excitation_Seed[0]) -
                (static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t)) <<
              16UL) + static_cast<uint32_T>(t)) + static_cast<uint32_T>(i);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    } else {
      // no actions
    }

    (void)rt_nrand_Upu32_Yd_f_pw_snf(&tseed);
    i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[1]) >> 16UL));
    t = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[1]) & 32768U));
    tseed = ((((static_cast<uint32_T>(rtP.excitation_Seed[1]) -
                (static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t)) <<
              16UL) + static_cast<uint32_T>(t)) + static_cast<uint32_T>(i);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    } else {
      // no actions
    }

    (void)rt_nrand_Upu32_Yd_f_pw_snf(&tseed);
    i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[2]) >> 16UL));
    t = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
      (rtP.excitation_Seed[2]) & 32768U));
    tseed = ((((static_cast<uint32_T>(rtP.excitation_Seed[2]) -
                (static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t)) <<
              16UL) + static_cast<uint32_T>(t)) + static_cast<uint32_T>(i);
    if (tseed < 1U) {
      tseed = 1144108930U;
    } else if (tseed > 2147483646U) {
      tseed = 2147483646U;
    } else {
      // no actions
    }

    (void)rt_nrand_Upu32_Yd_f_pw_snf(&tseed);

    // End of InitializeConditions for RandomNumber: '<S2>/excitation'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc1'

    // InitializeConditions for Memory: '<S90>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_d[0],
                      &rtP.Memory_InitialCondition_f[0], 166U * sizeof(boolean_T));

    // InitializeConditions for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE_j = rtP.DiscreteTimeIntegrator_IC;
    rtDW.DiscreteTimeIntegrator_PrevRe_b = 2;

    // InitializeConditions for UnitDelay: '<S90>/last_mv'
    rtDW.last_mv_DSTATE_n[0] = rtP.last_mv_InitialCondition_f[0];
    rtDW.last_mv_DSTATE_n[1] = rtP.last_mv_InitialCondition_f[1];
    rtDW.last_mv_DSTATE_n[2] = rtP.last_mv_InitialCondition_f[2];

    // InitializeConditions for Delay: '<S112>/MemoryX'
    rtDW.icLoad_n = true;

    // InitializeConditions for Delay: '<S112>/MemoryP'
    rtDW.icLoad_h = true;

    // SystemInitialize for Enabled SubSystem: '<S131>/MeasurementUpdate'
    // SystemInitialize for Product: '<S155>/Product3' incorporates:
    //   Outport: '<S155>/L*(y[k]-yhat[k|k-1])'

    rtDW.Product3_c[0] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3_c[1] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3_c[2] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3_c[3] = rtP.Lykyhatkk1_Y0_c;

    // End of SystemInitialize for SubSystem: '<S131>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc2'

    // InitializeConditions for Memory: '<S160>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_c[0],
                      &rtP.Memory_InitialCondition_j[0], 206U * sizeof(boolean_T));

    // InitializeConditions for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE_m[0] = rtP.DiscreteTimeIntegrator_IC_n[0];
    rtDW.DiscreteTimeIntegrator_DSTATE_m[1] = rtP.DiscreteTimeIntegrator_IC_n[1];
    rtDW.DiscreteTimeIntegrator_PrevRe_f = 2;

    // InitializeConditions for UnitDelay: '<S160>/last_mv'
    rtDW.last_mv_DSTATE_i[0] = rtP.last_mv_InitialCondition_b[0];
    rtDW.last_mv_DSTATE_i[1] = rtP.last_mv_InitialCondition_b[1];
    rtDW.last_mv_DSTATE_i[2] = rtP.last_mv_InitialCondition_b[2];

    // InitializeConditions for Delay: '<S182>/MemoryX'
    rtDW.icLoad_a = true;

    // InitializeConditions for Delay: '<S182>/MemoryP'
    rtDW.icLoad_p = true;

    // SystemInitialize for Enabled SubSystem: '<S201>/MeasurementUpdate'
    MeasurementUpdate_Init(rtDW.Product3_a, &rtP.MeasurementUpdate_j);

    // End of SystemInitialize for SubSystem: '<S201>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/wtMod'

    for (i = 0; i < 6; i++) {
      // InitializeConditions for Delay: '<S8>/Delay'
      rtDW.Delay_DSTATE[i] = rtP.ywt0[i];
    }

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc3'

    // InitializeConditions for DiscreteIntegrator: '<S5>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE[0] = rtP.DiscreteTimeIntegrator_IC_c[0];
    rtDW.DiscreteTimeIntegrator_DSTATE[1] = rtP.DiscreteTimeIntegrator_IC_c[1];
    rtDW.DiscreteTimeIntegrator_PrevRese = 2;

    // InitializeConditions for Memory: '<S230>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput[0],
                      &rtP.Memory_InitialCondition_b[0], 126U * sizeof(boolean_T));

    // InitializeConditions for UnitDelay: '<S230>/last_mv'
    rtDW.last_mv_DSTATE[0] = rtP.last_mv_InitialCondition_i[0];
    rtDW.last_mv_DSTATE[1] = rtP.last_mv_InitialCondition_i[1];
    rtDW.last_mv_DSTATE[2] = rtP.last_mv_InitialCondition_i[2];

    // InitializeConditions for Delay: '<S252>/MemoryX'
    rtDW.icLoad = true;

    // InitializeConditions for Delay: '<S252>/MemoryP'
    rtDW.icLoad_e = true;

    // SystemInitialize for Enabled SubSystem: '<S271>/MeasurementUpdate'
    MeasurementUpdate_Init(rtDW.Product3, &rtP.MeasurementUpdate_c);

    // End of SystemInitialize for SubSystem: '<S271>/MeasurementUpdate'
  }
}

// Model terminate function
void SupervisoryController::terminate()
{
  // (no terminate code required)
}

// Constructor
SupervisoryController::SupervisoryController() :
  rtU(),
  rtY(),
  rtDW(),
  rtPrevZCX(),
  rtM()
{
  // Currently there is no constructor body generated.
}

// Destructor
SupervisoryController::~SupervisoryController()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
SupervisoryController::RT_MODEL * SupervisoryController::getRTM()
{
  return (&rtM);
}

//
// File trailer for generated code.
//
// [EOF]
//
