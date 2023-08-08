//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2467
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon Aug  7 23:00:45 2023
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
  static const real_T b_a_0[1442]{ -0.0, -1.0001226651251065,
    0.0016804698402870472, -0.025, -0.0, -0.0, -1.0002454535095886,
    0.0033601438631711473, -0.050003066628127667, 4.201174600717618E-5, -0.0,
    -1.0003683651173199, 0.0050390227502789741, -0.075009202965867383,
    0.00012601534258645488, -0.0, -1.0004913999122134, 0.0067171071827700837,
    -0.10001841209380039, 0.00025199091134342927, -0.0, -1.0006145578582215,
    0.0083943978413372582, -0.12503069709160572, 0.00041991859091268134, -0.0,
    -1.0007378389193362, 0.010070895406206851, -0.15004606103806126,
    0.00062977853694611279, -0.0, -1.0008612430595885, 0.011746600557139133,
    -0.17506450701104467, 0.000881550922101284, -0.0, -1.000984770243049,
    0.013421513973428631, -0.20008603808753436, 0.0011752159360297622, -0.0,
    -1.0011084204338279, 0.015095636333904478, -0.22511065734361058,
    0.0015107537853654778, -0.0, -1.0012321935960742, 0.016768968316930755,
    -0.25013836785445626, 0.0018881446937130896, -0.0, -1.0013560896939764,
    0.01844151060040683, -0.27516917269435814, 0.0023073689016363582, -0.0,
    -1.0014801086917624, 0.020113263861767709, -0.30020307493670756,
    0.0027684066666465289, -0.0, -1.0016042505536991, 0.021784228777984374,
    -0.32524007765400159, 0.0032712382631907219, -0.0, -1.0017285152440929,
    0.023454406025564121, -0.3502801839178441, 0.0038158439826403317, -0.0,
    -1.0018529027272891, 0.025123796280550909, -0.37532339679894644,
    0.0044022041332794351, -0.0, -1.0019774129676724, 0.026792400218525705,
    -0.40036971936712867, 0.0050302990402932082, -0.0, -1.0021020459296663,
    0.028460218514606821, -0.42541915469132052, 0.005700109045756351, -0.0,
    -1.0022268015777338, 0.030127251843450248, -0.4504717058395622,
    0.006411614508621522, -0.0, -1.0023516798763767, 0.031793500879250015,
    -0.47552737587900556, 0.0071647958047077786, -0.0, -1.002476680790136,
    0.03345896629573851, -0.500586167875915, 0.0079596333266890289, 0.0,
    1.0001226651251065, -0.0016804698402870472, 0.025, 0.0, 0.0,
    1.0002454535095886, -0.0033601438631711473, 0.050003066628127667,
    -4.201174600717618E-5, 0.0, 1.0003683651173199, -0.0050390227502789741,
    0.075009202965867383, -0.00012601534258645488, 0.0, 1.0004913999122134,
    -0.0067171071827700837, 0.10001841209380039, -0.00025199091134342927, 0.0,
    1.0006145578582215, -0.0083943978413372582, 0.12503069709160572,
    -0.00041991859091268134, 0.0, 1.0007378389193362, -0.010070895406206851,
    0.15004606103806126, -0.00062977853694611279, 0.0, 1.0008612430595885,
    -0.011746600557139133, 0.17506450701104467, -0.000881550922101284, 0.0,
    1.000984770243049, -0.013421513973428631, 0.20008603808753436,
    -0.0011752159360297622, 0.0, 1.0011084204338279, -0.015095636333904478,
    0.22511065734361058, -0.0015107537853654778, 0.0, 1.0012321935960742,
    -0.016768968316930755, 0.25013836785445626, -0.0018881446937130896, 0.0,
    1.0013560896939764, -0.01844151060040683, 0.27516917269435814,
    -0.0023073689016363582, 0.0, 1.0014801086917624, -0.020113263861767709,
    0.30020307493670756, -0.0027684066666465289, 0.0, 1.0016042505536991,
    -0.021784228777984374, 0.32524007765400159, -0.0032712382631907219, 0.0,
    1.0017285152440929, -0.023454406025564121, 0.3502801839178441,
    -0.0038158439826403317, 0.0, 1.0018529027272891, -0.025123796280550909,
    0.37532339679894644, -0.0044022041332794351, 0.0, 1.0019774129676724,
    -0.026792400218525705, 0.40036971936712867, -0.0050302990402932082, 0.0,
    1.0021020459296663, -0.028460218514606821, 0.42541915469132052,
    -0.005700109045756351, 0.0, 1.0022268015777338, -0.030127251843450248,
    0.4504717058395622, -0.006411614508621522, 0.0, 1.0023516798763767,
    -0.031793500879250015, 0.47552737587900556, -0.0071647958047077786, 0.0,
    1.002476680790136, -0.03345896629573851, 0.500586167875915,
    -0.0079596333266890289, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    6.4394278387945068E-5, -0.999403766481181, -0.0, -0.025, -0.0,
    0.00012875806168090968, -0.99880799666941367, 1.6098569596986267E-6,
    -0.049985094162029529, -0.0, 0.00019309137599829017, -0.99821269023697423,
    4.8288085017213693E-6, -0.074955294078764875, -0.0, 0.00025739424744158326,
    -0.99761784685637822, 9.6560929016786243E-6, -0.099910611334689231, -0.0,
    0.00032166670209439935, -0.99702346620038029, 1.6090949087718206E-5,
    -0.12485105750609868, -0.0, 0.00038590876602247561, -0.996429547941974,
    2.4132616640078191E-5, -0.1497766441611082, -0.0, 0.00045012046527368919,
    -0.99583609175439181, 3.3780335790640084E-5, -0.17468738285965754, -0.0,
    0.00051430182587807046, -0.9952430973111045, 4.5033347422482317E-5,
    -0.19958328515351734, -0.0, 0.00057845287384781629, -0.99465056428582155,
    5.7890893069434082E-5, -0.22446436258629496, -0.0, 0.00064257363517730291,
    -0.99405849235249066, 7.2352214915629483E-5, -0.2493306266934405, -0.0,
    0.00070666413584309932, -0.99346688118529747, 8.8416555795062062E-5,
    -0.27418208900225277, -0.0, 0.00077072440180398047, -0.99287573045866573,
    0.00010608315919113954, -0.29901876103188524, -0.0, 0.00083475445900094006,
    -0.9922850398472568, 0.00012535126923623904, -0.3238406542933519, -0.0,
    0.00089875433335720413, -0.99169480902596963, 0.00014622013071126255,
    -0.34864778028953336, -0.0, 0.00096272405077824379, -0.99110503766994063,
    0.00016868898904519265, -0.37344015051518265, -0.0, 0.0010266636371517885,
    -0.99051572545454325, 0.00019275709031464874, -0.39821777645693118, -0.0,
    0.0010905731183478392, -0.98992687205538832, 0.00021842368124344344,
    -0.42298066959329478, -0.0, 0.0011544525202186815, -0.98933847714832324,
    0.00024568800920213941, -0.4477288413946795, -0.0, 0.0012183018685988985,
    -0.98875054040943233, 0.00027454932220760648, -0.47246230332338762, -0.0,
    0.0012821211893053839, -0.98816306151503641, 0.00030500686892257894,
    -0.49718106683362345, 0.0, -6.4394278387945068E-5, 0.999403766481181, 0.0,
    0.025, 0.0, -0.00012875806168090968, 0.99880799666941367,
    -1.6098569596986267E-6, 0.049985094162029529, 0.0, -0.00019309137599829017,
    0.99821269023697423, -4.8288085017213693E-6, 0.074955294078764875, 0.0,
    -0.00025739424744158326, 0.99761784685637822, -9.6560929016786243E-6,
    0.099910611334689231, 0.0, -0.00032166670209439935, 0.99702346620038029,
    -1.6090949087718206E-5, 0.12485105750609868, 0.0, -0.00038590876602247561,
    0.996429547941974, -2.4132616640078191E-5, 0.1497766441611082, 0.0,
    -0.00045012046527368919, 0.99583609175439181, -3.3780335790640084E-5,
    0.17468738285965754, 0.0, -0.00051430182587807046, 0.9952430973111045,
    -4.5033347422482317E-5, 0.19958328515351734, 0.0, -0.00057845287384781629,
    0.99465056428582155, -5.7890893069434082E-5, 0.22446436258629496, 0.0,
    -0.00064257363517730291, 0.99405849235249066, -7.2352214915629483E-5,
    0.2493306266934405, 0.0, -0.00070666413584309932, 0.99346688118529747,
    -8.8416555795062062E-5, 0.27418208900225277, 0.0, -0.00077072440180398047,
    0.99287573045866573, -0.00010608315919113954, 0.29901876103188524, 0.0,
    -0.00083475445900094006, 0.9922850398472568, -0.00012535126923623904,
    0.3238406542933519, 0.0, -0.00089875433335720413, 0.99169480902596963,
    -0.00014622013071126255, 0.34864778028953336, 0.0, -0.00096272405077824379,
    0.99110503766994063, -0.00016868898904519265, 0.37344015051518265, 0.0,
    -0.0010266636371517885, 0.99051572545454325, -0.00019275709031464874,
    0.39821777645693118, 0.0, -0.0010905731183478392, 0.98992687205538832,
    -0.00021842368124344344, 0.42298066959329478, 0.0, -0.0011544525202186815,
    0.98933847714832324, -0.00024568800920213941, 0.4477288413946795, 0.0,
    -0.0012183018685988985, 0.98875054040943233, -0.00027454932220760648,
    0.47246230332338762, 0.0, -0.0012821211893053839, 0.98816306151503641,
    -0.00030500686892257894, 0.49718106683362345, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

  static const real_T b_a_1[882]{ -0.98979340248168846, 0.036210018833518431,
    -0.0, -0.9796800787081813, 0.072076616087004344, -0.0, -0.96965918005766349,
    0.10760292110357451, -0.0, -0.9597298656277844, 0.14279203475849142, -0.0,
    -0.9498913021654386, 0.17764702971811902, -0.0, -0.94014266399718538,
    0.21217095069652278, -0.0, -0.93048313296030061, 0.24636681470973576, -0.0,
    -0.920911898334457, 0.28023761132771158, -0.0, -0.91142815677402456,
    0.31378630292398574, -0.0, -0.90203111224098764, 0.34701582492306576, -0.0,
    -0.89271997593847163, 0.37992908604557096, -0.0, -0.883493966244875,
    0.41252896855114257, -0.0, -0.87435230864860014, 0.44481832847914388, -0.0,
    -0.86529423568337782, 0.4767999958871712, -0.0, -0.85631898686418007,
    0.50847677508739519, -0.0, -0.847425808623716, 0.539851444880752, -0.0,
    -0.83861395424950524, 0.57092675878900434, -0.0, -0.82988268382152419,
    0.60170544528469228, -0.0, -0.82123126415041847, 0.63219020801899084, -0.0,
    -0.81265896871627885, 0.66238372604749651, -0.0, 0.98979340248168846,
    -0.036210018833518431, 0.0, 0.9796800787081813, -0.072076616087004344, 0.0,
    0.96965918005766349, -0.10760292110357451, 0.0, 0.9597298656277844,
    -0.14279203475849142, 0.0, 0.9498913021654386, -0.17764702971811902, 0.0,
    0.94014266399718538, -0.21217095069652278, 0.0, 0.93048313296030061,
    -0.24636681470973576, 0.0, 0.920911898334457, -0.28023761132771158, 0.0,
    0.91142815677402456, -0.31378630292398574, 0.0, 0.90203111224098764,
    -0.34701582492306576, 0.0, 0.89271997593847163, -0.37992908604557096, 0.0,
    0.883493966244875, -0.41252896855114257, 0.0, 0.87435230864860014,
    -0.44481832847914388, 0.0, 0.86529423568337782, -0.4767999958871712, 0.0,
    0.85631898686418007, -0.50847677508739519, 0.0, 0.847425808623716,
    -0.539851444880752, 0.0, 0.83861395424950524, -0.57092675878900434, 0.0,
    0.82988268382152419, -0.60170544528469228, 0.0, 0.82123126415041847,
    -0.63219020801899084, 0.0, 0.81265896871627885, -0.66238372604749651, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00030104618687214896, -1.0007224384072255,
    -0.0, -0.00059923720380821669, -1.0014344978436069, -0.0,
    -0.000894599067828991, -1.0021362741967341, -0.0, -0.0011871575592765928,
    -1.0028278624813931, -0.0, -0.0014769382239674068, -1.0035093568475058, -0.0,
    -0.0017639663753254257, -1.0041808505879968, -0.0, -0.002048267096496192,
    -1.0048424361465897, -0.0, -0.0023298652424415065, -1.005494205125532, -0.0,
    -0.0026087854420150866, -1.006136248293249, -0.0, -0.0028850521000193404,
    -1.00676865559193, -0.0, -0.0031586893992434336, -1.0073915161450431, -0.0,
    -0.0034297213024828167, -1.0080049182647837, -0.0, -0.0036981715545403807,
    -1.0086089494594537, -0.0, -0.0039640636842094108, -1.0092036964407751, -0.0,
    -0.0042274210062385031, -1.0097892451311354, -0.0, -0.0044882666232786052,
    -1.0103656806707684, -0.0, -0.0047466234278123473, -1.0109330874248688, -0.0,
    -0.0050025141040658241, -1.0114915489906429, -0.0, -0.0052559611299029833,
    -1.0120411482042946, -0.0, -0.0055069867787027873, -1.0125819671479479, -0.0,
    0.00030104618687214896, 1.0007224384072255, 0.0, 0.00059923720380821669,
    1.0014344978436069, 0.0, 0.000894599067828991, 1.0021362741967341, 0.0,
    0.0011871575592765928, 1.0028278624813931, 0.0, 0.0014769382239674068,
    1.0035093568475058, 0.0, 0.0017639663753254257, 1.0041808505879968, 0.0,
    0.002048267096496192, 1.0048424361465897, 0.0, 0.0023298652424415065,
    1.005494205125532, 0.0, 0.0026087854420150866, 1.006136248293249, 0.0,
    0.0028850521000193404, 1.00676865559193, 0.0, 0.0031586893992434336,
    1.0073915161450431, 0.0, 0.0034297213024828167, 1.0080049182647837, 0.0,
    0.0036981715545403807, 1.0086089494594537, 0.0, 0.0039640636842094108,
    1.0092036964407751, 0.0, 0.0042274210062385031, 1.0097892451311354, 0.0,
    0.0044882666232786052, 1.0103656806707684, 0.0, 0.0047466234278123473,
    1.0109330874248688, 0.0, 0.0050025141040658241, 1.0114915489906429, 0.0,
    0.0052559611299029833, 1.0120411482042946, 0.0, 0.0055069867787027873,
    1.0125819671479479, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
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

  static const real_T b_a[830]{ -0.9977021834978721, -0.0, -0.0, -0.025,
    -0.99540964695642165, -0.0, -0.0, -0.049942554587446807,
    -0.99312237824326788, -0.0, -0.0, -0.074827795761357341, -0.990840365253908,
    -0.0, -0.0, -0.099655855217439027, -0.9885635959116531, -0.0, -0.0,
    -0.12442686434878672, -0.98629205816756438, -0.0, -0.0, -0.14914095424657806,
    -0.98402574000038923, -0.0, -0.0, -0.17379825570076715, -0.98176462941649767,
    -0.0, -0.0, -0.19839889920077688, -0.97950871444981891, -0.0, -0.0,
    -0.22294301493618932, -0.977257983161778, -0.0, -0.0, -0.2474307327974348,
    -0.97501242364123264, -0.0, -0.0, -0.27186218237647924, -0.9727720240044101,
    -0.0, -0.0, -0.29623749296751006, -0.97053677239484437, -0.0, -0.0,
    -0.32055679356762035, -0.96830665698331353, -0.0, -0.0, -0.34482021287749148,
    -0.966081665967777, -0.0, -0.0, -0.36902787930207437, -0.963861787573313,
    -0.0, -0.0, -0.3931799209512688, -0.96164701005205655, -0.0, -0.0,
    -0.41727646564060167, -0.959437321683137, -0.0, -0.0, -0.44131764089190312,
    -0.95723271077261607, -0.0, -0.0, -0.46530357393398158, -0.95503316565342611,
    -0.0, -0.0, -0.489234391703297, 0.9977021834978721, 0.0, 0.0, 0.025,
    0.99540964695642165, 0.0, 0.0, 0.049942554587446807, 0.99312237824326788,
    0.0, 0.0, 0.074827795761357341, 0.990840365253908, 0.0, 0.0,
    0.099655855217439027, 0.9885635959116531, 0.0, 0.0, 0.12442686434878672,
    0.98629205816756438, 0.0, 0.0, 0.14914095424657806, 0.98402574000038923, 0.0,
    0.0, 0.17379825570076715, 0.98176462941649767, 0.0, 0.0, 0.19839889920077688,
    0.97950871444981891, 0.0, 0.0, 0.22294301493618932, 0.977257983161778, 0.0,
    0.0, 0.2474307327974348, 0.97501242364123264, 0.0, 0.0, 0.27186218237647924,
    0.9727720240044101, 0.0, 0.0, 0.29623749296751006, 0.97053677239484437, 0.0,
    0.0, 0.32055679356762035, 0.96830665698331353, 0.0, 0.0, 0.34482021287749148,
    0.966081665967777, 0.0, 0.0, 0.36902787930207437, 0.963861787573313, 0.0,
    0.0, 0.3931799209512688, 0.96164701005205655, 0.0, 0.0, 0.41727646564060167,
    0.959437321683137, 0.0, 0.0, 0.44131764089190312, 0.95723271077261607, 0.0,
    0.0, 0.46530357393398158, 0.95503316565342611, 0.0, 0.0, 0.489234391703297,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 };

  static const real_T b_Ac_0[824]{ -0.0, 0.0071318388640683669,
    0.012368697539837475, -0.0, -0.0, -0.0, 0.014263756082690558,
    0.024718035607500709, 0.00017829597160170919, 0.00030921743849593691, -0.0,
    0.021395752912117184, 0.037048025614083456, 0.00053488987366897319,
    0.00092716832868345472, -0.0, 0.028527830608018147, 0.0493586789617647,
    0.0010697836964719027, 0.0018533689690355412, -0.0, 0.035659990425483137,
    0.061650007043814951, 0.0017829794616723562, 0.0030873359430796592, -0.0,
    0.042792233619022145, 0.073922021244602523, 0.0026744792223094348,
    0.0046285861191750332, -0.0, 0.049924561442565955, 0.086174732939599832,
    0.0037442850627849886, 0.006476636650290097, -0.0, 0.057056975149466643,
    0.098408153495389644, 0.004992399098849137, 0.0086310049737800924, -0.0,
    0.0641894759924981, 0.11062229426967139, 0.0064188234775858031,
    0.011091208811164835, -0.0, 0.071322065223856493, 0.1228171666112674,
    0.0080235603773982558, 0.01385676616790662, -0.0, 0.078454744095160817,
    0.1349927818601292, 0.0098066120079946684, 0.016927195333188306, -0.0,
    0.085587513857453348, 0.14714915134734369, 0.011767980610373688,
    0.020302014879691535, -0.0, 0.09272037576120018, 0.1592862863951395,
    0.013907668456810022, 0.023980743663375128, -0.0, 0.099853331056291711,
    0.17140419831689321, 0.016225677850840024, 0.027962900823253617, -0.0,
    0.10698638099204312, 0.18350289841713557, 0.018722011127247317,
    0.032248005781175947, -0.0, 0.11411952681719491, 0.19558239799155772,
    0.021396670652048394, 0.036835578241604335, -0.0, 0.12125276977991335,
    0.20764270832701751, 0.024249658822478266, 0.041725138191393278, -0.0,
    0.12838611112779105, 0.21968384070154565, 0.0272809780669761,
    0.046916205899568721, -0.0, 0.13551955210784739, 0.23170580638435198,
    0.030490630845170875, 0.052408301917107367, -0.0, 0.14265309396652903,
    0.24370861663583168, 0.033878619647867057, 0.058200947076716171, 0.0,
    -0.0071318388640683669, -0.012368697539837475, 0.0, 0.0, 0.0,
    -0.014263756082690558, -0.024718035607500709, -0.00017829597160170919,
    -0.00030921743849593691, 0.0, -0.021395752912117184, -0.037048025614083456,
    -0.00053488987366897319, -0.00092716832868345472, 0.0, -0.028527830608018147,
    -0.0493586789617647, -0.0010697836964719027, -0.0018533689690355412, 0.0,
    -0.035659990425483137, -0.061650007043814951, -0.0017829794616723562,
    -0.0030873359430796592, 0.0, -0.042792233619022145, -0.073922021244602523,
    -0.0026744792223094348, -0.0046285861191750332, 0.0, -0.049924561442565955,
    -0.086174732939599832, -0.0037442850627849886, -0.006476636650290097, 0.0,
    -0.057056975149466643, -0.098408153495389644, -0.004992399098849137,
    -0.0086310049737800924, 0.0, -0.0641894759924981, -0.11062229426967139,
    -0.0064188234775858031, -0.011091208811164835, 0.0, -0.071322065223856493,
    -0.1228171666112674, -0.0080235603773982558, -0.01385676616790662, 0.0,
    -0.078454744095160817, -0.1349927818601292, -0.0098066120079946684,
    -0.016927195333188306, 0.0, -0.085587513857453348, -0.14714915134734369,
    -0.011767980610373688, -0.020302014879691535, 0.0, -0.09272037576120018,
    -0.1592862863951395, -0.013907668456810022, -0.023980743663375128, 0.0,
    -0.099853331056291711, -0.17140419831689321, -0.016225677850840024,
    -0.027962900823253617, 0.0, -0.10698638099204312, -0.18350289841713557,
    -0.018722011127247317, -0.032248005781175947, 0.0, -0.11411952681719491,
    -0.19558239799155772, -0.021396670652048394, -0.036835578241604335, 0.0,
    -0.12125276977991335, -0.20764270832701751, -0.024249658822478266,
    -0.041725138191393278, 0.0, -0.12838611112779105, -0.21968384070154565,
    -0.0272809780669761, -0.046916205899568721, 0.0, -0.13551955210784739,
    -0.23170580638435198, -0.030490630845170875, -0.052408301917107367, 0.0,
    -0.14265309396652903, -0.24370861663583168, -0.033878619647867057,
    -0.058200947076716171, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, -0.0,
    -0.13659152296751884, 0.063538603002101665, -0.0, -0.0, -0.0,
    -0.27320389247378107, 0.12726886009414035, -0.0034147880741879711,
    0.0015884650750525417, -0.0, -0.40983712341735679, 0.19119069203750322,
    -0.010244885386032498, 0.00477018657740605, -0.0, -0.54649123069354111,
    0.25530401966585881, -0.020490813471466418, 0.00954995387834363, -0.0,
    -0.68316622919435854, 0.31960876388510828, -0.03415309423880495,
    0.0159325543699901, -0.0, -0.81986213380856687, 0.38410484567333708,
    -0.051232249968663914, 0.023922773467117804, -0.0, -0.95657895942166182,
    0.44879218608076621, -0.071728803313878084, 0.03352539460895123, -0.0,
    -1.0933167209158812, 0.51367070622970379, -0.095643277299419627,
    0.044745199260970384, -0.0, -1.2300754331702088, 0.57874032731449654,
    -0.12297619532231666, 0.057586966916712973, -0.0, -1.3668551110603793,
    0.64400097060148154, -0.15372808115157188, 0.072055475099575383, -0.0,
    -1.5036557694588819, 0.70945255742893776, -0.18789945892808135,
    0.088155499364612422, -0.0, -1.6404774232349646, 0.77509500920703767,
    -0.22549085316455342, 0.10589181330033587, -0.0, -1.7773200872546391,
    0.84092824741779892, -0.26650278874542754, 0.1252691885305118, -0.0,
    -1.9141837763806839, 0.90695219361503621, -0.31093579092679352,
    0.14629239471595679, -0.0, -2.0510685054726494, 0.97316676942431268,
    -0.35879038533631064, 0.16896619955633269, -0.0, -2.1879742893868621,
    1.0395718965428922, -0.4100670979731269, 0.19329536879194051, -0.0,
    -2.3249011429764286, 1.1061674967396908, -0.46476645520779847,
    0.21928466620551282, -0.0, -2.4618490810912395, 1.1729534918552289,
    -0.52288898378220916, 0.24693885362400511, -0.0, -2.5988181185779742,
    1.2399298038015829, -0.58443521080949012, 0.27626269092038586, -0.0,
    -2.735808270280105, 1.3070963545623373, -0.64940566377393949,
    0.30726093601542542, 0.0, 0.13659152296751884, -0.063538603002101665, 0.0,
    0.0, 0.0, 0.27320389247378107, -0.12726886009414035, 0.0034147880741879711,
    -0.0015884650750525417, 0.0, 0.40983712341735679, -0.19119069203750322,
    0.010244885386032498, -0.00477018657740605, 0.0, 0.54649123069354111,
    -0.25530401966585881, 0.020490813471466418, -0.00954995387834363, 0.0,
    0.68316622919435854, -0.31960876388510828, 0.03415309423880495,
    -0.0159325543699901, 0.0, 0.81986213380856687, -0.38410484567333708,
    0.051232249968663914, -0.023922773467117804, 0.0, 0.95657895942166182,
    -0.44879218608076621, 0.071728803313878084, -0.03352539460895123, 0.0,
    1.0933167209158812, -0.51367070622970379, 0.095643277299419627,
    -0.044745199260970384, 0.0, 1.2300754331702088, -0.57874032731449654,
    0.12297619532231666, -0.057586966916712973, 0.0, 1.3668551110603793,
    -0.64400097060148154, 0.15372808115157188, -0.072055475099575383, 0.0,
    1.5036557694588819, -0.70945255742893776, 0.18789945892808135,
    -0.088155499364612422, 0.0, 1.6404774232349646, -0.77509500920703767,
    0.22549085316455342, -0.10589181330033587, 0.0, 1.7773200872546391,
    -0.84092824741779892, 0.26650278874542754, -0.1252691885305118, 0.0,
    1.9141837763806839, -0.90695219361503621, 0.31093579092679352,
    -0.14629239471595679, 0.0, 2.0510685054726494, -0.97316676942431268,
    0.35879038533631064, -0.16896619955633269, 0.0, 2.1879742893868621,
    -1.0395718965428922, 0.4100670979731269, -0.19329536879194051, 0.0,
    2.3249011429764286, -1.1061674967396908, 0.46476645520779847,
    -0.21928466620551282, 0.0, 2.4618490810912395, -1.1729534918552289,
    0.52288898378220916, -0.24693885362400511, 0.0, 2.5988181185779742,
    -1.2399298038015829, 0.58443521080949012, -0.27626269092038586, 0.0,
    2.735808270280105, -1.3070963545623373, 0.64940566377393949,
    -0.30726093601542542, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, -0.0,
    0.04588890657825119, -0.10453377294976968, -0.0, -0.0, -0.0,
    0.091790173501845068, -0.2090823342837666, 0.0011472226644562798,
    -0.002613344323744242, -0.0, 0.13770380323925227, -0.31364569595584796,
    0.0034419770020024066, -0.007840402680838408, -0.0, 0.183629798260016,
    -0.41822386991689187, 0.0068845720829837139, -0.015681545079734606, -0.0,
    0.22956816103475194, -0.52281686811480088, 0.011475317039484116,
    -0.026137141827656904, -0.0, 0.27551889403514823, -0.62742470249450566,
    0.017214521065352914, -0.039207563530526922, -0.0, 0.32148199973396535,
    -0.7320473849979684, 0.024102493416231621, -0.054893181092889556, -0.0,
    0.36745748060503614, -0.83668492756418655, 0.032139543409580758,
    -0.073194365717838758, -0.0, 0.41344533912326564, -0.9413373421291964,
    0.041325980424706657, -0.094111488906943416, -0.0, 0.45944557776463113,
    -1.0460046406260766, 0.0516621139027883, -0.11764492246017333, -0.0,
    0.505458199006182, -1.1506868349849515, 0.063148253346904076,
    -0.14379503847582525, -0.0, 0.55148320532603978, -1.2553839371329951,
    0.075784708322058622, -0.17256220935044903, -0.0, 0.597520599203398,
    -1.3600959589944344, 0.089571788455209619, -0.20394680777877389, -0.0,
    0.643570383118522, -1.4648229124905527, 0.10450980343529456,
    -0.23794920675363476, -0.0, 0.68963255955274927, -1.5695648095396941,
    0.12059906301325761, -0.27456977956589856, -0.0, 0.73570713098848894,
    -1.6743216620572658, 0.13783987700207634, -0.31380889980439092, -0.0,
    0.78179409990922211, -1.7790934819557425, 0.15623255527678856,
    -0.35566694135582255, -0.0, 0.82789346879950143, -1.8838802811446695,
    0.17577740777451911, -0.40014427840471611, -0.0, 0.8740052401449514,
    -1.9886820715306668, 0.19647474449450664, -0.44724128543333286, -0.0,
    0.920129416432268, -2.0934988650174318, 0.21832487549813043,
    -0.49695833722159954, 0.0, -0.04588890657825119, 0.10453377294976968, 0.0,
    0.0, 0.0, -0.091790173501845068, 0.2090823342837666, -0.0011472226644562798,
    0.002613344323744242, 0.0, -0.13770380323925227, 0.31364569595584796,
    -0.0034419770020024066, 0.007840402680838408, 0.0, -0.183629798260016,
    0.41822386991689187, -0.0068845720829837139, 0.015681545079734606, 0.0,
    -0.22956816103475194, 0.52281686811480088, -0.011475317039484116,
    0.026137141827656904, 0.0, -0.27551889403514823, 0.62742470249450566,
    -0.017214521065352914, 0.039207563530526922, 0.0, -0.32148199973396535,
    0.7320473849979684, -0.024102493416231621, 0.054893181092889556, 0.0,
    -0.36745748060503614, 0.83668492756418655, -0.032139543409580758,
    0.073194365717838758, 0.0, -0.41344533912326564, 0.9413373421291964,
    -0.041325980424706657, 0.094111488906943416, 0.0, -0.45944557776463113,
    1.0460046406260766, -0.0516621139027883, 0.11764492246017333, 0.0,
    -0.505458199006182, 1.1506868349849515, -0.063148253346904076,
    0.14379503847582525, 0.0, -0.55148320532603978, 1.2553839371329951,
    -0.075784708322058622, 0.17256220935044903, 0.0, -0.597520599203398,
    1.3600959589944344, -0.089571788455209619, 0.20394680777877389, 0.0,
    -0.643570383118522, 1.4648229124905527, -0.10450980343529456,
    0.23794920675363476, 0.0, -0.68963255955274927, 1.5695648095396941,
    -0.12059906301325761, 0.27456977956589856, 0.0, -0.73570713098848894,
    1.6743216620572658, -0.13783987700207634, 0.31380889980439092, 0.0,
    -0.78179409990922211, 1.7790934819557425, -0.15623255527678856,
    0.35566694135582255, 0.0, -0.82789346879950143, 1.8838802811446695,
    -0.17577740777451911, 0.40014427840471611, 0.0, -0.8740052401449514,
    1.9886820715306668, -0.19647474449450664, 0.44724128543333286, 0.0,
    -0.920129416432268, 2.0934988650174318, -0.21832487549813043,
    0.49695833722159954, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
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

  static const real_T b_Ac[664]{ -0.22599754912864506, -0.0, -0.0, -0.0,
    -0.45147579735946186, -0.0, -0.0, -0.0056499387282161265, -0.676435937950623,
    -0.0, -0.0, -0.016936833662202673, -0.90087916141841273, -0.0, -0.0,
    -0.033847732110968247, -1.1248066555435274, -0.0, -0.0,
    -0.056369711146428567, -1.3482196053773612, -0.0, -0.0,
    -0.084489877535016744, -1.5711191932482778, -0.0, -0.0, -0.11819536766945077,
    -1.7935065987678671, -0.0, -0.0, -0.15747334750065772, -2.0153829988371879,
    -0.0, -0.0, -0.20231101246985439, -2.2367495676529967, -0.0, -0.0,
    -0.25269558744078408, -2.4576074767139611, -0.0, -0.0, -0.30861432663210897,
    -2.67795789482686, -0.0, -0.0, -0.370054513549958, -2.8978019881127679, -0.0,
    -0.0, -0.43700346092062947, -3.1171409200132283, -0.0, -0.0,
    -0.50944851062344865, -3.3359758512964088, -0.0, -0.0, -0.58737703362377935,
    -3.5543079400632447, -0.0, -0.0, -0.67077642990618958, -3.7721383417535681,
    -0.0, -0.0, -0.75963412840777078, -3.9894682091522222, -0.0, -0.0,
    -0.85393758695161, -4.2062986923951629, -0.0, -0.0, -0.95367429218041566,
    -4.4226309389755434, -0.0, -0.0, -1.0588317594902947, 0.22599754912864506,
    0.0, 0.0, 0.0, 0.45147579735946186, 0.0, 0.0, 0.0056499387282161265,
    0.676435937950623, 0.0, 0.0, 0.016936833662202673, 0.90087916141841273, 0.0,
    0.0, 0.033847732110968247, 1.1248066555435274, 0.0, 0.0,
    0.056369711146428567, 1.3482196053773612, 0.0, 0.0, 0.084489877535016744,
    1.5711191932482778, 0.0, 0.0, 0.11819536766945077, 1.7935065987678671, 0.0,
    0.0, 0.15747334750065772, 2.0153829988371879, 0.0, 0.0, 0.20231101246985439,
    2.2367495676529967, 0.0, 0.0, 0.25269558744078408, 2.4576074767139611, 0.0,
    0.0, 0.30861432663210897, 2.67795789482686, 0.0, 0.0, 0.370054513549958,
    2.8978019881127679, 0.0, 0.0, 0.43700346092062947, 3.1171409200132283, 0.0,
    0.0, 0.50944851062344865, 3.3359758512964088, 0.0, 0.0, 0.58737703362377935,
    3.5543079400632447, 0.0, 0.0, 0.67077642990618958, 3.7721383417535681, 0.0,
    0.0, 0.75963412840777078, 3.9894682091522222, 0.0, 0.0, 0.85393758695161,
    4.2062986923951629, 0.0, 0.0, 0.95367429218041566, 4.4226309389755434, 0.0,
    0.0, 1.0588317594902947, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.087367622049086546, -0.0, -0.0, -0.0, 0.17453448933447702, -0.0, -0.0,
    0.0021841905512271637, 0.26150106315378036, -0.0, -0.0, 0.00654755278458909,
    0.34826780374462818, -0.0, -0.0, 0.013085079363433599, 0.4348351702871105,
    -0.0, -0.0, 0.0217917744570493, 0.52120362090620576, -0.0, -0.0,
    0.032662653714227059, 0.6073736126742052, -0.0, -0.0, 0.0456927442368822,
    0.69334560161313186, -0.0, -0.0, 0.060877084553737332, 0.779120042697154,
    -0.0, -0.0, 0.078210724594065625, 0.86469738985499234, -0.0, -0.0,
    0.097688725661494474, 0.95007809597232318, -0.0, -0.0, 0.11930616040786929,
    1.0352626128941742, -0.0, -0.0, 0.14305811280717737, 1.1202513914273164,
    -0.0, -0.0, 0.16893967812953173, 1.2050448813426495, -0.0, -0.0,
    0.19694596291521463, 1.2896435313775823, -0.0, -0.0, 0.22707208494878087,
    1.3740477892384069, -0.0, -0.0, 0.25931317323322045, 1.458258101602669, -0.0,
    -0.0, 0.29366436796418061, 1.5422749141215313, -0.0, -0.0,
    0.33012082050424735, 1.6260986714221315, -0.0, -0.0, 0.36867769335728562,
    1.709729817109936, -0.0, -0.0, 0.40933016014283891, -0.087367622049086546,
    0.0, 0.0, 0.0, -0.17453448933447702, 0.0, 0.0, -0.0021841905512271637,
    -0.26150106315378036, 0.0, 0.0, -0.00654755278458909, -0.34826780374462818,
    0.0, 0.0, -0.013085079363433599, -0.4348351702871105, 0.0, 0.0,
    -0.0217917744570493, -0.52120362090620576, 0.0, 0.0, -0.032662653714227059,
    -0.6073736126742052, 0.0, 0.0, -0.0456927442368822, -0.69334560161313186,
    0.0, 0.0, -0.060877084553737332, -0.779120042697154, 0.0, 0.0,
    -0.078210724594065625, -0.86469738985499234, 0.0, 0.0, -0.097688725661494474,
    -0.95007809597232318, 0.0, 0.0, -0.11930616040786929, -1.0352626128941742,
    0.0, 0.0, -0.14305811280717737, -1.1202513914273164, 0.0, 0.0,
    -0.16893967812953173, -1.2050448813426495, 0.0, 0.0, -0.19694596291521463,
    -1.2896435313775823, 0.0, 0.0, -0.22707208494878087, -1.3740477892384069,
    0.0, 0.0, -0.25931317323322045, -1.458258101602669, 0.0, 0.0,
    -0.29366436796418061, -1.5422749141215313, 0.0, 0.0, -0.33012082050424735,
    -1.6260986714221315, 0.0, 0.0, -0.36867769335728562, -1.709729817109936, 0.0,
    0.0, -0.40933016014283891, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.087367622049086546, -0.0, -0.0, -0.0, 0.17453448933447702, -0.0, -0.0,
    0.0021841905512271637, 0.26150106315378036, -0.0, -0.0, 0.00654755278458909,
    0.34826780374462818, -0.0, -0.0, 0.013085079363433599, 0.4348351702871105,
    -0.0, -0.0, 0.0217917744570493, 0.52120362090620576, -0.0, -0.0,
    0.032662653714227059, 0.6073736126742052, -0.0, -0.0, 0.0456927442368822,
    0.69334560161313186, -0.0, -0.0, 0.060877084553737332, 0.779120042697154,
    -0.0, -0.0, 0.078210724594065625, 0.86469738985499234, -0.0, -0.0,
    0.097688725661494474, 0.95007809597232318, -0.0, -0.0, 0.11930616040786929,
    1.0352626128941742, -0.0, -0.0, 0.14305811280717737, 1.1202513914273164,
    -0.0, -0.0, 0.16893967812953173, 1.2050448813426495, -0.0, -0.0,
    0.19694596291521463, 1.2896435313775823, -0.0, -0.0, 0.22707208494878087,
    1.3740477892384069, -0.0, -0.0, 0.25931317323322045, 1.458258101602669, -0.0,
    -0.0, 0.29366436796418061, 1.5422749141215313, -0.0, -0.0,
    0.33012082050424735, 1.6260986714221315, -0.0, -0.0, 0.36867769335728562,
    1.709729817109936, -0.0, -0.0, 0.40933016014283891, -0.087367622049086546,
    0.0, 0.0, 0.0, -0.17453448933447702, 0.0, 0.0, -0.0021841905512271637,
    -0.26150106315378036, 0.0, 0.0, -0.00654755278458909, -0.34826780374462818,
    0.0, 0.0, -0.013085079363433599, -0.4348351702871105, 0.0, 0.0,
    -0.0217917744570493, -0.52120362090620576, 0.0, 0.0, -0.032662653714227059,
    -0.6073736126742052, 0.0, 0.0, -0.0456927442368822, -0.69334560161313186,
    0.0, 0.0, -0.060877084553737332, -0.779120042697154, 0.0, 0.0,
    -0.078210724594065625, -0.86469738985499234, 0.0, 0.0, -0.097688725661494474,
    -0.95007809597232318, 0.0, 0.0, -0.11930616040786929, -1.0352626128941742,
    0.0, 0.0, -0.14305811280717737, -1.1202513914273164, 0.0, 0.0,
    -0.16893967812953173, -1.2050448813426495, 0.0, 0.0, -0.19694596291521463,
    -1.2896435313775823, 0.0, 0.0, -0.22707208494878087, -1.3740477892384069,
    0.0, 0.0, -0.25931317323322045, -1.458258101602669, 0.0, 0.0,
    -0.29366436796418061, -1.5422749141215313, 0.0, 0.0, -0.33012082050424735,
    -1.6260986714221315, 0.0, 0.0, -0.36867769335728562, -1.709729817109936, 0.0,
    0.0, -0.40933016014283891, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0,
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
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T a_0[618]{ -0.0, 0.0071318388640683669,
    0.012368697539837475, -0.0, -0.0, -0.0, 0.014263756082690558,
    0.024718035607500709, 0.00017829597160170919, 0.00030921743849593691, -0.0,
    0.021395752912117184, 0.037048025614083456, 0.00053488987366897319,
    0.00092716832868345472, -0.0, 0.028527830608018147, 0.0493586789617647,
    0.0010697836964719027, 0.0018533689690355412, -0.0, 0.035659990425483137,
    0.061650007043814951, 0.0017829794616723562, 0.0030873359430796592, -0.0,
    0.042792233619022145, 0.073922021244602523, 0.0026744792223094348,
    0.0046285861191750332, -0.0, 0.049924561442565955, 0.086174732939599832,
    0.0037442850627849886, 0.006476636650290097, -0.0, 0.057056975149466643,
    0.098408153495389644, 0.004992399098849137, 0.0086310049737800924, -0.0,
    0.0641894759924981, 0.11062229426967139, 0.0064188234775858031,
    0.011091208811164835, -0.0, 0.071322065223856493, 0.1228171666112674,
    0.0080235603773982558, 0.01385676616790662, -0.0, 0.078454744095160817,
    0.1349927818601292, 0.0098066120079946684, 0.016927195333188306, -0.0,
    0.085587513857453348, 0.14714915134734369, 0.011767980610373688,
    0.020302014879691535, -0.0, 0.09272037576120018, 0.1592862863951395,
    0.013907668456810022, 0.023980743663375128, -0.0, 0.099853331056291711,
    0.17140419831689321, 0.016225677850840024, 0.027962900823253617, -0.0,
    0.10698638099204312, 0.18350289841713557, 0.018722011127247317,
    0.032248005781175947, -0.0, 0.11411952681719491, 0.19558239799155772,
    0.021396670652048394, 0.036835578241604335, -0.0, 0.12125276977991335,
    0.20764270832701751, 0.024249658822478266, 0.041725138191393278, -0.0,
    0.12838611112779105, 0.21968384070154565, 0.0272809780669761,
    0.046916205899568721, -0.0, 0.13551955210784739, 0.23170580638435198,
    0.030490630845170875, 0.052408301917107367, -0.0, 0.14265309396652903,
    0.24370861663583168, 0.033878619647867057, 0.058200947076716171, 0.0,
    -0.0071318388640683669, -0.012368697539837475, 0.0, 0.0, 0.0,
    -0.014263756082690558, -0.024718035607500709, -0.00017829597160170919,
    -0.00030921743849593691, 0.0, -0.021395752912117184, -0.037048025614083456,
    -0.00053488987366897319, -0.00092716832868345472, 0.0, -0.028527830608018147,
    -0.0493586789617647, -0.0010697836964719027, -0.0018533689690355412, 0.0,
    -0.035659990425483137, -0.061650007043814951, -0.0017829794616723562,
    -0.0030873359430796592, 0.0, -0.042792233619022145, -0.073922021244602523,
    -0.0026744792223094348, -0.0046285861191750332, 0.0, -0.049924561442565955,
    -0.086174732939599832, -0.0037442850627849886, -0.006476636650290097, 0.0,
    -0.057056975149466643, -0.098408153495389644, -0.004992399098849137,
    -0.0086310049737800924, 0.0, -0.0641894759924981, -0.11062229426967139,
    -0.0064188234775858031, -0.011091208811164835, 0.0, -0.071322065223856493,
    -0.1228171666112674, -0.0080235603773982558, -0.01385676616790662, 0.0,
    -0.078454744095160817, -0.1349927818601292, -0.0098066120079946684,
    -0.016927195333188306, 0.0, -0.085587513857453348, -0.14714915134734369,
    -0.011767980610373688, -0.020302014879691535, 0.0, -0.09272037576120018,
    -0.1592862863951395, -0.013907668456810022, -0.023980743663375128, 0.0,
    -0.099853331056291711, -0.17140419831689321, -0.016225677850840024,
    -0.027962900823253617, 0.0, -0.10698638099204312, -0.18350289841713557,
    -0.018722011127247317, -0.032248005781175947, 0.0, -0.11411952681719491,
    -0.19558239799155772, -0.021396670652048394, -0.036835578241604335, 0.0,
    -0.12125276977991335, -0.20764270832701751, -0.024249658822478266,
    -0.041725138191393278, 0.0, -0.12838611112779105, -0.21968384070154565,
    -0.0272809780669761, -0.046916205899568721, 0.0, -0.13551955210784739,
    -0.23170580638435198, -0.030490630845170875, -0.052408301917107367, 0.0,
    -0.14265309396652903, -0.24370861663583168, -0.033878619647867057,
    -0.058200947076716171, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, -0.0,
    -0.13659152296751884, 0.063538603002101665, -0.0, -0.0, -0.0,
    -0.27320389247378107, 0.12726886009414035, -0.0034147880741879711,
    0.0015884650750525417, -0.0, -0.40983712341735679, 0.19119069203750322,
    -0.010244885386032498, 0.00477018657740605, -0.0, -0.54649123069354111,
    0.25530401966585881, -0.020490813471466418, 0.00954995387834363, -0.0,
    -0.68316622919435854, 0.31960876388510828, -0.03415309423880495,
    0.0159325543699901, -0.0, -0.81986213380856687, 0.38410484567333708,
    -0.051232249968663914, 0.023922773467117804, -0.0, -0.95657895942166182,
    0.44879218608076621, -0.071728803313878084, 0.03352539460895123, -0.0,
    -1.0933167209158812, 0.51367070622970379, -0.095643277299419627,
    0.044745199260970384, -0.0, -1.2300754331702088, 0.57874032731449654,
    -0.12297619532231666, 0.057586966916712973, -0.0, -1.3668551110603793,
    0.64400097060148154, -0.15372808115157188, 0.072055475099575383, -0.0,
    -1.5036557694588819, 0.70945255742893776, -0.18789945892808135,
    0.088155499364612422, -0.0, -1.6404774232349646, 0.77509500920703767,
    -0.22549085316455342, 0.10589181330033587, -0.0, -1.7773200872546391,
    0.84092824741779892, -0.26650278874542754, 0.1252691885305118, -0.0,
    -1.9141837763806839, 0.90695219361503621, -0.31093579092679352,
    0.14629239471595679, -0.0, -2.0510685054726494, 0.97316676942431268,
    -0.35879038533631064, 0.16896619955633269, -0.0, -2.1879742893868621,
    1.0395718965428922, -0.4100670979731269, 0.19329536879194051, -0.0,
    -2.3249011429764286, 1.1061674967396908, -0.46476645520779847,
    0.21928466620551282, -0.0, -2.4618490810912395, 1.1729534918552289,
    -0.52288898378220916, 0.24693885362400511, -0.0, -2.5988181185779742,
    1.2399298038015829, -0.58443521080949012, 0.27626269092038586, -0.0,
    -2.735808270280105, 1.3070963545623373, -0.64940566377393949,
    0.30726093601542542, 0.0, 0.13659152296751884, -0.063538603002101665, 0.0,
    0.0, 0.0, 0.27320389247378107, -0.12726886009414035, 0.0034147880741879711,
    -0.0015884650750525417, 0.0, 0.40983712341735679, -0.19119069203750322,
    0.010244885386032498, -0.00477018657740605, 0.0, 0.54649123069354111,
    -0.25530401966585881, 0.020490813471466418, -0.00954995387834363, 0.0,
    0.68316622919435854, -0.31960876388510828, 0.03415309423880495,
    -0.0159325543699901, 0.0, 0.81986213380856687, -0.38410484567333708,
    0.051232249968663914, -0.023922773467117804, 0.0, 0.95657895942166182,
    -0.44879218608076621, 0.071728803313878084, -0.03352539460895123, 0.0,
    1.0933167209158812, -0.51367070622970379, 0.095643277299419627,
    -0.044745199260970384, 0.0, 1.2300754331702088, -0.57874032731449654,
    0.12297619532231666, -0.057586966916712973, 0.0, 1.3668551110603793,
    -0.64400097060148154, 0.15372808115157188, -0.072055475099575383, 0.0,
    1.5036557694588819, -0.70945255742893776, 0.18789945892808135,
    -0.088155499364612422, 0.0, 1.6404774232349646, -0.77509500920703767,
    0.22549085316455342, -0.10589181330033587, 0.0, 1.7773200872546391,
    -0.84092824741779892, 0.26650278874542754, -0.1252691885305118, 0.0,
    1.9141837763806839, -0.90695219361503621, 0.31093579092679352,
    -0.14629239471595679, 0.0, 2.0510685054726494, -0.97316676942431268,
    0.35879038533631064, -0.16896619955633269, 0.0, 2.1879742893868621,
    -1.0395718965428922, 0.4100670979731269, -0.19329536879194051, 0.0,
    2.3249011429764286, -1.1061674967396908, 0.46476645520779847,
    -0.21928466620551282, 0.0, 2.4618490810912395, -1.1729534918552289,
    0.52288898378220916, -0.24693885362400511, 0.0, 2.5988181185779742,
    -1.2399298038015829, 0.58443521080949012, -0.27626269092038586, 0.0,
    2.735808270280105, -1.3070963545623373, 0.64940566377393949,
    -0.30726093601542542, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, -0.0,
    0.04588890657825119, -0.10453377294976968, -0.0, -0.0, -0.0,
    0.091790173501845068, -0.2090823342837666, 0.0011472226644562798,
    -0.002613344323744242, -0.0, 0.13770380323925227, -0.31364569595584796,
    0.0034419770020024066, -0.007840402680838408, -0.0, 0.183629798260016,
    -0.41822386991689187, 0.0068845720829837139, -0.015681545079734606, -0.0,
    0.22956816103475194, -0.52281686811480088, 0.011475317039484116,
    -0.026137141827656904, -0.0, 0.27551889403514823, -0.62742470249450566,
    0.017214521065352914, -0.039207563530526922, -0.0, 0.32148199973396535,
    -0.7320473849979684, 0.024102493416231621, -0.054893181092889556, -0.0,
    0.36745748060503614, -0.83668492756418655, 0.032139543409580758,
    -0.073194365717838758, -0.0, 0.41344533912326564, -0.9413373421291964,
    0.041325980424706657, -0.094111488906943416, -0.0, 0.45944557776463113,
    -1.0460046406260766, 0.0516621139027883, -0.11764492246017333, -0.0,
    0.505458199006182, -1.1506868349849515, 0.063148253346904076,
    -0.14379503847582525, -0.0, 0.55148320532603978, -1.2553839371329951,
    0.075784708322058622, -0.17256220935044903, -0.0, 0.597520599203398,
    -1.3600959589944344, 0.089571788455209619, -0.20394680777877389, -0.0,
    0.643570383118522, -1.4648229124905527, 0.10450980343529456,
    -0.23794920675363476, -0.0, 0.68963255955274927, -1.5695648095396941,
    0.12059906301325761, -0.27456977956589856, -0.0, 0.73570713098848894,
    -1.6743216620572658, 0.13783987700207634, -0.31380889980439092, -0.0,
    0.78179409990922211, -1.7790934819557425, 0.15623255527678856,
    -0.35566694135582255, -0.0, 0.82789346879950143, -1.8838802811446695,
    0.17577740777451911, -0.40014427840471611, -0.0, 0.8740052401449514,
    -1.9886820715306668, 0.19647474449450664, -0.44724128543333286, -0.0,
    0.920129416432268, -2.0934988650174318, 0.21832487549813043,
    -0.49695833722159954, 0.0, -0.04588890657825119, 0.10453377294976968, 0.0,
    0.0, 0.0, -0.091790173501845068, 0.2090823342837666, -0.0011472226644562798,
    0.002613344323744242, 0.0, -0.13770380323925227, 0.31364569595584796,
    -0.0034419770020024066, 0.007840402680838408, 0.0, -0.183629798260016,
    0.41822386991689187, -0.0068845720829837139, 0.015681545079734606, 0.0,
    -0.22956816103475194, 0.52281686811480088, -0.011475317039484116,
    0.026137141827656904, 0.0, -0.27551889403514823, 0.62742470249450566,
    -0.017214521065352914, 0.039207563530526922, 0.0, -0.32148199973396535,
    0.7320473849979684, -0.024102493416231621, 0.054893181092889556, 0.0,
    -0.36745748060503614, 0.83668492756418655, -0.032139543409580758,
    0.073194365717838758, 0.0, -0.41344533912326564, 0.9413373421291964,
    -0.041325980424706657, 0.094111488906943416, 0.0, -0.45944557776463113,
    1.0460046406260766, -0.0516621139027883, 0.11764492246017333, 0.0,
    -0.505458199006182, 1.1506868349849515, -0.063148253346904076,
    0.14379503847582525, 0.0, -0.55148320532603978, 1.2553839371329951,
    -0.075784708322058622, 0.17256220935044903, 0.0, -0.597520599203398,
    1.3600959589944344, -0.089571788455209619, 0.20394680777877389, 0.0,
    -0.643570383118522, 1.4648229124905527, -0.10450980343529456,
    0.23794920675363476, 0.0, -0.68963255955274927, 1.5695648095396941,
    -0.12059906301325761, 0.27456977956589856, 0.0, -0.73570713098848894,
    1.6743216620572658, -0.13783987700207634, 0.31380889980439092, 0.0,
    -0.78179409990922211, 1.7790934819557425, -0.15623255527678856,
    0.35566694135582255, 0.0, -0.82789346879950143, 1.8838802811446695,
    -0.17577740777451911, 0.40014427840471611, 0.0, -0.8740052401449514,
    1.9886820715306668, -0.19647474449450664, 0.44724128543333286, 0.0,
    -0.920129416432268, 2.0934988650174318, -0.21832487549813043,
    0.49695833722159954, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T b_Ac_1[504]{ -0.025118718490858984, -0.016819265413734148,
    -0.0, -0.0499861237076248, -0.032741132441156645, -0.0,
    -0.074604510545374814, -0.047774052825594666, -0.0, -0.0989761530204811,
    -0.0619264013180454, -0.0, -0.12310330446053273, -0.075206476377572148, -0.0,
    -0.14698819769253046, -0.087622500865329314, -0.0, -0.17063304522936948,
    -0.099182622732274267, -0.0, -0.19404003945462597, -0.1098949157006235, -0.0,
    -0.21721135280566262, -0.11976737993910992, -0.0, -0.24014913795506856,
    -0.12880794273209781, -0.0, -0.26285552799044903, -0.13702445914261108, -0.0,
    -0.28533263659257935, -0.14442471266933057, -0.0, -0.30758255821193842,
    -0.1510164158976148, -0.0, -0.32960736824363657, -0.15680721114459886, -0.0,
    -0.35140912320075196, -0.16180467109842528, -0.0, -0.37298986088609037,
    -0.16601629945166, -0.0, -0.39435160056238255, -0.16944953152894673, -0.0,
    -0.41549634312093342, -0.17211173490895176, -0.0, -0.43642607124873706,
    -0.17401021004065131, -0.0, -0.45714274959407175, -0.17515219085401298, -0.0,
    0.025118718490858984, 0.016819265413734148, 0.0, 0.0499861237076248,
    0.032741132441156645, 0.0, 0.074604510545374814, 0.047774052825594666, 0.0,
    0.0989761530204811, 0.0619264013180454, 0.0, 0.12310330446053273,
    0.075206476377572148, 0.0, 0.14698819769253046, 0.087622500865329314, 0.0,
    0.17063304522936948, 0.099182622732274267, 0.0, 0.19404003945462597,
    0.1098949157006235, 0.0, 0.21721135280566262, 0.11976737993910992, 0.0,
    0.24014913795506856, 0.12880794273209781, 0.0, 0.26285552799044903,
    0.13702445914261108, 0.0, 0.28533263659257935, 0.14442471266933057, 0.0,
    0.30758255821193842, 0.1510164158976148, 0.0, 0.32960736824363657,
    0.15680721114459886, 0.0, 0.35140912320075196, 0.16180467109842528, 0.0,
    0.37298986088609037, 0.16601629945166, 0.0, 0.39435160056238255,
    0.16944953152894673, 0.0, 0.41549634312093342, 0.17211173490895176, 0.0,
    0.43642607124873706, 0.17401021004065131, 0.0, 0.45714274959407175,
    0.17515219085401298, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.0189087824026894, -0.066906613140151408, -0.0, 0.03760442849307085,
    -0.13454624955425293, -0.0, 0.056088892993839096, -0.20291172113658829, -0.0,
    0.0743641128405729, -0.27199590536902313, -0.0, 0.092432007343536085,
    -0.34179174472431623, -0.0, 0.11029447834800676, -0.41229224607485848, -0.0,
    0.12795341039314814, -0.48349048010679008, -0.0, 0.14541067086943421,
    -0.55537958073944738, -0.0, 0.16266811017464336, -0.62795274455009054, -0.0,
    0.17972756186843311, -0.70120323020386432, -0.0, 0.19659084282550868,
    -0.775124357888944, -0.0, 0.21325975338739828, -0.84970950875681928, -0.0,
    0.22973607751284789, -0.9249521243676706, -0.0, 0.24602158292684784,
    -1.0008457061407889, -0.0, 0.26211802126830397, -1.0773838148099963, -0.0,
    0.27802712823636522, -1.1545600698840195, -0.0, 0.2937506237354206,
    -1.2323681491117719, -0.0, 0.30929021201877682, -1.3108017879525005, -0.0,
    0.32464758183102904, -1.3898547790507516, -0.0, 0.33982440654913676,
    -1.4695209717161131, -0.0, -0.0189087824026894, 0.066906613140151408, 0.0,
    -0.03760442849307085, 0.13454624955425293, 0.0, -0.056088892993839096,
    0.20291172113658829, 0.0, -0.0743641128405729, 0.27199590536902313, 0.0,
    -0.092432007343536085, 0.34179174472431623, 0.0, -0.11029447834800676,
    0.41229224607485848, 0.0, -0.12795341039314814, 0.48349048010679008, 0.0,
    -0.14541067086943421, 0.55537958073944738, 0.0, -0.16266811017464336,
    0.62795274455009054, 0.0, -0.17972756186843311, 0.70120323020386432, 0.0,
    -0.19659084282550868, 0.775124357888944, 0.0, -0.21325975338739828,
    0.84970950875681928, 0.0, -0.22973607751284789, 0.9249521243676706, 0.0,
    -0.24602158292684784, 1.0008457061407889, 0.0, -0.26211802126830397,
    1.0773838148099963, 0.0, -0.27802712823636522, 1.1545600698840195, 0.0,
    -0.2937506237354206, 1.2323681491117719, 0.0, -0.30929021201877682,
    1.3108017879525005, 0.0, -0.32464758183102904, 1.3898547790507516, 0.0,
    -0.33982440654913676, 1.4695209717161131, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0,
    0.0, -0.032190816623613963, 0.026334884546125613, -0.0, -0.0640451465215906,
    0.053854424500504594, -0.0, -0.095566067442515334, 0.082547291513030641,
    -0.0, -0.1267566291319891, 0.11240226049489649, -0.0, -0.15761985358735717,
    0.14340820867920523, -0.0, -0.18815873531012073, 0.17555411469012672, -0.0,
    -0.21837624155605259, 0.20882905762052142, -0.0, -0.24827531258303759,
    0.24322221611795444, -0.0, -0.27785886189665859, 0.27872286747902364, -0.0,
    -0.30712977649354861, 0.31532038675192614, -0.0, -0.3360909171025292,
    0.35300424584718826, -0.0, -0.36474511842355545, 0.39176401265648464, -0.0,
    -0.39309518936448734, 0.431589350179473, -0.0, -0.42114391327570738,
    0.47247001565857133, -0.0, -0.44889404818260403, 0.51439585972160617, -0.0,
    -0.4763483270159406, 0.557356825532259, -0.0, -0.50350945784012813,
    0.60134294794824061, -0.0, -0.53038012407942237, 0.646344352687124, -0.0,
    -0.55696298474206274, 0.69235125549976384, -0.0, -0.58326067464237263,
    0.739353961351236, -0.0, 0.032190816623613963, -0.026334884546125613, 0.0,
    0.0640451465215906, -0.053854424500504594, 0.0, 0.095566067442515334,
    -0.082547291513030641, 0.0, 0.1267566291319891, -0.11240226049489649, 0.0,
    0.15761985358735717, -0.14340820867920523, 0.0, 0.18815873531012073,
    -0.17555411469012672, 0.0, 0.21837624155605259, -0.20882905762052142, 0.0,
    0.24827531258303759, -0.24322221611795444, 0.0, 0.27785886189665859,
    -0.27872286747902364, 0.0, 0.30712977649354861, -0.31532038675192614, 0.0,
    0.3360909171025292, -0.35300424584718826, 0.0, 0.36474511842355545,
    -0.39176401265648464, 0.0, 0.39309518936448734, -0.431589350179473, 0.0,
    0.42114391327570738, -0.47247001565857133, 0.0, 0.44889404818260403,
    -0.51439585972160617, 0.0, 0.4763483270159406, -0.557356825532259, 0.0,
    0.50350945784012813, -0.60134294794824061, 0.0, 0.53038012407942237,
    -0.646344352687124, 0.0, 0.55696298474206274, -0.69235125549976384, 0.0,
    0.58326067464237263, -0.739353961351236, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T a[498]{ -0.22599754912864506, -0.0, -0.0, -0.0,
    -0.45147579735946186, -0.0, -0.0, -0.0056499387282161265, -0.676435937950623,
    -0.0, -0.0, -0.016936833662202673, -0.90087916141841273, -0.0, -0.0,
    -0.033847732110968247, -1.1248066555435274, -0.0, -0.0,
    -0.056369711146428567, -1.3482196053773612, -0.0, -0.0,
    -0.084489877535016744, -1.5711191932482778, -0.0, -0.0, -0.11819536766945077,
    -1.7935065987678671, -0.0, -0.0, -0.15747334750065772, -2.0153829988371879,
    -0.0, -0.0, -0.20231101246985439, -2.2367495676529967, -0.0, -0.0,
    -0.25269558744078408, -2.4576074767139611, -0.0, -0.0, -0.30861432663210897,
    -2.67795789482686, -0.0, -0.0, -0.370054513549958, -2.8978019881127679, -0.0,
    -0.0, -0.43700346092062947, -3.1171409200132283, -0.0, -0.0,
    -0.50944851062344865, -3.3359758512964088, -0.0, -0.0, -0.58737703362377935,
    -3.5543079400632447, -0.0, -0.0, -0.67077642990618958, -3.7721383417535681,
    -0.0, -0.0, -0.75963412840777078, -3.9894682091522222, -0.0, -0.0,
    -0.85393758695161, -4.2062986923951629, -0.0, -0.0, -0.95367429218041566,
    -4.4226309389755434, -0.0, -0.0, -1.0588317594902947, 0.22599754912864506,
    0.0, 0.0, 0.0, 0.45147579735946186, 0.0, 0.0, 0.0056499387282161265,
    0.676435937950623, 0.0, 0.0, 0.016936833662202673, 0.90087916141841273, 0.0,
    0.0, 0.033847732110968247, 1.1248066555435274, 0.0, 0.0,
    0.056369711146428567, 1.3482196053773612, 0.0, 0.0, 0.084489877535016744,
    1.5711191932482778, 0.0, 0.0, 0.11819536766945077, 1.7935065987678671, 0.0,
    0.0, 0.15747334750065772, 2.0153829988371879, 0.0, 0.0, 0.20231101246985439,
    2.2367495676529967, 0.0, 0.0, 0.25269558744078408, 2.4576074767139611, 0.0,
    0.0, 0.30861432663210897, 2.67795789482686, 0.0, 0.0, 0.370054513549958,
    2.8978019881127679, 0.0, 0.0, 0.43700346092062947, 3.1171409200132283, 0.0,
    0.0, 0.50944851062344865, 3.3359758512964088, 0.0, 0.0, 0.58737703362377935,
    3.5543079400632447, 0.0, 0.0, 0.67077642990618958, 3.7721383417535681, 0.0,
    0.0, 0.75963412840777078, 3.9894682091522222, 0.0, 0.0, 0.85393758695161,
    4.2062986923951629, 0.0, 0.0, 0.95367429218041566, 4.4226309389755434, 0.0,
    0.0, 1.0588317594902947, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.087367622049086546, -0.0, -0.0, -0.0, 0.17453448933447702, -0.0, -0.0,
    0.0021841905512271637, 0.26150106315378036, -0.0, -0.0, 0.00654755278458909,
    0.34826780374462818, -0.0, -0.0, 0.013085079363433599, 0.4348351702871105,
    -0.0, -0.0, 0.0217917744570493, 0.52120362090620576, -0.0, -0.0,
    0.032662653714227059, 0.6073736126742052, -0.0, -0.0, 0.0456927442368822,
    0.69334560161313186, -0.0, -0.0, 0.060877084553737332, 0.779120042697154,
    -0.0, -0.0, 0.078210724594065625, 0.86469738985499234, -0.0, -0.0,
    0.097688725661494474, 0.95007809597232318, -0.0, -0.0, 0.11930616040786929,
    1.0352626128941742, -0.0, -0.0, 0.14305811280717737, 1.1202513914273164,
    -0.0, -0.0, 0.16893967812953173, 1.2050448813426495, -0.0, -0.0,
    0.19694596291521463, 1.2896435313775823, -0.0, -0.0, 0.22707208494878087,
    1.3740477892384069, -0.0, -0.0, 0.25931317323322045, 1.458258101602669, -0.0,
    -0.0, 0.29366436796418061, 1.5422749141215313, -0.0, -0.0,
    0.33012082050424735, 1.6260986714221315, -0.0, -0.0, 0.36867769335728562,
    1.709729817109936, -0.0, -0.0, 0.40933016014283891, -0.087367622049086546,
    0.0, 0.0, 0.0, -0.17453448933447702, 0.0, 0.0, -0.0021841905512271637,
    -0.26150106315378036, 0.0, 0.0, -0.00654755278458909, -0.34826780374462818,
    0.0, 0.0, -0.013085079363433599, -0.4348351702871105, 0.0, 0.0,
    -0.0217917744570493, -0.52120362090620576, 0.0, 0.0, -0.032662653714227059,
    -0.6073736126742052, 0.0, 0.0, -0.0456927442368822, -0.69334560161313186,
    0.0, 0.0, -0.060877084553737332, -0.779120042697154, 0.0, 0.0,
    -0.078210724594065625, -0.86469738985499234, 0.0, 0.0, -0.097688725661494474,
    -0.95007809597232318, 0.0, 0.0, -0.11930616040786929, -1.0352626128941742,
    0.0, 0.0, -0.14305811280717737, -1.1202513914273164, 0.0, 0.0,
    -0.16893967812953173, -1.2050448813426495, 0.0, 0.0, -0.19694596291521463,
    -1.2896435313775823, 0.0, 0.0, -0.22707208494878087, -1.3740477892384069,
    0.0, 0.0, -0.25931317323322045, -1.458258101602669, 0.0, 0.0,
    -0.29366436796418061, -1.5422749141215313, 0.0, 0.0, -0.33012082050424735,
    -1.6260986714221315, 0.0, 0.0, -0.36867769335728562, -1.709729817109936, 0.0,
    0.0, -0.40933016014283891, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.087367622049086546, -0.0, -0.0, -0.0, 0.17453448933447702, -0.0, -0.0,
    0.0021841905512271637, 0.26150106315378036, -0.0, -0.0, 0.00654755278458909,
    0.34826780374462818, -0.0, -0.0, 0.013085079363433599, 0.4348351702871105,
    -0.0, -0.0, 0.0217917744570493, 0.52120362090620576, -0.0, -0.0,
    0.032662653714227059, 0.6073736126742052, -0.0, -0.0, 0.0456927442368822,
    0.69334560161313186, -0.0, -0.0, 0.060877084553737332, 0.779120042697154,
    -0.0, -0.0, 0.078210724594065625, 0.86469738985499234, -0.0, -0.0,
    0.097688725661494474, 0.95007809597232318, -0.0, -0.0, 0.11930616040786929,
    1.0352626128941742, -0.0, -0.0, 0.14305811280717737, 1.1202513914273164,
    -0.0, -0.0, 0.16893967812953173, 1.2050448813426495, -0.0, -0.0,
    0.19694596291521463, 1.2896435313775823, -0.0, -0.0, 0.22707208494878087,
    1.3740477892384069, -0.0, -0.0, 0.25931317323322045, 1.458258101602669, -0.0,
    -0.0, 0.29366436796418061, 1.5422749141215313, -0.0, -0.0,
    0.33012082050424735, 1.6260986714221315, -0.0, -0.0, 0.36867769335728562,
    1.709729817109936, -0.0, -0.0, 0.40933016014283891, -0.087367622049086546,
    0.0, 0.0, 0.0, -0.17453448933447702, 0.0, 0.0, -0.0021841905512271637,
    -0.26150106315378036, 0.0, 0.0, -0.00654755278458909, -0.34826780374462818,
    0.0, 0.0, -0.013085079363433599, -0.4348351702871105, 0.0, 0.0,
    -0.0217917744570493, -0.52120362090620576, 0.0, 0.0, -0.032662653714227059,
    -0.6073736126742052, 0.0, 0.0, -0.0456927442368822, -0.69334560161313186,
    0.0, 0.0, -0.060877084553737332, -0.779120042697154, 0.0, 0.0,
    -0.078210724594065625, -0.86469738985499234, 0.0, 0.0, -0.097688725661494474,
    -0.95007809597232318, 0.0, 0.0, -0.11930616040786929, -1.0352626128941742,
    0.0, 0.0, -0.14305811280717737, -1.1202513914273164, 0.0, 0.0,
    -0.16893967812953173, -1.2050448813426495, 0.0, 0.0, -0.19694596291521463,
    -1.2896435313775823, 0.0, 0.0, -0.22707208494878087, -1.3740477892384069,
    0.0, 0.0, -0.25931317323322045, -1.458258101602669, 0.0, 0.0,
    -0.29366436796418061, -1.5422749141215313, 0.0, 0.0, -0.33012082050424735,
    -1.6260986714221315, 0.0, 0.0, -0.36867769335728562, -1.709729817109936, 0.0,
    0.0, -0.40933016014283891, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T a_1[378]{ -0.025118718490858984, -0.016819265413734148,
    -0.0, -0.0499861237076248, -0.032741132441156645, -0.0,
    -0.074604510545374814, -0.047774052825594666, -0.0, -0.0989761530204811,
    -0.0619264013180454, -0.0, -0.12310330446053273, -0.075206476377572148, -0.0,
    -0.14698819769253046, -0.087622500865329314, -0.0, -0.17063304522936948,
    -0.099182622732274267, -0.0, -0.19404003945462597, -0.1098949157006235, -0.0,
    -0.21721135280566262, -0.11976737993910992, -0.0, -0.24014913795506856,
    -0.12880794273209781, -0.0, -0.26285552799044903, -0.13702445914261108, -0.0,
    -0.28533263659257935, -0.14442471266933057, -0.0, -0.30758255821193842,
    -0.1510164158976148, -0.0, -0.32960736824363657, -0.15680721114459886, -0.0,
    -0.35140912320075196, -0.16180467109842528, -0.0, -0.37298986088609037,
    -0.16601629945166, -0.0, -0.39435160056238255, -0.16944953152894673, -0.0,
    -0.41549634312093342, -0.17211173490895176, -0.0, -0.43642607124873706,
    -0.17401021004065131, -0.0, -0.45714274959407175, -0.17515219085401298, -0.0,
    0.025118718490858984, 0.016819265413734148, 0.0, 0.0499861237076248,
    0.032741132441156645, 0.0, 0.074604510545374814, 0.047774052825594666, 0.0,
    0.0989761530204811, 0.0619264013180454, 0.0, 0.12310330446053273,
    0.075206476377572148, 0.0, 0.14698819769253046, 0.087622500865329314, 0.0,
    0.17063304522936948, 0.099182622732274267, 0.0, 0.19404003945462597,
    0.1098949157006235, 0.0, 0.21721135280566262, 0.11976737993910992, 0.0,
    0.24014913795506856, 0.12880794273209781, 0.0, 0.26285552799044903,
    0.13702445914261108, 0.0, 0.28533263659257935, 0.14442471266933057, 0.0,
    0.30758255821193842, 0.1510164158976148, 0.0, 0.32960736824363657,
    0.15680721114459886, 0.0, 0.35140912320075196, 0.16180467109842528, 0.0,
    0.37298986088609037, 0.16601629945166, 0.0, 0.39435160056238255,
    0.16944953152894673, 0.0, 0.41549634312093342, 0.17211173490895176, 0.0,
    0.43642607124873706, 0.17401021004065131, 0.0, 0.45714274959407175,
    0.17515219085401298, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.0189087824026894, -0.066906613140151408, -0.0, 0.03760442849307085,
    -0.13454624955425293, -0.0, 0.056088892993839096, -0.20291172113658829, -0.0,
    0.0743641128405729, -0.27199590536902313, -0.0, 0.092432007343536085,
    -0.34179174472431623, -0.0, 0.11029447834800676, -0.41229224607485848, -0.0,
    0.12795341039314814, -0.48349048010679008, -0.0, 0.14541067086943421,
    -0.55537958073944738, -0.0, 0.16266811017464336, -0.62795274455009054, -0.0,
    0.17972756186843311, -0.70120323020386432, -0.0, 0.19659084282550868,
    -0.775124357888944, -0.0, 0.21325975338739828, -0.84970950875681928, -0.0,
    0.22973607751284789, -0.9249521243676706, -0.0, 0.24602158292684784,
    -1.0008457061407889, -0.0, 0.26211802126830397, -1.0773838148099963, -0.0,
    0.27802712823636522, -1.1545600698840195, -0.0, 0.2937506237354206,
    -1.2323681491117719, -0.0, 0.30929021201877682, -1.3108017879525005, -0.0,
    0.32464758183102904, -1.3898547790507516, -0.0, 0.33982440654913676,
    -1.4695209717161131, -0.0, -0.0189087824026894, 0.066906613140151408, 0.0,
    -0.03760442849307085, 0.13454624955425293, 0.0, -0.056088892993839096,
    0.20291172113658829, 0.0, -0.0743641128405729, 0.27199590536902313, 0.0,
    -0.092432007343536085, 0.34179174472431623, 0.0, -0.11029447834800676,
    0.41229224607485848, 0.0, -0.12795341039314814, 0.48349048010679008, 0.0,
    -0.14541067086943421, 0.55537958073944738, 0.0, -0.16266811017464336,
    0.62795274455009054, 0.0, -0.17972756186843311, 0.70120323020386432, 0.0,
    -0.19659084282550868, 0.775124357888944, 0.0, -0.21325975338739828,
    0.84970950875681928, 0.0, -0.22973607751284789, 0.9249521243676706, 0.0,
    -0.24602158292684784, 1.0008457061407889, 0.0, -0.26211802126830397,
    1.0773838148099963, 0.0, -0.27802712823636522, 1.1545600698840195, 0.0,
    -0.2937506237354206, 1.2323681491117719, 0.0, -0.30929021201877682,
    1.3108017879525005, 0.0, -0.32464758183102904, 1.3898547790507516, 0.0,
    -0.33982440654913676, 1.4695209717161131, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0,
    0.0, -0.032190816623613963, 0.026334884546125613, -0.0, -0.0640451465215906,
    0.053854424500504594, -0.0, -0.095566067442515334, 0.082547291513030641,
    -0.0, -0.1267566291319891, 0.11240226049489649, -0.0, -0.15761985358735717,
    0.14340820867920523, -0.0, -0.18815873531012073, 0.17555411469012672, -0.0,
    -0.21837624155605259, 0.20882905762052142, -0.0, -0.24827531258303759,
    0.24322221611795444, -0.0, -0.27785886189665859, 0.27872286747902364, -0.0,
    -0.30712977649354861, 0.31532038675192614, -0.0, -0.3360909171025292,
    0.35300424584718826, -0.0, -0.36474511842355545, 0.39176401265648464, -0.0,
    -0.39309518936448734, 0.431589350179473, -0.0, -0.42114391327570738,
    0.47247001565857133, -0.0, -0.44889404818260403, 0.51439585972160617, -0.0,
    -0.4763483270159406, 0.557356825532259, -0.0, -0.50350945784012813,
    0.60134294794824061, -0.0, -0.53038012407942237, 0.646344352687124, -0.0,
    -0.55696298474206274, 0.69235125549976384, -0.0, -0.58326067464237263,
    0.739353961351236, -0.0, 0.032190816623613963, -0.026334884546125613, 0.0,
    0.0640451465215906, -0.053854424500504594, 0.0, 0.095566067442515334,
    -0.082547291513030641, 0.0, 0.1267566291319891, -0.11240226049489649, 0.0,
    0.15761985358735717, -0.14340820867920523, 0.0, 0.18815873531012073,
    -0.17555411469012672, 0.0, 0.21837624155605259, -0.20882905762052142, 0.0,
    0.24827531258303759, -0.24322221611795444, 0.0, 0.27785886189665859,
    -0.27872286747902364, 0.0, 0.30712977649354861, -0.31532038675192614, 0.0,
    0.3360909171025292, -0.35300424584718826, 0.0, 0.36474511842355545,
    -0.39176401265648464, 0.0, 0.39309518936448734, -0.431589350179473, 0.0,
    0.42114391327570738, -0.47247001565857133, 0.0, 0.44889404818260403,
    -0.51439585972160617, 0.0, 0.4763483270159406, -0.557356825532259, 0.0,
    0.50350945784012813, -0.60134294794824061, 0.0, 0.53038012407942237,
    -0.646344352687124, 0.0, 0.55696298474206274, -0.69235125549976384, 0.0,
    0.58326067464237263, -0.739353961351236, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0,
    1.0 };

  static const real_T b_Kr_0[300]{ -0.0, 0.00013063329054594013,
    0.00022655638891913061, -0.0, -0.0, -0.0, 0.00026126801630562953,
    0.00045275817202038512, 0.0, 0.0, -0.0, 0.00039190420028970463,
    0.00067860555831980476, 0.0, 0.0, -0.0, 0.000522541865498165,
    0.00090409875667013951, 0.0, 0.0, -0.0, 0.00065318103492038263,
    0.0011292379757609637, 0.0, 0.0, -0.0, 0.00078382173153511111,
    0.0013540234241187911, 0.0, 0.0, -0.0, 0.0009144639783104947,
    0.0015784553101071898, 0.0, 0.0, -0.0, 0.0010451077982040777,
    0.0018025338419268967, 0.0, 0.0, -0.0, 0.0011757532141628137,
    0.0020262592276159342, 0.0, 0.0, -0.0, 0.0013064002491230743,
    0.0022496316750497225, 0.0, 0.0, -0.0, 0.0014370489260106588,
    0.0024726513919411971, 0.0, 0.0, -0.0, 0.0015676992677408031,
    0.00269531858584092, 0.0, 0.0, -0.0, 0.0016983512972181892,
    0.0029176334641371981, 0.0, 0.0, -0.0, 0.0018290050373369539,
    0.0031395962340561947, 0.0, 0.0, -0.0, 0.0019596605109806976,
    0.0033612071026620453, 0.0, 0.0, -0.0, 0.0020903177410224955,
    0.0035824662768569717, 0.0, 0.0, -0.0, 0.0022209767503249033,
    0.0038033739633813964, 0.0, 0.0, -0.0, 0.0023516375617399694,
    0.0040239303688140559, 0.0, 0.0, -0.0, 0.0024823001981092423,
    0.004244135699572116, 0.0, 0.0, -0.0, 0.0026129646822637812,
    0.004463990161911284, 0.0, 0.0, -0.0, -0.0025019353978715038,
    0.0011638312285314027, -0.0, -0.0, -0.0, -0.005004252640033722,
    0.0023311729688525766, -0.0, 0.0, -0.0, -0.0075069519993825069,
    0.0035020237695565378, -0.0, 0.0, -0.0, -0.01001003374875372,
    0.0046763821805602748, -0.0, 0.0, -0.0, -0.012513498160923319,
    0.0058542467531038558, -0.0, 0.0, -0.0, -0.015017345508607423,
    0.0070356160397495393, -0.0, 0.0, -0.0, -0.0175215760644624,
    0.008220488594380888, -0.0, 0.0, -0.0, -0.020026190101084949,
    0.0094088629722018768, -0.0, 0.0, -0.0, -0.022531187891012154,
    0.010600737729736007, -0.0, 0.0, -0.0, -0.025036569706721589,
    0.011796111424825416, -0.0, 0.0, -0.0, -0.027542335820631394,
    0.012994982616630004, -0.0, 0.0, -0.0, -0.030048486505100322,
    0.01419734986562653, -0.0, 0.0, -0.0, -0.032555022032427855,
    0.015403211733607739, -0.0, 0.0, -0.0, -0.035061942674854255,
    0.016612566783681476, -0.0, 0.0, -0.0, -0.037569248704560652,
    0.017825413580269792, -0.0, 0.0, -0.0, -0.040076940393669123,
    0.019041750689108086, -0.0, 0.0, -0.0, -0.042585018014242769,
    0.020261576677244188, -0.0, 0.0, -0.0, -0.045093481838285779,
    0.021484890113037513, -0.0, 0.0, -0.0, -0.047602332137743536,
    0.022711689566158148, -0.0, 0.0, -0.0, -0.050111569184502665,
    0.023941973607586004, -0.0, 0.0, -0.0, 0.00084054322877011171,
    -0.001914736296470494, -0.0, -0.0, -0.0, 0.0016813128609426524,
    -0.0038297434705267387, 0.0, -0.0, -0.0, 0.00252230894173239,
    -0.0057450217411265277, 0.0, -0.0, -0.0, 0.0033635315163737392,
    -0.0076605713271730861, 0.0, -0.0, -0.0, 0.004204980630120759,
    -0.0095763924475151373, 0.0, -0.0, -0.0, 0.0050466563282471526,
    -0.011492485320946967, 0.0, -0.0, -0.0, 0.0058885586560462645,
    -0.013408850166208491, 0.0, -0.0, -0.0, 0.0067306876588310825,
    -0.015325487201985316, 0.0, -0.0, -0.0, 0.007573043381934234,
    -0.017242396646908811, 0.0, -0.0, -0.0, 0.0084156258707079844,
    -0.019159578719556172, 0.0, -0.0, -0.0, 0.00925843517052424,
    -0.021077033638450482, 0.0, -0.0, -0.0, 0.01010147132677454,
    -0.022994761622060773, 0.0, -0.0, -0.0, 0.010944734384870066,
    -0.02491276288880211, 0.0, -0.0, -0.0, 0.011788224390241631,
    -0.026831037657035636, 0.0, -0.0, -0.0, 0.012631941388339681,
    -0.028749586145068649, 0.0, -0.0, -0.0, 0.013475885424634295,
    -0.030668408571154656, 0.0, -0.0, -0.0, 0.014320056544615187,
    -0.03258750515349345, 0.0, -0.0, -0.0, 0.015164454793791698,
    -0.034506876110231177, 0.0, -0.0, -0.0, 0.016009080217692803,
    -0.036426521659460379, 0.0, -0.0, -0.0, 0.016853932861867104,
    -0.038346442019220085, 0.0, -0.0 };

  static const real_T b_Kr_1[300]{ -0.0004600974465772233,
    -0.00030807706503736739, -0.0, -0.0, -0.0, -0.00091559160912372238,
    -0.00059971655957308809, -0.0, -1.1502436164430583E-5,
    -7.7019266259341861E-6, -0.0013665245230389403, -0.00087507329347635886,
    -0.0, -3.4392226392523637E-5, -2.269484061526139E-5, -0.0018129378412888371,
    -0.0011343006663543662, -0.0, -6.8555339468497146E-5, -4.4571672952170365E-5,
    -0.0022548728378846811, -0.0013775506803813825, -0.0,
    -0.00011387878550071806, -7.29291896110295E-5, -0.0026923704113301947,
    -0.0016049739530111637, -0.0, -0.00017025060644783508,
    -0.00010736795662056407, -0.0031254710880373429, -0.0018167197295737088,
    -0.0, -0.00023755986673108996, -0.00014749230544584317,
    -0.0035542150257110534, -0.0020129358957574352, -0.0, -0.0003156966439320235,
    -0.00019291029868518587, -0.0039786420167031451, -0.0021937689899778094,
    -0.0, -0.00040455201957479987, -0.00024323369607912176,
    -0.0043987914913357464, -0.0023593642156334686, -0.0,
    -0.00050401806999237839, -0.000298077920828567, -0.0048147025211944918,
    -0.0025098654532508553, -0.0, -0.000613987857275772, -0.00035706202621940371,
    -0.0052264138223917468, -0.0026454152725183784, -0.0,
    -0.00073435542030563428, -0.00041980866255067505, -0.0056339637588001621,
    -0.002766154944211108, -0.0, -0.00086501576586542795,
    -0.00048594404436363449, -0.00603739034525681, -0.0028722244520069961, -0.0,
    -0.001005864859835432, -0.0005550979179689122, -0.0064367312507381744,
    -0.0029637625041956145, -0.0, -0.0011567996184668524,
    -0.00062690352926908716, -0.0068320238015062576, -0.0030409065452803819,
    -0.0, -0.0013177178997353068, -0.00070099759187397753,
    -0.0072233049842260729, -0.0031037927674752558, -0.0, -0.0014885184947729634,
    -0.00077702025550598715, -0.0076106114490547766, -0.0031525561220968426,
    -0.0, -0.0016691011193786152, -0.00085461507469286851,
    -0.0079939795127027015, -0.0031873303308528824, -0.0, -0.0018593664056049846,
    -0.00093342897774528961, -0.0083734451614665452, -0.0032082478970280473,
    -0.0, -0.0020592158934225522, -0.0010131122360166117, 0.00034635057136882692,
    -0.0012255227859700041, -0.0, -0.0, -0.0, 0.00068879714289381386,
    -0.0024644722973817881, -0.0, 8.6587642842206746E-6, -3.0638069649250106E-5,
    0.0010273755190655818, -0.0037167168703096231, -0.0, 2.587869285656602E-5,
    -9.2249877083794811E-5, 0.00136212117856965, -0.0049821260421899823, -0.0,
    5.1563080833205561E-5, -0.00018516779884153537, 0.0016930692772501304,
    -0.0062605705408920447, -0.0, 8.56161102974468E-5, -0.000309720949896285,
    0.0020202546510464671, -0.0075519222738876131, -0.0, 0.00012794284222870006,
    -0.00046623521341858609, 0.0023437118189034568, -0.0088560543175195509, -0.0,
    0.00017844920850486172, -0.00065503327026577641, 0.0026634749856548046,
    -0.010172840906367841, -0.0, 0.00023704200397744817, -0.00087643462820376531,
    0.0029795780448804432, -0.011502157422712366, -0.0, 0.00030362887861881823,
    -0.0011307556508629614, 0.0032920545817378671, -0.012843880386091551, -0.0,
    0.0003781183297408293, -0.0014183095864307705, 0.0036009378757677073,
    -0.014197887442955979, -0.0, 0.00046041969428427596, -0.001739406596083059,
    0.0039062609036737877, -0.015564057356416118, -0.0, 0.00055044314117846862,
    -0.0020943537821569583, 0.004208056342077892, -0.016942269996083323, -0.0,
    0.00064809966377031331, -0.0024834552160673613, 0.0045063565702494722,
    -0.018332406328003231, -0.0, 0.00075330107232226065, -0.0029070119659694443,
    0.0048011936728105283, -0.019734348404680729, -0.0, 0.00086595998657849741,
    -0.0033653221241695254, 0.005092599442415878, -0.021147979355195683, -0.0,
    0.00098598982839876066, -0.0038586808342865437, 0.0053806053824090551,
    -0.022573183375408538, -0.0, 0.0011133048144591577, -0.0043873803181664357,
    0.00566524270945404, -0.024009845718255044, -0.0, 0.001247819949019384,
    -0.0049517099025516492, 0.0059465423561430517, -0.02545785268412926, -0.0,
    0.0013894510167557351, -0.0055519560455080256, 0.0062245349735806244,
    -0.026917091611354027, -0.0, 0.0015381145756593114, -0.006188402362611257,
    -0.00058963647118981388, 0.00048237385756712709, -0.0, -0.0, -0.0,
    -0.0011731095434256084, 0.00098644694826231457, -0.0, -1.4740911779745348E-5,
    1.2059346439178178E-5, -0.0017504755915684609, 0.0015120117716527783, -0.0,
    -4.4068650365385555E-5, 3.6720520145736047E-5, -0.0023217904775511451,
    0.002058862718734233, -0.0, -8.7830540154597067E-5, 7.4520814437055514E-5,
    -0.0028871095550439779, 0.0026267960547241894, -0.0, -0.00014587530209337569,
    0.00012599238240541136, -0.0034464876740782206, 0.0032156099020117707, -0.0,
    -0.00021805304096947513, 0.00019166228377351607, -0.0039999791856274274,
    0.0038251042232626273, -0.0, -0.00030421523282143063, 0.00027205253132381031,
    -0.0045476379461471166, 0.0044550808046775305, -0.0, -0.00040421471246211634,
    0.000367680136905376, -0.00508951732207315, 0.00510534323940326, -0.0,
    -0.0005179056611157942, 0.00047905715702231428, -0.0056256701942791927,
    0.0057756969110943886, -0.0, -0.00064514359416762291, 0.00060669073800739587,
    -0.0061561489624936232, 0.0064659489776245967, -0.0, -0.00078578534902460273,
    0.00075108316078475556, -0.0066810055496762692, 0.00717590835494616, -0.0,
    -0.00093968907308694335, 0.0009127318852253704, -0.0072002914063553315,
    0.00790538570109625, -0.0, -0.00110671421182885, 0.0010921295940990244,
    -0.0077140575149248509, 0.0086541934003487274, -0.0, -0.0012867214969877333,
    0.0012897642366264308, -0.00822235439390309, 0.0094221455475100978, -0.0,
    -0.0014795729348608545, 0.001506119071635149, -0.008725232102152182,
    0.010209057932358311, -0.0, -0.0016851317947084317, 0.0017416727103229015,
    -0.0092227402430593835, 0.011014748024223114, -0.0, -0.0019032625972622364,
    0.0019968991586318593, -0.009714927968680306, 0.011839034956706681, -0.0,
    -0.0021338311033387214, 0.0022722678592374374, -0.010201843983844449,
    0.012681739512543208, -0.0, -0.0023767043025557291, 0.0025682437331551043,
    -0.010683536550223399, 0.013542684108596249, -0.0, -0.0026317504021518404,
    0.0028852872209686846 };

  static const real_T b_Kr[240]{ -0.004139578033196244, -0.0, -0.0, -0.0,
    -0.0082696440756759651, -0.0, -0.0, -0.00010348945082990611,
    -0.012390219984248396, -0.0, -0.0, -0.00031023055272180522,
    -0.016501327565499841, -0.0, -0.0, -0.00061998605232801507,
    -0.020602988575909063, -0.0, -0.0, -0.0010325192414655112,
    -0.024695224721962426, -0.0, -0.0, -0.0015475939558632374,
    -0.02877805766026879, -0.0, -0.0, -0.002164974573912298,
    -0.032851508997674084, -0.0, -0.0, -0.0028844260154190179,
    -0.03691560029137566, -0.0, -0.0, -0.0037057137403608696,
    -0.040970353049036422, -0.0, -0.0, -0.004628603747645261,
    -0.045015788728898581, -0.0, -0.0, -0.0056528625738711718,
    -0.049051928739897263, -0.0, -0.0, -0.0067782572920936359,
    -0.053078794441773762, -0.0, -0.0, -0.0080045555105910677,
    -0.057096407145188648, -0.0, -0.0, -0.00933152537163541,
    -0.061104788111834459, -0.0, -0.0, -0.010758935550265128,
    -0.0651039585545483, -0.0, -0.0, -0.012286555253060988,
    -0.069093939637424051, -0.0, -0.0, -0.013914154216924698,
    -0.07307475247592439, -0.0, -0.0, -0.0156415027078603, -0.077046418136992548,
    -0.0, -0.0, -0.01746837151975841, -0.08100895763916377, -0.0, -0.0,
    -0.019394531973183225, 0.0016003053592458171, -0.0, -0.0, -0.0,
    0.0031969335104287152, -0.0, -0.0, 4.0007633981145425E-5,
    0.0047898929030980638, -0.0, -0.0, 0.00011993097174186333,
    0.006379191967387717, -0.0, -0.0, 0.00023967829431931491,
    0.00796483911406063, -0.0, -0.0, 0.0003991580935040078,
    0.0095468427345533641, -0.0, -0.0, 0.00059827907135552344,
    0.011125211201020505, -0.0, -0.0, 0.00083695013971935754,
    0.012699952866378959, -0.0, -0.0, 0.0011150804197448702,
    0.014271076064352163, -0.0, -0.0, 0.001432579241404344, 0.01583858910951419,
    -0.0, -0.0, 0.0017893561430131482, 0.017402500297333742, -0.0, -0.0,
    0.0021853208707510031, 0.018962817904218058, -0.0, -0.0,
    0.0026203833781843467, 0.020519550187556714, -0.0, -0.0,
    0.0030944538257897981, 0.022072705385765324, -0.0, -0.0,
    0.0036074425804787156, 0.023622291718329122, -0.0, -0.0,
    0.0041592602151228489, 0.025168317385846482, -0.0, -0.0,
    0.0047498175080810771, 0.02671079057007231, -0.0, -0.0,
    0.0053790254427272394, 0.028249719433961332, -0.0, -0.0,
    0.0060467952069790474, 0.02978511212171131, -0.0, -0.0, 0.00675303819282808,
    0.031316976758806125, -0.0, -0.0, 0.0074976659958708629,
    0.0016003053592458171, -0.0, -0.0, -0.0, 0.0031969335104287152, -0.0, -0.0,
    4.0007633981145425E-5, 0.0047898929030980638, -0.0, -0.0,
    0.00011993097174186333, 0.006379191967387717, -0.0, -0.0,
    0.00023967829431931491, 0.00796483911406063, -0.0, -0.0,
    0.0003991580935040078, 0.0095468427345533641, -0.0, -0.0,
    0.00059827907135552344, 0.011125211201020505, -0.0, -0.0,
    0.00083695013971935754, 0.012699952866378959, -0.0, -0.0,
    0.0011150804197448702, 0.014271076064352163, -0.0, -0.0,
    0.001432579241404344, 0.01583858910951419, -0.0, -0.0, 0.0017893561430131482,
    0.017402500297333742, -0.0, -0.0, 0.0021853208707510031,
    0.018962817904218058, -0.0, -0.0, 0.0026203833781843467,
    0.020519550187556714, -0.0, -0.0, 0.0030944538257897981,
    0.022072705385765324, -0.0, -0.0, 0.0036074425804787156,
    0.023622291718329122, -0.0, -0.0, 0.0041592602151228489,
    0.025168317385846482, -0.0, -0.0, 0.0047498175080810771, 0.02671079057007231,
    -0.0, -0.0, 0.0053790254427272394, 0.028249719433961332, -0.0, -0.0,
    0.0060467952069790474, 0.02978511212171131, -0.0, -0.0, 0.00675303819282808,
    0.031316976758806125, -0.0, -0.0, 0.0074976659958708629 };

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

  static const real_T b_Kx_0[21]{ -0.026404598946573047, -0.046701198280547959,
    0.0, 0.0, -0.027435037351643085, -0.047106463579682338, 0.0,
    0.53251906053640408, -0.24750295164334835, 0.0, 0.0, 0.52591591852487163,
    -0.24906487889444787, 0.0, -0.18632686544532984, 0.399349755453896, 0.0, 0.0,
    -0.17681546557705763, 0.40245619880339512, 0.0 };

  static const real_T b_Kx_1[21]{ 0.064401756957438827, 0.047455061923499582,
    0.014630349328867837, 0.0078609701236701246, 0.090742080898368632,
    0.043732737337692507, 0.0, -0.19527152302960851, 0.2933282397767637,
    -0.01094980165376826, 0.042692224192222988, -0.067749117999953076,
    0.27445318611580427, 0.0, 0.17088341788856595, -0.13723068591963278,
    0.018711210911735485, -0.01918263454114651, 0.115953552636297,
    -0.12895417294734363, 0.0 };

  static const real_T b_Hinv[16]{ 6.9623302122049981, 8.6381748287066245,
    8.6381748287066245, 0.0, 8.6381748287066245, 25.967662402148488,
    -3.3394025578514817, 0.0, 8.6381748287066245, -3.3394025578514817,
    25.967662402148491, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Hinv_0[16]{ 28.12258982300947, 2.763469259946183,
    4.5879287660733707, 0.0, 2.763469259946183, 1.9875333798155721,
    2.0875914381912626, 0.0, 4.5879287660733707, 2.0875914381912626,
    3.6968865774238209, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Hinv_1[16]{ 17.34236214914117, -4.5461861246542492,
    -7.1762044834945335, 0.0, -4.5461861246542492, 7.7917586452799936,
    9.183059203964179, 0.0, -7.1762044834945335, 9.183059203964179,
    17.834757857570036, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Linv[16]{ 0.60597727574674687, 1.9291506615478005,
    1.6951399583908207, 0.0, 0.0, 5.0535354185274359, -0.65531837746031507, 0.0,
    0.0, 0.0, 5.0958475646499171, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Linv_0[16]{ 4.7320151847238234, 0.19206168771609869,
    2.3861547356917727, 0.0, 0.0, 0.89927377084443239, 1.0857440144374766, 0.0,
    0.0, 0.0, 1.9227289401847107, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Linv_1[16]{ 3.7707236243550346, -0.48631557236454065,
    -1.6992653300611269, 0.0, 0.0, 1.7502660163283204, 2.1744717789307328, 0.0,
    0.0, 0.0, 4.2231218141997795, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Kx[15]{ 0.88219871842955933, 0.13622878929978996,
    0.85679023656649278, 0.0, 0.0, -0.34104619497474442, -0.0526642232255966,
    -0.33122361659364064, 0.0, 0.0, -0.34104619497474442, -0.0526642232255966,
    -0.33122361659364064, 0.0, 0.0 };

  static const real_T b_Ku1[9]{ 2.6891273495088517, -1.0395805743008937,
    -1.0395805743008937, -1.0395805743008937, 0.40188791009141395,
    0.40188791009141395, -1.0395805743008937, 0.40188791009141395,
    0.40188791009141395 };

  static const real_T b_Ku1_0[9]{ 0.0105374058219662, -0.00953798324465595,
    -0.0500367818799705, -0.0095379832446559482, 1.2044783155925893,
    -0.68758688074338681, -0.0500367818799705, -0.6875868807433867,
    0.68674611351228942 };

  static const real_T b_Ku1_1[9]{ 0.036210150692227161, 0.019541806654652709,
    0.018237449206281651, 0.019541806654652709, 0.29773963560129468,
    -0.16301113588039762, 0.018237449206281654, -0.16301113588039762,
    0.11322096555832134 };

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
      real_T Saturation_idx_0;
      real_T Saturation_idx_1;
      real_T Saturation_idx_2;
      real_T dwt;
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
      Saturation_idx_0 = 2.197 / (0.2 * rtU.k_2);

      // 'wtMod_:10' x0 = 0.5*k_2;
      Saturation_idx_1 = 0.5 * rtU.k_2;

      //  midpoint
      // 'wtMod_:11' sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));
      // 'wtMod_:13' for i = 1:2*no
      for (k = 0; k < 6; k++) {
        // Delay: '<S8>/Delay'
        Saturation_idx_2 = rtDW.Delay_DSTATE[k];

        // MATLAB Function: '<S8>/MATLAB Function' incorporates:
        //   Inport: '<Root>/k_2'
        //   Outport: '<Root>/currEv'

        // 'wtMod_:14' if yDest(i) ~= 0
        if (rtY.currEv.r[k] != 0.0) {
          //  drive ywt to 1
          // 'wtMod_:16' if (ywtT(i) <= 1)
          if (Saturation_idx_2 <= 1.0) {
            // 'wtMod_:17' ywtT(i) = ywtT(i) + dwt;
            Saturation_idx_2 += dwt;
          }

          // 'wtMod_:19' else
          //  drive ywt to 0
          // 'wtMod_:21' if (ywtT(i) > 0)
        } else if (Saturation_idx_2 > 0.0) {
          // 'wtMod_:22' ywtT(i) = ywtT(i) - dwt;
          Saturation_idx_2 -= dwt;
        } else {
          // no actions
        }

        // 'wtMod_:25' if ywtT(i) <= 0
        if (Saturation_idx_2 <= 0.0) {
          // 'wtMod_:26' ywt(i) = 0;
          rtb_ywt[k] = 0.0;
        } else {
          // 'wtMod_:27' else
          // 'wtMod_:28' ywt(i) = sigmoid(ywtT(i)*k_2);
          // 'wtMod_:11' @(x) 1/(1 + exp(-k_1*(x-x0)))
          rtb_ywt[k] = 1.0 / (std::exp((Saturation_idx_2 * rtU.k_2 -
            Saturation_idx_1) * -Saturation_idx_0) + 1.0);
        }

        // Delay: '<S8>/Delay'
        rtb_ywtT[k] = Saturation_idx_2;
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
      //   Constant: '<S3>/Constant1'
      //   Constant: '<S3>/Constant13'
      //   Constant: '<S4>/Constant1'
      //   Constant: '<S4>/Constant13'
      //   Constant: '<S5>/Constant1'
      //   Constant: '<S5>/Constant13'
      //   DataTypeConversion: '<S112>/DataTypeConversionEnable'
      //   DataTypeConversion: '<S182>/DataTypeConversionEnable'
      //   DataTypeConversion: '<S252>/DataTypeConversionEnable'
      //   Delay: '<S112>/MemoryX'
      //   Delay: '<S182>/MemoryP'
      //   Delay: '<S182>/MemoryX'
      //   Delay: '<S252>/MemoryP'
      //   Delay: '<S252>/MemoryX'
      //   Outport: '<Root>/yhat'
      //   Product: '<S155>/C[k]*xhat[k|k-1]'
      //   Product: '<S155>/D[k]*u[k]'
      //   Product: '<S155>/Product3'
      //   Product: '<S185>/Product'
      //   Product: '<S185>/Product1'
      //   Product: '<S255>/Product'
      //   Product: '<S255>/Product1'
      //   Sum: '<S155>/Sum'
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
        Saturation_idx_0 = rtDW.last_mv_DSTATE_n[0];
        rtb_TmpSignalConversionAtSFu_o4[2] = rtP.Constant1_Value_j[0] +
          rtDW.MemoryX_DSTATE_l[1];
        Saturation_idx_1 = rtDW.last_mv_DSTATE_n[1];
        rtb_TmpSignalConversionAtSFu_o4[3] = rtP.Constant1_Value_j[1] +
          rtDW.MemoryX_DSTATE_l[2];
        Saturation_idx_2 = rtDW.last_mv_DSTATE_n[2];
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

          Bc_2 = -(((a[k + 166] * Saturation_idx_1 + a[k] * Saturation_idx_0) +
                    a[k + 332] * Saturation_idx_2) + (dwt + Bc_2));
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

          rtb_Sum2_f[k] = ((b_Ku1[3 * k + 1] * Saturation_idx_1 + b_Ku1[3 * k] *
                            Saturation_idx_0) + b_Ku1[3 * k + 2] *
                           Saturation_idx_2) + (dwt + Bc_2);
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

        // Outputs for Atomic SubSystem: '<S112>/CalculatePL'
        // MATLAB Function: '<S114>/Discrete-Time KF - Calculate PLMZ' incorporates:
        //   Constant: '<S3>/Constant1'
        //   DataTypeConversion: '<S112>/DataTypeConversionEnable'
        //   Delay: '<S112>/MemoryP'
        //   Product: '<S132>/Product'
        //   Product: '<S132>/Product2'
        //   Sum: '<S132>/Add1'

        //  See help of ctrlKalmanFilterDTCalculatePL.m
        // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S152>:1' 
        //    Copyright 2014 The MathWorks, Inc.
        // '<S152>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
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
        } else {
          (void)std::memset(&rtb_N_f[0], 0, 12U * sizeof(real_T));
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
        }

        // End of MATLAB Function: '<S114>/Discrete-Time KF - Calculate PLMZ'
        // End of Outputs for SubSystem: '<S112>/CalculatePL'

        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S90>/umin_scale1'

        //  Determine if the Square-Root algorithm was used
        // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S154>:1' 
        // '<S154>:1:4' if isSqrtUsed
        Saturation_idx_0 = rtP.umin_scale1_Gain[0] * umax_incr[0];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = Saturation_idx_0 - rtU.u0[0];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S90>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain[1] * umax_incr[1];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = Saturation_idx_0 - rtU.u0[1];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S90>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain[2] * umax_incr[2];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S89>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = Saturation_idx_0 - rtU.u0[2];

        // Outputs for Enabled SubSystem: '<S131>/MeasurementUpdate' incorporates:
        //   EnablePort: '<S155>/Enable'

        if (rtP.Constant1_Value_e != 0.0) {
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
            // Product: '<S155>/C[k]*xhat[k|k-1]'
            rtb_C_0[k_0] = 0.0;
            rtb_C_0[k_0] += rtb_C_c[k_0] * rtDW.MemoryX_DSTATE_l[0];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 3] * rtDW.MemoryX_DSTATE_l[1];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 6] * rtDW.MemoryX_DSTATE_l[2];
            rtb_C_0[k_0] += rtb_C_c[k_0 + 9] * rtDW.MemoryX_DSTATE_l[3];

            // Product: '<S155>/D[k]*u[k]' incorporates:
            //   Delay: '<S112>/MemoryX'
            //   Product: '<S155>/C[k]*xhat[k|k-1]'

            tmp[k_0] = 0.0;
            tmp[k_0] += rtP.Constant13_Value_c[k_0] * rtb_Sum1[0];
            tmp[k_0] += rtP.Constant13_Value_c[k_0 + 3] * rtb_Sum1[1];
            tmp[k_0] += rtP.Constant13_Value_c[k_0 + 6] * rtb_Sum1[2];

            // Sum: '<S155>/Sum' incorporates:
            //   Constant: '<S3>/Constant13'
            //   Product: '<S155>/C[k]*xhat[k|k-1]'
            //   Product: '<S155>/D[k]*u[k]'
            //   Sum: '<S155>/Add1'

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
        } else if (rtDW.MeasurementUpdate_MODE) {
          // Disable for Product: '<S155>/Product3' incorporates:
          //   Outport: '<S155>/L*(y[k]-yhat[k|k-1])'
          //
          rtDW.Product3_c[0] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3_c[1] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3_c[2] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3_c[3] = rtP.Lykyhatkk1_Y0_c;
          rtDW.MeasurementUpdate_MODE = false;
        } else {
          // no actions
        }

        // End of Outputs for SubSystem: '<S131>/MeasurementUpdate'

        // Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
        //   Constant: '<S3>/Constant1'
        //   Constant: '<S3>/Constant13'
        //   DataTypeConversion: '<S112>/DataTypeConversionEnable'
        //   Delay: '<S112>/MemoryX'
        //   Inport: '<Root>/iRST'
        //   Product: '<S155>/C[k]*xhat[k|k-1]'
        //   Product: '<S155>/D[k]*u[k]'
        //   Product: '<S155>/Product3'
        //   Sum: '<S155>/Sum'
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

        rtY.u[2] = Saturation_idx_0;

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
        Saturation_idx_0 = rtDW.last_mv_DSTATE_i[0];
        Saturation_idx_1 = rtDW.last_mv_DSTATE_i[1];
        Saturation_idx_2 = rtDW.last_mv_DSTATE_i[2];

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

          Bc_2 = -(((a_0[k + 206] * Saturation_idx_1 + a_0[k] * Saturation_idx_0)
                    + a_0[k + 412] * Saturation_idx_2) + (dwt + Bc_2));
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

          rtb_Sum2_f[k] = ((b_Ku1_0[3 * k + 1] * Saturation_idx_1 + b_Ku1_0[3 *
                            k] * Saturation_idx_0) + b_Ku1_0[3 * k + 2] *
                           Saturation_idx_2) + (dwt + Bc_2);
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

        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S160>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_p[0] * umax_incr[0];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = Saturation_idx_0 - rtU.u0[0];

        // Sum: '<S159>/Sum6'
        rtb_Sum6[0] = y__m[0];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S160>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_p[1] * umax_incr[1];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = Saturation_idx_0 - rtU.u0[1];

        // Sum: '<S159>/Sum6'
        rtb_Sum6[1] = y__m[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S160>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_p[2] * umax_incr[2];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S159>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = Saturation_idx_0 - rtU.u0[2];

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

        rtY.u[2] = Saturation_idx_0;

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
        Saturation_idx_0 = rtDW.last_mv_DSTATE[0];
        Saturation_idx_1 = rtDW.last_mv_DSTATE[1];
        Saturation_idx_2 = rtDW.last_mv_DSTATE[2];

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

          Bc_2 = -(((a_1[k + 126] * Saturation_idx_1 + a_1[k] * Saturation_idx_0)
                    + a_1[k + 252] * Saturation_idx_2) + (static_cast<real_T>
                    (b_Mlim) + Bc_2));
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

          rtb_Sum2_f[k] = ((b_Ku1_1[3 * k + 1] * Saturation_idx_1 + b_Ku1_1[3 *
                            k] * Saturation_idx_0) + b_Ku1_1[3 * k + 2] *
                           Saturation_idx_2) + (dwt + Bc_2);
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

        // Saturate: '<S5>/Saturation' incorporates:
        //   Gain: '<S230>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_g[0] * umax_incr[0];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[0] = Saturation_idx_0 - rtU.u0[0];

        // Sum: '<S229>/Sum6'
        rtb_Sum6[0] = y__m[0];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Saturate: '<S5>/Saturation' incorporates:
        //   Gain: '<S230>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_g[1] * umax_incr[1];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[1] = Saturation_idx_0 - rtU.u0[1];

        // Sum: '<S229>/Sum6'
        rtb_Sum6[1] = y__m[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Saturate: '<S5>/Saturation' incorporates:
        //   Gain: '<S230>/umin_scale1'

        Saturation_idx_0 = rtP.umin_scale1_Gain_g[2] * umax_incr[2];
        if (Saturation_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S229>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Sum1[2] = Saturation_idx_0 - rtU.u0[2];

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

        rtY.u[2] = Saturation_idx_0;

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

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/wtMod'

    for (i = 0; i < 6; i++) {
      // InitializeConditions for Delay: '<S8>/Delay'
      rtDW.Delay_DSTATE[i] = rtP.ywt0[i];
    }
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
