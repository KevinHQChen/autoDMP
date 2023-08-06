//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2406
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sun Aug  6 01:50:56 2023
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

// Named constants for MATLAB Function: '<S38>/FixedHorizonOptimizer'
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

  // MATLAB Function: '<S157>/MATLAB Function1'
  in3_idx_0 = in3 - in2;
  stride_0_0 = (in7 - in6) + 1 != 1 ? static_cast<int32_T>(1) : static_cast<
    int32_T>(0);
  for (int32_T i{0}; i < in3_idx_0; i++) {
    in1[(in2 + i) + 12 * in4] = in5[i * stride_0_0 + in6] * in8[(in4 << 2UL) + i];
  }

  // End of MATLAB Function: '<S157>/MATLAB Function1'
}

//
// System initialize for function-call system:
//    '<S1>/paramEst1'
//    '<S1>/paramEst2'
//
void SupervisoryController::paramEst1_Init(real_T rty_theta[12], real_T rty_P
  [144], real_T rty_err[3], DW_paramEst1 *localDW, P_paramEst1 *localP)
{
  // InitializeConditions for Delay: '<S159>/Delay1'
  localDW->icLoad = true;

  // InitializeConditions for UnitDelay: '<S157>/Unit Delay3'
  localDW->UnitDelay3_DSTATE[0] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[1] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[2] = localP->UnitDelay3_InitialCondition;

  // InitializeConditions for Delay: '<S159>/Delay'
  localDW->icLoad_n = true;

  // SystemInitialize for Outport: '<S4>/theta'
  for (int32_T i{0}; i < 12; i++) {
    rty_theta[i] = localP->theta_Y0;
  }

  // End of SystemInitialize for Outport: '<S4>/theta'

  // SystemInitialize for Outport: '<S4>/P'
  for (int32_T i{0}; i < 144; i++) {
    rty_P[i] = localP->P_Y0;
  }

  // End of SystemInitialize for Outport: '<S4>/P'

  // SystemInitialize for Outport: '<S4>/err'
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
  real_T rtb_Add1_c[3];
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

  // Delay: '<S159>/Delay1'
  localDW->icLoad = ((rtu_rstP && (static_cast<uint32_T>
    (localZCE->Delay1_Reset_ZCE) != POS_ZCSIG)) || localDW->icLoad);
  localZCE->Delay1_Reset_ZCE = rtu_rstP ? static_cast<ZCSigState>(1) :
    static_cast<ZCSigState>(0);
  if (localDW->icLoad) {
    (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rtu_P0[0], 144U * sizeof
                      (real_T));
  }

  // Sum: '<S157>/Add1'
  rtb_Add1_c[0] = rtu_y[0] - rtu_y0[0];

  // Sum: '<S157>/Add3'
  rtb_Add3_idx_0 = rtu_u[0] - rtu_u0[0];

  // Sum: '<S157>/Add1'
  rtb_Add1_c[1] = rtu_y[1] - rtu_y0[1];

  // Sum: '<S157>/Add3'
  rtb_Add3_idx_1 = rtu_u[1] - rtu_u0[1];

  // Sum: '<S157>/Add1'
  rtb_Add1_c[2] = rtu_y[2] - rtu_y0[2];

  // Sum: '<S157>/Add3'
  rtb_Add3_idx_2 = rtu_u[2] - rtu_u0[2];

  // MATLAB Function: '<S157>/MATLAB Function1' incorporates:
  //   Sum: '<S157>/Add1'
  //   Sum: '<S157>/Add3'

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1': '<S158>:1' 
  // '<S158>:1:2' [z, phi] = getRegressors_(y, yPrev, u, sign_, no, ni, np, dt, mdlNum); 
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
    regs[i] = rtb_Add1_c[p1];
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

  // Delay: '<S159>/Delay'
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
    // Math: '<S159>/Transpose' incorporates:
    //   MATLAB Function: '<S159>/MATLAB Function'

    tmp[3 * i] = rtb_phi[i];
    tmp[3 * i + 1] = rtb_phi[i + 12];
    tmp[3 * i + 2] = rtb_phi[i + 24];
  }

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function': '<S160>:1' 
  // '<S160>:1:2' [dtheta, dP, L] = rls_(theta, phi, epsil, EN, p_, dPmod_, lambda, P, no, ni, np); 
  // 'rls_:3' dtheta = zeros(no*np, 1);
  // 'rls_:4' dP = zeros(no*np,no*np);
  // 'rls_:5' L = zeros(no*np, no);
  // 'rls_:7' L = P*phi*inv(lambda*eye(no) + phi'*P*phi);
  for (i = 0; i < 3; i++) {
    // Sum: '<S159>/Sum2' incorporates:
    //   Delay: '<S159>/Delay'
    //   MATLAB Function: '<S157>/MATLAB Function1'
    //   Math: '<S159>/Transpose'
    //   Sum: '<S157>/Add1'
    //   UnitDelay: '<S157>/Unit Delay3'

    rtb_Add3_idx_0 = 0.0;
    for (p2 = 0; p2 < 12; p2++) {
      // MATLAB Function: '<S159>/MATLAB Function'
      p1 = 3 * p2 + i;
      rtb_Add3_idx_0 += tmp[p1] * localDW->Delay_DSTATE[p2];

      // MATLAB Function: '<S159>/MATLAB Function' incorporates:
      //   Delay: '<S159>/Delay'
      //   Delay: '<S159>/Delay1'

      tmp_0[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        tmp_0[p1] += tmp[3 * ibtile + i] * localDW->Delay1_DSTATE[12 * p2 +
          ibtile];
      }
    }

    rty_err[i] = (rtb_Add1_c[i] - localDW->UnitDelay3_DSTATE[i]) -
      rtb_Add3_idx_0;

    // End of Sum: '<S159>/Sum2'

    // MATLAB Function: '<S159>/MATLAB Function'
    for (p2 = 0; p2 < 3; p2++) {
      rtb_Add3_idx_0 = 0.0;
      for (p1 = 0; p1 < 12; p1++) {
        rtb_Add3_idx_0 += tmp_0[3 * p1 + i] * rtb_phi[12 * p2 + p1];
      }

      p1 = 3 * p2 + i;
      b_b[p1] = static_cast<real_T>(b_b_0[p1]) * rtu_lambda + rtb_Add3_idx_0;
    }
  }

  // MATLAB Function: '<S159>/MATLAB Function' incorporates:
  //   Delay: '<S159>/Delay'
  //   Delay: '<S159>/Delay1'

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
    // Delay: '<S159>/Delay1'
    tmp_1 = _mm_loadu_pd(&localDW->Delay1_DSTATE[i]);

    // Product: '<S159>/Product1' incorporates:
    //   Delay: '<S159>/Delay1'

    tmp_2 = _mm_loadu_pd(&rty_P[i]);
    (void)_mm_storeu_pd(&rty_P[i], _mm_div_pd(_mm_sub_pd(tmp_1, tmp_2),
      _mm_set1_pd(rtu_lambda)));
  }

  for (i = 0; i <= 10; i += 2) {
    // Delay: '<S159>/Delay'
    tmp_1 = _mm_loadu_pd(&localDW->Delay_DSTATE[i]);

    // Sum: '<S159>/Sum' incorporates:
    //   Delay: '<S159>/Delay'

    tmp_2 = _mm_loadu_pd(&rty_theta[i]);
    (void)_mm_storeu_pd(&rty_theta[i], _mm_add_pd(tmp_1, tmp_2));
  }

  // Update for Delay: '<S159>/Delay1'
  localDW->icLoad = false;
  (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rty_P[0], 144U * sizeof(real_T));

  // Update for UnitDelay: '<S157>/Unit Delay3' incorporates:
  //   Sum: '<S157>/Add1'

  localDW->UnitDelay3_DSTATE[0] = rtb_Add1_c[0];
  localDW->UnitDelay3_DSTATE[1] = rtb_Add1_c[1];
  localDW->UnitDelay3_DSTATE[2] = rtb_Add1_c[2];

  // Update for Delay: '<S159>/Delay'
  localDW->icLoad_n = false;
  (void)std::memcpy(&localDW->Delay_DSTATE[0], &rty_theta[0], 12U * sizeof
                    (real_T));
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

// Function for MATLAB Function: '<S3>/MATLAB Function'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S109>/optimizer'
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

// Function for MATLAB Function: '<S113>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv(const real_T A[12], const real_T B_0[9],
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

  static const real_T b_Kr[240]{ -0.00384284772673063, -0.0, -0.0, -0.0,
    -0.0076717186107965888, -0.0, -0.0, -0.0, -0.011486663487449431, -0.0, -0.0,
    -0.0, -0.015287733007047537, -0.0, -0.0, -0.0, -0.019074977635728591, -0.0,
    -0.0, -0.0, -0.022848447656079605, -0.0, -0.0, -0.0, -0.026608193167804516,
    -0.0, -0.0, -0.0, -0.030354264088389361, -0.0, -0.0, -0.0,
    -0.034086710153765, -0.0, -0.0, -0.0, -0.037805580918967477, -0.0, -0.0,
    -0.0, -0.04151092575879594, -0.0, -0.0, -0.0, -0.04520279386846817, -0.0,
    -0.0, -0.0, -0.048881234264273758, -0.0, -0.0, -0.0, -0.052546295784224879,
    -0.0, -0.0, -0.0, -0.0561980270887047, -0.0, -0.0, -0.0,
    -0.059836476661113455, -0.0, -0.0, -0.0, -0.063461692808512118, -0.0, -0.0,
    -0.0, -0.067073723662263787, -0.0, -0.0, -0.0, -0.070672617178672723, -0.0,
    -0.0, -0.0, -0.074258421139621028, -0.0, -0.0, -0.0, 0.0041395750008033473,
    -0.0, -0.0, -0.0, 0.0082640939305366964, -0.0, -0.0, 0.0,
    0.01237361154971919, -0.0, -0.0, 0.0, 0.016468182419700179, -0.0, -0.0, 0.0,
    0.020547860903383636, -0.0, -0.0, 0.0, 0.024612701165949907, -0.0, -0.0, 0.0,
    0.028662757175574875, -0.0, -0.0, 0.0, 0.03269808270414646, -0.0, -0.0, 0.0,
    0.036718731327978558, -0.0, -0.0, 0.0, 0.040724756428522364, -0.0, -0.0, 0.0,
    0.044716211193075089, -0.0, -0.0, 0.0, 0.048693148615486133, -0.0, -0.0, 0.0,
    0.052655621496860691, -0.0, -0.0, 0.0, 0.056603682446260746, -0.0, -0.0, 0.0,
    0.060537383881403589, -0.0, -0.0, 0.0, 0.064456778029357731, -0.0, -0.0, 0.0,
    0.068361916927236324, -0.0, -0.0, 0.0, 0.072252852422888064, -0.0, -0.0, 0.0,
    0.076129636175585519, -0.0, -0.0, 0.0, 0.079992319656711067, -0.0, -0.0, 0.0,
    0.00056155944692921835, -0.0, -0.0, -0.0, 0.0011210764428963565, -0.0, -0.0,
    0.0, 0.0016785584165110745, -0.0, -0.0, 0.0, 0.0022340127693643965, -0.0,
    -0.0, 0.0, 0.0027874468761269785, -0.0, -0.0, 0.0, 0.0033388680846470208,
    -0.0, -0.0, 0.0, 0.0038882837160478245, -0.0, -0.0, 0.0,
    0.0044357010648249937, -0.0, -0.0, 0.0, 0.0049811273989432804, -0.0, -0.0,
    0.0, 0.0055245699599330853, -0.0, -0.0, 0.0, 0.0060660359629865959, -0.0,
    -0.0, 0.0, 0.0066055325970535853, -0.0, -0.0, 0.0, 0.0071430670249368576,
    -0.0, -0.0, 0.0, 0.007678646383387347, -0.0, -0.0, 0.0,
    0.0082122777831988725, -0.0, -0.0, 0.0, 0.0087439683093025449, -0.0, -0.0,
    0.0, 0.0092737250208608339, -0.0, -0.0, 0.0, 0.0098015549513612889, -0.0,
    -0.0, 0.0, 0.010327465108709924, -0.0, -0.0, 0.0, 0.010851462475324258, -0.0,
    -0.0, 0.0 };

  static const real_T b_Mlim[166]{ 612.0, 612.0, 612.0, 15.3, 612.0, 612.0,
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

  static const real_T b_Hinv[16]{ 15.969185026172692, 14.367770534785548,
    1.9490786550682766, 0.0, 14.367770534785548, 13.829880342623046,
    -2.0995776697049715, 0.0, 1.9490786550682766, -2.0995776697049715,
    29.022244012256408, 0.0, 0.0, 0.0, 0.0, 9.9999999999999974E-6 };

  static const real_T b_Linv[16]{ 0.66952848091349959, 3.9230116272599145,
    0.36179609146848607, 0.0, 0.0, 3.6983765413978609, -0.38973234489972963, 0.0,
    0.0, 0.0, 5.3872297159353071, 0.0, 0.0, 0.0, 0.0, 0.003162277660168379 };

  static const real_T b_Kx[15]{ 0.75062707300997611, 0.0, 0.78870934466740916,
    0.0, 0.0, -0.8085870914801655, 0.0, -0.84960990345118026, 0.0, 0.0,
    -0.10968993671997439, 0.0, -0.11525493979334633, 0.0, 0.0 };

  static const real_T b_Ku1[9]{ 2.1966849657348271, -2.3663030167819756,
    -0.32100392265224076, -2.3663030167819761, 2.5490182045099905,
    0.34579038980072829, -0.32100392265224081, 0.34579038980072824,
    0.046908646421976154 };

  static const uint8_T b_Mrows_0[166]{ 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U, 10U,
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

  real_T Bc[166];
  real_T a__1[166];
  real_T Product1_j[144];
  real_T rseq[80];
  real_T B_est[36];
  real_T D_est[27];
  real_T rtb_R_tmp[27];
  real_T y[24];
  real_T Abar[16];
  real_T rtb_A[16];
  real_T rtb_Q[16];
  real_T rtb_Transpose2_0[16];
  real_T rtb_Z[16];
  real_T rtb_y_g[16];
  real_T rtb_y_m[16];
  real_T Sum_h[12];
  real_T rtb_Add_e[12];
  real_T rtb_C[12];
  real_T rtb_N[12];
  real_T rtb_Product2[12];
  real_T rtb_Transpose2[12];
  real_T rtb_Add1[9];
  real_T rtb_N_0[9];
  real_T rtb_R[9];
  real_T tmp[9];
  real_T rtb_ywt[6];
  real_T rtb_ywtT[6];
  real_T rtb_xest[5];
  real_T rtb_Sum2[4];
  real_T rtb_TmpSignalConversionAtSFun_p[4];
  real_T Sum2_c[3];
  real_T rtb_Add1_n[3];
  real_T rtb_C_0[3];
  real_T rtb_y_[3];
  real_T tmp_0[3];
  real_T umax_incr[3];
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
        for (int32_T ii{0}; ii < 6; ii++) {
          // Inport: '<Root>/y'
          rtDW.traj[ii] = rtU.y[ii];
        }
      } else {
        // '<S1>:520:6' else
        //  hold last waypoint pos
        // '<S1>:520:7' traj(:,1) = traj(:, waypt);
        for (int32_T ii{0}; ii < 6; ii++) {
          rtDW.traj[ii] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6 +
            ii];
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
      __m128d tmp_1;
      __m128d tmp_2;
      __m128d tmp_3;
      real_T Bc_0;
      real_T dwt;
      real_T umin_scale1_idx_0;
      real_T umin_scale1_idx_1;
      real_T umin_scale1_idx_2;
      int32_T i;
      int32_T i_0;
      int32_T i_2;
      int32_T ii;
      int32_T rtb_Q_tmp;
      int32_T rtb_R_tmp_0;
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
            for (ii = 0; ii < 6; ii++) {
              // Inport: '<Root>/y'
              rtDW.traj[ii] = rtU.y[ii];
            }
          } else {
            // '<S1>:520:6' else
            //  hold last waypoint pos
            // '<S1>:520:7' traj(:,1) = traj(:, waypt);
            for (ii = 0; ii < 6; ii++) {
              rtDW.traj[ii] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                6 + ii];
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
          for (ii = 0; ii < 12; ii++) {
            P0_2_tmp[ii] = static_cast<int8_T>(ii + 1);
          }

          // '<S1>:59:5' theta0_2 = theta(1:np*no);
          for (k = 0; k < 12; k++) {
            for (ii = 0; ii < 12; ii++) {
              rtDW.P0_2[ii + 12 * k] = rtY.P_g[((P0_2_tmp[k] - 1) * 24 +
                P0_2_tmp[ii]) - 1];
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
          for (ii = 0; ii < 12; ii++) {
            P0_2_tmp[ii] = static_cast<int8_T>(ii + 13);
          }

          // '<S1>:59:11' theta0_1 = theta(np*no+1:2*np*no);
          for (k = 0; k < 12; k++) {
            for (ii = 0; ii < 12; ii++) {
              rtDW.P0_1[ii + 12 * k] = rtY.P_g[((P0_2_tmp[k] - 1) * 24 +
                P0_2_tmp[ii]) - 1];
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
      i_2 = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S159>/Product1'

        (void)std::memcpy(&rtY.P_g[k], &Product1_j[i_2], 12U * sizeof(real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i] = Sum_h[i];
        k += 24;
        i_2 += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S159>/Sum2'

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
      i_2 = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S163>/Product1'

        (void)std::memcpy(&rtY.P_g[k + 300], &Product1_j[i_2], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i + 12] = Sum_h[i];
        k += 24;
        i_2 += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S163>/Sum2'

      rtY.prmErr[3] = Sum2_c[0];
      rtY.prmErr[4] = Sum2_c[1];
      rtY.prmErr[5] = Sum2_c[2];

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      // MATLAB Function: '<S3>/MATLAB Function' incorporates:
      //   Inport: '<Root>/k_2'

      // [u, ywt, yhat, currTraj] = ...
      // ampc(traj(:,waypt), currEv.r, y, ymax, y0, x0, u0, umax, uwt,...
      // excitation, theta, thetaSgn, k_2);
      // [u, ywt, currTraj] = gmpc(traj(:,waypt), currEv.r, y, ymax, umax, uwt, k_2); 
      // '<S1>:59:32' [u, yhat(1:no), currTraj, ywt] = mpc1(traj(:,waypt), currEv.r, y, [0;0;0], 0, u0, umax, k_2, false); 
      // Simulink Function 'mpc1': '<S1>:882'
      // MATLAB Function 'SupervisoryController/mpc1/MATLAB Function': '<S86>:1' 
      // '<S86>:1:2' [ywt, ywtT, uwt, uwtT] = wtMod_(y, yDest, ywtT, uwtT, dt, no, ni, k_2); 
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
        // Delay: '<S3>/Delay'
        umin_scale1_idx_2 = rtDW.Delay_DSTATE[k];

        // MATLAB Function: '<S3>/MATLAB Function' incorporates:
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

        // Delay: '<S3>/Delay'
        rtb_ywtT[k] = umin_scale1_idx_2;
      }

      // MATLAB Function: '<S3>/MATLAB Function' incorporates:
      //   Inport: '<Root>/y'
      //   Outport: '<Root>/currEv'

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
      // '<S86>:1:3' r_ = zeros(no, 1);
      // '<S86>:1:4' if any(yDest(1:no))
      if (any(&rtY.currEv.r[0])) {
        // '<S86>:1:5' r_ = r(1:no);
        // '<S86>:1:6' y_ = y(1:no);
        ii = (static_cast<int32_T>(rtDW.waypt) - 1) * 6;
        Sum2_c[0] = rtDW.traj[ii];
        rtb_y_[0] = rtU.y[0];
        Sum2_c[1] = rtDW.traj[ii + 1];
        rtb_y_[1] = rtU.y[1];
        Sum2_c[2] = rtDW.traj[ii + 2];
        rtb_y_[2] = rtU.y[2];
      } else if (any(&rtY.currEv.r[3])) {
        // '<S86>:1:7' elseif any(yDest(no+1:2*no))
        // '<S86>:1:8' r_ = r(no+1:2*no);
        // '<S86>:1:9' y_ = y(no+1:2*no);
        ii = (static_cast<int32_T>(rtDW.waypt) - 1) * 6;
        Sum2_c[0] = rtDW.traj[ii + 3];
        rtb_y_[0] = rtU.y[3];
        Sum2_c[1] = rtDW.traj[ii + 4];
        rtb_y_[1] = rtU.y[4];
        Sum2_c[2] = rtDW.traj[ii + 5];
        rtb_y_[2] = rtU.y[5];
      } else {
        // '<S86>:1:10' else
        // '<S86>:1:11' r_ = zeros(no, 1);
        // '<S86>:1:12' y_ = y(1:no);
        Sum2_c[0] = 0.0;
        rtb_y_[0] = rtU.y[0];
        Sum2_c[1] = 0.0;
        rtb_y_[1] = rtU.y[1];
        Sum2_c[2] = 0.0;
        rtb_y_[2] = rtU.y[2];
      }

      // Delay: '<S111>/MemoryX' incorporates:
      //   Constant: '<S111>/X0'
      //   DataTypeConversion: '<S111>/DataTypeConversionReset'

      rtDW.icLoad = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE) ==
                      POS_ZCSIG) || rtDW.icLoad);
      rtPrevZCX.MemoryX_Reset_ZCE = 0U;
      if (rtDW.icLoad) {
        rtDW.MemoryX_DSTATE[0] = rtP.X0_Value_f[0];
        rtDW.MemoryX_DSTATE[1] = rtP.X0_Value_f[1];
        rtDW.MemoryX_DSTATE[2] = rtP.X0_Value_f[2];
        rtDW.MemoryX_DSTATE[3] = rtP.X0_Value_f[3];
      }

      // SignalConversion generated from: '<S110>/ SFunction ' incorporates:
      //   MATLAB Function: '<S109>/optimizer'

      rtb_TmpSignalConversionAtSFun_p[0] = Sum2_c[0];
      rtb_TmpSignalConversionAtSFun_p[1] = Sum2_c[1];
      rtb_TmpSignalConversionAtSFun_p[2] = Sum2_c[2];

      // MATLAB Function: '<S109>/optimizer' incorporates:
      //   Constant: '<S3>/Constant'
      //   Constant: '<S88>/Constant1'
      //   Delay: '<S111>/MemoryX'
      //   DiscreteIntegrator: '<S3>/Discrete-Time Integrator'
      //   Gain: '<S3>/Gain2'
      //   Inport: '<Root>/umax'
      //   SignalConversion generated from: '<S110>/ SFunction '
      //   Sum: '<S88>/Sum2'
      //   UnitDelay: '<S89>/last_mv'

      // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S110>:1'
      // '<S110>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
      // '<S110>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
      // '<S110>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
      // '<S110>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
      //  Inputs (in BlockDataType except iA)
      //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
      // '<S110>:1:24' xk = convertDataType(xk0,isDouble);
      // '<S110>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S110>:1:252' if isa(u,'double')
      // '<S110>:1:253' y = u;
      //    old_u:      last mv (calculated by MPC)
      // '<S110>:1:26' old_u = convertDataType(old_u0,isDouble);
      // '<S110>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S110>:1:252' if isa(u,'double')
      // '<S110>:1:253' y = u;
      //    ym:         current measured output (used only with built-in KF)
      // '<S110>:1:28' ym = convertDataType(ym0,isDouble);
      //    ref:        output reference
      // '<S110>:1:30' ref = convertDataType(ref0,isDouble);
      // '<S110>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S110>:1:252' if isa(u,'double')
      // '<S110>:1:253' y = u;
      //    md:         measured disturbance
      // '<S110>:1:32' md = convertDataType(md0,isDouble);
      //    umin:       run-time MV bound
      // '<S110>:1:34' umin = convertDataType(umin0,isDouble);
      //    umax:       run-time MV bound
      // '<S110>:1:36' umax = convertDataType(umax0,isDouble);
      // '<S110>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S110>:1:252' if isa(u,'double')
      // '<S110>:1:253' y = u;
      //    ymin:       run-time OV bound
      // '<S110>:1:38' ymin = convertDataType(ymin0,isDouble);
      //    ymax:       run-time OV bound
      // '<S110>:1:40' ymax = convertDataType(ymax0,isDouble);
      //    E:          run-time mixed constraints
      // '<S110>:1:42' E = convertDataType(E0,isDouble);
      //    F:          run-time mixed constraints
      // '<S110>:1:44' F = convertDataType(F0,isDouble);
      //    G:          run-time mixed constraints
      // '<S110>:1:46' G = convertDataType(G0,isDouble);
      //    S:          run-time mixed constraints
      // '<S110>:1:48' S = convertDataType(S0,isDouble);
      //    switch_in:  if it matches "enable_value", MPC is active in control
      // '<S110>:1:50' switch_in = int32(switch_in0);
      //    ext_mv:     external last mv (actual)
      // '<S110>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
      //    MVtarget:   MV reference
      // '<S110>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
      //    ywt:        run-time OV weights
      // '<S110>:1:56' ywt = convertDataType(ywt0,isDouble);
      //    uwt:        run-time MV weights
      // '<S110>:1:58' uwt = convertDataType(uwt0,isDouble);
      //    duwt:       run-time DMV weights
      // '<S110>:1:60' duwt = convertDataType(duwt0,isDouble);
      //    ewt:     run-time Slack weights
      // '<S110>:1:62' ewt = convertDataType(ewt0,isDouble);
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
      // '<S110>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
      // '<S110>:1:96' isAdaptive = false;
      // '<S110>:1:97' isLTV = false;
      // '<S110>:1:98' ZERO = zeros('like',ref);
      // '<S110>:1:99' ONE = ones('like',ref);
      // '<S110>:1:100' hasMD = nv>int32(1);
      //  Pre-allocate all the MEX block outputs for the simulation mode
      // '<S110>:1:105' if isSimulation
      //  Get reference and MD signals -- accounting for previewing
      // '<S110>:1:119' if isSimulation
      // '<S110>:1:126' else
      //  When doing code generation, use M code directly
      // '<S110>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
      ii = 0;
      for (k = 0; k < 20; k++) {
        rseq[ii] = rtb_TmpSignalConversionAtSFun_p[0];
        rseq[ii + 1] = rtb_TmpSignalConversionAtSFun_p[1];
        rseq[ii + 2] = rtb_TmpSignalConversionAtSFun_p[2];
        rseq[ii + 3] = rtP.Constant_Value;
        ii += 4;
      }

      //  External MV override.
      //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
      // '<S110>:1:133' old_u = old_u - uoff;
      // '<S110>:1:134' if no_mv
      // '<S110>:1:135' delmv = zeros(nu,1,'like',ref);
      //  Obtain x[k|k]
      // '<S110>:1:143' xk = xk - xoff;
      rtb_xest[0] = rtDW.MemoryX_DSTATE[0];
      rtb_xest[1] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE;
      umin_scale1_idx_0 = rtDW.last_mv_DSTATE[0];
      rtb_xest[2] = rtP.Constant1_Value_j[0] + rtDW.MemoryX_DSTATE[1];
      umin_scale1_idx_1 = rtDW.last_mv_DSTATE[1];
      rtb_xest[3] = rtP.Constant1_Value_j[1] + rtDW.MemoryX_DSTATE[2];
      umin_scale1_idx_2 = rtDW.last_mv_DSTATE[2];
      rtb_xest[4] = rtP.Constant1_Value_j[2] + rtDW.MemoryX_DSTATE[3];

      //  Remove offset
      // '<S110>:1:144' if CustomEstimation
      //  Input state is x(k|k)
      // '<S110>:1:146' xest = xk;
      //  Real-time MV target override
      //  Note: utargetValue is a vector length p*nu.
      // '<S110>:1:162' if no_uref
      //  no external utarget
      // '<S110>:1:164' utargetValue = utarget;
      //  Real-time custom constraint override (scaled E/F/S)
      // '<S110>:1:173' if ~no_cc
      // '<S110>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
      // '<S110>:1:183' if isSimulation
      // '<S110>:1:214' else
      //  When doing code generation, use M code directly
      // '<S110>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
      // '<S110>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
      // '<S110>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
      // '<S110>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
      // '<S110>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
      // '<S110>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
      // '<S110>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
      // '<S110>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
      // '<S110>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
      // '<S110>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
      // '<S110>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
      // '<S110>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
      // '<S110>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
      umax_incr_flag[0] = false;
      umax_incr[0] = 0.0;
      umax_incr_flag[1] = false;
      umax_incr[1] = 0.0;
      umax_incr_flag[2] = false;
      umax_incr[2] = 0.0;
      for (k = 0; k < 166; k++) {
        uint8_T b_Mrows;
        dwt = b_Mlim[k];
        Bc_0 = 0.0;
        for (ii = 0; ii < 5; ii++) {
          Bc_0 += b_a[166 * ii + k] * rtb_xest[ii];
        }

        Bc_0 = -(((a[k + 166] * umin_scale1_idx_1 + a[k] * umin_scale1_idx_0) +
                  a[k + 332] * umin_scale1_idx_2) + (dwt + Bc_0));
        b_Mrows = b_Mrows_0[k];
        if ((b_Mrows > 80UL) && (b_Mrows > 160UL) && (b_Mrows <= 220UL)) {
          ii = (static_cast<int32_T>(b_Mrows) - div_nde_s32_floor
                (static_cast<int32_T>(b_Mrows) - 161, static_cast<int32_T>(nu)) *
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
          Bc_0 += dwt;
        }

        Bc[k] = Bc_0;
      }

      rtb_Sum2[0] = 0.0;
      rtb_Sum2[1] = 0.0;
      rtb_Sum2[2] = 0.0;
      rtb_Sum2[3] = 0.0;
      for (k = 0; k < 3; k++) {
        dwt = 0.0;
        for (ii = 0; ii < 5; ii++) {
          dwt += b_Kx[5 * k + ii] * rtb_xest[ii];
        }

        Bc_0 = 0.0;
        for (ii = 0; ii < 80; ii++) {
          Bc_0 += b_Kr[80 * k + ii] * rseq[ii];
        }

        rtb_Sum2[k] = ((b_Ku1[3 * k + 1] * umin_scale1_idx_1 + b_Ku1[3 * k] *
                        umin_scale1_idx_0) + b_Ku1[3 * k + 2] *
                       umin_scale1_idx_2) + (dwt + Bc_0);
      }

      // Update for Memory: '<S89>/Memory' incorporates:
      //   MATLAB Function: '<S109>/optimizer'

      qpkwik(b_Linv, b_Hinv, rtb_Sum2, b_Ac, Bc, rtDW.Memory_PreviousInput, 680,
             1.0E-6, rtb_TmpSignalConversionAtSFun_p, a__1, &k);

      // MATLAB Function: '<S109>/optimizer' incorporates:
      //   UnitDelay: '<S89>/last_mv'

      if ((k < 0) || (k == 0)) {
        rtb_TmpSignalConversionAtSFun_p[0] = 0.0;
        rtb_TmpSignalConversionAtSFun_p[1] = 0.0;
        rtb_TmpSignalConversionAtSFun_p[2] = 0.0;
      }

      umax_incr[0] = rtDW.last_mv_DSTATE[0] + rtb_TmpSignalConversionAtSFun_p[0];
      umax_incr[1] = rtDW.last_mv_DSTATE[1] + rtb_TmpSignalConversionAtSFun_p[1];
      umax_incr[2] = rtDW.last_mv_DSTATE[2] + rtb_TmpSignalConversionAtSFun_p[2];

      // Delay: '<S111>/MemoryP' incorporates:
      //   Constant: '<S111>/P0'
      //   DataTypeConversion: '<S111>/DataTypeConversionReset'

      // '<S110>:1:231' if return_xseq || return_ovseq
      // '<S110>:1:233' else
      // '<S110>:1:234' yseq = zeros(p+1,ny,'like',rseq);
      // '<S110>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
      // '<S110>:1:238' if CustomEstimation
      // '<S110>:1:239' xk1 = xk;
      // '<S110>:1:244' xk1 = xk1 + xoff;
      //  Updated state must include offset
      //  return xest in original value
      // '<S110>:1:247' xest = xest + xoff;
      rtDW.icLoad_h = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE) ==
                        POS_ZCSIG) || rtDW.icLoad_h);
      rtPrevZCX.MemoryP_Reset_ZCE = 0U;
      if (rtDW.icLoad_h) {
        (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtP.P0_Value_a[0], sizeof
                          (real_T) << 4UL);
      }

      // MATLAB Function: '<S88>/MATLAB Function' incorporates:
      //   BusCreator: '<S3>/Bus Creator1'
      //   Constant: '<S3>/Constant12'
      //   Constant: '<S3>/Constant13'
      //   Constant: '<S3>/Constant3'
      //   Constant: '<S3>/Constant4'
      //   MATLAB Function: '<S133>/ScalarExpansion'

      // MATLAB Function 'SupervisoryController/mpc1/State Estimator OD (KF)/MATLAB Function': '<S112>:1' 
      // '<S112>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod1, Bod1, Cod1, Dod1, Dmn1, 3, 3, 1); 
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
      (void)std::memset(&rtb_A[0], 0, sizeof(real_T) << 4UL);
      rtb_A[0] = rtP.Constant3_Value;

      // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
      // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
      k = 0;
      i_2 = 0;
      for (i = 0; i < 3; i++) {
        rtb_A[k + 5] = rtP.Aod1[i_2];
        rtb_A[k + 6] = rtP.Aod1[i_2 + 1];
        rtb_A[k + 7] = rtP.Aod1[i_2 + 2];
        Sum_h[k] = rtP.Constant4_Value[i];
        rtb_C[i] = rtP.Constant12_Value_e[i];
        k += 4;
        i_2 += 3;
      }

      (void)std::memcpy(&rtb_C[3], &rtP.Cod1[0], 9U * sizeof(real_T));

      // 'stateEst_:22' D(1:no, 1:ni) = Dp;
      // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
      (void)std::memset(&B_est[0], 0, 36U * sizeof(real_T));

      // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
      (void)std::memset(&y[0], 0, 24U * sizeof(real_T));
      ii = 0;
      k = 0;
      for (i_2 = 0; i_2 < 3; i_2++) {
        y[ii] = rtP.Constant4_Value[i_2];
        y[ii + 13] = rtP.Bod1[k];
        y[ii + 14] = rtP.Bod1[k + 1];
        y[ii + 15] = rtP.Bod1[k + 2];
        ii += 4;
        k += 3;
      }

      ii = 0;
      for (k = 0; k < 6; k++) {
        B_est[ii] = y[ii];
        B_est[ii + 1] = y[ii + 1];
        B_est[ii + 2] = y[ii + 2];
        B_est[ii + 3] = y[ii + 3];
        ii += 4;
      }

      // 'stateEst_:26' D_est = [Dp Dod Dn];
      for (ii = 0; ii < 9; ii++) {
        D_est[ii] = rtP.Constant13_Value_c[ii];
        D_est[ii + 9] = rtP.Dod1[ii];
        D_est[ii + 18] = rtP.Dmn1[ii];
      }

      // 'stateEst_:27' Q = B_est * B_est';
      for (ii = 0; ii < 4; ii++) {
        k = 0;
        for (i_2 = 0; i_2 < 4; i_2++) {
          rtb_Q_tmp = k + ii;
          rtb_Q[rtb_Q_tmp] = 0.0;
          i = 0;
          for (i_0 = 0; i_0 < 9; i_0++) {
            rtb_Q[rtb_Q_tmp] += B_est[i + ii] * B_est[i + i_2];
            i += 4;
          }

          k += 4;
        }
      }

      // 'stateEst_:28' R = D_est * D_est';
      ii = 0;
      for (k = 0; k < 9; k++) {
        rtb_R_tmp[k] = D_est[ii];
        rtb_R_tmp[k + 9] = D_est[ii + 1];
        rtb_R_tmp[k + 18] = D_est[ii + 2];
        ii += 3;
      }

      ii = 0;
      k = 0;
      for (i_2 = 0; i_2 < 3; i_2++) {
        for (i = 0; i < 3; i++) {
          rtb_R_tmp_0 = i + ii;
          rtb_R[rtb_R_tmp_0] = 0.0;
          i_0 = 0;
          for (rtb_Q_tmp = 0; rtb_Q_tmp < 9; rtb_Q_tmp++) {
            rtb_R[rtb_R_tmp_0] += D_est[i_0 + i] * rtb_R_tmp[rtb_Q_tmp + k];
            i_0 += 3;
          }
        }

        ii += 3;
        k += 9;
      }

      // Outputs for Atomic SubSystem: '<S111>/ScalarExpansionQ'
      // 'stateEst_:29' N = B_est * D_est';
      //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
      //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
      //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
      //  yhat = C*xhat + D*u;
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
      // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S155>:1' 
      //    Copyright 2014-2015 The MathWorks, Inc.
      // '<S155>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
      ii = 0;
      for (k = 0; k < 4; k++) {
        i_2 = 0;
        i = 0;
        for (i_0 = 0; i_0 < 3; i_0++) {
          rtb_R_tmp_0 = i_2 + k;
          rtb_N[rtb_R_tmp_0] = 0.0;
          rtb_Q_tmp = 0;
          for (int32_T i_1{0}; i_1 < 9; i_1++) {
            rtb_N[rtb_R_tmp_0] += B_est[rtb_Q_tmp + k] * rtb_R_tmp[i_1 + i];
            rtb_Q_tmp += 4;
          }

          i_2 += 4;
          i += 9;
        }

        rtb_y_g[ii] = (rtb_Q[ii] + rtb_Q[k]) / 2.0;
        rtb_y_g[ii + 1] = (rtb_Q[ii + 1] + rtb_Q[k + 4]) / 2.0;
        rtb_y_g[ii + 2] = (rtb_Q[ii + 2] + rtb_Q[k + 8]) / 2.0;
        rtb_y_g[ii + 3] = (rtb_Q[ii + 3] + rtb_Q[k + 12]) / 2.0;
        ii += 4;
      }

      // End of MATLAB Function: '<S88>/MATLAB Function'
      // End of Outputs for SubSystem: '<S111>/ScalarExpansionQ'

      // Outputs for Atomic SubSystem: '<S111>/ReducedQRN'
      for (ii = 0; ii < 4; ii++) {
        // Product: '<S131>/Product' incorporates:
        //   Constant: '<S111>/G'
        //   Math: '<S131>/Transpose1'

        k = 0;
        for (i_2 = 0; i_2 < 4; i_2++) {
          i = k + ii;
          rtb_y_m[i] = 0.0;
          rtb_y_m[i] += rtb_y_g[ii] * rtP.G_Value_a[i_2];
          rtb_y_m[i] += rtb_y_g[ii + 4] * rtP.G_Value_a[i_2 + 4];
          rtb_y_m[i] += rtb_y_g[ii + 8] * rtP.G_Value_a[i_2 + 8];
          rtb_y_m[i] += rtb_y_g[ii + 12] * rtP.G_Value_a[i_2 + 12];
          k += 4;
        }
      }

      // Product: '<S131>/Product' incorporates:
      //   Constant: '<S111>/G'

      ii = 0;
      for (k = 0; k < 4; k++) {
        for (i_2 = 0; i_2 <= 2; i_2 += 2) {
          i = i_2 + ii;
          (void)_mm_storeu_pd(&rtb_Q[i], _mm_set1_pd(0.0));
          tmp_3 = _mm_loadu_pd(&rtb_Q[i]);
          (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(tmp_3, _mm_mul_pd
            (_mm_set1_pd(rtb_y_m[ii]), _mm_loadu_pd(&rtP.G_Value_a[i_2]))));
          tmp_3 = _mm_loadu_pd(&rtb_Q[i]);
          (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
            (rtb_y_m[ii + 1]), _mm_loadu_pd(&rtP.G_Value_a[i_2 + 4])), tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_Q[i]);
          (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
            (rtb_y_m[ii + 2]), _mm_loadu_pd(&rtP.G_Value_a[i_2 + 8])), tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_Q[i]);
          (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
            (rtb_y_m[ii + 3]), _mm_loadu_pd(&rtP.G_Value_a[i_2 + 12])), tmp_3));
        }

        ii += 4;
      }

      // Math: '<S131>/Transpose2' incorporates:
      //   Constant: '<S111>/H'

      ii = 0;
      for (k = 0; k < 3; k++) {
        rtb_Transpose2[ii] = rtP.H_Value_o[k];
        rtb_Transpose2[ii + 1] = rtP.H_Value_o[k + 3];
        rtb_Transpose2[ii + 2] = rtP.H_Value_o[k + 6];
        rtb_Transpose2[ii + 3] = rtP.H_Value_o[k + 9];
        ii += 4;
      }

      // End of Math: '<S131>/Transpose2'
      for (ii = 0; ii < 4; ii++) {
        // Sum: '<S131>/Add' incorporates:
        //   Math: '<S131>/Transpose2'
        //   Product: '<S131>/Product1'

        k = 0;
        for (i_2 = 0; i_2 < 3; i_2++) {
          i = k + ii;
          rtb_Add_e[i] = (((rtb_Transpose2[k + 1] * rtb_y_g[ii + 4] +
                            rtb_Transpose2[k] * rtb_y_g[ii]) + rtb_Transpose2[k
                           + 2] * rtb_y_g[ii + 8]) + rtb_Transpose2[k + 3] *
                          rtb_y_g[ii + 12]) + rtb_N[i];
          k += 4;
        }

        // End of Sum: '<S131>/Add'
      }

      // Product: '<S131>/Product2' incorporates:
      //   Constant: '<S111>/G'
      //   Product: '<S131>/Product3'
      //   Product: '<S131>/Product4'
      //   Sum: '<S131>/Add'

      ii = 0;
      for (k = 0; k < 3; k++) {
        for (i_2 = 0; i_2 <= 2; i_2 += 2) {
          i = i_2 + ii;
          (void)_mm_storeu_pd(&rtb_Product2[i], _mm_set1_pd(0.0));
          tmp_3 = _mm_loadu_pd(&rtb_Product2[i]);
          (void)_mm_storeu_pd(&rtb_Product2[i], _mm_add_pd(tmp_3, _mm_mul_pd
            (_mm_set1_pd(rtb_Add_e[ii]), _mm_loadu_pd(&rtP.G_Value_a[i_2]))));
          tmp_3 = _mm_loadu_pd(&rtb_Product2[i]);
          (void)_mm_storeu_pd(&rtb_Product2[i], _mm_add_pd(_mm_mul_pd
            (_mm_set1_pd(rtb_Add_e[ii + 1]), _mm_loadu_pd(&rtP.G_Value_a[i_2 + 4])),
            tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_Product2[i]);
          (void)_mm_storeu_pd(&rtb_Product2[i], _mm_add_pd(_mm_mul_pd
            (_mm_set1_pd(rtb_Add_e[ii + 2]), _mm_loadu_pd(&rtP.G_Value_a[i_2 + 8])),
            tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_Product2[i]);
          (void)_mm_storeu_pd(&rtb_Product2[i], _mm_add_pd(_mm_mul_pd
            (_mm_set1_pd(rtb_Add_e[ii + 3]), _mm_loadu_pd(&rtP.G_Value_a[i_2 +
            12])), tmp_3));
        }

        ii += 4;
      }

      // Product: '<S131>/Product3' incorporates:
      //   Product: '<S131>/Product4'

      ii = 0;
      i_2 = 0;
      for (i = 0; i < 3; i++) {
        // Product: '<S131>/Product2' incorporates:
        //   Product: '<S131>/Product3'
        //   Product: '<S131>/Product4'

        i_0 = 0;
        for (k = 0; k < 3; k++) {
          // Product: '<S131>/Product3' incorporates:
          //   Product: '<S131>/Product2'

          rtb_Q_tmp = k + ii;
          tmp[rtb_Q_tmp] = 0.0;

          // Product: '<S131>/Product4'
          rtb_N_0[rtb_Q_tmp] = 0.0;

          // Product: '<S131>/Product3' incorporates:
          //   Constant: '<S111>/H'
          //   Product: '<S131>/Product2'
          //   Sum: '<S131>/Add'

          tmp[rtb_Q_tmp] += rtb_Add_e[i_2] * rtP.H_Value_o[k];

          // Product: '<S131>/Product4' incorporates:
          //   Math: '<S131>/Transpose'
          //   Math: '<S131>/Transpose2'
          //   Product: '<S131>/Product2'
          //   Product: '<S131>/Product3'

          rtb_N_0[rtb_Q_tmp] += rtb_N[i_0] * rtb_Transpose2[i_2];

          // Product: '<S131>/Product3' incorporates:
          //   Constant: '<S111>/H'
          //   Product: '<S131>/Product2'
          //   Sum: '<S131>/Add'

          tmp[rtb_Q_tmp] += rtb_Add_e[i_2 + 1] * rtP.H_Value_o[k + 3];

          // Product: '<S131>/Product4' incorporates:
          //   Math: '<S131>/Transpose'
          //   Math: '<S131>/Transpose2'
          //   Product: '<S131>/Product2'
          //   Product: '<S131>/Product3'

          rtb_N_0[rtb_Q_tmp] += rtb_N[i_0 + 1] * rtb_Transpose2[i_2 + 1];

          // Product: '<S131>/Product3' incorporates:
          //   Constant: '<S111>/H'
          //   Product: '<S131>/Product2'
          //   Sum: '<S131>/Add'

          tmp[rtb_Q_tmp] += rtb_Add_e[i_2 + 2] * rtP.H_Value_o[k + 6];

          // Product: '<S131>/Product4' incorporates:
          //   Math: '<S131>/Transpose'
          //   Math: '<S131>/Transpose2'
          //   Product: '<S131>/Product2'
          //   Product: '<S131>/Product3'

          rtb_N_0[rtb_Q_tmp] += rtb_N[i_0 + 2] * rtb_Transpose2[i_2 + 2];

          // Product: '<S131>/Product3' incorporates:
          //   Constant: '<S111>/H'
          //   Product: '<S131>/Product2'
          //   Sum: '<S131>/Add'

          tmp[rtb_Q_tmp] += rtb_Add_e[i_2 + 3] * rtP.H_Value_o[k + 9];

          // Product: '<S131>/Product4' incorporates:
          //   Math: '<S131>/Transpose'
          //   Math: '<S131>/Transpose2'
          //   Product: '<S131>/Product2'
          //   Product: '<S131>/Product3'

          rtb_N_0[rtb_Q_tmp] += rtb_N[i_0 + 3] * rtb_Transpose2[i_2 + 3];

          // Product: '<S131>/Product3' incorporates:
          //   Product: '<S131>/Product2'
          //   Product: '<S131>/Product4'

          i_0 += 4;
        }

        // Product: '<S131>/Product2' incorporates:
        //   Product: '<S131>/Product3'
        //   Product: '<S131>/Product4'

        ii += 3;
        i_2 += 4;
      }

      // Sum: '<S131>/Add1' incorporates:
      //   MATLAB Function: '<S134>/ScalarExpansion'

      ii = 0;
      for (k = 0; k < 3; k++) {
        // Outputs for Atomic SubSystem: '<S111>/ScalarExpansionR'
        rtb_Add1[ii] = (rtb_R[ii] + rtb_R[k]) / 2.0 + (tmp[ii] + rtb_N_0[ii]);
        rtb_Add1[ii + 1] = (rtb_R[ii + 1] + rtb_R[k + 3]) / 2.0 + (tmp[ii + 1] +
          rtb_N_0[ii + 1]);
        rtb_Add1[ii + 2] = (rtb_R[ii + 2] + rtb_R[k + 6]) / 2.0 + (tmp[ii + 2] +
          rtb_N_0[ii + 2]);

        // End of Outputs for SubSystem: '<S111>/ScalarExpansionR'
        ii += 3;
      }

      // End of Sum: '<S131>/Add1'
      // End of Outputs for SubSystem: '<S111>/ReducedQRN'

      // Gain: '<S89>/umin_scale1'
      //  See help of ctrlKalmanFilterDTCalculatePL.m
      // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S151>:1' 
      //    Copyright 2014 The MathWorks, Inc.
      // '<S151>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
      //  Determine if the Square-Root algorithm was used
      // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S153>:1' 
      // '<S153>:1:4' if isSqrtUsed
      umin_scale1_idx_0 = rtP.umin_scale1_Gain[0] * umax_incr[0];

      // Sum: '<S88>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      umin_scale1_idx_1 = umin_scale1_idx_0 - rtU.u0[0];

      // End of Outputs for SubSystem: '<S1>/mpc1'

      // Outport: '<Root>/u'
      rtY.u[0] = umin_scale1_idx_0;

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      // Gain: '<S89>/umin_scale1'
      umin_scale1_idx_0 = rtP.umin_scale1_Gain[1] * umax_incr[1];

      // Sum: '<S88>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      umin_scale1_idx_2 = umin_scale1_idx_0 - rtU.u0[1];

      // End of Outputs for SubSystem: '<S1>/mpc1'

      // Outport: '<Root>/u'
      rtY.u[1] = umin_scale1_idx_0;

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      // Gain: '<S89>/umin_scale1'
      umin_scale1_idx_0 = rtP.umin_scale1_Gain[2] * umax_incr[2];

      // Sum: '<S88>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      dwt = umin_scale1_idx_0 - rtU.u0[2];

      // Outputs for Enabled SubSystem: '<S130>/MeasurementUpdate' incorporates:
      //   EnablePort: '<S154>/Enable'

      // Outputs for Atomic SubSystem: '<S111>/CalculatePL'
      // MATLAB Function: '<S113>/Discrete-Time KF - Calculate PLMZ' incorporates:
      //   Constant: '<S3>/Constant1'
      //   Constant: '<S3>/Constant13'
      //   DataTypeConversion: '<S111>/DataTypeConversionEnable'
      //   Delay: '<S111>/MemoryP'
      //   Delay: '<S111>/MemoryX'
      //   Product: '<S131>/Product'
      //   Product: '<S131>/Product2'
      //   Product: '<S154>/C[k]*xhat[k|k-1]'
      //   Product: '<S154>/D[k]*u[k]'
      //   Product: '<S154>/Product3'
      //   Sum: '<S131>/Add1'
      //   Sum: '<S154>/Add1'
      //   Sum: '<S154>/Sum'
      //   Sum: '<S88>/Sum6'

      if (rtP.Constant1_Value_e != 0.0) {
        ii = 0;
        for (k = 0; k < 3; k++) {
          i_2 = 0;
          i = 0;
          for (i_0 = 0; i_0 < 4; i_0++) {
            rtb_Q_tmp = i_2 + k;
            rtb_Transpose2[i_0 + ii] = rtb_C[rtb_Q_tmp];
            rtb_N[rtb_Q_tmp] = 0.0;
            rtb_N[rtb_Q_tmp] += rtDW.MemoryP_DSTATE[i] * rtb_C[k];
            rtb_N[rtb_Q_tmp] += rtDW.MemoryP_DSTATE[i + 1] * rtb_C[k + 3];
            rtb_N[rtb_Q_tmp] += rtDW.MemoryP_DSTATE[i + 2] * rtb_C[k + 6];
            rtb_N[rtb_Q_tmp] += rtDW.MemoryP_DSTATE[i + 3] * rtb_C[k + 9];
            i_2 += 3;
            i += 4;
          }

          ii += 4;
        }

        for (ii = 0; ii < 3; ii++) {
          k = 0;
          i_2 = 0;
          for (i = 0; i < 3; i++) {
            rtb_R_tmp_0 = k + ii;
            rtb_R[rtb_R_tmp_0] = (((rtb_Transpose2[i_2 + 1] * rtb_N[ii + 3] +
              rtb_Transpose2[i_2] * rtb_N[ii]) + rtb_Transpose2[i_2 + 2] *
              rtb_N[ii + 6]) + rtb_Transpose2[i_2 + 3] * rtb_N[ii + 9]) +
              rtb_Add1[rtb_R_tmp_0];
            k += 3;
            i_2 += 4;
          }
        }

        for (ii = 0; ii < 4; ii++) {
          for (k = 0; k < 4; k++) {
            i = k << 2UL;
            i_2 = ii + i;
            rtb_y_m[i_2] = 0.0;
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i] * rtb_A[ii];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 1] * rtb_A[ii + 4];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 2] * rtb_A[ii + 8];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 3] * rtb_A[ii + 12];
          }

          for (k = 0; k < 3; k++) {
            i_2 = k << 2UL;
            i = i_2 + ii;
            rtb_Add_e[i] = (((rtb_Transpose2[i_2 + 1] * rtb_y_m[ii + 4] +
                              rtb_Transpose2[i_2] * rtb_y_m[ii]) +
                             rtb_Transpose2[i_2 + 2] * rtb_y_m[ii + 8]) +
                            rtb_Transpose2[i_2 + 3] * rtb_y_m[ii + 12]) +
              rtb_Product2[i];
          }
        }

        mrdiv(rtb_Add_e, rtb_R, rtb_N);
        ii = 0;
        for (k = 0; k < 3; k++) {
          for (i_2 = 0; i_2 <= 2; i_2 += 2) {
            i = i_2 + ii;
            (void)_mm_storeu_pd(&rtb_Add_e[i], _mm_set1_pd(0.0));
            tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE[i_2]);
            tmp_2 = _mm_loadu_pd(&rtb_Add_e[i]);
            (void)_mm_storeu_pd(&rtb_Add_e[i], _mm_add_pd(tmp_2, _mm_mul_pd
              (_mm_set1_pd(rtb_Transpose2[ii]), tmp_3)));
            tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE[i_2 + 4]);
            tmp_2 = _mm_loadu_pd(&rtb_Add_e[i]);
            (void)_mm_storeu_pd(&rtb_Add_e[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_Transpose2[ii + 1]), tmp_3), tmp_2));
            tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE[i_2 + 8]);
            tmp_2 = _mm_loadu_pd(&rtb_Add_e[i]);
            (void)_mm_storeu_pd(&rtb_Add_e[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_Transpose2[ii + 2]), tmp_3), tmp_2));
            tmp_3 = _mm_loadu_pd(&rtDW.MemoryP_DSTATE[i_2 + 12]);
            tmp_2 = _mm_loadu_pd(&rtb_Add_e[i]);
            (void)_mm_storeu_pd(&rtb_Add_e[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
              (rtb_Transpose2[ii + 3]), tmp_3), tmp_2));
          }

          ii += 4;
        }

        mrdiv(rtb_Add_e, rtb_R, rtb_Transpose2);
        for (ii = 0; ii < 16; ii++) {
          b_I[ii] = 0;
        }

        b_I[0] = 1;
        b_I[5] = 1;
        b_I[10] = 1;
        b_I[15] = 1;
        for (ii = 0; ii < 4; ii++) {
          for (k = 0; k < 4; k++) {
            i = (k << 2UL) + ii;
            rtb_y_g[i] = static_cast<real_T>(b_I[i]) - ((rtb_C[3 * k + 1] *
              rtb_Transpose2[ii + 4] + rtb_C[3 * k] * rtb_Transpose2[ii]) +
              rtb_C[3 * k + 2] * rtb_Transpose2[ii + 8]);
          }

          for (k = 0; k < 4; k++) {
            i_2 = k << 2UL;
            i = ii + i_2;
            Abar[i] = 0.0;
            Abar[i] += rtDW.MemoryP_DSTATE[i_2] * rtb_y_g[ii];
            Abar[i] += rtDW.MemoryP_DSTATE[i_2 + 1] * rtb_y_g[ii + 4];
            Abar[i] += rtDW.MemoryP_DSTATE[i_2 + 2] * rtb_y_g[ii + 8];
            Abar[i] += rtDW.MemoryP_DSTATE[i_2 + 3] * rtb_y_g[ii + 12];
          }

          for (k = 0; k < 3; k++) {
            rtb_Q_tmp = (k << 2UL) + ii;
            rtb_Add_e[rtb_Q_tmp] = 0.0;
            rtb_Add_e[rtb_Q_tmp] += rtb_Add1[3 * k] * rtb_Transpose2[ii];
            rtb_Add_e[rtb_Q_tmp] += rtb_Add1[3 * k + 1] * rtb_Transpose2[ii + 4];
            rtb_Add_e[rtb_Q_tmp] += rtb_Add1[3 * k + 2] * rtb_Transpose2[ii + 8];
          }
        }

        for (ii = 0; ii < 4; ii++) {
          k = 0;
          for (i_2 = 0; i_2 < 4; i_2++) {
            i = k + ii;
            rtb_y_m[i] = 0.0;
            rtb_y_m[i] += Abar[ii] * rtb_y_g[i_2];
            rtb_y_m[i] += Abar[ii + 4] * rtb_y_g[i_2 + 4];
            rtb_y_m[i] += Abar[ii + 8] * rtb_y_g[i_2 + 8];
            rtb_y_m[i] += Abar[ii + 12] * rtb_y_g[i_2 + 12];
            rtb_Transpose2_0[i] = 0.0;
            rtb_Transpose2_0[i] += rtb_Add_e[ii] * rtb_Transpose2[i_2];
            rtb_Transpose2_0[i] += rtb_Add_e[ii + 4] * rtb_Transpose2[i_2 + 4];
            rtb_Transpose2_0[i] += rtb_Add_e[ii + 8] * rtb_Transpose2[i_2 + 8];
            k += 4;
          }
        }

        for (ii = 0; ii <= 14; ii += 2) {
          tmp_3 = _mm_loadu_pd(&rtb_y_m[ii]);
          tmp_2 = _mm_loadu_pd(&rtb_Transpose2_0[ii]);
          (void)_mm_storeu_pd(&rtb_Z[ii], _mm_add_pd(tmp_3, tmp_2));
        }

        mrdiv(rtb_Product2, rtb_Add1, rtb_Transpose2);
        for (ii = 0; ii < 4; ii++) {
          for (k = 0; k < 4; k++) {
            i = (k << 2UL) + ii;
            rtb_y_g[i] = rtb_A[i] - ((rtb_C[3 * k + 1] * rtb_Transpose2[ii + 4]
              + rtb_C[3 * k] * rtb_Transpose2[ii]) + rtb_C[3 * k + 2] *
              rtb_Transpose2[ii + 8]);
          }

          for (k = 0; k < 4; k++) {
            i_2 = k << 2UL;
            i = ii + i_2;
            Abar[i] = 0.0;
            Abar[i] += rtb_Z[i_2] * rtb_y_g[ii];
            Abar[i] += rtb_Z[i_2 + 1] * rtb_y_g[ii + 4];
            Abar[i] += rtb_Z[i_2 + 2] * rtb_y_g[ii + 8];
            Abar[i] += rtb_Z[i_2 + 3] * rtb_y_g[ii + 12];
          }
        }

        for (ii = 0; ii < 4; ii++) {
          k = 0;
          for (i_2 = 0; i_2 < 4; i_2++) {
            i = k + ii;
            rtb_y_m[i] = (((Abar[ii + 4] * rtb_y_g[i_2 + 4] + Abar[ii] *
                            rtb_y_g[i_2]) + Abar[ii + 8] * rtb_y_g[i_2 + 8]) +
                          Abar[ii + 12] * rtb_y_g[i_2 + 12]) + rtb_Q[i];
            rtb_Transpose2_0[i] = 0.0;
            rtb_Transpose2_0[i] += rtb_Transpose2[ii] * rtb_Product2[i_2];
            rtb_Transpose2_0[i] += rtb_Transpose2[ii + 4] * rtb_Product2[i_2 + 4];
            rtb_Transpose2_0[i] += rtb_Transpose2[ii + 8] * rtb_Product2[i_2 + 8];
            k += 4;
          }
        }

        for (ii = 0; ii <= 14; ii += 2) {
          tmp_3 = _mm_loadu_pd(&rtb_y_m[ii]);
          tmp_2 = _mm_loadu_pd(&rtb_Transpose2_0[ii]);
          (void)_mm_storeu_pd(&rtb_y_g[ii], _mm_sub_pd(tmp_3, tmp_2));
        }

        rtDW.MeasurementUpdate_MODE = true;
        for (ii = 0; ii <= 0; ii += 2) {
          tmp_3 = _mm_set1_pd(0.0);
          (void)_mm_storeu_pd(&rtb_C_0[ii], tmp_3);
          tmp_2 = _mm_loadu_pd(&rtb_C[ii]);
          tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
          (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
            _mm_set1_pd(rtDW.MemoryX_DSTATE[0])), tmp_1));
          tmp_2 = _mm_loadu_pd(&rtb_C[ii + 3]);
          tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
          (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
            _mm_set1_pd(rtDW.MemoryX_DSTATE[1])), tmp_1));
          tmp_2 = _mm_loadu_pd(&rtb_C[ii + 6]);
          tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
          (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
            _mm_set1_pd(rtDW.MemoryX_DSTATE[2])), tmp_1));
          tmp_2 = _mm_loadu_pd(&rtb_C[ii + 9]);
          tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
          (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
            _mm_set1_pd(rtDW.MemoryX_DSTATE[3])), tmp_1));
          (void)_mm_storeu_pd(&tmp_0[ii], tmp_3);
          tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
          (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[ii]), _mm_set1_pd(umin_scale1_idx_1)),
            tmp_3));
          tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
          (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[ii + 3]), _mm_set1_pd(umin_scale1_idx_2)),
            tmp_3));
          tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
          (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
            (&rtP.Constant13_Value_c[ii + 6]), _mm_set1_pd(dwt)), tmp_3));
          tmp_3 = _mm_loadu_pd(&rtb_C_0[ii]);
          tmp_2 = _mm_loadu_pd(&tmp_0[ii]);
          tmp_1 = _mm_loadu_pd(&rtb_y_[ii]);
          (void)_mm_storeu_pd(&rtb_Add1_n[ii], _mm_sub_pd(tmp_1, _mm_add_pd
            (tmp_3, tmp_2)));
        }

        for (ii = 2; ii < 3; ii++) {
          rtb_C_0[ii] = 0.0;
          rtb_C_0[ii] += rtb_C[ii] * rtDW.MemoryX_DSTATE[0];
          rtb_C_0[ii] += rtb_C[ii + 3] * rtDW.MemoryX_DSTATE[1];
          rtb_C_0[ii] += rtb_C[ii + 6] * rtDW.MemoryX_DSTATE[2];
          rtb_C_0[ii] += rtb_C[ii + 9] * rtDW.MemoryX_DSTATE[3];
          tmp_0[ii] = 0.0;
          tmp_0[ii] += rtP.Constant13_Value_c[ii] * umin_scale1_idx_1;
          tmp_0[ii] += rtP.Constant13_Value_c[ii + 3] * umin_scale1_idx_2;
          tmp_0[ii] += rtP.Constant13_Value_c[ii + 6] * dwt;
          rtb_Add1_n[ii] = rtb_y_[ii] - (rtb_C_0[ii] + tmp_0[ii]);
        }

        for (ii = 0; ii <= 2; ii += 2) {
          (void)_mm_storeu_pd(&rtDW.Product3[ii], _mm_set1_pd(0.0));
          tmp_3 = _mm_loadu_pd(&rtb_N[ii]);
          tmp_2 = _mm_loadu_pd(&rtDW.Product3[ii]);
          (void)_mm_storeu_pd(&rtDW.Product3[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
            _mm_set1_pd(rtb_Add1_n[0])), tmp_2));
          tmp_3 = _mm_loadu_pd(&rtb_N[ii + 4]);
          tmp_2 = _mm_loadu_pd(&rtDW.Product3[ii]);
          (void)_mm_storeu_pd(&rtDW.Product3[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
            _mm_set1_pd(rtb_Add1_n[1])), tmp_2));
          tmp_3 = _mm_loadu_pd(&rtb_N[ii + 8]);
          tmp_2 = _mm_loadu_pd(&rtDW.Product3[ii]);
          (void)_mm_storeu_pd(&rtDW.Product3[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
            _mm_set1_pd(rtb_Add1_n[2])), tmp_2));
        }
      } else {
        for (ii = 0; ii < 4; ii++) {
          for (k = 0; k < 4; k++) {
            i = k << 2UL;
            i_2 = ii + i;
            rtb_y_m[i_2] = 0.0;
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i] * rtb_A[ii];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 1] * rtb_A[ii + 4];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 2] * rtb_A[ii + 8];
            rtb_y_m[i_2] += rtDW.MemoryP_DSTATE[i + 3] * rtb_A[ii + 12];
          }

          for (k = 0; k < 4; k++) {
            i = (k << 2UL) + ii;
            rtb_y_g[i] = (((rtb_y_m[ii + 4] * rtb_A[k + 4] + rtb_y_m[ii] *
                            rtb_A[k]) + rtb_y_m[ii + 8] * rtb_A[k + 8]) +
                          rtb_y_m[ii + 12] * rtb_A[k + 12]) + rtb_Q[i];
          }
        }

        if (rtDW.MeasurementUpdate_MODE) {
          // Disable for Product: '<S154>/Product3' incorporates:
          //   Outport: '<S154>/L*(y[k]-yhat[k|k-1])'
          //
          rtDW.Product3[0] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3[1] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3[2] = rtP.Lykyhatkk1_Y0_c;
          rtDW.Product3[3] = rtP.Lykyhatkk1_Y0_c;
          rtDW.MeasurementUpdate_MODE = false;
        }
      }

      // End of MATLAB Function: '<S113>/Discrete-Time KF - Calculate PLMZ'
      // End of Outputs for SubSystem: '<S111>/CalculatePL'
      // End of Outputs for SubSystem: '<S130>/MeasurementUpdate'
      // End of Outputs for SubSystem: '<S1>/mpc1'
      for (ii = 0; ii <= 0; ii += 2) {
        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Product: '<S114>/Product' incorporates:
        //   Delay: '<S111>/MemoryX'

        tmp_3 = _mm_set1_pd(0.0);
        (void)_mm_storeu_pd(&rtb_C_0[ii], tmp_3);
        tmp_2 = _mm_loadu_pd(&rtb_C[ii]);
        tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
        (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[0])), tmp_1));
        tmp_2 = _mm_loadu_pd(&rtb_C[ii + 3]);
        tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
        (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[1])), tmp_1));
        tmp_2 = _mm_loadu_pd(&rtb_C[ii + 6]);
        tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
        (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[2])), tmp_1));
        tmp_2 = _mm_loadu_pd(&rtb_C[ii + 9]);
        tmp_1 = _mm_loadu_pd(&rtb_C_0[ii]);
        (void)_mm_storeu_pd(&rtb_C_0[ii], _mm_add_pd(_mm_mul_pd(tmp_2,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[3])), tmp_1));

        // Product: '<S114>/Product1' incorporates:
        //   Constant: '<S3>/Constant13'
        //   Product: '<S114>/Product'

        (void)_mm_storeu_pd(&tmp_0[ii], tmp_3);
        tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
        (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
          (&rtP.Constant13_Value_c[ii]), _mm_set1_pd(umin_scale1_idx_1)), tmp_3));
        tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
        (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
          (&rtP.Constant13_Value_c[ii + 3]), _mm_set1_pd(umin_scale1_idx_2)),
          tmp_3));
        tmp_3 = _mm_loadu_pd(&tmp_0[ii]);
        (void)_mm_storeu_pd(&tmp_0[ii], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
          (&rtP.Constant13_Value_c[ii + 6]), _mm_set1_pd(dwt)), tmp_3));

        // Product: '<S114>/Product'
        tmp_3 = _mm_loadu_pd(&rtb_C_0[ii]);

        // Product: '<S114>/Product1' incorporates:
        //   Product: '<S114>/Product'

        tmp_2 = _mm_loadu_pd(&tmp_0[ii]);

        // Sum: '<S114>/Add1' incorporates:
        //   Product: '<S114>/Product'

        (void)_mm_storeu_pd(&rtb_Add1_n[ii], _mm_add_pd(tmp_3, tmp_2));

        // End of Outputs for SubSystem: '<S1>/mpc1'
      }

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      for (ii = 2; ii < 3; ii++) {
        // Product: '<S114>/Product' incorporates:
        //   Delay: '<S111>/MemoryX'

        rtb_C_0[ii] = 0.0;
        rtb_C_0[ii] += rtb_C[ii] * rtDW.MemoryX_DSTATE[0];
        rtb_C_0[ii] += rtb_C[ii + 3] * rtDW.MemoryX_DSTATE[1];
        rtb_C_0[ii] += rtb_C[ii + 6] * rtDW.MemoryX_DSTATE[2];
        rtb_C_0[ii] += rtb_C[ii + 9] * rtDW.MemoryX_DSTATE[3];

        // Product: '<S114>/Product1' incorporates:
        //   Constant: '<S3>/Constant13'
        //   Product: '<S114>/Product'

        tmp_0[ii] = 0.0;
        tmp_0[ii] += rtP.Constant13_Value_c[ii] * umin_scale1_idx_1;
        tmp_0[ii] += rtP.Constant13_Value_c[ii + 3] * umin_scale1_idx_2;
        tmp_0[ii] += rtP.Constant13_Value_c[ii + 6] * dwt;

        // Sum: '<S114>/Add1' incorporates:
        //   Product: '<S114>/Product'

        rtb_Add1_n[ii] = rtb_C_0[ii] + tmp_0[ii];
      }

      for (k = 0; k < 6; k++) {
        // Update for Delay: '<S3>/Delay'
        rtDW.Delay_DSTATE[k] = rtb_ywtT[k];
      }

      // Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
      //   Sum: '<S3>/Sum'

      rtDW.DiscreteTimeIntegrator_DSTATE += (rtb_y_[0] - Sum2_c[0]) *
        rtP.DiscreteTimeIntegrator_gainval;

      // Update for UnitDelay: '<S89>/last_mv'
      rtDW.last_mv_DSTATE[0] = umax_incr[0];
      rtDW.last_mv_DSTATE[1] = umax_incr[1];
      rtDW.last_mv_DSTATE[2] = umax_incr[2];

      // Update for Delay: '<S111>/MemoryX'
      rtDW.icLoad = false;

      // End of Outputs for SubSystem: '<S1>/mpc1'
      for (ii = 0; ii <= 2; ii += 2) {
        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Product: '<S130>/B[k]*u[k]'
        tmp_3 = _mm_set1_pd(0.0);
        (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFun_p[ii], tmp_3);
        tmp_2 = _mm_loadu_pd(&Sum_h[ii]);
        tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFun_p[ii]);
        (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFun_p[ii], _mm_add_pd
                            (_mm_mul_pd(tmp_2, _mm_set1_pd(umin_scale1_idx_1)),
                             tmp_1));
        tmp_2 = _mm_loadu_pd(&Sum_h[ii + 4]);
        tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFun_p[ii]);
        (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFun_p[ii], _mm_add_pd
                            (_mm_mul_pd(tmp_2, _mm_set1_pd(umin_scale1_idx_2)),
                             tmp_1));
        tmp_2 = _mm_loadu_pd(&Sum_h[ii + 8]);
        tmp_1 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFun_p[ii]);
        (void)_mm_storeu_pd(&rtb_TmpSignalConversionAtSFun_p[ii], _mm_add_pd
                            (_mm_mul_pd(tmp_2, _mm_set1_pd(dwt)), tmp_1));

        // Product: '<S130>/A[k]*xhat[k|k-1]' incorporates:
        //   Delay: '<S111>/MemoryX'
        //   Product: '<S130>/B[k]*u[k]'

        (void)_mm_storeu_pd(&rtb_Sum2[ii], tmp_3);
        tmp_3 = _mm_loadu_pd(&rtb_A[ii]);
        tmp_2 = _mm_loadu_pd(&rtb_Sum2[ii]);
        (void)_mm_storeu_pd(&rtb_Sum2[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[0])), tmp_2));
        tmp_3 = _mm_loadu_pd(&rtb_A[ii + 4]);
        tmp_2 = _mm_loadu_pd(&rtb_Sum2[ii]);
        (void)_mm_storeu_pd(&rtb_Sum2[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[1])), tmp_2));
        tmp_3 = _mm_loadu_pd(&rtb_A[ii + 8]);
        tmp_2 = _mm_loadu_pd(&rtb_Sum2[ii]);
        (void)_mm_storeu_pd(&rtb_Sum2[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[2])), tmp_2));
        tmp_3 = _mm_loadu_pd(&rtb_A[ii + 12]);
        tmp_2 = _mm_loadu_pd(&rtb_Sum2[ii]);
        (void)_mm_storeu_pd(&rtb_Sum2[ii], _mm_add_pd(_mm_mul_pd(tmp_3,
          _mm_set1_pd(rtDW.MemoryX_DSTATE[3])), tmp_2));

        // End of Outputs for SubSystem: '<S1>/mpc1'
      }

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      // Update for Delay: '<S111>/MemoryX' incorporates:
      //   Sum: '<S130>/Add'

      rtDW.MemoryX_DSTATE[0] = (rtb_TmpSignalConversionAtSFun_p[0] + rtb_Sum2[0])
        + rtDW.Product3[0];
      rtDW.MemoryX_DSTATE[1] = (rtb_TmpSignalConversionAtSFun_p[1] + rtb_Sum2[1])
        + rtDW.Product3[1];
      rtDW.MemoryX_DSTATE[2] = (rtb_TmpSignalConversionAtSFun_p[2] + rtb_Sum2[2])
        + rtDW.Product3[2];
      rtDW.MemoryX_DSTATE[3] = (rtb_TmpSignalConversionAtSFun_p[3] + rtb_Sum2[3])
        + rtDW.Product3[3];

      // Update for Delay: '<S111>/MemoryP'
      rtDW.icLoad_h = false;
      (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtb_y_g[0], sizeof(real_T) <<
                        4UL);

      // Outport: '<Root>/yhat' incorporates:
      //   Sum: '<S88>/Sum3'

      rtY.yhat[0] = rtb_Add1_n[0];
      rtY.yhat[1] = rtb_Add1_n[1];

      // End of Outputs for SubSystem: '<S1>/mpc1'

      // Outport: '<Root>/u'
      rtY.u[2] = umin_scale1_idx_0;

      // Outputs for Function Call SubSystem: '<S1>/mpc1'
      // Outport: '<Root>/yhat' incorporates:
      //   Sum: '<S88>/Sum3'

      rtY.yhat[2] = rtb_Add1_n[2];

      // End of Outputs for SubSystem: '<S1>/mpc1'
      for (k = 0; k < 6; k++) {
        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Outport: '<Root>/ywt' incorporates:
        //   Gain: '<S3>/Gain'

        rtY.ywt[k] = rtP.beta * rtb_ywt[k];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/currTraj'
        rtY.currTraj[k] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6 +
          k];
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
      rtY.P_g[i] = static_cast<real_T>(tmp_3[i]);
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

    // InitializeConditions for Memory: '<S89>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput[0],
                      &rtP.Memory_InitialCondition_f[0], 166U * sizeof(boolean_T));
    for (i = 0; i < 6; i++) {
      // InitializeConditions for Delay: '<S3>/Delay'
      rtDW.Delay_DSTATE[i] = rtP.ywt0[i];
    }

    // InitializeConditions for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE = rtP.DiscreteTimeIntegrator_IC;

    // InitializeConditions for UnitDelay: '<S89>/last_mv'
    rtDW.last_mv_DSTATE[0] = rtP.last_mv_InitialCondition_f[0];
    rtDW.last_mv_DSTATE[1] = rtP.last_mv_InitialCondition_f[1];
    rtDW.last_mv_DSTATE[2] = rtP.last_mv_InitialCondition_f[2];

    // InitializeConditions for Delay: '<S111>/MemoryX'
    rtDW.icLoad = true;

    // InitializeConditions for Delay: '<S111>/MemoryP'
    rtDW.icLoad_h = true;

    // SystemInitialize for Enabled SubSystem: '<S130>/MeasurementUpdate'
    // SystemInitialize for Product: '<S154>/Product3' incorporates:
    //   Outport: '<S154>/L*(y[k]-yhat[k|k-1])'

    rtDW.Product3[0] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3[1] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3[2] = rtP.Lykyhatkk1_Y0_c;
    rtDW.Product3[3] = rtP.Lykyhatkk1_Y0_c;

    // End of SystemInitialize for SubSystem: '<S130>/MeasurementUpdate'
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
