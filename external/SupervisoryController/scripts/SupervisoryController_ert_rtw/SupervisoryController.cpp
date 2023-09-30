//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2516
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sat Sep 30 10:58:52 2023
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

// Named constants for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

  // MATLAB Function: '<S322>/MATLAB Function1'
  in3_idx_0 = in3 - in2;
  stride_0_0 = (in7 - in6) + 1 != 1 ? static_cast<int32_T>(1) : static_cast<
    int32_T>(0);
  for (int32_T i{0}; i < in3_idx_0; i++) {
    in1[(in2 + i) + 12 * in4] = in5[i * stride_0_0 + in6] * in8[(in4 << 2UL) + i];
  }

  // End of MATLAB Function: '<S322>/MATLAB Function1'
}

//
// System initialize for function-call system:
//    '<S1>/paramEst1'
//    '<S1>/paramEst2'
//
void SupervisoryController::paramEst1_Init(real_T rty_theta[12], real_T rty_P
  [144], real_T rty_err[3], DW_paramEst1 *localDW, P_paramEst1 *localP)
{
  // InitializeConditions for Delay: '<S324>/Delay1'
  localDW->icLoad = true;

  // InitializeConditions for UnitDelay: '<S322>/Unit Delay3'
  localDW->UnitDelay3_DSTATE[0] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[1] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[2] = localP->UnitDelay3_InitialCondition;

  // InitializeConditions for Delay: '<S324>/Delay'
  localDW->icLoad_n = true;

  // SystemInitialize for Outport: '<S7>/theta'
  for (int32_T i{0}; i < 12; i++) {
    rty_theta[i] = localP->theta_Y0;
  }

  // End of SystemInitialize for Outport: '<S7>/theta'

  // SystemInitialize for Outport: '<S7>/P'
  for (int32_T i{0}; i < 144; i++) {
    rty_P[i] = localP->P_Y0;
  }

  // End of SystemInitialize for Outport: '<S7>/P'

  // SystemInitialize for Outport: '<S7>/err'
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

  // Delay: '<S324>/Delay1'
  localDW->icLoad = ((rtu_rstP && (static_cast<uint32_T>
    (localZCE->Delay1_Reset_ZCE) != POS_ZCSIG)) || localDW->icLoad);
  localZCE->Delay1_Reset_ZCE = rtu_rstP ? static_cast<ZCSigState>(1) :
    static_cast<ZCSigState>(0);
  if (localDW->icLoad) {
    (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rtu_P0[0], 144U * sizeof
                      (real_T));
  }

  // Sum: '<S322>/Add1'
  rtb_Add1_pv[0] = rtu_y[0] - rtu_y0[0];

  // Sum: '<S322>/Add3'
  rtb_Add3_idx_0 = rtu_u[0] - rtu_u0[0];

  // Sum: '<S322>/Add1'
  rtb_Add1_pv[1] = rtu_y[1] - rtu_y0[1];

  // Sum: '<S322>/Add3'
  rtb_Add3_idx_1 = rtu_u[1] - rtu_u0[1];

  // Sum: '<S322>/Add1'
  rtb_Add1_pv[2] = rtu_y[2] - rtu_y0[2];

  // Sum: '<S322>/Add3'
  rtb_Add3_idx_2 = rtu_u[2] - rtu_u0[2];

  // MATLAB Function: '<S322>/MATLAB Function1' incorporates:
  //   Sum: '<S322>/Add1'
  //   Sum: '<S322>/Add3'

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1': '<S323>:1' 
  // '<S323>:1:2' [z, phi] = getRegressors_(y, yPrev, u, sign_, no, ni, np, dt, mdlNum); 
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

  // Delay: '<S324>/Delay'
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
    // Math: '<S324>/Transpose' incorporates:
    //   MATLAB Function: '<S324>/MATLAB Function'

    tmp[3 * i] = rtb_phi[i];
    tmp[3 * i + 1] = rtb_phi[i + 12];
    tmp[3 * i + 2] = rtb_phi[i + 24];
  }

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function': '<S325>:1' 
  // '<S325>:1:2' [dtheta, dP, L] = rls_(theta, phi, epsil, EN, p_, dPmod_, lambda, P, no, ni, np); 
  // 'rls_:3' dtheta = zeros(no*np, 1);
  // 'rls_:4' dP = zeros(no*np,no*np);
  // 'rls_:5' L = zeros(no*np, no);
  // 'rls_:7' L = P*phi*inv(lambda*eye(no) + phi'*P*phi);
  for (i = 0; i < 3; i++) {
    // Sum: '<S324>/Sum2' incorporates:
    //   Delay: '<S324>/Delay'
    //   MATLAB Function: '<S322>/MATLAB Function1'
    //   Math: '<S324>/Transpose'
    //   Sum: '<S322>/Add1'
    //   UnitDelay: '<S322>/Unit Delay3'

    rtb_Add3_idx_0 = 0.0;
    for (p2 = 0; p2 < 12; p2++) {
      // MATLAB Function: '<S324>/MATLAB Function'
      p1 = 3 * p2 + i;
      rtb_Add3_idx_0 += tmp[p1] * localDW->Delay_DSTATE[p2];

      // MATLAB Function: '<S324>/MATLAB Function' incorporates:
      //   Delay: '<S324>/Delay'
      //   Delay: '<S324>/Delay1'

      tmp_0[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        tmp_0[p1] += tmp[3 * ibtile + i] * localDW->Delay1_DSTATE[12 * p2 +
          ibtile];
      }
    }

    rty_err[i] = (rtb_Add1_pv[i] - localDW->UnitDelay3_DSTATE[i]) -
      rtb_Add3_idx_0;

    // End of Sum: '<S324>/Sum2'

    // MATLAB Function: '<S324>/MATLAB Function'
    for (p2 = 0; p2 < 3; p2++) {
      rtb_Add3_idx_0 = 0.0;
      for (p1 = 0; p1 < 12; p1++) {
        rtb_Add3_idx_0 += tmp_0[3 * p1 + i] * rtb_phi[12 * p2 + p1];
      }

      p1 = 3 * p2 + i;
      b_b[p1] = static_cast<real_T>(b_b_0[p1]) * rtu_lambda + rtb_Add3_idx_0;
    }
  }

  // MATLAB Function: '<S324>/MATLAB Function' incorporates:
  //   Delay: '<S324>/Delay'
  //   Delay: '<S324>/Delay1'

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
    // Delay: '<S324>/Delay1'
    tmp_1 = _mm_loadu_pd(&localDW->Delay1_DSTATE[i]);

    // Product: '<S324>/Product1' incorporates:
    //   Delay: '<S324>/Delay1'

    tmp_2 = _mm_loadu_pd(&rty_P[i]);
    (void)_mm_storeu_pd(&rty_P[i], _mm_div_pd(_mm_sub_pd(tmp_1, tmp_2),
      _mm_set1_pd(rtu_lambda)));
  }

  for (i = 0; i <= 10; i += 2) {
    // Delay: '<S324>/Delay'
    tmp_1 = _mm_loadu_pd(&localDW->Delay_DSTATE[i]);

    // Sum: '<S324>/Sum' incorporates:
    //   Delay: '<S324>/Delay'

    tmp_2 = _mm_loadu_pd(&rty_theta[i]);
    (void)_mm_storeu_pd(&rty_theta[i], _mm_add_pd(tmp_1, tmp_2));
  }

  // Update for Delay: '<S324>/Delay1'
  localDW->icLoad = false;
  (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rty_P[0], 144U * sizeof(real_T));

  // Update for UnitDelay: '<S322>/Unit Delay3' incorporates:
  //   Sum: '<S322>/Add1'

  localDW->UnitDelay3_DSTATE[0] = rtb_Add1_pv[0];
  localDW->UnitDelay3_DSTATE[1] = rtb_Add1_pv[1];
  localDW->UnitDelay3_DSTATE[2] = rtb_Add1_pv[2];

  // Update for Delay: '<S324>/Delay'
  localDW->icLoad_n = false;
  (void)std::memcpy(&localDW->Delay_DSTATE[0], &rty_theta[0], 12U * sizeof
                    (real_T));
}

//
// Output and update for atomic system:
//    '<S113>/ScalarExpansionR'
//    '<S183>/ScalarExpansionR'
//    '<S253>/ScalarExpansionR'
//
void SupervisoryController::ScalarExpansionR(const real_T rtu_u[9], real_T
  rty_y[9])
{
  int32_T tmp;

  // MATLAB Function: '<S136>/ScalarExpansion'
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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S158>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S158>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  tmp = 0;
  for (int32_T i{0}; i < 3; i++) {
    rty_y[tmp] = (rtu_u[tmp] + rtu_u[i]) / 2.0;
    rty_y[tmp + 1] = (rtu_u[tmp + 1] + rtu_u[i + 3]) / 2.0;
    rty_y[tmp + 2] = (rtu_u[tmp + 2] + rtu_u[i + 6]) / 2.0;
    tmp += 3;
  }

  // End of MATLAB Function: '<S136>/ScalarExpansion'
}

// Function for MATLAB Function: '<S185>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv(const real_T A[24], const real_T B_0[9],
  real_T Y[24])
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
  for (rtemp = 0; rtemp < 8; rtemp++) {
    int32_T Y_tmp;
    int32_T Y_tmp_0;
    int32_T Y_tmp_1;
    Y_tmp = (r1 << 3UL) + rtemp;
    Y[Y_tmp] = A[rtemp] / b_A[r1];
    Y_tmp_0 = (r2 << 3UL) + rtemp;
    Y[Y_tmp_0] = A[rtemp + 8] - b_A[r1 + 3] * Y[Y_tmp];
    Y_tmp_1 = (r3 << 3UL) + rtemp;
    Y[Y_tmp_1] = A[rtemp + 16] - b_A[r1 + 6] * Y[Y_tmp];
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
//    '<S183>/CalculatePL'
//    '<S253>/CalculatePL'
//
void SupervisoryController::CalculatePL(const real_T rtu_Ak[64], const real_T
  rtu_Ck[24], const real_T rtu_Qbark[64], const real_T rtu_Rbark[9], const
  real_T rtu_Nbark[24], boolean_T rtu_Enablek, const real_T rtu_Pk[64], real_T
  rty_Mk[24], real_T rty_Lk[24], real_T rty_Zk[64], real_T rty_Pk1[64])
{
  real_T Abar[64];
  real_T Abar_0[64];
  real_T Abar_1[64];
  real_T rty_Mk_0[64];
  real_T NRinv[24];
  real_T rtu_Ck_0[24];
  real_T yCov[9];
  int8_T b_I[64];

  // MATLAB Function: '<S185>/Discrete-Time KF - Calculate PLMZ'
  //  See help of ctrlKalmanFilterDTCalculatePL.m
  // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S223>:1' 
  //    Copyright 2014 The MathWorks, Inc.
  // '<S223>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
  if (rtu_Enablek) {
    __m128d tmp_0;
    __m128d tmp_1;
    real_T tmp;
    int32_T i;
    int32_T k;
    int32_T rtu_Ak_tmp;
    int32_T rtu_Ak_tmp_0;
    int32_T rtu_Pk_tmp;
    int32_T rty_Mk_tmp;
    i = 0;
    for (k = 0; k < 3; k++) {
      rty_Mk_tmp = 0;
      rtu_Ak_tmp = 0;
      for (rtu_Ak_tmp_0 = 0; rtu_Ak_tmp_0 < 8; rtu_Ak_tmp_0++) {
        int32_T NRinv_tmp;
        NRinv_tmp = rty_Mk_tmp + k;
        NRinv[rtu_Ak_tmp_0 + i] = rtu_Ck[NRinv_tmp];
        rtu_Ck_0[NRinv_tmp] = 0.0;
        rtu_Pk_tmp = 0;
        for (int32_T i_0{0}; i_0 < 8; i_0++) {
          rtu_Ck_0[NRinv_tmp] += rtu_Ck[rtu_Pk_tmp + k] * rtu_Pk[i_0 +
            rtu_Ak_tmp];
          rtu_Pk_tmp += 3;
        }

        rty_Mk_tmp += 3;
        rtu_Ak_tmp += 8;
      }

      i += 8;
    }

    for (i = 0; i < 3; i++) {
      k = 0;
      rty_Mk_tmp = 0;
      for (rtu_Ak_tmp = 0; rtu_Ak_tmp < 3; rtu_Ak_tmp++) {
        tmp = 0.0;
        rtu_Ak_tmp_0 = 0;
        for (rtu_Pk_tmp = 0; rtu_Pk_tmp < 8; rtu_Pk_tmp++) {
          tmp += rtu_Ck_0[rtu_Ak_tmp_0 + i] * NRinv[rtu_Pk_tmp + rty_Mk_tmp];
          rtu_Ak_tmp_0 += 3;
        }

        rtu_Ak_tmp_0 = k + i;
        yCov[rtu_Ak_tmp_0] = rtu_Rbark[rtu_Ak_tmp_0] + tmp;
        k += 3;
        rty_Mk_tmp += 8;
      }
    }

    for (i = 0; i < 8; i++) {
      for (k = 0; k < 8; k++) {
        rtu_Ak_tmp = k << 3UL;
        rtu_Ak_tmp_0 = i + rtu_Ak_tmp;
        Abar[rtu_Ak_tmp_0] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          Abar[rtu_Ak_tmp_0] += rtu_Ak[(rty_Mk_tmp << 3UL) + i] *
            rtu_Pk[rtu_Ak_tmp + rty_Mk_tmp];
        }
      }

      for (k = 0; k < 3; k++) {
        tmp = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          tmp += Abar[(rty_Mk_tmp << 3UL) + i] * NRinv[(k << 3UL) + rty_Mk_tmp];
        }

        rtu_Ak_tmp = (k << 3UL) + i;
        rtu_Ck_0[rtu_Ak_tmp] = rtu_Nbark[rtu_Ak_tmp] + tmp;
      }
    }

    mrdiv(rtu_Ck_0, yCov, rty_Lk);
    i = 0;
    for (k = 0; k < 3; k++) {
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
        rtu_Pk_tmp = rty_Mk_tmp + i;
        rtu_Ck_0[rtu_Pk_tmp] = 0.0;
        rtu_Ak_tmp = 0;
        for (rtu_Ak_tmp_0 = 0; rtu_Ak_tmp_0 < 8; rtu_Ak_tmp_0++) {
          rtu_Ck_0[rtu_Pk_tmp] += rtu_Pk[rtu_Ak_tmp + rty_Mk_tmp] *
            NRinv[rtu_Ak_tmp_0 + i];
          rtu_Ak_tmp += 8;
        }
      }

      i += 8;
    }

    mrdiv(rtu_Ck_0, yCov, rty_Mk);
    (void)std::memset(&b_I[0], 0, sizeof(int8_T) << 6UL);
    k = 0;
    for (i = 0; i < 8; i++) {
      b_I[k] = 1;
      k += 9;
    }

    for (i = 0; i < 8; i++) {
      for (k = 0; k < 8; k++) {
        rtu_Pk_tmp = (k << 3UL) + i;
        Abar[rtu_Pk_tmp] = static_cast<real_T>(b_I[rtu_Pk_tmp]) - ((rtu_Ck[3 * k
          + 1] * rty_Mk[i + 8] + rtu_Ck[3 * k] * rty_Mk[i]) + rtu_Ck[3 * k + 2] *
          rty_Mk[i + 16]);
      }

      for (k = 0; k < 8; k++) {
        rtu_Pk_tmp = k << 3UL;
        rtu_Ak_tmp = i + rtu_Pk_tmp;
        Abar_0[rtu_Ak_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          Abar_0[rtu_Ak_tmp] += Abar[(rty_Mk_tmp << 3UL) + i] *
            rtu_Pk[rtu_Pk_tmp + rty_Mk_tmp];
        }
      }

      for (k = 0; k < 3; k++) {
        rty_Mk_tmp = (k << 3UL) + i;
        NRinv[rty_Mk_tmp] = 0.0;
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k] * rty_Mk[i];
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k + 1] * rty_Mk[i + 8];
        NRinv[rty_Mk_tmp] += rtu_Rbark[3 * k + 2] * rty_Mk[i + 16];
      }
    }

    for (i = 0; i < 8; i++) {
      k = 0;
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
        rtu_Pk_tmp = k + i;
        Abar_1[rtu_Pk_tmp] = 0.0;
        rtu_Ak_tmp = 0;
        for (rtu_Ak_tmp_0 = 0; rtu_Ak_tmp_0 < 8; rtu_Ak_tmp_0++) {
          Abar_1[rtu_Pk_tmp] += Abar_0[rtu_Ak_tmp + i] * Abar[rtu_Ak_tmp +
            rty_Mk_tmp];
          rtu_Ak_tmp += 8;
        }

        rty_Mk_0[rtu_Pk_tmp] = 0.0;
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i] * rty_Mk[rty_Mk_tmp];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 8] * rty_Mk[rty_Mk_tmp + 8];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 16] * rty_Mk[rty_Mk_tmp + 16];
        k += 8;
      }
    }

    for (i = 0; i <= 62; i += 2) {
      tmp_0 = _mm_loadu_pd(&Abar_1[i]);
      tmp_1 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Zk[i], _mm_add_pd(tmp_0, tmp_1));
    }

    mrdiv(rtu_Nbark, rtu_Rbark, NRinv);
    for (i = 0; i < 8; i++) {
      for (k = 0; k < 8; k++) {
        rtu_Pk_tmp = (k << 3UL) + i;
        Abar[rtu_Pk_tmp] = rtu_Ak[rtu_Pk_tmp] - ((rtu_Ck[3 * k + 1] * NRinv[i +
          8] + rtu_Ck[3 * k] * NRinv[i]) + rtu_Ck[3 * k + 2] * NRinv[i + 16]);
      }

      for (k = 0; k < 8; k++) {
        rtu_Pk_tmp = k << 3UL;
        rtu_Ak_tmp = i + rtu_Pk_tmp;
        Abar_0[rtu_Ak_tmp] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          Abar_0[rtu_Ak_tmp] += Abar[(rty_Mk_tmp << 3UL) + i] *
            rty_Zk[rtu_Pk_tmp + rty_Mk_tmp];
        }
      }
    }

    for (i = 0; i < 8; i++) {
      k = 0;
      for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
        tmp = 0.0;
        rtu_Ak_tmp = 0;
        for (rtu_Ak_tmp_0 = 0; rtu_Ak_tmp_0 < 8; rtu_Ak_tmp_0++) {
          tmp += Abar_0[rtu_Ak_tmp + i] * Abar[rtu_Ak_tmp + rty_Mk_tmp];
          rtu_Ak_tmp += 8;
        }

        rtu_Pk_tmp = k + i;
        Abar_1[rtu_Pk_tmp] = rtu_Qbark[rtu_Pk_tmp] + tmp;
        rty_Mk_0[rtu_Pk_tmp] = 0.0;
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i] * rtu_Nbark[rty_Mk_tmp];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 8] * rtu_Nbark[rty_Mk_tmp + 8];
        rty_Mk_0[rtu_Pk_tmp] += NRinv[i + 16] * rtu_Nbark[rty_Mk_tmp + 16];
        k += 8;
      }
    }

    for (i = 0; i <= 62; i += 2) {
      tmp_0 = _mm_loadu_pd(&Abar_1[i]);
      tmp_1 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Pk1[i], _mm_sub_pd(tmp_0, tmp_1));
    }
  } else {
    (void)std::memset(&rty_Lk[0], 0, 24U * sizeof(real_T));
    (void)std::memset(&rty_Mk[0], 0, 24U * sizeof(real_T));
    (void)std::memcpy(&rty_Zk[0], &rtu_Pk[0], sizeof(real_T) << 6UL);
    for (int32_T i{0}; i < 8; i++) {
      int32_T rtu_Ak_tmp;
      int32_T rty_Mk_tmp;
      for (int32_T k{0}; k < 8; k++) {
        int32_T rtu_Ak_tmp_0;
        rtu_Ak_tmp = k << 3UL;
        rtu_Ak_tmp_0 = i + rtu_Ak_tmp;
        Abar[rtu_Ak_tmp_0] = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          Abar[rtu_Ak_tmp_0] += rtu_Ak[(rty_Mk_tmp << 3UL) + i] *
            rtu_Pk[rtu_Ak_tmp + rty_Mk_tmp];
        }
      }

      for (int32_T k{0}; k < 8; k++) {
        real_T tmp;
        tmp = 0.0;
        for (rty_Mk_tmp = 0; rty_Mk_tmp < 8; rty_Mk_tmp++) {
          rtu_Ak_tmp = rty_Mk_tmp << 3UL;
          tmp += Abar[rtu_Ak_tmp + i] * rtu_Ak[rtu_Ak_tmp + k];
        }

        rty_Mk_tmp = (k << 3UL) + i;
        rty_Pk1[rty_Mk_tmp] = rtu_Qbark[rty_Mk_tmp] + tmp;
      }
    }
  }

  // End of MATLAB Function: '<S185>/Discrete-Time KF - Calculate PLMZ'
}

//
// Output and update for atomic system:
//    '<S224>/SqrtUsedFcn'
//    '<S294>/SqrtUsedFcn'
//
void SupervisoryController::SqrtUsedFcn(const real_T rtu_u[64], boolean_T
  rtu_isSqrtUsed, real_T rty_P[64])
{
  //  Determine if the Square-Root algorithm was used
  // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S225>:1' 
  // '<S225>:1:4' if isSqrtUsed
  if (rtu_isSqrtUsed) {
    // '<S225>:1:5' P = u*u.';
    for (int32_T i{0}; i < 8; i++) {
      int32_T tmp;
      tmp = 0;
      for (int32_T i_0{0}; i_0 < 8; i_0++) {
        int32_T tmp_0;
        int32_T tmp_1;
        tmp_1 = tmp + i;
        rty_P[tmp_1] = 0.0;
        tmp_0 = 0;
        for (int32_T i_1{0}; i_1 < 8; i_1++) {
          rty_P[tmp_1] += rtu_u[tmp_0 + i] * rtu_u[tmp_0 + i_0];
          tmp_0 += 8;
        }

        tmp += 8;
      }
    }
  } else {
    // '<S225>:1:6' else
    // '<S225>:1:7' P = u;
    (void)std::memcpy(&rty_P[0], &rtu_u[0], sizeof(real_T) << 6UL);
  }
}

//
// System initialize for enable system:
//    '<S202>/MeasurementUpdate'
//    '<S272>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Init(real_T rty_Lykyhatkk1[8],
  P_MeasurementUpdate *localP)
{
  // SystemInitialize for Outport: '<S226>/L*(y[k]-yhat[k|k-1])'
  for (int32_T i{0}; i < 8; i++) {
    rty_Lykyhatkk1[i] = localP->Lykyhatkk1_Y0;
  }

  // End of SystemInitialize for Outport: '<S226>/L*(y[k]-yhat[k|k-1])'
}

//
// Disable for enable system:
//    '<S202>/MeasurementUpdate'
//    '<S272>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[8],
  DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  // Disable for Outport: '<S226>/L*(y[k]-yhat[k|k-1])'
  for (int32_T i{0}; i < 8; i++) {
    // Outputs for Enabled SubSystem: '<S202>/MeasurementUpdate' incorporates:
    //   EnablePort: '<S226>/Enable'

    rty_Lykyhatkk1[i] = localP->Lykyhatkk1_Y0;

    // End of Outputs for SubSystem: '<S202>/MeasurementUpdate'
  }

  // End of Disable for Outport: '<S226>/L*(y[k]-yhat[k|k-1])'
  localDW->MeasurementUpdate_MODE = false;
}

//
// Output and update for enable system:
//    '<S202>/MeasurementUpdate'
//    '<S272>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate(boolean_T rtu_Enable, const real_T
  rtu_Lk[24], const real_T rtu_yk[3], const real_T rtu_Ck[24], const real_T
  rtu_xhatkk1[8], const real_T rtu_Dk[9], const real_T rtu_uk[3], real_T
  rty_Lykyhatkk1[8], DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  real_T rtu_Ck_0[3];
  real_T rtu_Dk_0[3];
  real_T rtu_yk_0[3];

  // Outputs for Enabled SubSystem: '<S202>/MeasurementUpdate' incorporates:
  //   EnablePort: '<S226>/Enable'

  if (rtu_Enable) {
    localDW->MeasurementUpdate_MODE = true;
    for (int32_T i{0}; i < 3; i++) {
      int32_T tmp;

      // Product: '<S226>/C[k]*xhat[k|k-1]'
      rtu_Ck_0[i] = 0.0;
      tmp = 0;
      for (int32_T i_0{0}; i_0 < 8; i_0++) {
        rtu_Ck_0[i] += rtu_Ck[tmp + i] * rtu_xhatkk1[i_0];
        tmp += 3;
      }

      // Product: '<S226>/D[k]*u[k]' incorporates:
      //   Product: '<S226>/C[k]*xhat[k|k-1]'

      rtu_Dk_0[i] = 0.0;
      rtu_Dk_0[i] += rtu_Dk[i] * rtu_uk[0];
      rtu_Dk_0[i] += rtu_Dk[i + 3] * rtu_uk[1];
      rtu_Dk_0[i] += rtu_Dk[i + 6] * rtu_uk[2];

      // Sum: '<S226>/Sum' incorporates:
      //   Product: '<S226>/C[k]*xhat[k|k-1]'
      //   Sum: '<S226>/Add1'

      rtu_yk_0[i] = rtu_yk[i] - (rtu_Ck_0[i] + rtu_Dk_0[i]);
    }

    for (int32_T i{0}; i <= 6; i += 2) {
      __m128d tmp_0;
      __m128d tmp_1;

      // Product: '<S226>/Product3'
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[0])), tmp_1));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i + 8]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[1])), tmp_1));
      tmp_0 = _mm_loadu_pd(&rtu_Lk[i + 16]);
      tmp_1 = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_0,
        _mm_set1_pd(rtu_yk_0[2])), tmp_1));
    }
  } else if (localDW->MeasurementUpdate_MODE) {
    MeasurementUpdate_Disable(rty_Lykyhatkk1, localDW, localP);
  } else {
    // no actions
  }

  // End of Outputs for SubSystem: '<S202>/MeasurementUpdate'
}

//
// Output and update for atomic system:
//    '<S183>/ReducedQRN'
//    '<S253>/ReducedQRN'
//
void SupervisoryController::ReducedQRN(const real_T rtu_G[64], const real_T
  rtu_H[24], const real_T rtu_Q[64], const real_T rtu_R[9], const real_T rtu_N
  [24], real_T rty_Qbar[64], real_T rty_Rbar[9], real_T rty_Nbar[24])
{
  real_T rtu_Q_0[64];
  real_T rtb_Add_a[24];
  real_T rtb_Transpose2_k[24];
  real_T rtu_H_0[9];
  real_T rtu_N_0[9];
  int32_T i;
  int32_T i_0;
  int32_T i_1;
  int32_T rtb_Add_ih_tmp;
  int32_T rtu_H_tmp;
  int32_T rtu_Q_tmp;

  // Product: '<S203>/Product' incorporates:
  //   Math: '<S203>/Transpose1'

  for (i = 0; i < 8; i++) {
    i_1 = 0;
    for (i_0 = 0; i_0 < 8; i_0++) {
      rtu_Q_tmp = i_1 + i;
      rtu_Q_0[rtu_Q_tmp] = 0.0;
      rtb_Add_ih_tmp = 0;
      for (rtu_H_tmp = 0; rtu_H_tmp < 8; rtu_H_tmp++) {
        rtu_Q_0[rtu_Q_tmp] += rtu_Q[rtb_Add_ih_tmp + i] * rtu_G[rtb_Add_ih_tmp +
          i_0];
        rtb_Add_ih_tmp += 8;
      }

      i_1 += 8;
    }
  }

  i = 0;
  for (i_1 = 0; i_1 < 8; i_1++) {
    for (i_0 = 0; i_0 < 8; i_0++) {
      rtb_Add_ih_tmp = i_0 + i;
      rty_Qbar[rtb_Add_ih_tmp] = 0.0;
      rtu_H_tmp = 0;
      for (rtu_Q_tmp = 0; rtu_Q_tmp < 8; rtu_Q_tmp++) {
        rty_Qbar[rtb_Add_ih_tmp] += rtu_G[rtu_H_tmp + i_0] * rtu_Q_0[rtu_Q_tmp +
          i];
        rtu_H_tmp += 8;
      }
    }

    i += 8;
  }

  // End of Product: '<S203>/Product'

  // Math: '<S203>/Transpose2'
  i = 0;
  for (i_1 = 0; i_1 < 3; i_1++) {
    i_0 = 0;
    for (rtb_Add_ih_tmp = 0; rtb_Add_ih_tmp < 8; rtb_Add_ih_tmp++) {
      rtb_Transpose2_k[rtb_Add_ih_tmp + i] = rtu_H[i_0 + i_1];
      i_0 += 3;
    }

    i += 8;
  }

  // End of Math: '<S203>/Transpose2'

  // Sum: '<S203>/Add' incorporates:
  //   Math: '<S203>/Transpose2'
  //   Product: '<S203>/Product1'

  for (i = 0; i < 8; i++) {
    i_1 = 0;
    for (i_0 = 0; i_0 < 3; i_0++) {
      real_T tmp;
      tmp = 0.0;
      rtb_Add_ih_tmp = 0;
      for (rtu_H_tmp = 0; rtu_H_tmp < 8; rtu_H_tmp++) {
        tmp += rtu_Q[rtb_Add_ih_tmp + i] * rtb_Transpose2_k[rtu_H_tmp + i_1];
        rtb_Add_ih_tmp += 8;
      }

      rtb_Add_ih_tmp = i_1 + i;
      rtb_Add_a[rtb_Add_ih_tmp] = rtu_N[rtb_Add_ih_tmp] + tmp;
      i_1 += 8;
    }
  }

  // End of Sum: '<S203>/Add'
  for (i = 0; i < 3; i++) {
    // Product: '<S203>/Product2' incorporates:
    //   Sum: '<S203>/Add'

    for (i_1 = 0; i_1 < 8; i_1++) {
      i_0 = i << 3UL;
      rtb_Add_ih_tmp = i_1 + i_0;
      rty_Nbar[rtb_Add_ih_tmp] = 0.0;
      for (rtu_H_tmp = 0; rtu_H_tmp < 8; rtu_H_tmp++) {
        rty_Nbar[rtb_Add_ih_tmp] += rtu_G[(rtu_H_tmp << 3UL) + i_1] *
          rtb_Add_a[i_0 + rtu_H_tmp];
      }
    }

    // End of Product: '<S203>/Product2'
    for (i_1 = 0; i_1 < 3; i_1++) {
      // Product: '<S203>/Product3' incorporates:
      //   Product: '<S203>/Product4'

      rtb_Add_ih_tmp = 3 * i_1 + i;
      rtu_H_0[rtb_Add_ih_tmp] = 0.0;

      // Product: '<S203>/Product4'
      rtu_N_0[rtb_Add_ih_tmp] = 0.0;
      for (i_0 = 0; i_0 < 8; i_0++) {
        // Product: '<S203>/Product3' incorporates:
        //   Product: '<S203>/Product4'
        //   Sum: '<S203>/Add'

        rtu_H_tmp = (i_1 << 3UL) + i_0;
        rtu_H_0[rtb_Add_ih_tmp] += rtu_H[3 * i_0 + i] * rtb_Add_a[rtu_H_tmp];

        // Product: '<S203>/Product4' incorporates:
        //   Math: '<S203>/Transpose'
        //   Math: '<S203>/Transpose2'

        rtu_N_0[rtb_Add_ih_tmp] += rtu_N[(i << 3UL) + i_0] *
          rtb_Transpose2_k[rtu_H_tmp];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
    __m128d tmp_0;
    __m128d tmp_1;
    __m128d tmp_2;

    // Sum: '<S203>/Add1'
    tmp_0 = _mm_loadu_pd(&rtu_H_0[i]);
    tmp_1 = _mm_loadu_pd(&rtu_N_0[i]);
    tmp_2 = _mm_loadu_pd(&rtu_R[i]);
    (void)_mm_storeu_pd(&rty_Rbar[i], _mm_add_pd(_mm_add_pd(tmp_0, tmp_1), tmp_2));
  }

  // Sum: '<S203>/Add1'
  for (i = 8; i < 9; i++) {
    rty_Rbar[i] = (rtu_H_0[i] + rtu_N_0[i]) + rtu_R[i];
  }
}

//
// Output and update for atomic system:
//    '<S183>/ScalarExpansionQ'
//    '<S253>/ScalarExpansionQ'
//
void SupervisoryController::ScalarExpansionQ(const real_T rtu_u[64], real_T
  rty_y[64])
{
  int32_T tmp;

  // MATLAB Function: '<S205>/ScalarExpansion'
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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S227>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S227>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  tmp = 0;
  for (int32_T i_0{0}; i_0 < 8; i_0++) {
    int32_T tmp_0;
    tmp_0 = 0;
    for (int32_T i{0}; i < 8; i++) {
      int32_T tmp_1;
      tmp_1 = i + tmp;
      rty_y[tmp_1] = (rtu_u[tmp_0 + i_0] + rtu_u[tmp_1]) / 2.0;
      tmp_0 += 8;
    }

    tmp += 8;
  }

  // End of MATLAB Function: '<S205>/ScalarExpansion'
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

// Function for MATLAB Function: '<S9>/MATLAB Function'
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

// Function for MATLAB Function: '<S111>/optimizer'
int32_T SupervisoryController::xpotrf(real_T b_A[16])
{
  int32_T info;
  int32_T j;
  boolean_T exitg1;
  info = 0;
  j = 0;
  exitg1 = false;
  while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
          false) && (j < 4)) {
    real_T c;
    real_T ssq;
    int32_T idxAjj;
    idxAjj = (j << 2UL) + j;
    ssq = 0.0;
    if (j >= 1) {
      for (int32_T b_k{0}; b_k < j; b_k++) {
        c = b_A[(b_k << 2UL) + j];
        ssq += c * c;
      }
    }

    ssq = b_A[idxAjj] - ssq;
    if (ssq > 0.0) {
      ssq = std::sqrt(ssq);
      b_A[idxAjj] = ssq;
      if (j + 1 < 4) {
        int32_T b_iy;
        int32_T d;
        int32_T jm1;
        if (j != 0) {
          b_iy = (((j - 1) << 2UL) + j) + 2;
          for (int32_T b_k{j + 2}; b_k <= b_iy; b_k += 4) {
            jm1 = b_k - j;
            c = -b_A[(((jm1 - 2) >> 2UL) << 2UL) + j];
            d = jm1 + 2;
            for (jm1 = b_k; jm1 <= d; jm1++) {
              int32_T tmp_0;
              tmp_0 = ((idxAjj + jm1) - b_k) + 1;
              b_A[tmp_0] += b_A[jm1 - 1] * c;
            }
          }
        }

        ssq = 1.0 / ssq;
        jm1 = (idxAjj - j) + 4;
        b_iy = (((((jm1 - idxAjj) - 1) / 2) << 1UL) + idxAjj) + 2;
        d = b_iy - 2;
        for (int32_T b_k{idxAjj + 2}; b_k <= d; b_k += 2) {
          __m128d tmp;
          tmp = _mm_loadu_pd(&b_A[b_k - 1]);
          (void)_mm_storeu_pd(&b_A[b_k - 1], _mm_mul_pd(tmp, _mm_set1_pd(ssq)));
        }

        for (int32_T b_k{b_iy}; b_k <= jm1; b_k++) {
          b_A[b_k - 1] *= ssq;
        }
      }

      j++;
    } else {
      b_A[idxAjj] = ssq;
      info = j + 1;
      exitg1 = true;
    }
  }

  return info;
}

// Function for MATLAB Function: '<S111>/optimizer'
real_T SupervisoryController::minimum(const real_T x[4])
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
      if (ex > tmp) {
        ex = tmp;
      }
    }
  }

  return ex;
}

// Function for MATLAB Function: '<S111>/optimizer'
void SupervisoryController::trisolve(const real_T b_A[16], real_T B_1[16])
{
  for (int32_T j{0}; j < 4; j++) {
    int32_T jBcol;
    jBcol = j << 2UL;
    for (int32_T b_k{0}; b_k < 4; b_k++) {
      real_T tmp_0;
      int32_T kAcol;
      int32_T tmp;
      kAcol = b_k << 2UL;
      tmp = b_k + jBcol;
      tmp_0 = B_1[tmp];
      if (tmp_0 != 0.0) {
        B_1[tmp] = tmp_0 / b_A[b_k + kAcol];
        for (int32_T i{b_k + 2}; i < 5; i++) {
          int32_T tmp_1;
          tmp_1 = (i + jBcol) - 1;
          B_1[tmp_1] -= b_A[(i + kAcol) - 1] * B_1[tmp];
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S251>/optimizer'
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

// Function for MATLAB Function: '<S251>/optimizer'
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

// Function for MATLAB Function: '<S251>/optimizer'
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

// Function for MATLAB Function: '<S181>/optimizer'
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

// Function for MATLAB Function: '<S181>/optimizer'
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

// Function for MATLAB Function: '<S181>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S111>/optimizer'
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

// Function for MATLAB Function: '<S115>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv_c(const real_T A[21], const real_T B_2[9],
  real_T Y[21])
{
  real_T b_A[9];
  real_T a21;
  real_T maxval;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  int32_T rtemp;
  (void)std::memcpy(&b_A[0], &B_2[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = std::abs(B_2[0]);
  a21 = std::abs(B_2[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (std::abs(B_2[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  b_A[r2] = B_2[r2] / B_2[r1];
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
  for (rtemp = 0; rtemp < 7; rtemp++) {
    int32_T Y_tmp;
    int32_T Y_tmp_0;
    int32_T Y_tmp_1;
    Y_tmp = 7 * r1 + rtemp;
    Y[Y_tmp] = A[rtemp] / b_A[r1];
    Y_tmp_0 = 7 * r2 + rtemp;
    Y[Y_tmp_0] = A[rtemp + 7] - b_A[r1 + 3] * Y[Y_tmp];
    Y_tmp_1 = 7 * r3 + rtemp;
    Y[Y_tmp_1] = A[rtemp + 14] - b_A[r1 + 6] * Y[Y_tmp];
    Y[Y_tmp_0] /= b_A[r2 + 3];
    Y[Y_tmp_1] -= b_A[r2 + 6] * Y[Y_tmp_0];
    Y[Y_tmp_1] /= b_A[r3 + 6];
    Y[Y_tmp_0] -= b_A[r3 + 3] * Y[Y_tmp_1];
    Y[Y_tmp] -= Y[Y_tmp_1] * b_A[r3];
    Y[Y_tmp] -= Y[Y_tmp_0] * b_A[r2];
  }
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
  static const real_T b_Mx_1[2060]{ -0.0, -1.0001226651251065,
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
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0125, 0.0125, 0.0125, -0.0, -0.0,
    -0.025, 0.025, 0.025, -0.0, -0.0, -0.037500000000000006,
    0.037500000000000006, 0.037500000000000006, -0.0, -0.0, -0.05, 0.05, 0.05,
    -0.0, -0.0, -0.0625, 0.0625, 0.0625, -0.0, -0.0, -0.075, 0.075, 0.075, -0.0,
    -0.0, -0.0875, 0.0875, 0.0875, -0.0, -0.0, -0.099999999999999992,
    0.099999999999999992, 0.099999999999999992, -0.0, -0.0, -0.11249999999999999,
    0.11249999999999999, 0.11249999999999999, -0.0, -0.0, -0.12499999999999999,
    0.12499999999999999, 0.12499999999999999, -0.0, -0.0, -0.13749999999999998,
    0.13749999999999998, 0.13749999999999998, -0.0, -0.0, -0.15, 0.15, 0.15,
    -0.0, -0.0, -0.1625, 0.1625, 0.1625, -0.0, -0.0, -0.17500000000000002,
    0.17500000000000002, 0.17500000000000002, -0.0, -0.0, -0.18750000000000003,
    0.18750000000000003, 0.18750000000000003, -0.0, -0.0, -0.20000000000000004,
    0.20000000000000004, 0.20000000000000004, -0.0, -0.0, -0.21250000000000005,
    0.21250000000000005, 0.21250000000000005, -0.0, -0.0, -0.22500000000000006,
    0.22500000000000006, 0.22500000000000006, -0.0, -0.0, -0.23750000000000007,
    0.23750000000000007, 0.23750000000000007, -0.0, -0.0, -0.25000000000000006,
    0.25000000000000006, 0.25000000000000006, -0.0, -0.0, 0.0125, -0.0125,
    -0.0125, 0.0, 0.0, 0.025, -0.025, -0.025, 0.0, 0.0, 0.037500000000000006,
    -0.037500000000000006, -0.037500000000000006, 0.0, 0.0, 0.05, -0.05, -0.05,
    0.0, 0.0, 0.0625, -0.0625, -0.0625, 0.0, 0.0, 0.075, -0.075, -0.075, 0.0,
    0.0, 0.0875, -0.0875, -0.0875, 0.0, 0.0, 0.099999999999999992,
    -0.099999999999999992, -0.099999999999999992, 0.0, 0.0, 0.11249999999999999,
    -0.11249999999999999, -0.11249999999999999, 0.0, 0.0, 0.12499999999999999,
    -0.12499999999999999, -0.12499999999999999, 0.0, 0.0, 0.13749999999999998,
    -0.13749999999999998, -0.13749999999999998, 0.0, 0.0, 0.15, -0.15, -0.15,
    0.0, 0.0, 0.1625, -0.1625, -0.1625, 0.0, 0.0, 0.17500000000000002,
    -0.17500000000000002, -0.17500000000000002, 0.0, 0.0, 0.18750000000000003,
    -0.18750000000000003, -0.18750000000000003, 0.0, 0.0, 0.20000000000000004,
    -0.20000000000000004, -0.20000000000000004, 0.0, 0.0, 0.21250000000000005,
    -0.21250000000000005, -0.21250000000000005, 0.0, 0.0, 0.22500000000000006,
    -0.22500000000000006, -0.22500000000000006, 0.0, 0.0, 0.23750000000000007,
    -0.23750000000000007, -0.23750000000000007, 0.0, 0.0, 0.25000000000000006,
    -0.25000000000000006, -0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5,
    0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5,
    -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.5, -0.5,
    -0.5, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0125, -0.0125, 0.0125, -0.0, -0.0, 0.025, -0.025, 0.025, -0.0, -0.0,
    0.037500000000000006, -0.037500000000000006, 0.037500000000000006, -0.0,
    -0.0, 0.05, -0.05, 0.05, -0.0, -0.0, 0.0625, -0.0625, 0.0625, -0.0, -0.0,
    0.075, -0.075, 0.075, -0.0, -0.0, 0.0875, -0.0875, 0.0875, -0.0, -0.0,
    0.099999999999999992, -0.099999999999999992, 0.099999999999999992, -0.0,
    -0.0, 0.11249999999999999, -0.11249999999999999, 0.11249999999999999, -0.0,
    -0.0, 0.12499999999999999, -0.12499999999999999, 0.12499999999999999, -0.0,
    -0.0, 0.13749999999999998, -0.13749999999999998, 0.13749999999999998, -0.0,
    -0.0, 0.15, -0.15, 0.15, -0.0, -0.0, 0.1625, -0.1625, 0.1625, -0.0, -0.0,
    0.17500000000000002, -0.17500000000000002, 0.17500000000000002, -0.0, -0.0,
    0.18750000000000003, -0.18750000000000003, 0.18750000000000003, -0.0, -0.0,
    0.20000000000000004, -0.20000000000000004, 0.20000000000000004, -0.0, -0.0,
    0.21250000000000005, -0.21250000000000005, 0.21250000000000005, -0.0, -0.0,
    0.22500000000000006, -0.22500000000000006, 0.22500000000000006, -0.0, -0.0,
    0.23750000000000007, -0.23750000000000007, 0.23750000000000007, -0.0, -0.0,
    0.25000000000000006, -0.25000000000000006, 0.25000000000000006, -0.0, -0.0,
    -0.0125, 0.0125, -0.0125, 0.0, 0.0, -0.025, 0.025, -0.025, 0.0, 0.0,
    -0.037500000000000006, 0.037500000000000006, -0.037500000000000006, 0.0, 0.0,
    -0.05, 0.05, -0.05, 0.0, 0.0, -0.0625, 0.0625, -0.0625, 0.0, 0.0, -0.075,
    0.075, -0.075, 0.0, 0.0, -0.0875, 0.0875, -0.0875, 0.0, 0.0,
    -0.099999999999999992, 0.099999999999999992, -0.099999999999999992, 0.0, 0.0,
    -0.11249999999999999, 0.11249999999999999, -0.11249999999999999, 0.0, 0.0,
    -0.12499999999999999, 0.12499999999999999, -0.12499999999999999, 0.0, 0.0,
    -0.13749999999999998, 0.13749999999999998, -0.13749999999999998, 0.0, 0.0,
    -0.15, 0.15, -0.15, 0.0, 0.0, -0.1625, 0.1625, -0.1625, 0.0, 0.0,
    -0.17500000000000002, 0.17500000000000002, -0.17500000000000002, 0.0, 0.0,
    -0.18750000000000003, 0.18750000000000003, -0.18750000000000003, 0.0, 0.0,
    -0.20000000000000004, 0.20000000000000004, -0.20000000000000004, 0.0, 0.0,
    -0.21250000000000005, 0.21250000000000005, -0.21250000000000005, 0.0, 0.0,
    -0.22500000000000006, 0.22500000000000006, -0.22500000000000006, 0.0, 0.0,
    -0.23750000000000007, 0.23750000000000007, -0.23750000000000007, 0.0, 0.0,
    -0.25000000000000006, 0.25000000000000006, -0.25000000000000006, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5,
    -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5,
    -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5,
    -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0,
    0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5,
    0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, 0.0125, -0.0125, -0.0, -0.0,
    0.025, 0.025, -0.025, -0.0, -0.0, 0.037500000000000006, 0.037500000000000006,
    -0.037500000000000006, -0.0, -0.0, 0.05, 0.05, -0.05, -0.0, -0.0, 0.0625,
    0.0625, -0.0625, -0.0, -0.0, 0.075, 0.075, -0.075, -0.0, -0.0, 0.0875,
    0.0875, -0.0875, -0.0, -0.0, 0.099999999999999992, 0.099999999999999992,
    -0.099999999999999992, -0.0, -0.0, 0.11249999999999999, 0.11249999999999999,
    -0.11249999999999999, -0.0, -0.0, 0.12499999999999999, 0.12499999999999999,
    -0.12499999999999999, -0.0, -0.0, 0.13749999999999998, 0.13749999999999998,
    -0.13749999999999998, -0.0, -0.0, 0.15, 0.15, -0.15, -0.0, -0.0, 0.1625,
    0.1625, -0.1625, -0.0, -0.0, 0.17500000000000002, 0.17500000000000002,
    -0.17500000000000002, -0.0, -0.0, 0.18750000000000003, 0.18750000000000003,
    -0.18750000000000003, -0.0, -0.0, 0.20000000000000004, 0.20000000000000004,
    -0.20000000000000004, -0.0, -0.0, 0.21250000000000005, 0.21250000000000005,
    -0.21250000000000005, -0.0, -0.0, 0.22500000000000006, 0.22500000000000006,
    -0.22500000000000006, -0.0, -0.0, 0.23750000000000007, 0.23750000000000007,
    -0.23750000000000007, -0.0, -0.0, 0.25000000000000006, 0.25000000000000006,
    -0.25000000000000006, -0.0, -0.0, -0.0125, -0.0125, 0.0125, 0.0, 0.0, -0.025,
    -0.025, 0.025, 0.0, 0.0, -0.037500000000000006, -0.037500000000000006,
    0.037500000000000006, 0.0, 0.0, -0.05, -0.05, 0.05, 0.0, 0.0, -0.0625,
    -0.0625, 0.0625, 0.0, 0.0, -0.075, -0.075, 0.075, 0.0, 0.0, -0.0875, -0.0875,
    0.0875, 0.0, 0.0, -0.099999999999999992, -0.099999999999999992,
    0.099999999999999992, 0.0, 0.0, -0.11249999999999999, -0.11249999999999999,
    0.11249999999999999, 0.0, 0.0, -0.12499999999999999, -0.12499999999999999,
    0.12499999999999999, 0.0, 0.0, -0.13749999999999998, -0.13749999999999998,
    0.13749999999999998, 0.0, 0.0, -0.15, -0.15, 0.15, 0.0, 0.0, -0.1625,
    -0.1625, 0.1625, 0.0, 0.0, -0.17500000000000002, -0.17500000000000002,
    0.17500000000000002, 0.0, 0.0, -0.18750000000000003, -0.18750000000000003,
    0.18750000000000003, 0.0, 0.0, -0.20000000000000004, -0.20000000000000004,
    0.20000000000000004, 0.0, 0.0, -0.21250000000000005, -0.21250000000000005,
    0.21250000000000005, 0.0, 0.0, -0.22500000000000006, -0.22500000000000006,
    0.22500000000000006, 0.0, 0.0, -0.23750000000000007, -0.23750000000000007,
    0.23750000000000007, 0.0, 0.0, -0.25000000000000006, -0.25000000000000006,
    0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0,
    -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0,
    -0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5,
    0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5,
    0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5,
    -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5,
    0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, -0.5, -0.5,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mx_0[1328]{ -0.9977021834978721, -0.0, -0.0, -0.025,
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
    0.0, 0.0, 0.0, 0.0, -0.0125, 0.0125, 0.0125, -0.0, -0.025, 0.025, 0.025,
    -0.0, -0.037500000000000006, 0.037500000000000006, 0.037500000000000006,
    -0.0, -0.05, 0.05, 0.05, -0.0, -0.0625, 0.0625, 0.0625, -0.0, -0.075, 0.075,
    0.075, -0.0, -0.0875, 0.0875, 0.0875, -0.0, -0.099999999999999992,
    0.099999999999999992, 0.099999999999999992, -0.0, -0.11249999999999999,
    0.11249999999999999, 0.11249999999999999, -0.0, -0.12499999999999999,
    0.12499999999999999, 0.12499999999999999, -0.0, -0.13749999999999998,
    0.13749999999999998, 0.13749999999999998, -0.0, -0.15, 0.15, 0.15, -0.0,
    -0.1625, 0.1625, 0.1625, -0.0, -0.17500000000000002, 0.17500000000000002,
    0.17500000000000002, -0.0, -0.18750000000000003, 0.18750000000000003,
    0.18750000000000003, -0.0, -0.20000000000000004, 0.20000000000000004,
    0.20000000000000004, -0.0, -0.21250000000000005, 0.21250000000000005,
    0.21250000000000005, -0.0, -0.22500000000000006, 0.22500000000000006,
    0.22500000000000006, -0.0, -0.23750000000000007, 0.23750000000000007,
    0.23750000000000007, -0.0, -0.25000000000000006, 0.25000000000000006,
    0.25000000000000006, -0.0, 0.0125, -0.0125, -0.0125, 0.0, 0.025, -0.025,
    -0.025, 0.0, 0.037500000000000006, -0.037500000000000006,
    -0.037500000000000006, 0.0, 0.05, -0.05, -0.05, 0.0, 0.0625, -0.0625,
    -0.0625, 0.0, 0.075, -0.075, -0.075, 0.0, 0.0875, -0.0875, -0.0875, 0.0,
    0.099999999999999992, -0.099999999999999992, -0.099999999999999992, 0.0,
    0.11249999999999999, -0.11249999999999999, -0.11249999999999999, 0.0,
    0.12499999999999999, -0.12499999999999999, -0.12499999999999999, 0.0,
    0.13749999999999998, -0.13749999999999998, -0.13749999999999998, 0.0, 0.15,
    -0.15, -0.15, 0.0, 0.1625, -0.1625, -0.1625, 0.0, 0.17500000000000002,
    -0.17500000000000002, -0.17500000000000002, 0.0, 0.18750000000000003,
    -0.18750000000000003, -0.18750000000000003, 0.0, 0.20000000000000004,
    -0.20000000000000004, -0.20000000000000004, 0.0, 0.21250000000000005,
    -0.21250000000000005, -0.21250000000000005, 0.0, 0.22500000000000006,
    -0.22500000000000006, -0.22500000000000006, 0.0, 0.23750000000000007,
    -0.23750000000000007, -0.23750000000000007, 0.0, 0.25000000000000006,
    -0.25000000000000006, -0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, -0.5, 0.5, 0.5, -0.0, 0.5, -0.5, -0.5,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125,
    -0.0125, 0.0125, -0.0, 0.025, -0.025, 0.025, -0.0, 0.037500000000000006,
    -0.037500000000000006, 0.037500000000000006, -0.0, 0.05, -0.05, 0.05, -0.0,
    0.0625, -0.0625, 0.0625, -0.0, 0.075, -0.075, 0.075, -0.0, 0.0875, -0.0875,
    0.0875, -0.0, 0.099999999999999992, -0.099999999999999992,
    0.099999999999999992, -0.0, 0.11249999999999999, -0.11249999999999999,
    0.11249999999999999, -0.0, 0.12499999999999999, -0.12499999999999999,
    0.12499999999999999, -0.0, 0.13749999999999998, -0.13749999999999998,
    0.13749999999999998, -0.0, 0.15, -0.15, 0.15, -0.0, 0.1625, -0.1625, 0.1625,
    -0.0, 0.17500000000000002, -0.17500000000000002, 0.17500000000000002, -0.0,
    0.18750000000000003, -0.18750000000000003, 0.18750000000000003, -0.0,
    0.20000000000000004, -0.20000000000000004, 0.20000000000000004, -0.0,
    0.21250000000000005, -0.21250000000000005, 0.21250000000000005, -0.0,
    0.22500000000000006, -0.22500000000000006, 0.22500000000000006, -0.0,
    0.23750000000000007, -0.23750000000000007, 0.23750000000000007, -0.0,
    0.25000000000000006, -0.25000000000000006, 0.25000000000000006, -0.0,
    -0.0125, 0.0125, -0.0125, 0.0, -0.025, 0.025, -0.025, 0.0,
    -0.037500000000000006, 0.037500000000000006, -0.037500000000000006, 0.0,
    -0.05, 0.05, -0.05, 0.0, -0.0625, 0.0625, -0.0625, 0.0, -0.075, 0.075,
    -0.075, 0.0, -0.0875, 0.0875, -0.0875, 0.0, -0.099999999999999992,
    0.099999999999999992, -0.099999999999999992, 0.0, -0.11249999999999999,
    0.11249999999999999, -0.11249999999999999, 0.0, -0.12499999999999999,
    0.12499999999999999, -0.12499999999999999, 0.0, -0.13749999999999998,
    0.13749999999999998, -0.13749999999999998, 0.0, -0.15, 0.15, -0.15, 0.0,
    -0.1625, 0.1625, -0.1625, 0.0, -0.17500000000000002, 0.17500000000000002,
    -0.17500000000000002, 0.0, -0.18750000000000003, 0.18750000000000003,
    -0.18750000000000003, 0.0, -0.20000000000000004, 0.20000000000000004,
    -0.20000000000000004, 0.0, -0.21250000000000005, 0.21250000000000005,
    -0.21250000000000005, 0.0, -0.22500000000000006, 0.22500000000000006,
    -0.22500000000000006, 0.0, -0.23750000000000007, 0.23750000000000007,
    -0.23750000000000007, 0.0, -0.25000000000000006, 0.25000000000000006,
    -0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.5,
    -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5,
    -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5,
    -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5,
    -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5,
    -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, 0.5,
    -0.5, 0.5, -0.0, 0.5, -0.5, 0.5, -0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5,
    0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5,
    0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5,
    0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5,
    0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5,
    0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5, 0.5, -0.5, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, 0.0125, -0.0125, -0.0,
    0.025, 0.025, -0.025, -0.0, 0.037500000000000006, 0.037500000000000006,
    -0.037500000000000006, -0.0, 0.05, 0.05, -0.05, -0.0, 0.0625, 0.0625,
    -0.0625, -0.0, 0.075, 0.075, -0.075, -0.0, 0.0875, 0.0875, -0.0875, -0.0,
    0.099999999999999992, 0.099999999999999992, -0.099999999999999992, -0.0,
    0.11249999999999999, 0.11249999999999999, -0.11249999999999999, -0.0,
    0.12499999999999999, 0.12499999999999999, -0.12499999999999999, -0.0,
    0.13749999999999998, 0.13749999999999998, -0.13749999999999998, -0.0, 0.15,
    0.15, -0.15, -0.0, 0.1625, 0.1625, -0.1625, -0.0, 0.17500000000000002,
    0.17500000000000002, -0.17500000000000002, -0.0, 0.18750000000000003,
    0.18750000000000003, -0.18750000000000003, -0.0, 0.20000000000000004,
    0.20000000000000004, -0.20000000000000004, -0.0, 0.21250000000000005,
    0.21250000000000005, -0.21250000000000005, -0.0, 0.22500000000000006,
    0.22500000000000006, -0.22500000000000006, -0.0, 0.23750000000000007,
    0.23750000000000007, -0.23750000000000007, -0.0, 0.25000000000000006,
    0.25000000000000006, -0.25000000000000006, -0.0, -0.0125, -0.0125, 0.0125,
    0.0, -0.025, -0.025, 0.025, 0.0, -0.037500000000000006,
    -0.037500000000000006, 0.037500000000000006, 0.0, -0.05, -0.05, 0.05, 0.0,
    -0.0625, -0.0625, 0.0625, 0.0, -0.075, -0.075, 0.075, 0.0, -0.0875, -0.0875,
    0.0875, 0.0, -0.099999999999999992, -0.099999999999999992,
    0.099999999999999992, 0.0, -0.11249999999999999, -0.11249999999999999,
    0.11249999999999999, 0.0, -0.12499999999999999, -0.12499999999999999,
    0.12499999999999999, 0.0, -0.13749999999999998, -0.13749999999999998,
    0.13749999999999998, 0.0, -0.15, -0.15, 0.15, 0.0, -0.1625, -0.1625, 0.1625,
    0.0, -0.17500000000000002, -0.17500000000000002, 0.17500000000000002, 0.0,
    -0.18750000000000003, -0.18750000000000003, 0.18750000000000003, 0.0,
    -0.20000000000000004, -0.20000000000000004, 0.20000000000000004, 0.0,
    -0.21250000000000005, -0.21250000000000005, 0.21250000000000005, 0.0,
    -0.22500000000000006, -0.22500000000000006, 0.22500000000000006, 0.0,
    -0.23750000000000007, -0.23750000000000007, 0.23750000000000007, 0.0,
    -0.25000000000000006, -0.25000000000000006, 0.25000000000000006, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5,
    0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5,
    -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5,
    0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5,
    -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5,
    0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5,
    -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5,
    0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5,
    -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5,
    0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5,
    -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T b_Mx[1008]{ -1.00107782, -0.0, -0.0, -1.0021568016959521,
    -0.0, -0.0, -1.003236946339956, -0.0, -0.0, -1.00431825518546, -0.0, -0.0,
    -1.0054007294872638, -0.0, -0.0, -1.0064843705015196, -0.0, -0.0,
    -1.0075691794857335, -0.0, -0.0, -1.0086551576987668, -0.0, -0.0,
    -1.0097423064008375, -0.0, -0.0, -1.0108306268535223, -0.0, -0.0,
    -1.0119201203197574, -0.0, -0.0, -1.0130107880638404, -0.0, -0.0,
    -1.0141026313514312, -0.0, -0.0, -1.0151956514495544, -0.0, -0.0,
    -1.0162898496265995, -0.0, -0.0, -1.017385227152324, -0.0, -0.0,
    -1.0184817852978532, -0.0, -0.0, -1.0195795253356827, -0.0, -0.0,
    -1.02067844853968, -0.0, -0.0, -1.021778556185085, -0.0, -0.0, 1.00107782,
    0.0, 0.0, 1.0021568016959521, 0.0, 0.0, 1.003236946339956, 0.0, 0.0,
    1.00431825518546, 0.0, 0.0, 1.0054007294872638, 0.0, 0.0, 1.0064843705015196,
    0.0, 0.0, 1.0075691794857335, 0.0, 0.0, 1.0086551576987668, 0.0, 0.0,
    1.0097423064008375, 0.0, 0.0, 1.0108306268535223, 0.0, 0.0,
    1.0119201203197574, 0.0, 0.0, 1.0130107880638404, 0.0, 0.0,
    1.0141026313514312, 0.0, 0.0, 1.0151956514495544, 0.0, 0.0,
    1.0162898496265995, 0.0, 0.0, 1.017385227152324, 0.0, 0.0,
    1.0184817852978532, 0.0, 0.0, 1.0195795253356827, 0.0, 0.0, 1.02067844853968,
    0.0, 0.0, 1.021778556185085, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    -0.99930767, -0.0, -0.0, -0.99861581932082888, -0.0, -0.0,
    -0.99792444763063848, -0.0, -0.0, -0.99723355459781038, -0.0, -0.0,
    -0.99654313989095566, -0.0, -0.0, -0.9958532031789149, -0.0, -0.0,
    -0.995163744130758, -0.0, -0.0, -0.99447476241578392, -0.0, -0.0,
    -0.99378625770352058, -0.0, -0.0, -0.99309822966372463, -0.0, -0.0,
    -0.99241067796638149, -0.0, -0.0, -0.991723602281705, -0.0, -0.0,
    -0.99103700228013725, -0.0, -0.0, -0.99035087763234864, -0.0, -0.0,
    -0.9896652280092374, -0.0, -0.0, -0.98898005308192971, -0.0, -0.0,
    -0.98829535252177947, -0.0, -0.0, -0.987611126000368, -0.0, -0.0,
    -0.98692737318950419, -0.0, -0.0, -0.98624409376122391, -0.0, 0.0,
    0.99930767, 0.0, 0.0, 0.99861581932082888, 0.0, 0.0, 0.99792444763063848,
    0.0, 0.0, 0.99723355459781038, 0.0, 0.0, 0.99654313989095566, 0.0, 0.0,
    0.9958532031789149, 0.0, 0.0, 0.995163744130758, 0.0, 0.0,
    0.99447476241578392, 0.0, 0.0, 0.99378625770352058, 0.0, 0.0,
    0.99309822966372463, 0.0, 0.0, 0.99241067796638149, 0.0, 0.0,
    0.991723602281705, 0.0, 0.0, 0.99103700228013725, 0.0, 0.0,
    0.99035087763234864, 0.0, 0.0, 0.9896652280092374, 0.0, 0.0,
    0.98898005308192971, 0.0, 0.0, 0.98829535252177947, 0.0, 0.0,
    0.987611126000368, 0.0, 0.0, 0.98692737318950419, 0.0, 0.0,
    0.98624409376122391, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0125, 0.0125,
    0.0125, -0.025, 0.025, 0.025, -0.037500000000000006, 0.037500000000000006,
    0.037500000000000006, -0.05, 0.05, 0.05, -0.0625, 0.0625, 0.0625, -0.075,
    0.075, 0.075, -0.0875, 0.0875, 0.0875, -0.099999999999999992,
    0.099999999999999992, 0.099999999999999992, -0.11249999999999999,
    0.11249999999999999, 0.11249999999999999, -0.12499999999999999,
    0.12499999999999999, 0.12499999999999999, -0.13749999999999998,
    0.13749999999999998, 0.13749999999999998, -0.15, 0.15, 0.15, -0.1625, 0.1625,
    0.1625, -0.17500000000000002, 0.17500000000000002, 0.17500000000000002,
    -0.18750000000000003, 0.18750000000000003, 0.18750000000000003,
    -0.20000000000000004, 0.20000000000000004, 0.20000000000000004,
    -0.21250000000000005, 0.21250000000000005, 0.21250000000000005,
    -0.22500000000000006, 0.22500000000000006, 0.22500000000000006,
    -0.23750000000000007, 0.23750000000000007, 0.23750000000000007,
    -0.25000000000000006, 0.25000000000000006, 0.25000000000000006, 0.0125,
    -0.0125, -0.0125, 0.025, -0.025, -0.025, 0.037500000000000006,
    -0.037500000000000006, -0.037500000000000006, 0.05, -0.05, -0.05, 0.0625,
    -0.0625, -0.0625, 0.075, -0.075, -0.075, 0.0875, -0.0875, -0.0875,
    0.099999999999999992, -0.099999999999999992, -0.099999999999999992,
    0.11249999999999999, -0.11249999999999999, -0.11249999999999999,
    0.12499999999999999, -0.12499999999999999, -0.12499999999999999,
    0.13749999999999998, -0.13749999999999998, -0.13749999999999998, 0.15, -0.15,
    -0.15, 0.1625, -0.1625, -0.1625, 0.17500000000000002, -0.17500000000000002,
    -0.17500000000000002, 0.18750000000000003, -0.18750000000000003,
    -0.18750000000000003, 0.20000000000000004, -0.20000000000000004,
    -0.20000000000000004, 0.21250000000000005, -0.21250000000000005,
    -0.21250000000000005, 0.22500000000000006, -0.22500000000000006,
    -0.22500000000000006, 0.23750000000000007, -0.23750000000000007,
    -0.23750000000000007, 0.25000000000000006, -0.25000000000000006,
    -0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.5, -0.5,
    0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
    -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5,
    0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
    0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
    0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5,
    -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
    -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
    0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5,
    -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0125, -0.0125, 0.0125, 0.025, -0.025, 0.025, 0.037500000000000006,
    -0.037500000000000006, 0.037500000000000006, 0.05, -0.05, 0.05, 0.0625,
    -0.0625, 0.0625, 0.075, -0.075, 0.075, 0.0875, -0.0875, 0.0875,
    0.099999999999999992, -0.099999999999999992, 0.099999999999999992,
    0.11249999999999999, -0.11249999999999999, 0.11249999999999999,
    0.12499999999999999, -0.12499999999999999, 0.12499999999999999,
    0.13749999999999998, -0.13749999999999998, 0.13749999999999998, 0.15, -0.15,
    0.15, 0.1625, -0.1625, 0.1625, 0.17500000000000002, -0.17500000000000002,
    0.17500000000000002, 0.18750000000000003, -0.18750000000000003,
    0.18750000000000003, 0.20000000000000004, -0.20000000000000004,
    0.20000000000000004, 0.21250000000000005, -0.21250000000000005,
    0.21250000000000005, 0.22500000000000006, -0.22500000000000006,
    0.22500000000000006, 0.23750000000000007, -0.23750000000000007,
    0.23750000000000007, 0.25000000000000006, -0.25000000000000006,
    0.25000000000000006, -0.0125, 0.0125, -0.0125, -0.025, 0.025, -0.025,
    -0.037500000000000006, 0.037500000000000006, -0.037500000000000006, -0.05,
    0.05, -0.05, -0.0625, 0.0625, -0.0625, -0.075, 0.075, -0.075, -0.0875,
    0.0875, -0.0875, -0.099999999999999992, 0.099999999999999992,
    -0.099999999999999992, -0.11249999999999999, 0.11249999999999999,
    -0.11249999999999999, -0.12499999999999999, 0.12499999999999999,
    -0.12499999999999999, -0.13749999999999998, 0.13749999999999998,
    -0.13749999999999998, -0.15, 0.15, -0.15, -0.1625, 0.1625, -0.1625,
    -0.17500000000000002, 0.17500000000000002, -0.17500000000000002,
    -0.18750000000000003, 0.18750000000000003, -0.18750000000000003,
    -0.20000000000000004, 0.20000000000000004, -0.20000000000000004,
    -0.21250000000000005, 0.21250000000000005, -0.21250000000000005,
    -0.22500000000000006, 0.22500000000000006, -0.22500000000000006,
    -0.23750000000000007, 0.23750000000000007, -0.23750000000000007,
    -0.25000000000000006, 0.25000000000000006, -0.25000000000000006, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
    -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5,
    0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
    0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
    -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
    -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
    0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5,
    -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
    -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, 0.0125, -0.0125, 0.025,
    0.025, -0.025, 0.037500000000000006, 0.037500000000000006,
    -0.037500000000000006, 0.05, 0.05, -0.05, 0.0625, 0.0625, -0.0625, 0.075,
    0.075, -0.075, 0.0875, 0.0875, -0.0875, 0.099999999999999992,
    0.099999999999999992, -0.099999999999999992, 0.11249999999999999,
    0.11249999999999999, -0.11249999999999999, 0.12499999999999999,
    0.12499999999999999, -0.12499999999999999, 0.13749999999999998,
    0.13749999999999998, -0.13749999999999998, 0.15, 0.15, -0.15, 0.1625, 0.1625,
    -0.1625, 0.17500000000000002, 0.17500000000000002, -0.17500000000000002,
    0.18750000000000003, 0.18750000000000003, -0.18750000000000003,
    0.20000000000000004, 0.20000000000000004, -0.20000000000000004,
    0.21250000000000005, 0.21250000000000005, -0.21250000000000005,
    0.22500000000000006, 0.22500000000000006, -0.22500000000000006,
    0.23750000000000007, 0.23750000000000007, -0.23750000000000007,
    0.25000000000000006, 0.25000000000000006, -0.25000000000000006, -0.0125,
    -0.0125, 0.0125, -0.025, -0.025, 0.025, -0.037500000000000006,
    -0.037500000000000006, 0.037500000000000006, -0.05, -0.05, 0.05, -0.0625,
    -0.0625, 0.0625, -0.075, -0.075, 0.075, -0.0875, -0.0875, 0.0875,
    -0.099999999999999992, -0.099999999999999992, 0.099999999999999992,
    -0.11249999999999999, -0.11249999999999999, 0.11249999999999999,
    -0.12499999999999999, -0.12499999999999999, 0.12499999999999999,
    -0.13749999999999998, -0.13749999999999998, 0.13749999999999998, -0.15,
    -0.15, 0.15, -0.1625, -0.1625, 0.1625, -0.17500000000000002,
    -0.17500000000000002, 0.17500000000000002, -0.18750000000000003,
    -0.18750000000000003, 0.18750000000000003, -0.20000000000000004,
    -0.20000000000000004, 0.20000000000000004, -0.21250000000000005,
    -0.21250000000000005, 0.21250000000000005, -0.22500000000000006,
    -0.22500000000000006, 0.22500000000000006, -0.23750000000000007,
    -0.23750000000000007, 0.23750000000000007, -0.25000000000000006,
    -0.25000000000000006, 0.25000000000000006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5,
    0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
    0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
    -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5,
    0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
    0.5, 0.5, -0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
    -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
    0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5,
    -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5,
    -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T c_1[1000]{ 0.0, 0.0, 0.0, 0.0, 0.0125, 0.5, -0.0125, -0.5,
    -0.0125, -0.5, 1.0001226651251065, -6.4394278387945068E-5, 0.0, 0.0, -0.0125,
    -0.5, 0.0125, 0.5, -0.0125, -0.5, -0.0016804698402870472, 0.999403766481181,
    0.0, 0.0, -0.0125, -0.5, -0.0125, -0.5, 0.0125, 0.5, 0.025, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.5, -0.025, -0.5, -0.025, -0.5,
    1.0002454535095886, -0.00012875806168090968, 0.0, 0.0, -0.025, -0.5, 0.025,
    0.5, -0.025, -0.5, -0.0033601438631711473, 0.99880799666941367, 0.0, 0.0,
    -0.025, -0.5, -0.025, -0.5, 0.025, 0.5, 0.050003066628127667,
    -1.6098569596986267E-6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -4.201174600717618E-5, 0.049985094162029529, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.037500000000000006, 0.5,
    -0.037500000000000006, -0.5, -0.037500000000000006, -0.5, 1.0003683651173199,
    -0.00019309137599829017, 0.0, 0.0, -0.037500000000000006, -0.5,
    0.037500000000000006, 0.5, -0.037500000000000006, -0.5,
    -0.0050390227502789741, 0.99821269023697423, 0.0, 0.0, -0.037500000000000006,
    -0.5, -0.037500000000000006, -0.5, 0.037500000000000006, 0.5,
    0.075009202965867383, -4.8288085017213693E-6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.00012601534258645488, 0.074955294078764875, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.5, -0.05, -0.5, -0.05, -0.5,
    1.0004913999122134, -0.00025739424744158326, 0.0, 0.0, -0.05, -0.5, 0.05,
    0.5, -0.05, -0.5, -0.0067171071827700837, 0.99761784685637822, 0.0, 0.0,
    -0.05, -0.5, -0.05, -0.5, 0.05, 0.5, 0.10001841209380039,
    -9.6560929016786243E-6, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.00025199091134342927, 0.099910611334689231, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0625, 0.5, -0.0625, -0.5, -0.0625, -0.5,
    1.0006145578582215, -0.00032166670209439935, 0.0, 0.0, -0.0625, -0.5, 0.0625,
    0.5, -0.0625, -0.5, -0.0083943978413372582, 0.99702346620038029, 0.0, 0.0,
    -0.0625, -0.5, -0.0625, -0.5, 0.0625, 0.5, 0.12503069709160572,
    -1.6090949087718206E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.00041991859091268134, 0.12485105750609868, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.075, 0.5, -0.075, -0.5, -0.075, -0.5,
    1.0007378389193362, -0.00038590876602247561, 0.0, 0.0, -0.075, -0.5, 0.075,
    0.5, -0.075, -0.5, -0.010070895406206851, 0.996429547941974, 0.0, 0.0,
    -0.075, -0.5, -0.075, -0.5, 0.075, 0.5, 0.15004606103806126,
    -2.4132616640078191E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.00062977853694611279, 0.1497766441611082, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0875, 0.5, -0.0875, -0.5, -0.0875, -0.5,
    1.0008612430595885, -0.00045012046527368919, 0.0, 0.0, -0.0875, -0.5, 0.0875,
    0.5, -0.0875, -0.5, -0.011746600557139133, 0.99583609175439181, 0.0, 0.0,
    -0.0875, -0.5, -0.0875, -0.5, 0.0875, 0.5, 0.17506450701104467,
    -3.3780335790640084E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.000881550922101284, 0.17468738285965754, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.099999999999999992, 0.5,
    -0.099999999999999992, -0.5, -0.099999999999999992, -0.5, 1.000984770243049,
    -0.00051430182587807046, 0.0, 0.0, -0.099999999999999992, -0.5,
    0.099999999999999992, 0.5, -0.099999999999999992, -0.5,
    -0.013421513973428631, 0.9952430973111045, 0.0, 0.0, -0.099999999999999992,
    -0.5, -0.099999999999999992, -0.5, 0.099999999999999992, 0.5,
    0.20008603808753436, -4.5033347422482317E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0011752159360297622, 0.19958328515351734, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.11249999999999999, 0.5,
    -0.11249999999999999, -0.5, -0.11249999999999999, -0.5, 1.0011084204338279,
    -0.00057845287384781629, 0.0, 0.0, -0.11249999999999999, -0.5,
    0.11249999999999999, 0.5, -0.11249999999999999, -0.5, -0.015095636333904478,
    0.99465056428582155, 0.0, 0.0, -0.11249999999999999, -0.5,
    -0.11249999999999999, -0.5, 0.11249999999999999, 0.5, 0.22511065734361058,
    -5.7890893069434082E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0015107537853654778, 0.22446436258629496, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.12499999999999999, 0.5, -0.12499999999999999,
    -0.5, -0.12499999999999999, -0.5, 1.0012321935960742,
    -0.00064257363517730291, 0.0, 0.0, -0.12499999999999999, -0.5,
    0.12499999999999999, 0.5, -0.12499999999999999, -0.5, -0.016768968316930755,
    0.99405849235249066, 0.0, 0.0, -0.12499999999999999, -0.5,
    -0.12499999999999999, -0.5, 0.12499999999999999, 0.5, 0.25013836785445626,
    -7.2352214915629483E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0018881446937130896, 0.2493306266934405, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.13749999999999998, 0.5, -0.13749999999999998,
    -0.5, -0.13749999999999998, -0.5, 1.0013560896939764,
    -0.00070666413584309932, 0.0, 0.0, -0.13749999999999998, -0.5,
    0.13749999999999998, 0.5, -0.13749999999999998, -0.5, -0.01844151060040683,
    0.99346688118529747, 0.0, 0.0, -0.13749999999999998, -0.5,
    -0.13749999999999998, -0.5, 0.13749999999999998, 0.5, 0.27516917269435814,
    -8.8416555795062062E-5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0023073689016363582, 0.27418208900225277, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.15, 0.5, -0.15, -0.5, -0.15, -0.5,
    1.0014801086917624, -0.00077072440180398047, 0.0, 0.0, -0.15, -0.5, 0.15,
    0.5, -0.15, -0.5, -0.020113263861767709, 0.99287573045866573, 0.0, 0.0,
    -0.15, -0.5, -0.15, -0.5, 0.15, 0.5, 0.30020307493670756,
    -0.00010608315919113954, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0027684066666465289, 0.29901876103188524, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1625, 0.5, -0.1625, -0.5, -0.1625, -0.5,
    1.0016042505536991, -0.00083475445900094006, 0.0, 0.0, -0.1625, -0.5, 0.1625,
    0.5, -0.1625, -0.5, -0.021784228777984374, 0.9922850398472568, 0.0, 0.0,
    -0.1625, -0.5, -0.1625, -0.5, 0.1625, 0.5, 0.32524007765400159,
    -0.00012535126923623904, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0032712382631907219, 0.3238406542933519, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17500000000000002, 0.5, -0.17500000000000002,
    -0.5, -0.17500000000000002, -0.5, 1.0017285152440929,
    -0.00089875433335720413, 0.0, 0.0, -0.17500000000000002, -0.5,
    0.17500000000000002, 0.5, -0.17500000000000002, -0.5, -0.023454406025564121,
    0.99169480902596963, 0.0, 0.0, -0.17500000000000002, -0.5,
    -0.17500000000000002, -0.5, 0.17500000000000002, 0.5, 0.3502801839178441,
    -0.00014622013071126255, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0038158439826403317, 0.34864778028953336, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.18750000000000003, 0.5, -0.18750000000000003,
    -0.5, -0.18750000000000003, -0.5, 1.0018529027272891,
    -0.00096272405077824379, 0.0, 0.0, -0.18750000000000003, -0.5,
    0.18750000000000003, 0.5, -0.18750000000000003, -0.5, -0.025123796280550909,
    0.99110503766994063, 0.0, 0.0, -0.18750000000000003, -0.5,
    -0.18750000000000003, -0.5, 0.18750000000000003, 0.5, 0.37532339679894644,
    -0.00016868898904519265, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0044022041332794351, 0.37344015051518265, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.20000000000000004, 0.5, -0.20000000000000004,
    -0.5, -0.20000000000000004, -0.5, 1.0019774129676724, -0.0010266636371517885,
    0.0, 0.0, -0.20000000000000004, -0.5, 0.20000000000000004, 0.5,
    -0.20000000000000004, -0.5, -0.026792400218525705, 0.99051572545454325, 0.0,
    0.0, -0.20000000000000004, -0.5, -0.20000000000000004, -0.5,
    0.20000000000000004, 0.5, 0.40036971936712867, -0.00019275709031464874, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0050302990402932082,
    0.39821777645693118, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.21250000000000005, 0.5, -0.21250000000000005, -0.5,
    -0.21250000000000005, -0.5, 1.0021020459296663, -0.0010905731183478392, 0.0,
    0.0, -0.21250000000000005, -0.5, 0.21250000000000005, 0.5,
    -0.21250000000000005, -0.5, -0.028460218514606821, 0.98992687205538832, 0.0,
    0.0, -0.21250000000000005, -0.5, -0.21250000000000005, -0.5,
    0.21250000000000005, 0.5, 0.42541915469132052, -0.00021842368124344344, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.005700109045756351,
    0.42298066959329478, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.22500000000000006, 0.5, -0.22500000000000006, -0.5,
    -0.22500000000000006, -0.5, 1.0022268015777338, -0.0011544525202186815, 0.0,
    0.0, -0.22500000000000006, -0.5, 0.22500000000000006, 0.5,
    -0.22500000000000006, -0.5, -0.030127251843450248, 0.98933847714832324, 0.0,
    0.0, -0.22500000000000006, -0.5, -0.22500000000000006, -0.5,
    0.22500000000000006, 0.5, 0.4504717058395622, -0.00024568800920213941, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.006411614508621522, 0.4477288413946795,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.23750000000000007, 0.5, -0.23750000000000007, -0.5, -0.23750000000000007,
    -0.5, 1.0023516798763767, -0.0012183018685988985, 0.0, 0.0,
    -0.23750000000000007, -0.5, 0.23750000000000007, 0.5, -0.23750000000000007,
    -0.5, -0.031793500879250015, 0.98875054040943233, 0.0, 0.0,
    -0.23750000000000007, -0.5, -0.23750000000000007, -0.5, 0.23750000000000007,
    0.5, 0.47552737587900556, -0.00027454932220760648, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0071647958047077786, 0.47246230332338762, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25000000000000006, 0.5,
    -0.25000000000000006, -0.5, -0.25000000000000006, -0.5, 1.002476680790136,
    -0.0012821211893053839, 0.0, 0.0, -0.25000000000000006, -0.5,
    0.25000000000000006, 0.5, -0.25000000000000006, -0.5, -0.03345896629573851,
    0.98816306151503641, 0.0, 0.0, -0.25000000000000006, -0.5,
    -0.25000000000000006, -0.5, 0.25000000000000006, 0.5, 0.500586167875915,
    -0.00030500686892257894, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0079596333266890289, 0.49718106683362345, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 };

  static const real_T b_Ac_1[824]{ -0.0, 0.0071318388640683669,
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

  static const real_T b_Ac_0[664]{ -0.22599754912864506, -0.0, -0.0, -0.0,
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

  static const real_T c_0[640]{ 0.9977021834978721, 0.0, 0.0125, 0.5, -0.0125,
    -0.5, -0.0125, -0.5, 0.0, 0.0, -0.0125, -0.5, 0.0125, 0.5, -0.0125, -0.5,
    0.0, 0.0, -0.0125, -0.5, -0.0125, -0.5, 0.0125, 0.5, 0.025, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.99540964695642165, 0.0, 0.025, 0.5, -0.025, -0.5,
    -0.025, -0.5, 0.0, 0.0, -0.025, -0.5, 0.025, 0.5, -0.025, -0.5, 0.0, 0.0,
    -0.025, -0.5, -0.025, -0.5, 0.025, 0.5, 0.049942554587446807, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.99312237824326788, 0.0, 0.037500000000000006, 0.5,
    -0.037500000000000006, -0.5, -0.037500000000000006, -0.5, 0.0, 0.0,
    -0.037500000000000006, -0.5, 0.037500000000000006, 0.5,
    -0.037500000000000006, -0.5, 0.0, 0.0, -0.037500000000000006, -0.5,
    -0.037500000000000006, -0.5, 0.037500000000000006, 0.5, 0.074827795761357341,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.990840365253908, 0.0, 0.05, 0.5, -0.05,
    -0.5, -0.05, -0.5, 0.0, 0.0, -0.05, -0.5, 0.05, 0.5, -0.05, -0.5, 0.0, 0.0,
    -0.05, -0.5, -0.05, -0.5, 0.05, 0.5, 0.099655855217439027, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.9885635959116531, 0.0, 0.0625, 0.5, -0.0625, -0.5,
    -0.0625, -0.5, 0.0, 0.0, -0.0625, -0.5, 0.0625, 0.5, -0.0625, -0.5, 0.0, 0.0,
    -0.0625, -0.5, -0.0625, -0.5, 0.0625, 0.5, 0.12442686434878672, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.98629205816756438, 0.0, 0.075, 0.5, -0.075, -0.5,
    -0.075, -0.5, 0.0, 0.0, -0.075, -0.5, 0.075, 0.5, -0.075, -0.5, 0.0, 0.0,
    -0.075, -0.5, -0.075, -0.5, 0.075, 0.5, 0.14914095424657806, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.98402574000038923, 0.0, 0.0875, 0.5, -0.0875, -0.5,
    -0.0875, -0.5, 0.0, 0.0, -0.0875, -0.5, 0.0875, 0.5, -0.0875, -0.5, 0.0, 0.0,
    -0.0875, -0.5, -0.0875, -0.5, 0.0875, 0.5, 0.17379825570076715, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.98176462941649767, 0.0, 0.099999999999999992, 0.5,
    -0.099999999999999992, -0.5, -0.099999999999999992, -0.5, 0.0, 0.0,
    -0.099999999999999992, -0.5, 0.099999999999999992, 0.5,
    -0.099999999999999992, -0.5, 0.0, 0.0, -0.099999999999999992, -0.5,
    -0.099999999999999992, -0.5, 0.099999999999999992, 0.5, 0.19839889920077688,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97950871444981891, 0.0,
    0.11249999999999999, 0.5, -0.11249999999999999, -0.5, -0.11249999999999999,
    -0.5, 0.0, 0.0, -0.11249999999999999, -0.5, 0.11249999999999999, 0.5,
    -0.11249999999999999, -0.5, 0.0, 0.0, -0.11249999999999999, -0.5,
    -0.11249999999999999, -0.5, 0.11249999999999999, 0.5, 0.22294301493618932,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.977257983161778, 0.0,
    0.12499999999999999, 0.5, -0.12499999999999999, -0.5, -0.12499999999999999,
    -0.5, 0.0, 0.0, -0.12499999999999999, -0.5, 0.12499999999999999, 0.5,
    -0.12499999999999999, -0.5, 0.0, 0.0, -0.12499999999999999, -0.5,
    -0.12499999999999999, -0.5, 0.12499999999999999, 0.5, 0.2474307327974348,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97501242364123264, 0.0,
    0.13749999999999998, 0.5, -0.13749999999999998, -0.5, -0.13749999999999998,
    -0.5, 0.0, 0.0, -0.13749999999999998, -0.5, 0.13749999999999998, 0.5,
    -0.13749999999999998, -0.5, 0.0, 0.0, -0.13749999999999998, -0.5,
    -0.13749999999999998, -0.5, 0.13749999999999998, 0.5, 0.27186218237647924,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9727720240044101, 0.0, 0.15, 0.5, -0.15,
    -0.5, -0.15, -0.5, 0.0, 0.0, -0.15, -0.5, 0.15, 0.5, -0.15, -0.5, 0.0, 0.0,
    -0.15, -0.5, -0.15, -0.5, 0.15, 0.5, 0.29623749296751006, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97053677239484437, 0.0, 0.1625, 0.5, -0.1625, -0.5, -0.1625,
    -0.5, 0.0, 0.0, -0.1625, -0.5, 0.1625, 0.5, -0.1625, -0.5, 0.0, 0.0, -0.1625,
    -0.5, -0.1625, -0.5, 0.1625, 0.5, 0.32055679356762035, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.96830665698331353, 0.0, 0.17500000000000002, 0.5,
    -0.17500000000000002, -0.5, -0.17500000000000002, -0.5, 0.0, 0.0,
    -0.17500000000000002, -0.5, 0.17500000000000002, 0.5, -0.17500000000000002,
    -0.5, 0.0, 0.0, -0.17500000000000002, -0.5, -0.17500000000000002, -0.5,
    0.17500000000000002, 0.5, 0.34482021287749148, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.966081665967777, 0.0, 0.18750000000000003, 0.5, -0.18750000000000003,
    -0.5, -0.18750000000000003, -0.5, 0.0, 0.0, -0.18750000000000003, -0.5,
    0.18750000000000003, 0.5, -0.18750000000000003, -0.5, 0.0, 0.0,
    -0.18750000000000003, -0.5, -0.18750000000000003, -0.5, 0.18750000000000003,
    0.5, 0.36902787930207437, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.963861787573313, 0.0, 0.20000000000000004, 0.5, -0.20000000000000004, -0.5,
    -0.20000000000000004, -0.5, 0.0, 0.0, -0.20000000000000004, -0.5,
    0.20000000000000004, 0.5, -0.20000000000000004, -0.5, 0.0, 0.0,
    -0.20000000000000004, -0.5, -0.20000000000000004, -0.5, 0.20000000000000004,
    0.5, 0.3931799209512688, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.96164701005205655, 0.0, 0.21250000000000005, 0.5, -0.21250000000000005,
    -0.5, -0.21250000000000005, -0.5, 0.0, 0.0, -0.21250000000000005, -0.5,
    0.21250000000000005, 0.5, -0.21250000000000005, -0.5, 0.0, 0.0,
    -0.21250000000000005, -0.5, -0.21250000000000005, -0.5, 0.21250000000000005,
    0.5, 0.41727646564060167, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.959437321683137, 0.0, 0.22500000000000006, 0.5, -0.22500000000000006, -0.5,
    -0.22500000000000006, -0.5, 0.0, 0.0, -0.22500000000000006, -0.5,
    0.22500000000000006, 0.5, -0.22500000000000006, -0.5, 0.0, 0.0,
    -0.22500000000000006, -0.5, -0.22500000000000006, -0.5, 0.22500000000000006,
    0.5, 0.44131764089190312, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.95723271077261607, 0.0, 0.23750000000000007, 0.5, -0.23750000000000007,
    -0.5, -0.23750000000000007, -0.5, 0.0, 0.0, -0.23750000000000007, -0.5,
    0.23750000000000007, 0.5, -0.23750000000000007, -0.5, 0.0, 0.0,
    -0.23750000000000007, -0.5, -0.23750000000000007, -0.5, 0.23750000000000007,
    0.5, 0.46530357393398158, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.95503316565342611, 0.0, 0.25000000000000006, 0.5, -0.25000000000000006,
    -0.5, -0.25000000000000006, -0.5, 0.0, 0.0, -0.25000000000000006, -0.5,
    0.25000000000000006, 0.5, -0.25000000000000006, -0.5, 0.0, 0.0,
    -0.25000000000000006, -0.5, -0.25000000000000006, -0.5, 0.25000000000000006,
    0.5, 0.489234391703297, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mu1_1[618]{ -0.0, 0.0071318388640683669,
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

  static const real_T b_Ac[504]{ -0.144117, 0.0930431, -0.0, -0.28838933218494,
    0.18602178347057702, -0.0, -0.43281716397495551, 0.27893609500922684, -0.0,
    -0.577400662970631, 0.37178607918256912, -0.0, -0.72213999695319386,
    0.46457178052636866, -0.0, -0.86703533388470988, 0.55729324354555687, -0.0,
    -1.0120868419082774, 0.649950512714253, -0.0, -1.1572946893482228,
    0.74254363247578548, -0.0, -1.3026590447102959, 0.83507264724271346, -0.0,
    -1.4481800766818653, 0.92753760139684793, -0.0, -1.5938579541321145,
    1.0199385392892728, -0.0, -1.739692846112237, 1.1122755052403666, -0.0,
    -1.8856849218556335, 1.2045485435398235, -0.0, -2.0318343507781078,
    1.2967576984466747, -0.0, -2.1781413024780631, 1.3889030141893091, -0.0,
    -2.3246059467366997, 1.4809845349654953, -0.0, -2.4712284535182114,
    1.5730023049424027, -0.0, -2.6180089929699819, 1.6649563682566217, -0.0,
    -2.7649477354227843, 1.7568467690141867, -0.0, -2.9120448513909771,
    1.848673551290595, -0.0, 0.144117, -0.0930431, 0.0, 0.28838933218494,
    -0.18602178347057702, 0.0, 0.43281716397495551, -0.27893609500922684, 0.0,
    0.577400662970631, -0.37178607918256912, 0.0, 0.72213999695319386,
    -0.46457178052636866, 0.0, 0.86703533388470988, -0.55729324354555687, 0.0,
    1.0120868419082774, -0.649950512714253, 0.0, 1.1572946893482228,
    -0.74254363247578548, 0.0, 1.3026590447102959, -0.83507264724271346, 0.0,
    1.4481800766818653, -0.92753760139684793, 0.0, 1.5938579541321145,
    -1.0199385392892728, 0.0, 1.739692846112237, -1.1122755052403666, 0.0,
    1.8856849218556335, -1.2045485435398235, 0.0, 2.0318343507781078,
    -1.2967576984466747, 0.0, 2.1781413024780631, -1.3889030141893091, 0.0,
    2.3246059467366997, -1.4809845349654953, 0.0, 2.4712284535182114,
    -1.5730023049424027, 0.0, 2.6180089929699819, -1.6649563682566217, 0.0,
    2.7649477354227843, -1.7568467690141867, 0.0, 2.9120448513909771,
    -1.848673551290595, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.101003,
    -0.119576, -0.0, 0.20211486305346, -0.23906921394792, -0.0,
    0.30333570649515623, -0.35847969915902744, -0.0, 0.4046656477863308,
    -0.47780751290890866, -0.0, 0.50610480451482776, -0.59705271243349645, -0.0,
    0.60765329439522986, -0.71621535492909738, -0.0, 0.70931123526899487,
    -0.83529549755241927, -0.0, 0.81107874510459244, -0.95429319742059882, -0.0,
    0.912955941997641, -1.0732085116112287, -0.0, 1.0149429441710447,
    -1.1920414971623847, -0.0, 1.1170398699751309, -1.3107922110726542, -0.0,
    1.2192468378877874, -1.4294607103011623, -0.0, 1.3215639665145995,
    -1.5480470517675995, -0.0, 1.4239913745889881, -1.6665512923522492, -0.0,
    1.5265291809723474, -1.7849734888960149, -0.0, 1.6291775046541828,
    -1.9033136982004475, -0.0, 1.7319364647522488, -2.0215719770277722, -0.0,
    1.8348061805126878, -2.1397483821009167, -0.0, 1.9377867713101677,
    -2.2578429701035367, -0.0, 2.0408783566480211, -2.3758557976800447, -0.0,
    -0.101003, 0.119576, 0.0, -0.20211486305346, 0.23906921394792, 0.0,
    -0.30333570649515623, 0.35847969915902744, 0.0, -0.4046656477863308,
    0.47780751290890866, 0.0, -0.50610480451482776, 0.59705271243349645, 0.0,
    -0.60765329439522986, 0.71621535492909738, 0.0, -0.70931123526899487,
    0.83529549755241927, 0.0, -0.81107874510459244, 0.95429319742059882, 0.0,
    -0.912955941997641, 1.0732085116112287, 0.0, -1.0149429441710447,
    1.1920414971623847, 0.0, -1.1170398699751309, 1.3107922110726542, 0.0,
    -1.2192468378877874, 1.4294607103011623, 0.0, -1.3215639665145995,
    1.5480470517675995, 0.0, -1.4239913745889881, 1.6665512923522492, 0.0,
    -1.5265291809723474, 1.7849734888960149, 0.0, -1.6291775046541828,
    1.9033136982004475, 0.0, -1.7319364647522488, 2.0215719770277722, 0.0,
    -1.8348061805126878, 2.1397483821009167, 0.0, -1.9377867713101677,
    2.2578429701035367, 0.0, -2.0408783566480211, 2.3758557976800447, 0.0, -0.0,
    -1.0, -0.0, 0.0, 1.0, 0.0, 0.0275178, 0.0102807, -0.0, 0.055065259235195993,
    0.020554282362969002, -0.0, 0.082642409672904871, 0.030820752016660648, -0.0,
    0.11024928331489851, 0.04108011388541695, -0.0, 0.13788591219744095,
    0.051332372890170661, -0.0, 0.16555232839132558, 0.061577533948447609, -0.0,
    0.19324856400191229, 0.071815601974369078, -0.0, 0.22097465116916482,
    0.082046581878654157, -0.0, 0.24873062206768792, 0.092270478568622111, -0.0,
    0.27651650890676488, 0.10248729694819469, -0.0, 0.30433234393039471,
    0.11269704191789855, -0.0, 0.33217815941732975, 0.12289971837486753, -0.0,
    0.36005398768111291, 0.13309533121284506, -0.0, 0.38795986107011532,
    0.14328388532218647, -0.0, 0.41589581196757386, 0.15346538558986136, -0.0,
    0.4438618727916287, 0.16363983689945594, -0.0, 0.47185807599536089,
    0.17380724413117532, -0.0, 0.49988445406683013, 0.18396761216184598, -0.0,
    0.52794103952911242, 0.19412094586491796, -0.0, 0.55602786494033762,
    0.20426725011046729, -0.0, -0.0275178, -0.0102807, 0.0,
    -0.055065259235195993, -0.020554282362969002, 0.0, -0.082642409672904871,
    -0.030820752016660648, 0.0, -0.11024928331489851, -0.04108011388541695, 0.0,
    -0.13788591219744095, -0.051332372890170661, 0.0, -0.16555232839132558,
    -0.061577533948447609, 0.0, -0.19324856400191229, -0.071815601974369078, 0.0,
    -0.22097465116916482, -0.082046581878654157, 0.0, -0.24873062206768792,
    -0.092270478568622111, 0.0, -0.27651650890676488, -0.10248729694819469, 0.0,
    -0.30433234393039471, -0.11269704191789855, 0.0, -0.33217815941732975,
    -0.12289971837486753, 0.0, -0.36005398768111291, -0.13309533121284506, 0.0,
    -0.38795986107011532, -0.14328388532218647, 0.0, -0.41589581196757386,
    -0.15346538558986136, 0.0, -0.4438618727916287, -0.16363983689945594, 0.0,
    -0.47185807599536089, -0.17380724413117532, 0.0, -0.49988445406683013,
    -0.18396761216184598, 0.0, -0.52794103952911242, -0.19412094586491796, 0.0,
    -0.55602786494033762, -0.20426725011046729, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b_Mu1_0[498]{ -0.22599754912864506, -0.0, -0.0, -0.0,
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

  static const real_T c[480]{ 1.00107782, 0.0, 0.0125, 0.5, -0.0125, -0.5,
    -0.0125, -0.5, 0.0, 0.99930767, -0.0125, -0.5, 0.0125, 0.5, -0.0125, -0.5,
    0.0, 0.0, -0.0125, -0.5, -0.0125, -0.5, 0.0125, 0.5, 1.0021568016959521, 0.0,
    0.025, 0.5, -0.025, -0.5, -0.025, -0.5, 0.0, 0.99861581932082888, -0.025,
    -0.5, 0.025, 0.5, -0.025, -0.5, 0.0, 0.0, -0.025, -0.5, -0.025, -0.5, 0.025,
    0.5, 1.003236946339956, 0.0, 0.037500000000000006, 0.5,
    -0.037500000000000006, -0.5, -0.037500000000000006, -0.5, 0.0,
    0.99792444763063848, -0.037500000000000006, -0.5, 0.037500000000000006, 0.5,
    -0.037500000000000006, -0.5, 0.0, 0.0, -0.037500000000000006, -0.5,
    -0.037500000000000006, -0.5, 0.037500000000000006, 0.5, 1.00431825518546,
    0.0, 0.05, 0.5, -0.05, -0.5, -0.05, -0.5, 0.0, 0.99723355459781038, -0.05,
    -0.5, 0.05, 0.5, -0.05, -0.5, 0.0, 0.0, -0.05, -0.5, -0.05, -0.5, 0.05, 0.5,
    1.0054007294872638, 0.0, 0.0625, 0.5, -0.0625, -0.5, -0.0625, -0.5, 0.0,
    0.99654313989095566, -0.0625, -0.5, 0.0625, 0.5, -0.0625, -0.5, 0.0, 0.0,
    -0.0625, -0.5, -0.0625, -0.5, 0.0625, 0.5, 1.0064843705015196, 0.0, 0.075,
    0.5, -0.075, -0.5, -0.075, -0.5, 0.0, 0.9958532031789149, -0.075, -0.5,
    0.075, 0.5, -0.075, -0.5, 0.0, 0.0, -0.075, -0.5, -0.075, -0.5, 0.075, 0.5,
    1.0075691794857335, 0.0, 0.0875, 0.5, -0.0875, -0.5, -0.0875, -0.5, 0.0,
    0.995163744130758, -0.0875, -0.5, 0.0875, 0.5, -0.0875, -0.5, 0.0, 0.0,
    -0.0875, -0.5, -0.0875, -0.5, 0.0875, 0.5, 1.0086551576987668, 0.0,
    0.099999999999999992, 0.5, -0.099999999999999992, -0.5,
    -0.099999999999999992, -0.5, 0.0, 0.99447476241578392, -0.099999999999999992,
    -0.5, 0.099999999999999992, 0.5, -0.099999999999999992, -0.5, 0.0, 0.0,
    -0.099999999999999992, -0.5, -0.099999999999999992, -0.5,
    0.099999999999999992, 0.5, 1.0097423064008375, 0.0, 0.11249999999999999, 0.5,
    -0.11249999999999999, -0.5, -0.11249999999999999, -0.5, 0.0,
    0.99378625770352058, -0.11249999999999999, -0.5, 0.11249999999999999, 0.5,
    -0.11249999999999999, -0.5, 0.0, 0.0, -0.11249999999999999, -0.5,
    -0.11249999999999999, -0.5, 0.11249999999999999, 0.5, 1.0108306268535223,
    0.0, 0.12499999999999999, 0.5, -0.12499999999999999, -0.5,
    -0.12499999999999999, -0.5, 0.0, 0.99309822966372463, -0.12499999999999999,
    -0.5, 0.12499999999999999, 0.5, -0.12499999999999999, -0.5, 0.0, 0.0,
    -0.12499999999999999, -0.5, -0.12499999999999999, -0.5, 0.12499999999999999,
    0.5, 1.0119201203197574, 0.0, 0.13749999999999998, 0.5, -0.13749999999999998,
    -0.5, -0.13749999999999998, -0.5, 0.0, 0.99241067796638149,
    -0.13749999999999998, -0.5, 0.13749999999999998, 0.5, -0.13749999999999998,
    -0.5, 0.0, 0.0, -0.13749999999999998, -0.5, -0.13749999999999998, -0.5,
    0.13749999999999998, 0.5, 1.0130107880638404, 0.0, 0.15, 0.5, -0.15, -0.5,
    -0.15, -0.5, 0.0, 0.991723602281705, -0.15, -0.5, 0.15, 0.5, -0.15, -0.5,
    0.0, 0.0, -0.15, -0.5, -0.15, -0.5, 0.15, 0.5, 1.0141026313514312, 0.0,
    0.1625, 0.5, -0.1625, -0.5, -0.1625, -0.5, 0.0, 0.99103700228013725, -0.1625,
    -0.5, 0.1625, 0.5, -0.1625, -0.5, 0.0, 0.0, -0.1625, -0.5, -0.1625, -0.5,
    0.1625, 0.5, 1.0151956514495544, 0.0, 0.17500000000000002, 0.5,
    -0.17500000000000002, -0.5, -0.17500000000000002, -0.5, 0.0,
    0.99035087763234864, -0.17500000000000002, -0.5, 0.17500000000000002, 0.5,
    -0.17500000000000002, -0.5, 0.0, 0.0, -0.17500000000000002, -0.5,
    -0.17500000000000002, -0.5, 0.17500000000000002, 0.5, 1.0162898496265995,
    0.0, 0.18750000000000003, 0.5, -0.18750000000000003, -0.5,
    -0.18750000000000003, -0.5, 0.0, 0.9896652280092374, -0.18750000000000003,
    -0.5, 0.18750000000000003, 0.5, -0.18750000000000003, -0.5, 0.0, 0.0,
    -0.18750000000000003, -0.5, -0.18750000000000003, -0.5, 0.18750000000000003,
    0.5, 1.017385227152324, 0.0, 0.20000000000000004, 0.5, -0.20000000000000004,
    -0.5, -0.20000000000000004, -0.5, 0.0, 0.98898005308192971,
    -0.20000000000000004, -0.5, 0.20000000000000004, 0.5, -0.20000000000000004,
    -0.5, 0.0, 0.0, -0.20000000000000004, -0.5, -0.20000000000000004, -0.5,
    0.20000000000000004, 0.5, 1.0184817852978532, 0.0, 0.21250000000000005, 0.5,
    -0.21250000000000005, -0.5, -0.21250000000000005, -0.5, 0.0,
    0.98829535252177947, -0.21250000000000005, -0.5, 0.21250000000000005, 0.5,
    -0.21250000000000005, -0.5, 0.0, 0.0, -0.21250000000000005, -0.5,
    -0.21250000000000005, -0.5, 0.21250000000000005, 0.5, 1.0195795253356827,
    0.0, 0.22500000000000006, 0.5, -0.22500000000000006, -0.5,
    -0.22500000000000006, -0.5, 0.0, 0.987611126000368, -0.22500000000000006,
    -0.5, 0.22500000000000006, 0.5, -0.22500000000000006, -0.5, 0.0, 0.0,
    -0.22500000000000006, -0.5, -0.22500000000000006, -0.5, 0.22500000000000006,
    0.5, 1.02067844853968, 0.0, 0.23750000000000007, 0.5, -0.23750000000000007,
    -0.5, -0.23750000000000007, -0.5, 0.0, 0.98692737318950419,
    -0.23750000000000007, -0.5, 0.23750000000000007, 0.5, -0.23750000000000007,
    -0.5, 0.0, 0.0, -0.23750000000000007, -0.5, -0.23750000000000007, -0.5,
    0.23750000000000007, 0.5, 1.021778556185085, 0.0, 0.25000000000000006, 0.5,
    -0.25000000000000006, -0.5, -0.25000000000000006, -0.5, 0.0,
    0.98624409376122391, -0.25000000000000006, -0.5, 0.25000000000000006, 0.5,
    -0.25000000000000006, -0.5, 0.0, 0.0, -0.25000000000000006, -0.5,
    -0.25000000000000006, -0.5, 0.25000000000000006, 0.5 };

  static const real_T b_Mu1[378]{ -0.144117, 0.0930431, -0.0, -0.28838933218494,
    0.18602178347057702, -0.0, -0.43281716397495551, 0.27893609500922684, -0.0,
    -0.577400662970631, 0.37178607918256912, -0.0, -0.72213999695319386,
    0.46457178052636866, -0.0, -0.86703533388470988, 0.55729324354555687, -0.0,
    -1.0120868419082774, 0.649950512714253, -0.0, -1.1572946893482228,
    0.74254363247578548, -0.0, -1.3026590447102959, 0.83507264724271346, -0.0,
    -1.4481800766818653, 0.92753760139684793, -0.0, -1.5938579541321145,
    1.0199385392892728, -0.0, -1.739692846112237, 1.1122755052403666, -0.0,
    -1.8856849218556335, 1.2045485435398235, -0.0, -2.0318343507781078,
    1.2967576984466747, -0.0, -2.1781413024780631, 1.3889030141893091, -0.0,
    -2.3246059467366997, 1.4809845349654953, -0.0, -2.4712284535182114,
    1.5730023049424027, -0.0, -2.6180089929699819, 1.6649563682566217, -0.0,
    -2.7649477354227843, 1.7568467690141867, -0.0, -2.9120448513909771,
    1.848673551290595, -0.0, 0.144117, -0.0930431, 0.0, 0.28838933218494,
    -0.18602178347057702, 0.0, 0.43281716397495551, -0.27893609500922684, 0.0,
    0.577400662970631, -0.37178607918256912, 0.0, 0.72213999695319386,
    -0.46457178052636866, 0.0, 0.86703533388470988, -0.55729324354555687, 0.0,
    1.0120868419082774, -0.649950512714253, 0.0, 1.1572946893482228,
    -0.74254363247578548, 0.0, 1.3026590447102959, -0.83507264724271346, 0.0,
    1.4481800766818653, -0.92753760139684793, 0.0, 1.5938579541321145,
    -1.0199385392892728, 0.0, 1.739692846112237, -1.1122755052403666, 0.0,
    1.8856849218556335, -1.2045485435398235, 0.0, 2.0318343507781078,
    -1.2967576984466747, 0.0, 2.1781413024780631, -1.3889030141893091, 0.0,
    2.3246059467366997, -1.4809845349654953, 0.0, 2.4712284535182114,
    -1.5730023049424027, 0.0, 2.6180089929699819, -1.6649563682566217, 0.0,
    2.7649477354227843, -1.7568467690141867, 0.0, 2.9120448513909771,
    -1.848673551290595, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.101003,
    -0.119576, -0.0, 0.20211486305346, -0.23906921394792, -0.0,
    0.30333570649515623, -0.35847969915902744, -0.0, 0.4046656477863308,
    -0.47780751290890866, -0.0, 0.50610480451482776, -0.59705271243349645, -0.0,
    0.60765329439522986, -0.71621535492909738, -0.0, 0.70931123526899487,
    -0.83529549755241927, -0.0, 0.81107874510459244, -0.95429319742059882, -0.0,
    0.912955941997641, -1.0732085116112287, -0.0, 1.0149429441710447,
    -1.1920414971623847, -0.0, 1.1170398699751309, -1.3107922110726542, -0.0,
    1.2192468378877874, -1.4294607103011623, -0.0, 1.3215639665145995,
    -1.5480470517675995, -0.0, 1.4239913745889881, -1.6665512923522492, -0.0,
    1.5265291809723474, -1.7849734888960149, -0.0, 1.6291775046541828,
    -1.9033136982004475, -0.0, 1.7319364647522488, -2.0215719770277722, -0.0,
    1.8348061805126878, -2.1397483821009167, -0.0, 1.9377867713101677,
    -2.2578429701035367, -0.0, 2.0408783566480211, -2.3758557976800447, -0.0,
    -0.101003, 0.119576, 0.0, -0.20211486305346, 0.23906921394792, 0.0,
    -0.30333570649515623, 0.35847969915902744, 0.0, -0.4046656477863308,
    0.47780751290890866, 0.0, -0.50610480451482776, 0.59705271243349645, 0.0,
    -0.60765329439522986, 0.71621535492909738, 0.0, -0.70931123526899487,
    0.83529549755241927, 0.0, -0.81107874510459244, 0.95429319742059882, 0.0,
    -0.912955941997641, 1.0732085116112287, 0.0, -1.0149429441710447,
    1.1920414971623847, 0.0, -1.1170398699751309, 1.3107922110726542, 0.0,
    -1.2192468378877874, 1.4294607103011623, 0.0, -1.3215639665145995,
    1.5480470517675995, 0.0, -1.4239913745889881, 1.6665512923522492, 0.0,
    -1.5265291809723474, 1.7849734888960149, 0.0, -1.6291775046541828,
    1.9033136982004475, 0.0, -1.7319364647522488, 2.0215719770277722, 0.0,
    -1.8348061805126878, 2.1397483821009167, 0.0, -1.9377867713101677,
    2.2578429701035367, 0.0, -2.0408783566480211, 2.3758557976800447, 0.0, -0.0,
    -1.0, -0.0, 0.0, 1.0, 0.0, 0.0275178, 0.0102807, -0.0, 0.055065259235195993,
    0.020554282362969002, -0.0, 0.082642409672904871, 0.030820752016660648, -0.0,
    0.11024928331489851, 0.04108011388541695, -0.0, 0.13788591219744095,
    0.051332372890170661, -0.0, 0.16555232839132558, 0.061577533948447609, -0.0,
    0.19324856400191229, 0.071815601974369078, -0.0, 0.22097465116916482,
    0.082046581878654157, -0.0, 0.24873062206768792, 0.092270478568622111, -0.0,
    0.27651650890676488, 0.10248729694819469, -0.0, 0.30433234393039471,
    0.11269704191789855, -0.0, 0.33217815941732975, 0.12289971837486753, -0.0,
    0.36005398768111291, 0.13309533121284506, -0.0, 0.38795986107011532,
    0.14328388532218647, -0.0, 0.41589581196757386, 0.15346538558986136, -0.0,
    0.4438618727916287, 0.16363983689945594, -0.0, 0.47185807599536089,
    0.17380724413117532, -0.0, 0.49988445406683013, 0.18396761216184598, -0.0,
    0.52794103952911242, 0.19412094586491796, -0.0, 0.55602786494033762,
    0.20426725011046729, -0.0, -0.0275178, -0.0102807, 0.0,
    -0.055065259235195993, -0.020554282362969002, 0.0, -0.082642409672904871,
    -0.030820752016660648, 0.0, -0.11024928331489851, -0.04108011388541695, 0.0,
    -0.13788591219744095, -0.051332372890170661, 0.0, -0.16555232839132558,
    -0.061577533948447609, 0.0, -0.19324856400191229, -0.071815601974369078, 0.0,
    -0.22097465116916482, -0.082046581878654157, 0.0, -0.24873062206768792,
    -0.092270478568622111, 0.0, -0.27651650890676488, -0.10248729694819469, 0.0,
    -0.30433234393039471, -0.11269704191789855, 0.0, -0.33217815941732975,
    -0.12289971837486753, 0.0, -0.36005398768111291, -0.13309533121284506, 0.0,
    -0.38795986107011532, -0.14328388532218647, 0.0, -0.41589581196757386,
    -0.15346538558986136, 0.0, -0.4438618727916287, -0.16363983689945594, 0.0,
    -0.47185807599536089, -0.17380724413117532, 0.0, -0.49988445406683013,
    -0.18396761216184598, 0.0, -0.52794103952911242, -0.19412094586491796, 0.0,
    -0.55602786494033762, -0.20426725011046729, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0,
    1.0 };

  static const real_T b_SuJm_1[300]{ 0.0, -0.0071318388640683669,
    -0.012368697539837475, 0.0, 0.0, 0.0, -0.014263756082690558,
    -0.024718035607500709, -0.00017829597160170919, -0.00030921743849593691, 0.0,
    -0.021395752912117184, -0.037048025614083456, -0.00053488987366897319,
    -0.00092716832868345472, 0.0, -0.028527830608018147, -0.0493586789617647,
    -0.0010697836964719027, -0.0018533689690355412, 0.0, -0.035659990425483137,
    -0.061650007043814951, -0.0017829794616723562, -0.0030873359430796592, 0.0,
    -0.042792233619022145, -0.073922021244602523, -0.0026744792223094348,
    -0.0046285861191750332, 0.0, -0.049924561442565955, -0.086174732939599832,
    -0.0037442850627849886, -0.006476636650290097, 0.0, -0.057056975149466643,
    -0.098408153495389644, -0.004992399098849137, -0.0086310049737800924, 0.0,
    -0.0641894759924981, -0.11062229426967139, -0.0064188234775858031,
    -0.011091208811164835, 0.0, -0.071322065223856493, -0.1228171666112674,
    -0.0080235603773982558, -0.01385676616790662, 0.0, -0.078454744095160817,
    -0.1349927818601292, -0.0098066120079946684, -0.016927195333188306, 0.0,
    -0.085587513857453348, -0.14714915134734369, -0.011767980610373688,
    -0.020302014879691535, 0.0, -0.09272037576120018, -0.1592862863951395,
    -0.013907668456810022, -0.023980743663375128, 0.0, -0.099853331056291711,
    -0.17140419831689321, -0.016225677850840024, -0.027962900823253617, 0.0,
    -0.10698638099204312, -0.18350289841713557, -0.018722011127247317,
    -0.032248005781175947, 0.0, -0.11411952681719491, -0.19558239799155772,
    -0.021396670652048394, -0.036835578241604335, 0.0, -0.12125276977991335,
    -0.20764270832701751, -0.024249658822478266, -0.041725138191393278, 0.0,
    -0.12838611112779105, -0.21968384070154565, -0.0272809780669761,
    -0.046916205899568721, 0.0, -0.13551955210784739, -0.23170580638435198,
    -0.030490630845170875, -0.052408301917107367, 0.0, -0.14265309396652903,
    -0.24370861663583168, -0.033878619647867057, -0.058200947076716171, 0.0,
    0.13659152296751884, -0.063538603002101665, 0.0, 0.0, 0.0,
    0.27320389247378107, -0.12726886009414035, 0.0034147880741879711,
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
    -0.30726093601542542, 0.0, -0.04588890657825119, 0.10453377294976968, 0.0,
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
    0.49695833722159954 };

  static const real_T e_1[300]{ 0.0, 0.0, 0.0, -0.0071318388640683669,
    0.13659152296751884, -0.04588890657825119, -0.012368697539837475,
    -0.063538603002101665, 0.10453377294976968, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.014263756082690558, 0.27320389247378107,
    -0.091790173501845068, -0.024718035607500709, -0.12726886009414035,
    0.2090823342837666, -0.00017829597160170919, 0.0034147880741879711,
    -0.0011472226644562798, -0.00030921743849593691, -0.0015884650750525417,
    0.002613344323744242, 0.0, 0.0, 0.0, -0.021395752912117184,
    0.40983712341735679, -0.13770380323925227, -0.037048025614083456,
    -0.19119069203750322, 0.31364569595584796, -0.00053488987366897319,
    0.010244885386032498, -0.0034419770020024066, -0.00092716832868345472,
    -0.00477018657740605, 0.007840402680838408, 0.0, 0.0, 0.0,
    -0.028527830608018147, 0.54649123069354111, -0.183629798260016,
    -0.0493586789617647, -0.25530401966585881, 0.41822386991689187,
    -0.0010697836964719027, 0.020490813471466418, -0.0068845720829837139,
    -0.0018533689690355412, -0.00954995387834363, 0.015681545079734606, 0.0, 0.0,
    0.0, -0.035659990425483137, 0.68316622919435854, -0.22956816103475194,
    -0.061650007043814951, -0.31960876388510828, 0.52281686811480088,
    -0.0017829794616723562, 0.03415309423880495, -0.011475317039484116,
    -0.0030873359430796592, -0.0159325543699901, 0.026137141827656904, 0.0, 0.0,
    0.0, -0.042792233619022145, 0.81986213380856687, -0.27551889403514823,
    -0.073922021244602523, -0.38410484567333708, 0.62742470249450566,
    -0.0026744792223094348, 0.051232249968663914, -0.017214521065352914,
    -0.0046285861191750332, -0.023922773467117804, 0.039207563530526922, 0.0,
    0.0, 0.0, -0.049924561442565955, 0.95657895942166182, -0.32148199973396535,
    -0.086174732939599832, -0.44879218608076621, 0.7320473849979684,
    -0.0037442850627849886, 0.071728803313878084, -0.024102493416231621,
    -0.006476636650290097, -0.03352539460895123, 0.054893181092889556, 0.0, 0.0,
    0.0, -0.057056975149466643, 1.0933167209158812, -0.36745748060503614,
    -0.098408153495389644, -0.51367070622970379, 0.83668492756418655,
    -0.004992399098849137, 0.095643277299419627, -0.032139543409580758,
    -0.0086310049737800924, -0.044745199260970384, 0.073194365717838758, 0.0,
    0.0, 0.0, -0.0641894759924981, 1.2300754331702088, -0.41344533912326564,
    -0.11062229426967139, -0.57874032731449654, 0.9413373421291964,
    -0.0064188234775858031, 0.12297619532231666, -0.041325980424706657,
    -0.011091208811164835, -0.057586966916712973, 0.094111488906943416, 0.0, 0.0,
    0.0, -0.071322065223856493, 1.3668551110603793, -0.45944557776463113,
    -0.1228171666112674, -0.64400097060148154, 1.0460046406260766,
    -0.0080235603773982558, 0.15372808115157188, -0.0516621139027883,
    -0.01385676616790662, -0.072055475099575383, 0.11764492246017333, 0.0, 0.0,
    0.0, -0.078454744095160817, 1.5036557694588819, -0.505458199006182,
    -0.1349927818601292, -0.70945255742893776, 1.1506868349849515,
    -0.0098066120079946684, 0.18789945892808135, -0.063148253346904076,
    -0.016927195333188306, -0.088155499364612422, 0.14379503847582525, 0.0, 0.0,
    0.0, -0.085587513857453348, 1.6404774232349646, -0.55148320532603978,
    -0.14714915134734369, -0.77509500920703767, 1.2553839371329951,
    -0.011767980610373688, 0.22549085316455342, -0.075784708322058622,
    -0.020302014879691535, -0.10589181330033587, 0.17256220935044903, 0.0, 0.0,
    0.0, -0.09272037576120018, 1.7773200872546391, -0.597520599203398,
    -0.1592862863951395, -0.84092824741779892, 1.3600959589944344,
    -0.013907668456810022, 0.26650278874542754, -0.089571788455209619,
    -0.023980743663375128, -0.1252691885305118, 0.20394680777877389, 0.0, 0.0,
    0.0, -0.099853331056291711, 1.9141837763806839, -0.643570383118522,
    -0.17140419831689321, -0.90695219361503621, 1.4648229124905527,
    -0.016225677850840024, 0.31093579092679352, -0.10450980343529456,
    -0.027962900823253617, -0.14629239471595679, 0.23794920675363476, 0.0, 0.0,
    0.0, -0.10698638099204312, 2.0510685054726494, -0.68963255955274927,
    -0.18350289841713557, -0.97316676942431268, 1.5695648095396941,
    -0.018722011127247317, 0.35879038533631064, -0.12059906301325761,
    -0.032248005781175947, -0.16896619955633269, 0.27456977956589856, 0.0, 0.0,
    0.0, -0.11411952681719491, 2.1879742893868621, -0.73570713098848894,
    -0.19558239799155772, -1.0395718965428922, 1.6743216620572658,
    -0.021396670652048394, 0.4100670979731269, -0.13783987700207634,
    -0.036835578241604335, -0.19329536879194051, 0.31380889980439092, 0.0, 0.0,
    0.0, -0.12125276977991335, 2.3249011429764286, -0.78179409990922211,
    -0.20764270832701751, -1.1061674967396908, 1.7790934819557425,
    -0.024249658822478266, 0.46476645520779847, -0.15623255527678856,
    -0.041725138191393278, -0.21928466620551282, 0.35566694135582255, 0.0, 0.0,
    0.0, -0.12838611112779105, 2.4618490810912395, -0.82789346879950143,
    -0.21968384070154565, -1.1729534918552289, 1.8838802811446695,
    -0.0272809780669761, 0.52288898378220916, -0.17577740777451911,
    -0.046916205899568721, -0.24693885362400511, 0.40014427840471611, 0.0, 0.0,
    0.0, -0.13551955210784739, 2.5988181185779742, -0.8740052401449514,
    -0.23170580638435198, -1.2399298038015829, 1.9886820715306668,
    -0.030490630845170875, 0.58443521080949012, -0.19647474449450664,
    -0.052408301917107367, -0.27626269092038586, 0.44724128543333286, 0.0, 0.0,
    0.0, -0.14265309396652903, 2.735808270280105, -0.920129416432268,
    -0.24370861663583168, -1.3070963545623373, 2.0934988650174318,
    -0.033878619647867057, 0.64940566377393949, -0.21832487549813043,
    -0.058200947076716171, -0.30726093601542542, 0.49695833722159954 };

  static const real_T b_SuJm[240]{ 0.22599754912864506, 0.0, 0.0, 0.0,
    0.45147579735946186, 0.0, 0.0, 0.0056499387282161265, 0.676435937950623, 0.0,
    0.0, 0.016936833662202673, 0.90087916141841273, 0.0, 0.0,
    0.033847732110968247, 1.1248066555435274, 0.0, 0.0, 0.056369711146428567,
    1.3482196053773612, 0.0, 0.0, 0.084489877535016744, 1.5711191932482778, 0.0,
    0.0, 0.11819536766945077, 1.7935065987678671, 0.0, 0.0, 0.15747334750065772,
    2.0153829988371879, 0.0, 0.0, 0.20231101246985439, 2.2367495676529967, 0.0,
    0.0, 0.25269558744078408, 2.4576074767139611, 0.0, 0.0, 0.30861432663210897,
    2.67795789482686, 0.0, 0.0, 0.370054513549958, 2.8978019881127679, 0.0, 0.0,
    0.43700346092062947, 3.1171409200132283, 0.0, 0.0, 0.50944851062344865,
    3.3359758512964088, 0.0, 0.0, 0.58737703362377935, 3.5543079400632447, 0.0,
    0.0, 0.67077642990618958, 3.7721383417535681, 0.0, 0.0, 0.75963412840777078,
    3.9894682091522222, 0.0, 0.0, 0.85393758695161, 4.2062986923951629, 0.0, 0.0,
    0.95367429218041566, 4.4226309389755434, 0.0, 0.0, 1.0588317594902947,
    -0.087367622049086546, 0.0, 0.0, 0.0, -0.17453448933447702, 0.0, 0.0,
    -0.0021841905512271637, -0.26150106315378036, 0.0, 0.0, -0.00654755278458909,
    -0.34826780374462818, 0.0, 0.0, -0.013085079363433599, -0.4348351702871105,
    0.0, 0.0, -0.0217917744570493, -0.52120362090620576, 0.0, 0.0,
    -0.032662653714227059, -0.6073736126742052, 0.0, 0.0, -0.0456927442368822,
    -0.69334560161313186, 0.0, 0.0, -0.060877084553737332, -0.779120042697154,
    0.0, 0.0, -0.078210724594065625, -0.86469738985499234, 0.0, 0.0,
    -0.097688725661494474, -0.95007809597232318, 0.0, 0.0, -0.11930616040786929,
    -1.0352626128941742, 0.0, 0.0, -0.14305811280717737, -1.1202513914273164,
    0.0, 0.0, -0.16893967812953173, -1.2050448813426495, 0.0, 0.0,
    -0.19694596291521463, -1.2896435313775823, 0.0, 0.0, -0.22707208494878087,
    -1.3740477892384069, 0.0, 0.0, -0.25931317323322045, -1.458258101602669, 0.0,
    0.0, -0.29366436796418061, -1.5422749141215313, 0.0, 0.0,
    -0.33012082050424735, -1.6260986714221315, 0.0, 0.0, -0.36867769335728562,
    -1.709729817109936, 0.0, 0.0, -0.40933016014283891, -0.087367622049086546,
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
    0.0, -0.40933016014283891 };

  static const real_T e[240]{ 0.22599754912864506, -0.087367622049086546,
    -0.087367622049086546, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.45147579735946186, -0.17453448933447702, -0.17453448933447702, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0056499387282161265, -0.0021841905512271637,
    -0.0021841905512271637, 0.676435937950623, -0.26150106315378036,
    -0.26150106315378036, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.016936833662202673,
    -0.00654755278458909, -0.00654755278458909, 0.90087916141841273,
    -0.34826780374462818, -0.34826780374462818, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.033847732110968247, -0.013085079363433599, -0.013085079363433599,
    1.1248066555435274, -0.4348351702871105, -0.4348351702871105, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.056369711146428567, -0.0217917744570493,
    -0.0217917744570493, 1.3482196053773612, -0.52120362090620576,
    -0.52120362090620576, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.084489877535016744,
    -0.032662653714227059, -0.032662653714227059, 1.5711191932482778,
    -0.6073736126742052, -0.6073736126742052, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.11819536766945077, -0.0456927442368822, -0.0456927442368822,
    1.7935065987678671, -0.69334560161313186, -0.69334560161313186, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.15747334750065772, -0.060877084553737332,
    -0.060877084553737332, 2.0153829988371879, -0.779120042697154,
    -0.779120042697154, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.20231101246985439,
    -0.078210724594065625, -0.078210724594065625, 2.2367495676529967,
    -0.86469738985499234, -0.86469738985499234, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.25269558744078408, -0.097688725661494474, -0.097688725661494474,
    2.4576074767139611, -0.95007809597232318, -0.95007809597232318, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.30861432663210897, -0.11930616040786929,
    -0.11930616040786929, 2.67795789482686, -1.0352626128941742,
    -1.0352626128941742, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.370054513549958,
    -0.14305811280717737, -0.14305811280717737, 2.8978019881127679,
    -1.1202513914273164, -1.1202513914273164, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.43700346092062947, -0.16893967812953173, -0.16893967812953173,
    3.1171409200132283, -1.2050448813426495, -1.2050448813426495, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.50944851062344865, -0.19694596291521463,
    -0.19694596291521463, 3.3359758512964088, -1.2896435313775823,
    -1.2896435313775823, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.58737703362377935,
    -0.22707208494878087, -0.22707208494878087, 3.5543079400632447,
    -1.3740477892384069, -1.3740477892384069, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.67077642990618958, -0.25931317323322045, -0.25931317323322045,
    3.7721383417535681, -1.458258101602669, -1.458258101602669, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.75963412840777078, -0.29366436796418061,
    -0.29366436796418061, 3.9894682091522222, -1.5422749141215313,
    -1.5422749141215313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.85393758695161,
    -0.33012082050424735, -0.33012082050424735, 4.2062986923951629,
    -1.6260986714221315, -1.6260986714221315, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.95367429218041566, -0.36867769335728562, -0.36867769335728562,
    4.4226309389755434, -1.709729817109936, -1.709729817109936, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0588317594902947, -0.40933016014283891,
    -0.40933016014283891 };

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

  static const real_T b_SuJm_0[180]{ 0.144117, -0.0930431, 0.0, 0.28838933218494,
    -0.18602178347057702, 0.0, 0.43281716397495551, -0.27893609500922684, 0.0,
    0.577400662970631, -0.37178607918256912, 0.0, 0.72213999695319386,
    -0.46457178052636866, 0.0, 0.86703533388470988, -0.55729324354555687, 0.0,
    1.0120868419082774, -0.649950512714253, 0.0, 1.1572946893482228,
    -0.74254363247578548, 0.0, 1.3026590447102959, -0.83507264724271346, 0.0,
    1.4481800766818653, -0.92753760139684793, 0.0, 1.5938579541321145,
    -1.0199385392892728, 0.0, 1.739692846112237, -1.1122755052403666, 0.0,
    1.8856849218556335, -1.2045485435398235, 0.0, 2.0318343507781078,
    -1.2967576984466747, 0.0, 2.1781413024780631, -1.3889030141893091, 0.0,
    2.3246059467366997, -1.4809845349654953, 0.0, 2.4712284535182114,
    -1.5730023049424027, 0.0, 2.6180089929699819, -1.6649563682566217, 0.0,
    2.7649477354227843, -1.7568467690141867, 0.0, 2.9120448513909771,
    -1.848673551290595, 0.0, -0.101003, 0.119576, 0.0, -0.20211486305346,
    0.23906921394792, 0.0, -0.30333570649515623, 0.35847969915902744, 0.0,
    -0.4046656477863308, 0.47780751290890866, 0.0, -0.50610480451482776,
    0.59705271243349645, 0.0, -0.60765329439522986, 0.71621535492909738, 0.0,
    -0.70931123526899487, 0.83529549755241927, 0.0, -0.81107874510459244,
    0.95429319742059882, 0.0, -0.912955941997641, 1.0732085116112287, 0.0,
    -1.0149429441710447, 1.1920414971623847, 0.0, -1.1170398699751309,
    1.3107922110726542, 0.0, -1.2192468378877874, 1.4294607103011623, 0.0,
    -1.3215639665145995, 1.5480470517675995, 0.0, -1.4239913745889881,
    1.6665512923522492, 0.0, -1.5265291809723474, 1.7849734888960149, 0.0,
    -1.6291775046541828, 1.9033136982004475, 0.0, -1.7319364647522488,
    2.0215719770277722, 0.0, -1.8348061805126878, 2.1397483821009167, 0.0,
    -1.9377867713101677, 2.2578429701035367, 0.0, -2.0408783566480211,
    2.3758557976800447, 0.0, -0.0275178, -0.0102807, 0.0, -0.055065259235195993,
    -0.020554282362969002, 0.0, -0.082642409672904871, -0.030820752016660648,
    0.0, -0.11024928331489851, -0.04108011388541695, 0.0, -0.13788591219744095,
    -0.051332372890170661, 0.0, -0.16555232839132558, -0.061577533948447609, 0.0,
    -0.19324856400191229, -0.071815601974369078, 0.0, -0.22097465116916482,
    -0.082046581878654157, 0.0, -0.24873062206768792, -0.092270478568622111, 0.0,
    -0.27651650890676488, -0.10248729694819469, 0.0, -0.30433234393039471,
    -0.11269704191789855, 0.0, -0.33217815941732975, -0.12289971837486753, 0.0,
    -0.36005398768111291, -0.13309533121284506, 0.0, -0.38795986107011532,
    -0.14328388532218647, 0.0, -0.41589581196757386, -0.15346538558986136, 0.0,
    -0.4438618727916287, -0.16363983689945594, 0.0, -0.47185807599536089,
    -0.17380724413117532, 0.0, -0.49988445406683013, -0.18396761216184598, 0.0,
    -0.52794103952911242, -0.19412094586491796, 0.0, -0.55602786494033762,
    -0.20426725011046729, 0.0 };

  static const real_T e_0[180]{ 0.144117, -0.101003, -0.0275178, -0.0930431,
    0.119576, -0.0102807, 0.0, 0.0, 0.0, 0.28838933218494, -0.20211486305346,
    -0.055065259235195993, -0.18602178347057702, 0.23906921394792,
    -0.020554282362969002, 0.0, 0.0, 0.0, 0.43281716397495551,
    -0.30333570649515623, -0.082642409672904871, -0.27893609500922684,
    0.35847969915902744, -0.030820752016660648, 0.0, 0.0, 0.0, 0.577400662970631,
    -0.4046656477863308, -0.11024928331489851, -0.37178607918256912,
    0.47780751290890866, -0.04108011388541695, 0.0, 0.0, 0.0,
    0.72213999695319386, -0.50610480451482776, -0.13788591219744095,
    -0.46457178052636866, 0.59705271243349645, -0.051332372890170661, 0.0, 0.0,
    0.0, 0.86703533388470988, -0.60765329439522986, -0.16555232839132558,
    -0.55729324354555687, 0.71621535492909738, -0.061577533948447609, 0.0, 0.0,
    0.0, 1.0120868419082774, -0.70931123526899487, -0.19324856400191229,
    -0.649950512714253, 0.83529549755241927, -0.071815601974369078, 0.0, 0.0,
    0.0, 1.1572946893482228, -0.81107874510459244, -0.22097465116916482,
    -0.74254363247578548, 0.95429319742059882, -0.082046581878654157, 0.0, 0.0,
    0.0, 1.3026590447102959, -0.912955941997641, -0.24873062206768792,
    -0.83507264724271346, 1.0732085116112287, -0.092270478568622111, 0.0, 0.0,
    0.0, 1.4481800766818653, -1.0149429441710447, -0.27651650890676488,
    -0.92753760139684793, 1.1920414971623847, -0.10248729694819469, 0.0, 0.0,
    0.0, 1.5938579541321145, -1.1170398699751309, -0.30433234393039471,
    -1.0199385392892728, 1.3107922110726542, -0.11269704191789855, 0.0, 0.0, 0.0,
    1.739692846112237, -1.2192468378877874, -0.33217815941732975,
    -1.1122755052403666, 1.4294607103011623, -0.12289971837486753, 0.0, 0.0, 0.0,
    1.8856849218556335, -1.3215639665145995, -0.36005398768111291,
    -1.2045485435398235, 1.5480470517675995, -0.13309533121284506, 0.0, 0.0, 0.0,
    2.0318343507781078, -1.4239913745889881, -0.38795986107011532,
    -1.2967576984466747, 1.6665512923522492, -0.14328388532218647, 0.0, 0.0, 0.0,
    2.1781413024780631, -1.5265291809723474, -0.41589581196757386,
    -1.3889030141893091, 1.7849734888960149, -0.15346538558986136, 0.0, 0.0, 0.0,
    2.3246059467366997, -1.6291775046541828, -0.4438618727916287,
    -1.4809845349654953, 1.9033136982004475, -0.16363983689945594, 0.0, 0.0, 0.0,
    2.4712284535182114, -1.7319364647522488, -0.47185807599536089,
    -1.5730023049424027, 2.0215719770277722, -0.17380724413117532, 0.0, 0.0, 0.0,
    2.6180089929699819, -1.8348061805126878, -0.49988445406683013,
    -1.6649563682566217, 2.1397483821009167, -0.18396761216184598, 0.0, 0.0, 0.0,
    2.7649477354227843, -1.9377867713101677, -0.52794103952911242,
    -1.7568467690141867, 2.2578429701035367, -0.19412094586491796, 0.0, 0.0, 0.0,
    2.9120448513909771, -2.0408783566480211, -0.55602786494033762,
    -1.848673551290595, 2.3758557976800447, -0.20426725011046729, 0.0, 0.0, 0.0
  };

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

  static const real_T g[16]{ 1.5936374091129624, -1.3562674366201879,
    -0.16194972057192442, 0.0, -1.3562674366201877, 1.3230276631452211,
    0.0844082916306713, 0.0, -0.16194972057192439, 0.084408291630671317,
    0.080052135759209747, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const real_T g_0[16]{ 2.7232488148062077, -1.0395805743008937,
    -1.0395805743008937, 0.0, -1.0395805743008937, 0.43600937538877005,
    0.40188791009141395, 0.0, -1.0395805743008937, 0.40188791009141395,
    0.43600937538877005, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const real_T g_1[16]{ 0.044658871119322277, -0.00953798324465595,
    -0.0500367818799705, 0.0, -0.0095379832446559482, 1.2385997808899454,
    -0.68758688074338681, 0.0, -0.0500367818799705, -0.6875868807433867,
    0.72086757880964547, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const real_T W_1[5]{ 0.0, 0.018316915599999997, 0.018316915599999997,
    0.0, 0.0 };

  static const real_T W_0[4]{ 0.018316915599999997, 0.0, 0.0,
    0.018316915599999997 };

  static const real_T W[3]{ 0.018316915599999997, 0.018316915599999997, 0.0 };

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

  static const int16_T b_Mlim[126]{ 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 80, 80, 80, 0, 0, 0 };

  static const int8_T b_A[400]{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
    1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T b_Jm[180]{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0 };

  static const int8_T d[180]{ 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0,
    1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
    0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1,
    0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
    1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0,
    0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0,
    0, 0, 1, 0, 0, 0, 1 };

  static const int8_T f[180]{ 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0 };

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

  static const uint8_T b_Mrows_0[126]{ 1U, 2U, 3U, 4U, 5U, 6U, 7U, 8U, 9U, 10U,
    11U, 12U, 13U, 14U, 15U, 16U, 17U, 18U, 19U, 20U, 21U, 22U, 23U, 24U, 25U,
    26U, 27U, 28U, 29U, 30U, 31U, 32U, 33U, 34U, 35U, 36U, 37U, 38U, 39U, 40U,
    41U, 42U, 43U, 44U, 45U, 46U, 47U, 48U, 49U, 50U, 51U, 52U, 53U, 54U, 55U,
    56U, 57U, 58U, 59U, 60U, 61U, 62U, 63U, 64U, 65U, 66U, 67U, 68U, 69U, 70U,
    71U, 72U, 73U, 74U, 75U, 76U, 77U, 78U, 79U, 80U, 81U, 82U, 83U, 84U, 85U,
    86U, 87U, 88U, 89U, 90U, 91U, 92U, 93U, 94U, 95U, 96U, 97U, 98U, 99U, 100U,
    101U, 102U, 103U, 104U, 105U, 106U, 107U, 108U, 109U, 110U, 111U, 112U, 113U,
    114U, 115U, 116U, 117U, 118U, 119U, 120U, 121U, 122U, 123U, 181U, 182U, 183U
  };

  real_T WySuJm_0[300];
  real_T WySuJm[240];
  real_T Bc_0[206];
  real_T a__1_0[206];
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T WySuJm_1[180];
  real_T Bc[166];
  real_T a__1[166];
  real_T Product1_j[144];
  real_T Bc_1[126];
  real_T a__1_1[126];
  real_T rseq_0[100];
  real_T rseq[80];
  real_T B_est_0[72];
  real_T rtb_A[64];
  real_T rtb_Product[64];
  real_T rtb_Q[64];
  real_T rtb_Z_e[64];
  real_T B_est[63];
  real_T rseq_1[60];
  real_T Abar[49];
  real_T rtb_A_e[49];
  real_T rtb_Q_j[49];
  real_T rtb_Transpose2_0[49];
  real_T rtb_Z[49];
  real_T rtb_y_g[49];
  real_T rtb_y_m[49];
  real_T y_0[48];
  real_T y[42];
  real_T b_Kx[30];
  real_T D_est[27];
  real_T rtb_R_tmp[27];
  real_T rtb_B[24];
  real_T rtb_C[24];
  real_T rtb_L[24];
  real_T rtb_N[24];
  real_T rtb_Product2[24];
  real_T rtb_Add_k[21];
  real_T rtb_B_o[21];
  real_T rtb_C_c[21];
  real_T rtb_N_f[21];
  real_T rtb_Product2_bg[21];
  real_T rtb_Transpose2[21];
  real_T b_L[16];
  real_T b_Linv[16];
  real_T Sum_h[12];
  real_T rtb_xest_p[10];
  real_T I2Jm_0[9];
  real_T rtb_R[9];
  real_T rtb_R_0[9];
  real_T rtb_y[9];
  real_T rtb_A_1[8];
  real_T rtb_xest[8];
  real_T rtb_A_j[7];
  real_T rtb_Sum2_f[7];
  real_T rtb_ywt[6];
  real_T rtb_ywtT[6];
  real_T rtb_TmpSignalConversionAtSFun_h[5];
  real_T rtb_TmpSignalConversionAtSFun_p[4];
  real_T zopt[4];
  real_T Sum2_c[3];
  real_T b_Wu[3];
  real_T rtb_C_1[3];
  real_T rtb_Product1_nb[3];
  real_T rtb_Sum6_c[3];
  real_T tmp[3];
  int32_T kidx;
  uint16_T waypt;
  int8_T a[3600];
  int8_T b_I[49];
  int8_T b[16];
  int8_T P0_2_tmp[12];
  boolean_T rtb_iAout_a[206];
  boolean_T rtb_iAout_l[166];
  boolean_T rtb_iAout_k[126];
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
        for (int32_T kidx_0{0}; kidx_0 < 6; kidx_0++) {
          // Inport: '<Root>/y'
          rtDW.traj[kidx_0] = rtU.y[kidx_0];
        }
      } else {
        // '<S1>:520:6' else
        //  hold last waypoint pos
        // '<S1>:520:7' traj(:,1) = traj(:, waypt);
        for (int32_T kidx_0{0}; kidx_0 < 6; kidx_0++) {
          rtDW.traj[kidx_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
            6 + kidx_0];
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
      __m128d tmp_2;
      real_T dwt;
      real_T r__h_idx_0;
      real_T r__h_idx_1;
      real_T r__h_idx_2;
      real_T s;
      int32_T Tries;
      int32_T i;
      int32_T kidx_0;
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
            for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
              // Inport: '<Root>/y'
              rtDW.traj[kidx_0] = rtU.y[kidx_0];
            }
          } else {
            // '<S1>:520:6' else
            //  hold last waypoint pos
            // '<S1>:520:7' traj(:,1) = traj(:, waypt);
            for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
              rtDW.traj[kidx_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) -
                1) * 6 + kidx_0];
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
              kidx = 0;
              exitg1 = false;
              while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                       (0U)) == false) && (kidx < 6)) {
                if (!(rtU.nextEv.r[kidx] == rtP.nullEv.r[kidx])) {
                  rstTheta2 = false;
                  exitg1 = true;
                } else {
                  kidx++;
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
        kidx = 3;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (kidx - 3 < 3)) {
          if (rtU.yo[kidx]) {
            rstP2 = true;
            exitg1 = true;
          } else {
            kidx++;
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
          for (kidx_0 = 0; kidx_0 < 12; kidx_0++) {
            P0_2_tmp[kidx_0] = static_cast<int8_T>(kidx_0 + 1);
          }

          // '<S1>:59:5' theta0_2 = theta(1:np*no);
          for (Tries = 0; Tries < 12; Tries++) {
            for (kidx_0 = 0; kidx_0 < 12; kidx_0++) {
              rtDW.P0_2[kidx_0 + 12 * Tries] = rtY.P_p[((P0_2_tmp[Tries] - 1) *
                24 + P0_2_tmp[kidx_0]) - 1];
            }

            rtDW.theta0_2[Tries] = rtY.theta[Tries];
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
        Tries = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (Tries < 3)) {
          if (rtU.yo[Tries]) {
            rstP1 = true;
            exitg1 = true;
          } else {
            Tries++;
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
          for (kidx_0 = 0; kidx_0 < 12; kidx_0++) {
            P0_2_tmp[kidx_0] = static_cast<int8_T>(kidx_0 + 13);
          }

          // '<S1>:59:11' theta0_1 = theta(np*no+1:2*np*no);
          for (Tries = 0; Tries < 12; Tries++) {
            for (kidx_0 = 0; kidx_0 < 12; kidx_0++) {
              rtDW.P0_1[kidx_0 + 12 * Tries] = rtY.P_p[((P0_2_tmp[Tries] - 1) *
                24 + P0_2_tmp[kidx_0]) - 1];
            }

            rtDW.theta0_1[Tries] = rtY.theta[Tries + 12];
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
      kidx = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (kidx < 6)) {
        if (!rtU.enAdapt[kidx]) {
          c_y = false;
          exitg1 = true;
        } else {
          kidx++;
        }
      }

      if (!c_y) {
        // '<S1>:59:16' enAdapt_(:) = false;
        for (Tries = 0; Tries < 6; Tries++) {
          rtDW.enAdapt_[Tries] = false;
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
      Tries = 0;
      kidx = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S324>/Product1'

        (void)std::memcpy(&rtY.P_p[Tries], &Product1_j[kidx], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i] = Sum_h[i];
        Tries += 24;
        kidx += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S324>/Sum2'

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
      Tries = 0;
      kidx = 0;
      for (i = 0; i < 12; i++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S328>/Product1'

        (void)std::memcpy(&rtY.P_p[Tries + 300], &Product1_j[kidx], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i + 12] = Sum_h[i];
        Tries += 24;
        kidx += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S328>/Sum2'

      rtY.prmErr[3] = Sum2_c[0];
      rtY.prmErr[4] = Sum2_c[1];
      rtY.prmErr[5] = Sum2_c[2];

      // '<S1>:59:28' currTraj = traj(:,waypt);
      for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
        // Outport: '<Root>/currTraj'
        rtY.currTraj[kidx_0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
          6 + kidx_0];
      }

      // Outport: '<Root>/sig' incorporates:
      //   Outport: '<Root>/ywt'

      // '<S1>:59:29' sig = gainSchSig(ywt);
      rtY.sig = gainSchSig(rtY.ywt);

      // Outputs for Function Call SubSystem: '<S1>/wtMod'
      // MATLAB Function: '<S9>/MATLAB Function' incorporates:
      //   Inport: '<Root>/k_2'

      // '<S1>:59:30' [ywt, y__, r_] = wtMod(currEv.r, y, k_2, currTraj);
      // Simulink Function 'wtMod': '<S1>:911'
      // MATLAB Function 'SupervisoryController/wtMod/MATLAB Function': '<S330>:1' 
      // '<S330>:1:2' [ywt, ywtT, uwt, uwtT] = wtMod_(y, yDest, ywtT, uwtT, dt, no, ni, k_2); 
      // 'wtMod_:3' ywt = zeros(1,2*no);
      // 'wtMod_:4' uwt = dt*ones(1,ni);
      //  time-scaled sigmoid (around 0->1 in 0->k_2 seconds, k_2 being a time scale factor) 
      // 'wtMod_:7' r = 0.2;
      //  10% to 90% rise time
      // 'wtMod_:8' dwt = dt/k_2;
      dwt = rtP.dt / rtU.k_2;

      // 'wtMod_:9' k_1 = 2.197/(r*k_2);
      s = 2.197 / (0.2 * rtU.k_2);

      // 'wtMod_:10' x0 = 0.5*k_2;
      r__h_idx_0 = 0.5 * rtU.k_2;

      //  midpoint
      // 'wtMod_:11' sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));
      // 'wtMod_:13' for i = 1:2*no
      for (kidx = 0; kidx < 6; kidx++) {
        // Delay: '<S9>/Delay'
        r__h_idx_1 = rtDW.Delay_DSTATE[kidx];

        // MATLAB Function: '<S9>/MATLAB Function' incorporates:
        //   Inport: '<Root>/k_2'
        //   Outport: '<Root>/currEv'

        // 'wtMod_:14' if yDest(i) ~= 0
        if (rtY.currEv.r[kidx] != 0.0) {
          //  drive ywt to 1
          // 'wtMod_:16' if (ywtT(i) <= 1)
          if (r__h_idx_1 <= 1.0) {
            // 'wtMod_:17' ywtT(i) = ywtT(i) + dwt;
            r__h_idx_1 += dwt;
          }

          // 'wtMod_:19' else
          //  drive ywt to 0
          // 'wtMod_:21' if (ywtT(i) > 0)
        } else if (r__h_idx_1 > 0.0) {
          // 'wtMod_:22' ywtT(i) = ywtT(i) - dwt;
          r__h_idx_1 -= dwt;
        } else {
          // no actions
        }

        // 'wtMod_:25' if ywtT(i) <= 0
        if (r__h_idx_1 <= 0.0) {
          // 'wtMod_:26' ywt(i) = 0;
          rtb_ywt[kidx] = 0.0;
        } else {
          // 'wtMod_:27' else
          // 'wtMod_:28' ywt(i) = sigmoid(ywtT(i)*k_2);
          // 'wtMod_:11' @(x) 1/(1 + exp(-k_1*(x-x0)))
          rtb_ywt[kidx] = 1.0 / (std::exp((r__h_idx_1 * rtU.k_2 - r__h_idx_0) *
            -s) + 1.0);
        }

        // Delay: '<S9>/Delay'
        rtb_ywtT[kidx] = r__h_idx_1;
      }

      // MATLAB Function: '<S9>/MATLAB Function' incorporates:
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
      // '<S330>:1:3' r_ = zeros(no, 1);
      // '<S330>:1:4' if any(yDest(1:no))
      if (any(&rtY.currEv.r[0])) {
        // '<S330>:1:5' r_ = r(1:no);
        // '<S330>:1:6' y_ = y(1:no);
        r__h_idx_0 = rtY.currTraj[0];
        Sum2_c[0] = rtU.y[0];
        r__h_idx_1 = rtY.currTraj[1];
        Sum2_c[1] = rtU.y[1];
        r__h_idx_2 = rtY.currTraj[2];
        Sum2_c[2] = rtU.y[2];
      } else if (any(&rtY.currEv.r[3])) {
        // '<S330>:1:7' elseif any(yDest(no+1:2*no))
        // '<S330>:1:8' r_ = r(no+1:2*no);
        // '<S330>:1:9' y_ = y(no+1:2*no);
        r__h_idx_0 = rtY.currTraj[3];
        Sum2_c[0] = rtU.y[3];
        r__h_idx_1 = rtY.currTraj[4];
        Sum2_c[1] = rtU.y[4];
        r__h_idx_2 = rtY.currTraj[5];
        Sum2_c[2] = rtU.y[5];
      } else {
        // '<S330>:1:10' else
        // '<S330>:1:11' r_ = zeros(no, 1);
        // '<S330>:1:12' y_ = y(1:no);
        r__h_idx_0 = 0.0;
        Sum2_c[0] = rtU.y[0];
        r__h_idx_1 = 0.0;
        Sum2_c[1] = rtU.y[1];
        r__h_idx_2 = 0.0;
        Sum2_c[2] = rtU.y[2];
      }

      // End of Outputs for SubSystem: '<S1>/wtMod'
      for (Tries = 0; Tries <= 4; Tries += 2) {
        // Outputs for Function Call SubSystem: '<S1>/wtMod'
        // Update for Delay: '<S9>/Delay'
        tmp_2 = _mm_loadu_pd(&rtb_ywtT[Tries]);
        (void)_mm_storeu_pd(&rtDW.Delay_DSTATE[Tries], tmp_2);

        // Gain: '<S9>/Gain' incorporates:
        //   Delay: '<S9>/Delay'

        tmp_2 = _mm_loadu_pd(&rtb_ywt[Tries]);

        // Outport: '<Root>/ywt' incorporates:
        //   Delay: '<S9>/Delay'
        //   Gain: '<S9>/Gain'

        (void)_mm_storeu_pd(&rtY.ywt[Tries], _mm_mul_pd(_mm_set1_pd(rtP.beta),
          tmp_2));

        // End of Outputs for SubSystem: '<S1>/wtMod'
      }

      // Outport: '<Root>/sig' incorporates:
      //   Constant: '<S183>/G'
      //   Constant: '<S183>/H'
      //   Constant: '<S253>/G'
      //   Constant: '<S253>/H'
      //   Constant: '<S3>/Constant1'
      //   Constant: '<S3>/Constant13'
      //   Constant: '<S4>/Constant1'
      //   Constant: '<S4>/Constant13'
      //   Constant: '<S5>/Constant1'
      //   Constant: '<S5>/Constant13'
      //   DataTypeConversion: '<S113>/DataTypeConversionEnable'
      //   DataTypeConversion: '<S183>/DataTypeConversionEnable'
      //   DataTypeConversion: '<S253>/DataTypeConversionEnable'
      //   Delay: '<S183>/MemoryP'
      //   Delay: '<S183>/MemoryX'
      //   Delay: '<S253>/MemoryP'
      //   Delay: '<S253>/MemoryX'
      //   Outport: '<Root>/yhat'
      //   Product: '<S116>/Product'
      //   Product: '<S116>/Product1'
      //   Product: '<S156>/C[k]*xhat[k|k-1]'
      //   Product: '<S156>/D[k]*u[k]'
      //   Product: '<S156>/Product3'
      //   Product: '<S186>/Product'
      //   Product: '<S186>/Product1'
      //   Product: '<S256>/Product'
      //   Product: '<S256>/Product1'
      //   Sum: '<S160>/Sum3'
      //   Sum: '<S230>/Sum3'
      //   Sum: '<S90>/Sum3'

      // [u, yhat] = ampc(currTraj, y, ymax, ywt, y0, x0, u0, umax, uwt, excitation, theta, thetaSgn); 
      // [u, ywt, currTraj] = gmpc(traj(:,waypt), currEv.r, y, ymax, umax, uwt, k_2); 
      // '<S1>:59:33' if sig == 1
      if (rtY.sig == 1.0) {
        __m128d tmp_0;
        __m128d tmp_1;
        real_T Saturation_l_idx_0;
        real_T Saturation_l_idx_1;
        real_T Saturation_l_idx_2;
        int32_T I2Jm_tmp;
        int32_T a_tmp;
        int32_T b_Linv_tmp;
        int32_T rtb_R_tmp_0;
        int16_T ixw;
        boolean_T guard11;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'

        // '<S1>:59:34' [u, yhat(1:no)] = mpc1(r_, y__, [0;0;0], 0, u0, umax, uwt, iRST); 
        // Simulink Function 'mpc1': '<S1>:882'
        if (rtU.iRST && (rtDW.DiscreteTimeIntegrator_PrevRe_b <= 0)) {
          rtDW.DiscreteTimeIntegrator_DSTATE_j = rtP.DiscreteTimeIntegrator_IC;
        }

        // Gain: '<S3>/Gain1' incorporates:
        //   Inport: '<Root>/uwt'
        //   Product: '<S116>/Product1'

        rtb_Product1_nb[0] = rtP.beta * rtU.uwt[0];
        rtb_Product1_nb[1] = rtP.beta * rtU.uwt[1];
        rtb_Product1_nb[2] = rtP.beta * rtU.uwt[2];

        // Delay: '<S113>/MemoryX' incorporates:
        //   Constant: '<S113>/X0'
        //   DataTypeConversion: '<S113>/DataTypeConversionReset'

        rtDW.icLoad_n = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE_l) ==
                          POS_ZCSIG) || rtDW.icLoad_n);
        rtPrevZCX.MemoryX_Reset_ZCE_l = 0U;
        if (rtDW.icLoad_n) {
          for (Tries = 0; Tries < 7; Tries++) {
            rtDW.MemoryX_DSTATE_l[Tries] = rtP.X0_Value_f[Tries];
          }
        }

        // Sum: '<S90>/Sum2' incorporates:
        //   Delay: '<S113>/MemoryX'

        rtb_Sum2_f[0] = rtDW.MemoryX_DSTATE_l[0];

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (Tries = 0; Tries <= 4; Tries += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_2 = _mm_loadu_pd(&rtDW.MemoryX_DSTATE_l[Tries + 1]);
          (void)_mm_storeu_pd(&rtb_Sum2_f[Tries + 1], _mm_add_pd(tmp_2,
            _mm_loadu_pd(&rtP.Constant1_Value_j[Tries])));

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // MATLAB Function: '<S111>/optimizer' incorporates:
        //   Constant: '<S3>/Constant'
        //   Constant: '<S90>/Constant1'
        //   Delay: '<S113>/MemoryX'
        //   DiscreteIntegrator: '<S3>/Discrete-Time Integrator'
        //   Gain: '<S3>/Gain2'
        //   Math: '<S91>/Math Function1'
        //   Product: '<S116>/Product1'
        //   SignalConversion generated from: '<S112>/ SFunction '
        //   Sum: '<S90>/Sum2'
        //   UnitDelay: '<S91>/last_mv'

        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S112>:1'
        // '<S112>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S112>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S112>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S112>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S112>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S112>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S112>:1:252' if isa(u,'double')
        // '<S112>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S112>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S112>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S112>:1:252' if isa(u,'double')
        // '<S112>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S112>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S112>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S112>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S112>:1:252' if isa(u,'double')
        // '<S112>:1:253' y = u;
        //    md:         measured disturbance
        // '<S112>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S112>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S112>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S112>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S112>:1:252' if isa(u,'double')
        // '<S112>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S112>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S112>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S112>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S112>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S112>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S112>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S112>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S112>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S112>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S112>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S112>:1:58' uwt = convertDataType(uwt0,isDouble);
        // '<S112>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S112>:1:252' if isa(u,'double')
        // '<S112>:1:253' y = u;
        //    duwt:       run-time DMV weights
        // '<S112>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S112>:1:62' ewt = convertDataType(ewt0,isDouble);
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
        // '<S112>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S112>:1:96' isAdaptive = false;
        // '<S112>:1:97' isLTV = false;
        // '<S112>:1:98' ZERO = zeros('like',ref);
        // '<S112>:1:99' ONE = ones('like',ref);
        // '<S112>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S112>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S112>:1:119' if isSimulation
        // '<S112>:1:126' else
        //  When doing code generation, use M code directly
        // '<S112>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        kidx_0 = 0;
        for (kidx = 0; kidx < 20; kidx++) {
          rseq[kidx_0] = r__h_idx_0;
          rseq[kidx_0 + 1] = r__h_idx_1;
          rseq[kidx_0 + 2] = r__h_idx_2;
          rseq[kidx_0 + 3] = rtP.Constant_Value;
          kidx_0 += 4;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S112>:1:133' old_u = old_u - uoff;
        Saturation_l_idx_0 = rtDW.last_mv_DSTATE_n[0];
        Saturation_l_idx_1 = rtDW.last_mv_DSTATE_n[1];
        Saturation_l_idx_2 = rtDW.last_mv_DSTATE_n[2];

        // '<S112>:1:134' if no_mv
        // '<S112>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S112>:1:143' xk = xk - xoff;
        rtb_xest[0] = rtb_Sum2_f[0];
        rtb_xest[1] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE_j;
        for (Tries = 0; Tries < 6; Tries++) {
          rtb_xest[Tries + 2] = rtb_Sum2_f[Tries + 1];
        }

        //  Remove offset
        // '<S112>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S112>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S112>:1:162' if no_uref
        //  no external utarget
        // '<S112>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S112>:1:173' if ~no_cc
        // '<S112>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S112>:1:183' if isSimulation
        // '<S112>:1:214' else
        //  When doing code generation, use M code directly
        // '<S112>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S112>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S112>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S112>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S112>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S112>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S112>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S112>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S112>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S112>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S112>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S112>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S112>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        (void)std::memcpy(&b_Linv[0], &g_0[0], sizeof(real_T) << 4UL);
        (void)std::memset(&rtb_iAout_l[0], 0, 166U * sizeof(boolean_T));
        if (rtb_Product1_nb[0] < 0.0) {
          b_Wu[0] = 0.0;
        } else {
          b_Wu[0] = rtb_Product1_nb[0] * rtb_Product1_nb[0];
        }

        if (rtb_Product1_nb[1] < 0.0) {
          b_Wu[1] = 0.0;
        } else {
          b_Wu[1] = rtb_Product1_nb[1] * rtb_Product1_nb[1];
        }

        if (rtb_Product1_nb[2] < 0.0) {
          b_Wu[2] = 0.0;
        } else {
          b_Wu[2] = rtb_Product1_nb[2] * rtb_Product1_nb[2];
        }

        (void)std::memset(&rtb_y[0], 0, 9U * sizeof(real_T));
        rtb_y[0] = 1.0;
        rtb_y[4] = 1.0;
        rtb_y[8] = 1.0;
        kidx = -1;
        for (Tries = 0; Tries < 20; Tries++) {
          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            for (i = 0; i < 20; i++) {
              a_tmp = static_cast<int32_T>(b_A[20 * Tries + i]);
              a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0]) * a_tmp);
              a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 1]) * a_tmp);
              a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 2]) * a_tmp);
              kidx += 3;
            }
          }
        }

        kidx_0 = 0;
        for (Tries = 0; Tries < 3; Tries++) {
          for (kidx = 0; kidx < 60; kidx++) {
            I2Jm_tmp = kidx + kidx_0;
            I2Jm[I2Jm_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 60; a_tmp++) {
              I2Jm[I2Jm_tmp] += static_cast<real_T>(static_cast<int32_T>(
                static_cast<int32_T>(a[i + kidx]) * static_cast<int32_T>
                (b_Jm[a_tmp + kidx_0])));
              i += 60;
            }
          }

          kidx_0 += 60;
        }

        ixw = 1;
        for (kidx = 0; kidx < 80; kidx++) {
          dwt = W_0[ixw - 1];
          WySuJm[kidx] = dwt * b_SuJm[kidx];
          WySuJm[kidx + 80] = b_SuJm[kidx + 80] * dwt;
          WySuJm[kidx + 160] = b_SuJm[kidx + 160] * dwt;
          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 4) {
            ixw = 1;
          }
        }

        ixw = 1;
        for (kidx = 0; kidx < 60; kidx++) {
          dwt = b_Wu[ixw - 1];
          WuI2Jm[kidx] = dwt * I2Jm[kidx];
          WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * dwt;
          WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * dwt;
          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 3) {
            ixw = 1;
          }

          WduJm[kidx] = 0.034121465297356074 * static_cast<real_T>(b_Jm[kidx]);
          WduJm[kidx + 60] = static_cast<real_T>(b_Jm[kidx + 60]) *
            0.034121465297356074;
          WduJm[kidx + 120] = static_cast<real_T>(b_Jm[kidx + 120]) *
            0.034121465297356074;
        }

        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 3; Tries++) {
            rtb_R_tmp_0 = 3 * Tries + kidx_0;
            rtb_R[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 80; kidx++) {
              rtb_R[rtb_R_tmp_0] += e[3 * kidx + kidx_0] * WySuJm[80 * Tries +
                kidx];
            }

            s = 0.0;
            I2Jm_0[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 60; kidx++) {
              i = 60 * Tries + kidx;
              s += static_cast<real_T>(f[3 * kidx + kidx_0]) * WduJm[i];
              I2Jm_0[rtb_R_tmp_0] += I2Jm[60 * kidx_0 + kidx] * WuI2Jm[i];
            }

            rtb_R_0[rtb_R_tmp_0] = rtb_R[rtb_R_tmp_0] + s;
          }
        }

        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          i = 0;
          a_tmp = 0;
          for (I2Jm_tmp = 0; I2Jm_tmp < 3; I2Jm_tmp++) {
            b_Linv_tmp = I2Jm_tmp + Tries;
            b_Linv[I2Jm_tmp + kidx_0] = rtb_R_0[b_Linv_tmp] + I2Jm_0[b_Linv_tmp];
            s = 0.0;
            b_Linv_tmp = 0;
            for (rtb_R_tmp_0 = 0; rtb_R_tmp_0 < 60; rtb_R_tmp_0++) {
              s += static_cast<real_T>(d[b_Linv_tmp + kidx]) *
                WuI2Jm[rtb_R_tmp_0 + a_tmp];
              b_Linv_tmp += 3;
            }

            b_Linv_tmp = i + kidx;
            rtb_y[b_Linv_tmp] = rtb_R[b_Linv_tmp] + s;
            i += 3;
            a_tmp += 60;
          }

          kidx_0 += 4;
          Tries += 3;
        }

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (kidx_0 = 0; kidx_0 <= 178; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_2 = _mm_loadu_pd(&WuI2Jm[kidx_0]);
          (void)_mm_storeu_pd(&WuI2Jm[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 8; Tries++) {
            // MATLAB Function: '<S111>/optimizer'
            i = (kidx_0 << 3UL) + Tries;
            rtb_B[i] = 0.0;
            for (kidx = 0; kidx < 80; kidx++) {
              rtb_B[i] += c_0[(kidx << 3UL) + Tries] * WySuJm[80 * kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 21; Tries++) {
            // MATLAB Function: '<S111>/optimizer'
            i = 21 * kidx_0 + Tries;
            B_est[i] = 0.0;
            for (kidx = 0; kidx < 80; kidx++) {
              B_est[i] += WySuJm[80 * kidx_0 + kidx] * 0.0;
            }
          }
        }

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (kidx_0 = 0; kidx_0 <= 238; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_2 = _mm_loadu_pd(&WySuJm[kidx_0]);
          (void)_mm_storeu_pd(&WySuJm[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // MATLAB Function: '<S111>/optimizer' incorporates:
        //   Inport: '<Root>/umax'
        //   Memory: '<S91>/Memory'
        //   UnitDelay: '<S91>/last_mv'

        kidx = 0;
        (void)std::memcpy(&b_L[0], &b_Linv[0], sizeof(real_T) << 4UL);
        Tries = xpotrf(b_L);
        guard11 = false;
        if (Tries == 0) {
          rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
          rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
          rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
          rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
          if (minimum(rtb_TmpSignalConversionAtSFun_p) > 1.4901161193847656E-7)
          {
          } else {
            guard11 = true;
          }
        } else {
          guard11 = true;
        }

        if (guard11) {
          boolean_T exitg2;
          dwt = 0.0;
          Tries = 0;
          exitg2 = false;
          while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (Tries < 4)) {
            s = ((std::abs(b_Linv[Tries + 4]) + std::abs(b_Linv[Tries])) + std::
                 abs(b_Linv[Tries + 8])) + std::abs(b_Linv[Tries + 12]);
            if (std::isnan(s)) {
              dwt = (rtNaN);
              exitg2 = true;
            } else {
              if (s > dwt) {
                dwt = s;
              }

              Tries++;
            }
          }

          if (dwt >= 1.0E+10) {
            kidx = 2;
          } else {
            Tries = 0;
            exitg1 = false;
            while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                     (0U)) == false) && (Tries <= 4)) {
              boolean_T guard2;
              dwt = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
                1.4901161193847656E-7;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b[kidx_0] = 0;
              }

              b[0] = 1;
              b[5] = 1;
              b[10] = 1;
              b[15] = 1;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b_Linv[kidx_0] += dwt * static_cast<real_T>(b[kidx_0]);
                b_L[kidx_0] = b_Linv[kidx_0];
              }

              kidx = xpotrf(b_L);
              guard2 = false;
              if (kidx == 0) {
                rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
                rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
                rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
                rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
                if (minimum(rtb_TmpSignalConversionAtSFun_p) >
                    1.4901161193847656E-7) {
                  kidx = 1;
                  exitg1 = true;
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                kidx = 3;
                Tries++;
              }
            }
          }
        }

        if (kidx > 1) {
          b_Wu[0] = rtDW.last_mv_DSTATE_n[0];
          b_Wu[1] = rtDW.last_mv_DSTATE_n[1];
          b_Wu[2] = rtDW.last_mv_DSTATE_n[2];
        } else {
          for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
            b[kidx_0] = 0;
          }

          b[0] = 1;
          b[5] = 1;
          b[10] = 1;
          b[15] = 1;
          kidx_0 = 0;
          for (kidx = 0; kidx < 4; kidx++) {
            b_Linv[kidx_0] = static_cast<real_T>(b[kidx_0]);
            b_Linv[kidx_0 + 1] = static_cast<real_T>(b[kidx_0 + 1]);
            b_Linv[kidx_0 + 2] = static_cast<real_T>(b[kidx_0 + 2]);
            b_Linv[kidx_0 + 3] = static_cast<real_T>(b[kidx_0 + 3]);
            kidx_0 += 4;
          }

          trisolve(b_L, b_Linv);
          umax_incr_flag[0] = false;
          b_Wu[0] = 0.0;
          umax_incr_flag[1] = false;
          b_Wu[1] = 0.0;
          umax_incr_flag[2] = false;
          b_Wu[2] = 0.0;
          for (kidx = 0; kidx < 166; kidx++) {
            uint8_T b_Mrows;
            dwt = b_Mlim_0[kidx];
            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
              s += b_Mx_0[166 * kidx_0 + kidx] * rtb_xest[kidx_0];
            }

            s = -(((b_Mu1_0[kidx + 166] * Saturation_l_idx_1 + b_Mu1_0[kidx] *
                    Saturation_l_idx_0) + b_Mu1_0[kidx + 332] *
                   Saturation_l_idx_2) + (dwt + s));
            b_Mrows = b_Mrows_1[kidx];
            if ((b_Mrows > 80UL) && (b_Mrows > 160UL) && (b_Mrows <= 220UL)) {
              Tries = (static_cast<int32_T>(b_Mrows) - div_nde_s32_floor(
                        static_cast<int32_T>(b_Mrows) - 161, static_cast<int32_T>
                        (nu)) * static_cast<int32_T>(nu)) - 161;
              rstP2 = umax_incr_flag[Tries];
              if (!umax_incr_flag[Tries]) {
                dwt = -rtU.umax[Tries] - (-dwt);
                rstP2 = true;
              } else {
                dwt = b_Wu[Tries];
              }

              b_Wu[Tries] = dwt;
              umax_incr_flag[Tries] = rstP2;
              s += dwt;
            }

            Bc[kidx] = s;
          }

          rtb_TmpSignalConversionAtSFun_p[0] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[1] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[2] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[3] = 0.0;
          for (kidx = 0; kidx < 3; kidx++) {
            real_T WuI2Jm_0;
            r__h_idx_1 = 0.0;
            for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
              r__h_idx_1 += rtb_B[(kidx << 3UL) + kidx_0] * rtb_xest[kidx_0];
            }

            dwt = 0.0;
            for (kidx_0 = 0; kidx_0 < 80; kidx_0++) {
              dwt += WySuJm[80 * kidx + kidx_0] * rseq[kidx_0];
            }

            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 21; kidx_0++) {
              s += B_est[21 * kidx + kidx_0];
            }

            WuI2Jm_0 = 0.0;
            for (kidx_0 = 0; kidx_0 < 60; kidx_0++) {
              WuI2Jm_0 += WuI2Jm[60 * kidx + kidx_0] * 0.0;
            }

            rtb_TmpSignalConversionAtSFun_p[kidx] = ((((rtb_y[3 * kidx + 1] *
              Saturation_l_idx_1 + rtb_y[3 * kidx] * Saturation_l_idx_0) +
              rtb_y[3 * kidx + 2] * Saturation_l_idx_2) + (r__h_idx_1 + dwt)) +
              s) + WuI2Jm_0;
          }

          (void)std::memcpy(&rtb_iAout_l[0], &rtDW.Memory_PreviousInput_d[0],
                            166U * sizeof(boolean_T));
          kidx_0 = 0;
          for (Tries = 0; Tries < 4; Tries++) {
            kidx = 0;
            for (i = 0; i < 4; i++) {
              b_Linv_tmp = kidx + Tries;
              b_L[b_Linv_tmp] = 0.0;
              b_L[b_Linv_tmp] += b_Linv[kidx_0] * b_Linv[kidx];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 1] * b_Linv[kidx + 1];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 2] * b_Linv[kidx + 2];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 3] * b_Linv[kidx + 3];
              kidx += 4;
            }

            kidx_0 += 4;
          }

          qpkwik(b_Linv, b_L, rtb_TmpSignalConversionAtSFun_p, b_Ac_0, Bc,
                 rtb_iAout_l, 680, 1.0E-6, zopt, a__1, &kidx);
          if ((kidx < 0) || (kidx == 0)) {
            zopt[0] = 0.0;
            zopt[1] = 0.0;
            zopt[2] = 0.0;
          }

          b_Wu[0] = rtDW.last_mv_DSTATE_n[0] + zopt[0];
          b_Wu[1] = rtDW.last_mv_DSTATE_n[1] + zopt[1];
          b_Wu[2] = rtDW.last_mv_DSTATE_n[2] + zopt[2];
        }

        // Delay: '<S113>/MemoryP' incorporates:
        //   Constant: '<S113>/P0'
        //   DataTypeConversion: '<S113>/DataTypeConversionReset'

        // '<S112>:1:231' if return_xseq || return_ovseq
        // '<S112>:1:233' else
        // '<S112>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S112>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S112>:1:238' if CustomEstimation
        // '<S112>:1:239' xk1 = xk;
        // '<S112>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S112>:1:247' xest = xest + xoff;
        rtDW.icLoad_h = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE_n) ==
                          POS_ZCSIG) || rtDW.icLoad_h);
        rtPrevZCX.MemoryP_Reset_ZCE_n = 0U;
        if (rtDW.icLoad_h) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE_e[0], &rtP.P0_Value_a[0], 49U *
                            sizeof(real_T));
        }

        // MATLAB Function: '<S90>/MATLAB Function' incorporates:
        //   BusCreator: '<S3>/Bus Creator1'
        //   Constant: '<S3>/Constant12'
        //   Constant: '<S3>/Constant13'
        //   Constant: '<S3>/Constant3'
        //   Constant: '<S3>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc1/State Estimator OD (KF)/MATLAB Function': '<S114>:1' 
        // '<S114>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod1, Bod1, Cod1(1:3,:), Dod1(1:3,:), Dmn1, 3, 3, 1); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&rtb_B_o[0], 0, 21U * sizeof(real_T));

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
        (void)std::memset(&rtb_A_e[0], 0, 49U * sizeof(real_T));
        rtb_A_e[0] = rtP.Constant3_Value;

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        rtb_B_o[0] = rtP.Constant4_Value[0];
        rtb_C_c[0] = rtP.Constant12_Value_e[0];
        rtb_B_o[7] = rtP.Constant4_Value[1];
        rtb_C_c[1] = rtP.Constant12_Value_e[1];
        rtb_B_o[14] = rtP.Constant4_Value[2];
        rtb_C_c[2] = rtP.Constant12_Value_e[2];
        kidx_0 = 0;
        Tries = 0;
        kidx = 0;
        i = 0;
        for (a_tmp = 0; a_tmp < 6; a_tmp++) {
          for (I2Jm_tmp = 0; I2Jm_tmp < 6; I2Jm_tmp++) {
            rtb_A_e[(I2Jm_tmp + kidx) + 8] = rtP.Aod1[I2Jm_tmp + i];
          }

          rtb_C_c[kidx_0 + 3] = rtP.Cod1[Tries];
          rtb_C_c[kidx_0 + 4] = rtP.Cod1[Tries + 1];
          rtb_C_c[kidx_0 + 5] = rtP.Cod1[Tries + 2];
          kidx_0 += 3;
          Tries += 4;
          kidx += 7;
          i += 6;
        }

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est[0], 0, 63U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y[0], 0, 42U * sizeof(real_T));
        y[0] = rtP.Constant4_Value[0];
        y[7] = rtP.Constant4_Value[1];
        y[14] = rtP.Constant4_Value[2];
        for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
          y[kidx_0 + 22] = rtP.Bod1[kidx_0];
          y[kidx_0 + 29] = rtP.Bod1[kidx_0 + 6];
          y[kidx_0 + 36] = rtP.Bod1[kidx_0 + 12];
        }

        (void)std::memcpy(&B_est[0], &y[0], 42U * sizeof(real_T));

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          D_est[kidx_0] = rtP.Constant13_Value_c[kidx_0];
          D_est[kidx_0 + 9] = rtP.Dod1[Tries];
          D_est[kidx_0 + 18] = rtP.Dmn1[kidx_0];
          D_est[kidx_0 + 1] = rtP.Constant13_Value_c[kidx_0 + 1];
          D_est[kidx_0 + 10] = rtP.Dod1[Tries + 1];
          D_est[kidx_0 + 19] = rtP.Dmn1[kidx_0 + 1];
          D_est[kidx_0 + 2] = rtP.Constant13_Value_c[kidx_0 + 2];
          D_est[kidx_0 + 11] = rtP.Dod1[Tries + 2];
          D_est[kidx_0 + 20] = rtP.Dmn1[kidx_0 + 2];
          kidx_0 += 3;
          Tries += 4;
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
          Tries = 0;
          for (kidx = 0; kidx < 7; kidx++) {
            b_Linv_tmp = Tries + kidx_0;
            rtb_Q_j[b_Linv_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 9; a_tmp++) {
              rtb_Q_j[b_Linv_tmp] += B_est[i + kidx_0] * B_est[i + kidx];
              i += 7;
            }

            Tries += 7;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        kidx_0 = 0;
        for (Tries = 0; Tries < 9; Tries++) {
          rtb_R_tmp[Tries] = D_est[kidx_0];
          rtb_R_tmp[Tries + 9] = D_est[kidx_0 + 1];
          rtb_R_tmp[Tries + 18] = D_est[kidx_0 + 2];
          kidx_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 3; Tries++) {
            rtb_R_tmp_0 = 3 * kidx_0 + Tries;
            rtb_R[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_R[rtb_R_tmp_0] += D_est[3 * kidx + Tries] * rtb_R_tmp[9 *
                kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 7; Tries++) {
            i = 7 * kidx_0 + Tries;
            rtb_N_f[i] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_N_f[i] += B_est[7 * kidx + Tries] * rtb_R_tmp[9 * kidx_0 +
                kidx];
            }
          }
        }

        // End of MATLAB Function: '<S90>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S113>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S113>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S113>/ScalarExpansionQ'
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
        kidx_0 = 0;
        for (Tries = 0; Tries < 7; Tries++) {
          kidx = 0;
          for (i = 0; i < 7; i++) {
            I2Jm_tmp = i + kidx_0;
            rtb_y_g[I2Jm_tmp] = (rtb_Q_j[kidx + Tries] + rtb_Q_j[I2Jm_tmp]) /
              2.0;
            kidx += 7;
          }

          kidx_0 += 7;
        }

        // End of MATLAB Function: '<S135>/ScalarExpansion'
        // End of Outputs for SubSystem: '<S113>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S113>/ReducedQRN'
        for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
          // Product: '<S133>/Product' incorporates:
          //   Constant: '<S113>/G'
          //   Math: '<S133>/Transpose1'

          Tries = 0;
          for (kidx = 0; kidx < 7; kidx++) {
            I2Jm_tmp = Tries + kidx_0;
            rtb_y_m[I2Jm_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 7; a_tmp++) {
              rtb_y_m[I2Jm_tmp] += rtb_y_g[i + kidx_0] * rtP.G_Value_a[i + kidx];
              i += 7;
            }

            Tries += 7;
          }
        }

        // Product: '<S133>/Product' incorporates:
        //   Constant: '<S113>/G'
        //   Constant: '<S113>/H'
        //   Math: '<S133>/Transpose2'

        kidx_0 = 0;
        for (Tries = 0; Tries < 7; Tries++) {
          kidx = 0;
          for (i = 0; i < 7; i++) {
            b_Linv_tmp = kidx + Tries;
            rtb_Q_j[b_Linv_tmp] = 0.0;
            a_tmp = 0;
            for (I2Jm_tmp = 0; I2Jm_tmp < 7; I2Jm_tmp++) {
              rtb_Q_j[b_Linv_tmp] += rtP.G_Value_a[a_tmp + Tries] *
                rtb_y_m[I2Jm_tmp + kidx];
              a_tmp += 7;
            }

            kidx += 7;
          }

          rtb_Transpose2[Tries] = rtP.H_Value_o[kidx_0];
          rtb_Transpose2[Tries + 7] = rtP.H_Value_o[kidx_0 + 1];
          rtb_Transpose2[Tries + 14] = rtP.H_Value_o[kidx_0 + 2];
          kidx_0 += 3;
        }

        for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
          // Sum: '<S133>/Add' incorporates:
          //   Math: '<S133>/Transpose2'
          //   Product: '<S133>/Product1'

          Tries = 0;
          for (kidx = 0; kidx < 3; kidx++) {
            s = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 7; a_tmp++) {
              s += rtb_y_g[i + kidx_0] * rtb_Transpose2[a_tmp + Tries];
              i += 7;
            }

            i = Tries + kidx_0;
            rtb_Add_k[i] = rtb_N_f[i] + s;
            Tries += 7;
          }

          // End of Sum: '<S133>/Add'
        }

        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 7; Tries++) {
            // Product: '<S133>/Product2' incorporates:
            //   Constant: '<S113>/G'
            //   Sum: '<S133>/Add'

            i = 7 * kidx_0 + Tries;
            rtb_Product2_bg[i] = 0.0;
            for (kidx = 0; kidx < 7; kidx++) {
              rtb_Product2_bg[i] += rtP.G_Value_a[7 * kidx + Tries] * rtb_Add_k
                [7 * kidx_0 + kidx];
            }

            // End of Product: '<S133>/Product2'
          }

          for (Tries = 0; Tries < 3; Tries++) {
            // Product: '<S133>/Product3' incorporates:
            //   Constant: '<S113>/H'
            //   Product: '<S133>/Product4'

            kidx = 3 * Tries + kidx_0;
            rtb_R_0[kidx] = 0.0;

            // Product: '<S133>/Product4'
            I2Jm_0[kidx] = 0.0;
            for (i = 0; i < 7; i++) {
              // Product: '<S133>/Product3' incorporates:
              //   Constant: '<S113>/H'
              //   Product: '<S133>/Product4'

              a_tmp = 7 * Tries + i;
              rtb_R_0[kidx] += rtP.H_Value_o[3 * i + kidx_0] * rtb_Add_k[a_tmp];
              I2Jm_0[kidx] += rtb_N_f[7 * kidx_0 + i] * rtb_Transpose2[a_tmp];
            }
          }
        }

        // End of Outputs for SubSystem: '<S113>/ReducedQRN'
        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (kidx_0 = 0; kidx_0 <= 6; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          // Outputs for Atomic SubSystem: '<S113>/ReducedQRN'
          tmp_2 = _mm_loadu_pd(&rtb_R_0[kidx_0]);
          tmp_0 = _mm_loadu_pd(&I2Jm_0[kidx_0]);
          tmp_1 = _mm_loadu_pd(&rtb_y[kidx_0]);
          (void)_mm_storeu_pd(&rtb_R[kidx_0], _mm_add_pd(_mm_add_pd(tmp_2, tmp_0),
            tmp_1));

          // End of Outputs for SubSystem: '<S113>/ReducedQRN'
          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Outputs for Atomic SubSystem: '<S113>/ReducedQRN'
        for (kidx_0 = 8; kidx_0 < 9; kidx_0++) {
          // Sum: '<S133>/Add1'
          rtb_R[kidx_0] = (rtb_R_0[kidx_0] + I2Jm_0[kidx_0]) + rtb_y[kidx_0];
        }

        // End of Outputs for SubSystem: '<S113>/ReducedQRN'

        // Outputs for Atomic SubSystem: '<S113>/CalculatePL'
        // MATLAB Function: '<S115>/Discrete-Time KF - Calculate PLMZ' incorporates:
        //   Constant: '<S113>/H'
        //   Constant: '<S3>/Constant1'
        //   DataTypeConversion: '<S113>/DataTypeConversionEnable'
        //   Delay: '<S113>/MemoryP'
        //   Math: '<S133>/Transpose'
        //   Math: '<S133>/Transpose2'
        //   Product: '<S133>/Product'
        //   Product: '<S133>/Product2'
        //   Product: '<S133>/Product3'
        //   Product: '<S133>/Product4'
        //   Sum: '<S133>/Add'
        //   Sum: '<S133>/Add1'

        //  See help of ctrlKalmanFilterDTCalculatePL.m
        // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S153>:1' 
        //    Copyright 2014 The MathWorks, Inc.
        // '<S153>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
        if (rtP.Constant1_Value_e != 0.0) {
          kidx_0 = 0;
          for (Tries = 0; Tries < 7; Tries++) {
            rtb_Transpose2[Tries] = rtb_C_c[kidx_0];
            rtb_Transpose2[Tries + 7] = rtb_C_c[kidx_0 + 1];
            rtb_Transpose2[Tries + 14] = rtb_C_c[kidx_0 + 2];
            kidx_0 += 3;
          }

          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            for (Tries = 0; Tries < 7; Tries++) {
              i = 3 * Tries + kidx_0;
              rtb_N_f[i] = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                rtb_N_f[i] += rtb_C_c[3 * kidx + kidx_0] *
                  rtDW.MemoryP_DSTATE_e[7 * Tries + kidx];
              }
            }

            for (Tries = 0; Tries < 3; Tries++) {
              s = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                s += rtb_N_f[3 * kidx + kidx_0] * rtb_Transpose2[7 * Tries +
                  kidx];
              }

              b_Linv_tmp = 3 * Tries + kidx_0;
              rtb_y[b_Linv_tmp] = rtb_R[b_Linv_tmp] + s;
            }
          }

          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            for (Tries = 0; Tries < 7; Tries++) {
              i = 7 * Tries + kidx_0;
              rtb_y_m[i] = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                rtb_y_m[i] += rtb_A_e[7 * kidx + kidx_0] *
                  rtDW.MemoryP_DSTATE_e[7 * Tries + kidx];
              }
            }

            for (Tries = 0; Tries < 3; Tries++) {
              s = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                s += rtb_y_m[7 * kidx + kidx_0] * rtb_Transpose2[7 * Tries +
                  kidx];
              }

              i = 7 * Tries + kidx_0;
              rtb_Add_k[i] = rtb_Product2_bg[i] + s;
            }
          }

          mrdiv_c(rtb_Add_k, rtb_y, rtb_N_f);
          kidx_0 = 0;
          for (Tries = 0; Tries < 3; Tries++) {
            for (kidx = 0; kidx < 7; kidx++) {
              i = kidx + kidx_0;
              rtb_Add_k[i] = 0.0;
              a_tmp = 0;
              for (I2Jm_tmp = 0; I2Jm_tmp < 7; I2Jm_tmp++) {
                rtb_Add_k[i] += rtDW.MemoryP_DSTATE_e[a_tmp + kidx] *
                  rtb_Transpose2[I2Jm_tmp + kidx_0];
                a_tmp += 7;
              }
            }

            kidx_0 += 7;
          }

          mrdiv_c(rtb_Add_k, rtb_y, rtb_Transpose2);
          for (kidx_0 = 0; kidx_0 < 49; kidx_0++) {
            b_I[kidx_0] = 0;
          }

          kidx_0 = 0;
          for (kidx = 0; kidx < 7; kidx++) {
            b_I[kidx_0] = 1;
            kidx_0 += 8;
          }

          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            for (Tries = 0; Tries < 7; Tries++) {
              I2Jm_tmp = 7 * Tries + kidx_0;
              rtb_y_g[I2Jm_tmp] = static_cast<real_T>(b_I[I2Jm_tmp]) -
                ((rtb_C_c[3 * Tries + 1] * rtb_Transpose2[kidx_0 + 7] + rtb_C_c
                  [3 * Tries] * rtb_Transpose2[kidx_0]) + rtb_C_c[3 * Tries + 2]
                 * rtb_Transpose2[kidx_0 + 14]);
            }

            for (Tries = 0; Tries < 7; Tries++) {
              I2Jm_tmp = 7 * Tries + kidx_0;
              Abar[I2Jm_tmp] = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                Abar[I2Jm_tmp] += rtb_y_g[7 * kidx + kidx_0] *
                  rtDW.MemoryP_DSTATE_e[7 * Tries + kidx];
              }
            }

            for (Tries = 0; Tries < 3; Tries++) {
              kidx = 7 * Tries + kidx_0;
              rtb_Add_k[kidx] = 0.0;
              rtb_Add_k[kidx] += rtb_R[3 * Tries] * rtb_Transpose2[kidx_0];
              rtb_Add_k[kidx] += rtb_R[3 * Tries + 1] * rtb_Transpose2[kidx_0 +
                7];
              rtb_Add_k[kidx] += rtb_R[3 * Tries + 2] * rtb_Transpose2[kidx_0 +
                14];
            }
          }

          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            Tries = 0;
            for (kidx = 0; kidx < 7; kidx++) {
              I2Jm_tmp = Tries + kidx_0;
              rtb_y_m[I2Jm_tmp] = 0.0;
              i = 0;
              for (a_tmp = 0; a_tmp < 7; a_tmp++) {
                rtb_y_m[I2Jm_tmp] += Abar[i + kidx_0] * rtb_y_g[i + kidx];
                i += 7;
              }

              rtb_Transpose2_0[I2Jm_tmp] = 0.0;
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Add_k[kidx_0] *
                rtb_Transpose2[kidx];
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Add_k[kidx_0 + 7] *
                rtb_Transpose2[kidx + 7];
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Add_k[kidx_0 + 14] *
                rtb_Transpose2[kidx + 14];
              Tries += 7;
            }
          }

          for (kidx_0 = 0; kidx_0 <= 46; kidx_0 += 2) {
            tmp_2 = _mm_loadu_pd(&rtb_y_m[kidx_0]);
            tmp_0 = _mm_loadu_pd(&rtb_Transpose2_0[kidx_0]);
            (void)_mm_storeu_pd(&rtb_Z[kidx_0], _mm_add_pd(tmp_2, tmp_0));
          }

          for (kidx_0 = 48; kidx_0 < 49; kidx_0++) {
            rtb_Z[kidx_0] = rtb_y_m[kidx_0] + rtb_Transpose2_0[kidx_0];
          }

          mrdiv_c(rtb_Product2_bg, rtb_R, rtb_Transpose2);
          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            for (Tries = 0; Tries < 7; Tries++) {
              I2Jm_tmp = 7 * Tries + kidx_0;
              rtb_y_g[I2Jm_tmp] = rtb_A_e[I2Jm_tmp] - ((rtb_C_c[3 * Tries + 1] *
                rtb_Transpose2[kidx_0 + 7] + rtb_C_c[3 * Tries] *
                rtb_Transpose2[kidx_0]) + rtb_C_c[3 * Tries + 2] *
                rtb_Transpose2[kidx_0 + 14]);
            }

            for (Tries = 0; Tries < 7; Tries++) {
              I2Jm_tmp = 7 * Tries + kidx_0;
              Abar[I2Jm_tmp] = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                Abar[I2Jm_tmp] += rtb_y_g[7 * kidx + kidx_0] * rtb_Z[7 * Tries +
                  kidx];
              }
            }
          }

          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            Tries = 0;
            for (kidx = 0; kidx < 7; kidx++) {
              s = 0.0;
              i = 0;
              for (a_tmp = 0; a_tmp < 7; a_tmp++) {
                s += Abar[i + kidx_0] * rtb_y_g[i + kidx];
                i += 7;
              }

              I2Jm_tmp = Tries + kidx_0;
              rtb_y_m[I2Jm_tmp] = rtb_Q_j[I2Jm_tmp] + s;
              rtb_Transpose2_0[I2Jm_tmp] = 0.0;
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Transpose2[kidx_0] *
                rtb_Product2_bg[kidx];
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Transpose2[kidx_0 + 7] *
                rtb_Product2_bg[kidx + 7];
              rtb_Transpose2_0[I2Jm_tmp] += rtb_Transpose2[kidx_0 + 14] *
                rtb_Product2_bg[kidx + 14];
              Tries += 7;
            }
          }

          for (kidx_0 = 0; kidx_0 <= 46; kidx_0 += 2) {
            tmp_2 = _mm_loadu_pd(&rtb_y_m[kidx_0]);
            tmp_0 = _mm_loadu_pd(&rtb_Transpose2_0[kidx_0]);
            (void)_mm_storeu_pd(&rtb_y_g[kidx_0], _mm_sub_pd(tmp_2, tmp_0));
          }

          for (kidx_0 = 48; kidx_0 < 49; kidx_0++) {
            rtb_y_g[kidx_0] = rtb_y_m[kidx_0] - rtb_Transpose2_0[kidx_0];
          }
        } else {
          (void)std::memset(&rtb_N_f[0], 0, 21U * sizeof(real_T));
          for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
            for (Tries = 0; Tries < 7; Tries++) {
              i = 7 * Tries + kidx_0;
              rtb_y_m[i] = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                rtb_y_m[i] += rtb_A_e[7 * kidx + kidx_0] *
                  rtDW.MemoryP_DSTATE_e[7 * Tries + kidx];
              }
            }

            for (Tries = 0; Tries < 7; Tries++) {
              s = 0.0;
              for (kidx = 0; kidx < 7; kidx++) {
                s += rtb_y_m[7 * kidx + kidx_0] * rtb_A_e[7 * kidx + Tries];
              }

              I2Jm_tmp = 7 * Tries + kidx_0;
              rtb_y_g[I2Jm_tmp] = rtb_Q_j[I2Jm_tmp] + s;
            }
          }
        }

        // End of MATLAB Function: '<S115>/Discrete-Time KF - Calculate PLMZ'
        // End of Outputs for SubSystem: '<S113>/CalculatePL'

        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S91>/umin_scale1'

        //  Determine if the Square-Root algorithm was used
        // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S155>:1' 
        // '<S155>:1:4' if isSqrtUsed
        Saturation_l_idx_0 = rtP.umin_scale1_Gain[0] * b_Wu[0];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S90>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[0] = Saturation_l_idx_0 - rtU.u0[0];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S91>/umin_scale1'

        Saturation_l_idx_0 = rtP.umin_scale1_Gain[1] * b_Wu[1];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S90>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[1] = Saturation_l_idx_0 - rtU.u0[1];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        // Saturate: '<S3>/Saturation' incorporates:
        //   Gain: '<S91>/umin_scale1'

        Saturation_l_idx_0 = rtP.umin_scale1_Gain[2] * b_Wu[2];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_k) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_k;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_h) {
          // Saturate: '<S3>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_h;
        } else {
          // no actions
        }

        // Sum: '<S90>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[2] = Saturation_l_idx_0 - rtU.u0[2];

        // Outputs for Enabled SubSystem: '<S132>/MeasurementUpdate' incorporates:
        //   EnablePort: '<S156>/Enable'

        if (rtP.Constant1_Value_e != 0.0) {
          rtDW.MeasurementUpdate_MODE = true;
          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            // Product: '<S156>/C[k]*xhat[k|k-1]' incorporates:
            //   Delay: '<S113>/MemoryX'

            rtb_C_1[kidx_0] = 0.0;
            Tries = 0;
            for (kidx = 0; kidx < 7; kidx++) {
              rtb_C_1[kidx_0] += rtb_C_c[Tries + kidx_0] *
                rtDW.MemoryX_DSTATE_l[kidx];
              Tries += 3;
            }

            // Product: '<S156>/D[k]*u[k]' incorporates:
            //   Product: '<S156>/C[k]*xhat[k|k-1]'

            tmp[kidx_0] = 0.0;
            tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0] * rtb_Product1_nb[0];
            tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0 + 3] * rtb_Product1_nb
              [1];
            tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0 + 6] * rtb_Product1_nb
              [2];

            // Sum: '<S156>/Sum' incorporates:
            //   Constant: '<S3>/Constant13'
            //   Product: '<S156>/C[k]*xhat[k|k-1]'
            //   Product: '<S156>/D[k]*u[k]'
            //   Sum: '<S156>/Add1'

            rtb_Sum6_c[kidx_0] = Sum2_c[kidx_0] - (rtb_C_1[kidx_0] + tmp[kidx_0]);
          }

          for (kidx_0 = 0; kidx_0 <= 4; kidx_0 += 2) {
            (void)_mm_storeu_pd(&rtDW.Product3_c[kidx_0], _mm_set1_pd(0.0));
            tmp_2 = _mm_loadu_pd(&rtb_N_f[kidx_0]);
            tmp_0 = _mm_loadu_pd(&rtDW.Product3_c[kidx_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[kidx_0], _mm_add_pd(_mm_mul_pd
              (tmp_2, _mm_set1_pd(rtb_Sum6_c[0])), tmp_0));
            tmp_2 = _mm_loadu_pd(&rtb_N_f[kidx_0 + 7]);
            tmp_0 = _mm_loadu_pd(&rtDW.Product3_c[kidx_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[kidx_0], _mm_add_pd(_mm_mul_pd
              (tmp_2, _mm_set1_pd(rtb_Sum6_c[1])), tmp_0));
            tmp_2 = _mm_loadu_pd(&rtb_N_f[kidx_0 + 14]);
            tmp_0 = _mm_loadu_pd(&rtDW.Product3_c[kidx_0]);
            (void)_mm_storeu_pd(&rtDW.Product3_c[kidx_0], _mm_add_pd(_mm_mul_pd
              (tmp_2, _mm_set1_pd(rtb_Sum6_c[2])), tmp_0));
          }

          for (kidx_0 = 6; kidx_0 < 7; kidx_0++) {
            // Product: '<S156>/Product3'
            rtDW.Product3_c[kidx_0] = 0.0;
            rtDW.Product3_c[kidx_0] += rtb_N_f[kidx_0] * rtb_Sum6_c[0];
            rtDW.Product3_c[kidx_0] += rtb_N_f[kidx_0 + 7] * rtb_Sum6_c[1];
            rtDW.Product3_c[kidx_0] += rtb_N_f[kidx_0 + 14] * rtb_Sum6_c[2];
          }
        } else if (rtDW.MeasurementUpdate_MODE) {
          for (Tries = 0; Tries < 7; Tries++) {
            // Disable for Product: '<S156>/Product3' incorporates:
            //   Outport: '<S156>/L*(y[k]-yhat[k|k-1])'
            //
            rtDW.Product3_c[Tries] = rtP.Lykyhatkk1_Y0_c;
          }

          rtDW.MeasurementUpdate_MODE = false;
        } else {
          // no actions
        }

        // End of Outputs for SubSystem: '<S132>/MeasurementUpdate'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          // Product: '<S116>/Product' incorporates:
          //   Delay: '<S113>/MemoryX'

          rtb_C_1[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 7; kidx++) {
            rtb_C_1[kidx_0] += rtb_C_c[Tries + kidx_0] *
              rtDW.MemoryX_DSTATE_l[kidx];
            Tries += 3;
          }

          // Product: '<S116>/Product1' incorporates:
          //   Product: '<S116>/Product'

          tmp[kidx_0] = 0.0;
          tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0] * rtb_Product1_nb[0];
          tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0 + 3] * rtb_Product1_nb[1];
          tmp[kidx_0] += rtP.Constant13_Value_c[kidx_0 + 6] * rtb_Product1_nb[2];

          // Sum: '<S116>/Add1' incorporates:
          //   Product: '<S116>/Product'
          //   Product: '<S116>/Product1'

          rtb_Sum6_c[kidx_0] = rtb_C_1[kidx_0] + tmp[kidx_0];
        }

        // Update for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' incorporates:
        //   Constant: '<S3>/Constant1'
        //   Constant: '<S3>/Constant13'
        //   DataTypeConversion: '<S113>/DataTypeConversionEnable'
        //   Inport: '<Root>/iRST'
        //   Product: '<S116>/Product'
        //   Product: '<S116>/Product1'
        //   Product: '<S156>/C[k]*xhat[k|k-1]'
        //   Product: '<S156>/D[k]*u[k]'
        //   Product: '<S156>/Product3'
        //   Sum: '<S3>/Sum'

        rtDW.DiscreteTimeIntegrator_DSTATE_j += (Sum2_c[0] - r__h_idx_0) *
          rtP.DiscreteTimeIntegrator_gainval;
        rtDW.DiscreteTimeIntegrator_PrevRe_b = static_cast<int8_T>(rtU.iRST ? 1 :
          0);

        // Update for Memory: '<S91>/Memory'
        (void)std::memcpy(&rtDW.Memory_PreviousInput_d[0], &rtb_iAout_l[0], 166U
                          * sizeof(boolean_T));

        // Update for UnitDelay: '<S91>/last_mv'
        rtDW.last_mv_DSTATE_n[0] = b_Wu[0];
        rtDW.last_mv_DSTATE_n[1] = b_Wu[1];
        rtDW.last_mv_DSTATE_n[2] = b_Wu[2];

        // Update for Delay: '<S113>/MemoryX'
        rtDW.icLoad_n = false;
        for (kidx_0 = 0; kidx_0 < 7; kidx_0++) {
          // Product: '<S132>/B[k]*u[k]'
          rtb_Sum2_f[kidx_0] = 0.0;
          rtb_Sum2_f[kidx_0] += rtb_B_o[kidx_0] * rtb_Product1_nb[0];
          rtb_Sum2_f[kidx_0] += rtb_B_o[kidx_0 + 7] * rtb_Product1_nb[1];
          rtb_Sum2_f[kidx_0] += rtb_B_o[kidx_0 + 14] * rtb_Product1_nb[2];

          // Product: '<S132>/A[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S113>/MemoryX'
          //   Product: '<S132>/B[k]*u[k]'

          rtb_A_j[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 7; kidx++) {
            rtb_A_j[kidx_0] += rtb_A_e[Tries + kidx_0] *
              rtDW.MemoryX_DSTATE_l[kidx];
            Tries += 7;
          }

          // End of Product: '<S132>/A[k]*xhat[k|k-1]'
        }

        // End of Outputs for SubSystem: '<S1>/mpc1'
        for (kidx_0 = 0; kidx_0 <= 4; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc1'
          tmp_2 = _mm_loadu_pd(&rtb_Sum2_f[kidx_0]);
          tmp_0 = _mm_loadu_pd(&rtb_A_j[kidx_0]);
          tmp_1 = _mm_loadu_pd(&rtDW.Product3_c[kidx_0]);
          (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE_l[kidx_0], _mm_add_pd
                              (_mm_add_pd(tmp_2, tmp_0), tmp_1));

          // End of Outputs for SubSystem: '<S1>/mpc1'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        for (kidx_0 = 6; kidx_0 < 7; kidx_0++) {
          // Update for Delay: '<S113>/MemoryX' incorporates:
          //   Sum: '<S132>/Add'

          rtDW.MemoryX_DSTATE_l[kidx_0] = (rtb_Sum2_f[kidx_0] + rtb_A_j[kidx_0])
            + rtDW.Product3_c[kidx_0];
        }

        // Update for Delay: '<S113>/MemoryP' incorporates:
        //   Delay: '<S113>/MemoryX'
        //   Product: '<S132>/B[k]*u[k]'
        //   Sum: '<S132>/Add'

        rtDW.icLoad_h = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE_e[0], &rtb_y_g[0], 49U * sizeof
                          (real_T));
        rtY.yhat[0] = rtb_Sum6_c[0];
        rtY.yhat[1] = rtb_Sum6_c[1];

        // End of Outputs for SubSystem: '<S1>/mpc1'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S90>/Sum3'

        rtY.u[2] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc1'
        rtY.yhat[2] = rtb_Sum6_c[2];

        // End of Outputs for SubSystem: '<S1>/mpc1'
      } else if (rtY.sig == 2.0) {
        real_T Saturation_l_idx_0;
        real_T Saturation_l_idx_1;
        real_T Saturation_l_idx_2;
        int32_T I2Jm_tmp;
        int32_T a_tmp;
        int32_T b_Linv_tmp;
        int32_T rtb_R_tmp_0;
        int16_T ixw;
        boolean_T guard11;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // DiscreteIntegrator: '<S4>/Discrete-Time Integrator' incorporates:
        //   Inport: '<Root>/iRST'

        // '<S1>:59:35' elseif sig == 2
        // '<S1>:59:36' [u, yhat(1:no)] = mpc2(r_, y__, [0;0;0], [0;0], u0, umax, uwt, iRST); 
        // Simulink Function 'mpc2': '<S1>:902'
        if (rtU.iRST && (rtDW.DiscreteTimeIntegrator_PrevRe_f <= 0)) {
          rtDW.DiscreteTimeIntegrator_DSTATE_m[0] =
            rtP.DiscreteTimeIntegrator_IC_n[0];
          rtDW.DiscreteTimeIntegrator_DSTATE_m[1] =
            rtP.DiscreteTimeIntegrator_IC_n[1];
        }

        // Gain: '<S4>/Gain1' incorporates:
        //   Inport: '<Root>/uwt'
        //   Product: '<S186>/Product1'

        rtb_Product1_nb[0] = rtP.beta * rtU.uwt[0];
        rtb_Product1_nb[1] = rtP.beta * rtU.uwt[1];
        rtb_Product1_nb[2] = rtP.beta * rtU.uwt[2];

        // Delay: '<S183>/MemoryX' incorporates:
        //   Constant: '<S183>/X0'
        //   DataTypeConversion: '<S183>/DataTypeConversionReset'

        rtDW.icLoad_a = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE_g) ==
                          POS_ZCSIG) || rtDW.icLoad_a);
        rtPrevZCX.MemoryX_Reset_ZCE_g = 0U;
        if (rtDW.icLoad_a) {
          (void)std::memcpy(&rtDW.MemoryX_DSTATE_c[0], &rtP.X0_Value_k[0],
                            sizeof(real_T) << 3UL);
        }

        // Sum: '<S160>/Sum2' incorporates:
        //   Delay: '<S183>/MemoryX'

        rtb_xest[0] = rtDW.MemoryX_DSTATE_c[0];
        rtb_xest[1] = rtDW.MemoryX_DSTATE_c[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'
        for (Tries = 0; Tries <= 4; Tries += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc2'
          tmp_2 = _mm_loadu_pd(&rtDW.MemoryX_DSTATE_c[Tries + 2]);
          (void)_mm_storeu_pd(&rtb_xest[Tries + 2], _mm_add_pd(tmp_2,
            _mm_loadu_pd(&rtP.Constant1_Value_p[Tries])));

          // End of Outputs for SubSystem: '<S1>/mpc2'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // SignalConversion generated from: '<S182>/ SFunction ' incorporates:
        //   Constant: '<S160>/Constant1'
        //   Constant: '<S4>/Constant'
        //   Delay: '<S183>/MemoryX'
        //   MATLAB Function: '<S181>/optimizer'
        //   Sum: '<S160>/Sum2'

        rtb_TmpSignalConversionAtSFun_h[0] = r__h_idx_0;
        rtb_TmpSignalConversionAtSFun_h[1] = r__h_idx_1;
        rtb_TmpSignalConversionAtSFun_h[2] = r__h_idx_2;
        rtb_TmpSignalConversionAtSFun_h[3] = rtP.Constant_Value_p[0];
        rtb_TmpSignalConversionAtSFun_h[4] = rtP.Constant_Value_p[1];

        // MATLAB Function: '<S181>/optimizer' incorporates:
        //   DiscreteIntegrator: '<S4>/Discrete-Time Integrator'
        //   Gain: '<S4>/Gain2'
        //   Math: '<S161>/Math Function1'
        //   Product: '<S186>/Product1'
        //   SignalConversion generated from: '<S182>/ SFunction '
        //   UnitDelay: '<S161>/last_mv'
        //
        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S182>:1'
        // '<S182>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S182>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S182>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S182>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S182>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S182>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S182>:1:252' if isa(u,'double')
        // '<S182>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S182>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S182>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S182>:1:252' if isa(u,'double')
        // '<S182>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S182>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S182>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S182>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S182>:1:252' if isa(u,'double')
        // '<S182>:1:253' y = u;
        //    md:         measured disturbance
        // '<S182>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S182>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S182>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S182>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S182>:1:252' if isa(u,'double')
        // '<S182>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S182>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S182>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S182>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S182>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S182>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S182>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S182>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S182>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S182>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S182>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S182>:1:58' uwt = convertDataType(uwt0,isDouble);
        // '<S182>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S182>:1:252' if isa(u,'double')
        // '<S182>:1:253' y = u;
        //    duwt:       run-time DMV weights
        // '<S182>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S182>:1:62' ewt = convertDataType(ewt0,isDouble);
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
        // '<S182>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S182>:1:96' isAdaptive = false;
        // '<S182>:1:97' isLTV = false;
        // '<S182>:1:98' ZERO = zeros('like',ref);
        // '<S182>:1:99' ONE = ones('like',ref);
        // '<S182>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S182>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S182>:1:119' if isSimulation
        // '<S182>:1:126' else
        //  When doing code generation, use M code directly
        // '<S182>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        kidx_0 = 0;
        for (kidx = 0; kidx < 20; kidx++) {
          for (Tries = 0; Tries < 5; Tries++) {
            rseq_0[Tries + kidx_0] = rtb_TmpSignalConversionAtSFun_h[Tries];
          }

          kidx_0 += 5;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S182>:1:133' old_u = old_u - uoff;
        Saturation_l_idx_0 = rtDW.last_mv_DSTATE_i[0];
        Saturation_l_idx_1 = rtDW.last_mv_DSTATE_i[1];
        Saturation_l_idx_2 = rtDW.last_mv_DSTATE_i[2];

        // '<S182>:1:134' if no_mv
        // '<S182>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S182>:1:143' xk = xk - xoff;
        rtb_xest_p[0] = rtb_xest[0];
        rtb_xest_p[2] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE_m[0];
        rtb_xest_p[1] = rtb_xest[1];
        rtb_xest_p[3] = rtP.dt * rtDW.DiscreteTimeIntegrator_DSTATE_m[1];
        for (Tries = 0; Tries < 6; Tries++) {
          rtb_xest_p[Tries + 4] = rtb_xest[Tries + 2];
        }

        //  Remove offset
        // '<S182>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S182>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S182>:1:162' if no_uref
        //  no external utarget
        // '<S182>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S182>:1:173' if ~no_cc
        // '<S182>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S182>:1:183' if isSimulation
        // '<S182>:1:214' else
        //  When doing code generation, use M code directly
        // '<S182>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S182>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S182>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S182>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S182>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S182>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S182>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S182>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S182>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S182>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S182>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S182>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S182>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        (void)std::memcpy(&b_Linv[0], &g_1[0], sizeof(real_T) << 4UL);
        (void)std::memset(&rtb_iAout_a[0], 0, 206U * sizeof(boolean_T));
        if (rtb_Product1_nb[0] < 0.0) {
          b_Wu[0] = 0.0;
        } else {
          b_Wu[0] = rtb_Product1_nb[0] * rtb_Product1_nb[0];
        }

        if (rtb_Product1_nb[1] < 0.0) {
          b_Wu[1] = 0.0;
        } else {
          b_Wu[1] = rtb_Product1_nb[1] * rtb_Product1_nb[1];
        }

        if (rtb_Product1_nb[2] < 0.0) {
          b_Wu[2] = 0.0;
        } else {
          b_Wu[2] = rtb_Product1_nb[2] * rtb_Product1_nb[2];
        }

        (void)std::memset(&rtb_y[0], 0, 9U * sizeof(real_T));
        rtb_y[0] = 1.0;
        rtb_y[4] = 1.0;
        rtb_y[8] = 1.0;
        kidx = -1;
        for (Tries = 0; Tries < 20; Tries++) {
          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            for (i = 0; i < 20; i++) {
              a_tmp = static_cast<int32_T>(b_A[20 * Tries + i]);
              a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0]) * a_tmp);
              a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 1]) * a_tmp);
              a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 2]) * a_tmp);
              kidx += 3;
            }
          }
        }

        kidx_0 = 0;
        for (Tries = 0; Tries < 3; Tries++) {
          for (kidx = 0; kidx < 60; kidx++) {
            I2Jm_tmp = kidx + kidx_0;
            I2Jm[I2Jm_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 60; a_tmp++) {
              I2Jm[I2Jm_tmp] += static_cast<real_T>(static_cast<int32_T>(
                static_cast<int32_T>(a[i + kidx]) * static_cast<int32_T>
                (b_Jm[a_tmp + kidx_0])));
              i += 60;
            }
          }

          kidx_0 += 60;
        }

        ixw = 1;
        for (kidx = 0; kidx < 100; kidx++) {
          dwt = W_1[ixw - 1];
          WySuJm_0[kidx] = dwt * b_SuJm_1[kidx];
          WySuJm_0[kidx + 100] = b_SuJm_1[kidx + 100] * dwt;
          WySuJm_0[kidx + 200] = b_SuJm_1[kidx + 200] * dwt;
          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 5) {
            ixw = 1;
          }
        }

        ixw = 1;
        for (kidx = 0; kidx < 60; kidx++) {
          dwt = b_Wu[ixw - 1];
          WuI2Jm[kidx] = dwt * I2Jm[kidx];
          WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * dwt;
          WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * dwt;
          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 3) {
            ixw = 1;
          }

          WduJm[kidx] = 0.034121465297356074 * static_cast<real_T>(b_Jm[kidx]);
          WduJm[kidx + 60] = static_cast<real_T>(b_Jm[kidx + 60]) *
            0.034121465297356074;
          WduJm[kidx + 120] = static_cast<real_T>(b_Jm[kidx + 120]) *
            0.034121465297356074;
        }

        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 3; Tries++) {
            rtb_R_tmp_0 = 3 * Tries + kidx_0;
            rtb_R[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 100; kidx++) {
              rtb_R[rtb_R_tmp_0] += e_1[3 * kidx + kidx_0] * WySuJm_0[100 *
                Tries + kidx];
            }

            s = 0.0;
            I2Jm_0[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 60; kidx++) {
              i = 60 * Tries + kidx;
              s += static_cast<real_T>(f[3 * kidx + kidx_0]) * WduJm[i];
              I2Jm_0[rtb_R_tmp_0] += I2Jm[60 * kidx_0 + kidx] * WuI2Jm[i];
            }

            rtb_R_0[rtb_R_tmp_0] = rtb_R[rtb_R_tmp_0] + s;
          }
        }

        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          i = 0;
          a_tmp = 0;
          for (I2Jm_tmp = 0; I2Jm_tmp < 3; I2Jm_tmp++) {
            b_Linv_tmp = I2Jm_tmp + Tries;
            b_Linv[I2Jm_tmp + kidx_0] = rtb_R_0[b_Linv_tmp] + I2Jm_0[b_Linv_tmp];
            s = 0.0;
            b_Linv_tmp = 0;
            for (rtb_R_tmp_0 = 0; rtb_R_tmp_0 < 60; rtb_R_tmp_0++) {
              s += static_cast<real_T>(d[b_Linv_tmp + kidx]) *
                WuI2Jm[rtb_R_tmp_0 + a_tmp];
              b_Linv_tmp += 3;
            }

            b_Linv_tmp = i + kidx;
            rtb_y[b_Linv_tmp] = rtb_R[b_Linv_tmp] + s;
            i += 3;
            a_tmp += 60;
          }

          kidx_0 += 4;
          Tries += 3;
        }

        // End of Outputs for SubSystem: '<S1>/mpc2'
        for (kidx_0 = 0; kidx_0 <= 178; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc2'
          tmp_2 = _mm_loadu_pd(&WuI2Jm[kidx_0]);
          (void)_mm_storeu_pd(&WuI2Jm[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc2'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 10; Tries++) {
            // MATLAB Function: '<S181>/optimizer'
            i = 10 * kidx_0 + Tries;
            b_Kx[i] = 0.0;
            for (kidx = 0; kidx < 100; kidx++) {
              b_Kx[i] += c_1[10 * kidx + Tries] * WySuJm_0[100 * kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 21; Tries++) {
            // MATLAB Function: '<S181>/optimizer'
            i = 21 * kidx_0 + Tries;
            B_est[i] = 0.0;
            for (kidx = 0; kidx < 100; kidx++) {
              B_est[i] += WySuJm_0[100 * kidx_0 + kidx] * 0.0;
            }
          }
        }

        // End of Outputs for SubSystem: '<S1>/mpc2'
        for (kidx_0 = 0; kidx_0 <= 298; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc2'
          tmp_2 = _mm_loadu_pd(&WySuJm_0[kidx_0]);
          (void)_mm_storeu_pd(&WySuJm_0[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc2'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // MATLAB Function: '<S181>/optimizer' incorporates:
        //   Inport: '<Root>/umax'
        //   Memory: '<S161>/Memory'
        //   UnitDelay: '<S161>/last_mv'

        kidx = 0;
        (void)std::memcpy(&b_L[0], &b_Linv[0], sizeof(real_T) << 4UL);
        Tries = xpotrf(b_L);
        guard11 = false;
        if (Tries == 0) {
          rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
          rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
          rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
          rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
          if (minimum(rtb_TmpSignalConversionAtSFun_p) > 1.4901161193847656E-7)
          {
          } else {
            guard11 = true;
          }
        } else {
          guard11 = true;
        }

        if (guard11) {
          boolean_T exitg2;
          dwt = 0.0;
          Tries = 0;
          exitg2 = false;
          while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (Tries < 4)) {
            s = ((std::abs(b_Linv[Tries + 4]) + std::abs(b_Linv[Tries])) + std::
                 abs(b_Linv[Tries + 8])) + std::abs(b_Linv[Tries + 12]);
            if (std::isnan(s)) {
              dwt = (rtNaN);
              exitg2 = true;
            } else {
              if (s > dwt) {
                dwt = s;
              }

              Tries++;
            }
          }

          if (dwt >= 1.0E+10) {
            kidx = 2;
          } else {
            Tries = 0;
            exitg1 = false;
            while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                     (0U)) == false) && (Tries <= 4)) {
              boolean_T guard2;
              dwt = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
                1.4901161193847656E-7;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b[kidx_0] = 0;
              }

              b[0] = 1;
              b[5] = 1;
              b[10] = 1;
              b[15] = 1;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b_Linv[kidx_0] += dwt * static_cast<real_T>(b[kidx_0]);
                b_L[kidx_0] = b_Linv[kidx_0];
              }

              kidx = xpotrf(b_L);
              guard2 = false;
              if (kidx == 0) {
                rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
                rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
                rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
                rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
                if (minimum(rtb_TmpSignalConversionAtSFun_p) >
                    1.4901161193847656E-7) {
                  kidx = 1;
                  exitg1 = true;
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                kidx = 3;
                Tries++;
              }
            }
          }
        }

        if (kidx > 1) {
          b_Wu[0] = rtDW.last_mv_DSTATE_i[0];
          b_Wu[1] = rtDW.last_mv_DSTATE_i[1];
          b_Wu[2] = rtDW.last_mv_DSTATE_i[2];
        } else {
          for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
            b[kidx_0] = 0;
          }

          b[0] = 1;
          b[5] = 1;
          b[10] = 1;
          b[15] = 1;
          kidx_0 = 0;
          for (kidx = 0; kidx < 4; kidx++) {
            b_Linv[kidx_0] = static_cast<real_T>(b[kidx_0]);
            b_Linv[kidx_0 + 1] = static_cast<real_T>(b[kidx_0 + 1]);
            b_Linv[kidx_0 + 2] = static_cast<real_T>(b[kidx_0 + 2]);
            b_Linv[kidx_0 + 3] = static_cast<real_T>(b[kidx_0 + 3]);
            kidx_0 += 4;
          }

          trisolve(b_L, b_Linv);
          umax_incr_flag[0] = false;
          b_Wu[0] = 0.0;
          umax_incr_flag[1] = false;
          b_Wu[1] = 0.0;
          umax_incr_flag[2] = false;
          b_Wu[2] = 0.0;
          for (kidx = 0; kidx < 206; kidx++) {
            dwt = b_Mlim_1[kidx];
            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 10; kidx_0++) {
              s += b_Mx_1[206 * kidx_0 + kidx] * rtb_xest_p[kidx_0];
            }

            s = -(((b_Mu1_1[kidx + 206] * Saturation_l_idx_1 + b_Mu1_1[kidx] *
                    Saturation_l_idx_0) + b_Mu1_1[kidx + 412] *
                   Saturation_l_idx_2) + (dwt + s));
            ixw = b_Mrows_2[kidx];
            if ((ixw > 100) && (ixw > 200) && (ixw <= 260)) {
              Tries = (static_cast<int32_T>(ixw) - div_nde_s32_floor(
                        static_cast<int32_T>(ixw) - 201, static_cast<int32_T>(nu))
                       * static_cast<int32_T>(nu)) - 201;
              rstP2 = umax_incr_flag[Tries];
              if (!umax_incr_flag[Tries]) {
                dwt = -rtU.umax[Tries] - (-dwt);
                rstP2 = true;
              } else {
                dwt = b_Wu[Tries];
              }

              b_Wu[Tries] = dwt;
              umax_incr_flag[Tries] = rstP2;
              s += dwt;
            }

            Bc_0[kidx] = s;
          }

          rtb_TmpSignalConversionAtSFun_p[0] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[1] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[2] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[3] = 0.0;
          for (kidx = 0; kidx < 3; kidx++) {
            real_T WuI2Jm_0;
            r__h_idx_0 = 0.0;
            for (kidx_0 = 0; kidx_0 < 10; kidx_0++) {
              r__h_idx_0 += b_Kx[10 * kidx + kidx_0] * rtb_xest_p[kidx_0];
            }

            dwt = 0.0;
            for (kidx_0 = 0; kidx_0 < 100; kidx_0++) {
              dwt += WySuJm_0[100 * kidx + kidx_0] * rseq_0[kidx_0];
            }

            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 21; kidx_0++) {
              s += B_est[21 * kidx + kidx_0];
            }

            WuI2Jm_0 = 0.0;
            for (kidx_0 = 0; kidx_0 < 60; kidx_0++) {
              WuI2Jm_0 += WuI2Jm[60 * kidx + kidx_0] * 0.0;
            }

            rtb_TmpSignalConversionAtSFun_p[kidx] = ((((rtb_y[3 * kidx + 1] *
              Saturation_l_idx_1 + rtb_y[3 * kidx] * Saturation_l_idx_0) +
              rtb_y[3 * kidx + 2] * Saturation_l_idx_2) + (r__h_idx_0 + dwt)) +
              s) + WuI2Jm_0;
          }

          (void)std::memcpy(&rtb_iAout_a[0], &rtDW.Memory_PreviousInput_c[0],
                            206U * sizeof(boolean_T));
          kidx_0 = 0;
          for (Tries = 0; Tries < 4; Tries++) {
            kidx = 0;
            for (i = 0; i < 4; i++) {
              b_Linv_tmp = kidx + Tries;
              b_L[b_Linv_tmp] = 0.0;
              b_L[b_Linv_tmp] += b_Linv[kidx_0] * b_Linv[kidx];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 1] * b_Linv[kidx + 1];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 2] * b_Linv[kidx + 2];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 3] * b_Linv[kidx + 3];
              kidx += 4;
            }

            kidx_0 += 4;
          }

          qpkwik_o(b_Linv, b_L, rtb_TmpSignalConversionAtSFun_p, b_Ac_1, Bc_0,
                   rtb_iAout_a, 840, 1.0E-6, zopt, a__1_0, &kidx);
          if ((kidx < 0) || (kidx == 0)) {
            zopt[0] = 0.0;
            zopt[1] = 0.0;
            zopt[2] = 0.0;
          }

          b_Wu[0] = rtDW.last_mv_DSTATE_i[0] + zopt[0];
          b_Wu[1] = rtDW.last_mv_DSTATE_i[1] + zopt[1];
          b_Wu[2] = rtDW.last_mv_DSTATE_i[2] + zopt[2];
        }

        // Delay: '<S183>/MemoryP' incorporates:
        //   Constant: '<S183>/P0'
        //   DataTypeConversion: '<S183>/DataTypeConversionReset'

        // '<S182>:1:231' if return_xseq || return_ovseq
        // '<S182>:1:233' else
        // '<S182>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S182>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S182>:1:238' if CustomEstimation
        // '<S182>:1:239' xk1 = xk;
        // '<S182>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S182>:1:247' xest = xest + xoff;
        rtDW.icLoad_p = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE_i) ==
                          POS_ZCSIG) || rtDW.icLoad_p);
        rtPrevZCX.MemoryP_Reset_ZCE_i = 0U;
        if (rtDW.icLoad_p) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE_h[0], &rtP.P0_Value_c[0],
                            sizeof(real_T) << 6UL);
        }

        // MATLAB Function: '<S160>/MATLAB Function' incorporates:
        //   Constant: '<S4>/Constant12'
        //   Constant: '<S4>/Constant13'
        //   Constant: '<S4>/Constant3'
        //   Constant: '<S4>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc2/State Estimator OD (KF)/MATLAB Function': '<S184>:1' 
        // '<S184>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod2, Bod2, Cod2(1:3,:), Dod2(1:3,:), Dmn1, 3, 3, 2); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&rtb_B[0], 0, 24U * sizeof(real_T));

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
        (void)std::memset(&rtb_A[0], 0, sizeof(real_T) << 6UL);
        rtb_A[0] = rtP.Constant3_Value_d[0];
        rtb_A[1] = rtP.Constant3_Value_d[1];
        rtb_A[8] = rtP.Constant3_Value_d[2];
        rtb_A[9] = rtP.Constant3_Value_d[3];
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 6; kidx++) {
          for (i = 0; i < 6; i++) {
            rtb_A[(i + kidx_0) + 18] = rtP.Aod2[i + Tries];
          }

          kidx_0 += 8;
          Tries += 6;
        }

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          rtb_B[kidx_0] = rtP.Constant4_Value_n[Tries];
          rtb_B[kidx_0 + 1] = rtP.Constant4_Value_n[Tries + 1];
          kidx_0 += 8;
          Tries += 2;
        }

        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        kidx_0 = 0;
        for (Tries = 0; Tries < 2; Tries++) {
          rtb_C[kidx_0] = rtP.Constant12_Value_i[kidx_0];
          rtb_C[kidx_0 + 1] = rtP.Constant12_Value_i[kidx_0 + 1];
          rtb_C[kidx_0 + 2] = rtP.Constant12_Value_i[kidx_0 + 2];
          kidx_0 += 3;
        }

        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 6; kidx++) {
          rtb_C[kidx_0 + 6] = rtP.Cod2[Tries];
          rtb_C[kidx_0 + 7] = rtP.Cod2[Tries + 1];
          rtb_C[kidx_0 + 8] = rtP.Cod2[Tries + 2];
          kidx_0 += 3;
          Tries += 5;
        }

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est_0[0], 0, 72U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y_0[0], 0, 48U * sizeof(real_T));
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          y_0[kidx_0] = rtP.Constant4_Value_n[Tries];
          y_0[kidx_0 + 1] = rtP.Constant4_Value_n[Tries + 1];
          kidx_0 += 8;
          Tries += 2;
        }

        for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
          y_0[kidx_0 + 26] = rtP.Bod2[kidx_0];
          y_0[kidx_0 + 34] = rtP.Bod2[kidx_0 + 6];
          y_0[kidx_0 + 42] = rtP.Bod2[kidx_0 + 12];
        }

        (void)std::memcpy(&B_est_0[0], &y_0[0], 48U * sizeof(real_T));

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          D_est[kidx_0] = rtP.Constant13_Value_g[kidx_0];
          D_est[kidx_0 + 9] = rtP.Dod2[Tries];
          D_est[kidx_0 + 18] = rtP.Dmn1[kidx_0];
          D_est[kidx_0 + 1] = rtP.Constant13_Value_g[kidx_0 + 1];
          D_est[kidx_0 + 10] = rtP.Dod2[Tries + 1];
          D_est[kidx_0 + 19] = rtP.Dmn1[kidx_0 + 1];
          D_est[kidx_0 + 2] = rtP.Constant13_Value_g[kidx_0 + 2];
          D_est[kidx_0 + 11] = rtP.Dod2[Tries + 2];
          D_est[kidx_0 + 20] = rtP.Dmn1[kidx_0 + 2];
          kidx_0 += 3;
          Tries += 5;
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            I2Jm_tmp = Tries + kidx_0;
            rtb_Q[I2Jm_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 9; a_tmp++) {
              rtb_Q[I2Jm_tmp] += B_est_0[i + kidx_0] * B_est_0[i + kidx];
              i += 8;
            }

            Tries += 8;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        kidx_0 = 0;
        for (Tries = 0; Tries < 9; Tries++) {
          rtb_R_tmp[Tries] = D_est[kidx_0];
          rtb_R_tmp[Tries + 9] = D_est[kidx_0 + 1];
          rtb_R_tmp[Tries + 18] = D_est[kidx_0 + 2];
          kidx_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 3; Tries++) {
            rtb_R_tmp_0 = 3 * kidx_0 + Tries;
            rtb_R[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_R[rtb_R_tmp_0] += D_est[3 * kidx + Tries] * rtb_R_tmp[9 *
                kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 8; Tries++) {
            i = (kidx_0 << 3UL) + Tries;
            rtb_N[i] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_N[i] += B_est_0[(kidx << 3UL) + Tries] * rtb_R_tmp[9 * kidx_0
                + kidx];
            }
          }
        }

        // End of MATLAB Function: '<S160>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S183>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S183>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S183>/ScalarExpansionQ'
        ScalarExpansionQ(rtb_Q, rtb_Z_e);

        // End of Outputs for SubSystem: '<S183>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S183>/ReducedQRN'
        ReducedQRN(rtP.G_Value_g, rtP.H_Value_k, rtb_Z_e, rtb_y, rtb_N,
                   rtb_Product, rtb_R, rtb_Product2);

        // End of Outputs for SubSystem: '<S183>/ReducedQRN'

        // Outputs for Atomic SubSystem: '<S183>/CalculatePL'
        CalculatePL(rtb_A, rtb_C, rtb_Product, rtb_R, rtb_Product2,
                    rtP.Constant1_Value_pe != 0.0, rtDW.MemoryP_DSTATE_h, rtb_N,
                    rtb_L, rtb_Z_e, rtb_Q);

        // End of Outputs for SubSystem: '<S183>/CalculatePL'

        // MATLAB Function: '<S224>/SqrtUsedFcn' incorporates:
        //   Constant: '<S183>/G'
        //   Constant: '<S183>/H'
        //   Constant: '<S224>/isSqrtUsed'
        //   Constant: '<S4>/Constant1'
        //   DataTypeConversion: '<S183>/DataTypeConversionEnable'
        //   Delay: '<S183>/MemoryP'

        SqrtUsedFcn(rtb_Z_e, rtP.isSqrtUsed_Value_a, rtb_Product);

        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S161>/umin_scale1'

        Saturation_l_idx_0 = rtP.umin_scale1_Gain_p[0] * b_Wu[0];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S160>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[0] = Saturation_l_idx_0 - rtU.u0[0];

        // Sum: '<S160>/Sum6'
        rtb_Sum6_c[0] = Sum2_c[0];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S161>/umin_scale1'

        Saturation_l_idx_0 = rtP.umin_scale1_Gain_p[1] * b_Wu[1];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S160>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[1] = Saturation_l_idx_0 - rtU.u0[1];

        // Sum: '<S160>/Sum6'
        rtb_Sum6_c[1] = Sum2_c[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Saturate: '<S4>/Saturation' incorporates:
        //   Gain: '<S161>/umin_scale1'

        Saturation_l_idx_0 = rtP.umin_scale1_Gain_p[2] * b_Wu[2];
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_h) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_h;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_o) {
          // Saturate: '<S4>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_o;
        } else {
          // no actions
        }

        // Sum: '<S160>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        rtb_Product1_nb[2] = Saturation_l_idx_0 - rtU.u0[2];

        // Sum: '<S160>/Sum6'
        rtb_Sum6_c[2] = Sum2_c[2];

        // Outputs for Enabled SubSystem: '<S202>/MeasurementUpdate'
        MeasurementUpdate(rtP.Constant1_Value_pe != 0.0, rtb_L, rtb_Sum6_c,
                          rtb_C, rtDW.MemoryX_DSTATE_c, rtP.Constant13_Value_g,
                          rtb_Product1_nb, rtDW.Product3_a,
                          &rtDW.MeasurementUpdate_j, &rtP.MeasurementUpdate_j);

        // End of Outputs for SubSystem: '<S202>/MeasurementUpdate'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          // Product: '<S186>/Product'
          rtb_C_1[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            rtb_C_1[kidx_0] += rtb_C[Tries + kidx_0] *
              rtDW.MemoryX_DSTATE_c[kidx];
            Tries += 3;
          }

          // Product: '<S186>/Product1' incorporates:
          //   Product: '<S186>/Product'

          tmp[kidx_0] = 0.0;
          tmp[kidx_0] += rtP.Constant13_Value_g[kidx_0] * rtb_Product1_nb[0];
          tmp[kidx_0] += rtP.Constant13_Value_g[kidx_0 + 3] * rtb_Product1_nb[1];
          tmp[kidx_0] += rtP.Constant13_Value_g[kidx_0 + 6] * rtb_Product1_nb[2];

          // Sum: '<S186>/Add1' incorporates:
          //   Product: '<S186>/Product'
          //   Product: '<S186>/Product1'

          rtb_Sum6_c[kidx_0] = rtb_C_1[kidx_0] + tmp[kidx_0];
        }

        // Update for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' incorporates:
        //   Constant: '<S4>/Constant1'
        //   Constant: '<S4>/Constant13'
        //   DataTypeConversion: '<S183>/DataTypeConversionEnable'
        //   Delay: '<S183>/MemoryX'
        //   Inport: '<Root>/iRST'
        //   Product: '<S186>/Product'
        //   Product: '<S186>/Product1'
        //   Sum: '<S4>/Sum'

        rtDW.DiscreteTimeIntegrator_DSTATE_m[0] += (Sum2_c[1] - r__h_idx_1) *
          rtP.DiscreteTimeIntegrator_gainva_b;
        rtDW.DiscreteTimeIntegrator_DSTATE_m[1] += (Sum2_c[2] - r__h_idx_2) *
          rtP.DiscreteTimeIntegrator_gainva_b;
        rtDW.DiscreteTimeIntegrator_PrevRe_f = static_cast<int8_T>(rtU.iRST ? 1 :
          0);

        // Update for Memory: '<S161>/Memory'
        (void)std::memcpy(&rtDW.Memory_PreviousInput_c[0], &rtb_iAout_a[0], 206U
                          * sizeof(boolean_T));

        // Update for UnitDelay: '<S161>/last_mv'
        rtDW.last_mv_DSTATE_i[0] = b_Wu[0];
        rtDW.last_mv_DSTATE_i[1] = b_Wu[1];
        rtDW.last_mv_DSTATE_i[2] = b_Wu[2];

        // Update for Delay: '<S183>/MemoryX'
        rtDW.icLoad_a = false;
        for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
          // Product: '<S202>/B[k]*u[k]'
          rtb_xest[kidx_0] = 0.0;
          rtb_xest[kidx_0] += rtb_B[kidx_0] * rtb_Product1_nb[0];
          rtb_xest[kidx_0] += rtb_B[kidx_0 + 8] * rtb_Product1_nb[1];
          rtb_xest[kidx_0] += rtb_B[kidx_0 + 16] * rtb_Product1_nb[2];

          // Product: '<S202>/A[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S183>/MemoryX'
          //   Product: '<S202>/B[k]*u[k]'

          rtb_A_1[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            rtb_A_1[kidx_0] += rtb_A[Tries + kidx_0] *
              rtDW.MemoryX_DSTATE_c[kidx];
            Tries += 8;
          }

          // End of Product: '<S202>/A[k]*xhat[k|k-1]'
        }

        // End of Outputs for SubSystem: '<S1>/mpc2'
        for (kidx_0 = 0; kidx_0 <= 6; kidx_0 += 2) {
          __m128d tmp_0;
          __m128d tmp_1;

          // Outputs for Function Call SubSystem: '<S1>/mpc2'
          tmp_2 = _mm_loadu_pd(&rtb_xest[kidx_0]);
          tmp_0 = _mm_loadu_pd(&rtb_A_1[kidx_0]);
          tmp_1 = _mm_loadu_pd(&rtDW.Product3_a[kidx_0]);
          (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE_c[kidx_0], _mm_add_pd
                              (_mm_add_pd(tmp_2, tmp_0), tmp_1));

          // End of Outputs for SubSystem: '<S1>/mpc2'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        // Update for Delay: '<S183>/MemoryP' incorporates:
        //   Delay: '<S183>/MemoryX'
        //   Product: '<S202>/B[k]*u[k]'
        //   Sum: '<S202>/Add'

        rtDW.icLoad_p = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE_h[0], &rtb_Q[0], sizeof(real_T) <<
                          6UL);
        rtY.yhat[0] = rtb_Sum6_c[0];
        rtY.yhat[1] = rtb_Sum6_c[1];

        // End of Outputs for SubSystem: '<S1>/mpc2'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S160>/Sum3'

        rtY.u[2] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc2'
        rtY.yhat[2] = rtb_Sum6_c[2];

        // End of Outputs for SubSystem: '<S1>/mpc2'
      } else if (rtY.sig == 3.0) {
        real_T Saturation_l_idx_0;
        real_T Saturation_l_idx_1;
        real_T Saturation_l_idx_2;
        int32_T I2Jm_tmp;
        int32_T a_tmp;
        int32_T b_Linv_tmp;
        int32_T rtb_R_tmp_0;
        int16_T ixw;
        boolean_T guard11{ false };

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Gain: '<S5>/Gain1' incorporates:
        //   Inport: '<Root>/uwt'
        //   Product: '<S256>/Product1'

        // '<S1>:59:37' elseif sig == 3
        // '<S1>:59:38' [u, yhat(1:no)] = mpc3(r_, y__, [0;0;0], [0;0], u0, umax, uwt, iRST); 
        // Simulink Function 'mpc3': '<S1>:936'
        rtb_Product1_nb[0] = rtP.beta * rtU.uwt[0];
        rtb_Product1_nb[1] = rtP.beta * rtU.uwt[1];
        rtb_Product1_nb[2] = rtP.beta * rtU.uwt[2];

        // Delay: '<S253>/MemoryX' incorporates:
        //   Constant: '<S253>/X0'
        //   DataTypeConversion: '<S253>/DataTypeConversionReset'

        rtDW.icLoad = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE) ==
                        POS_ZCSIG) || rtDW.icLoad);
        rtPrevZCX.MemoryX_Reset_ZCE = 0U;
        if (rtDW.icLoad) {
          (void)std::memcpy(&rtDW.MemoryX_DSTATE[0], &rtP.X0_Value_a[0], sizeof
                            (real_T) << 3UL);
        }

        // MATLAB Function: '<S251>/optimizer' incorporates:
        //   Delay: '<S253>/MemoryX'
        //   UnitDelay: '<S231>/last_mv'

        // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S252>:1'
        // '<S252>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
        // '<S252>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
        // '<S252>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
        // '<S252>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
        //  Inputs (in BlockDataType except iA)
        //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
        // '<S252>:1:24' xk = convertDataType(xk0,isDouble);
        // '<S252>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S252>:1:252' if isa(u,'double')
        // '<S252>:1:253' y = u;
        //    old_u:      last mv (calculated by MPC)
        // '<S252>:1:26' old_u = convertDataType(old_u0,isDouble);
        // '<S252>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S252>:1:252' if isa(u,'double')
        // '<S252>:1:253' y = u;
        //    ym:         current measured output (used only with built-in KF)
        // '<S252>:1:28' ym = convertDataType(ym0,isDouble);
        //    ref:        output reference
        // '<S252>:1:30' ref = convertDataType(ref0,isDouble);
        // '<S252>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S252>:1:252' if isa(u,'double')
        // '<S252>:1:253' y = u;
        //    md:         measured disturbance
        // '<S252>:1:32' md = convertDataType(md0,isDouble);
        //    umin:       run-time MV bound
        // '<S252>:1:34' umin = convertDataType(umin0,isDouble);
        //    umax:       run-time MV bound
        // '<S252>:1:36' umax = convertDataType(umax0,isDouble);
        // '<S252>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S252>:1:252' if isa(u,'double')
        // '<S252>:1:253' y = u;
        //    ymin:       run-time OV bound
        // '<S252>:1:38' ymin = convertDataType(ymin0,isDouble);
        //    ymax:       run-time OV bound
        // '<S252>:1:40' ymax = convertDataType(ymax0,isDouble);
        //    E:          run-time mixed constraints
        // '<S252>:1:42' E = convertDataType(E0,isDouble);
        //    F:          run-time mixed constraints
        // '<S252>:1:44' F = convertDataType(F0,isDouble);
        //    G:          run-time mixed constraints
        // '<S252>:1:46' G = convertDataType(G0,isDouble);
        //    S:          run-time mixed constraints
        // '<S252>:1:48' S = convertDataType(S0,isDouble);
        //    switch_in:  if it matches "enable_value", MPC is active in control 
        // '<S252>:1:50' switch_in = int32(switch_in0);
        //    ext_mv:     external last mv (actual)
        // '<S252>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
        //    MVtarget:   MV reference
        // '<S252>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
        //    ywt:        run-time OV weights
        // '<S252>:1:56' ywt = convertDataType(ywt0,isDouble);
        //    uwt:        run-time MV weights
        // '<S252>:1:58' uwt = convertDataType(uwt0,isDouble);
        // '<S252>:1:250' if isDouble
        //  convert an input signal to double precision when necessary
        // '<S252>:1:252' if isa(u,'double')
        // '<S252>:1:253' y = u;
        //    duwt:       run-time DMV weights
        // '<S252>:1:60' duwt = convertDataType(duwt0,isDouble);
        //    ewt:     run-time Slack weights
        // '<S252>:1:62' ewt = convertDataType(ewt0,isDouble);
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
        // '<S252>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
        // '<S252>:1:96' isAdaptive = false;
        // '<S252>:1:97' isLTV = false;
        // '<S252>:1:98' ZERO = zeros('like',ref);
        // '<S252>:1:99' ONE = ones('like',ref);
        // '<S252>:1:100' hasMD = nv>int32(1);
        //  Pre-allocate all the MEX block outputs for the simulation mode
        // '<S252>:1:105' if isSimulation
        //  Get reference and MD signals -- accounting for previewing
        // '<S252>:1:119' if isSimulation
        // '<S252>:1:126' else
        //  When doing code generation, use M code directly
        // '<S252>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
        kidx_0 = 0;
        for (kidx = 0; kidx < 20; kidx++) {
          rseq_1[kidx_0] = r__h_idx_0;
          rseq_1[kidx_0 + 1] = r__h_idx_1;
          rseq_1[kidx_0 + 2] = r__h_idx_2;
          kidx_0 += 3;
        }

        //  External MV override.
        //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
        // '<S252>:1:133' old_u = old_u - uoff;
        Saturation_l_idx_0 = rtDW.last_mv_DSTATE_k[0];
        Saturation_l_idx_1 = rtDW.last_mv_DSTATE_k[1];
        Saturation_l_idx_2 = rtDW.last_mv_DSTATE_k[2];

        // '<S252>:1:134' if no_mv
        // '<S252>:1:135' delmv = zeros(nu,1,'like',ref);
        //  Obtain x[k|k]
        // '<S252>:1:143' xk = xk - xoff;
        rtb_xest[0] = rtDW.MemoryX_DSTATE[0];
        rtb_xest[1] = rtDW.MemoryX_DSTATE[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'
        for (Tries = 0; Tries <= 4; Tries += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc3'
          tmp_2 = _mm_loadu_pd(&rtDW.MemoryX_DSTATE[Tries + 2]);
          (void)_mm_storeu_pd(&rtb_xest[Tries + 2], _mm_add_pd(tmp_2,
            _mm_loadu_pd(&rtP.Constant1_Value_h[Tries])));

          // End of Outputs for SubSystem: '<S1>/mpc3'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // MATLAB Function: '<S251>/optimizer' incorporates:
        //   Constant: '<S230>/Constant1'
        //   Delay: '<S253>/MemoryX'
        //   Math: '<S231>/Math Function1'
        //   Product: '<S256>/Product1'

        //  Remove offset
        // '<S252>:1:144' if CustomEstimation
        //  Input state is x(k|k)
        // '<S252>:1:146' xest = xk;
        //  Real-time MV target override
        //  Note: utargetValue is a vector length p*nu.
        // '<S252>:1:162' if no_uref
        //  no external utarget
        // '<S252>:1:164' utargetValue = utarget;
        //  Real-time custom constraint override (scaled E/F/S)
        // '<S252>:1:173' if ~no_cc
        // '<S252>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
        // '<S252>:1:183' if isSimulation
        // '<S252>:1:214' else
        //  When doing code generation, use M code directly
        // '<S252>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(... 
        // '<S252>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
        // '<S252>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
        // '<S252>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
        // '<S252>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
        // '<S252>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
        // '<S252>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
        // '<S252>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
        // '<S252>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
        // '<S252>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ... 
        // '<S252>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
        // '<S252>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
        // '<S252>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
        (void)std::memcpy(&b_Linv[0], &g[0], sizeof(real_T) << 4UL);
        (void)std::memset(&rtb_iAout_k[0], 0, 126U * sizeof(boolean_T));
        if (rtb_Product1_nb[0] < 0.0) {
          b_Wu[0] = 0.0;
        } else {
          b_Wu[0] = rtb_Product1_nb[0] * rtb_Product1_nb[0];
        }

        if (rtb_Product1_nb[1] < 0.0) {
          b_Wu[1] = 0.0;
        } else {
          b_Wu[1] = rtb_Product1_nb[1] * rtb_Product1_nb[1];
        }

        if (rtb_Product1_nb[2] < 0.0) {
          b_Wu[2] = 0.0;
        } else {
          b_Wu[2] = rtb_Product1_nb[2] * rtb_Product1_nb[2];
        }

        (void)std::memset(&rtb_y[0], 0, 9U * sizeof(real_T));
        rtb_y[0] = 1.0;
        rtb_y[4] = 1.0;
        rtb_y[8] = 1.0;
        kidx = -1;
        for (Tries = 0; Tries < 20; Tries++) {
          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            for (i = 0; i < 20; i++) {
              a_tmp = static_cast<int32_T>(b_A[20 * Tries + i]);
              a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0]) * a_tmp);
              a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 1]) * a_tmp);
              a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(rtb_y[3 *
                kidx_0 + 2]) * a_tmp);
              kidx += 3;
            }
          }
        }

        ixw = 1;
        for (kidx = 0; kidx < 60; kidx++) {
          dwt = W[ixw - 1];
          for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
            I2Jm_tmp = 60 * kidx_0 + kidx;
            I2Jm[I2Jm_tmp] = 0.0;
            for (Tries = 0; Tries < 60; Tries++) {
              I2Jm[I2Jm_tmp] += static_cast<real_T>(static_cast<int32_T>(
                static_cast<int32_T>(a[60 * Tries + kidx]) * static_cast<int32_T>
                (b_Jm[60 * kidx_0 + Tries])));
            }

            WySuJm_1[I2Jm_tmp] = b_SuJm_0[I2Jm_tmp] * dwt;
          }

          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 3) {
            ixw = 1;
          }
        }

        ixw = 1;
        for (kidx = 0; kidx < 60; kidx++) {
          dwt = b_Wu[ixw - 1];
          WuI2Jm[kidx] = dwt * I2Jm[kidx];
          WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * dwt;
          WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * dwt;
          ixw = static_cast<int16_T>(ixw + 1);
          if (ixw > 3) {
            ixw = 1;
          }

          WduJm[kidx] = 0.034121465297356074 * static_cast<real_T>(b_Jm[kidx]);
          WduJm[kidx + 60] = static_cast<real_T>(b_Jm[kidx + 60]) *
            0.034121465297356074;
          WduJm[kidx + 120] = static_cast<real_T>(b_Jm[kidx + 120]) *
            0.034121465297356074;
        }

        kidx_0 = 0;
        for (Tries = 0; Tries < 3; Tries++) {
          kidx = 0;
          i = 0;
          for (a_tmp = 0; a_tmp < 3; a_tmp++) {
            rtb_R_tmp_0 = kidx + Tries;
            rtb_R[rtb_R_tmp_0] = 0.0;
            s = 0.0;
            I2Jm_0[rtb_R_tmp_0] = 0.0;
            I2Jm_tmp = 0;
            for (b_Linv_tmp = 0; b_Linv_tmp < 60; b_Linv_tmp++) {
              int32_T rtb_R_tmp_1;
              int32_T rtb_R_tmp_2;
              rtb_R_tmp_1 = I2Jm_tmp + Tries;
              rtb_R_tmp_2 = b_Linv_tmp + i;
              rtb_R[rtb_R_tmp_0] += e_0[rtb_R_tmp_1] * WySuJm_1[rtb_R_tmp_2];
              s += static_cast<real_T>(f[rtb_R_tmp_1]) * WduJm[rtb_R_tmp_2];
              I2Jm_0[rtb_R_tmp_0] += I2Jm[b_Linv_tmp + kidx_0] *
                WuI2Jm[rtb_R_tmp_2];
              I2Jm_tmp += 3;
            }

            rtb_R_0[rtb_R_tmp_0] = rtb_R[rtb_R_tmp_0] + s;
            kidx += 3;
            i += 60;
          }

          kidx_0 += 60;
        }

        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          i = 0;
          a_tmp = 0;
          for (I2Jm_tmp = 0; I2Jm_tmp < 3; I2Jm_tmp++) {
            b_Linv_tmp = I2Jm_tmp + Tries;
            b_Linv[I2Jm_tmp + kidx_0] = rtb_R_0[b_Linv_tmp] + I2Jm_0[b_Linv_tmp];
            s = 0.0;
            b_Linv_tmp = 0;
            for (rtb_R_tmp_0 = 0; rtb_R_tmp_0 < 60; rtb_R_tmp_0++) {
              s += static_cast<real_T>(d[b_Linv_tmp + kidx]) *
                WuI2Jm[rtb_R_tmp_0 + a_tmp];
              b_Linv_tmp += 3;
            }

            b_Linv_tmp = i + kidx;
            rtb_y[b_Linv_tmp] = rtb_R[b_Linv_tmp] + s;
            i += 3;
            a_tmp += 60;
          }

          kidx_0 += 4;
          Tries += 3;
        }

        // End of Outputs for SubSystem: '<S1>/mpc3'
        for (kidx_0 = 0; kidx_0 <= 178; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc3'
          tmp_2 = _mm_loadu_pd(&WuI2Jm[kidx_0]);
          (void)_mm_storeu_pd(&WuI2Jm[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc3'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 8; Tries++) {
            // MATLAB Function: '<S251>/optimizer'
            i = (kidx_0 << 3UL) + Tries;
            rtb_B[i] = 0.0;
            for (kidx = 0; kidx < 60; kidx++) {
              rtb_B[i] += c[(kidx << 3UL) + Tries] * WySuJm_1[60 * kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 21; Tries++) {
            // MATLAB Function: '<S251>/optimizer'
            i = 21 * kidx_0 + Tries;
            B_est[i] = 0.0;
            for (kidx = 0; kidx < 60; kidx++) {
              B_est[i] += WySuJm_1[60 * kidx_0 + kidx] * 0.0;
            }
          }
        }

        // End of Outputs for SubSystem: '<S1>/mpc3'
        for (kidx_0 = 0; kidx_0 <= 178; kidx_0 += 2) {
          // Outputs for Function Call SubSystem: '<S1>/mpc3'
          tmp_2 = _mm_loadu_pd(&WySuJm_1[kidx_0]);
          (void)_mm_storeu_pd(&WySuJm_1[kidx_0], _mm_mul_pd(tmp_2, _mm_set1_pd
            (-1.0)));

          // End of Outputs for SubSystem: '<S1>/mpc3'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // MATLAB Function: '<S251>/optimizer' incorporates:
        //   Inport: '<Root>/umax'
        //   Memory: '<S231>/Memory'
        //   UnitDelay: '<S231>/last_mv'

        kidx = 0;
        (void)std::memcpy(&b_L[0], &b_Linv[0], sizeof(real_T) << 4UL);
        Tries = xpotrf(b_L);
        guard11 = false;
        if (Tries == 0) {
          rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
          rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
          rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
          rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
          if (minimum(rtb_TmpSignalConversionAtSFun_p) > 1.4901161193847656E-7)
          {
          } else {
            guard11 = true;
          }
        } else {
          guard11 = true;
        }

        if (guard11) {
          boolean_T exitg2;
          dwt = 0.0;
          Tries = 0;
          exitg2 = false;
          while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (Tries < 4)) {
            s = ((std::abs(b_Linv[Tries + 4]) + std::abs(b_Linv[Tries])) + std::
                 abs(b_Linv[Tries + 8])) + std::abs(b_Linv[Tries + 12]);
            if (std::isnan(s)) {
              dwt = (rtNaN);
              exitg2 = true;
            } else {
              if (s > dwt) {
                dwt = s;
              }

              Tries++;
            }
          }

          if (dwt >= 1.0E+10) {
            kidx = 2;
          } else {
            Tries = 0;
            exitg1 = false;
            while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                     (0U)) == false) && (Tries <= 4)) {
              boolean_T guard2{ false };

              dwt = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
                1.4901161193847656E-7;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b[kidx_0] = 0;
              }

              b[0] = 1;
              b[5] = 1;
              b[10] = 1;
              b[15] = 1;
              for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
                b_Linv[kidx_0] += dwt * static_cast<real_T>(b[kidx_0]);
                b_L[kidx_0] = b_Linv[kidx_0];
              }

              kidx = xpotrf(b_L);
              guard2 = false;
              if (kidx == 0) {
                rtb_TmpSignalConversionAtSFun_p[0] = b_L[0];
                rtb_TmpSignalConversionAtSFun_p[1] = b_L[5];
                rtb_TmpSignalConversionAtSFun_p[2] = b_L[10];
                rtb_TmpSignalConversionAtSFun_p[3] = b_L[15];
                if (minimum(rtb_TmpSignalConversionAtSFun_p) >
                    1.4901161193847656E-7) {
                  kidx = 1;
                  exitg1 = true;
                } else {
                  guard2 = true;
                }
              } else {
                guard2 = true;
              }

              if (guard2) {
                kidx = 3;
                Tries++;
              }
            }
          }
        }

        if (kidx > 1) {
          r__h_idx_0 = rtDW.last_mv_DSTATE_k[0];
          r__h_idx_1 = rtDW.last_mv_DSTATE_k[1];
          r__h_idx_2 = rtDW.last_mv_DSTATE_k[2];
        } else {
          for (kidx_0 = 0; kidx_0 < 16; kidx_0++) {
            b[kidx_0] = 0;
          }

          b[0] = 1;
          b[5] = 1;
          b[10] = 1;
          b[15] = 1;
          kidx_0 = 0;
          for (kidx = 0; kidx < 4; kidx++) {
            b_Linv[kidx_0] = static_cast<real_T>(b[kidx_0]);
            b_Linv[kidx_0 + 1] = static_cast<real_T>(b[kidx_0 + 1]);
            b_Linv[kidx_0 + 2] = static_cast<real_T>(b[kidx_0 + 2]);
            b_Linv[kidx_0 + 3] = static_cast<real_T>(b[kidx_0 + 3]);
            kidx_0 += 4;
          }

          trisolve(b_L, b_Linv);
          umax_incr_flag[0] = false;
          b_Wu[0] = 0.0;
          umax_incr_flag[1] = false;
          b_Wu[1] = 0.0;
          umax_incr_flag[2] = false;
          b_Wu[2] = 0.0;
          for (kidx = 0; kidx < 126; kidx++) {
            uint8_T b_Mrows;
            ixw = b_Mlim[kidx];
            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
              s += b_Mx[126 * kidx_0 + kidx] * rtb_xest[kidx_0];
            }

            s = -(((b_Mu1[kidx + 126] * Saturation_l_idx_1 + b_Mu1[kidx] *
                    Saturation_l_idx_0) + b_Mu1[kidx + 252] * Saturation_l_idx_2)
                  + (static_cast<real_T>(ixw) + s));
            b_Mrows = b_Mrows_0[kidx];
            if ((b_Mrows > 60UL) && (b_Mrows > 120UL) && (b_Mrows <= 180UL)) {
              Tries = (static_cast<int32_T>(b_Mrows) - div_nde_s32_floor(
                        static_cast<int32_T>(b_Mrows) - 121, static_cast<int32_T>
                        (nu)) * static_cast<int32_T>(nu)) - 121;
              rstP2 = umax_incr_flag[Tries];
              if (!umax_incr_flag[Tries]) {
                dwt = -rtU.umax[Tries] - (-static_cast<real_T>(ixw));
                rstP2 = true;
              } else {
                dwt = b_Wu[Tries];
              }

              b_Wu[Tries] = dwt;
              umax_incr_flag[Tries] = rstP2;
              s += dwt;
            }

            Bc_1[kidx] = s;
          }

          rtb_TmpSignalConversionAtSFun_p[0] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[1] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[2] = 0.0;
          rtb_TmpSignalConversionAtSFun_p[3] = 0.0;
          for (kidx = 0; kidx < 3; kidx++) {
            real_T WuI2Jm_0;
            r__h_idx_1 = 0.0;
            for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
              r__h_idx_1 += rtb_B[(kidx << 3UL) + kidx_0] * rtb_xest[kidx_0];
            }

            dwt = 0.0;
            for (kidx_0 = 0; kidx_0 < 60; kidx_0++) {
              dwt += WySuJm_1[60 * kidx + kidx_0] * rseq_1[kidx_0];
            }

            s = 0.0;
            for (kidx_0 = 0; kidx_0 < 21; kidx_0++) {
              s += B_est[21 * kidx + kidx_0];
            }

            WuI2Jm_0 = 0.0;
            for (kidx_0 = 0; kidx_0 < 60; kidx_0++) {
              WuI2Jm_0 += WuI2Jm[60 * kidx + kidx_0] * 0.0;
            }

            rtb_TmpSignalConversionAtSFun_p[kidx] = ((((rtb_y[3 * kidx + 1] *
              Saturation_l_idx_1 + rtb_y[3 * kidx] * Saturation_l_idx_0) +
              rtb_y[3 * kidx + 2] * Saturation_l_idx_2) + (r__h_idx_1 + dwt)) +
              s) + WuI2Jm_0;
          }

          (void)std::memcpy(&rtb_iAout_k[0], &rtDW.Memory_PreviousInput_n[0],
                            126U * sizeof(boolean_T));
          kidx_0 = 0;
          for (Tries = 0; Tries < 4; Tries++) {
            kidx = 0;
            for (i = 0; i < 4; i++) {
              b_Linv_tmp = kidx + Tries;
              b_L[b_Linv_tmp] = 0.0;
              b_L[b_Linv_tmp] += b_Linv[kidx_0] * b_Linv[kidx];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 1] * b_Linv[kidx + 1];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 2] * b_Linv[kidx + 2];
              b_L[b_Linv_tmp] += b_Linv[kidx_0 + 3] * b_Linv[kidx + 3];
              kidx += 4;
            }

            kidx_0 += 4;
          }

          qpkwik_f(b_Linv, b_L, rtb_TmpSignalConversionAtSFun_p, b_Ac, Bc_1,
                   rtb_iAout_k, 520, 1.0E-6, zopt, a__1_1, &kidx);
          if ((kidx < 0) || (kidx == 0)) {
            zopt[0] = 0.0;
            zopt[1] = 0.0;
            zopt[2] = 0.0;
          }

          r__h_idx_0 = rtDW.last_mv_DSTATE_k[0] + zopt[0];
          r__h_idx_1 = rtDW.last_mv_DSTATE_k[1] + zopt[1];
          r__h_idx_2 = rtDW.last_mv_DSTATE_k[2] + zopt[2];
        }

        // Delay: '<S253>/MemoryP' incorporates:
        //   Constant: '<S253>/P0'
        //   DataTypeConversion: '<S253>/DataTypeConversionReset'

        // '<S252>:1:231' if return_xseq || return_ovseq
        // '<S252>:1:233' else
        // '<S252>:1:234' yseq = zeros(p+1,ny,'like',rseq);
        // '<S252>:1:235' xseq = zeros(p+1,nxQP,'like',rseq);
        // '<S252>:1:238' if CustomEstimation
        // '<S252>:1:239' xk1 = xk;
        // '<S252>:1:244' xk1 = xk1 + xoff;
        //  Updated state must include offset
        //  return xest in original value
        // '<S252>:1:247' xest = xest + xoff;
        rtDW.icLoad_e = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE) ==
                          POS_ZCSIG) || rtDW.icLoad_e);
        rtPrevZCX.MemoryP_Reset_ZCE = 0U;
        if (rtDW.icLoad_e) {
          (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtP.P0_Value_m[0], sizeof
                            (real_T) << 6UL);
        }

        // MATLAB Function: '<S230>/MATLAB Function' incorporates:
        //   BusCreator: '<S5>/Bus Creator1'
        //   Constant: '<S5>/Constant12'
        //   Constant: '<S5>/Constant13'
        //   Constant: '<S5>/Constant3'
        //   Constant: '<S5>/Constant4'

        // MATLAB Function 'SupervisoryController/mpc3/State Estimator OD (KF)/MATLAB Function': '<S254>:1' 
        // '<S254>:1:2' [A, B, C, D, Q, R, N] = stateEst_(Ap, Bp, Cp, Dp, Aod3, Bod3, Cod3(1:3,:), Dod3(1:3,:), Dmn1, 3, 3, 2); 
        // 'stateEst_:3' nsp = ns_;
        //  n_plant_states
        // 'stateEst_:4' nsod = size(Aod,1);
        //  n_od_states
        // 'stateEst_:5' ns = nsp + nsod;
        //  n_states = n_plant_states + n_od_states
        // 'stateEst_:7' A = zeros(ns);
        //  n_states x n_states
        // 'stateEst_:8' B = zeros(ns,ni);
        (void)std::memset(&rtb_B[0], 0, 24U * sizeof(real_T));

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
        (void)std::memset(&rtb_A[0], 0, sizeof(real_T) << 6UL);
        rtb_A[0] = rtP.Constant3_Value_g[0];
        rtb_A[1] = rtP.Constant3_Value_g[1];
        rtb_A[8] = rtP.Constant3_Value_g[2];
        rtb_A[9] = rtP.Constant3_Value_g[3];
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 6; kidx++) {
          for (i = 0; i < 6; i++) {
            rtb_A[(i + kidx_0) + 18] = rtP.Aod3[i + Tries];
          }

          kidx_0 += 8;
          Tries += 6;
        }

        // 'stateEst_:20' B(1:nsp, 1:ni) = Bp;
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          rtb_B[kidx_0] = rtP.Constant4_Value_f[Tries];
          rtb_B[kidx_0 + 1] = rtP.Constant4_Value_f[Tries + 1];
          kidx_0 += 8;
          Tries += 2;
        }

        // 'stateEst_:21' C(1:no, 1:ns) = [Cp Cod];
        for (Tries = 0; Tries < 6; Tries++) {
          rtb_C[Tries] = rtP.Constant12_Value_f[Tries];
        }

        (void)std::memcpy(&rtb_C[6], &rtP.Cod3[0], 18U * sizeof(real_T));

        // 'stateEst_:22' D(1:no, 1:ni) = Dp;
        // 'stateEst_:24' B_est = zeros(ns, ni + no + no);
        (void)std::memset(&B_est_0[0], 0, 72U * sizeof(real_T));

        // 'stateEst_:25' B_est(1:ns, 1:ni+no) = blkdiag(Bp, Bod);
        (void)std::memset(&y_0[0], 0, 48U * sizeof(real_T));
        kidx_0 = 0;
        Tries = 0;
        for (kidx = 0; kidx < 3; kidx++) {
          y_0[kidx_0] = rtP.Constant4_Value_f[Tries];
          y_0[kidx_0 + 1] = rtP.Constant4_Value_f[Tries + 1];
          kidx_0 += 8;
          Tries += 2;
        }

        for (kidx_0 = 0; kidx_0 < 6; kidx_0++) {
          y_0[kidx_0 + 26] = rtP.Bod3[kidx_0];
          y_0[kidx_0 + 34] = rtP.Bod3[kidx_0 + 6];
          y_0[kidx_0 + 42] = rtP.Bod3[kidx_0 + 12];
        }

        (void)std::memcpy(&B_est_0[0], &y_0[0], 48U * sizeof(real_T));

        // 'stateEst_:26' D_est = [Dp Dod Dn];
        for (kidx_0 = 0; kidx_0 < 9; kidx_0++) {
          D_est[kidx_0] = rtP.Constant13_Value_a[kidx_0];
          D_est[kidx_0 + 9] = rtP.Dod3[kidx_0];
          D_est[kidx_0 + 18] = rtP.Dmn1[kidx_0];
        }

        // 'stateEst_:27' Q = B_est * B_est';
        for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            I2Jm_tmp = Tries + kidx_0;
            rtb_Q[I2Jm_tmp] = 0.0;
            i = 0;
            for (a_tmp = 0; a_tmp < 9; a_tmp++) {
              rtb_Q[I2Jm_tmp] += B_est_0[i + kidx_0] * B_est_0[i + kidx];
              i += 8;
            }

            Tries += 8;
          }
        }

        // 'stateEst_:28' R = D_est * D_est';
        kidx_0 = 0;
        for (Tries = 0; Tries < 9; Tries++) {
          rtb_R_tmp[Tries] = D_est[kidx_0];
          rtb_R_tmp[Tries + 9] = D_est[kidx_0 + 1];
          rtb_R_tmp[Tries + 18] = D_est[kidx_0 + 2];
          kidx_0 += 3;
        }

        // 'stateEst_:29' N = B_est * D_est';
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          for (Tries = 0; Tries < 3; Tries++) {
            rtb_R_tmp_0 = 3 * kidx_0 + Tries;
            rtb_R[rtb_R_tmp_0] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_R[rtb_R_tmp_0] += D_est[3 * kidx + Tries] * rtb_R_tmp[9 *
                kidx_0 + kidx];
            }
          }

          for (Tries = 0; Tries < 8; Tries++) {
            i = (kidx_0 << 3UL) + Tries;
            rtb_N[i] = 0.0;
            for (kidx = 0; kidx < 9; kidx++) {
              rtb_N[i] += B_est_0[(kidx << 3UL) + Tries] * rtb_R_tmp[9 * kidx_0
                + kidx];
            }
          }
        }

        // End of MATLAB Function: '<S230>/MATLAB Function'

        // Outputs for Atomic SubSystem: '<S253>/ScalarExpansionR'
        //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
        //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
        //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
        //  yhat = C*xhat + D*u;
        ScalarExpansionR(rtb_R, rtb_y);

        // End of Outputs for SubSystem: '<S253>/ScalarExpansionR'

        // Outputs for Atomic SubSystem: '<S253>/ScalarExpansionQ'
        ScalarExpansionQ(rtb_Q, rtb_Z_e);

        // End of Outputs for SubSystem: '<S253>/ScalarExpansionQ'

        // Outputs for Atomic SubSystem: '<S253>/ReducedQRN'
        ReducedQRN(rtP.G_Value_h, rtP.H_Value_oa, rtb_Z_e, rtb_y, rtb_N,
                   rtb_Product, rtb_R, rtb_Product2);

        // End of Outputs for SubSystem: '<S253>/ReducedQRN'

        // Outputs for Atomic SubSystem: '<S253>/CalculatePL'
        CalculatePL(rtb_A, rtb_C, rtb_Product, rtb_R, rtb_Product2,
                    rtP.Constant1_Value_n != 0.0, rtDW.MemoryP_DSTATE, rtb_N,
                    rtb_L, rtb_Z_e, rtb_Q);

        // End of Outputs for SubSystem: '<S253>/CalculatePL'

        // MATLAB Function: '<S294>/SqrtUsedFcn' incorporates:
        //   Constant: '<S253>/G'
        //   Constant: '<S253>/H'
        //   Constant: '<S294>/isSqrtUsed'
        //   Constant: '<S5>/Constant1'
        //   DataTypeConversion: '<S253>/DataTypeConversionEnable'
        //   Delay: '<S253>/MemoryP'

        SqrtUsedFcn(rtb_Z_e, rtP.isSqrtUsed_Value_p, rtb_Product);

        // Gain: '<S231>/umin_scale1'
        Saturation_l_idx_0 = rtP.umin_scale1_Gain_g[0] * r__h_idx_0;

        // Saturate: '<S5>/Saturation'
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S230>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        b_Wu[0] = Saturation_l_idx_0 - rtU.u0[0];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[0] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Gain: '<S231>/umin_scale1'
        Saturation_l_idx_0 = rtP.umin_scale1_Gain_g[1] * r__h_idx_1;

        // Saturate: '<S5>/Saturation'
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S230>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        b_Wu[1] = Saturation_l_idx_0 - rtU.u0[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u'
        rtY.u[1] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Gain: '<S231>/umin_scale1'
        Saturation_l_idx_0 = rtP.umin_scale1_Gain_g[2] * r__h_idx_2;

        // Saturate: '<S5>/Saturation'
        if (Saturation_l_idx_0 > rtP.Saturation_UpperSat_c) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_UpperSat_c;
        } else if (Saturation_l_idx_0 < rtP.Saturation_LowerSat_b) {
          // Saturate: '<S5>/Saturation'
          Saturation_l_idx_0 = rtP.Saturation_LowerSat_b;
        } else {
          // no actions
        }

        // Sum: '<S230>/Sum1' incorporates:
        //   Inport: '<Root>/u0'

        b_Wu[2] = Saturation_l_idx_0 - rtU.u0[2];

        // Outputs for Enabled SubSystem: '<S272>/MeasurementUpdate'
        MeasurementUpdate(rtP.Constant1_Value_n != 0.0, rtb_L, Sum2_c, rtb_C,
                          rtDW.MemoryX_DSTATE, rtP.Constant13_Value_a, b_Wu,
                          rtDW.Product3, &rtDW.MeasurementUpdate_c,
                          &rtP.MeasurementUpdate_c);

        // End of Outputs for SubSystem: '<S272>/MeasurementUpdate'
        for (kidx_0 = 0; kidx_0 < 3; kidx_0++) {
          // Product: '<S256>/Product'
          rtb_C_1[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            rtb_C_1[kidx_0] += rtb_C[Tries + kidx_0] * rtDW.MemoryX_DSTATE[kidx];
            Tries += 3;
          }

          // Product: '<S256>/Product1' incorporates:
          //   Product: '<S256>/Product'

          tmp[kidx_0] = 0.0;
          tmp[kidx_0] += rtP.Constant13_Value_a[kidx_0] * b_Wu[0];
          tmp[kidx_0] += rtP.Constant13_Value_a[kidx_0 + 3] * b_Wu[1];
          tmp[kidx_0] += rtP.Constant13_Value_a[kidx_0 + 6] * b_Wu[2];

          // Sum: '<S256>/Add1' incorporates:
          //   Product: '<S256>/Product'
          //   Product: '<S256>/Product1'

          Sum2_c[kidx_0] = rtb_C_1[kidx_0] + tmp[kidx_0];
        }

        // Update for Memory: '<S231>/Memory' incorporates:
        //   Constant: '<S5>/Constant1'
        //   Constant: '<S5>/Constant13'
        //   DataTypeConversion: '<S253>/DataTypeConversionEnable'
        //   Delay: '<S253>/MemoryX'
        //   Product: '<S256>/Product'
        //   Product: '<S256>/Product1'

        (void)std::memcpy(&rtDW.Memory_PreviousInput_n[0], &rtb_iAout_k[0], 126U
                          * sizeof(boolean_T));

        // Update for UnitDelay: '<S231>/last_mv'
        rtDW.last_mv_DSTATE_k[0] = r__h_idx_0;
        rtDW.last_mv_DSTATE_k[1] = r__h_idx_1;
        rtDW.last_mv_DSTATE_k[2] = r__h_idx_2;

        // Update for Delay: '<S253>/MemoryX'
        rtDW.icLoad = false;
        for (kidx_0 = 0; kidx_0 < 8; kidx_0++) {
          // Product: '<S272>/B[k]*u[k]'
          rtb_xest[kidx_0] = 0.0;
          rtb_xest[kidx_0] += rtb_B[kidx_0] * b_Wu[0];
          rtb_xest[kidx_0] += rtb_B[kidx_0 + 8] * b_Wu[1];
          rtb_xest[kidx_0] += rtb_B[kidx_0 + 16] * b_Wu[2];

          // Product: '<S272>/A[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S253>/MemoryX'
          //   Product: '<S272>/B[k]*u[k]'

          rtb_A_1[kidx_0] = 0.0;
          Tries = 0;
          for (kidx = 0; kidx < 8; kidx++) {
            rtb_A_1[kidx_0] += rtb_A[Tries + kidx_0] * rtDW.MemoryX_DSTATE[kidx];
            Tries += 8;
          }

          // End of Product: '<S272>/A[k]*xhat[k|k-1]'
        }

        // End of Outputs for SubSystem: '<S1>/mpc3'
        for (kidx_0 = 0; kidx_0 <= 6; kidx_0 += 2) {
          __m128d tmp_0;
          __m128d tmp_1;

          // Outputs for Function Call SubSystem: '<S1>/mpc3'
          tmp_2 = _mm_loadu_pd(&rtb_xest[kidx_0]);
          tmp_0 = _mm_loadu_pd(&rtb_A_1[kidx_0]);
          tmp_1 = _mm_loadu_pd(&rtDW.Product3[kidx_0]);
          (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE[kidx_0], _mm_add_pd
                              (_mm_add_pd(tmp_2, tmp_0), tmp_1));

          // End of Outputs for SubSystem: '<S1>/mpc3'
        }

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        // Update for Delay: '<S253>/MemoryP' incorporates:
        //   Delay: '<S253>/MemoryX'
        //   Product: '<S272>/B[k]*u[k]'
        //   Sum: '<S272>/Add'

        rtDW.icLoad_e = false;
        (void)std::memcpy(&rtDW.MemoryP_DSTATE[0], &rtb_Q[0], sizeof(real_T) <<
                          6UL);
        rtY.yhat[0] = Sum2_c[0];
        rtY.yhat[1] = Sum2_c[1];

        // End of Outputs for SubSystem: '<S1>/mpc3'

        // Outport: '<Root>/u' incorporates:
        //   Outport: '<Root>/yhat'
        //   Sum: '<S230>/Sum3'

        rtY.u[2] = Saturation_l_idx_0;

        // Outputs for Function Call SubSystem: '<S1>/mpc3'
        rtY.yhat[2] = Sum2_c[2];

        // End of Outputs for SubSystem: '<S1>/mpc3'
        // u = mpcg(r_, y__, umax, uwt, iRST);
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
      rtY.P_p[i] = static_cast<real_T>(tmp_3[i]);
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

    // InitializeConditions for DiscreteIntegrator: '<S3>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE_j = rtP.DiscreteTimeIntegrator_IC;
    rtDW.DiscreteTimeIntegrator_PrevRe_b = 2;

    // InitializeConditions for Memory: '<S91>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_d[0],
                      &rtP.Memory_InitialCondition_f[0], 166U * sizeof(boolean_T));

    // InitializeConditions for UnitDelay: '<S91>/last_mv'
    rtDW.last_mv_DSTATE_n[0] = rtP.last_mv_InitialCondition_f[0];
    rtDW.last_mv_DSTATE_n[1] = rtP.last_mv_InitialCondition_f[1];
    rtDW.last_mv_DSTATE_n[2] = rtP.last_mv_InitialCondition_f[2];

    // InitializeConditions for Delay: '<S113>/MemoryX'
    rtDW.icLoad_n = true;

    // InitializeConditions for Delay: '<S113>/MemoryP'
    rtDW.icLoad_h = true;

    // SystemInitialize for Enabled SubSystem: '<S132>/MeasurementUpdate'
    for (i = 0; i < 7; i++) {
      // SystemInitialize for Product: '<S156>/Product3' incorporates:
      //   Outport: '<S156>/L*(y[k]-yhat[k|k-1])'

      rtDW.Product3_c[i] = rtP.Lykyhatkk1_Y0_c;
    }

    // End of SystemInitialize for SubSystem: '<S132>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc2'

    // InitializeConditions for DiscreteIntegrator: '<S4>/Discrete-Time Integrator' 
    rtDW.DiscreteTimeIntegrator_DSTATE_m[0] = rtP.DiscreteTimeIntegrator_IC_n[0];
    rtDW.DiscreteTimeIntegrator_DSTATE_m[1] = rtP.DiscreteTimeIntegrator_IC_n[1];
    rtDW.DiscreteTimeIntegrator_PrevRe_f = 2;

    // InitializeConditions for Memory: '<S161>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_c[0],
                      &rtP.Memory_InitialCondition_j[0], 206U * sizeof(boolean_T));

    // InitializeConditions for UnitDelay: '<S161>/last_mv'
    rtDW.last_mv_DSTATE_i[0] = rtP.last_mv_InitialCondition_b[0];
    rtDW.last_mv_DSTATE_i[1] = rtP.last_mv_InitialCondition_b[1];
    rtDW.last_mv_DSTATE_i[2] = rtP.last_mv_InitialCondition_b[2];

    // InitializeConditions for Delay: '<S183>/MemoryX'
    rtDW.icLoad_a = true;

    // InitializeConditions for Delay: '<S183>/MemoryP'
    rtDW.icLoad_p = true;

    // SystemInitialize for Enabled SubSystem: '<S202>/MeasurementUpdate'
    MeasurementUpdate_Init(rtDW.Product3_a, &rtP.MeasurementUpdate_j);

    // End of SystemInitialize for SubSystem: '<S202>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc3'

    // InitializeConditions for Memory: '<S231>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_n[0],
                      &rtP.Memory_InitialCondition_b[0], 126U * sizeof(boolean_T));

    // InitializeConditions for UnitDelay: '<S231>/last_mv'
    rtDW.last_mv_DSTATE_k[0] = rtP.last_mv_InitialCondition_i[0];
    rtDW.last_mv_DSTATE_k[1] = rtP.last_mv_InitialCondition_i[1];
    rtDW.last_mv_DSTATE_k[2] = rtP.last_mv_InitialCondition_i[2];

    // InitializeConditions for Delay: '<S253>/MemoryX'
    rtDW.icLoad = true;

    // InitializeConditions for Delay: '<S253>/MemoryP'
    rtDW.icLoad_e = true;

    // SystemInitialize for Enabled SubSystem: '<S272>/MeasurementUpdate'
    MeasurementUpdate_Init(rtDW.Product3, &rtP.MeasurementUpdate_c);

    // End of SystemInitialize for SubSystem: '<S272>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/wtMod'

    for (i = 0; i < 6; i++) {
      // InitializeConditions for Delay: '<S9>/Delay'
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
