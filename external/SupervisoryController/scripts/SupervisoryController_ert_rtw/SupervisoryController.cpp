//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
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

const int32_T ny{ 6 };

const int32_T p{ 20 };

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

  // MATLAB Function: '<S158>/MATLAB Function1'
  in3_idx_0 = in3 - in2;
  stride_0_0 = (in7 - in6) + 1 != 1 ? static_cast<int32_T>(1) : static_cast<
    int32_T>(0);
  for (int32_T i{0}; i < in3_idx_0; i++) {
    in1[(in2 + i) + 12 * in4] = in5[i * stride_0_0 + in6] * in8[(in4 << 2UL) + i];
  }

  // End of MATLAB Function: '<S158>/MATLAB Function1'
}

//
// System initialize for function-call system:
//    '<S1>/paramEst1'
//    '<S1>/paramEst2'
//
void SupervisoryController::paramEst1_Init(real_T rty_theta[12], real_T rty_P
  [144], real_T rty_err[3], DW_paramEst1 *localDW, P_paramEst1 *localP)
{
  // InitializeConditions for Delay: '<S160>/Delay1'
  localDW->icLoad = true;

  // InitializeConditions for UnitDelay: '<S158>/Unit Delay3'
  localDW->UnitDelay3_DSTATE[0] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[1] = localP->UnitDelay3_InitialCondition;
  localDW->UnitDelay3_DSTATE[2] = localP->UnitDelay3_InitialCondition;

  // InitializeConditions for Delay: '<S160>/Delay'
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

  // Delay: '<S160>/Delay1'
  localDW->icLoad = ((rtu_rstP && (static_cast<uint32_T>
    (localZCE->Delay1_Reset_ZCE) != POS_ZCSIG)) || localDW->icLoad);
  localZCE->Delay1_Reset_ZCE = rtu_rstP ? static_cast<ZCSigState>(1) :
    static_cast<ZCSigState>(0);
  if (localDW->icLoad) {
    (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rtu_P0[0], 144U * sizeof
                      (real_T));
  }

  // Sum: '<S158>/Add1'
  rtb_Add1_c[0] = rtu_y[0] - rtu_y0[0];

  // Sum: '<S158>/Add3'
  rtb_Add3_idx_0 = rtu_u[0] - rtu_u0[0];

  // Sum: '<S158>/Add1'
  rtb_Add1_c[1] = rtu_y[1] - rtu_y0[1];

  // Sum: '<S158>/Add3'
  rtb_Add3_idx_1 = rtu_u[1] - rtu_u0[1];

  // Sum: '<S158>/Add1'
  rtb_Add1_c[2] = rtu_y[2] - rtu_y0[2];

  // Sum: '<S158>/Add3'
  rtb_Add3_idx_2 = rtu_u[2] - rtu_u0[2];

  // MATLAB Function: '<S158>/MATLAB Function1' incorporates:
  //   Sum: '<S158>/Add1'
  //   Sum: '<S158>/Add3'

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/MATLAB Function1': '<S159>:1' 
  // '<S159>:1:2' [z, phi] = getRegressors_(y, yPrev, u, sign_, no, ni, np, dt, mdlNum); 
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

  // Delay: '<S160>/Delay'
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
    // Math: '<S160>/Transpose' incorporates:
    //   MATLAB Function: '<S160>/MATLAB Function'

    tmp[3 * i] = rtb_phi[i];
    tmp[3 * i + 1] = rtb_phi[i + 12];
    tmp[3 * i + 2] = rtb_phi[i + 24];
  }

  // MATLAB Function 'SupervisoryController/paramEst1/Param Estimator (RLS)/RLS/MATLAB Function': '<S161>:1' 
  // '<S161>:1:2' [dtheta, dP, L] = rls_(theta, phi, epsil, EN, p_, dPmod_, lambda, P, no, ni, np); 
  // 'rls_:3' dtheta = zeros(no*np, 1);
  // 'rls_:4' dP = zeros(no*np,no*np);
  // 'rls_:5' L = zeros(no*np, no);
  // 'rls_:7' L = P*phi*inv(lambda*eye(no) + phi'*P*phi);
  for (i = 0; i < 3; i++) {
    // Sum: '<S160>/Sum2' incorporates:
    //   Delay: '<S160>/Delay'
    //   MATLAB Function: '<S158>/MATLAB Function1'
    //   Math: '<S160>/Transpose'
    //   Sum: '<S158>/Add1'
    //   UnitDelay: '<S158>/Unit Delay3'

    rtb_Add3_idx_0 = 0.0;
    for (p2 = 0; p2 < 12; p2++) {
      // MATLAB Function: '<S160>/MATLAB Function'
      p1 = 3 * p2 + i;
      rtb_Add3_idx_0 += tmp[p1] * localDW->Delay_DSTATE[p2];

      // MATLAB Function: '<S160>/MATLAB Function' incorporates:
      //   Delay: '<S160>/Delay'
      //   Delay: '<S160>/Delay1'

      tmp_0[p1] = 0.0;
      for (ibtile = 0; ibtile < 12; ibtile++) {
        tmp_0[p1] += tmp[3 * ibtile + i] * localDW->Delay1_DSTATE[12 * p2 +
          ibtile];
      }
    }

    rty_err[i] = (rtb_Add1_c[i] - localDW->UnitDelay3_DSTATE[i]) -
      rtb_Add3_idx_0;

    // End of Sum: '<S160>/Sum2'

    // MATLAB Function: '<S160>/MATLAB Function'
    for (p2 = 0; p2 < 3; p2++) {
      rtb_Add3_idx_0 = 0.0;
      for (p1 = 0; p1 < 12; p1++) {
        rtb_Add3_idx_0 += tmp_0[3 * p1 + i] * rtb_phi[12 * p2 + p1];
      }

      p1 = 3 * p2 + i;
      b_b[p1] = static_cast<real_T>(b_b_0[p1]) * rtu_lambda + rtb_Add3_idx_0;
    }
  }

  // MATLAB Function: '<S160>/MATLAB Function' incorporates:
  //   Delay: '<S160>/Delay'
  //   Delay: '<S160>/Delay1'

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
    // Delay: '<S160>/Delay1'
    tmp_1 = _mm_loadu_pd(&localDW->Delay1_DSTATE[i]);

    // Product: '<S160>/Product1' incorporates:
    //   Delay: '<S160>/Delay1'

    tmp_2 = _mm_loadu_pd(&rty_P[i]);
    (void)_mm_storeu_pd(&rty_P[i], _mm_div_pd(_mm_sub_pd(tmp_1, tmp_2),
      _mm_set1_pd(rtu_lambda)));
  }

  for (i = 0; i <= 10; i += 2) {
    // Delay: '<S160>/Delay'
    tmp_1 = _mm_loadu_pd(&localDW->Delay_DSTATE[i]);

    // Sum: '<S160>/Sum' incorporates:
    //   Delay: '<S160>/Delay'

    tmp_2 = _mm_loadu_pd(&rty_theta[i]);
    (void)_mm_storeu_pd(&rty_theta[i], _mm_add_pd(tmp_1, tmp_2));
  }

  // Update for Delay: '<S160>/Delay1'
  localDW->icLoad = false;
  (void)std::memcpy(&localDW->Delay1_DSTATE[0], &rty_P[0], 144U * sizeof(real_T));

  // Update for UnitDelay: '<S158>/Unit Delay3' incorporates:
  //   Sum: '<S158>/Add1'

  localDW->UnitDelay3_DSTATE[0] = rtb_Add1_c[0];
  localDW->UnitDelay3_DSTATE[1] = rtb_Add1_c[1];
  localDW->UnitDelay3_DSTATE[2] = rtb_Add1_c[2];

  // Update for Delay: '<S160>/Delay'
  localDW->icLoad_n = false;
  (void)std::memcpy(&localDW->Delay_DSTATE[0], &rty_theta[0], 12U * sizeof
                    (real_T));
}

//
// Output and update for atomic system:
//    '<S2>/MATLAB Function'
//    '<S3>/MATLAB Function'
//
void SupervisoryController::MATLABFunction(const real_T rtu_yDest[6], real_T
  rtu_k_2, real_T rty_ywtT[6], real_T rty_ywt[6], real_T rty_uwt[3], P *rtP)
{
  real_T dwt;
  real_T sigmoid_workspace_k_1;
  real_T sigmoid_workspace_x0;

  // MATLAB Function 'SupervisoryController/ampc/MATLAB Function': '<S7>:1'
  // '<S7>:1:2' [ywt, ywtT, uwt, uwtT] = wtMod_(y, yDest, ywtT, uwtT, dt, no, ni, k_2); 
  // 'wtMod_:3' ywt = zeros(1,2*no);
  // 'wtMod_:4' uwt = dt*ones(1,ni);
  rty_uwt[0] = rtP->dt;
  rty_uwt[1] = rtP->dt;
  rty_uwt[2] = rtP->dt;

  //  time-scaled sigmoid (around 0->1 in 0->k_2 seconds, k_2 being a time scale factor) 
  // 'wtMod_:7' r = 0.2;
  //  10% to 90% rise time
  // 'wtMod_:8' dwt = dt/k_2;
  dwt = rtP->dt / rtu_k_2;

  // 'wtMod_:9' k_1 = 2.197/(r*k_2);
  sigmoid_workspace_k_1 = 2.197 / (0.2 * rtu_k_2);

  // 'wtMod_:10' x0 = 0.5*k_2;
  sigmoid_workspace_x0 = 0.5 * rtu_k_2;

  //  midpoint
  // 'wtMod_:11' sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));
  // 'wtMod_:13' for i = 1:2*no
  for (int32_T i{0}; i < 6; i++) {
    // 'wtMod_:14' if yDest(i) ~= 0
    if (rtu_yDest[i] != 0.0) {
      //  drive ywt to 1
      // 'wtMod_:16' if (ywtT(i) <= 1)
      if (rty_ywtT[i] <= 1.0) {
        // 'wtMod_:17' ywtT(i) = ywtT(i) + dwt;
        rty_ywtT[i] += dwt;
      }

      // 'wtMod_:19' else
      //  drive ywt to 0
      // 'wtMod_:21' if (ywtT(i) > 0)
    } else if (rty_ywtT[i] > 0.0) {
      // 'wtMod_:22' ywtT(i) = ywtT(i) - dwt;
      rty_ywtT[i] -= dwt;
    } else {
      // no actions
    }

    // 'wtMod_:25' if ywtT(i) <= 0
    if (rty_ywtT[i] <= 0.0) {
      // 'wtMod_:26' ywt(i) = 0;
      rty_ywt[i] = 0.0;
    } else {
      // 'wtMod_:27' else
      // 'wtMod_:28' ywt(i) = sigmoid(ywtT(i)*k_2);
      // 'wtMod_:11' @(x) 1/(1 + exp(-k_1*(x-x0)))
      rty_ywt[i] = 1.0 / (std::exp((rty_ywtT[i] * rtu_k_2 - sigmoid_workspace_x0)
        * -sigmoid_workspace_k_1) + 1.0);
    }
  }

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

// Function for Chart: '<Root>/SupervisoryController'
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

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
void SupervisoryController::mpc_checkhessian(real_T b_H[16], real_T b_L[16],
  real_T *BadH)
{
  real_T varargin_1[4];
  int32_T Tries;
  int8_T b[16];
  boolean_T guard1{ false };

  *BadH = 0.0;
  (void)std::memcpy(&b_L[0], &b_H[0], sizeof(real_T) << 4UL);
  Tries = xpotrf(b_L);
  guard1 = false;
  if (Tries == 0) {
    varargin_1[0] = b_L[0];
    varargin_1[1] = b_L[5];
    varargin_1[2] = b_L[10];
    varargin_1[3] = b_L[15];
    if (minimum(varargin_1) > 1.4901161193847656E-7) {
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    real_T normH;
    boolean_T exitg2;
    normH = 0.0;
    Tries = 0;
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (Tries < 4)) {
      real_T s;
      s = ((std::abs(b_H[Tries + 4]) + std::abs(b_H[Tries])) + std::abs
           (b_H[Tries + 8])) + std::abs(b_H[Tries + 12]);
      if (std::isnan(s)) {
        normH = (rtNaN);
        exitg2 = true;
      } else {
        if (s > normH) {
          normH = s;
        }

        Tries++;
      }
    }

    if (normH >= 1.0E+10) {
      *BadH = 2.0;
    } else {
      boolean_T exitg1;
      Tries = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (Tries <= 4)) {
        int32_T j;
        boolean_T guard2{ false };

        normH = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
          1.4901161193847656E-7;
        for (j = 0; j < 16; j++) {
          b[j] = 0;
        }

        b[0] = 1;
        b[5] = 1;
        b[10] = 1;
        b[15] = 1;
        for (j = 0; j < 16; j++) {
          b_H[j] += normH * static_cast<real_T>(b[j]);
          b_L[j] = b_H[j];
        }

        j = xpotrf(b_L);
        guard2 = false;
        if (j == 0) {
          varargin_1[0] = b_L[0];
          varargin_1[1] = b_L[5];
          varargin_1[2] = b_L[10];
          varargin_1[3] = b_L[15];
          if (minimum(varargin_1) > 1.4901161193847656E-7) {
            *BadH = 1.0;
            exitg1 = true;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          *BadH = 3.0;
          Tries++;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
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
  real_T a_c;
  real_T b;
  real_T y;
  a_c = std::abs(u0);
  b = std::abs(u1);
  if (a_c < b) {
    a_c /= b;
    y = std::sqrt(a_c * a_c + 1.0) * b;
  } else if (a_c > b) {
    b /= a_c;
    y = std::sqrt(b * b + 1.0) * a_c;
  } else if (std::isnan(b)) {
    y = (rtNaN);
  } else {
    y = a_c * 1.4142135623730951;
  }

  return y;
}

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
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

// Function for MATLAB Function: '<S112>/optimizer'
void SupervisoryController::KWIKfactor(const real_T b_Ac[984], const int32_T iC
  [246], int32_T nA, const real_T b_Linv[16], real_T D[16], real_T b_H[16],
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
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 245];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 491];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 737];
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

// Function for MATLAB Function: '<S112>/optimizer'
void SupervisoryController::DropConstraint(int32_T kDrop, boolean_T iA[246],
  int32_T *nA, int32_T iC[246])
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

// Function for MATLAB Function: '<S112>/optimizer'
void SupervisoryController::qpkwik(const real_T b_Linv[16], const real_T b_Hinv
  [16], const real_T f[4], const real_T b_Ac[984], const real_T b[246],
  boolean_T iA[246], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
  lambda[246], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[246];
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
  int32_T iC[246];
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
  for (i = 0; i < 246; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 246; tmp++) {
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
            (void)std::memset(&iA[0], 0, 246U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 246U * sizeof(int32_T));
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
              (void)std::memset(&iA[0], 0, 246U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 246U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 246U * sizeof(real_T));
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
      for (i = 0; i < 246; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 246] * x[1]);
          z[2] = std::abs(b_Ac[i + 492] * x[2]);
          z[3] = std::abs(b_Ac[i + 738] * x[3]);
          t = std::fmax(t, maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 246] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 492] * x
                    [2]) + b_Ac[i + 738] * x[3]) - b[i]) / t;
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
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 738])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 492])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 246])),
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
                    _mm_set1_pd(b_Ac[tmp + 738])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 492])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 246])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((D[iSave + 1] * b_Ac[tmp + 246] + D[iSave] * b_Ac[tmp])
                          + D[iSave + 2] * b_Ac[tmp + 492]) + D[iSave + 3] *
                    b_Ac[tmp + 738];
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

              t = b_Ac[tmp + 246];
              cVal_tmp = b_Ac[tmp + 492];
              cVal_tmp_0 = b_Ac[tmp + 738];
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
                  if ((iC_0 <= 246) && (lambda[iC_0 - 1] < 0.0)) {
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
              for (tmp = 0; tmp < 246; tmp++) {
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

// Function for MATLAB Function: '<S112>/optimizer'
void SupervisoryController::mpcblock_optimizer(const real_T rseq[120], const
  real_T vseq[21], const real_T umax[3], const real_T ymin[6], const real_T
  ymax[6], int32_T switch_in, const real_T x[8], const real_T old_u[3], const
  boolean_T iA[246], const real_T b_Mlim[246], const real_T b_Mx[1968], const
  real_T b_Mu1[738], const real_T b_Mv[5166], const real_T b_utarget[60], const
  real_T b_uoff[3], const real_T b_yoff[6], int32_T b_enable_value, real_T b_H
  [16], const real_T b_Ac[984], const real_T ywt[6], const real_T uwt[3], const
  real_T b_Wdu[3], const real_T b_Jm[180], const real_T b_SuJm[360], const
  real_T b_Su1[360], const real_T b_Sx[960], const real_T b_Hv[2520], const
  real_T b_I1[180], const int32_T b_Mrows[246], const real_T b_RYscale[6], const
  real_T b_RMVscale[3], real_T u[3], real_T *cost, real_T useq[63], real_T
  *status, boolean_T iAout[246])
{
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

  real_T WySuJm[360];
  real_T Bc[246];
  real_T a__1[246];
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T aux2[120];
  real_T b_Sx_0[120];
  real_T b_Kv[63];
  real_T K[60];
  real_T a_1[60];
  real_T aux3[60];
  real_T b_Kx[24];
  real_T b_Linv[16];
  real_T c_Linv[16];
  real_T B_0[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T b_Wy[6];
  real_T ymax_incr[6];
  real_T ymin_incr[6];
  real_T aux[4];
  real_T f[4];
  real_T zopt[4];
  real_T b_Wu[3];
  real_T umax_incr[3];
  real_T ywt_0;
  int32_T kidx;
  int8_T a_0[3600];
  int8_T b_B[16];
  boolean_T ymax_incr_flag[6];
  boolean_T ymin_incr_flag[6];
  boolean_T umax_incr_flag[3];
  *cost = 0.0;
  *status = 1.0;
  (void)std::memset(&iAout[0], 0, 246U * sizeof(boolean_T));
  if (switch_in != b_enable_value) {
    u[0] = old_u[0] + b_uoff[0];
    u[1] = old_u[1] + b_uoff[1];
    u[2] = old_u[2] + b_uoff[2];
    for (int32_T i{0}; i < 21; i++) {
      useq[i] = u[0];
      useq[i + 21] = u[1];
      useq[i + 42] = u[2];
    }
  } else {
    __m128d tmp;
    __m128d tmp_0;
    int32_T a_tmp;
    int32_T b_Kx_tmp;
    int32_T b_Linv_tmp;
    int32_T i;
    int32_T i1;
    int32_T j2;
    int16_T ixw;
    for (kidx = 0; kidx < 6; kidx++) {
      ywt_0 = ywt[kidx];
      if (ywt_0 < 0.0) {
        b_Wy[kidx] = 0.0;
      } else {
        b_Wy[kidx] = ywt_0 * ywt_0;
      }
    }

    if (uwt[0] < 0.0) {
      b_Wu[0] = 0.0;
    } else {
      b_Wu[0] = uwt[0] * uwt[0];
    }

    if (uwt[1] < 0.0) {
      b_Wu[1] = 0.0;
    } else {
      b_Wu[1] = uwt[1] * uwt[1];
    }

    if (uwt[2] < 0.0) {
      b_Wu[2] = 0.0;
    } else {
      b_Wu[2] = uwt[2] * uwt[2];
    }

    (void)std::memset(&B_0[0], 0, 9U * sizeof(real_T));
    B_0[0] = 1.0;
    B_0[4] = 1.0;
    B_0[8] = 1.0;
    kidx = -1;
    for (i = 0; i < 20; i++) {
      for (j2 = 0; j2 < 3; j2++) {
        for (i1 = 0; i1 < 20; i1++) {
          a_tmp = static_cast<int32_T>(b_A[20 * i + i1]);
          a_0[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2]) *
            a_tmp);
          a_0[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2 +
            1]) * a_tmp);
          a_0[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2 +
            2]) * a_tmp);
          kidx += 3;
        }
      }
    }

    i = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      for (j2 = 0; j2 < 60; j2++) {
        a_tmp = j2 + i;
        I2Jm[a_tmp] = 0.0;
        i1 = 0;
        for (b_Linv_tmp = 0; b_Linv_tmp < 60; b_Linv_tmp++) {
          I2Jm[a_tmp] += static_cast<real_T>(a_0[i1 + j2]) * b_Jm[b_Linv_tmp + i];
          i1 += 60;
        }
      }

      i += 60;
    }

    ixw = 1;
    for (kidx = 0; kidx < 120; kidx++) {
      ywt_0 = b_Wy[ixw - 1];
      WySuJm[kidx] = ywt_0 * b_SuJm[kidx];
      WySuJm[kidx + 120] = b_SuJm[kidx + 120] * ywt_0;
      WySuJm[kidx + 240] = b_SuJm[kidx + 240] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 6) {
        ixw = 1;
      }
    }

    ixw = 1;
    for (kidx = 0; kidx < 60; kidx++) {
      ywt_0 = b_Wu[ixw - 1];
      WuI2Jm[kidx] = ywt_0 * I2Jm[kidx];
      WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * ywt_0;
      WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 3) {
        ixw = 1;
      }
    }

    ixw = 1;
    for (kidx = 0; kidx < 60; kidx++) {
      ywt_0 = b_Wdu[ixw - 1];
      WduJm[kidx] = ywt_0 * b_Jm[kidx];
      WduJm[kidx + 60] = b_Jm[kidx + 60] * ywt_0;
      WduJm[kidx + 120] = b_Jm[kidx + 120] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 3) {
        ixw = 1;
      }
    }

    for (i = 0; i < 3; i++) {
      for (kidx = 0; kidx < 3; kidx++) {
        a_tmp = 3 * kidx + i;
        B_0[a_tmp] = 0.0;
        for (j2 = 0; j2 < 120; j2++) {
          B_0[a_tmp] += b_SuJm[120 * i + j2] * WySuJm[120 * kidx + j2];
        }

        b_Jm_0[a_tmp] = 0.0;
        ywt_0 = 0.0;
        for (j2 = 0; j2 < 60; j2++) {
          i1 = 60 * i + j2;
          b_Linv_tmp = 60 * kidx + j2;
          ywt_0 += I2Jm[i1] * WuI2Jm[b_Linv_tmp];
          b_Jm_0[a_tmp] += b_Jm[i1] * WduJm[b_Linv_tmp];
        }

        b_H[i + (kidx << 2UL)] = (B_0[a_tmp] + b_Jm_0[a_tmp]) + ywt_0;
      }
    }

    for (i = 0; i < 3; i++) {
      for (kidx = 0; kidx < 3; kidx++) {
        i1 = 3 * kidx + i;
        b_Jm_0[i1] = 0.0;
        for (j2 = 0; j2 < 120; j2++) {
          b_Jm_0[i1] += b_Su1[120 * i + j2] * WySuJm[120 * kidx + j2];
        }

        b_I1_0[i1] = 0.0;
        for (j2 = 0; j2 < 60; j2++) {
          b_I1_0[i1] += b_I1[60 * i + j2] * WuI2Jm[60 * kidx + j2];
        }
      }
    }

    for (i = 0; i <= 6; i += 2) {
      tmp = _mm_loadu_pd(&b_Jm_0[i]);
      tmp_0 = _mm_loadu_pd(&b_I1_0[i]);
      (void)_mm_storeu_pd(&B_0[i], _mm_add_pd(tmp, tmp_0));
    }

    for (i = 8; i < 9; i++) {
      B_0[i] = b_Jm_0[i] + b_I1_0[i];
    }

    for (i = 0; i <= 178; i += 2) {
      tmp = _mm_loadu_pd(&WuI2Jm[i]);
      (void)_mm_storeu_pd(&WuI2Jm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    i = 0;
    for (kidx = 0; kidx < 8; kidx++) {
      j2 = 0;
      i1 = 0;
      for (b_Linv_tmp = 0; b_Linv_tmp < 3; b_Linv_tmp++) {
        b_Kx_tmp = j2 + kidx;
        b_Kx[b_Kx_tmp] = 0.0;
        for (a_tmp = 0; a_tmp < 120; a_tmp++) {
          b_Kx[b_Kx_tmp] += b_Sx[a_tmp + i] * WySuJm[a_tmp + i1];
        }

        j2 += 8;
        i1 += 120;
      }

      i += 120;
    }

    i = 0;
    for (kidx = 0; kidx < 21; kidx++) {
      j2 = 0;
      i1 = 0;
      for (b_Linv_tmp = 0; b_Linv_tmp < 3; b_Linv_tmp++) {
        b_Kx_tmp = j2 + kidx;
        b_Kv[b_Kx_tmp] = 0.0;
        for (a_tmp = 0; a_tmp < 120; a_tmp++) {
          b_Kv[b_Kx_tmp] += b_Hv[a_tmp + i] * WySuJm[a_tmp + i1];
        }

        j2 += 21;
        i1 += 120;
      }

      i += 120;
    }

    for (i = 0; i <= 358; i += 2) {
      tmp = _mm_loadu_pd(&WySuJm[i]);
      (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    (void)std::memcpy(&b_Linv[0], &b_H[0], sizeof(real_T) << 4UL);
    mpc_checkhessian(b_Linv, c_Linv, &ywt_0);
    if (ywt_0 > 1.0) {
      u[0] = old_u[0] + b_uoff[0];
      u[1] = old_u[1] + b_uoff[1];
      u[2] = old_u[2] + b_uoff[2];
      for (i = 0; i < 21; i++) {
        useq[i] = u[0];
        useq[i + 21] = u[1];
        useq[i + 42] = u[2];
      }

      *status = -2.0;
    } else {
      real_T WySuJm_0;
      real_T b_Kv_0;
      for (i = 0; i < 16; i++) {
        b_B[i] = 0;
      }

      b_B[0] = 1;
      b_B[5] = 1;
      b_B[10] = 1;
      b_B[15] = 1;
      i = 0;
      for (kidx = 0; kidx < 4; kidx++) {
        b_Linv[i] = static_cast<real_T>(b_B[i]);
        b_Linv[i + 1] = static_cast<real_T>(b_B[i + 1]);
        b_Linv[i + 2] = static_cast<real_T>(b_B[i + 2]);
        b_Linv[i + 3] = static_cast<real_T>(b_B[i + 3]);
        i += 4;
      }

      trisolve(c_Linv, b_Linv);
      for (i = 0; i < 246; i++) {
        ywt_0 = 0.0;
        for (kidx = 0; kidx < 8; kidx++) {
          ywt_0 += b_Mx[246 * kidx + i] * x[kidx];
        }

        WySuJm_0 = 0.0;
        for (kidx = 0; kidx < 21; kidx++) {
          WySuJm_0 += b_Mv[246 * kidx + i] * vseq[kidx];
        }

        Bc[i] = -((((b_Mu1[i + 246] * old_u[1] + b_Mu1[i] * old_u[0]) + b_Mu1[i
                    + 492] * old_u[2]) + (b_Mlim[i] + ywt_0)) + WySuJm_0);
      }

      for (i = 0; i < 6; i++) {
        ymax_incr_flag[i] = false;
        ymax_incr[i] = 0.0;
        ymin_incr_flag[i] = false;
        ymin_incr[i] = 0.0;
      }

      umax_incr_flag[0] = false;
      umax_incr[0] = 0.0;
      umax_incr_flag[1] = false;
      umax_incr[1] = 0.0;
      umax_incr_flag[2] = false;
      umax_incr[2] = 0.0;
      if (b_Mrows[0] > 0) {
        boolean_T exitg1;
        kidx = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (kidx < 246)) {
          if (b_Mrows[kidx] <= 120) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 1,
                  static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 1;
            c_Del_Save_Flag0 = ymax_incr_flag[i];
            if (!ymax_incr_flag[i]) {
              ywt_0 = -(b_RYscale[i] * ymax[i] - b_yoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = ymax_incr[i];
            }

            ymax_incr[i] = ywt_0;
            ymax_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 240) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 121,
                  static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 121;
            c_Del_Save_Flag0 = ymin_incr_flag[i];
            if (!ymin_incr_flag[i]) {
              ywt_0 = (b_RYscale[i] * ymin[i] - b_yoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = ymin_incr[i];
            }

            ymin_incr[i] = ywt_0;
            ymin_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 300) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 241,
                  static_cast<int32_T>(nu)) * static_cast<int32_T>(nu)) - 241;
            c_Del_Save_Flag0 = umax_incr_flag[i];
            if (!umax_incr_flag[i]) {
              ywt_0 = -(b_RMVscale[i] * umax[i] - b_uoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = umax_incr[i];
            }

            umax_incr[i] = ywt_0;
            umax_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 360) {
            kidx++;
          } else {
            exitg1 = true;
          }
        }
      }

      f[0] = 0.0;
      f[1] = 0.0;
      f[2] = 0.0;
      f[3] = 0.0;
      for (kidx = 0; kidx < 3; kidx++) {
        real_T WuI2Jm_0;
        ywt_0 = 0.0;
        for (i = 0; i < 8; i++) {
          ywt_0 += b_Kx[(kidx << 3UL) + i] * x[i];
        }

        WySuJm_0 = 0.0;
        for (i = 0; i < 120; i++) {
          WySuJm_0 += WySuJm[120 * kidx + i] * rseq[i];
        }

        b_Kv_0 = 0.0;
        for (i = 0; i < 21; i++) {
          b_Kv_0 += b_Kv[21 * kidx + i] * vseq[i];
        }

        WuI2Jm_0 = 0.0;
        for (i = 0; i < 60; i++) {
          WuI2Jm_0 += WuI2Jm[60 * kidx + i] * b_utarget[i];
        }

        f[kidx] = ((((B_0[3 * kidx + 1] * old_u[1] + B_0[3 * kidx] * old_u[0]) +
                     B_0[3 * kidx + 2] * old_u[2]) + (ywt_0 + WySuJm_0)) +
                   b_Kv_0) + WuI2Jm_0;
      }

      (void)std::memcpy(&iAout[0], &iA[0], 246U * sizeof(boolean_T));
      i = 0;
      for (kidx = 0; kidx < 4; kidx++) {
        j2 = 0;
        for (i1 = 0; i1 < 4; i1++) {
          b_Linv_tmp = j2 + kidx;
          c_Linv[b_Linv_tmp] = 0.0;
          c_Linv[b_Linv_tmp] += b_Linv[i] * b_Linv[j2];
          c_Linv[b_Linv_tmp] += b_Linv[i + 1] * b_Linv[j2 + 1];
          c_Linv[b_Linv_tmp] += b_Linv[i + 2] * b_Linv[j2 + 2];
          c_Linv[b_Linv_tmp] += b_Linv[i + 3] * b_Linv[j2 + 3];
          j2 += 4;
        }

        i += 4;
      }

      qpkwik(b_Linv, c_Linv, f, b_Ac, Bc, iAout, 1000, 1.0E-6, zopt, a__1, &kidx);
      if ((kidx < 0) || (kidx == 0)) {
        zopt[0] = 0.0;
        zopt[1] = 0.0;
        zopt[2] = 0.0;
        zopt[3] = 0.0;
      }

      *status = static_cast<real_T>(kidx);
      u[0] = (old_u[0] + zopt[0]) + b_uoff[0];
      u[1] = (old_u[1] + zopt[1]) + b_uoff[1];
      u[2] = (old_u[2] + zopt[2]) + b_uoff[2];
      if (kidx > 0) {
        for (i = 0; i < 120; i++) {
          b_Sx_0[i] = 0.0;
          for (kidx = 0; kidx < 8; kidx++) {
            b_Sx_0[i] += b_Sx[120 * kidx + i] * x[kidx];
          }

          ywt_0 = 0.0;
          for (kidx = 0; kidx < 21; kidx++) {
            ywt_0 += b_Hv[120 * kidx + i] * vseq[kidx];
          }

          aux2[i] = ((((b_Su1[i + 120] * old_u[1] + b_Su1[i] * old_u[0]) +
                       b_Su1[i + 240] * old_u[2]) + b_Sx_0[i]) + ywt_0) - rseq[i];
        }

        for (i = 0; i <= 58; i += 2) {
          (void)_mm_storeu_pd(&aux3[i], _mm_sub_pd(_mm_add_pd(_mm_add_pd
            (_mm_mul_pd(_mm_loadu_pd(&b_I1[i + 60]), _mm_set1_pd(old_u[1])),
             _mm_mul_pd(_mm_loadu_pd(&b_I1[i]), _mm_set1_pd(old_u[0]))),
            _mm_mul_pd(_mm_loadu_pd(&b_I1[i + 120]), _mm_set1_pd(old_u[2]))),
            _mm_loadu_pd(&b_utarget[i])));
        }

        kidx = -1;
        for (i = 0; i < 20; i++) {
          K[kidx + 1] = b_Wu[0];
          K[kidx + 2] = b_Wu[1];
          K[kidx + 3] = b_Wu[2];
          kidx += 3;
        }

        kidx = -1;
        for (i = 0; i < 20; i++) {
          for (j2 = 0; j2 < 6; j2++) {
            b_Sx_0[(kidx + j2) + 1] = b_Wy[j2];
          }

          kidx += 6;
        }

        aux[0] = zopt[0];
        aux[1] = zopt[1];
        aux[2] = zopt[2];
        aux[3] = zopt[3];
        ywt_0 = 0.0;
        for (i = 0; i < 60; i++) {
          WySuJm_0 = aux3[i];
          ywt_0 += K[i] * WySuJm_0 * WySuJm_0;
        }

        WySuJm_0 = 0.0;
        for (i = 0; i < 120; i++) {
          b_Kv_0 = aux2[i];
          WySuJm_0 += b_Sx_0[i] * b_Kv_0 * b_Kv_0;
        }

        b_Kv_0 = 0.0;
        for (i = 0; i < 4; i++) {
          b_Kv_0 += ((((b_H[i + 4] * aux[1] + b_H[i] * aux[0]) + b_H[i + 8] *
                       aux[2]) + b_H[i + 12] * aux[3]) + 2.0 * f[i]) * aux[i];
        }

        *cost = (ywt_0 + WySuJm_0) + b_Kv_0;
      }

      (void)std::memset(&B_0[0], 0, 9U * sizeof(real_T));
      B_0[0] = 1.0;
      B_0[4] = 1.0;
      B_0[8] = 1.0;
      kidx = -1;
      for (i = 0; i < 20; i++) {
        for (j2 = 0; j2 < 3; j2++) {
          for (i1 = 0; i1 < 20; i1++) {
            a_tmp = static_cast<int32_T>(b_A[20 * i + i1]);
            a_0[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2])
              * a_tmp);
            a_0[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2
              + 1]) * a_tmp);
            a_0[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(B_0[3 * j2
              + 2]) * a_tmp);
            kidx += 3;
          }
        }
      }

      i = 0;
      for (kidx = 0; kidx < 3; kidx++) {
        for (j2 = 0; j2 < 60; j2++) {
          a_tmp = j2 + i;
          I2Jm[a_tmp] = 0.0;
          i1 = 0;
          for (b_Linv_tmp = 0; b_Linv_tmp < 60; b_Linv_tmp++) {
            I2Jm[a_tmp] += static_cast<real_T>(a_0[i1 + j2]) * b_Jm[b_Linv_tmp +
              i];
            i1 += 60;
          }
        }

        i += 60;
      }

      ywt_0 = old_u[0] + b_uoff[0];
      WySuJm_0 = old_u[1] + b_uoff[1];
      b_Kv_0 = old_u[2] + b_uoff[2];
      for (i = 0; i <= 58; i += 2) {
        tmp = _mm_set1_pd(0.0);
        (void)_mm_storeu_pd(&aux3[i], tmp);
        (void)_mm_storeu_pd(&K[i], tmp);
        tmp = _mm_loadu_pd(&I2Jm[i]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[0])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i]),
          _mm_set1_pd(ywt_0)), tmp));
        tmp = _mm_loadu_pd(&I2Jm[i + 60]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[1])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i +
          60]), _mm_set1_pd(WySuJm_0)), tmp));
        tmp = _mm_loadu_pd(&I2Jm[i + 120]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[2])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i +
          120]), _mm_set1_pd(b_Kv_0)), tmp));
        tmp = _mm_loadu_pd(&aux3[i]);
        tmp_0 = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&a_1[i], _mm_add_pd(tmp, tmp_0));
      }

      i = 0;
      for (kidx = 0; kidx < 20; kidx++) {
        useq[kidx] = a_1[i];
        useq[kidx + 21] = a_1[i + 1];
        useq[kidx + 42] = a_1[i + 2];
        i += 3;
      }

      useq[20] = useq[19];
      useq[41] = useq[40];
      useq[62] = useq[61];
    }
  }
}

// Function for MATLAB Function: '<S134>/optimizer'
void SupervisoryController::mpcblock_optimizer_h(const real_T rseq[120], const
  real_T vseq[21], const real_T umax[3], const real_T ymin[6], const real_T
  ymax[6], int32_T switch_in, const real_T x[10], const real_T old_u[3], const
  boolean_T iA[246], const real_T b_Mlim[246], const real_T b_Mx[2460], const
  real_T b_Mu1[738], const real_T b_Mv[5166], const real_T b_utarget[60], const
  real_T b_uoff[3], const real_T b_yoff[6], int32_T b_enable_value, real_T b_H
  [16], const real_T b_Ac[984], const real_T ywt[6], const real_T uwt[3], const
  real_T b_Wdu[3], const real_T b_Jm[180], const real_T b_SuJm[360], const
  real_T b_Su1[360], const real_T b_Sx[1200], const real_T b_Hv[2520], const
  real_T b_I1[180], const int32_T b_Mrows[246], const real_T b_RYscale[6], const
  real_T b_RMVscale[3], real_T u[3], real_T *cost, real_T useq[63], real_T
  *status, boolean_T iAout[246])
{
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

  real_T WySuJm[360];
  real_T Bc[246];
  real_T a__1[246];
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T aux2[120];
  real_T b_Sx_0[120];
  real_T b_Kv[63];
  real_T K[60];
  real_T a_3[60];
  real_T aux3[60];
  real_T b_Kx[30];
  real_T b_Linv[16];
  real_T c_Linv[16];
  real_T B_2[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T b_Wy[6];
  real_T ymax_incr[6];
  real_T ymin_incr[6];
  real_T aux[4];
  real_T f[4];
  real_T zopt[4];
  real_T b_Wu[3];
  real_T umax_incr[3];
  real_T ywt_0;
  int32_T kidx;
  int8_T a_2[3600];
  int8_T b_B[16];
  boolean_T ymax_incr_flag[6];
  boolean_T ymin_incr_flag[6];
  boolean_T umax_incr_flag[3];
  *cost = 0.0;
  *status = 1.0;
  (void)std::memset(&iAout[0], 0, 246U * sizeof(boolean_T));
  if (switch_in != b_enable_value) {
    u[0] = old_u[0] + b_uoff[0];
    u[1] = old_u[1] + b_uoff[1];
    u[2] = old_u[2] + b_uoff[2];
    for (int32_T i{0}; i < 21; i++) {
      useq[i] = u[0];
      useq[i + 21] = u[1];
      useq[i + 42] = u[2];
    }
  } else {
    __m128d tmp;
    __m128d tmp_0;
    int32_T a_tmp;
    int32_T b_Kx_tmp;
    int32_T b_Linv_tmp;
    int32_T i;
    int32_T i1;
    int32_T j2;
    int16_T ixw;
    for (kidx = 0; kidx < 6; kidx++) {
      ywt_0 = ywt[kidx];
      if (ywt_0 < 0.0) {
        b_Wy[kidx] = 0.0;
      } else {
        b_Wy[kidx] = ywt_0 * ywt_0;
      }
    }

    if (uwt[0] < 0.0) {
      b_Wu[0] = 0.0;
    } else {
      b_Wu[0] = uwt[0] * uwt[0];
    }

    if (uwt[1] < 0.0) {
      b_Wu[1] = 0.0;
    } else {
      b_Wu[1] = uwt[1] * uwt[1];
    }

    if (uwt[2] < 0.0) {
      b_Wu[2] = 0.0;
    } else {
      b_Wu[2] = uwt[2] * uwt[2];
    }

    (void)std::memset(&B_2[0], 0, 9U * sizeof(real_T));
    B_2[0] = 1.0;
    B_2[4] = 1.0;
    B_2[8] = 1.0;
    kidx = -1;
    for (i = 0; i < 20; i++) {
      for (j2 = 0; j2 < 3; j2++) {
        for (i1 = 0; i1 < 20; i1++) {
          a_tmp = static_cast<int32_T>(b_A[20 * i + i1]);
          a_2[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2]) *
            a_tmp);
          a_2[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2 +
            1]) * a_tmp);
          a_2[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2 +
            2]) * a_tmp);
          kidx += 3;
        }
      }
    }

    i = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      for (j2 = 0; j2 < 60; j2++) {
        a_tmp = j2 + i;
        I2Jm[a_tmp] = 0.0;
        i1 = 0;
        for (b_Linv_tmp = 0; b_Linv_tmp < 60; b_Linv_tmp++) {
          I2Jm[a_tmp] += static_cast<real_T>(a_2[i1 + j2]) * b_Jm[b_Linv_tmp + i];
          i1 += 60;
        }
      }

      i += 60;
    }

    ixw = 1;
    for (kidx = 0; kidx < 120; kidx++) {
      ywt_0 = b_Wy[ixw - 1];
      WySuJm[kidx] = ywt_0 * b_SuJm[kidx];
      WySuJm[kidx + 120] = b_SuJm[kidx + 120] * ywt_0;
      WySuJm[kidx + 240] = b_SuJm[kidx + 240] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 6) {
        ixw = 1;
      }
    }

    ixw = 1;
    for (kidx = 0; kidx < 60; kidx++) {
      ywt_0 = b_Wu[ixw - 1];
      WuI2Jm[kidx] = ywt_0 * I2Jm[kidx];
      WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * ywt_0;
      WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 3) {
        ixw = 1;
      }
    }

    ixw = 1;
    for (kidx = 0; kidx < 60; kidx++) {
      ywt_0 = b_Wdu[ixw - 1];
      WduJm[kidx] = ywt_0 * b_Jm[kidx];
      WduJm[kidx + 60] = b_Jm[kidx + 60] * ywt_0;
      WduJm[kidx + 120] = b_Jm[kidx + 120] * ywt_0;
      ixw = static_cast<int16_T>(ixw + 1);
      if (ixw > 3) {
        ixw = 1;
      }
    }

    for (i = 0; i < 3; i++) {
      for (kidx = 0; kidx < 3; kidx++) {
        a_tmp = 3 * kidx + i;
        B_2[a_tmp] = 0.0;
        for (j2 = 0; j2 < 120; j2++) {
          B_2[a_tmp] += b_SuJm[120 * i + j2] * WySuJm[120 * kidx + j2];
        }

        b_Jm_0[a_tmp] = 0.0;
        ywt_0 = 0.0;
        for (j2 = 0; j2 < 60; j2++) {
          i1 = 60 * i + j2;
          b_Linv_tmp = 60 * kidx + j2;
          ywt_0 += I2Jm[i1] * WuI2Jm[b_Linv_tmp];
          b_Jm_0[a_tmp] += b_Jm[i1] * WduJm[b_Linv_tmp];
        }

        b_H[i + (kidx << 2UL)] = (B_2[a_tmp] + b_Jm_0[a_tmp]) + ywt_0;
      }
    }

    for (i = 0; i < 3; i++) {
      for (kidx = 0; kidx < 3; kidx++) {
        i1 = 3 * kidx + i;
        b_Jm_0[i1] = 0.0;
        for (j2 = 0; j2 < 120; j2++) {
          b_Jm_0[i1] += b_Su1[120 * i + j2] * WySuJm[120 * kidx + j2];
        }

        b_I1_0[i1] = 0.0;
        for (j2 = 0; j2 < 60; j2++) {
          b_I1_0[i1] += b_I1[60 * i + j2] * WuI2Jm[60 * kidx + j2];
        }
      }
    }

    for (i = 0; i <= 6; i += 2) {
      tmp = _mm_loadu_pd(&b_Jm_0[i]);
      tmp_0 = _mm_loadu_pd(&b_I1_0[i]);
      (void)_mm_storeu_pd(&B_2[i], _mm_add_pd(tmp, tmp_0));
    }

    for (i = 8; i < 9; i++) {
      B_2[i] = b_Jm_0[i] + b_I1_0[i];
    }

    for (i = 0; i <= 178; i += 2) {
      tmp = _mm_loadu_pd(&WuI2Jm[i]);
      (void)_mm_storeu_pd(&WuI2Jm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    i = 0;
    for (kidx = 0; kidx < 10; kidx++) {
      j2 = 0;
      i1 = 0;
      for (b_Linv_tmp = 0; b_Linv_tmp < 3; b_Linv_tmp++) {
        b_Kx_tmp = j2 + kidx;
        b_Kx[b_Kx_tmp] = 0.0;
        for (a_tmp = 0; a_tmp < 120; a_tmp++) {
          b_Kx[b_Kx_tmp] += b_Sx[a_tmp + i] * WySuJm[a_tmp + i1];
        }

        j2 += 10;
        i1 += 120;
      }

      i += 120;
    }

    i = 0;
    for (kidx = 0; kidx < 21; kidx++) {
      j2 = 0;
      i1 = 0;
      for (b_Linv_tmp = 0; b_Linv_tmp < 3; b_Linv_tmp++) {
        b_Kx_tmp = j2 + kidx;
        b_Kv[b_Kx_tmp] = 0.0;
        for (a_tmp = 0; a_tmp < 120; a_tmp++) {
          b_Kv[b_Kx_tmp] += b_Hv[a_tmp + i] * WySuJm[a_tmp + i1];
        }

        j2 += 21;
        i1 += 120;
      }

      i += 120;
    }

    for (i = 0; i <= 358; i += 2) {
      tmp = _mm_loadu_pd(&WySuJm[i]);
      (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    (void)std::memcpy(&b_Linv[0], &b_H[0], sizeof(real_T) << 4UL);
    mpc_checkhessian(b_Linv, c_Linv, &ywt_0);
    if (ywt_0 > 1.0) {
      u[0] = old_u[0] + b_uoff[0];
      u[1] = old_u[1] + b_uoff[1];
      u[2] = old_u[2] + b_uoff[2];
      for (i = 0; i < 21; i++) {
        useq[i] = u[0];
        useq[i + 21] = u[1];
        useq[i + 42] = u[2];
      }

      *status = -2.0;
    } else {
      real_T WySuJm_0;
      real_T b_Kv_0;
      for (i = 0; i < 16; i++) {
        b_B[i] = 0;
      }

      b_B[0] = 1;
      b_B[5] = 1;
      b_B[10] = 1;
      b_B[15] = 1;
      i = 0;
      for (kidx = 0; kidx < 4; kidx++) {
        b_Linv[i] = static_cast<real_T>(b_B[i]);
        b_Linv[i + 1] = static_cast<real_T>(b_B[i + 1]);
        b_Linv[i + 2] = static_cast<real_T>(b_B[i + 2]);
        b_Linv[i + 3] = static_cast<real_T>(b_B[i + 3]);
        i += 4;
      }

      trisolve(c_Linv, b_Linv);
      for (i = 0; i < 246; i++) {
        ywt_0 = 0.0;
        for (kidx = 0; kidx < 10; kidx++) {
          ywt_0 += b_Mx[246 * kidx + i] * x[kidx];
        }

        WySuJm_0 = 0.0;
        for (kidx = 0; kidx < 21; kidx++) {
          WySuJm_0 += b_Mv[246 * kidx + i] * vseq[kidx];
        }

        Bc[i] = -((((b_Mu1[i + 246] * old_u[1] + b_Mu1[i] * old_u[0]) + b_Mu1[i
                    + 492] * old_u[2]) + (b_Mlim[i] + ywt_0)) + WySuJm_0);
      }

      for (i = 0; i < 6; i++) {
        ymax_incr_flag[i] = false;
        ymax_incr[i] = 0.0;
        ymin_incr_flag[i] = false;
        ymin_incr[i] = 0.0;
      }

      umax_incr_flag[0] = false;
      umax_incr[0] = 0.0;
      umax_incr_flag[1] = false;
      umax_incr[1] = 0.0;
      umax_incr_flag[2] = false;
      umax_incr[2] = 0.0;
      if (b_Mrows[0] > 0) {
        boolean_T exitg1;
        kidx = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (kidx < 246)) {
          if (b_Mrows[kidx] <= 120) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 1,
                  static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 1;
            c_Del_Save_Flag0 = ymax_incr_flag[i];
            if (!ymax_incr_flag[i]) {
              ywt_0 = -(b_RYscale[i] * ymax[i] - b_yoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = ymax_incr[i];
            }

            ymax_incr[i] = ywt_0;
            ymax_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 240) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 121,
                  static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 121;
            c_Del_Save_Flag0 = ymin_incr_flag[i];
            if (!ymin_incr_flag[i]) {
              ywt_0 = (b_RYscale[i] * ymin[i] - b_yoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = ymin_incr[i];
            }

            ymin_incr[i] = ywt_0;
            ymin_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 300) {
            boolean_T c_Del_Save_Flag0;
            i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 241,
                  static_cast<int32_T>(nu)) * static_cast<int32_T>(nu)) - 241;
            c_Del_Save_Flag0 = umax_incr_flag[i];
            if (!umax_incr_flag[i]) {
              ywt_0 = -(b_RMVscale[i] * umax[i] - b_uoff[i]) - (-b_Mlim[kidx]);
              c_Del_Save_Flag0 = true;
            } else {
              ywt_0 = umax_incr[i];
            }

            umax_incr[i] = ywt_0;
            umax_incr_flag[i] = c_Del_Save_Flag0;
            Bc[kidx] += ywt_0;
            kidx++;
          } else if (b_Mrows[kidx] <= 360) {
            kidx++;
          } else {
            exitg1 = true;
          }
        }
      }

      f[0] = 0.0;
      f[1] = 0.0;
      f[2] = 0.0;
      f[3] = 0.0;
      for (kidx = 0; kidx < 3; kidx++) {
        real_T WuI2Jm_0;
        ywt_0 = 0.0;
        for (i = 0; i < 10; i++) {
          ywt_0 += b_Kx[10 * kidx + i] * x[i];
        }

        WySuJm_0 = 0.0;
        for (i = 0; i < 120; i++) {
          WySuJm_0 += WySuJm[120 * kidx + i] * rseq[i];
        }

        b_Kv_0 = 0.0;
        for (i = 0; i < 21; i++) {
          b_Kv_0 += b_Kv[21 * kidx + i] * vseq[i];
        }

        WuI2Jm_0 = 0.0;
        for (i = 0; i < 60; i++) {
          WuI2Jm_0 += WuI2Jm[60 * kidx + i] * b_utarget[i];
        }

        f[kidx] = ((((B_2[3 * kidx + 1] * old_u[1] + B_2[3 * kidx] * old_u[0]) +
                     B_2[3 * kidx + 2] * old_u[2]) + (ywt_0 + WySuJm_0)) +
                   b_Kv_0) + WuI2Jm_0;
      }

      (void)std::memcpy(&iAout[0], &iA[0], 246U * sizeof(boolean_T));
      i = 0;
      for (kidx = 0; kidx < 4; kidx++) {
        j2 = 0;
        for (i1 = 0; i1 < 4; i1++) {
          b_Linv_tmp = j2 + kidx;
          c_Linv[b_Linv_tmp] = 0.0;
          c_Linv[b_Linv_tmp] += b_Linv[i] * b_Linv[j2];
          c_Linv[b_Linv_tmp] += b_Linv[i + 1] * b_Linv[j2 + 1];
          c_Linv[b_Linv_tmp] += b_Linv[i + 2] * b_Linv[j2 + 2];
          c_Linv[b_Linv_tmp] += b_Linv[i + 3] * b_Linv[j2 + 3];
          j2 += 4;
        }

        i += 4;
      }

      qpkwik(b_Linv, c_Linv, f, b_Ac, Bc, iAout, 1000, 1.0E-6, zopt, a__1, &kidx);
      if ((kidx < 0) || (kidx == 0)) {
        zopt[0] = 0.0;
        zopt[1] = 0.0;
        zopt[2] = 0.0;
        zopt[3] = 0.0;
      }

      *status = static_cast<real_T>(kidx);
      u[0] = (old_u[0] + zopt[0]) + b_uoff[0];
      u[1] = (old_u[1] + zopt[1]) + b_uoff[1];
      u[2] = (old_u[2] + zopt[2]) + b_uoff[2];
      if (kidx > 0) {
        for (i = 0; i < 120; i++) {
          b_Sx_0[i] = 0.0;
          for (kidx = 0; kidx < 10; kidx++) {
            b_Sx_0[i] += b_Sx[120 * kidx + i] * x[kidx];
          }

          ywt_0 = 0.0;
          for (kidx = 0; kidx < 21; kidx++) {
            ywt_0 += b_Hv[120 * kidx + i] * vseq[kidx];
          }

          aux2[i] = ((((b_Su1[i + 120] * old_u[1] + b_Su1[i] * old_u[0]) +
                       b_Su1[i + 240] * old_u[2]) + b_Sx_0[i]) + ywt_0) - rseq[i];
        }

        for (i = 0; i <= 58; i += 2) {
          (void)_mm_storeu_pd(&aux3[i], _mm_sub_pd(_mm_add_pd(_mm_add_pd
            (_mm_mul_pd(_mm_loadu_pd(&b_I1[i + 60]), _mm_set1_pd(old_u[1])),
             _mm_mul_pd(_mm_loadu_pd(&b_I1[i]), _mm_set1_pd(old_u[0]))),
            _mm_mul_pd(_mm_loadu_pd(&b_I1[i + 120]), _mm_set1_pd(old_u[2]))),
            _mm_loadu_pd(&b_utarget[i])));
        }

        kidx = -1;
        for (i = 0; i < 20; i++) {
          K[kidx + 1] = b_Wu[0];
          K[kidx + 2] = b_Wu[1];
          K[kidx + 3] = b_Wu[2];
          kidx += 3;
        }

        kidx = -1;
        for (i = 0; i < 20; i++) {
          for (j2 = 0; j2 < 6; j2++) {
            b_Sx_0[(kidx + j2) + 1] = b_Wy[j2];
          }

          kidx += 6;
        }

        aux[0] = zopt[0];
        aux[1] = zopt[1];
        aux[2] = zopt[2];
        aux[3] = zopt[3];
        ywt_0 = 0.0;
        for (i = 0; i < 60; i++) {
          WySuJm_0 = aux3[i];
          ywt_0 += K[i] * WySuJm_0 * WySuJm_0;
        }

        WySuJm_0 = 0.0;
        for (i = 0; i < 120; i++) {
          b_Kv_0 = aux2[i];
          WySuJm_0 += b_Sx_0[i] * b_Kv_0 * b_Kv_0;
        }

        b_Kv_0 = 0.0;
        for (i = 0; i < 4; i++) {
          b_Kv_0 += ((((b_H[i + 4] * aux[1] + b_H[i] * aux[0]) + b_H[i + 8] *
                       aux[2]) + b_H[i + 12] * aux[3]) + 2.0 * f[i]) * aux[i];
        }

        *cost = (ywt_0 + WySuJm_0) + b_Kv_0;
      }

      (void)std::memset(&B_2[0], 0, 9U * sizeof(real_T));
      B_2[0] = 1.0;
      B_2[4] = 1.0;
      B_2[8] = 1.0;
      kidx = -1;
      for (i = 0; i < 20; i++) {
        for (j2 = 0; j2 < 3; j2++) {
          for (i1 = 0; i1 < 20; i1++) {
            a_tmp = static_cast<int32_T>(b_A[20 * i + i1]);
            a_2[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2])
              * a_tmp);
            a_2[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2
              + 1]) * a_tmp);
            a_2[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(B_2[3 * j2
              + 2]) * a_tmp);
            kidx += 3;
          }
        }
      }

      i = 0;
      for (kidx = 0; kidx < 3; kidx++) {
        for (j2 = 0; j2 < 60; j2++) {
          a_tmp = j2 + i;
          I2Jm[a_tmp] = 0.0;
          i1 = 0;
          for (b_Linv_tmp = 0; b_Linv_tmp < 60; b_Linv_tmp++) {
            I2Jm[a_tmp] += static_cast<real_T>(a_2[i1 + j2]) * b_Jm[b_Linv_tmp +
              i];
            i1 += 60;
          }
        }

        i += 60;
      }

      ywt_0 = old_u[0] + b_uoff[0];
      WySuJm_0 = old_u[1] + b_uoff[1];
      b_Kv_0 = old_u[2] + b_uoff[2];
      for (i = 0; i <= 58; i += 2) {
        tmp = _mm_set1_pd(0.0);
        (void)_mm_storeu_pd(&aux3[i], tmp);
        (void)_mm_storeu_pd(&K[i], tmp);
        tmp = _mm_loadu_pd(&I2Jm[i]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[0])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i]),
          _mm_set1_pd(ywt_0)), tmp));
        tmp = _mm_loadu_pd(&I2Jm[i + 60]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[1])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i +
          60]), _mm_set1_pd(WySuJm_0)), tmp));
        tmp = _mm_loadu_pd(&I2Jm[i + 120]);
        tmp_0 = _mm_loadu_pd(&aux3[i]);
        (void)_mm_storeu_pd(&aux3[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
          (zopt[2])), tmp_0));
        tmp = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&K[i], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_I1[i +
          120]), _mm_set1_pd(b_Kv_0)), tmp));
        tmp = _mm_loadu_pd(&aux3[i]);
        tmp_0 = _mm_loadu_pd(&K[i]);
        (void)_mm_storeu_pd(&a_3[i], _mm_add_pd(tmp, tmp_0));
      }

      i = 0;
      for (kidx = 0; kidx < 20; kidx++) {
        useq[kidx] = a_3[i];
        useq[kidx + 21] = a_3[i + 1];
        useq[kidx + 42] = a_3[i + 2];
        i += 3;
      }

      useq[20] = useq[19];
      useq[41] = useq[40];
      useq[62] = useq[61];
    }
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
  static const real_T d_0[2460]{ -0.0, -0.99801496201563555,
    -0.00012982817451900192, -0.0, -0.0, -0.0, -0.0, -0.99603367687094413,
    -0.00025924037690453866, -0.0, -0.0, -0.0, -0.0, -0.99405613771726609,
    -0.00038823760144369459, -0.0, -0.0, -0.0, -0.0, -0.99208233771810028,
    -0.00051682084032238652, -0.0, -0.0, -0.0, -0.0, -0.99011227004908309,
    -0.00064499108362950342, -0.0, -0.0, -0.0, -0.0, -0.9881459278979674,
    -0.00077274931936103868, -0.0, -0.0, -0.0, -0.0, -0.98618330446460145,
    -0.00090009653342421417, -0.0, -0.0, -0.0, -0.0, -0.984224392960908,
    -0.0010270337096415967, -0.0, -0.0, -0.0, -0.0, -0.98226918661086315,
    -0.001153561829755207, -0.0, -0.0, -0.0, -0.0, -0.98031767865047559,
    -0.0012796818734306208, -0.0, -0.0, -0.0, -0.0, -0.97836986232776568,
    -0.0014053948182610615, -0.0, -0.0, -0.0, -0.0, -0.97642573090274454,
    -0.0015307016397714866, -0.0, -0.0, -0.0, -0.0, -0.97448527764739334,
    -0.0016556033114226655, -0.0, -0.0, -0.0, -0.0, -0.97254849584564251,
    -0.0017801008046152487, -0.0, -0.0, -0.0, -0.0, -0.97061537879335091,
    -0.0019041950886938322, -0.0, -0.0, -0.0, -0.0, -0.968685919798285,
    -0.0020278871309510108, -0.0, -0.0, -0.0, -0.0, -0.96676011218009861,
    -0.0021511778966314264, -0.0, -0.0, -0.0, -0.0, -0.96483794927031175,
    -0.0022740683489358075, -0.0, -0.0, -0.0, -0.0, -0.9629194244122905,
    -0.0023965594490250011, -0.0, -0.0, -0.0, -0.0, -0.96100453096122618,
    -0.0025186521560239968, -0.0, -0.0, -0.0, 0.0, 0.99801496201563555,
    0.00012982817451900192, 0.0, 0.0, 0.0, 0.0, 0.99603367687094413,
    0.00025924037690453866, 0.0, 0.0, 0.0, 0.0, 0.99405613771726609,
    0.00038823760144369459, 0.0, 0.0, 0.0, 0.0, 0.99208233771810028,
    0.00051682084032238652, 0.0, 0.0, 0.0, 0.0, 0.99011227004908309,
    0.00064499108362950342, 0.0, 0.0, 0.0, 0.0, 0.9881459278979674,
    0.00077274931936103868, 0.0, 0.0, 0.0, 0.0, 0.98618330446460145,
    0.00090009653342421417, 0.0, 0.0, 0.0, 0.0, 0.984224392960908,
    0.0010270337096415967, 0.0, 0.0, 0.0, 0.0, 0.98226918661086315,
    0.001153561829755207, 0.0, 0.0, 0.0, 0.0, 0.98031767865047559,
    0.0012796818734306208, 0.0, 0.0, 0.0, 0.0, 0.97836986232776568,
    0.0014053948182610615, 0.0, 0.0, 0.0, 0.0, 0.97642573090274454,
    0.0015307016397714866, 0.0, 0.0, 0.0, 0.0, 0.97448527764739334,
    0.0016556033114226655, 0.0, 0.0, 0.0, 0.0, 0.97254849584564251,
    0.0017801008046152487, 0.0, 0.0, 0.0, 0.0, 0.97061537879335091,
    0.0019041950886938322, 0.0, 0.0, 0.0, 0.0, 0.968685919798285,
    0.0020278871309510108, 0.0, 0.0, 0.0, 0.0, 0.96676011218009861,
    0.0021511778966314264, 0.0, 0.0, 0.0, 0.0, 0.96483794927031175,
    0.0022740683489358075, 0.0, 0.0, 0.0, 0.0, 0.9629194244122905,
    0.0023965594490250011, 0.0, 0.0, 0.0, 0.0, 0.96100453096122618,
    0.0025186521560239968, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    0.001444494825713849, -0.99878101747798131, -0.0, -0.0, -0.0, -0.0,
    0.0028843614603847466, -0.99756333333822533, -0.0, -0.0, -0.0, -0.0,
    0.0043196109666541723, -0.99634694659889633, -0.0, -0.0, -0.0, -0.0,
    0.0057502543837855847, -0.99513185627791911, -0.0, -0.0, -0.0, -0.0,
    0.0071763027277104837, -0.99391806139298267, -0.0, -0.0, -0.0, -0.0,
    0.0085977669910743838, -0.99270556096154294, -0.0, -0.0, -0.0, -0.0,
    0.0100146581432827, -0.99149435400082664, -0.0, -0.0, -0.0, -0.0,
    0.011426987130546547, -0.99028443952783429, -0.0, -0.0, -0.0, -0.0,
    0.012834764875928464, -0.9890758165593434, -0.0, -0.0, -0.0, -0.0,
    0.014238002279388031, -0.987868484111912, -0.0, -0.0, -0.0, -0.0,
    0.015636710217827417, -0.98666244120188173, -0.0, -0.0, -0.0, -0.0,
    0.017030899545136844, -0.98545768684538126, -0.0, -0.0, -0.0, -0.0,
    0.018420581092239943, -0.98425422005832941, -0.0, -0.0, -0.0, -0.0,
    0.019805765667139059, -0.98305203985643841, -0.0, -0.0, -0.0, -0.0,
    0.021186464054960434, -0.98185114525521711, -0.0, -0.0, -0.0, -0.0,
    0.022562687017999343, -0.98065153526997417, -0.0, -0.0, -0.0, -0.0,
    0.023934445295765105, -0.97945320891582144, -0.0, -0.0, -0.0, -0.0,
    0.025301749605026048, -0.97825616520767711, -0.0, -0.0, -0.0, -0.0,
    0.026664610639854355, -0.97706040316026854, -0.0, -0.0, -0.0, -0.0,
    0.028023039071670849, -0.975865921788136, -0.0, -0.0, -0.0, 0.0,
    -0.001444494825713849, 0.99878101747798131, 0.0, 0.0, 0.0, 0.0,
    -0.0028843614603847466, 0.99756333333822533, 0.0, 0.0, 0.0, 0.0,
    -0.0043196109666541723, 0.99634694659889633, 0.0, 0.0, 0.0, 0.0,
    -0.0057502543837855847, 0.99513185627791911, 0.0, 0.0, 0.0, 0.0,
    -0.0071763027277104837, 0.99391806139298267, 0.0, 0.0, 0.0, 0.0,
    -0.0085977669910743838, 0.99270556096154294, 0.0, 0.0, 0.0, 0.0,
    -0.0100146581432827, 0.99149435400082664, 0.0, 0.0, 0.0, 0.0,
    -0.011426987130546547, 0.99028443952783429, 0.0, 0.0, 0.0, 0.0,
    -0.012834764875928464, 0.9890758165593434, 0.0, 0.0, 0.0, 0.0,
    -0.014238002279388031, 0.987868484111912, 0.0, 0.0, 0.0, 0.0,
    -0.015636710217827417, 0.98666244120188173, 0.0, 0.0, 0.0, 0.0,
    -0.017030899545136844, 0.98545768684538126, 0.0, 0.0, 0.0, 0.0,
    -0.018420581092239943, 0.98425422005832941, 0.0, 0.0, 0.0, 0.0,
    -0.019805765667139059, 0.98305203985643841, 0.0, 0.0, 0.0, 0.0,
    -0.021186464054960434, 0.98185114525521711, 0.0, 0.0, 0.0, 0.0,
    -0.022562687017999343, 0.98065153526997417, 0.0, 0.0, 0.0, 0.0,
    -0.023934445295765105, 0.97945320891582144, 0.0, 0.0, 0.0, 0.0,
    -0.025301749605026048, 0.97825616520767711, 0.0, 0.0, 0.0, 0.0,
    -0.026664610639854355, 0.97706040316026854, 0.0, 0.0, 0.0, 0.0,
    -0.028023039071670849, 0.975865921788136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.99801496201563555,
    -0.00012982817451900192, -0.0, -0.0, -0.0, -0.0, -0.99603367687094413,
    -0.00025924037690453866, -0.0, -0.0, -0.0, -0.0, -0.99405613771726609,
    -0.00038823760144369459, -0.0, -0.0, -0.0, -0.0, -0.99208233771810028,
    -0.00051682084032238652, -0.0, -0.0, -0.0, -0.0, -0.99011227004908309,
    -0.00064499108362950342, -0.0, -0.0, -0.0, -0.0, -0.9881459278979674,
    -0.00077274931936103868, -0.0, -0.0, -0.0, -0.0, -0.98618330446460145,
    -0.00090009653342421417, -0.0, -0.0, -0.0, -0.0, -0.984224392960908,
    -0.0010270337096415967, -0.0, -0.0, -0.0, -0.0, -0.98226918661086315,
    -0.001153561829755207, -0.0, -0.0, -0.0, -0.0, -0.98031767865047559,
    -0.0012796818734306208, -0.0, -0.0, -0.0, -0.0, -0.97836986232776568,
    -0.0014053948182610615, -0.0, -0.0, -0.0, -0.0, -0.97642573090274454,
    -0.0015307016397714866, -0.0, -0.0, -0.0, -0.0, -0.97448527764739334,
    -0.0016556033114226655, -0.0, -0.0, -0.0, -0.0, -0.97254849584564251,
    -0.0017801008046152487, -0.0, -0.0, -0.0, -0.0, -0.97061537879335091,
    -0.0019041950886938322, -0.0, -0.0, -0.0, -0.0, -0.968685919798285,
    -0.0020278871309510108, -0.0, -0.0, -0.0, -0.0, -0.96676011218009861,
    -0.0021511778966314264, -0.0, -0.0, -0.0, -0.0, -0.96483794927031175,
    -0.0022740683489358075, -0.0, -0.0, -0.0, -0.0, -0.9629194244122905,
    -0.0023965594490250011, -0.0, -0.0, -0.0, -0.0, -0.96100453096122618,
    -0.0025186521560239968, 0.0, 0.0, 0.0, 0.0, 0.99801496201563555,
    0.00012982817451900192, 0.0, 0.0, 0.0, 0.0, 0.99603367687094413,
    0.00025924037690453866, 0.0, 0.0, 0.0, 0.0, 0.99405613771726609,
    0.00038823760144369459, 0.0, 0.0, 0.0, 0.0, 0.99208233771810028,
    0.00051682084032238652, 0.0, 0.0, 0.0, 0.0, 0.99011227004908309,
    0.00064499108362950342, 0.0, 0.0, 0.0, 0.0, 0.9881459278979674,
    0.00077274931936103868, 0.0, 0.0, 0.0, 0.0, 0.98618330446460145,
    0.00090009653342421417, 0.0, 0.0, 0.0, 0.0, 0.984224392960908,
    0.0010270337096415967, 0.0, 0.0, 0.0, 0.0, 0.98226918661086315,
    0.001153561829755207, 0.0, 0.0, 0.0, 0.0, 0.98031767865047559,
    0.0012796818734306208, 0.0, 0.0, 0.0, 0.0, 0.97836986232776568,
    0.0014053948182610615, 0.0, 0.0, 0.0, 0.0, 0.97642573090274454,
    0.0015307016397714866, 0.0, 0.0, 0.0, 0.0, 0.97448527764739334,
    0.0016556033114226655, 0.0, 0.0, 0.0, 0.0, 0.97254849584564251,
    0.0017801008046152487, 0.0, 0.0, 0.0, 0.0, 0.97061537879335091,
    0.0019041950886938322, 0.0, 0.0, 0.0, 0.0, 0.968685919798285,
    0.0020278871309510108, 0.0, 0.0, 0.0, 0.0, 0.96676011218009861,
    0.0021511778966314264, 0.0, 0.0, 0.0, 0.0, 0.96483794927031175,
    0.0022740683489358075, 0.0, 0.0, 0.0, 0.0, 0.9629194244122905,
    0.0023965594490250011, 0.0, 0.0, 0.0, 0.0, 0.96100453096122618,
    0.0025186521560239968, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0,
    0.001444494825713849, -0.99878101747798131, -0.0, -0.0, -0.0, -0.0,
    0.0028843614603847466, -0.99756333333822533, -0.0, -0.0, -0.0, -0.0,
    0.0043196109666541723, -0.99634694659889633, -0.0, -0.0, -0.0, -0.0,
    0.0057502543837855847, -0.99513185627791911, -0.0, -0.0, -0.0, -0.0,
    0.0071763027277104837, -0.99391806139298267, -0.0, -0.0, -0.0, -0.0,
    0.0085977669910743838, -0.99270556096154294, -0.0, -0.0, -0.0, -0.0,
    0.0100146581432827, -0.99149435400082664, -0.0, -0.0, -0.0, -0.0,
    0.011426987130546547, -0.99028443952783429, -0.0, -0.0, -0.0, -0.0,
    0.012834764875928464, -0.9890758165593434, -0.0, -0.0, -0.0, -0.0,
    0.014238002279388031, -0.987868484111912, -0.0, -0.0, -0.0, -0.0,
    0.015636710217827417, -0.98666244120188173, -0.0, -0.0, -0.0, -0.0,
    0.017030899545136844, -0.98545768684538126, -0.0, -0.0, -0.0, -0.0,
    0.018420581092239943, -0.98425422005832941, -0.0, -0.0, -0.0, -0.0,
    0.019805765667139059, -0.98305203985643841, -0.0, -0.0, -0.0, -0.0,
    0.021186464054960434, -0.98185114525521711, -0.0, -0.0, -0.0, -0.0,
    0.022562687017999343, -0.98065153526997417, -0.0, -0.0, -0.0, -0.0,
    0.023934445295765105, -0.97945320891582144, -0.0, -0.0, -0.0, -0.0,
    0.025301749605026048, -0.97825616520767711, -0.0, -0.0, -0.0, -0.0,
    0.026664610639854355, -0.97706040316026854, -0.0, -0.0, -0.0, -0.0,
    0.028023039071670849, -0.975865921788136, 0.0, 0.0, 0.0, 0.0,
    -0.001444494825713849, 0.99878101747798131, 0.0, 0.0, 0.0, 0.0,
    -0.0028843614603847466, 0.99756333333822533, 0.0, 0.0, 0.0, 0.0,
    -0.0043196109666541723, 0.99634694659889633, 0.0, 0.0, 0.0, 0.0,
    -0.0057502543837855847, 0.99513185627791911, 0.0, 0.0, 0.0, 0.0,
    -0.0071763027277104837, 0.99391806139298267, 0.0, 0.0, 0.0, 0.0,
    -0.0085977669910743838, 0.99270556096154294, 0.0, 0.0, 0.0, 0.0,
    -0.0100146581432827, 0.99149435400082664, 0.0, 0.0, 0.0, 0.0,
    -0.011426987130546547, 0.99028443952783429, 0.0, 0.0, 0.0, 0.0,
    -0.012834764875928464, 0.9890758165593434, 0.0, 0.0, 0.0, 0.0,
    -0.014238002279388031, 0.987868484111912, 0.0, 0.0, 0.0, 0.0,
    -0.015636710217827417, 0.98666244120188173, 0.0, 0.0, 0.0, 0.0,
    -0.017030899545136844, 0.98545768684538126, 0.0, 0.0, 0.0, 0.0,
    -0.018420581092239943, 0.98425422005832941, 0.0, 0.0, 0.0, 0.0,
    -0.019805765667139059, 0.98305203985643841, 0.0, 0.0, 0.0, 0.0,
    -0.021186464054960434, 0.98185114525521711, 0.0, 0.0, 0.0, 0.0,
    -0.022562687017999343, 0.98065153526997417, 0.0, 0.0, 0.0, 0.0,
    -0.023934445295765105, 0.97945320891582144, 0.0, 0.0, 0.0, 0.0,
    -0.025301749605026048, 0.97825616520767711, 0.0, 0.0, 0.0, 0.0,
    -0.026664610639854355, 0.97706040316026854, 0.0, 0.0, 0.0, 0.0,
    -0.028023039071670849, 0.975865921788136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T d_1[2460]{ -0.99747143658641013, -0.0011849245467832059,
    -0.0, -0.0, -0.0, -0.0, -0.9949489107331112, -0.0023674262193844361, -0.0,
    -0.0, -0.0, -0.0, -0.99243240790182863, -0.0035475109996374703, -0.0, -0.0,
    -0.0, -0.0, -0.9899219135892513, -0.0047251848550434227, -0.0, -0.0, -0.0,
    -0.0, -0.98741741332694732, -0.0059004537388052374, -0.0, -0.0, -0.0, -0.0,
    -0.98491889268128, -0.0070733235898620978, -0.0, -0.0, -0.0, -0.0,
    -0.98242633725332451, -0.0082438003329237571, -0.0, -0.0, -0.0, -0.0,
    -0.97993973267878431, -0.009411889878504786, -0.0, -0.0, -0.0, -0.0,
    -0.97745906462790788, -0.010577598122958734, -0.0, -0.0, -0.0, -0.0,
    -0.97498431880540515, -0.011740930948512213, -0.0, -0.0, -0.0, -0.0,
    -0.97251548095036522, -0.012901894223298897, -0.0, -0.0, -0.0, -0.0,
    -0.97005253683617321, -0.014060493801393438, -0.0, -0.0, -0.0, -0.0,
    -0.9675954722704283, -0.015216735522845306, -0.0, -0.0, -0.0, -0.0,
    -0.96514427309486084, -0.01637062521371254, -0.0, -0.0, -0.0, -0.0,
    -0.96269892518525058, -0.017522168686095428, -0.0, -0.0, -0.0, -0.0,
    -0.96025941445134477, -0.018671371738170094, -0.0, -0.0, -0.0, -0.0,
    -0.95782572683677636, -0.019818240154222011, -0.0, -0.0, -0.0, -0.0,
    -0.95539784831898256, -0.020962779704679434, -0.0, -0.0, -0.0, -0.0,
    -0.95297576490912339, -0.022104996146146753, -0.0, -0.0, -0.0, -0.0,
    -0.95055946265200075, -0.023244895221437765, -0.0, -0.0, -0.0, -0.0,
    0.99747143658641013, 0.0011849245467832059, 0.0, 0.0, 0.0, 0.0,
    0.9949489107331112, 0.0023674262193844361, 0.0, 0.0, 0.0, 0.0,
    0.99243240790182863, 0.0035475109996374703, 0.0, 0.0, 0.0, 0.0,
    0.9899219135892513, 0.0047251848550434227, 0.0, 0.0, 0.0, 0.0,
    0.98741741332694732, 0.0059004537388052374, 0.0, 0.0, 0.0, 0.0,
    0.98491889268128, 0.0070733235898620978, 0.0, 0.0, 0.0, 0.0,
    0.98242633725332451, 0.0082438003329237571, 0.0, 0.0, 0.0, 0.0,
    0.97993973267878431, 0.009411889878504786, 0.0, 0.0, 0.0, 0.0,
    0.97745906462790788, 0.010577598122958734, 0.0, 0.0, 0.0, 0.0,
    0.97498431880540515, 0.011740930948512213, 0.0, 0.0, 0.0, 0.0,
    0.97251548095036522, 0.012901894223298897, 0.0, 0.0, 0.0, 0.0,
    0.97005253683617321, 0.014060493801393438, 0.0, 0.0, 0.0, 0.0,
    0.9675954722704283, 0.015216735522845306, 0.0, 0.0, 0.0, 0.0,
    0.96514427309486084, 0.01637062521371254, 0.0, 0.0, 0.0, 0.0,
    0.96269892518525058, 0.017522168686095428, 0.0, 0.0, 0.0, 0.0,
    0.96025941445134477, 0.018671371738170094, 0.0, 0.0, 0.0, 0.0,
    0.95782572683677636, 0.019818240154222011, 0.0, 0.0, 0.0, 0.0,
    0.95539784831898256, 0.020962779704679434, 0.0, 0.0, 0.0, 0.0,
    0.95297576490912339, 0.022104996146146753, 0.0, 0.0, 0.0, 0.0,
    0.95055946265200075, 0.023244895221437765, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.00030050237926615431, -1.0004838136541614, -0.0, -0.0, -0.0,
    -0.0, 0.00060039030636459647, -1.0009675053113292, -0.0, -0.0, -0.0, -0.0,
    0.00089966529831621029, -1.001451075640559, -0.0, -0.0, -0.0, -0.0,
    0.0011983288685070489, -1.0019345253094327, -0.0, -0.0, -0.0, -0.0,
    0.0014963825266970827, -1.0024178549840617, -0.0, -0.0, -0.0, -0.0,
    0.0017938277790289267, -1.0029010653290908, -0.0, -0.0, -0.0, -0.0,
    0.0020906661280365455, -1.003384157007702, -0.0, -0.0, -0.0, -0.0,
    0.0023868990726539398, -1.0038671306816174, -0.0, -0.0, -0.0, -0.0,
    0.0026825281082238088, -1.004349987011103, -0.0, -0.0, -0.0, -0.0,
    0.0029775547265061954, -1.004832726654973, -0.0, -0.0, -0.0, -0.0,
    0.003271980415687107, -1.005315350270592, -0.0, -0.0, -0.0, -0.0,
    0.0035658066603871182, -1.0057978585138798, -0.0, -0.0, -0.0, -0.0,
    0.003859034941669952, -1.0062802520393141, -0.0, -0.0, -0.0, -0.0,
    0.00415166673705104, -1.0067625314999344, -0.0, -0.0, -0.0, -0.0,
    0.0044437035205060612, -1.0072446975473455, -0.0, -0.0, -0.0, -0.0,
    0.0047351467624794641, -1.0077267508317205, -0.0, -0.0, -0.0, -0.0,
    0.0050259979298929629, -1.008208692001805, -0.0, -0.0, -0.0, -0.0,
    0.0053162584861540171, -1.0086905217049205, -0.0, -0.0, -0.0, -0.0,
    0.0056059298911642881, -1.0091722405869672, -0.0, -0.0, -0.0, -0.0,
    0.0058950136013280795, -1.0096538492924283, -0.0, -0.0, -0.0, -0.0,
    -0.00030050237926615431, 1.0004838136541614, 0.0, 0.0, 0.0, 0.0,
    -0.00060039030636459647, 1.0009675053113292, 0.0, 0.0, 0.0, 0.0,
    -0.00089966529831621029, 1.001451075640559, 0.0, 0.0, 0.0, 0.0,
    -0.0011983288685070489, 1.0019345253094327, 0.0, 0.0, 0.0, 0.0,
    -0.0014963825266970827, 1.0024178549840617, 0.0, 0.0, 0.0, 0.0,
    -0.0017938277790289267, 1.0029010653290908, 0.0, 0.0, 0.0, 0.0,
    -0.0020906661280365455, 1.003384157007702, 0.0, 0.0, 0.0, 0.0,
    -0.0023868990726539398, 1.0038671306816174, 0.0, 0.0, 0.0, 0.0,
    -0.0026825281082238088, 1.004349987011103, 0.0, 0.0, 0.0, 0.0,
    -0.0029775547265061954, 1.004832726654973, 0.0, 0.0, 0.0, 0.0,
    -0.003271980415687107, 1.005315350270592, 0.0, 0.0, 0.0, 0.0,
    -0.0035658066603871182, 1.0057978585138798, 0.0, 0.0, 0.0, 0.0,
    -0.003859034941669952, 1.0062802520393141, 0.0, 0.0, 0.0, 0.0,
    -0.00415166673705104, 1.0067625314999344, 0.0, 0.0, 0.0, 0.0,
    -0.0044437035205060612, 1.0072446975473455, 0.0, 0.0, 0.0, 0.0,
    -0.0047351467624794641, 1.0077267508317205, 0.0, 0.0, 0.0, 0.0,
    -0.0050259979298929629, 1.008208692001805, 0.0, 0.0, 0.0, 0.0,
    -0.0053162584861540171, 1.0086905217049205, 0.0, 0.0, 0.0, 0.0,
    -0.0056059298911642881, 1.0091722405869672, 0.0, 0.0, 0.0, 0.0,
    -0.0058950136013280795, 1.0096538492924283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.99747143658641013,
    -0.0011849245467832059, -0.0, -0.0, -0.0, -0.0, -0.9949489107331112,
    -0.0023674262193844361, -0.0, -0.0, -0.0, -0.0, -0.99243240790182863,
    -0.0035475109996374703, -0.0, -0.0, -0.0, -0.0, -0.9899219135892513,
    -0.0047251848550434227, -0.0, -0.0, -0.0, -0.0, -0.98741741332694732,
    -0.0059004537388052374, -0.0, -0.0, -0.0, -0.0, -0.98491889268128,
    -0.0070733235898620978, -0.0, -0.0, -0.0, -0.0, -0.98242633725332451,
    -0.0082438003329237571, -0.0, -0.0, -0.0, -0.0, -0.97993973267878431,
    -0.009411889878504786, -0.0, -0.0, -0.0, -0.0, -0.97745906462790788,
    -0.010577598122958734, -0.0, -0.0, -0.0, -0.0, -0.97498431880540515,
    -0.011740930948512213, -0.0, -0.0, -0.0, -0.0, -0.97251548095036522,
    -0.012901894223298897, -0.0, -0.0, -0.0, -0.0, -0.97005253683617321,
    -0.014060493801393438, -0.0, -0.0, -0.0, -0.0, -0.9675954722704283,
    -0.015216735522845306, -0.0, -0.0, -0.0, -0.0, -0.96514427309486084,
    -0.01637062521371254, -0.0, -0.0, -0.0, -0.0, -0.96269892518525058,
    -0.017522168686095428, -0.0, -0.0, -0.0, -0.0, -0.96025941445134477,
    -0.018671371738170094, -0.0, -0.0, -0.0, -0.0, -0.95782572683677636,
    -0.019818240154222011, -0.0, -0.0, -0.0, -0.0, -0.95539784831898256,
    -0.020962779704679434, -0.0, -0.0, -0.0, -0.0, -0.95297576490912339,
    -0.022104996146146753, -0.0, -0.0, -0.0, -0.0, -0.95055946265200075,
    -0.023244895221437765, -0.0, 0.0, 0.0, 0.0, 0.99747143658641013,
    0.0011849245467832059, 0.0, 0.0, 0.0, 0.0, 0.9949489107331112,
    0.0023674262193844361, 0.0, 0.0, 0.0, 0.0, 0.99243240790182863,
    0.0035475109996374703, 0.0, 0.0, 0.0, 0.0, 0.9899219135892513,
    0.0047251848550434227, 0.0, 0.0, 0.0, 0.0, 0.98741741332694732,
    0.0059004537388052374, 0.0, 0.0, 0.0, 0.0, 0.98491889268128,
    0.0070733235898620978, 0.0, 0.0, 0.0, 0.0, 0.98242633725332451,
    0.0082438003329237571, 0.0, 0.0, 0.0, 0.0, 0.97993973267878431,
    0.009411889878504786, 0.0, 0.0, 0.0, 0.0, 0.97745906462790788,
    0.010577598122958734, 0.0, 0.0, 0.0, 0.0, 0.97498431880540515,
    0.011740930948512213, 0.0, 0.0, 0.0, 0.0, 0.97251548095036522,
    0.012901894223298897, 0.0, 0.0, 0.0, 0.0, 0.97005253683617321,
    0.014060493801393438, 0.0, 0.0, 0.0, 0.0, 0.9675954722704283,
    0.015216735522845306, 0.0, 0.0, 0.0, 0.0, 0.96514427309486084,
    0.01637062521371254, 0.0, 0.0, 0.0, 0.0, 0.96269892518525058,
    0.017522168686095428, 0.0, 0.0, 0.0, 0.0, 0.96025941445134477,
    0.018671371738170094, 0.0, 0.0, 0.0, 0.0, 0.95782572683677636,
    0.019818240154222011, 0.0, 0.0, 0.0, 0.0, 0.95539784831898256,
    0.020962779704679434, 0.0, 0.0, 0.0, 0.0, 0.95297576490912339,
    0.022104996146146753, 0.0, 0.0, 0.0, 0.0, 0.95055946265200075,
    0.023244895221437765, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    0.00030050237926615431, -1.0004838136541614, -0.0, -0.0, -0.0, -0.0,
    0.00060039030636459647, -1.0009675053113292, -0.0, -0.0, -0.0, -0.0,
    0.00089966529831621029, -1.001451075640559, -0.0, -0.0, -0.0, -0.0,
    0.0011983288685070489, -1.0019345253094327, -0.0, -0.0, -0.0, -0.0,
    0.0014963825266970827, -1.0024178549840617, -0.0, -0.0, -0.0, -0.0,
    0.0017938277790289267, -1.0029010653290908, -0.0, -0.0, -0.0, -0.0,
    0.0020906661280365455, -1.003384157007702, -0.0, -0.0, -0.0, -0.0,
    0.0023868990726539398, -1.0038671306816174, -0.0, -0.0, -0.0, -0.0,
    0.0026825281082238088, -1.004349987011103, -0.0, -0.0, -0.0, -0.0,
    0.0029775547265061954, -1.004832726654973, -0.0, -0.0, -0.0, -0.0,
    0.003271980415687107, -1.005315350270592, -0.0, -0.0, -0.0, -0.0,
    0.0035658066603871182, -1.0057978585138798, -0.0, -0.0, -0.0, -0.0,
    0.003859034941669952, -1.0062802520393141, -0.0, -0.0, -0.0, -0.0,
    0.00415166673705104, -1.0067625314999344, -0.0, -0.0, -0.0, -0.0,
    0.0044437035205060612, -1.0072446975473455, -0.0, -0.0, -0.0, -0.0,
    0.0047351467624794641, -1.0077267508317205, -0.0, -0.0, -0.0, -0.0,
    0.0050259979298929629, -1.008208692001805, -0.0, -0.0, -0.0, -0.0,
    0.0053162584861540171, -1.0086905217049205, -0.0, -0.0, -0.0, -0.0,
    0.0056059298911642881, -1.0091722405869672, -0.0, -0.0, -0.0, -0.0,
    0.0058950136013280795, -1.0096538492924283, -0.0, 0.0, 0.0, 0.0,
    -0.00030050237926615431, 1.0004838136541614, 0.0, 0.0, 0.0, 0.0,
    -0.00060039030636459647, 1.0009675053113292, 0.0, 0.0, 0.0, 0.0,
    -0.00089966529831621029, 1.001451075640559, 0.0, 0.0, 0.0, 0.0,
    -0.0011983288685070489, 1.0019345253094327, 0.0, 0.0, 0.0, 0.0,
    -0.0014963825266970827, 1.0024178549840617, 0.0, 0.0, 0.0, 0.0,
    -0.0017938277790289267, 1.0029010653290908, 0.0, 0.0, 0.0, 0.0,
    -0.0020906661280365455, 1.003384157007702, 0.0, 0.0, 0.0, 0.0,
    -0.0023868990726539398, 1.0038671306816174, 0.0, 0.0, 0.0, 0.0,
    -0.0026825281082238088, 1.004349987011103, 0.0, 0.0, 0.0, 0.0,
    -0.0029775547265061954, 1.004832726654973, 0.0, 0.0, 0.0, 0.0,
    -0.003271980415687107, 1.005315350270592, 0.0, 0.0, 0.0, 0.0,
    -0.0035658066603871182, 1.0057978585138798, 0.0, 0.0, 0.0, 0.0,
    -0.003859034941669952, 1.0062802520393141, 0.0, 0.0, 0.0, 0.0,
    -0.00415166673705104, 1.0067625314999344, 0.0, 0.0, 0.0, 0.0,
    -0.0044437035205060612, 1.0072446975473455, 0.0, 0.0, 0.0, 0.0,
    -0.0047351467624794641, 1.0077267508317205, 0.0, 0.0, 0.0, 0.0,
    -0.0050259979298929629, 1.008208692001805, 0.0, 0.0, 0.0, 0.0,
    -0.0053162584861540171, 1.0086905217049205, 0.0, 0.0, 0.0, 0.0,
    -0.0056059298911642881, 1.0091722405869672, 0.0, 0.0, 0.0, 0.0,
    -0.0058950136013280795, 1.0096538492924283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T d[1968]{ -0.99867923888102927, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.99736022217199194, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.99604294756893919, -0.0, -0.0, -0.0, -0.0, -0.0, -0.9947274127709651,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.99341361548020291, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.99210155340182049, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.99079122424401689, -0.0, -0.0, -0.0, -0.0, -0.0, -0.989482625718018, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.98817575553807258, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.98687061142144838, -0.0, -0.0, -0.0, -0.0, -0.0, -0.985567191088428, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.98426549226230531, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.9829655126693807, -0.0, -0.0, -0.0, -0.0, -0.0, -0.98166725003895783,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.98037070210333943, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.97907586659782331, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.97778274126069831, -0.0, -0.0, -0.0, -0.0, -0.0, -0.97649132383324055,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.97520161205970934, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.97391360368734325, -0.0, -0.0, -0.0, -0.0, -0.0,
    0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99736022217199194, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99604294756893919, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.9947274127709651, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99341361548020291, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99210155340182049, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99079122424401689, 0.0, 0.0, 0.0, 0.0, 0.0, 0.989482625718018, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98817575553807258, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.98687061142144838, 0.0, 0.0, 0.0, 0.0, 0.0, 0.985567191088428, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98426549226230531, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.9829655126693807, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98166725003895783, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98037070210333943, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97907586659782331, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97778274126069831, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97649132383324055, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97520161205970934, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97391360368734325, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -0.99867923888102927, -0.0, -0.0, -0.0, -0.0, -0.0, -0.99736022217199194,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.99604294756893919, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.9947274127709651, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.99341361548020291, -0.0, -0.0, -0.0, -0.0, -0.0, -0.99210155340182049,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.99079122424401689, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.989482625718018, -0.0, -0.0, -0.0, -0.0, -0.0, -0.98817575553807258,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.98687061142144838, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.985567191088428, -0.0, -0.0, -0.0, -0.0, -0.0, -0.98426549226230531,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.9829655126693807, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.98166725003895783, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.98037070210333943, -0.0, -0.0, -0.0, -0.0, -0.0, -0.97907586659782331,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.97778274126069831, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.97649132383324055, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.97520161205970934, -0.0, -0.0, -0.0, -0.0, -0.0, -0.97391360368734325,
    -0.0, -0.0, 0.0, 0.0, 0.0, 0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99736022217199194, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99604294756893919, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.9947274127709651, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99341361548020291, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99210155340182049, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99079122424401689, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.989482625718018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98817575553807258, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98687061142144838, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.985567191088428, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98426549226230531, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.9829655126693807, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.98166725003895783, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98037070210333943, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97907586659782331, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97778274126069831, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97649132383324055, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97520161205970934, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97391360368734325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0,
    -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0,
    -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0,
    -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0,
    -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0,
    -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0,
    -0.0, -0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T n_0[1200]{ 0.0, 0.99801496201563555,
    0.00012982817451900192, 0.0, 0.0, 0.0, 0.0, 0.99603367687094413,
    0.00025924037690453866, 0.0, 0.0, 0.0, 0.0, 0.99405613771726609,
    0.00038823760144369459, 0.0, 0.0, 0.0, 0.0, 0.99208233771810028,
    0.00051682084032238652, 0.0, 0.0, 0.0, 0.0, 0.99011227004908309,
    0.00064499108362950342, 0.0, 0.0, 0.0, 0.0, 0.9881459278979674,
    0.00077274931936103868, 0.0, 0.0, 0.0, 0.0, 0.98618330446460145,
    0.00090009653342421417, 0.0, 0.0, 0.0, 0.0, 0.984224392960908,
    0.0010270337096415967, 0.0, 0.0, 0.0, 0.0, 0.98226918661086315,
    0.001153561829755207, 0.0, 0.0, 0.0, 0.0, 0.98031767865047559,
    0.0012796818734306208, 0.0, 0.0, 0.0, 0.0, 0.97836986232776568,
    0.0014053948182610615, 0.0, 0.0, 0.0, 0.0, 0.97642573090274454,
    0.0015307016397714866, 0.0, 0.0, 0.0, 0.0, 0.97448527764739334,
    0.0016556033114226655, 0.0, 0.0, 0.0, 0.0, 0.97254849584564251,
    0.0017801008046152487, 0.0, 0.0, 0.0, 0.0, 0.97061537879335091,
    0.0019041950886938322, 0.0, 0.0, 0.0, 0.0, 0.968685919798285,
    0.0020278871309510108, 0.0, 0.0, 0.0, 0.0, 0.96676011218009861,
    0.0021511778966314264, 0.0, 0.0, 0.0, 0.0, 0.96483794927031175,
    0.0022740683489358075, 0.0, 0.0, 0.0, 0.0, 0.9629194244122905,
    0.0023965594490250011, 0.0, 0.0, 0.0, 0.0, 0.96100453096122618,
    0.0025186521560239968, 0.0, 0.0, 0.0, 0.0, -0.001444494825713849,
    0.99878101747798131, 0.0, 0.0, 0.0, 0.0, -0.0028843614603847466,
    0.99756333333822533, 0.0, 0.0, 0.0, 0.0, -0.0043196109666541723,
    0.99634694659889633, 0.0, 0.0, 0.0, 0.0, -0.0057502543837855847,
    0.99513185627791911, 0.0, 0.0, 0.0, 0.0, -0.0071763027277104837,
    0.99391806139298267, 0.0, 0.0, 0.0, 0.0, -0.0085977669910743838,
    0.99270556096154294, 0.0, 0.0, 0.0, 0.0, -0.0100146581432827,
    0.99149435400082664, 0.0, 0.0, 0.0, 0.0, -0.011426987130546547,
    0.99028443952783429, 0.0, 0.0, 0.0, 0.0, -0.012834764875928464,
    0.9890758165593434, 0.0, 0.0, 0.0, 0.0, -0.014238002279388031,
    0.987868484111912, 0.0, 0.0, 0.0, 0.0, -0.015636710217827417,
    0.98666244120188173, 0.0, 0.0, 0.0, 0.0, -0.017030899545136844,
    0.98545768684538126, 0.0, 0.0, 0.0, 0.0, -0.018420581092239943,
    0.98425422005832941, 0.0, 0.0, 0.0, 0.0, -0.019805765667139059,
    0.98305203985643841, 0.0, 0.0, 0.0, 0.0, -0.021186464054960434,
    0.98185114525521711, 0.0, 0.0, 0.0, 0.0, -0.022562687017999343,
    0.98065153526997417, 0.0, 0.0, 0.0, 0.0, -0.023934445295765105,
    0.97945320891582144, 0.0, 0.0, 0.0, 0.0, -0.025301749605026048,
    0.97825616520767711, 0.0, 0.0, 0.0, 0.0, -0.026664610639854355,
    0.97706040316026854, 0.0, 0.0, 0.0, 0.0, -0.028023039071670849,
    0.975865921788136, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99801496201563555,
    0.00012982817451900192, 0.0, 0.0, 0.0, 0.0, 0.99603367687094413,
    0.00025924037690453866, 0.0, 0.0, 0.0, 0.0, 0.99405613771726609,
    0.00038823760144369459, 0.0, 0.0, 0.0, 0.0, 0.99208233771810028,
    0.00051682084032238652, 0.0, 0.0, 0.0, 0.0, 0.99011227004908309,
    0.00064499108362950342, 0.0, 0.0, 0.0, 0.0, 0.9881459278979674,
    0.00077274931936103868, 0.0, 0.0, 0.0, 0.0, 0.98618330446460145,
    0.00090009653342421417, 0.0, 0.0, 0.0, 0.0, 0.984224392960908,
    0.0010270337096415967, 0.0, 0.0, 0.0, 0.0, 0.98226918661086315,
    0.001153561829755207, 0.0, 0.0, 0.0, 0.0, 0.98031767865047559,
    0.0012796818734306208, 0.0, 0.0, 0.0, 0.0, 0.97836986232776568,
    0.0014053948182610615, 0.0, 0.0, 0.0, 0.0, 0.97642573090274454,
    0.0015307016397714866, 0.0, 0.0, 0.0, 0.0, 0.97448527764739334,
    0.0016556033114226655, 0.0, 0.0, 0.0, 0.0, 0.97254849584564251,
    0.0017801008046152487, 0.0, 0.0, 0.0, 0.0, 0.97061537879335091,
    0.0019041950886938322, 0.0, 0.0, 0.0, 0.0, 0.968685919798285,
    0.0020278871309510108, 0.0, 0.0, 0.0, 0.0, 0.96676011218009861,
    0.0021511778966314264, 0.0, 0.0, 0.0, 0.0, 0.96483794927031175,
    0.0022740683489358075, 0.0, 0.0, 0.0, 0.0, 0.9629194244122905,
    0.0023965594490250011, 0.0, 0.0, 0.0, 0.0, 0.96100453096122618,
    0.0025186521560239968, 0.0, 0.0, 0.0, 0.0, -0.001444494825713849,
    0.99878101747798131, 0.0, 0.0, 0.0, 0.0, -0.0028843614603847466,
    0.99756333333822533, 0.0, 0.0, 0.0, 0.0, -0.0043196109666541723,
    0.99634694659889633, 0.0, 0.0, 0.0, 0.0, -0.0057502543837855847,
    0.99513185627791911, 0.0, 0.0, 0.0, 0.0, -0.0071763027277104837,
    0.99391806139298267, 0.0, 0.0, 0.0, 0.0, -0.0085977669910743838,
    0.99270556096154294, 0.0, 0.0, 0.0, 0.0, -0.0100146581432827,
    0.99149435400082664, 0.0, 0.0, 0.0, 0.0, -0.011426987130546547,
    0.99028443952783429, 0.0, 0.0, 0.0, 0.0, -0.012834764875928464,
    0.9890758165593434, 0.0, 0.0, 0.0, 0.0, -0.014238002279388031,
    0.987868484111912, 0.0, 0.0, 0.0, 0.0, -0.015636710217827417,
    0.98666244120188173, 0.0, 0.0, 0.0, 0.0, -0.017030899545136844,
    0.98545768684538126, 0.0, 0.0, 0.0, 0.0, -0.018420581092239943,
    0.98425422005832941, 0.0, 0.0, 0.0, 0.0, -0.019805765667139059,
    0.98305203985643841, 0.0, 0.0, 0.0, 0.0, -0.021186464054960434,
    0.98185114525521711, 0.0, 0.0, 0.0, 0.0, -0.022562687017999343,
    0.98065153526997417, 0.0, 0.0, 0.0, 0.0, -0.023934445295765105,
    0.97945320891582144, 0.0, 0.0, 0.0, 0.0, -0.025301749605026048,
    0.97825616520767711, 0.0, 0.0, 0.0, 0.0, -0.026664610639854355,
    0.97706040316026854, 0.0, 0.0, 0.0, 0.0, -0.028023039071670849,
    0.975865921788136, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0 };

  static const real_T n_1[1200]{ 0.99747143658641013, 0.0011849245467832059, 0.0,
    0.0, 0.0, 0.0, 0.9949489107331112, 0.0023674262193844361, 0.0, 0.0, 0.0, 0.0,
    0.99243240790182863, 0.0035475109996374703, 0.0, 0.0, 0.0, 0.0,
    0.9899219135892513, 0.0047251848550434227, 0.0, 0.0, 0.0, 0.0,
    0.98741741332694732, 0.0059004537388052374, 0.0, 0.0, 0.0, 0.0,
    0.98491889268128, 0.0070733235898620978, 0.0, 0.0, 0.0, 0.0,
    0.98242633725332451, 0.0082438003329237571, 0.0, 0.0, 0.0, 0.0,
    0.97993973267878431, 0.009411889878504786, 0.0, 0.0, 0.0, 0.0,
    0.97745906462790788, 0.010577598122958734, 0.0, 0.0, 0.0, 0.0,
    0.97498431880540515, 0.011740930948512213, 0.0, 0.0, 0.0, 0.0,
    0.97251548095036522, 0.012901894223298897, 0.0, 0.0, 0.0, 0.0,
    0.97005253683617321, 0.014060493801393438, 0.0, 0.0, 0.0, 0.0,
    0.9675954722704283, 0.015216735522845306, 0.0, 0.0, 0.0, 0.0,
    0.96514427309486084, 0.01637062521371254, 0.0, 0.0, 0.0, 0.0,
    0.96269892518525058, 0.017522168686095428, 0.0, 0.0, 0.0, 0.0,
    0.96025941445134477, 0.018671371738170094, 0.0, 0.0, 0.0, 0.0,
    0.95782572683677636, 0.019818240154222011, 0.0, 0.0, 0.0, 0.0,
    0.95539784831898256, 0.020962779704679434, 0.0, 0.0, 0.0, 0.0,
    0.95297576490912339, 0.022104996146146753, 0.0, 0.0, 0.0, 0.0,
    0.95055946265200075, 0.023244895221437765, 0.0, 0.0, 0.0, 0.0,
    -0.00030050237926615431, 1.0004838136541614, 0.0, 0.0, 0.0, 0.0,
    -0.00060039030636459647, 1.0009675053113292, 0.0, 0.0, 0.0, 0.0,
    -0.00089966529831621029, 1.001451075640559, 0.0, 0.0, 0.0, 0.0,
    -0.0011983288685070489, 1.0019345253094327, 0.0, 0.0, 0.0, 0.0,
    -0.0014963825266970827, 1.0024178549840617, 0.0, 0.0, 0.0, 0.0,
    -0.0017938277790289267, 1.0029010653290908, 0.0, 0.0, 0.0, 0.0,
    -0.0020906661280365455, 1.003384157007702, 0.0, 0.0, 0.0, 0.0,
    -0.0023868990726539398, 1.0038671306816174, 0.0, 0.0, 0.0, 0.0,
    -0.0026825281082238088, 1.004349987011103, 0.0, 0.0, 0.0, 0.0,
    -0.0029775547265061954, 1.004832726654973, 0.0, 0.0, 0.0, 0.0,
    -0.003271980415687107, 1.005315350270592, 0.0, 0.0, 0.0, 0.0,
    -0.0035658066603871182, 1.0057978585138798, 0.0, 0.0, 0.0, 0.0,
    -0.003859034941669952, 1.0062802520393141, 0.0, 0.0, 0.0, 0.0,
    -0.00415166673705104, 1.0067625314999344, 0.0, 0.0, 0.0, 0.0,
    -0.0044437035205060612, 1.0072446975473455, 0.0, 0.0, 0.0, 0.0,
    -0.0047351467624794641, 1.0077267508317205, 0.0, 0.0, 0.0, 0.0,
    -0.0050259979298929629, 1.008208692001805, 0.0, 0.0, 0.0, 0.0,
    -0.0053162584861540171, 1.0086905217049205, 0.0, 0.0, 0.0, 0.0,
    -0.0056059298911642881, 1.0091722405869672, 0.0, 0.0, 0.0, 0.0,
    -0.0058950136013280795, 1.0096538492924283, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.99747143658641013, 0.0011849245467832059, 0.0, 0.0, 0.0, 0.0,
    0.9949489107331112, 0.0023674262193844361, 0.0, 0.0, 0.0, 0.0,
    0.99243240790182863, 0.0035475109996374703, 0.0, 0.0, 0.0, 0.0,
    0.9899219135892513, 0.0047251848550434227, 0.0, 0.0, 0.0, 0.0,
    0.98741741332694732, 0.0059004537388052374, 0.0, 0.0, 0.0, 0.0,
    0.98491889268128, 0.0070733235898620978, 0.0, 0.0, 0.0, 0.0,
    0.98242633725332451, 0.0082438003329237571, 0.0, 0.0, 0.0, 0.0,
    0.97993973267878431, 0.009411889878504786, 0.0, 0.0, 0.0, 0.0,
    0.97745906462790788, 0.010577598122958734, 0.0, 0.0, 0.0, 0.0,
    0.97498431880540515, 0.011740930948512213, 0.0, 0.0, 0.0, 0.0,
    0.97251548095036522, 0.012901894223298897, 0.0, 0.0, 0.0, 0.0,
    0.97005253683617321, 0.014060493801393438, 0.0, 0.0, 0.0, 0.0,
    0.9675954722704283, 0.015216735522845306, 0.0, 0.0, 0.0, 0.0,
    0.96514427309486084, 0.01637062521371254, 0.0, 0.0, 0.0, 0.0,
    0.96269892518525058, 0.017522168686095428, 0.0, 0.0, 0.0, 0.0,
    0.96025941445134477, 0.018671371738170094, 0.0, 0.0, 0.0, 0.0,
    0.95782572683677636, 0.019818240154222011, 0.0, 0.0, 0.0, 0.0,
    0.95539784831898256, 0.020962779704679434, 0.0, 0.0, 0.0, 0.0,
    0.95297576490912339, 0.022104996146146753, 0.0, 0.0, 0.0, 0.0,
    0.95055946265200075, 0.023244895221437765, 0.0, 0.0, 0.0, 0.0,
    -0.00030050237926615431, 1.0004838136541614, 0.0, 0.0, 0.0, 0.0,
    -0.00060039030636459647, 1.0009675053113292, 0.0, 0.0, 0.0, 0.0,
    -0.00089966529831621029, 1.001451075640559, 0.0, 0.0, 0.0, 0.0,
    -0.0011983288685070489, 1.0019345253094327, 0.0, 0.0, 0.0, 0.0,
    -0.0014963825266970827, 1.0024178549840617, 0.0, 0.0, 0.0, 0.0,
    -0.0017938277790289267, 1.0029010653290908, 0.0, 0.0, 0.0, 0.0,
    -0.0020906661280365455, 1.003384157007702, 0.0, 0.0, 0.0, 0.0,
    -0.0023868990726539398, 1.0038671306816174, 0.0, 0.0, 0.0, 0.0,
    -0.0026825281082238088, 1.004349987011103, 0.0, 0.0, 0.0, 0.0,
    -0.0029775547265061954, 1.004832726654973, 0.0, 0.0, 0.0, 0.0,
    -0.003271980415687107, 1.005315350270592, 0.0, 0.0, 0.0, 0.0,
    -0.0035658066603871182, 1.0057978585138798, 0.0, 0.0, 0.0, 0.0,
    -0.003859034941669952, 1.0062802520393141, 0.0, 0.0, 0.0, 0.0,
    -0.00415166673705104, 1.0067625314999344, 0.0, 0.0, 0.0, 0.0,
    -0.0044437035205060612, 1.0072446975473455, 0.0, 0.0, 0.0, 0.0,
    -0.0047351467624794641, 1.0077267508317205, 0.0, 0.0, 0.0, 0.0,
    -0.0050259979298929629, 1.008208692001805, 0.0, 0.0, 0.0, 0.0,
    -0.0053162584861540171, 1.0086905217049205, 0.0, 0.0, 0.0, 0.0,
    -0.0056059298911642881, 1.0091722405869672, 0.0, 0.0, 0.0, 0.0,
    -0.0058950136013280795, 1.0096538492924283, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T g[984]{ -0.3579347593572087, -0.0, -0.0,
    -0.3579347593572087, -0.0, -0.0, -0.71539677240113031, -0.0, -0.0,
    -0.71539677240113031, -0.0, -0.0, -1.0723866635167145, -0.0, -0.0,
    -1.0723866635167145, -0.0, -0.0, -1.4289050562642476, -0.0, -0.0,
    -1.4289050562642476, -0.0, -0.0, -1.7849525733804419, -0.0, -0.0,
    -1.7849525733804419, -0.0, -0.0, -2.1405298367795229, -0.0, -0.0,
    -2.1405298367795229, -0.0, -0.0, -2.4956374675543165, -0.0, -0.0,
    -2.4956374675543165, -0.0, -0.0, -2.8502760859773328, -0.0, -0.0,
    -2.8502760859773328, -0.0, -0.0, -3.2044463115018504, -0.0, -0.0,
    -3.2044463115018504, -0.0, -0.0, -3.5581487627629982, -0.0, -0.0,
    -3.5581487627629982, -0.0, -0.0, -3.9113840575788359, -0.0, -0.0,
    -3.9113840575788359, -0.0, -0.0, -4.2641528129514326, -0.0, -0.0,
    -4.2641528129514326, -0.0, -0.0, -4.616455645067945, -0.0, -0.0,
    -4.616455645067945, -0.0, -0.0, -4.9682931693016954, -0.0, -0.0,
    -4.9682931693016954, -0.0, -0.0, -5.3196660002132425, -0.0, -0.0,
    -5.3196660002132425, -0.0, -0.0, -5.6705747515514595, -0.0, -0.0,
    -5.6705747515514595, -0.0, -0.0, -6.021020036254602, -0.0, -0.0,
    -6.021020036254602, -0.0, -0.0, -6.3710024664513822, -0.0, -0.0,
    -6.3710024664513822, -0.0, -0.0, -6.7205226534620355, -0.0, -0.0,
    -6.7205226534620355, -0.0, -0.0, -7.06958120779939, -0.0, -0.0,
    -7.06958120779939, -0.0, -0.0, 0.3579347593572087, 0.0, 0.0,
    0.3579347593572087, 0.0, 0.0, 0.71539677240113031, 0.0, 0.0,
    0.71539677240113031, 0.0, 0.0, 1.0723866635167145, 0.0, 0.0,
    1.0723866635167145, 0.0, 0.0, 1.4289050562642476, 0.0, 0.0,
    1.4289050562642476, 0.0, 0.0, 1.7849525733804419, 0.0, 0.0,
    1.7849525733804419, 0.0, 0.0, 2.1405298367795229, 0.0, 0.0,
    2.1405298367795229, 0.0, 0.0, 2.4956374675543165, 0.0, 0.0,
    2.4956374675543165, 0.0, 0.0, 2.8502760859773328, 0.0, 0.0,
    2.8502760859773328, 0.0, 0.0, 3.2044463115018504, 0.0, 0.0,
    3.2044463115018504, 0.0, 0.0, 3.5581487627629982, 0.0, 0.0,
    3.5581487627629982, 0.0, 0.0, 3.9113840575788359, 0.0, 0.0,
    3.9113840575788359, 0.0, 0.0, 4.2641528129514326, 0.0, 0.0,
    4.2641528129514326, 0.0, 0.0, 4.616455645067945, 0.0, 0.0, 4.616455645067945,
    0.0, 0.0, 4.9682931693016954, 0.0, 0.0, 4.9682931693016954, 0.0, 0.0,
    5.3196660002132425, 0.0, 0.0, 5.3196660002132425, 0.0, 0.0,
    5.6705747515514595, 0.0, 0.0, 5.6705747515514595, 0.0, 0.0,
    6.021020036254602, 0.0, 0.0, 6.021020036254602, 0.0, 0.0, 6.3710024664513822,
    0.0, 0.0, 6.3710024664513822, 0.0, 0.0, 6.7205226534620355, 0.0, 0.0,
    6.7205226534620355, 0.0, 0.0, 7.06958120779939, 0.0, 0.0, 7.06958120779939,
    0.0, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.042560663779893933, -0.0, -0.0,
    0.042560663779893933, -0.0, -0.0, 0.0850651150898698, -0.0, -0.0,
    0.0850651150898698, -0.0, -0.0, 0.12751342817317227, -0.0, -0.0,
    0.12751342817317227, -0.0, -0.0, 0.16990567717498842, -0.0, -0.0,
    0.16990567717498842, -0.0, -0.0, 0.21224193614257725, -0.0, -0.0,
    0.21224193614257725, -0.0, -0.0, 0.254522279025399, -0.0, -0.0,
    0.254522279025399, -0.0, -0.0, 0.29674677967524432, -0.0, -0.0,
    0.29674677967524432, -0.0, -0.0, 0.33891551184636343, -0.0, -0.0,
    0.33891551184636343, -0.0, -0.0, 0.38102854919559465, -0.0, -0.0,
    0.38102854919559465, -0.0, -0.0, 0.4230859652824932, -0.0, -0.0,
    0.4230859652824932, -0.0, -0.0, 0.46508783356945982, -0.0, -0.0,
    0.46508783356945982, -0.0, -0.0, 0.50703422742186888, -0.0, -0.0,
    0.50703422742186888, -0.0, -0.0, 0.54892522010819667, -0.0, -0.0,
    0.54892522010819667, -0.0, -0.0, 0.59076088480014921, -0.0, -0.0,
    0.59076088480014921, -0.0, -0.0, 0.6325412945727904, -0.0, -0.0,
    0.6325412945727904, -0.0, -0.0, 0.67426652240466922, -0.0, -0.0,
    0.67426652240466922, -0.0, -0.0, 0.71593664117794742, -0.0, -0.0,
    0.71593664117794742, -0.0, -0.0, 0.757551723678527, -0.0, -0.0,
    0.757551723678527, -0.0, -0.0, 0.79911184259617707, -0.0, -0.0,
    0.79911184259617707, -0.0, -0.0, 0.840617070524661, -0.0, -0.0,
    0.840617070524661, -0.0, -0.0, -0.042560663779893933, 0.0, 0.0,
    -0.042560663779893933, 0.0, 0.0, -0.0850651150898698, 0.0, 0.0,
    -0.0850651150898698, 0.0, 0.0, -0.12751342817317227, 0.0, 0.0,
    -0.12751342817317227, 0.0, 0.0, -0.16990567717498842, 0.0, 0.0,
    -0.16990567717498842, 0.0, 0.0, -0.21224193614257725, 0.0, 0.0,
    -0.21224193614257725, 0.0, 0.0, -0.254522279025399, 0.0, 0.0,
    -0.254522279025399, 0.0, 0.0, -0.29674677967524432, 0.0, 0.0,
    -0.29674677967524432, 0.0, 0.0, -0.33891551184636343, 0.0, 0.0,
    -0.33891551184636343, 0.0, 0.0, -0.38102854919559465, 0.0, 0.0,
    -0.38102854919559465, 0.0, 0.0, -0.4230859652824932, 0.0, 0.0,
    -0.4230859652824932, 0.0, 0.0, -0.46508783356945982, 0.0, 0.0,
    -0.46508783356945982, 0.0, 0.0, -0.50703422742186888, 0.0, 0.0,
    -0.50703422742186888, 0.0, 0.0, -0.54892522010819667, 0.0, 0.0,
    -0.54892522010819667, 0.0, 0.0, -0.59076088480014921, 0.0, 0.0,
    -0.59076088480014921, 0.0, 0.0, -0.6325412945727904, 0.0, 0.0,
    -0.6325412945727904, 0.0, 0.0, -0.67426652240466922, 0.0, 0.0,
    -0.67426652240466922, 0.0, 0.0, -0.71593664117794742, 0.0, 0.0,
    -0.71593664117794742, 0.0, 0.0, -0.757551723678527, 0.0, 0.0,
    -0.757551723678527, 0.0, 0.0, -0.79911184259617707, 0.0, 0.0,
    -0.79911184259617707, 0.0, 0.0, -0.840617070524661, 0.0, 0.0,
    -0.840617070524661, 0.0, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.0099627081344734868, -0.0, -0.0, 0.0099627081344734868, -0.0, -0.0,
    0.01991225791140331, -0.0, -0.0, 0.01991225791140331, -0.0, -0.0,
    0.029848666709836498, -0.0, -0.0, 0.029848666709836498, -0.0, -0.0,
    0.039771951885866519, -0.0, -0.0, 0.039771951885866519, -0.0, -0.0,
    0.049682130772663577, -0.0, -0.0, 0.049682130772663577, -0.0, -0.0,
    0.059579220680504912, -0.0, -0.0, 0.059579220680504912, -0.0, -0.0,
    0.069463238896805016, -0.0, -0.0, 0.069463238896805016, -0.0, -0.0,
    0.079334202686145827, -0.0, -0.0, 0.079334202686145827, -0.0, -0.0,
    0.089192129290306912, -0.0, -0.0, 0.089192129290306912, -0.0, -0.0,
    0.099037035928295547, -0.0, -0.0, 0.099037035928295547, -0.0, -0.0,
    0.10886893979637684, -0.0, -0.0, 0.10886893979637684, -0.0, -0.0,
    0.11868785806810371, -0.0, -0.0, 0.11868785806810371, -0.0, -0.0,
    0.12849380789434692, -0.0, -0.0, 0.12849380789434692, -0.0, -0.0,
    0.13828680640332505, -0.0, -0.0, 0.13828680640332505, -0.0, -0.0,
    0.1480668707006344, -0.0, -0.0, 0.1480668707006344, -0.0, -0.0,
    0.15783401786927884, -0.0, -0.0, 0.15783401786927884, -0.0, -0.0,
    0.16758826496969964, -0.0, -0.0, 0.16758826496969964, -0.0, -0.0,
    0.17732962903980537, -0.0, -0.0, 0.17732962903980537, -0.0, -0.0,
    0.18705812709500158, -0.0, -0.0, 0.18705812709500158, -0.0, -0.0,
    0.1967737761282205, -0.0, -0.0, 0.1967737761282205, -0.0, -0.0,
    -0.0099627081344734868, 0.0, 0.0, -0.0099627081344734868, 0.0, 0.0,
    -0.01991225791140331, 0.0, 0.0, -0.01991225791140331, 0.0, 0.0,
    -0.029848666709836498, 0.0, 0.0, -0.029848666709836498, 0.0, 0.0,
    -0.039771951885866519, 0.0, 0.0, -0.039771951885866519, 0.0, 0.0,
    -0.049682130772663577, 0.0, 0.0, -0.049682130772663577, 0.0, 0.0,
    -0.059579220680504912, 0.0, 0.0, -0.059579220680504912, 0.0, 0.0,
    -0.069463238896805016, 0.0, 0.0, -0.069463238896805016, 0.0, 0.0,
    -0.079334202686145827, 0.0, 0.0, -0.079334202686145827, 0.0, 0.0,
    -0.089192129290306912, 0.0, 0.0, -0.089192129290306912, 0.0, 0.0,
    -0.099037035928295547, 0.0, 0.0, -0.099037035928295547, 0.0, 0.0,
    -0.10886893979637684, 0.0, 0.0, -0.10886893979637684, 0.0, 0.0,
    -0.11868785806810371, 0.0, 0.0, -0.11868785806810371, 0.0, 0.0,
    -0.12849380789434692, 0.0, 0.0, -0.12849380789434692, 0.0, 0.0,
    -0.13828680640332505, 0.0, 0.0, -0.13828680640332505, 0.0, 0.0,
    -0.1480668707006344, 0.0, 0.0, -0.1480668707006344, 0.0, 0.0,
    -0.15783401786927884, 0.0, 0.0, -0.15783401786927884, 0.0, 0.0,
    -0.16758826496969964, 0.0, 0.0, -0.16758826496969964, 0.0, 0.0,
    -0.17732962903980537, 0.0, 0.0, -0.17732962903980537, 0.0, 0.0,
    -0.18705812709500158, 0.0, 0.0, -0.18705812709500158, 0.0, 0.0,
    -0.1967737761282205, 0.0, 0.0, -0.1967737761282205, 0.0, 0.0, -0.0, -0.0,
    -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
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
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T g_0[984]{ -0.0, 0.1888968412469153, 0.060492001339316988,
    -0.0, 0.1888968412469153, 0.060492001339316988, -0.0, 0.3773313347058973,
    0.1209347881183509, -0.0, 0.3773313347058973, 0.1209347881183509, -0.0,
    0.56530446924504474, 0.18132836030302127, -0.0, 0.56530446924504474,
    0.18132836030302127, -0.0, 0.75281723176956472, 0.24167271798767209, -0.0,
    0.75281723176956472, 0.24167271798767209, -0.0, 0.939870607225484,
    0.30196786139466048, -0.0, 0.939870607225484, 0.30196786139466048, -0.0,
    1.1264655786033526, 0.36221379087394617, -0.0, 1.1264655786033526,
    0.36221379087394617, -0.0, 1.3126031269419416, 0.42241050690268223, -0.0,
    1.3126031269419416, 0.42241050690268223, -0.0, 1.4982842313319336,
    0.48255801008480659, -0.0, 1.4982842313319336, 0.48255801008480659, -0.0,
    1.6835098689196064, 0.54265630115063446, -0.0, 1.6835098689196064,
    0.54265630115063446, -0.0, 1.8682810149105107, 0.6027053809564521, -0.0,
    1.8682810149105107, 0.6027053809564521, -0.0, 2.05259864257314,
    0.66270525048411089, -0.0, 2.05259864257314, 0.66270525048411089, -0.0,
    2.2364637232425952, 0.72265591084062319, -0.0, 2.2364637232425952,
    0.72265591084062319, -0.0, 2.41987722632424, 0.78255736325775849, -0.0,
    2.41987722632424, 0.78255736325775849, -0.0, 2.6028401192973529,
    0.842409609091641, -0.0, 2.6028401192973529, 0.842409609091641, -0.0,
    2.7853533677187707, 0.90221264982234783, -0.0, 2.7853533677187707,
    0.90221264982234783, -0.0, 2.9674179352265249, 0.96196648705350851, -0.0,
    2.9674179352265249, 0.96196648705350851, -0.0, 3.1490347835434722,
    1.0216711225119051, -0.0, 3.1490347835434722, 1.0216711225119051, -0.0,
    3.3302048724809192, 1.081326558047073, -0.0, 3.3302048724809192,
    1.081326558047073, -0.0, 3.5109291599422385, 1.1409327956309046, -0.0,
    3.5109291599422385, 1.1409327956309046, -0.0, 3.69120860192648,
    1.2004898373572503, -0.0, 3.69120860192648, 1.2004898373572503, 0.0,
    -0.1888968412469153, -0.060492001339316988, 0.0, -0.1888968412469153,
    -0.060492001339316988, 0.0, -0.3773313347058973, -0.1209347881183509, 0.0,
    -0.3773313347058973, -0.1209347881183509, 0.0, -0.56530446924504474,
    -0.18132836030302127, 0.0, -0.56530446924504474, -0.18132836030302127, 0.0,
    -0.75281723176956472, -0.24167271798767209, 0.0, -0.75281723176956472,
    -0.24167271798767209, 0.0, -0.939870607225484, -0.30196786139466048, 0.0,
    -0.939870607225484, -0.30196786139466048, 0.0, -1.1264655786033526,
    -0.36221379087394617, 0.0, -1.1264655786033526, -0.36221379087394617, 0.0,
    -1.3126031269419416, -0.42241050690268223, 0.0, -1.3126031269419416,
    -0.42241050690268223, 0.0, -1.4982842313319336, -0.48255801008480659, 0.0,
    -1.4982842313319336, -0.48255801008480659, 0.0, -1.6835098689196064,
    -0.54265630115063446, 0.0, -1.6835098689196064, -0.54265630115063446, 0.0,
    -1.8682810149105107, -0.6027053809564521, 0.0, -1.8682810149105107,
    -0.6027053809564521, 0.0, -2.05259864257314, -0.66270525048411089, 0.0,
    -2.05259864257314, -0.66270525048411089, 0.0, -2.2364637232425952,
    -0.72265591084062319, 0.0, -2.2364637232425952, -0.72265591084062319, 0.0,
    -2.41987722632424, -0.78255736325775849, 0.0, -2.41987722632424,
    -0.78255736325775849, 0.0, -2.6028401192973529, -0.842409609091641, 0.0,
    -2.6028401192973529, -0.842409609091641, 0.0, -2.7853533677187707,
    -0.90221264982234783, 0.0, -2.7853533677187707, -0.90221264982234783, 0.0,
    -2.9674179352265249, -0.96196648705350851, 0.0, -2.9674179352265249,
    -0.96196648705350851, 0.0, -3.1490347835434722, -1.0216711225119051, 0.0,
    -3.1490347835434722, -1.0216711225119051, 0.0, -3.3302048724809192,
    -1.081326558047073, 0.0, -3.3302048724809192, -1.081326558047073, 0.0,
    -3.5109291599422385, -1.1409327956309046, 0.0, -3.5109291599422385,
    -1.1409327956309046, 0.0, -3.69120860192648, -1.2004898373572503, 0.0,
    -3.69120860192648, -1.2004898373572503, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    -0.0, -0.1599184353223182, 0.089801925435736424, -0.0, -0.1599184353223182,
    0.089801925435736424, -0.0, -0.31964914489275253, 0.17947362197539324, -0.0,
    -0.31964914489275253, 0.17947362197539324, -0.0, -0.47919231323908462,
    0.26901527273781028, -0.0, -0.47919231323908462, 0.26901527273781028, -0.0,
    -0.63854812478731549, 0.35842706059465185, -0.0, -0.63854812478731549,
    0.35842706059465185, -0.0, -0.797716763861511, 0.44770916817072121, -0.0,
    -0.797716763861511, 0.44770916817072121, -0.0, -0.956698414683647,
    0.53686177784427469, -0.0, -0.956698414683647, 0.53686177784427469, -0.0,
    -1.1154932613734565, 0.62588507174733543, -0.0, -1.1154932613734565,
    0.62588507174733543, -0.0, -1.2741014879482764, 0.71477923176600711, -0.0,
    -1.2741014879482764, 0.71477923176600711, -0.0, -1.4325232783228956,
    0.80354443954078658, -0.0, -1.4325232783228956, 0.80354443954078658, -0.0,
    -1.5907588163094044, 0.89218087646687683, -0.0, -1.5907588163094044,
    0.89218087646687683, -0.0, -1.748808285617043, 0.98068872369449933, -0.0,
    -1.748808285617043, 0.98068872369449933, -0.0, -1.9066718698520528,
    1.069068162129206, -0.0, -1.9066718698520528, 1.069068162129206, -0.0,
    -2.0643497525175265, 1.1573193724321906, -0.0, -2.0643497525175265,
    1.1573193724321906, -0.0, -2.221842117013261, 1.2454425350206004, -0.0,
    -2.221842117013261, 1.2454425350206004, -0.0, -2.3791491466356085,
    1.3334378300678469, -0.0, -2.3791491466356085, 1.3334378300678469, -0.0,
    -2.536271024577331, 1.421305437503916, -0.0, -2.536271024577331,
    1.421305437503916, -0.0, -2.6932079339274537, 1.5090455370156788, -0.0,
    -2.6932079339274537, 1.5090455370156788, -0.0, -2.8499600576711197,
    1.596658308047201, -0.0, -2.8499600576711197, 1.596658308047201, -0.0,
    -3.0065275786894463, 1.6841439298000525, -0.0, -3.0065275786894463,
    1.6841439298000525, -0.0, -3.16291067975938, 1.7715025812336165, -0.0,
    -3.16291067975938, 1.7715025812336165, 0.0, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.1599184353223182, -0.089801925435736424, 0.0,
    0.31964914489275253, -0.17947362197539324, 0.0, 0.31964914489275253,
    -0.17947362197539324, 0.0, 0.47919231323908462, -0.26901527273781028, 0.0,
    0.47919231323908462, -0.26901527273781028, 0.0, 0.63854812478731549,
    -0.35842706059465185, 0.0, 0.63854812478731549, -0.35842706059465185, 0.0,
    0.797716763861511, -0.44770916817072121, 0.0, 0.797716763861511,
    -0.44770916817072121, 0.0, 0.956698414683647, -0.53686177784427469, 0.0,
    0.956698414683647, -0.53686177784427469, 0.0, 1.1154932613734565,
    -0.62588507174733543, 0.0, 1.1154932613734565, -0.62588507174733543, 0.0,
    1.2741014879482764, -0.71477923176600711, 0.0, 1.2741014879482764,
    -0.71477923176600711, 0.0, 1.4325232783228956, -0.80354443954078658, 0.0,
    1.4325232783228956, -0.80354443954078658, 0.0, 1.5907588163094044,
    -0.89218087646687683, 0.0, 1.5907588163094044, -0.89218087646687683, 0.0,
    1.748808285617043, -0.98068872369449933, 0.0, 1.748808285617043,
    -0.98068872369449933, 0.0, 1.9066718698520528, -1.069068162129206, 0.0,
    1.9066718698520528, -1.069068162129206, 0.0, 2.0643497525175265,
    -1.1573193724321906, 0.0, 2.0643497525175265, -1.1573193724321906, 0.0,
    2.221842117013261, -1.2454425350206004, 0.0, 2.221842117013261,
    -1.2454425350206004, 0.0, 2.3791491466356085, -1.3334378300678469, 0.0,
    2.3791491466356085, -1.3334378300678469, 0.0, 2.536271024577331,
    -1.421305437503916, 0.0, 2.536271024577331, -1.421305437503916, 0.0,
    2.6932079339274537, -1.5090455370156788, 0.0, 2.6932079339274537,
    -1.5090455370156788, 0.0, 2.8499600576711197, -1.596658308047201, 0.0,
    2.8499600576711197, -1.596658308047201, 0.0, 3.0065275786894463,
    -1.6841439298000525, 0.0, 3.0065275786894463, -1.6841439298000525, 0.0,
    3.16291067975938, -1.7715025812336165, 0.0, 3.16291067975938,
    -1.7715025812336165, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, -0.0,
    0.10536597514810161, -0.1119360731420971, -0.0, 0.10536597514810161,
    -0.1119360731420971, -0.0, 0.21068448591173913, -0.22372201869524067, -0.0,
    0.21068448591173913, -0.22372201869524067, -0.0, 0.31595540965099317,
    -0.3353580258245521, -0.0, 0.31595540965099317, -0.3353580258245521, -0.0,
    0.4211786242426373, -0.446844283480486, -0.0, 0.4211786242426373,
    -0.446844283480486, -0.0, 0.52635400807880228, -0.55818098039902453, -0.0,
    0.52635400807880228, -0.55818098039902453, -0.0, 0.63148144006564322,
    -0.6693683051018724, -0.0, 0.63148144006564322, -0.6693683051018724, -0.0,
    0.73656079962200982, -0.78040644589665065, -0.0, 0.73656079962200982,
    -0.78040644589665065, -0.0, 0.84159196667811931, -0.89129559087709187, -0.0,
    0.84159196667811931, -0.89129559087709187, -0.0, 0.94657482167423235,
    -1.0020359279232338, -0.0, 0.94657482167423235, -1.0020359279232338, -0.0,
    1.051509245559332, -1.1126276447016141, -0.0, 1.051509245559332,
    -1.1126276447016141, -0.0, 1.1563951197898057, -1.2230709286654644, -0.0,
    1.1563951197898057, -1.2230709286654644, -0.0, 1.2612323263281291,
    -1.3333659670549041, -0.0, 1.2612323263281291, -1.3333659670549041, -0.0,
    1.3660207476415547, -1.4435129468971351, -0.0, 1.3660207476415547,
    -1.4435129468971351, -0.0, 1.4707602667008017, -1.5535120550066355, -0.0,
    1.4707602667008017, -1.5535120550066355, -0.0, 1.5754507669787494,
    -1.6633634779853532, -0.0, 1.5754507669787494, -1.6633634779853532, -0.0,
    1.6800921324491334, -1.7730674022229005, -0.0, 1.6800921324491334,
    -1.7730674022229005, -0.0, 1.7846842475852445, -1.8826240138967472, -0.0,
    1.7846842475852445, -1.8826240138967472, -0.0, 1.8892269973586309,
    -1.9920334989724147, -0.0, 1.8892269973586309, -1.9920334989724147, -0.0,
    1.9937202672378032, -2.1012960432036696, -0.0, 1.9937202672378032,
    -2.1012960432036696, -0.0, 2.098163943186941, -2.2104118321327175, -0.0,
    2.098163943186941, -2.2104118321327175, 0.0, -0.10536597514810161,
    0.1119360731420971, 0.0, -0.10536597514810161, 0.1119360731420971, 0.0,
    -0.21068448591173913, 0.22372201869524067, 0.0, -0.21068448591173913,
    0.22372201869524067, 0.0, -0.31595540965099317, 0.3353580258245521, 0.0,
    -0.31595540965099317, 0.3353580258245521, 0.0, -0.4211786242426373,
    0.446844283480486, 0.0, -0.4211786242426373, 0.446844283480486, 0.0,
    -0.52635400807880228, 0.55818098039902453, 0.0, -0.52635400807880228,
    0.55818098039902453, 0.0, -0.63148144006564322, 0.6693683051018724, 0.0,
    -0.63148144006564322, 0.6693683051018724, 0.0, -0.73656079962200982,
    0.78040644589665065, 0.0, -0.73656079962200982, 0.78040644589665065, 0.0,
    -0.84159196667811931, 0.89129559087709187, 0.0, -0.84159196667811931,
    0.89129559087709187, 0.0, -0.94657482167423235, 1.0020359279232338, 0.0,
    -0.94657482167423235, 1.0020359279232338, 0.0, -1.051509245559332,
    1.1126276447016141, 0.0, -1.051509245559332, 1.1126276447016141, 0.0,
    -1.1563951197898057, 1.2230709286654644, 0.0, -1.1563951197898057,
    1.2230709286654644, 0.0, -1.2612323263281291, 1.3333659670549041, 0.0,
    -1.2612323263281291, 1.3333659670549041, 0.0, -1.3660207476415547,
    1.4435129468971351, 0.0, -1.3660207476415547, 1.4435129468971351, 0.0,
    -1.4707602667008017, 1.5535120550066355, 0.0, -1.4707602667008017,
    1.5535120550066355, 0.0, -1.5754507669787494, 1.6633634779853532, 0.0,
    -1.5754507669787494, 1.6633634779853532, 0.0, -1.6800921324491334,
    1.7730674022229005, 0.0, -1.6800921324491334, 1.7730674022229005, 0.0,
    -1.7846842475852445, 1.8826240138967472, 0.0, -1.7846842475852445,
    1.8826240138967472, 0.0, -1.8892269973586309, 1.9920334989724147, 0.0,
    -1.8892269973586309, 1.9920334989724147, 0.0, -1.9937202672378032,
    2.1012960432036696, 0.0, -1.9937202672378032, 2.1012960432036696, 0.0,
    -2.098163943186941, 2.2104118321327175, 0.0, -2.098163943186941,
    2.2104118321327175, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
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
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 };

  static const real_T g_1[984]{ -0.12647183082737606, 0.14396327290074137, -0.0,
    -0.12647183082737606, 0.14396327290074137, -0.0, -0.25266713091650589,
    0.28784633762178591, -0.0, -0.25266713091650589, 0.28784633762178591, -0.0,
    -0.37858657539014817, 0.43164948302538075, -0.0, -0.37858657539014817,
    0.43164948302538075, -0.0, -0.50423083775077449, 0.57537299731355918, -0.0,
    -0.50423083775077449, 0.57537299731355918, -0.0, -0.62960058988446832,
    0.71901716802974069, -0.0, -0.62960058988446832, 0.71901716802974069, -0.0,
    -0.75469650206481409, 0.862582282060328, -0.0, -0.75469650206481409,
    0.862582282060328, -0.0, -0.87951924295677675, 1.0060686256362998, -0.0,
    -0.87951924295677675, 1.0060686256362998, -0.0, -1.0040694796205727,
    1.1494764843347998, -0.0, -1.0040694796205727, 1.1494764843347998, -0.0,
    -1.1283478775155311, 1.2928061430807216, -0.0, -1.1283478775155311,
    1.2928061430807216, -0.0, -1.2523551005039453, 1.4360578861482902, -0.0,
    -1.2523551005039453, 1.4360578861482902, -0.0, -1.3760918108549158,
    1.57923199716264, -0.0, -1.3760918108549158, 1.57923199716264, -0.0,
    -1.4995586692481842, 1.7223287591013881, -0.0, -1.4995586692481842,
    1.7223287591013881, -0.0, -1.6227563347779563, 1.8653484542962038, -0.0,
    -1.6227563347779563, 1.8653484542962038, -0.0, -1.7456854649567182,
    2.0082913644343763, -0.0, -1.7456854649567182, 2.0082913644343763, -0.0,
    -1.8683467157190412, 2.1511577705603755, -0.0, -1.8683467157190412,
    2.1511577705603755, -0.0, -1.9907407414253795, 2.2939479530774123, -0.0,
    -1.9907407414253795, 2.2939479530774123, -0.0, -2.112868194865857,
    2.4366621917489919, -0.0, -2.112868194865857, 2.4366621917489919, -0.0,
    -2.2347297272640461, 2.5793007657004661, -0.0, -2.2347297272640461,
    2.5793007657004661, -0.0, -2.3563259882807368, 2.7218639534205811, -0.0,
    -2.3563259882807368, 2.7218639534205811, -0.0, -2.477657626017697,
    2.86435203276302, -0.0, -2.477657626017697, 2.86435203276302, -0.0,
    0.12647183082737606, -0.14396327290074137, 0.0, 0.12647183082737606,
    -0.14396327290074137, 0.0, 0.25266713091650589, -0.28784633762178591, 0.0,
    0.25266713091650589, -0.28784633762178591, 0.0, 0.37858657539014817,
    -0.43164948302538075, 0.0, 0.37858657539014817, -0.43164948302538075, 0.0,
    0.50423083775077449, -0.57537299731355918, 0.0, 0.50423083775077449,
    -0.57537299731355918, 0.0, 0.62960058988446832, -0.71901716802974069, 0.0,
    0.62960058988446832, -0.71901716802974069, 0.0, 0.75469650206481409,
    -0.862582282060328, 0.0, 0.75469650206481409, -0.862582282060328, 0.0,
    0.87951924295677675, -1.0060686256362998, 0.0, 0.87951924295677675,
    -1.0060686256362998, 0.0, 1.0040694796205727, -1.1494764843347998, 0.0,
    1.0040694796205727, -1.1494764843347998, 0.0, 1.1283478775155311,
    -1.2928061430807216, 0.0, 1.1283478775155311, -1.2928061430807216, 0.0,
    1.2523551005039453, -1.4360578861482902, 0.0, 1.2523551005039453,
    -1.4360578861482902, 0.0, 1.3760918108549158, -1.57923199716264, 0.0,
    1.3760918108549158, -1.57923199716264, 0.0, 1.4995586692481842,
    -1.7223287591013881, 0.0, 1.4995586692481842, -1.7223287591013881, 0.0,
    1.6227563347779563, -1.8653484542962038, 0.0, 1.6227563347779563,
    -1.8653484542962038, 0.0, 1.7456854649567182, -2.0082913644343763, 0.0,
    1.7456854649567182, -2.0082913644343763, 0.0, 1.8683467157190412,
    -2.1511577705603755, 0.0, 1.8683467157190412, -2.1511577705603755, 0.0,
    1.9907407414253795, -2.2939479530774123, 0.0, 1.9907407414253795,
    -2.2939479530774123, 0.0, 2.112868194865857, -2.4366621917489919, 0.0,
    2.112868194865857, -2.4366621917489919, 0.0, 2.2347297272640461,
    -2.5793007657004661, 0.0, 2.2347297272640461, -2.5793007657004661, 0.0,
    2.3563259882807368, -2.7218639534205811, 0.0, 2.3563259882807368,
    -2.7218639534205811, 0.0, 2.477657626017697, -2.86435203276302, 0.0,
    2.477657626017697, -2.86435203276302, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.18891678043810189, -0.2603991729373682, -0.0, 0.18891678043810189,
    -0.2603991729373682, -0.0, 0.37743412338800142, -0.52070047841969547, -0.0,
    0.37743412338800142, -0.52070047841969547, -0.0, 0.56555300944331421,
    -0.78090434240065543, -0.0, 0.56555300944331421, -0.78090434240065543, -0.0,
    0.7532744168461627, -1.0410111898780745, -0.0, 0.7532744168461627,
    -1.0410111898780745, -0.0, 0.94059932149283521, -1.3010214448962558, -0.0,
    0.94059932149283521, -1.3010214448962558, -0.0, 1.1275286969394305,
    -1.560935530548297, -0.0, 1.1275286969394305, -1.560935530548297, -0.0,
    1.3140635144074895, -1.8207538689784042, -0.0, 1.3140635144074895,
    -1.8207538689784042, -0.0, 1.5002047427896132, -2.0804768813841976, -0.0,
    1.5002047427896132, -2.0804768813841976, -0.0, 1.6859533486550675,
    -2.3401049880190148, -0.0, 1.6859533486550675, -2.3401049880190148, -0.0,
    1.8713102962553732, -2.5996386081942049, -0.0, 1.8713102962553732,
    -2.5996386081942049, -0.0, 2.0562765475298841, -2.8590781602814213, -0.0,
    2.0562765475298841, -2.8590781602814213, -0.0, 2.2408530621113516,
    -3.1184240617149053, -0.0, 2.2408530621113516, -3.1184240617149053, -0.0,
    2.4250407973314743, -3.3776767289937664, -0.0, 2.4250407973314743,
    -3.3776767289937664, -0.0, 2.6088407082264355, -3.6368365776842562, -0.0,
    2.6088407082264355, -3.6368365776842562, -0.0, 2.7922537475424285,
    -3.8959040224220369, -0.0, 2.7922537475424285, -3.8959040224220369, -0.0,
    2.975280865741166, -4.154879476914445, -0.0, 2.975280865741166,
    -4.154879476914445, -0.0, 3.1579230110053778, -4.4137633539427483, -0.0,
    3.1579230110053778, -4.4137633539427483, -0.0, 3.3401811292442956,
    -4.6725560653643994, -0.0, 3.3401811292442956, -4.6725560653643994, -0.0,
    3.5220561640991237, -4.9312580221152817, -0.0, 3.5220561640991237,
    -4.9312580221152817, -0.0, 3.7035490569484968, -5.1898696342119521, -0.0,
    3.7035490569484968, -5.1898696342119521, -0.0, -0.18891678043810189,
    0.2603991729373682, 0.0, -0.18891678043810189, 0.2603991729373682, 0.0,
    -0.37743412338800142, 0.52070047841969547, 0.0, -0.37743412338800142,
    0.52070047841969547, 0.0, -0.56555300944331421, 0.78090434240065543, 0.0,
    -0.56555300944331421, 0.78090434240065543, 0.0, -0.7532744168461627,
    1.0410111898780745, 0.0, -0.7532744168461627, 1.0410111898780745, 0.0,
    -0.94059932149283521, 1.3010214448962558, 0.0, -0.94059932149283521,
    1.3010214448962558, 0.0, -1.1275286969394305, 1.560935530548297, 0.0,
    -1.1275286969394305, 1.560935530548297, 0.0, -1.3140635144074895,
    1.8207538689784042, 0.0, -1.3140635144074895, 1.8207538689784042, 0.0,
    -1.5002047427896132, 2.0804768813841976, 0.0, -1.5002047427896132,
    2.0804768813841976, 0.0, -1.6859533486550675, 2.3401049880190148, 0.0,
    -1.6859533486550675, 2.3401049880190148, 0.0, -1.8713102962553732,
    2.5996386081942049, 0.0, -1.8713102962553732, 2.5996386081942049, 0.0,
    -2.0562765475298841, 2.8590781602814213, 0.0, -2.0562765475298841,
    2.8590781602814213, 0.0, -2.2408530621113516, 3.1184240617149053, 0.0,
    -2.2408530621113516, 3.1184240617149053, 0.0, -2.4250407973314743,
    3.3776767289937664, 0.0, -2.4250407973314743, 3.3776767289937664, 0.0,
    -2.6088407082264355, 3.6368365776842562, 0.0, -2.6088407082264355,
    3.6368365776842562, 0.0, -2.7922537475424285, 3.8959040224220369, 0.0,
    -2.7922537475424285, 3.8959040224220369, 0.0, -2.975280865741166,
    4.154879476914445, 0.0, -2.975280865741166, 4.154879476914445, 0.0,
    -3.1579230110053778, 4.4137633539427483, 0.0, -3.1579230110053778,
    4.4137633539427483, 0.0, -3.3401811292442956, 4.6725560653643994, 0.0,
    -3.3401811292442956, 4.6725560653643994, 0.0, -3.5220561640991237,
    4.9312580221152817, 0.0, -3.5220561640991237, 4.9312580221152817, 0.0,
    -3.7035490569484968, 5.1898696342119521, 0.0, -3.7035490569484968,
    5.1898696342119521, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.13919179616088076, 0.07490340326547, -0.0, 0.13919179616088076,
    0.07490340326547, -0.0, 0.27800912838762087, 0.15000797759616485, -0.0,
    0.27800912838762087, 0.15000797759616485, -0.0, 0.41645288308364015,
    0.22531337660988532, -0.0, 0.41645288308364015, 0.22531337660988532, -0.0,
    0.55452394451511977, 0.30081925480716887, -0.0, 0.55452394451511977,
    0.30081925480716887, -0.0, 0.69222319481614147, 0.37652526756918425, -0.0,
    0.69222319481614147, 0.37652526756918425, -0.0, 0.82955151399381344,
    0.4524310711556313, -0.0, 0.82955151399381344, 0.4524310711556313, -0.0,
    0.96650977993338516, 0.52853632270264572, -0.0, 0.96650977993338516,
    0.52853632270264572, -0.0, 1.1030988684033489, 0.60484068022070858, -0.0,
    1.1030988684033489, 0.60484068022070858, -0.0, 1.2393196530605293,
    0.68134380259256155, -0.0, 1.2393196530605293, 0.68134380259256155, -0.0,
    1.375173005455161, 0.75804534957112635, -0.0, 1.375173005455161,
    0.75804534957112635, -0.0, 1.5106597950359535, 0.83494498177742982, -0.0,
    1.5106597950359535, 0.83494498177742982, -0.0, 1.6457808891551446,
    0.91204236069853384, -0.0, 1.6457808891551446, 0.91204236069853384, -0.0,
    1.7805371530735408, 0.98933714868547, -0.0, 1.7805371530735408,
    0.98933714868547, -0.0, 1.914929449965546, 1.0668290089511794, -0.0,
    1.914929449965546, 1.0668290089511794, -0.0, 2.048958640924178,
    1.1445176055684576, -0.0, 2.048958640924178, 1.1445176055684576, -0.0,
    2.1826255849660736, 1.2224026034679043, -0.0, 2.1826255849660736,
    1.2224026034679043, -0.0, 2.315931139036481, 1.3004836684358778, -0.0,
    2.315931139036481, 1.3004836684358778, -0.0, 2.4488761580142389,
    1.3787604671124549, -0.0, 2.4488761580142389, 1.3787604671124549, -0.0,
    2.5814614947167467, 1.4572326669893951, -0.0, 2.5814614947167467,
    1.4572326669893951, -0.0, 2.7136879999049208, 1.53589993640811, -0.0,
    2.7136879999049208, 1.53589993640811, -0.0, -0.13919179616088076,
    -0.07490340326547, 0.0, -0.13919179616088076, -0.07490340326547, 0.0,
    -0.27800912838762087, -0.15000797759616485, 0.0, -0.27800912838762087,
    -0.15000797759616485, 0.0, -0.41645288308364015, -0.22531337660988532, 0.0,
    -0.41645288308364015, -0.22531337660988532, 0.0, -0.55452394451511977,
    -0.30081925480716887, 0.0, -0.55452394451511977, -0.30081925480716887, 0.0,
    -0.69222319481614147, -0.37652526756918425, 0.0, -0.69222319481614147,
    -0.37652526756918425, 0.0, -0.82955151399381344, -0.4524310711556313, 0.0,
    -0.82955151399381344, -0.4524310711556313, 0.0, -0.96650977993338516,
    -0.52853632270264572, 0.0, -0.96650977993338516, -0.52853632270264572, 0.0,
    -1.1030988684033489, -0.60484068022070858, 0.0, -1.1030988684033489,
    -0.60484068022070858, 0.0, -1.2393196530605293, -0.68134380259256155, 0.0,
    -1.2393196530605293, -0.68134380259256155, 0.0, -1.375173005455161,
    -0.75804534957112635, 0.0, -1.375173005455161, -0.75804534957112635, 0.0,
    -1.5106597950359535, -0.83494498177742982, 0.0, -1.5106597950359535,
    -0.83494498177742982, 0.0, -1.6457808891551446, -0.91204236069853384, 0.0,
    -1.6457808891551446, -0.91204236069853384, 0.0, -1.7805371530735408,
    -0.98933714868547, 0.0, -1.7805371530735408, -0.98933714868547, 0.0,
    -1.914929449965546, -1.0668290089511794, 0.0, -1.914929449965546,
    -1.0668290089511794, 0.0, -2.048958640924178, -1.1445176055684576, 0.0,
    -2.048958640924178, -1.1445176055684576, 0.0, -2.1826255849660736,
    -1.2224026034679043, 0.0, -2.1826255849660736, -1.2224026034679043, 0.0,
    -2.315931139036481, -1.3004836684358778, 0.0, -2.315931139036481,
    -1.3004836684358778, 0.0, -2.4488761580142389, -1.3787604671124549, 0.0,
    -2.4488761580142389, -1.3787604671124549, 0.0, -2.5814614947167467,
    -1.4572326669893951, 0.0, -2.5814614947167467, -1.4572326669893951, 0.0,
    -2.7136879999049208, -1.53589993640811, 0.0, -2.7136879999049208,
    -1.53589993640811, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0,
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
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0 };

  static const real_T n[960]{ 0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99736022217199194, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99604294756893919, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.9947274127709651, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99341361548020291, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99210155340182049, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99079122424401689, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.989482625718018, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98817575553807258, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98687061142144838, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.985567191088428, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98426549226230531, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.9829655126693807, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.98166725003895783, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98037070210333943, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97907586659782331, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97778274126069831, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97649132383324055, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97520161205970934, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97391360368734325, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99736022217199194, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99604294756893919, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.9947274127709651, 0.0, 0.0, 0.0, 0.0, 0.0, 0.99341361548020291, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.99210155340182049, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99079122424401689, 0.0, 0.0, 0.0, 0.0, 0.0, 0.989482625718018, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98817575553807258, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.98687061142144838, 0.0, 0.0, 0.0, 0.0, 0.0, 0.985567191088428, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98426549226230531, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.9829655126693807, 0.0, 0.0, 0.0, 0.0, 0.0, 0.98166725003895783, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.98037070210333943, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97907586659782331, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97778274126069831, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.97649132383324055, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.97520161205970934, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97391360368734325, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T e[738]{ -0.3579347593572087, -0.0, -0.0,
    -0.3579347593572087, -0.0, -0.0, -0.71539677240113031, -0.0, -0.0,
    -0.71539677240113031, -0.0, -0.0, -1.0723866635167145, -0.0, -0.0,
    -1.0723866635167145, -0.0, -0.0, -1.4289050562642476, -0.0, -0.0,
    -1.4289050562642476, -0.0, -0.0, -1.7849525733804419, -0.0, -0.0,
    -1.7849525733804419, -0.0, -0.0, -2.1405298367795229, -0.0, -0.0,
    -2.1405298367795229, -0.0, -0.0, -2.4956374675543165, -0.0, -0.0,
    -2.4956374675543165, -0.0, -0.0, -2.8502760859773328, -0.0, -0.0,
    -2.8502760859773328, -0.0, -0.0, -3.2044463115018504, -0.0, -0.0,
    -3.2044463115018504, -0.0, -0.0, -3.5581487627629982, -0.0, -0.0,
    -3.5581487627629982, -0.0, -0.0, -3.9113840575788359, -0.0, -0.0,
    -3.9113840575788359, -0.0, -0.0, -4.2641528129514326, -0.0, -0.0,
    -4.2641528129514326, -0.0, -0.0, -4.616455645067945, -0.0, -0.0,
    -4.616455645067945, -0.0, -0.0, -4.9682931693016954, -0.0, -0.0,
    -4.9682931693016954, -0.0, -0.0, -5.3196660002132425, -0.0, -0.0,
    -5.3196660002132425, -0.0, -0.0, -5.6705747515514595, -0.0, -0.0,
    -5.6705747515514595, -0.0, -0.0, -6.021020036254602, -0.0, -0.0,
    -6.021020036254602, -0.0, -0.0, -6.3710024664513822, -0.0, -0.0,
    -6.3710024664513822, -0.0, -0.0, -6.7205226534620355, -0.0, -0.0,
    -6.7205226534620355, -0.0, -0.0, -7.06958120779939, -0.0, -0.0,
    -7.06958120779939, -0.0, -0.0, 0.3579347593572087, 0.0, 0.0,
    0.3579347593572087, 0.0, 0.0, 0.71539677240113031, 0.0, 0.0,
    0.71539677240113031, 0.0, 0.0, 1.0723866635167145, 0.0, 0.0,
    1.0723866635167145, 0.0, 0.0, 1.4289050562642476, 0.0, 0.0,
    1.4289050562642476, 0.0, 0.0, 1.7849525733804419, 0.0, 0.0,
    1.7849525733804419, 0.0, 0.0, 2.1405298367795229, 0.0, 0.0,
    2.1405298367795229, 0.0, 0.0, 2.4956374675543165, 0.0, 0.0,
    2.4956374675543165, 0.0, 0.0, 2.8502760859773328, 0.0, 0.0,
    2.8502760859773328, 0.0, 0.0, 3.2044463115018504, 0.0, 0.0,
    3.2044463115018504, 0.0, 0.0, 3.5581487627629982, 0.0, 0.0,
    3.5581487627629982, 0.0, 0.0, 3.9113840575788359, 0.0, 0.0,
    3.9113840575788359, 0.0, 0.0, 4.2641528129514326, 0.0, 0.0,
    4.2641528129514326, 0.0, 0.0, 4.616455645067945, 0.0, 0.0, 4.616455645067945,
    0.0, 0.0, 4.9682931693016954, 0.0, 0.0, 4.9682931693016954, 0.0, 0.0,
    5.3196660002132425, 0.0, 0.0, 5.3196660002132425, 0.0, 0.0,
    5.6705747515514595, 0.0, 0.0, 5.6705747515514595, 0.0, 0.0,
    6.021020036254602, 0.0, 0.0, 6.021020036254602, 0.0, 0.0, 6.3710024664513822,
    0.0, 0.0, 6.3710024664513822, 0.0, 0.0, 6.7205226534620355, 0.0, 0.0,
    6.7205226534620355, 0.0, 0.0, 7.06958120779939, 0.0, 0.0, 7.06958120779939,
    0.0, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.042560663779893933, -0.0, -0.0,
    0.042560663779893933, -0.0, -0.0, 0.0850651150898698, -0.0, -0.0,
    0.0850651150898698, -0.0, -0.0, 0.12751342817317227, -0.0, -0.0,
    0.12751342817317227, -0.0, -0.0, 0.16990567717498842, -0.0, -0.0,
    0.16990567717498842, -0.0, -0.0, 0.21224193614257725, -0.0, -0.0,
    0.21224193614257725, -0.0, -0.0, 0.254522279025399, -0.0, -0.0,
    0.254522279025399, -0.0, -0.0, 0.29674677967524432, -0.0, -0.0,
    0.29674677967524432, -0.0, -0.0, 0.33891551184636343, -0.0, -0.0,
    0.33891551184636343, -0.0, -0.0, 0.38102854919559465, -0.0, -0.0,
    0.38102854919559465, -0.0, -0.0, 0.4230859652824932, -0.0, -0.0,
    0.4230859652824932, -0.0, -0.0, 0.46508783356945982, -0.0, -0.0,
    0.46508783356945982, -0.0, -0.0, 0.50703422742186888, -0.0, -0.0,
    0.50703422742186888, -0.0, -0.0, 0.54892522010819667, -0.0, -0.0,
    0.54892522010819667, -0.0, -0.0, 0.59076088480014921, -0.0, -0.0,
    0.59076088480014921, -0.0, -0.0, 0.6325412945727904, -0.0, -0.0,
    0.6325412945727904, -0.0, -0.0, 0.67426652240466922, -0.0, -0.0,
    0.67426652240466922, -0.0, -0.0, 0.71593664117794742, -0.0, -0.0,
    0.71593664117794742, -0.0, -0.0, 0.757551723678527, -0.0, -0.0,
    0.757551723678527, -0.0, -0.0, 0.79911184259617707, -0.0, -0.0,
    0.79911184259617707, -0.0, -0.0, 0.840617070524661, -0.0, -0.0,
    0.840617070524661, -0.0, -0.0, -0.042560663779893933, 0.0, 0.0,
    -0.042560663779893933, 0.0, 0.0, -0.0850651150898698, 0.0, 0.0,
    -0.0850651150898698, 0.0, 0.0, -0.12751342817317227, 0.0, 0.0,
    -0.12751342817317227, 0.0, 0.0, -0.16990567717498842, 0.0, 0.0,
    -0.16990567717498842, 0.0, 0.0, -0.21224193614257725, 0.0, 0.0,
    -0.21224193614257725, 0.0, 0.0, -0.254522279025399, 0.0, 0.0,
    -0.254522279025399, 0.0, 0.0, -0.29674677967524432, 0.0, 0.0,
    -0.29674677967524432, 0.0, 0.0, -0.33891551184636343, 0.0, 0.0,
    -0.33891551184636343, 0.0, 0.0, -0.38102854919559465, 0.0, 0.0,
    -0.38102854919559465, 0.0, 0.0, -0.4230859652824932, 0.0, 0.0,
    -0.4230859652824932, 0.0, 0.0, -0.46508783356945982, 0.0, 0.0,
    -0.46508783356945982, 0.0, 0.0, -0.50703422742186888, 0.0, 0.0,
    -0.50703422742186888, 0.0, 0.0, -0.54892522010819667, 0.0, 0.0,
    -0.54892522010819667, 0.0, 0.0, -0.59076088480014921, 0.0, 0.0,
    -0.59076088480014921, 0.0, 0.0, -0.6325412945727904, 0.0, 0.0,
    -0.6325412945727904, 0.0, 0.0, -0.67426652240466922, 0.0, 0.0,
    -0.67426652240466922, 0.0, 0.0, -0.71593664117794742, 0.0, 0.0,
    -0.71593664117794742, 0.0, 0.0, -0.757551723678527, 0.0, 0.0,
    -0.757551723678527, 0.0, 0.0, -0.79911184259617707, 0.0, 0.0,
    -0.79911184259617707, 0.0, 0.0, -0.840617070524661, 0.0, 0.0,
    -0.840617070524661, 0.0, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.0099627081344734868, -0.0, -0.0, 0.0099627081344734868, -0.0, -0.0,
    0.01991225791140331, -0.0, -0.0, 0.01991225791140331, -0.0, -0.0,
    0.029848666709836498, -0.0, -0.0, 0.029848666709836498, -0.0, -0.0,
    0.039771951885866519, -0.0, -0.0, 0.039771951885866519, -0.0, -0.0,
    0.049682130772663577, -0.0, -0.0, 0.049682130772663577, -0.0, -0.0,
    0.059579220680504912, -0.0, -0.0, 0.059579220680504912, -0.0, -0.0,
    0.069463238896805016, -0.0, -0.0, 0.069463238896805016, -0.0, -0.0,
    0.079334202686145827, -0.0, -0.0, 0.079334202686145827, -0.0, -0.0,
    0.089192129290306912, -0.0, -0.0, 0.089192129290306912, -0.0, -0.0,
    0.099037035928295547, -0.0, -0.0, 0.099037035928295547, -0.0, -0.0,
    0.10886893979637684, -0.0, -0.0, 0.10886893979637684, -0.0, -0.0,
    0.11868785806810371, -0.0, -0.0, 0.11868785806810371, -0.0, -0.0,
    0.12849380789434692, -0.0, -0.0, 0.12849380789434692, -0.0, -0.0,
    0.13828680640332505, -0.0, -0.0, 0.13828680640332505, -0.0, -0.0,
    0.1480668707006344, -0.0, -0.0, 0.1480668707006344, -0.0, -0.0,
    0.15783401786927884, -0.0, -0.0, 0.15783401786927884, -0.0, -0.0,
    0.16758826496969964, -0.0, -0.0, 0.16758826496969964, -0.0, -0.0,
    0.17732962903980537, -0.0, -0.0, 0.17732962903980537, -0.0, -0.0,
    0.18705812709500158, -0.0, -0.0, 0.18705812709500158, -0.0, -0.0,
    0.1967737761282205, -0.0, -0.0, 0.1967737761282205, -0.0, -0.0,
    -0.0099627081344734868, 0.0, 0.0, -0.0099627081344734868, 0.0, 0.0,
    -0.01991225791140331, 0.0, 0.0, -0.01991225791140331, 0.0, 0.0,
    -0.029848666709836498, 0.0, 0.0, -0.029848666709836498, 0.0, 0.0,
    -0.039771951885866519, 0.0, 0.0, -0.039771951885866519, 0.0, 0.0,
    -0.049682130772663577, 0.0, 0.0, -0.049682130772663577, 0.0, 0.0,
    -0.059579220680504912, 0.0, 0.0, -0.059579220680504912, 0.0, 0.0,
    -0.069463238896805016, 0.0, 0.0, -0.069463238896805016, 0.0, 0.0,
    -0.079334202686145827, 0.0, 0.0, -0.079334202686145827, 0.0, 0.0,
    -0.089192129290306912, 0.0, 0.0, -0.089192129290306912, 0.0, 0.0,
    -0.099037035928295547, 0.0, 0.0, -0.099037035928295547, 0.0, 0.0,
    -0.10886893979637684, 0.0, 0.0, -0.10886893979637684, 0.0, 0.0,
    -0.11868785806810371, 0.0, 0.0, -0.11868785806810371, 0.0, 0.0,
    -0.12849380789434692, 0.0, 0.0, -0.12849380789434692, 0.0, 0.0,
    -0.13828680640332505, 0.0, 0.0, -0.13828680640332505, 0.0, 0.0,
    -0.1480668707006344, 0.0, 0.0, -0.1480668707006344, 0.0, 0.0,
    -0.15783401786927884, 0.0, 0.0, -0.15783401786927884, 0.0, 0.0,
    -0.16758826496969964, 0.0, 0.0, -0.16758826496969964, 0.0, 0.0,
    -0.17732962903980537, 0.0, 0.0, -0.17732962903980537, 0.0, 0.0,
    -0.18705812709500158, 0.0, 0.0, -0.18705812709500158, 0.0, 0.0,
    -0.1967737761282205, 0.0, 0.0, -0.1967737761282205, 0.0, 0.0, -0.0, -0.0,
    -1.0, 0.0, 0.0, 1.0 };

  static const real_T e_0[738]{ -0.0, 0.1888968412469153, 0.060492001339316988,
    -0.0, 0.1888968412469153, 0.060492001339316988, -0.0, 0.3773313347058973,
    0.1209347881183509, -0.0, 0.3773313347058973, 0.1209347881183509, -0.0,
    0.56530446924504474, 0.18132836030302127, -0.0, 0.56530446924504474,
    0.18132836030302127, -0.0, 0.75281723176956472, 0.24167271798767209, -0.0,
    0.75281723176956472, 0.24167271798767209, -0.0, 0.939870607225484,
    0.30196786139466048, -0.0, 0.939870607225484, 0.30196786139466048, -0.0,
    1.1264655786033526, 0.36221379087394617, -0.0, 1.1264655786033526,
    0.36221379087394617, -0.0, 1.3126031269419416, 0.42241050690268223, -0.0,
    1.3126031269419416, 0.42241050690268223, -0.0, 1.4982842313319336,
    0.48255801008480659, -0.0, 1.4982842313319336, 0.48255801008480659, -0.0,
    1.6835098689196064, 0.54265630115063446, -0.0, 1.6835098689196064,
    0.54265630115063446, -0.0, 1.8682810149105107, 0.6027053809564521, -0.0,
    1.8682810149105107, 0.6027053809564521, -0.0, 2.05259864257314,
    0.66270525048411089, -0.0, 2.05259864257314, 0.66270525048411089, -0.0,
    2.2364637232425952, 0.72265591084062319, -0.0, 2.2364637232425952,
    0.72265591084062319, -0.0, 2.41987722632424, 0.78255736325775849, -0.0,
    2.41987722632424, 0.78255736325775849, -0.0, 2.6028401192973529,
    0.842409609091641, -0.0, 2.6028401192973529, 0.842409609091641, -0.0,
    2.7853533677187707, 0.90221264982234783, -0.0, 2.7853533677187707,
    0.90221264982234783, -0.0, 2.9674179352265249, 0.96196648705350851, -0.0,
    2.9674179352265249, 0.96196648705350851, -0.0, 3.1490347835434722,
    1.0216711225119051, -0.0, 3.1490347835434722, 1.0216711225119051, -0.0,
    3.3302048724809192, 1.081326558047073, -0.0, 3.3302048724809192,
    1.081326558047073, -0.0, 3.5109291599422385, 1.1409327956309046, -0.0,
    3.5109291599422385, 1.1409327956309046, -0.0, 3.69120860192648,
    1.2004898373572503, -0.0, 3.69120860192648, 1.2004898373572503, 0.0,
    -0.1888968412469153, -0.060492001339316988, 0.0, -0.1888968412469153,
    -0.060492001339316988, 0.0, -0.3773313347058973, -0.1209347881183509, 0.0,
    -0.3773313347058973, -0.1209347881183509, 0.0, -0.56530446924504474,
    -0.18132836030302127, 0.0, -0.56530446924504474, -0.18132836030302127, 0.0,
    -0.75281723176956472, -0.24167271798767209, 0.0, -0.75281723176956472,
    -0.24167271798767209, 0.0, -0.939870607225484, -0.30196786139466048, 0.0,
    -0.939870607225484, -0.30196786139466048, 0.0, -1.1264655786033526,
    -0.36221379087394617, 0.0, -1.1264655786033526, -0.36221379087394617, 0.0,
    -1.3126031269419416, -0.42241050690268223, 0.0, -1.3126031269419416,
    -0.42241050690268223, 0.0, -1.4982842313319336, -0.48255801008480659, 0.0,
    -1.4982842313319336, -0.48255801008480659, 0.0, -1.6835098689196064,
    -0.54265630115063446, 0.0, -1.6835098689196064, -0.54265630115063446, 0.0,
    -1.8682810149105107, -0.6027053809564521, 0.0, -1.8682810149105107,
    -0.6027053809564521, 0.0, -2.05259864257314, -0.66270525048411089, 0.0,
    -2.05259864257314, -0.66270525048411089, 0.0, -2.2364637232425952,
    -0.72265591084062319, 0.0, -2.2364637232425952, -0.72265591084062319, 0.0,
    -2.41987722632424, -0.78255736325775849, 0.0, -2.41987722632424,
    -0.78255736325775849, 0.0, -2.6028401192973529, -0.842409609091641, 0.0,
    -2.6028401192973529, -0.842409609091641, 0.0, -2.7853533677187707,
    -0.90221264982234783, 0.0, -2.7853533677187707, -0.90221264982234783, 0.0,
    -2.9674179352265249, -0.96196648705350851, 0.0, -2.9674179352265249,
    -0.96196648705350851, 0.0, -3.1490347835434722, -1.0216711225119051, 0.0,
    -3.1490347835434722, -1.0216711225119051, 0.0, -3.3302048724809192,
    -1.081326558047073, 0.0, -3.3302048724809192, -1.081326558047073, 0.0,
    -3.5109291599422385, -1.1409327956309046, 0.0, -3.5109291599422385,
    -1.1409327956309046, 0.0, -3.69120860192648, -1.2004898373572503, 0.0,
    -3.69120860192648, -1.2004898373572503, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    -0.0, -0.1599184353223182, 0.089801925435736424, -0.0, -0.1599184353223182,
    0.089801925435736424, -0.0, -0.31964914489275253, 0.17947362197539324, -0.0,
    -0.31964914489275253, 0.17947362197539324, -0.0, -0.47919231323908462,
    0.26901527273781028, -0.0, -0.47919231323908462, 0.26901527273781028, -0.0,
    -0.63854812478731549, 0.35842706059465185, -0.0, -0.63854812478731549,
    0.35842706059465185, -0.0, -0.797716763861511, 0.44770916817072121, -0.0,
    -0.797716763861511, 0.44770916817072121, -0.0, -0.956698414683647,
    0.53686177784427469, -0.0, -0.956698414683647, 0.53686177784427469, -0.0,
    -1.1154932613734565, 0.62588507174733543, -0.0, -1.1154932613734565,
    0.62588507174733543, -0.0, -1.2741014879482764, 0.71477923176600711, -0.0,
    -1.2741014879482764, 0.71477923176600711, -0.0, -1.4325232783228956,
    0.80354443954078658, -0.0, -1.4325232783228956, 0.80354443954078658, -0.0,
    -1.5907588163094044, 0.89218087646687683, -0.0, -1.5907588163094044,
    0.89218087646687683, -0.0, -1.748808285617043, 0.98068872369449933, -0.0,
    -1.748808285617043, 0.98068872369449933, -0.0, -1.9066718698520528,
    1.069068162129206, -0.0, -1.9066718698520528, 1.069068162129206, -0.0,
    -2.0643497525175265, 1.1573193724321906, -0.0, -2.0643497525175265,
    1.1573193724321906, -0.0, -2.221842117013261, 1.2454425350206004, -0.0,
    -2.221842117013261, 1.2454425350206004, -0.0, -2.3791491466356085,
    1.3334378300678469, -0.0, -2.3791491466356085, 1.3334378300678469, -0.0,
    -2.536271024577331, 1.421305437503916, -0.0, -2.536271024577331,
    1.421305437503916, -0.0, -2.6932079339274537, 1.5090455370156788, -0.0,
    -2.6932079339274537, 1.5090455370156788, -0.0, -2.8499600576711197,
    1.596658308047201, -0.0, -2.8499600576711197, 1.596658308047201, -0.0,
    -3.0065275786894463, 1.6841439298000525, -0.0, -3.0065275786894463,
    1.6841439298000525, -0.0, -3.16291067975938, 1.7715025812336165, -0.0,
    -3.16291067975938, 1.7715025812336165, 0.0, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.1599184353223182, -0.089801925435736424, 0.0,
    0.31964914489275253, -0.17947362197539324, 0.0, 0.31964914489275253,
    -0.17947362197539324, 0.0, 0.47919231323908462, -0.26901527273781028, 0.0,
    0.47919231323908462, -0.26901527273781028, 0.0, 0.63854812478731549,
    -0.35842706059465185, 0.0, 0.63854812478731549, -0.35842706059465185, 0.0,
    0.797716763861511, -0.44770916817072121, 0.0, 0.797716763861511,
    -0.44770916817072121, 0.0, 0.956698414683647, -0.53686177784427469, 0.0,
    0.956698414683647, -0.53686177784427469, 0.0, 1.1154932613734565,
    -0.62588507174733543, 0.0, 1.1154932613734565, -0.62588507174733543, 0.0,
    1.2741014879482764, -0.71477923176600711, 0.0, 1.2741014879482764,
    -0.71477923176600711, 0.0, 1.4325232783228956, -0.80354443954078658, 0.0,
    1.4325232783228956, -0.80354443954078658, 0.0, 1.5907588163094044,
    -0.89218087646687683, 0.0, 1.5907588163094044, -0.89218087646687683, 0.0,
    1.748808285617043, -0.98068872369449933, 0.0, 1.748808285617043,
    -0.98068872369449933, 0.0, 1.9066718698520528, -1.069068162129206, 0.0,
    1.9066718698520528, -1.069068162129206, 0.0, 2.0643497525175265,
    -1.1573193724321906, 0.0, 2.0643497525175265, -1.1573193724321906, 0.0,
    2.221842117013261, -1.2454425350206004, 0.0, 2.221842117013261,
    -1.2454425350206004, 0.0, 2.3791491466356085, -1.3334378300678469, 0.0,
    2.3791491466356085, -1.3334378300678469, 0.0, 2.536271024577331,
    -1.421305437503916, 0.0, 2.536271024577331, -1.421305437503916, 0.0,
    2.6932079339274537, -1.5090455370156788, 0.0, 2.6932079339274537,
    -1.5090455370156788, 0.0, 2.8499600576711197, -1.596658308047201, 0.0,
    2.8499600576711197, -1.596658308047201, 0.0, 3.0065275786894463,
    -1.6841439298000525, 0.0, 3.0065275786894463, -1.6841439298000525, 0.0,
    3.16291067975938, -1.7715025812336165, 0.0, 3.16291067975938,
    -1.7715025812336165, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, -0.0,
    0.10536597514810161, -0.1119360731420971, -0.0, 0.10536597514810161,
    -0.1119360731420971, -0.0, 0.21068448591173913, -0.22372201869524067, -0.0,
    0.21068448591173913, -0.22372201869524067, -0.0, 0.31595540965099317,
    -0.3353580258245521, -0.0, 0.31595540965099317, -0.3353580258245521, -0.0,
    0.4211786242426373, -0.446844283480486, -0.0, 0.4211786242426373,
    -0.446844283480486, -0.0, 0.52635400807880228, -0.55818098039902453, -0.0,
    0.52635400807880228, -0.55818098039902453, -0.0, 0.63148144006564322,
    -0.6693683051018724, -0.0, 0.63148144006564322, -0.6693683051018724, -0.0,
    0.73656079962200982, -0.78040644589665065, -0.0, 0.73656079962200982,
    -0.78040644589665065, -0.0, 0.84159196667811931, -0.89129559087709187, -0.0,
    0.84159196667811931, -0.89129559087709187, -0.0, 0.94657482167423235,
    -1.0020359279232338, -0.0, 0.94657482167423235, -1.0020359279232338, -0.0,
    1.051509245559332, -1.1126276447016141, -0.0, 1.051509245559332,
    -1.1126276447016141, -0.0, 1.1563951197898057, -1.2230709286654644, -0.0,
    1.1563951197898057, -1.2230709286654644, -0.0, 1.2612323263281291,
    -1.3333659670549041, -0.0, 1.2612323263281291, -1.3333659670549041, -0.0,
    1.3660207476415547, -1.4435129468971351, -0.0, 1.3660207476415547,
    -1.4435129468971351, -0.0, 1.4707602667008017, -1.5535120550066355, -0.0,
    1.4707602667008017, -1.5535120550066355, -0.0, 1.5754507669787494,
    -1.6633634779853532, -0.0, 1.5754507669787494, -1.6633634779853532, -0.0,
    1.6800921324491334, -1.7730674022229005, -0.0, 1.6800921324491334,
    -1.7730674022229005, -0.0, 1.7846842475852445, -1.8826240138967472, -0.0,
    1.7846842475852445, -1.8826240138967472, -0.0, 1.8892269973586309,
    -1.9920334989724147, -0.0, 1.8892269973586309, -1.9920334989724147, -0.0,
    1.9937202672378032, -2.1012960432036696, -0.0, 1.9937202672378032,
    -2.1012960432036696, -0.0, 2.098163943186941, -2.2104118321327175, -0.0,
    2.098163943186941, -2.2104118321327175, 0.0, -0.10536597514810161,
    0.1119360731420971, 0.0, -0.10536597514810161, 0.1119360731420971, 0.0,
    -0.21068448591173913, 0.22372201869524067, 0.0, -0.21068448591173913,
    0.22372201869524067, 0.0, -0.31595540965099317, 0.3353580258245521, 0.0,
    -0.31595540965099317, 0.3353580258245521, 0.0, -0.4211786242426373,
    0.446844283480486, 0.0, -0.4211786242426373, 0.446844283480486, 0.0,
    -0.52635400807880228, 0.55818098039902453, 0.0, -0.52635400807880228,
    0.55818098039902453, 0.0, -0.63148144006564322, 0.6693683051018724, 0.0,
    -0.63148144006564322, 0.6693683051018724, 0.0, -0.73656079962200982,
    0.78040644589665065, 0.0, -0.73656079962200982, 0.78040644589665065, 0.0,
    -0.84159196667811931, 0.89129559087709187, 0.0, -0.84159196667811931,
    0.89129559087709187, 0.0, -0.94657482167423235, 1.0020359279232338, 0.0,
    -0.94657482167423235, 1.0020359279232338, 0.0, -1.051509245559332,
    1.1126276447016141, 0.0, -1.051509245559332, 1.1126276447016141, 0.0,
    -1.1563951197898057, 1.2230709286654644, 0.0, -1.1563951197898057,
    1.2230709286654644, 0.0, -1.2612323263281291, 1.3333659670549041, 0.0,
    -1.2612323263281291, 1.3333659670549041, 0.0, -1.3660207476415547,
    1.4435129468971351, 0.0, -1.3660207476415547, 1.4435129468971351, 0.0,
    -1.4707602667008017, 1.5535120550066355, 0.0, -1.4707602667008017,
    1.5535120550066355, 0.0, -1.5754507669787494, 1.6633634779853532, 0.0,
    -1.5754507669787494, 1.6633634779853532, 0.0, -1.6800921324491334,
    1.7730674022229005, 0.0, -1.6800921324491334, 1.7730674022229005, 0.0,
    -1.7846842475852445, 1.8826240138967472, 0.0, -1.7846842475852445,
    1.8826240138967472, 0.0, -1.8892269973586309, 1.9920334989724147, 0.0,
    -1.8892269973586309, 1.9920334989724147, 0.0, -1.9937202672378032,
    2.1012960432036696, 0.0, -1.9937202672378032, 2.1012960432036696, 0.0,
    -2.098163943186941, 2.2104118321327175, 0.0, -2.098163943186941,
    2.2104118321327175, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T e_1[738]{ -0.12647183082737606, 0.14396327290074137, -0.0,
    -0.12647183082737606, 0.14396327290074137, -0.0, -0.25266713091650589,
    0.28784633762178591, -0.0, -0.25266713091650589, 0.28784633762178591, -0.0,
    -0.37858657539014817, 0.43164948302538075, -0.0, -0.37858657539014817,
    0.43164948302538075, -0.0, -0.50423083775077449, 0.57537299731355918, -0.0,
    -0.50423083775077449, 0.57537299731355918, -0.0, -0.62960058988446832,
    0.71901716802974069, -0.0, -0.62960058988446832, 0.71901716802974069, -0.0,
    -0.75469650206481409, 0.862582282060328, -0.0, -0.75469650206481409,
    0.862582282060328, -0.0, -0.87951924295677675, 1.0060686256362998, -0.0,
    -0.87951924295677675, 1.0060686256362998, -0.0, -1.0040694796205727,
    1.1494764843347998, -0.0, -1.0040694796205727, 1.1494764843347998, -0.0,
    -1.1283478775155311, 1.2928061430807216, -0.0, -1.1283478775155311,
    1.2928061430807216, -0.0, -1.2523551005039453, 1.4360578861482902, -0.0,
    -1.2523551005039453, 1.4360578861482902, -0.0, -1.3760918108549158,
    1.57923199716264, -0.0, -1.3760918108549158, 1.57923199716264, -0.0,
    -1.4995586692481842, 1.7223287591013881, -0.0, -1.4995586692481842,
    1.7223287591013881, -0.0, -1.6227563347779563, 1.8653484542962038, -0.0,
    -1.6227563347779563, 1.8653484542962038, -0.0, -1.7456854649567182,
    2.0082913644343763, -0.0, -1.7456854649567182, 2.0082913644343763, -0.0,
    -1.8683467157190412, 2.1511577705603755, -0.0, -1.8683467157190412,
    2.1511577705603755, -0.0, -1.9907407414253795, 2.2939479530774123, -0.0,
    -1.9907407414253795, 2.2939479530774123, -0.0, -2.112868194865857,
    2.4366621917489919, -0.0, -2.112868194865857, 2.4366621917489919, -0.0,
    -2.2347297272640461, 2.5793007657004661, -0.0, -2.2347297272640461,
    2.5793007657004661, -0.0, -2.3563259882807368, 2.7218639534205811, -0.0,
    -2.3563259882807368, 2.7218639534205811, -0.0, -2.477657626017697,
    2.86435203276302, -0.0, -2.477657626017697, 2.86435203276302, -0.0,
    0.12647183082737606, -0.14396327290074137, 0.0, 0.12647183082737606,
    -0.14396327290074137, 0.0, 0.25266713091650589, -0.28784633762178591, 0.0,
    0.25266713091650589, -0.28784633762178591, 0.0, 0.37858657539014817,
    -0.43164948302538075, 0.0, 0.37858657539014817, -0.43164948302538075, 0.0,
    0.50423083775077449, -0.57537299731355918, 0.0, 0.50423083775077449,
    -0.57537299731355918, 0.0, 0.62960058988446832, -0.71901716802974069, 0.0,
    0.62960058988446832, -0.71901716802974069, 0.0, 0.75469650206481409,
    -0.862582282060328, 0.0, 0.75469650206481409, -0.862582282060328, 0.0,
    0.87951924295677675, -1.0060686256362998, 0.0, 0.87951924295677675,
    -1.0060686256362998, 0.0, 1.0040694796205727, -1.1494764843347998, 0.0,
    1.0040694796205727, -1.1494764843347998, 0.0, 1.1283478775155311,
    -1.2928061430807216, 0.0, 1.1283478775155311, -1.2928061430807216, 0.0,
    1.2523551005039453, -1.4360578861482902, 0.0, 1.2523551005039453,
    -1.4360578861482902, 0.0, 1.3760918108549158, -1.57923199716264, 0.0,
    1.3760918108549158, -1.57923199716264, 0.0, 1.4995586692481842,
    -1.7223287591013881, 0.0, 1.4995586692481842, -1.7223287591013881, 0.0,
    1.6227563347779563, -1.8653484542962038, 0.0, 1.6227563347779563,
    -1.8653484542962038, 0.0, 1.7456854649567182, -2.0082913644343763, 0.0,
    1.7456854649567182, -2.0082913644343763, 0.0, 1.8683467157190412,
    -2.1511577705603755, 0.0, 1.8683467157190412, -2.1511577705603755, 0.0,
    1.9907407414253795, -2.2939479530774123, 0.0, 1.9907407414253795,
    -2.2939479530774123, 0.0, 2.112868194865857, -2.4366621917489919, 0.0,
    2.112868194865857, -2.4366621917489919, 0.0, 2.2347297272640461,
    -2.5793007657004661, 0.0, 2.2347297272640461, -2.5793007657004661, 0.0,
    2.3563259882807368, -2.7218639534205811, 0.0, 2.3563259882807368,
    -2.7218639534205811, 0.0, 2.477657626017697, -2.86435203276302, 0.0,
    2.477657626017697, -2.86435203276302, 0.0, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.18891678043810189, -0.2603991729373682, -0.0, 0.18891678043810189,
    -0.2603991729373682, -0.0, 0.37743412338800142, -0.52070047841969547, -0.0,
    0.37743412338800142, -0.52070047841969547, -0.0, 0.56555300944331421,
    -0.78090434240065543, -0.0, 0.56555300944331421, -0.78090434240065543, -0.0,
    0.7532744168461627, -1.0410111898780745, -0.0, 0.7532744168461627,
    -1.0410111898780745, -0.0, 0.94059932149283521, -1.3010214448962558, -0.0,
    0.94059932149283521, -1.3010214448962558, -0.0, 1.1275286969394305,
    -1.560935530548297, -0.0, 1.1275286969394305, -1.560935530548297, -0.0,
    1.3140635144074895, -1.8207538689784042, -0.0, 1.3140635144074895,
    -1.8207538689784042, -0.0, 1.5002047427896132, -2.0804768813841976, -0.0,
    1.5002047427896132, -2.0804768813841976, -0.0, 1.6859533486550675,
    -2.3401049880190148, -0.0, 1.6859533486550675, -2.3401049880190148, -0.0,
    1.8713102962553732, -2.5996386081942049, -0.0, 1.8713102962553732,
    -2.5996386081942049, -0.0, 2.0562765475298841, -2.8590781602814213, -0.0,
    2.0562765475298841, -2.8590781602814213, -0.0, 2.2408530621113516,
    -3.1184240617149053, -0.0, 2.2408530621113516, -3.1184240617149053, -0.0,
    2.4250407973314743, -3.3776767289937664, -0.0, 2.4250407973314743,
    -3.3776767289937664, -0.0, 2.6088407082264355, -3.6368365776842562, -0.0,
    2.6088407082264355, -3.6368365776842562, -0.0, 2.7922537475424285,
    -3.8959040224220369, -0.0, 2.7922537475424285, -3.8959040224220369, -0.0,
    2.975280865741166, -4.154879476914445, -0.0, 2.975280865741166,
    -4.154879476914445, -0.0, 3.1579230110053778, -4.4137633539427483, -0.0,
    3.1579230110053778, -4.4137633539427483, -0.0, 3.3401811292442956,
    -4.6725560653643994, -0.0, 3.3401811292442956, -4.6725560653643994, -0.0,
    3.5220561640991237, -4.9312580221152817, -0.0, 3.5220561640991237,
    -4.9312580221152817, -0.0, 3.7035490569484968, -5.1898696342119521, -0.0,
    3.7035490569484968, -5.1898696342119521, -0.0, -0.18891678043810189,
    0.2603991729373682, 0.0, -0.18891678043810189, 0.2603991729373682, 0.0,
    -0.37743412338800142, 0.52070047841969547, 0.0, -0.37743412338800142,
    0.52070047841969547, 0.0, -0.56555300944331421, 0.78090434240065543, 0.0,
    -0.56555300944331421, 0.78090434240065543, 0.0, -0.7532744168461627,
    1.0410111898780745, 0.0, -0.7532744168461627, 1.0410111898780745, 0.0,
    -0.94059932149283521, 1.3010214448962558, 0.0, -0.94059932149283521,
    1.3010214448962558, 0.0, -1.1275286969394305, 1.560935530548297, 0.0,
    -1.1275286969394305, 1.560935530548297, 0.0, -1.3140635144074895,
    1.8207538689784042, 0.0, -1.3140635144074895, 1.8207538689784042, 0.0,
    -1.5002047427896132, 2.0804768813841976, 0.0, -1.5002047427896132,
    2.0804768813841976, 0.0, -1.6859533486550675, 2.3401049880190148, 0.0,
    -1.6859533486550675, 2.3401049880190148, 0.0, -1.8713102962553732,
    2.5996386081942049, 0.0, -1.8713102962553732, 2.5996386081942049, 0.0,
    -2.0562765475298841, 2.8590781602814213, 0.0, -2.0562765475298841,
    2.8590781602814213, 0.0, -2.2408530621113516, 3.1184240617149053, 0.0,
    -2.2408530621113516, 3.1184240617149053, 0.0, -2.4250407973314743,
    3.3776767289937664, 0.0, -2.4250407973314743, 3.3776767289937664, 0.0,
    -2.6088407082264355, 3.6368365776842562, 0.0, -2.6088407082264355,
    3.6368365776842562, 0.0, -2.7922537475424285, 3.8959040224220369, 0.0,
    -2.7922537475424285, 3.8959040224220369, 0.0, -2.975280865741166,
    4.154879476914445, 0.0, -2.975280865741166, 4.154879476914445, 0.0,
    -3.1579230110053778, 4.4137633539427483, 0.0, -3.1579230110053778,
    4.4137633539427483, 0.0, -3.3401811292442956, 4.6725560653643994, 0.0,
    -3.3401811292442956, 4.6725560653643994, 0.0, -3.5220561640991237,
    4.9312580221152817, 0.0, -3.5220561640991237, 4.9312580221152817, 0.0,
    -3.7035490569484968, 5.1898696342119521, 0.0, -3.7035490569484968,
    5.1898696342119521, 0.0, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0,
    0.13919179616088076, 0.07490340326547, -0.0, 0.13919179616088076,
    0.07490340326547, -0.0, 0.27800912838762087, 0.15000797759616485, -0.0,
    0.27800912838762087, 0.15000797759616485, -0.0, 0.41645288308364015,
    0.22531337660988532, -0.0, 0.41645288308364015, 0.22531337660988532, -0.0,
    0.55452394451511977, 0.30081925480716887, -0.0, 0.55452394451511977,
    0.30081925480716887, -0.0, 0.69222319481614147, 0.37652526756918425, -0.0,
    0.69222319481614147, 0.37652526756918425, -0.0, 0.82955151399381344,
    0.4524310711556313, -0.0, 0.82955151399381344, 0.4524310711556313, -0.0,
    0.96650977993338516, 0.52853632270264572, -0.0, 0.96650977993338516,
    0.52853632270264572, -0.0, 1.1030988684033489, 0.60484068022070858, -0.0,
    1.1030988684033489, 0.60484068022070858, -0.0, 1.2393196530605293,
    0.68134380259256155, -0.0, 1.2393196530605293, 0.68134380259256155, -0.0,
    1.375173005455161, 0.75804534957112635, -0.0, 1.375173005455161,
    0.75804534957112635, -0.0, 1.5106597950359535, 0.83494498177742982, -0.0,
    1.5106597950359535, 0.83494498177742982, -0.0, 1.6457808891551446,
    0.91204236069853384, -0.0, 1.6457808891551446, 0.91204236069853384, -0.0,
    1.7805371530735408, 0.98933714868547, -0.0, 1.7805371530735408,
    0.98933714868547, -0.0, 1.914929449965546, 1.0668290089511794, -0.0,
    1.914929449965546, 1.0668290089511794, -0.0, 2.048958640924178,
    1.1445176055684576, -0.0, 2.048958640924178, 1.1445176055684576, -0.0,
    2.1826255849660736, 1.2224026034679043, -0.0, 2.1826255849660736,
    1.2224026034679043, -0.0, 2.315931139036481, 1.3004836684358778, -0.0,
    2.315931139036481, 1.3004836684358778, -0.0, 2.4488761580142389,
    1.3787604671124549, -0.0, 2.4488761580142389, 1.3787604671124549, -0.0,
    2.5814614947167467, 1.4572326669893951, -0.0, 2.5814614947167467,
    1.4572326669893951, -0.0, 2.7136879999049208, 1.53589993640811, -0.0,
    2.7136879999049208, 1.53589993640811, -0.0, -0.13919179616088076,
    -0.07490340326547, 0.0, -0.13919179616088076, -0.07490340326547, 0.0,
    -0.27800912838762087, -0.15000797759616485, 0.0, -0.27800912838762087,
    -0.15000797759616485, 0.0, -0.41645288308364015, -0.22531337660988532, 0.0,
    -0.41645288308364015, -0.22531337660988532, 0.0, -0.55452394451511977,
    -0.30081925480716887, 0.0, -0.55452394451511977, -0.30081925480716887, 0.0,
    -0.69222319481614147, -0.37652526756918425, 0.0, -0.69222319481614147,
    -0.37652526756918425, 0.0, -0.82955151399381344, -0.4524310711556313, 0.0,
    -0.82955151399381344, -0.4524310711556313, 0.0, -0.96650977993338516,
    -0.52853632270264572, 0.0, -0.96650977993338516, -0.52853632270264572, 0.0,
    -1.1030988684033489, -0.60484068022070858, 0.0, -1.1030988684033489,
    -0.60484068022070858, 0.0, -1.2393196530605293, -0.68134380259256155, 0.0,
    -1.2393196530605293, -0.68134380259256155, 0.0, -1.375173005455161,
    -0.75804534957112635, 0.0, -1.375173005455161, -0.75804534957112635, 0.0,
    -1.5106597950359535, -0.83494498177742982, 0.0, -1.5106597950359535,
    -0.83494498177742982, 0.0, -1.6457808891551446, -0.91204236069853384, 0.0,
    -1.6457808891551446, -0.91204236069853384, 0.0, -1.7805371530735408,
    -0.98933714868547, 0.0, -1.7805371530735408, -0.98933714868547, 0.0,
    -1.914929449965546, -1.0668290089511794, 0.0, -1.914929449965546,
    -1.0668290089511794, 0.0, -2.048958640924178, -1.1445176055684576, 0.0,
    -2.048958640924178, -1.1445176055684576, 0.0, -2.1826255849660736,
    -1.2224026034679043, 0.0, -2.1826255849660736, -1.2224026034679043, 0.0,
    -2.315931139036481, -1.3004836684358778, 0.0, -2.315931139036481,
    -1.3004836684358778, 0.0, -2.4488761580142389, -1.3787604671124549, 0.0,
    -2.4488761580142389, -1.3787604671124549, 0.0, -2.5814614947167467,
    -1.4572326669893951, 0.0, -2.5814614947167467, -1.4572326669893951, 0.0,
    -2.7136879999049208, -1.53589993640811, 0.0, -2.7136879999049208,
    -1.53589993640811, 0.0, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T l[360]{ 0.3579347593572087, 0.0, 0.0, 0.3579347593572087,
    0.0, 0.0, 0.71539677240113031, 0.0, 0.0, 0.71539677240113031, 0.0, 0.0,
    1.0723866635167145, 0.0, 0.0, 1.0723866635167145, 0.0, 0.0,
    1.4289050562642476, 0.0, 0.0, 1.4289050562642476, 0.0, 0.0,
    1.7849525733804419, 0.0, 0.0, 1.7849525733804419, 0.0, 0.0,
    2.1405298367795229, 0.0, 0.0, 2.1405298367795229, 0.0, 0.0,
    2.4956374675543165, 0.0, 0.0, 2.4956374675543165, 0.0, 0.0,
    2.8502760859773328, 0.0, 0.0, 2.8502760859773328, 0.0, 0.0,
    3.2044463115018504, 0.0, 0.0, 3.2044463115018504, 0.0, 0.0,
    3.5581487627629982, 0.0, 0.0, 3.5581487627629982, 0.0, 0.0,
    3.9113840575788359, 0.0, 0.0, 3.9113840575788359, 0.0, 0.0,
    4.2641528129514326, 0.0, 0.0, 4.2641528129514326, 0.0, 0.0,
    4.616455645067945, 0.0, 0.0, 4.616455645067945, 0.0, 0.0, 4.9682931693016954,
    0.0, 0.0, 4.9682931693016954, 0.0, 0.0, 5.3196660002132425, 0.0, 0.0,
    5.3196660002132425, 0.0, 0.0, 5.6705747515514595, 0.0, 0.0,
    5.6705747515514595, 0.0, 0.0, 6.021020036254602, 0.0, 0.0, 6.021020036254602,
    0.0, 0.0, 6.3710024664513822, 0.0, 0.0, 6.3710024664513822, 0.0, 0.0,
    6.7205226534620355, 0.0, 0.0, 6.7205226534620355, 0.0, 0.0, 7.06958120779939,
    0.0, 0.0, 7.06958120779939, 0.0, 0.0, -0.042560663779893933, 0.0, 0.0,
    -0.042560663779893933, 0.0, 0.0, -0.0850651150898698, 0.0, 0.0,
    -0.0850651150898698, 0.0, 0.0, -0.12751342817317227, 0.0, 0.0,
    -0.12751342817317227, 0.0, 0.0, -0.16990567717498842, 0.0, 0.0,
    -0.16990567717498842, 0.0, 0.0, -0.21224193614257725, 0.0, 0.0,
    -0.21224193614257725, 0.0, 0.0, -0.254522279025399, 0.0, 0.0,
    -0.254522279025399, 0.0, 0.0, -0.29674677967524432, 0.0, 0.0,
    -0.29674677967524432, 0.0, 0.0, -0.33891551184636343, 0.0, 0.0,
    -0.33891551184636343, 0.0, 0.0, -0.38102854919559465, 0.0, 0.0,
    -0.38102854919559465, 0.0, 0.0, -0.4230859652824932, 0.0, 0.0,
    -0.4230859652824932, 0.0, 0.0, -0.46508783356945982, 0.0, 0.0,
    -0.46508783356945982, 0.0, 0.0, -0.50703422742186888, 0.0, 0.0,
    -0.50703422742186888, 0.0, 0.0, -0.54892522010819667, 0.0, 0.0,
    -0.54892522010819667, 0.0, 0.0, -0.59076088480014921, 0.0, 0.0,
    -0.59076088480014921, 0.0, 0.0, -0.6325412945727904, 0.0, 0.0,
    -0.6325412945727904, 0.0, 0.0, -0.67426652240466922, 0.0, 0.0,
    -0.67426652240466922, 0.0, 0.0, -0.71593664117794742, 0.0, 0.0,
    -0.71593664117794742, 0.0, 0.0, -0.757551723678527, 0.0, 0.0,
    -0.757551723678527, 0.0, 0.0, -0.79911184259617707, 0.0, 0.0,
    -0.79911184259617707, 0.0, 0.0, -0.840617070524661, 0.0, 0.0,
    -0.840617070524661, 0.0, 0.0, -0.0099627081344734868, 0.0, 0.0,
    -0.0099627081344734868, 0.0, 0.0, -0.01991225791140331, 0.0, 0.0,
    -0.01991225791140331, 0.0, 0.0, -0.029848666709836498, 0.0, 0.0,
    -0.029848666709836498, 0.0, 0.0, -0.039771951885866519, 0.0, 0.0,
    -0.039771951885866519, 0.0, 0.0, -0.049682130772663577, 0.0, 0.0,
    -0.049682130772663577, 0.0, 0.0, -0.059579220680504912, 0.0, 0.0,
    -0.059579220680504912, 0.0, 0.0, -0.069463238896805016, 0.0, 0.0,
    -0.069463238896805016, 0.0, 0.0, -0.079334202686145827, 0.0, 0.0,
    -0.079334202686145827, 0.0, 0.0, -0.089192129290306912, 0.0, 0.0,
    -0.089192129290306912, 0.0, 0.0, -0.099037035928295547, 0.0, 0.0,
    -0.099037035928295547, 0.0, 0.0, -0.10886893979637684, 0.0, 0.0,
    -0.10886893979637684, 0.0, 0.0, -0.11868785806810371, 0.0, 0.0,
    -0.11868785806810371, 0.0, 0.0, -0.12849380789434692, 0.0, 0.0,
    -0.12849380789434692, 0.0, 0.0, -0.13828680640332505, 0.0, 0.0,
    -0.13828680640332505, 0.0, 0.0, -0.1480668707006344, 0.0, 0.0,
    -0.1480668707006344, 0.0, 0.0, -0.15783401786927884, 0.0, 0.0,
    -0.15783401786927884, 0.0, 0.0, -0.16758826496969964, 0.0, 0.0,
    -0.16758826496969964, 0.0, 0.0, -0.17732962903980537, 0.0, 0.0,
    -0.17732962903980537, 0.0, 0.0, -0.18705812709500158, 0.0, 0.0,
    -0.18705812709500158, 0.0, 0.0, -0.1967737761282205, 0.0, 0.0,
    -0.1967737761282205, 0.0, 0.0 };

  static const real_T l_0[360]{ 0.0, -0.1888968412469153, -0.060492001339316988,
    0.0, -0.1888968412469153, -0.060492001339316988, 0.0, -0.3773313347058973,
    -0.1209347881183509, 0.0, -0.3773313347058973, -0.1209347881183509, 0.0,
    -0.56530446924504474, -0.18132836030302127, 0.0, -0.56530446924504474,
    -0.18132836030302127, 0.0, -0.75281723176956472, -0.24167271798767209, 0.0,
    -0.75281723176956472, -0.24167271798767209, 0.0, -0.939870607225484,
    -0.30196786139466048, 0.0, -0.939870607225484, -0.30196786139466048, 0.0,
    -1.1264655786033526, -0.36221379087394617, 0.0, -1.1264655786033526,
    -0.36221379087394617, 0.0, -1.3126031269419416, -0.42241050690268223, 0.0,
    -1.3126031269419416, -0.42241050690268223, 0.0, -1.4982842313319336,
    -0.48255801008480659, 0.0, -1.4982842313319336, -0.48255801008480659, 0.0,
    -1.6835098689196064, -0.54265630115063446, 0.0, -1.6835098689196064,
    -0.54265630115063446, 0.0, -1.8682810149105107, -0.6027053809564521, 0.0,
    -1.8682810149105107, -0.6027053809564521, 0.0, -2.05259864257314,
    -0.66270525048411089, 0.0, -2.05259864257314, -0.66270525048411089, 0.0,
    -2.2364637232425952, -0.72265591084062319, 0.0, -2.2364637232425952,
    -0.72265591084062319, 0.0, -2.41987722632424, -0.78255736325775849, 0.0,
    -2.41987722632424, -0.78255736325775849, 0.0, -2.6028401192973529,
    -0.842409609091641, 0.0, -2.6028401192973529, -0.842409609091641, 0.0,
    -2.7853533677187707, -0.90221264982234783, 0.0, -2.7853533677187707,
    -0.90221264982234783, 0.0, -2.9674179352265249, -0.96196648705350851, 0.0,
    -2.9674179352265249, -0.96196648705350851, 0.0, -3.1490347835434722,
    -1.0216711225119051, 0.0, -3.1490347835434722, -1.0216711225119051, 0.0,
    -3.3302048724809192, -1.081326558047073, 0.0, -3.3302048724809192,
    -1.081326558047073, 0.0, -3.5109291599422385, -1.1409327956309046, 0.0,
    -3.5109291599422385, -1.1409327956309046, 0.0, -3.69120860192648,
    -1.2004898373572503, 0.0, -3.69120860192648, -1.2004898373572503, 0.0,
    0.1599184353223182, -0.089801925435736424, 0.0, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.31964914489275253, -0.17947362197539324, 0.0,
    0.31964914489275253, -0.17947362197539324, 0.0, 0.47919231323908462,
    -0.26901527273781028, 0.0, 0.47919231323908462, -0.26901527273781028, 0.0,
    0.63854812478731549, -0.35842706059465185, 0.0, 0.63854812478731549,
    -0.35842706059465185, 0.0, 0.797716763861511, -0.44770916817072121, 0.0,
    0.797716763861511, -0.44770916817072121, 0.0, 0.956698414683647,
    -0.53686177784427469, 0.0, 0.956698414683647, -0.53686177784427469, 0.0,
    1.1154932613734565, -0.62588507174733543, 0.0, 1.1154932613734565,
    -0.62588507174733543, 0.0, 1.2741014879482764, -0.71477923176600711, 0.0,
    1.2741014879482764, -0.71477923176600711, 0.0, 1.4325232783228956,
    -0.80354443954078658, 0.0, 1.4325232783228956, -0.80354443954078658, 0.0,
    1.5907588163094044, -0.89218087646687683, 0.0, 1.5907588163094044,
    -0.89218087646687683, 0.0, 1.748808285617043, -0.98068872369449933, 0.0,
    1.748808285617043, -0.98068872369449933, 0.0, 1.9066718698520528,
    -1.069068162129206, 0.0, 1.9066718698520528, -1.069068162129206, 0.0,
    2.0643497525175265, -1.1573193724321906, 0.0, 2.0643497525175265,
    -1.1573193724321906, 0.0, 2.221842117013261, -1.2454425350206004, 0.0,
    2.221842117013261, -1.2454425350206004, 0.0, 2.3791491466356085,
    -1.3334378300678469, 0.0, 2.3791491466356085, -1.3334378300678469, 0.0,
    2.536271024577331, -1.421305437503916, 0.0, 2.536271024577331,
    -1.421305437503916, 0.0, 2.6932079339274537, -1.5090455370156788, 0.0,
    2.6932079339274537, -1.5090455370156788, 0.0, 2.8499600576711197,
    -1.596658308047201, 0.0, 2.8499600576711197, -1.596658308047201, 0.0,
    3.0065275786894463, -1.6841439298000525, 0.0, 3.0065275786894463,
    -1.6841439298000525, 0.0, 3.16291067975938, -1.7715025812336165, 0.0,
    3.16291067975938, -1.7715025812336165, 0.0, -0.10536597514810161,
    0.1119360731420971, 0.0, -0.10536597514810161, 0.1119360731420971, 0.0,
    -0.21068448591173913, 0.22372201869524067, 0.0, -0.21068448591173913,
    0.22372201869524067, 0.0, -0.31595540965099317, 0.3353580258245521, 0.0,
    -0.31595540965099317, 0.3353580258245521, 0.0, -0.4211786242426373,
    0.446844283480486, 0.0, -0.4211786242426373, 0.446844283480486, 0.0,
    -0.52635400807880228, 0.55818098039902453, 0.0, -0.52635400807880228,
    0.55818098039902453, 0.0, -0.63148144006564322, 0.6693683051018724, 0.0,
    -0.63148144006564322, 0.6693683051018724, 0.0, -0.73656079962200982,
    0.78040644589665065, 0.0, -0.73656079962200982, 0.78040644589665065, 0.0,
    -0.84159196667811931, 0.89129559087709187, 0.0, -0.84159196667811931,
    0.89129559087709187, 0.0, -0.94657482167423235, 1.0020359279232338, 0.0,
    -0.94657482167423235, 1.0020359279232338, 0.0, -1.051509245559332,
    1.1126276447016141, 0.0, -1.051509245559332, 1.1126276447016141, 0.0,
    -1.1563951197898057, 1.2230709286654644, 0.0, -1.1563951197898057,
    1.2230709286654644, 0.0, -1.2612323263281291, 1.3333659670549041, 0.0,
    -1.2612323263281291, 1.3333659670549041, 0.0, -1.3660207476415547,
    1.4435129468971351, 0.0, -1.3660207476415547, 1.4435129468971351, 0.0,
    -1.4707602667008017, 1.5535120550066355, 0.0, -1.4707602667008017,
    1.5535120550066355, 0.0, -1.5754507669787494, 1.6633634779853532, 0.0,
    -1.5754507669787494, 1.6633634779853532, 0.0, -1.6800921324491334,
    1.7730674022229005, 0.0, -1.6800921324491334, 1.7730674022229005, 0.0,
    -1.7846842475852445, 1.8826240138967472, 0.0, -1.7846842475852445,
    1.8826240138967472, 0.0, -1.8892269973586309, 1.9920334989724147, 0.0,
    -1.8892269973586309, 1.9920334989724147, 0.0, -1.9937202672378032,
    2.1012960432036696, 0.0, -1.9937202672378032, 2.1012960432036696, 0.0,
    -2.098163943186941, 2.2104118321327175, 0.0, -2.098163943186941,
    2.2104118321327175 };

  static const real_T l_1[360]{ 0.12647183082737606, -0.14396327290074137, 0.0,
    0.12647183082737606, -0.14396327290074137, 0.0, 0.25266713091650589,
    -0.28784633762178591, 0.0, 0.25266713091650589, -0.28784633762178591, 0.0,
    0.37858657539014817, -0.43164948302538075, 0.0, 0.37858657539014817,
    -0.43164948302538075, 0.0, 0.50423083775077449, -0.57537299731355918, 0.0,
    0.50423083775077449, -0.57537299731355918, 0.0, 0.62960058988446832,
    -0.71901716802974069, 0.0, 0.62960058988446832, -0.71901716802974069, 0.0,
    0.75469650206481409, -0.862582282060328, 0.0, 0.75469650206481409,
    -0.862582282060328, 0.0, 0.87951924295677675, -1.0060686256362998, 0.0,
    0.87951924295677675, -1.0060686256362998, 0.0, 1.0040694796205727,
    -1.1494764843347998, 0.0, 1.0040694796205727, -1.1494764843347998, 0.0,
    1.1283478775155311, -1.2928061430807216, 0.0, 1.1283478775155311,
    -1.2928061430807216, 0.0, 1.2523551005039453, -1.4360578861482902, 0.0,
    1.2523551005039453, -1.4360578861482902, 0.0, 1.3760918108549158,
    -1.57923199716264, 0.0, 1.3760918108549158, -1.57923199716264, 0.0,
    1.4995586692481842, -1.7223287591013881, 0.0, 1.4995586692481842,
    -1.7223287591013881, 0.0, 1.6227563347779563, -1.8653484542962038, 0.0,
    1.6227563347779563, -1.8653484542962038, 0.0, 1.7456854649567182,
    -2.0082913644343763, 0.0, 1.7456854649567182, -2.0082913644343763, 0.0,
    1.8683467157190412, -2.1511577705603755, 0.0, 1.8683467157190412,
    -2.1511577705603755, 0.0, 1.9907407414253795, -2.2939479530774123, 0.0,
    1.9907407414253795, -2.2939479530774123, 0.0, 2.112868194865857,
    -2.4366621917489919, 0.0, 2.112868194865857, -2.4366621917489919, 0.0,
    2.2347297272640461, -2.5793007657004661, 0.0, 2.2347297272640461,
    -2.5793007657004661, 0.0, 2.3563259882807368, -2.7218639534205811, 0.0,
    2.3563259882807368, -2.7218639534205811, 0.0, 2.477657626017697,
    -2.86435203276302, 0.0, 2.477657626017697, -2.86435203276302, 0.0,
    -0.18891678043810189, 0.2603991729373682, 0.0, -0.18891678043810189,
    0.2603991729373682, 0.0, -0.37743412338800142, 0.52070047841969547, 0.0,
    -0.37743412338800142, 0.52070047841969547, 0.0, -0.56555300944331421,
    0.78090434240065543, 0.0, -0.56555300944331421, 0.78090434240065543, 0.0,
    -0.7532744168461627, 1.0410111898780745, 0.0, -0.7532744168461627,
    1.0410111898780745, 0.0, -0.94059932149283521, 1.3010214448962558, 0.0,
    -0.94059932149283521, 1.3010214448962558, 0.0, -1.1275286969394305,
    1.560935530548297, 0.0, -1.1275286969394305, 1.560935530548297, 0.0,
    -1.3140635144074895, 1.8207538689784042, 0.0, -1.3140635144074895,
    1.8207538689784042, 0.0, -1.5002047427896132, 2.0804768813841976, 0.0,
    -1.5002047427896132, 2.0804768813841976, 0.0, -1.6859533486550675,
    2.3401049880190148, 0.0, -1.6859533486550675, 2.3401049880190148, 0.0,
    -1.8713102962553732, 2.5996386081942049, 0.0, -1.8713102962553732,
    2.5996386081942049, 0.0, -2.0562765475298841, 2.8590781602814213, 0.0,
    -2.0562765475298841, 2.8590781602814213, 0.0, -2.2408530621113516,
    3.1184240617149053, 0.0, -2.2408530621113516, 3.1184240617149053, 0.0,
    -2.4250407973314743, 3.3776767289937664, 0.0, -2.4250407973314743,
    3.3776767289937664, 0.0, -2.6088407082264355, 3.6368365776842562, 0.0,
    -2.6088407082264355, 3.6368365776842562, 0.0, -2.7922537475424285,
    3.8959040224220369, 0.0, -2.7922537475424285, 3.8959040224220369, 0.0,
    -2.975280865741166, 4.154879476914445, 0.0, -2.975280865741166,
    4.154879476914445, 0.0, -3.1579230110053778, 4.4137633539427483, 0.0,
    -3.1579230110053778, 4.4137633539427483, 0.0, -3.3401811292442956,
    4.6725560653643994, 0.0, -3.3401811292442956, 4.6725560653643994, 0.0,
    -3.5220561640991237, 4.9312580221152817, 0.0, -3.5220561640991237,
    4.9312580221152817, 0.0, -3.7035490569484968, 5.1898696342119521, 0.0,
    -3.7035490569484968, 5.1898696342119521, 0.0, -0.13919179616088076,
    -0.07490340326547, 0.0, -0.13919179616088076, -0.07490340326547, 0.0,
    -0.27800912838762087, -0.15000797759616485, 0.0, -0.27800912838762087,
    -0.15000797759616485, 0.0, -0.41645288308364015, -0.22531337660988532, 0.0,
    -0.41645288308364015, -0.22531337660988532, 0.0, -0.55452394451511977,
    -0.30081925480716887, 0.0, -0.55452394451511977, -0.30081925480716887, 0.0,
    -0.69222319481614147, -0.37652526756918425, 0.0, -0.69222319481614147,
    -0.37652526756918425, 0.0, -0.82955151399381344, -0.4524310711556313, 0.0,
    -0.82955151399381344, -0.4524310711556313, 0.0, -0.96650977993338516,
    -0.52853632270264572, 0.0, -0.96650977993338516, -0.52853632270264572, 0.0,
    -1.1030988684033489, -0.60484068022070858, 0.0, -1.1030988684033489,
    -0.60484068022070858, 0.0, -1.2393196530605293, -0.68134380259256155, 0.0,
    -1.2393196530605293, -0.68134380259256155, 0.0, -1.375173005455161,
    -0.75804534957112635, 0.0, -1.375173005455161, -0.75804534957112635, 0.0,
    -1.5106597950359535, -0.83494498177742982, 0.0, -1.5106597950359535,
    -0.83494498177742982, 0.0, -1.6457808891551446, -0.91204236069853384, 0.0,
    -1.6457808891551446, -0.91204236069853384, 0.0, -1.7805371530735408,
    -0.98933714868547, 0.0, -1.7805371530735408, -0.98933714868547, 0.0,
    -1.914929449965546, -1.0668290089511794, 0.0, -1.914929449965546,
    -1.0668290089511794, 0.0, -2.048958640924178, -1.1445176055684576, 0.0,
    -2.048958640924178, -1.1445176055684576, 0.0, -2.1826255849660736,
    -1.2224026034679043, 0.0, -2.1826255849660736, -1.2224026034679043, 0.0,
    -2.315931139036481, -1.3004836684358778, 0.0, -2.315931139036481,
    -1.3004836684358778, 0.0, -2.4488761580142389, -1.3787604671124549, 0.0,
    -2.4488761580142389, -1.3787604671124549, 0.0, -2.5814614947167467,
    -1.4572326669893951, 0.0, -2.5814614947167467, -1.4572326669893951, 0.0,
    -2.7136879999049208, -1.53589993640811, 0.0, -2.7136879999049208,
    -1.53589993640811, 0.0 };

  static const real_T c[246]{ 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0, 612.0,
    612.0, 612.0, 80.0, 80.0, 80.0, -0.0, -0.0, -0.0 };

  static const real_T h[180]{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T o[180]{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 1.0 };

  static const real_T b_A_0[100]{ 0.99801496201563555, 0.00012982817451900192,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.001444494825713849,
    0.99878101747798131, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99801496201563555, 0.00012982817451900192, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.001444494825713849, 0.99878101747798131, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0 };

  static const real_T b_A_1[100]{ 0.99747143658641013, 0.0011849245467832059,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.00030050237926615431,
    1.0004838136541614, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.99747143658641013, 0.0011849245467832059, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.00030050237926615431, 1.0004838136541614, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0 };

  static const real_T b_A[64]{ 0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.99867923888102927, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T c_a_0[60]{ -1.3493010369352151E-16,
    -1.1675558559895421E-15, 2.2450629369655467E-16, -9.2500152910841517E-16,
    0.024689453048712479, -7.943777873772287E-16, 2.8823938352524109E-16,
    3.5570950392137461E-16, 8.7943241351772457E-16, 1.9285310195954591E-15,
    0.14580276217308596, -0.021598327103105264, 0.14580276217309096,
    -0.021598327103112134, -1.0054774783728353E-15, 0.02222088852315177,
    0.0031720277995458445, 1.5008670922687785E-16, -0.0024685645256298031,
    0.003172027799553423, -0.015523318927034492, 0.085203655514006313,
    -0.015523318927008747, 0.0852036555139457, -7.0589066818986132E-16,
    -0.0029029803765247804, 0.023198516335863745, -3.1464216881124013E-16,
    -0.0029029803765413028, -0.0014909367129008872, 1.7463654942675622E-15,
    5.2715088122423186E-17, 2.3678313046590209E-16, 2.760023978062948E-15,
    3.5851626435326736E-16, -1.5415647015806093E-15, -3.5518690565450184E-16,
    0.02468945304871167, -4.6913613292291945E-16, -2.8471103718024361E-15,
    0.14580276217267166, -0.021598327102856273, 0.14580276217266042,
    -0.021598327102846406, 1.2378015318971689E-15, -0.0024685645252104874,
    0.003172027799306045, -1.3815255191351443E-16, 0.022220888523571687,
    0.0031720277992777178, -0.015523318926905885, 0.085203655513833534,
    -0.015523318926929857, 0.085203655513889545, 8.926914482010702E-16,
    -0.0029029803766526282, -0.0014909367127887497, -5.4423394757560848E-17,
    -0.0029029803766385991, 0.023198516335942394 };

  static const real_T c_a_1[60]{ 0.13063875476843981, -0.053690379918034888,
    0.13010284540699737, -0.048537966075219825, 0.022672305893804586,
    0.0024236170740209533, 3.10900862053873E-15, -0.0014887611813181533,
    -0.0026564643729729276, 7.43199707040387E-16, -0.055206463682881686,
    0.19733828508082932, -0.050054049837865343, 0.14780123586116181,
    0.0039397008388678181, -0.022251681085583872, 5.6867795365447866E-17,
    -0.0011403806103201641, 0.0019004869693013892, -4.7316906414455058E-16,
    3.6351433545115045E-15, 1.4040431858075756E-14, -3.2789053806757555E-15,
    1.4046092584410489E-14, -4.72902566177774E-16, -1.3774131006650049E-14,
    0.024689453048711677, 1.2811698620675945E-16, -1.4277844130261156E-14,
    -3.0407729517469778E-16, 0.13010284540701042, -0.048537966075151331,
    0.13063875476844858, -0.053690379917926627, -0.0014887611813312498,
    -0.0026564643730341941, -3.1420550328853386E-15, 0.022672305893797325,
    0.0024236170739114441, -8.5716857860795223E-16, -0.050054049838186253,
    0.14780123586192084, -0.055206463683200487, 0.19733828508159162,
    -0.0011403806100064946, 0.0019004869685422226, -1.3321387737354083E-16,
    0.0039397008391853314, -0.022251681086325373, 3.2521052678624355E-16,
    7.8573329095329507E-16, 1.7729695040513638E-16, -3.1569984595464808E-15,
    -1.6399356126597423E-15, 1.3442269461789139E-16, -4.8148531004784216E-16,
    -3.0505811159304027E-16, 2.1140431622661049E-15, 1.6917293461581915E-15,
    0.024689453048712246 };

  static const real_T d_a_0[60]{ -1.3286135634905549E-15,
    -5.2719344895617172E-16, -8.875729112587828E-16, -6.3371986351093964E-16,
    0.024689453048712038, -2.8609304354609092E-16, 3.487319325503805E-16,
    9.8894977402454656E-17, 9.55001448329591E-16, 7.161099874107621E-16,
    0.14554453682362364, -0.021553069813391, 0.14554453682362928,
    -0.021553069813398929, -1.0180638245063373E-15, 0.022220888523201865,
    0.0031720277995270292, 1.3458197354116882E-16, -0.0024685645255787584,
    0.0031720277995358984, -0.015615580788843781, 0.085097778382963646,
    -0.015615580788818249, 0.085097778382903513, -7.03315603900164E-16,
    -0.0029029803765225669, 0.023198516335863107, -3.1093226788899258E-16,
    -0.0029029803765393014, -0.0014909367129007042, 1.4809047485715802E-15,
    5.6341244474062862E-16, 1.1072821416542014E-15, 3.1356200306111844E-15,
    1.0442222682848953E-16, 1.5578150588167511E-15, -2.1249548531906002E-16,
    0.024689453048711903, -4.9179393494252262E-15, -7.4317177983426109E-16,
    0.14554453682321053, -0.0215530698131427, 0.14554453682319984,
    -0.02155306981313385, 1.2215716358558572E-15, -0.0024685645251605244,
    0.0031720277992875533, -1.5444858114192032E-16, 0.022220888523622011,
    0.003172027799260489, -0.01561558078871641, 0.085097778382792283,
    -0.015615580788740587, 0.085097778382848141, 8.9228929914245286E-16,
    -0.0029029803766484259, -0.0014909367127899484, -5.1559145272375194E-17,
    -0.0029029803766346973, 0.023198516335941916 };

  static const real_T d_a_1[60]{ 0.13032456047964627, -0.053561558989666934,
    0.12978845788638857, -0.048407287350832322, 0.02267230589380282,
    0.00242361707403653, 3.0965418352169834E-15, -0.0014887611813189075,
    -0.0026564643729636412, 7.4307922647746408E-16, -0.055126171262810331,
    0.19736834454368832, -0.04997189962178173, 0.14781343384486076,
    0.003939700838866628, -0.022251681085588802, 4.964896767976845E-17,
    -0.0011403806103175548, 0.0019004869692955454, -4.7569554365069832E-16,
    1.9227121156520315E-15, 1.2072502537598143E-14, -1.501053003416049E-15,
    8.23344142778545E-15, -1.0236402466968285E-15, -1.1446044762745138E-14,
    0.024689453048712146, 8.9626519721957567E-16, -8.1482188373998946E-15,
    4.9749100862038536E-17, 0.12978845788640161, -0.048407287350765014,
    0.13032456047965524, -0.053561558989558521, -0.0014887611813321185,
    -0.0026564643730245946, -3.1357416027813126E-15, 0.022672305893795569,
    0.0024236170739276846, -8.5660274826988225E-16, -0.049971899622101336,
    0.14781343384561868, -0.055126171263127528, 0.19736834454444707,
    -0.0011403806100045892, 0.001900486968539325, -1.3682064082528937E-16,
    0.0039397008391829306, -0.022251681086328183, 3.275458957324977E-16,
    -1.2628574529309357E-17, 1.2008069891760323E-15, -1.9618966754900136E-15,
    -2.3814872173872742E-15, -1.6953886056165856E-16, -1.3350979816451865E-15,
    -2.7270192393292958E-17, 1.7396574001077273E-15, 2.704285886690011E-15,
    0.024689453048712118 };

  static const real_T c_a[48]{ 0.18830679432649863, 0.18830679432650269,
    0.022059066214098054, 2.5230328845122769E-17, 8.9692687745952415E-17,
    -0.0026303868346078179, -1.2571160900117593E-16, -5.71467306939525E-17,
    -4.1204398628128764E-15, -2.816212382747728E-15, 4.1692819362639037E-15,
    0.02468945304871217, -2.1028121777614864E-17, 2.7309599930514198E-15,
    2.4125046479575912E-16, -1.0204125760738813E-16, -4.6450947826779657E-16,
    -2.761872039395948E-15, 6.2425094920799318E-16, -1.976925546372601E-17,
    0.024689453048712048, 2.9356060741096639E-15, 2.4332183317588693E-17,
    5.7062545372552586E-17, 0.18830679432521083, 0.18830679432520178,
    -0.0026303868333161242, 2.8773390387444934E-18, 2.0233643636645766E-16,
    0.022059066215401157, -4.2901040066494445E-18, 2.8716963918091621E-17,
    -2.0337178809473559E-15, -1.7540017816779046E-15, 2.02407570490301E-15,
    2.4141424086370343E-16, 2.3030997484589184E-17, 1.5571067788478613E-15,
    0.024689453048712475, -3.3566203160325031E-18, -1.136476839024345E-15,
    3.0764737487189331E-15, 1.0347365692837263E-15, -1.0058861543012247E-16,
    5.4361235871880759E-17, -3.0427365314654311E-15, -4.1801243508450244E-18,
    0.024689453048711962 };

  static const real_T d_a[48]{ 0.1880580860341125, 0.18805808603411617,
    0.022059066214095539, 4.3831064811932122E-17, 6.7337364339229245E-17,
    -0.0026303868346073834, -1.1320233898242198E-16, -4.3175671180392858E-17,
    -2.8947889241124266E-15, -1.9732620600075495E-15, 3.9540462567412853E-15,
    0.024689453048712139, -1.7582053102690953E-17, 1.2530704858616323E-15,
    1.2935134781352392E-16, 1.7434853494517857E-17, -1.1846528643462433E-15,
    -3.1504850347149165E-15, 5.1527402493887339E-16, -4.8432915240592263E-17,
    0.024689453048712073, 2.0666863622775897E-15, -1.0939020792469785E-16,
    3.0001112140270809E-17, 0.18805808603282764, 0.18805808603281798,
    -0.0026303868333193889, 1.9821674517527381E-17, 1.7797214899950445E-16,
    0.022059066215400495, 9.2515026921393948E-18, 4.29739330190316E-17,
    -1.4701369148339633E-16, 1.6096813127316034E-16, 1.2236226095001052E-15,
    -1.5956560300803347E-16, -3.1124234652486909E-16, 1.0258195202799193E-15,
    0.024689453048712163, 1.595109096193717E-17, -5.4774339592263138E-16,
    2.0496139574284967E-15, -1.0942770620180964E-16, 5.2962691312108752E-18,
    1.6388199070448957E-17, -1.134213914809288E-15, -2.63722928475333E-17,
    0.024689453048712107 };

  static const real_T a_7[30]{ -0.1888968412469153, -0.060492001339316988,
    -0.1888968412469153, -0.060492001339316988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1599184353223182, -0.089801925435736424, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.10536597514810161,
    0.1119360731420971, -0.10536597514810161, 0.1119360731420971, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T a_8[30]{ 0.12647183082737606, -0.14396327290074137,
    0.12647183082737606, -0.14396327290074137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.18891678043810189, 0.2603991729373682, -0.18891678043810189,
    0.2603991729373682, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.13919179616088076,
    -0.07490340326547, -0.13919179616088076, -0.07490340326547, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T a_a[30]{ -0.1888968412469153, -0.060492001339316988,
    -0.1888968412469153, -0.060492001339316988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1599184353223182, -0.089801925435736424, 0.1599184353223182,
    -0.089801925435736424, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.10536597514810161,
    0.1119360731420971, -0.10536597514810161, 0.1119360731420971, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T a_b[30]{ 0.12647183082737606, -0.14396327290074137,
    0.12647183082737606, -0.14396327290074137, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.18891678043810189, 0.2603991729373682, -0.18891678043810189,
    0.2603991729373682, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.13919179616088076,
    -0.07490340326547, -0.13919179616088076, -0.07490340326547, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T a_6[24]{ 0.3579347593572087, 0.3579347593572087, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.042560663779893933, -0.042560663779893933, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.0099627081344734868, -0.0099627081344734868, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0 };

  static const real_T a_9[24]{ 0.3579347593572087, 0.3579347593572087, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.042560663779893933, -0.042560663779893933, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.0099627081344734868, -0.0099627081344734868, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0 };

  static const real_T f[16]{ 0.034121465297356074, 0.0, 0.0, 0.0, 0.0,
    0.034121465297356074, 0.0, 0.0, 0.0, 0.0, 0.034121465297356074, 0.0, 0.0,
    0.0, 0.0, 100000.0 };

  static const int32_T q[246]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
    54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
    73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91,
    92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
    109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123,
    124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138,
    139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
    154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
    169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183,
    184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198,
    199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213,
    214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228,
    229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243,
    301, 302, 303 };

  static const int8_T b_a_0[60]{ 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T b_a_1[60]{ 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  static const int8_T b_a[48]{ 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1 };

  static const int8_T b[5]{ 0, 1, 2, 4, 5 };

  real_T Product1_j[144];
  real_T rseq[120];
  real_T rtb_useq_f[63];
  real_T tmp[60];
  real_T vseq[21];
  real_T f_0[16];
  real_T Sum_h[12];
  real_T xk_0[10];
  real_T xk_2[10];
  real_T xk_3[10];
  real_T xk[8];
  real_T xk_1[8];
  real_T rtb_Gain2[6];
  real_T rtb_Gain3[6];
  real_T rtb_ywt[6];
  real_T tmp_1[6];
  real_T tmp_3[6];
  real_T y_innov[6];
  real_T y_innov_0[6];
  real_T y_innov_1[6];
  real_T Sum2_c[3];
  real_T ext_mv[3];
  real_T rtb_Delay2[3];
  real_T rtb_Gain1[3];
  real_T rtb_u_l[3];
  real_T tmp_0[3];
  real_T tmp_2[3];
  real_T tmp_4[3];
  real_T ext_mv_idx_2;
  real_T rtb_sig;
  real_T status;
  uint16_T waypt;
  int8_T P0_2_tmp[12];
  boolean_T tmp_6[246];
  boolean_T x[5];
  boolean_T b_x[4];
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
        for (int32_T k{0}; k < 6; k++) {
          // Inport: '<Root>/y'
          rtDW.traj[k] = rtU.y[k];
        }
      } else {
        // '<S1>:520:6' else
        //  hold last waypoint pos
        // '<S1>:520:7' traj(:,1) = traj(:, waypt);
        for (int32_T k{0}; k < 6; k++) {
          rtDW.traj[k] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6 +
            k];
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
      __m128d tmp_5;
      real_T Sum2_c_0;
      real_T ext_mv_idx_0;
      real_T ext_mv_idx_1;
      real_T rtb_Delay2_tmp;
      int32_T i;
      int32_T i_0;
      int32_T k;
      boolean_T c_y;
      boolean_T exitg1;
      boolean_T guard1{ false };

      boolean_T guard2{ false };

      boolean_T guard3{ false };

      boolean_T guard4{ false };

      boolean_T guard5{ false };

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
            for (k = 0; k < 6; k++) {
              // Inport: '<Root>/y'
              rtDW.traj[k] = rtU.y[k];
            }
          } else {
            // '<S1>:520:6' else
            //  hold last waypoint pos
            // '<S1>:520:7' traj(:,1) = traj(:, waypt);
            for (k = 0; k < 6; k++) {
              rtDW.traj[k] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                6 + k];
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
          for (k = 0; k < 12; k++) {
            P0_2_tmp[k] = static_cast<int8_T>(k + 1);
          }

          // '<S1>:59:5' theta0_2 = theta(1:np*no);
          for (i = 0; i < 12; i++) {
            for (k = 0; k < 12; k++) {
              rtDW.P0_2[k + 12 * i] = rtY.P_e[((P0_2_tmp[i] - 1) * 24 +
                P0_2_tmp[k]) - 1];
            }

            rtDW.theta0_2[i] = rtY.theta[i];
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
          for (k = 0; k < 12; k++) {
            P0_2_tmp[k] = static_cast<int8_T>(k + 13);
          }

          // '<S1>:59:11' theta0_1 = theta(np*no+1:2*np*no);
          for (i = 0; i < 12; i++) {
            for (k = 0; k < 12; k++) {
              rtDW.P0_1[k + 12 * i] = rtY.P_e[((P0_2_tmp[i] - 1) * 24 +
                P0_2_tmp[k]) - 1];
            }

            rtDW.theta0_1[i] = rtY.theta[i + 12];
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
        for (i = 0; i < 6; i++) {
          rtDW.enAdapt_[i] = false;
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
      i = 0;
      i_0 = 0;
      for (k = 0; k < 12; k++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S160>/Product1'

        (void)std::memcpy(&rtY.P_e[i], &Product1_j[i_0], 12U * sizeof(real_T));

        // Outport: '<Root>/theta'
        rtY.theta[k] = Sum_h[k];
        i += 24;
        i_0 += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S160>/Sum2'

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
      i = 0;
      i_0 = 0;
      for (k = 0; k < 12; k++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S164>/Product1'

        (void)std::memcpy(&rtY.P_e[i + 300], &Product1_j[i_0], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[k + 12] = Sum_h[k];
        i += 24;
        i_0 += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S164>/Sum2'

      rtY.prmErr[3] = Sum2_c[0];
      rtY.prmErr[4] = Sum2_c[1];
      rtY.prmErr[5] = Sum2_c[2];

      // [u, ywt, yhat, currTraj] = ...
      // ampc(traj(:,waypt), currEv.r, y, ymax, y0, x0, u0, umax, uwt,...
      // excitation, theta, thetaSgn, k_2);
      // '<S1>:59:31' [u, ywt, currTraj] = gmpc(traj(:,waypt), currEv.r, y, ymax, umax, uwt, k_2); 
      // Simulink Function 'gmpc': '<S1>:863'
      for (i = 0; i <= 4; i += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Gain: '<S3>/Gain2' incorporates:
        //   Inport: '<Root>/ymax'

        (void)_mm_storeu_pd(&rtb_Gain2[i], _mm_mul_pd(_mm_set1_pd
          (rtP.Gain2_Gain_b), _mm_loadu_pd(&rtU.ymax[i])));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      for (i = 0; i <= 4; i += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Gain: '<S3>/Gain3' incorporates:
        //   Inport: '<Root>/ymax'

        (void)_mm_storeu_pd(&rtb_Gain3[i], _mm_mul_pd(_mm_set1_pd
          (rtP.Gain3_Gain_d), _mm_loadu_pd(&rtU.ymax[i])));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      // Outputs for Function Call SubSystem: '<S1>/gmpc'
      // Update for Delay: '<S3>/Delay' incorporates:
      //   Delay: '<S3>/Delay1'
      //   Inport: '<Root>/k_2'
      //   MATLAB Function: '<S3>/MATLAB Function'
      //   Outport: '<Root>/currEv'

      MATLABFunction(rtY.currEv.r, rtU.k_2, rtDW.Delay_DSTATE, rtb_ywt, Sum2_c,
                     &rtP);

      // MATLAB Function: '<S3>/MATLAB Function1'
      // MATLAB Function 'SupervisoryController/gmpc/MATLAB Function1': '<S87>:1' 
      // '<S87>:1:2' sig = gainSchSig_(ywt);
      // 'gainSchSig_:3' if isempty(sigPrev)
      // 'gainSchSig_:7' if ( ywt(1) > 0.5 && all(ywt(2:end) < 0.5) )...
      // 'gainSchSig_:8'         || ( ywt(4) > 0.5 && all(ywt([1:3,5:end]) < 0.5) ) 
      guard1 = false;
      guard2 = false;
      guard3 = false;
      guard4 = false;
      guard5 = false;
      if (rtb_ywt[0] > 0.5) {
        for (k = 0; k < 5; k++) {
          x[k] = (rtb_ywt[k + 1] < 0.5);
        }

        rstP2 = true;
        k = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (k <= 4)) {
          if (!x[k]) {
            rstP2 = false;
            exitg1 = true;
          } else {
            k++;
          }
        }

        if (rstP2) {
          // 'gainSchSig_:9' sig = 1;
          rtb_sig = 1.0;
        } else {
          guard5 = true;
        }
      } else {
        guard5 = true;
      }

      if (guard5) {
        if (rtb_ywt[3] > 0.5) {
          for (k = 0; k < 5; k++) {
            x[k] = (rtb_ywt[b[k]] < 0.5);
          }

          rstP1 = true;
          k = 0;
          exitg1 = false;
          while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (k <= 4)) {
            if (!x[k]) {
              rstP1 = false;
              exitg1 = true;
            } else {
              k++;
            }
          }

          if (rstP1) {
            // 'gainSchSig_:9' sig = 1;
            rtb_sig = 1.0;
          } else {
            guard4 = true;
          }
        } else {
          guard4 = true;
        }
      }

      if (guard4) {
        if ((rtb_ywt[1] > 0.5) && (rtb_ywt[2] > 0.5)) {
          b_x[0] = (rtb_ywt[0] < 0.5);
          b_x[1] = (rtb_ywt[3] < 0.5);
          b_x[2] = (rtb_ywt[4] < 0.5);
          b_x[3] = (rtb_ywt[5] < 0.5);
          c_y = true;
          k = 0;
          exitg1 = false;
          while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (k <= 3)) {
            if (!b_x[k]) {
              c_y = false;
              exitg1 = true;
            } else {
              k++;
            }
          }

          if (c_y) {
            // 'gainSchSig_:10' elseif ( ywt(2) > 0.5 && ywt(3) > 0.5 && all(ywt([1,4:end]) < 0.5) )... 
            // 'gainSchSig_:11'         || ( ywt(5) > 0.5 && ywt(6) > 0.5 && all(ywt([1:4]) < 0.5) ) 
            // 'gainSchSig_:12' sig = 2;
            rtb_sig = 2.0;
          } else {
            guard3 = true;
          }
        } else {
          guard3 = true;
        }
      }

      if (guard3) {
        if ((rtb_ywt[4] > 0.5) && (rtb_ywt[5] > 0.5)) {
          b_x[0] = (rtb_ywt[0] < 0.5);
          b_x[1] = (rtb_ywt[1] < 0.5);
          b_x[2] = (rtb_ywt[2] < 0.5);
          b_x[3] = (rtb_ywt[3] < 0.5);
          rstP2 = true;
          k = 0;
          exitg1 = false;
          while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (k <= 3)) {
            if (!b_x[k]) {
              rstP2 = false;
              exitg1 = true;
            } else {
              k++;
            }
          }

          if (rstP2) {
            // 'gainSchSig_:10' elseif ( ywt(2) > 0.5 && ywt(3) > 0.5 && all(ywt([1,4:end]) < 0.5) )... 
            // 'gainSchSig_:11'         || ( ywt(5) > 0.5 && ywt(6) > 0.5 && all(ywt([1:4]) < 0.5) ) 
            // 'gainSchSig_:12' sig = 2;
            rtb_sig = 2.0;
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }
      }

      if (guard2) {
        if ((rtb_ywt[0] > 0.5) && (rtb_ywt[2] > 0.5)) {
          b_x[0] = (rtb_ywt[1] < 0.5);
          b_x[1] = (rtb_ywt[3] < 0.5);
          b_x[2] = (rtb_ywt[4] < 0.5);
          b_x[3] = (rtb_ywt[5] < 0.5);
          rstP2 = true;
          k = 0;
          exitg1 = false;
          while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (k <= 3)) {
            if (!b_x[k]) {
              rstP2 = false;
              exitg1 = true;
            } else {
              k++;
            }
          }

          if (rstP2) {
            // 'gainSchSig_:13' elseif ( ywt(1) > 0.5 && ywt(3) > 0.5 && all(ywt([2,4:end]) < 0.5) )... 
            // 'gainSchSig_:14'         || ( ywt(4) > 0.5 && ywt(6) > 0.5 && all(ywt([1:3,5]) < 0.5) ) 
            // 'gainSchSig_:15' sig = 3;
            rtb_sig = 3.0;
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }
      }

      if (guard1) {
        if (rtb_ywt[3] > 0.5) {
          if (rtb_ywt[5] > 0.5) {
            b_x[0] = (rtb_ywt[0] < 0.5);
            b_x[1] = (rtb_ywt[1] < 0.5);
            b_x[2] = (rtb_ywt[2] < 0.5);
            b_x[3] = (rtb_ywt[4] < 0.5);
            rstP2 = true;
            k = 0;
            exitg1 = false;
            while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>
                     (0U)) == false) && (k <= 3)) {
              if (!b_x[k]) {
                rstP2 = false;
                exitg1 = true;
              } else {
                k++;
              }
            }

            if (rstP2) {
              // 'gainSchSig_:13' elseif ( ywt(1) > 0.5 && ywt(3) > 0.5 && all(ywt([2,4:end]) < 0.5) )... 
              // 'gainSchSig_:14'         || ( ywt(4) > 0.5 && ywt(6) > 0.5 && all(ywt([1:3,5]) < 0.5) ) 
              // 'gainSchSig_:15' sig = 3;
              rtb_sig = 3.0;
            } else {
              // 'gainSchSig_:16' else
              // 'gainSchSig_:17' sig = sigPrev;
              rtb_sig = rtDW.sigPrev;
            }
          } else {
            // 'gainSchSig_:16' else
            // 'gainSchSig_:17' sig = sigPrev;
            rtb_sig = rtDW.sigPrev;
          }
        } else {
          // 'gainSchSig_:16' else
          // 'gainSchSig_:17' sig = sigPrev;
          rtb_sig = rtDW.sigPrev;
        }
      }

      // 'gainSchSig_:19' sigPrev = sig;
      rtDW.sigPrev = rtb_sig;

      // End of MATLAB Function: '<S3>/MATLAB Function1'
      // End of Outputs for SubSystem: '<S1>/gmpc'
      for (k = 0; k <= 4; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Gain: '<S3>/Gain'
        tmp_5 = _mm_loadu_pd(&rtb_ywt[k]);

        // Gain: '<S3>/Gain'
        (void)_mm_storeu_pd(&rtb_ywt[k], _mm_mul_pd(_mm_set1_pd(rtP.beta), tmp_5));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      // Outputs for Function Call SubSystem: '<S1>/gmpc'
      // Gain: '<S3>/Gain1' incorporates:
      //   Inport: '<Root>/uwt'

      rtb_Gain1[0] = rtP.beta * rtU.uwt[0];
      rtb_Gain1[1] = rtP.beta * rtU.uwt[1];
      rtb_Gain1[2] = rtP.beta * rtU.uwt[2];

      // Delay: '<S3>/Delay2'
      rtb_Delay2[0] = rtDW.Delay2_DSTATE[0];
      rtb_Delay2[1] = rtDW.Delay2_DSTATE[1];
      rtb_Delay2[2] = rtDW.Delay2_DSTATE[2];

      // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S113>:1'
      // '<S113>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
      // '<S113>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
      // '<S113>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
      // '<S113>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
      //  Inputs (in BlockDataType except iA)
      //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
      // '<S113>:1:24' xk = convertDataType(xk0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    old_u:      last mv (calculated by MPC)
      // '<S113>:1:26' old_u = convertDataType(old_u0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    ym:         current measured output (used only with built-in KF)
      // '<S113>:1:28' ym = convertDataType(ym0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    ref:        output reference
      // '<S113>:1:30' ref = convertDataType(ref0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    md:         measured disturbance
      // '<S113>:1:32' md = convertDataType(md0,isDouble);
      //    umin:       run-time MV bound
      // '<S113>:1:34' umin = convertDataType(umin0,isDouble);
      //    umax:       run-time MV bound
      // '<S113>:1:36' umax = convertDataType(umax0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    ymin:       run-time OV bound
      // '<S113>:1:38' ymin = convertDataType(ymin0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    ymax:       run-time OV bound
      // '<S113>:1:40' ymax = convertDataType(ymax0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    E:          run-time mixed constraints
      // '<S113>:1:42' E = convertDataType(E0,isDouble);
      //    F:          run-time mixed constraints
      // '<S113>:1:44' F = convertDataType(F0,isDouble);
      //    G:          run-time mixed constraints
      // '<S113>:1:46' G = convertDataType(G0,isDouble);
      //    S:          run-time mixed constraints
      // '<S113>:1:48' S = convertDataType(S0,isDouble);
      //    switch_in:  if it matches "enable_value", MPC is active in control
      // '<S113>:1:50' switch_in = int32(switch_in0);
      //    ext_mv:     external last mv (actual)
      // '<S113>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    MVtarget:   MV reference
      // '<S113>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
      //    ywt:        run-time OV weights
      // '<S113>:1:56' ywt = convertDataType(ywt0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    uwt:        run-time MV weights
      // '<S113>:1:58' uwt = convertDataType(uwt0,isDouble);
      // '<S113>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S113>:1:252' if isa(u,'double')
      // '<S113>:1:253' y = u;
      //    duwt:       run-time DMV weights
      // '<S113>:1:60' duwt = convertDataType(duwt0,isDouble);
      //    ewt:     run-time Slack weights
      // '<S113>:1:62' ewt = convertDataType(ewt0,isDouble);
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
      // '<S113>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
      // '<S113>:1:96' isAdaptive = false;
      // '<S113>:1:97' isLTV = false;
      // '<S113>:1:98' ZERO = zeros('like',ref);
      // '<S113>:1:99' ONE = ones('like',ref);
      // '<S113>:1:100' hasMD = nv>int32(1);
      //  Pre-allocate all the MEX block outputs for the simulation mode
      // '<S113>:1:105' if isSimulation
      //  Get reference and MD signals -- accounting for previewing
      // '<S113>:1:119' if isSimulation
      // '<S113>:1:126' else
      //  When doing code generation, use M code directly
      // '<S113>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
      for (k = 0; k < 21; k++) {
        // MATLAB Function: '<S112>/optimizer'
        vseq[k] = 1.0;
      }

      for (k = 0; k < 20; k++) {
        for (i = 0; i < 6; i++) {
          // MATLAB Function: '<S112>/optimizer'
          rseq[i + k * static_cast<int32_T>(ny)] = rtDW.traj
            [(static_cast<int32_T>(rtDW.waypt) - 1) * 6 + i];
        }
      }

      // MATLAB Function: '<S112>/optimizer' incorporates:
      //   Gain: '<S92>/ext.mv_scale'
      //   UnitDelay: '<S92>/last_mv'

      //  External MV override.
      //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
      // '<S113>:1:133' old_u = old_u - uoff;
      // '<S113>:1:134' if no_mv
      // '<S113>:1:136' else
      // '<S113>:1:137' ext_mv = ext_mv - uoff;
      ext_mv[0] = rtP.extmv_scale_Gain_g[0] * rtb_Delay2[0];
      ext_mv[1] = rtP.extmv_scale_Gain_g[1] * rtb_Delay2[1];
      ext_mv[2] = rtP.extmv_scale_Gain_g[2] * rtb_Delay2[2];

      //  Bias correction
      // '<S113>:1:138' delmv = ext_mv - old_u;
      // '<S113>:1:139' old_u = ext_mv;
      //  Obtain x[k|k]
      // '<S113>:1:143' xk = xk - xoff;
      //  Remove offset
      // '<S113>:1:144' if CustomEstimation
      // '<S113>:1:147' else
      //  Default state estimation.
      //  Scale measured output and remove offset.
      // '<S113>:1:150' ym = ym.*RYscale(myindex) - myoff;
      //  Correct x(k|k-1) for possible external mv override.
      //  NOTE:  Offset was removed from x[k|k-1] at k=0.
      // '<S113>:1:153' xk = xk + Bu*delmv;
      ext_mv_idx_0 = ext_mv[0] - rtDW.last_mv_DSTATE[0];
      ext_mv_idx_1 = ext_mv[1] - rtDW.last_mv_DSTATE[1];
      ext_mv_idx_2 = ext_mv[2] - rtDW.last_mv_DSTATE[2];

      // End of Outputs for SubSystem: '<S1>/gmpc'
      for (k = 0; k <= 6; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Memory: '<S92>/last_x' incorporates:
        //   MATLAB Function: '<S112>/optimizer'

        tmp_5 = _mm_loadu_pd(&rtDW.last_x_PreviousInput[k]);

        // MATLAB Function: '<S112>/optimizer'
        (void)_mm_storeu_pd(&xk[k], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
          (_mm_loadu_pd(&a_6[k + 8]), _mm_set1_pd(ext_mv_idx_1)), _mm_mul_pd
          (_mm_loadu_pd(&a_6[k]), _mm_set1_pd(ext_mv_idx_0))), _mm_mul_pd
          (_mm_loadu_pd(&a_6[k + 16]), _mm_set1_pd(ext_mv_idx_2))), tmp_5));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      // Outputs for Function Call SubSystem: '<S1>/gmpc'
      //  Measurement upate to x(k|k)
      // '<S113>:1:155' ym_est = C(myindex,:)*xk + Dv(myindex,:)*v;
      // '<S113>:1:156' y_innov = ym - ym_est;
      for (k = 0; k < 6; k++) {
        // MATLAB Function: '<S112>/optimizer' incorporates:
        //   Inport: '<Root>/y'

        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 8; i_0++) {
          ext_mv_idx_2 += static_cast<real_T>(b_a[i + k]) * xk[i_0];
          i += 6;
        }

        y_innov[k] = rtU.y[k] - ext_mv_idx_2;
      }

      // '<S113>:1:157' xest = xk + M*y_innov;
      //  Real-time MV target override
      //  Note: utargetValue is a vector length p*nu.
      // '<S113>:1:162' if no_uref
      //  no external utarget
      // '<S113>:1:164' utargetValue = utarget;
      //  Real-time custom constraint override (scaled E/F/S)
      // '<S113>:1:173' if ~no_cc
      // '<S113>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
      // '<S113>:1:183' if isSimulation
      // '<S113>:1:214' else
      //  When doing code generation, use M code directly
      // '<S113>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
      // '<S113>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
      // '<S113>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
      // '<S113>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
      // '<S113>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
      // '<S113>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
      // '<S113>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
      // '<S113>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
      // '<S113>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
      // '<S113>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
      // '<S113>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
      // '<S113>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
      // '<S113>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
      for (k = 0; k < 8; k++) {
        // MATLAB Function: '<S112>/optimizer'
        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 6; i_0++) {
          ext_mv_idx_2 += c_a[i + k] * y_innov[i_0];
          i += 8;
        }

        xk_1[k] = xk[k] + ext_mv_idx_2;
      }

      // MATLAB Function: '<S112>/optimizer'
      (void)std::memset(&rtDW.dv[0], 0, 5166U * sizeof(real_T));
      (void)std::memset(&tmp[0], 0, 60U * sizeof(real_T));
      tmp_0[0] = 0.0;
      tmp_0[1] = 0.0;
      tmp_0[2] = 0.0;
      for (k = 0; k < 6; k++) {
        tmp_1[k] = 0.0;
      }

      tmp_2[0] = 0.034121465297356074;
      tmp_2[1] = 0.034121465297356074;
      tmp_2[2] = 0.034121465297356074;
      (void)std::memset(&rtDW.dv1[0], 0, 2520U * sizeof(real_T));
      for (k = 0; k < 6; k++) {
        tmp_3[k] = 1.0;
      }

      tmp_4[0] = 1.0;
      tmp_4[1] = 1.0;
      tmp_4[2] = 1.0;

      // Memory: '<S92>/Memory'
      (void)std::memcpy(&tmp_6[0], &rtDW.Memory_PreviousInput[0], 246U * sizeof
                        (boolean_T));

      // MATLAB Function: '<S112>/optimizer'
      (void)std::memcpy(&f_0[0], &f[0], sizeof(real_T) << 4UL);

      // Update for Memory: '<S92>/Memory' incorporates:
      //   Gain: '<S3>/Gain'
      //   Gain: '<S3>/Gain1'
      //   Gain: '<S3>/Gain2'
      //   Gain: '<S3>/Gain3'
      //   Inport: '<Root>/umax'
      //   MATLAB Function: '<S112>/optimizer'
      //   Math: '<S92>/Math Function'
      //   Math: '<S92>/Math Function1'

      mpcblock_optimizer(rseq, vseq, rtU.umax, rtb_Gain2, rtb_Gain3,
                         static_cast<int32_T>(std::round(rtb_sig)), xk_1, ext_mv,
                         tmp_6, c, d, e, rtDW.dv, tmp, tmp_0, tmp_1, 1, f_0, g,
                         rtb_ywt, rtb_Gain1, tmp_2, h, l, l, n, rtDW.dv1, o, q,
                         tmp_3, tmp_4, Sum2_c, &ext_mv_idx_2, rtb_useq_f,
                         &status, rtDW.Memory_PreviousInput);

      // '<S113>:1:231' if return_xseq || return_ovseq
      // '<S113>:1:232' [yseq, xseq] = mpc_computeSequence(isLTV, xest, useq, vseq, uoff, yoff, xoff, p, ny, nxQP, nv, A, Bu, Bv, C, Dv); 
      // '<S113>:1:238' if CustomEstimation
      // '<S113>:1:240' else
      //  update x[k+1|k], assuming that above u and v will be applied.
      // '<S113>:1:242' xk1 = A*xk + Bu*(u - uoff) + Bv*v + L*y_innov;
      // '<S113>:1:244' xk1 = xk1 + xoff;
      //  Updated state must include offset
      //  return xest in original value
      // '<S113>:1:247' xest = xest + xoff;
      // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S135>:1'
      // '<S135>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
      // '<S135>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
      // '<S135>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
      // '<S135>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
      //  Inputs (in BlockDataType except iA)
      //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
      // '<S135>:1:24' xk = convertDataType(xk0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    old_u:      last mv (calculated by MPC)
      // '<S135>:1:26' old_u = convertDataType(old_u0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    ym:         current measured output (used only with built-in KF)
      // '<S135>:1:28' ym = convertDataType(ym0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    ref:        output reference
      // '<S135>:1:30' ref = convertDataType(ref0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    md:         measured disturbance
      // '<S135>:1:32' md = convertDataType(md0,isDouble);
      //    umin:       run-time MV bound
      // '<S135>:1:34' umin = convertDataType(umin0,isDouble);
      //    umax:       run-time MV bound
      // '<S135>:1:36' umax = convertDataType(umax0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    ymin:       run-time OV bound
      // '<S135>:1:38' ymin = convertDataType(ymin0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    ymax:       run-time OV bound
      // '<S135>:1:40' ymax = convertDataType(ymax0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    E:          run-time mixed constraints
      // '<S135>:1:42' E = convertDataType(E0,isDouble);
      //    F:          run-time mixed constraints
      // '<S135>:1:44' F = convertDataType(F0,isDouble);
      //    G:          run-time mixed constraints
      // '<S135>:1:46' G = convertDataType(G0,isDouble);
      //    S:          run-time mixed constraints
      // '<S135>:1:48' S = convertDataType(S0,isDouble);
      //    switch_in:  if it matches "enable_value", MPC is active in control
      // '<S135>:1:50' switch_in = int32(switch_in0);
      //    ext_mv:     external last mv (actual)
      // '<S135>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    MVtarget:   MV reference
      // '<S135>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
      //    ywt:        run-time OV weights
      // '<S135>:1:56' ywt = convertDataType(ywt0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    uwt:        run-time MV weights
      // '<S135>:1:58' uwt = convertDataType(uwt0,isDouble);
      // '<S135>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S135>:1:252' if isa(u,'double')
      // '<S135>:1:253' y = u;
      //    duwt:       run-time DMV weights
      // '<S135>:1:60' duwt = convertDataType(duwt0,isDouble);
      //    ewt:     run-time Slack weights
      // '<S135>:1:62' ewt = convertDataType(ewt0,isDouble);
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
      // '<S135>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
      // '<S135>:1:96' isAdaptive = false;
      // '<S135>:1:97' isLTV = false;
      // '<S135>:1:98' ZERO = zeros('like',ref);
      // '<S135>:1:99' ONE = ones('like',ref);
      // '<S135>:1:100' hasMD = nv>int32(1);
      //  Pre-allocate all the MEX block outputs for the simulation mode
      // '<S135>:1:105' if isSimulation
      //  Get reference and MD signals -- accounting for previewing
      // '<S135>:1:119' if isSimulation
      // '<S135>:1:126' else
      //  When doing code generation, use M code directly
      // '<S135>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
      for (k = 0; k < 21; k++) {
        // MATLAB Function: '<S134>/optimizer'
        vseq[k] = 1.0;
      }

      for (k = 0; k < 20; k++) {
        for (i = 0; i < 6; i++) {
          // MATLAB Function: '<S134>/optimizer'
          rseq[i + k * static_cast<int32_T>(ny)] = rtDW.traj
            [(static_cast<int32_T>(rtDW.waypt) - 1) * 6 + i];
        }
      }

      // MATLAB Function: '<S134>/optimizer' incorporates:
      //   Gain: '<S114>/ext.mv_scale'
      //   UnitDelay: '<S114>/last_mv'

      //  External MV override.
      //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
      // '<S135>:1:133' old_u = old_u - uoff;
      // '<S135>:1:134' if no_mv
      // '<S135>:1:136' else
      // '<S135>:1:137' ext_mv = ext_mv - uoff;
      ext_mv[0] = rtP.extmv_scale_Gain_o[0] * rtb_Delay2[0];
      ext_mv[1] = rtP.extmv_scale_Gain_o[1] * rtb_Delay2[1];
      ext_mv[2] = rtP.extmv_scale_Gain_o[2] * rtb_Delay2[2];

      //  Bias correction
      // '<S135>:1:138' delmv = ext_mv - old_u;
      // '<S135>:1:139' old_u = ext_mv;
      //  Obtain x[k|k]
      // '<S135>:1:143' xk = xk - xoff;
      //  Remove offset
      // '<S135>:1:144' if CustomEstimation
      // '<S135>:1:147' else
      //  Default state estimation.
      //  Scale measured output and remove offset.
      // '<S135>:1:150' ym = ym.*RYscale(myindex) - myoff;
      //  Correct x(k|k-1) for possible external mv override.
      //  NOTE:  Offset was removed from x[k|k-1] at k=0.
      // '<S135>:1:153' xk = xk + Bu*delmv;
      ext_mv_idx_0 = ext_mv[0] - rtDW.last_mv_DSTATE_i[0];
      ext_mv_idx_1 = ext_mv[1] - rtDW.last_mv_DSTATE_i[1];
      ext_mv_idx_2 = ext_mv[2] - rtDW.last_mv_DSTATE_i[2];

      // End of Outputs for SubSystem: '<S1>/gmpc'
      for (k = 0; k <= 8; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Memory: '<S114>/last_x' incorporates:
        //   MATLAB Function: '<S134>/optimizer'

        tmp_5 = _mm_loadu_pd(&rtDW.last_x_PreviousInput_h[k]);

        // MATLAB Function: '<S134>/optimizer'
        (void)_mm_storeu_pd(&xk_0[k], _mm_add_pd(_mm_add_pd(_mm_add_pd
          (_mm_mul_pd(_mm_loadu_pd(&a_7[k + 10]), _mm_set1_pd(ext_mv_idx_1)),
           _mm_mul_pd(_mm_loadu_pd(&a_7[k]), _mm_set1_pd(ext_mv_idx_0))),
          _mm_mul_pd(_mm_loadu_pd(&a_7[k + 20]), _mm_set1_pd(ext_mv_idx_2))),
          tmp_5));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      // Outputs for Function Call SubSystem: '<S1>/gmpc'
      //  Measurement upate to x(k|k)
      // '<S135>:1:155' ym_est = C(myindex,:)*xk + Dv(myindex,:)*v;
      // '<S135>:1:156' y_innov = ym - ym_est;
      for (k = 0; k < 6; k++) {
        // MATLAB Function: '<S134>/optimizer' incorporates:
        //   Inport: '<Root>/y'

        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 10; i_0++) {
          ext_mv_idx_2 += static_cast<real_T>(b_a_0[i + k]) * xk_0[i_0];
          i += 6;
        }

        y_innov_0[k] = rtU.y[k] - ext_mv_idx_2;
      }

      // '<S135>:1:157' xest = xk + M*y_innov;
      //  Real-time MV target override
      //  Note: utargetValue is a vector length p*nu.
      // '<S135>:1:162' if no_uref
      //  no external utarget
      // '<S135>:1:164' utargetValue = utarget;
      //  Real-time custom constraint override (scaled E/F/S)
      // '<S135>:1:173' if ~no_cc
      // '<S135>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
      // '<S135>:1:183' if isSimulation
      // '<S135>:1:214' else
      //  When doing code generation, use M code directly
      // '<S135>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
      // '<S135>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
      // '<S135>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
      // '<S135>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
      // '<S135>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
      // '<S135>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
      // '<S135>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
      // '<S135>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
      // '<S135>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
      // '<S135>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
      // '<S135>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
      // '<S135>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
      // '<S135>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
      for (k = 0; k < 10; k++) {
        // MATLAB Function: '<S134>/optimizer'
        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 6; i_0++) {
          ext_mv_idx_2 += c_a_0[i + k] * y_innov_0[i_0];
          i += 10;
        }

        xk_2[k] = xk_0[k] + ext_mv_idx_2;
      }

      // MATLAB Function: '<S134>/optimizer'
      (void)std::memset(&rtDW.dv[0], 0, 5166U * sizeof(real_T));
      (void)std::memset(&tmp[0], 0, 60U * sizeof(real_T));
      tmp_0[0] = 0.0;
      tmp_0[1] = 0.0;
      tmp_0[2] = 0.0;
      for (k = 0; k < 6; k++) {
        tmp_1[k] = 0.0;
      }

      tmp_2[0] = 0.034121465297356074;
      tmp_2[1] = 0.034121465297356074;
      tmp_2[2] = 0.034121465297356074;
      (void)std::memset(&rtDW.dv1[0], 0, 2520U * sizeof(real_T));
      for (k = 0; k < 6; k++) {
        tmp_3[k] = 1.0;
      }

      tmp_4[0] = 1.0;
      tmp_4[1] = 1.0;
      tmp_4[2] = 1.0;

      // Memory: '<S114>/Memory'
      (void)std::memcpy(&tmp_6[0], &rtDW.Memory_PreviousInput_g[0], 246U *
                        sizeof(boolean_T));

      // MATLAB Function: '<S134>/optimizer'
      (void)std::memcpy(&f_0[0], &f[0], sizeof(real_T) << 4UL);

      // Update for Memory: '<S114>/Memory' incorporates:
      //   Gain: '<S3>/Gain'
      //   Gain: '<S3>/Gain1'
      //   Gain: '<S3>/Gain2'
      //   Gain: '<S3>/Gain3'
      //   Inport: '<Root>/umax'
      //   MATLAB Function: '<S134>/optimizer'
      //   Math: '<S114>/Math Function'
      //   Math: '<S114>/Math Function1'

      mpcblock_optimizer_h(rseq, vseq, rtU.umax, rtb_Gain2, rtb_Gain3,
                           static_cast<int32_T>(std::round(rtb_sig)), xk_2,
                           ext_mv, tmp_6, c, d_0, e_0, rtDW.dv, tmp, tmp_0,
                           tmp_1, 2, f_0, g_0, rtb_ywt, rtb_Gain1, tmp_2, h, l_0,
                           l_0, n_0, rtDW.dv1, o, q, tmp_3, tmp_4, rtb_u_l,
                           &ext_mv_idx_2, rtb_useq_f, &status,
                           rtDW.Memory_PreviousInput_g);

      // '<S135>:1:231' if return_xseq || return_ovseq
      // '<S135>:1:232' [yseq, xseq] = mpc_computeSequence(isLTV, xest, useq, vseq, uoff, yoff, xoff, p, ny, nxQP, nv, A, Bu, Bv, C, Dv); 
      // '<S135>:1:238' if CustomEstimation
      // '<S135>:1:240' else
      //  update x[k+1|k], assuming that above u and v will be applied.
      // '<S135>:1:242' xk1 = A*xk + Bu*(u - uoff) + Bv*v + L*y_innov;
      // '<S135>:1:244' xk1 = xk1 + xoff;
      //  Updated state must include offset
      //  return xest in original value
      // '<S135>:1:247' xest = xest + xoff;
      // MATLAB Function 'MPC Controller/MPC/optimizer/optimizer': '<S157>:1'
      // '<S157>:1:17' coder.extrinsic('mpcblock_optimizer_double_mex');
      // '<S157>:1:18' coder.extrinsic('mpcblock_optimizer_single_mex');
      // '<S157>:1:19' coder.extrinsic('mpcblock_refmd_double_mex');
      // '<S157>:1:20' coder.extrinsic('mpcblock_refmd_single_mex');
      //  Inputs (in BlockDataType except iA)
      //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
      // '<S157>:1:24' xk = convertDataType(xk0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    old_u:      last mv (calculated by MPC)
      // '<S157>:1:26' old_u = convertDataType(old_u0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    ym:         current measured output (used only with built-in KF)
      // '<S157>:1:28' ym = convertDataType(ym0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    ref:        output reference
      // '<S157>:1:30' ref = convertDataType(ref0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    md:         measured disturbance
      // '<S157>:1:32' md = convertDataType(md0,isDouble);
      //    umin:       run-time MV bound
      // '<S157>:1:34' umin = convertDataType(umin0,isDouble);
      //    umax:       run-time MV bound
      // '<S157>:1:36' umax = convertDataType(umax0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    ymin:       run-time OV bound
      // '<S157>:1:38' ymin = convertDataType(ymin0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    ymax:       run-time OV bound
      // '<S157>:1:40' ymax = convertDataType(ymax0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    E:          run-time mixed constraints
      // '<S157>:1:42' E = convertDataType(E0,isDouble);
      //    F:          run-time mixed constraints
      // '<S157>:1:44' F = convertDataType(F0,isDouble);
      //    G:          run-time mixed constraints
      // '<S157>:1:46' G = convertDataType(G0,isDouble);
      //    S:          run-time mixed constraints
      // '<S157>:1:48' S = convertDataType(S0,isDouble);
      //    switch_in:  if it matches "enable_value", MPC is active in control
      // '<S157>:1:50' switch_in = int32(switch_in0);
      //    ext_mv:     external last mv (actual)
      // '<S157>:1:52' ext_mv = convertDataType(ext_mv0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    MVtarget:   MV reference
      // '<S157>:1:54' MVtarget = convertDataType(MVtarget0,isDouble);
      //    ywt:        run-time OV weights
      // '<S157>:1:56' ywt = convertDataType(ywt0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    uwt:        run-time MV weights
      // '<S157>:1:58' uwt = convertDataType(uwt0,isDouble);
      // '<S157>:1:250' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S157>:1:252' if isa(u,'double')
      // '<S157>:1:253' y = u;
      //    duwt:       run-time DMV weights
      // '<S157>:1:60' duwt = convertDataType(duwt0,isDouble);
      //    ewt:     run-time Slack weights
      // '<S157>:1:62' ewt = convertDataType(ewt0,isDouble);
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
      // '<S157>:1:95' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
      // '<S157>:1:96' isAdaptive = false;
      // '<S157>:1:97' isLTV = false;
      // '<S157>:1:98' ZERO = zeros('like',ref);
      // '<S157>:1:99' ONE = ones('like',ref);
      // '<S157>:1:100' hasMD = nv>int32(1);
      //  Pre-allocate all the MEX block outputs for the simulation mode
      // '<S157>:1:105' if isSimulation
      //  Get reference and MD signals -- accounting for previewing
      // '<S157>:1:119' if isSimulation
      // '<S157>:1:126' else
      //  When doing code generation, use M code directly
      // '<S157>:1:128' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
      for (k = 0; k < 21; k++) {
        // MATLAB Function: '<S156>/optimizer'
        vseq[k] = 1.0;
      }

      for (k = 0; k < 20; k++) {
        for (i = 0; i < 6; i++) {
          // MATLAB Function: '<S156>/optimizer'
          rseq[i + k * static_cast<int32_T>(ny)] = rtDW.traj
            [(static_cast<int32_T>(rtDW.waypt) - 1) * 6 + i];
        }
      }

      // Gain: '<S136>/ext.mv_scale'
      //  External MV override.
      //  NOTE: old_u and ext_mv input signals are dimensionless but include offset 
      // '<S157>:1:133' old_u = old_u - uoff;
      // '<S157>:1:134' if no_mv
      // '<S157>:1:136' else
      // '<S157>:1:137' ext_mv = ext_mv - uoff;
      status = rtP.extmv_scale_Gain_a[0] * rtb_Delay2[0];

      // MATLAB Function: '<S156>/optimizer' incorporates:
      //   Gain: '<S136>/ext.mv_scale'

      rtb_Delay2[0] = status;

      // Gain: '<S136>/ext.mv_scale'
      rtb_Delay2_tmp = rtP.extmv_scale_Gain_a[1] * rtb_Delay2[1];

      // MATLAB Function: '<S156>/optimizer' incorporates:
      //   Gain: '<S136>/ext.mv_scale'

      rtb_Delay2[1] = rtb_Delay2_tmp;

      // Gain: '<S136>/ext.mv_scale'
      ext_mv_idx_2 = rtP.extmv_scale_Gain_a[2] * rtb_Delay2[2];

      // MATLAB Function: '<S156>/optimizer' incorporates:
      //   Gain: '<S136>/ext.mv_scale'
      //   UnitDelay: '<S136>/last_mv'

      rtb_Delay2[2] = ext_mv_idx_2;

      //  Bias correction
      // '<S157>:1:138' delmv = ext_mv - old_u;
      // '<S157>:1:139' old_u = ext_mv;
      //  Obtain x[k|k]
      // '<S157>:1:143' xk = xk - xoff;
      //  Remove offset
      // '<S157>:1:144' if CustomEstimation
      // '<S157>:1:147' else
      //  Default state estimation.
      //  Scale measured output and remove offset.
      // '<S157>:1:150' ym = ym.*RYscale(myindex) - myoff;
      //  Correct x(k|k-1) for possible external mv override.
      //  NOTE:  Offset was removed from x[k|k-1] at k=0.
      // '<S157>:1:153' xk = xk + Bu*delmv;
      status -= rtDW.last_mv_DSTATE_g[0];
      rtb_Delay2_tmp -= rtDW.last_mv_DSTATE_g[1];
      ext_mv_idx_2 -= rtDW.last_mv_DSTATE_g[2];

      // End of Outputs for SubSystem: '<S1>/gmpc'
      for (k = 0; k <= 8; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/gmpc'
        // Memory: '<S136>/last_x' incorporates:
        //   MATLAB Function: '<S156>/optimizer'

        tmp_5 = _mm_loadu_pd(&rtDW.last_x_PreviousInput_c[k]);

        // MATLAB Function: '<S156>/optimizer'
        (void)_mm_storeu_pd(&xk_2[k], _mm_add_pd(_mm_add_pd(_mm_add_pd
          (_mm_mul_pd(_mm_loadu_pd(&a_8[k + 10]), _mm_set1_pd(rtb_Delay2_tmp)),
           _mm_mul_pd(_mm_loadu_pd(&a_8[k]), _mm_set1_pd(status))), _mm_mul_pd
          (_mm_loadu_pd(&a_8[k + 20]), _mm_set1_pd(ext_mv_idx_2))), tmp_5));

        // End of Outputs for SubSystem: '<S1>/gmpc'
      }

      // Outputs for Function Call SubSystem: '<S1>/gmpc'
      //  Measurement upate to x(k|k)
      // '<S157>:1:155' ym_est = C(myindex,:)*xk + Dv(myindex,:)*v;
      // '<S157>:1:156' y_innov = ym - ym_est;
      for (k = 0; k < 6; k++) {
        // MATLAB Function: '<S156>/optimizer' incorporates:
        //   Inport: '<Root>/y'

        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 10; i_0++) {
          ext_mv_idx_2 += static_cast<real_T>(b_a_1[i + k]) * xk_2[i_0];
          i += 6;
        }

        y_innov_1[k] = rtU.y[k] - ext_mv_idx_2;
      }

      // '<S157>:1:157' xest = xk + M*y_innov;
      //  Real-time MV target override
      //  Note: utargetValue is a vector length p*nu.
      // '<S157>:1:162' if no_uref
      //  no external utarget
      // '<S157>:1:164' utargetValue = utarget;
      //  Real-time custom constraint override (scaled E/F/S)
      // '<S157>:1:173' if ~no_cc
      // '<S157>:1:182' return_sequence = return_mvseq || return_xseq || return_ovseq; 
      // '<S157>:1:183' if isSimulation
      // '<S157>:1:214' else
      //  When doing code generation, use M code directly
      // '<S157>:1:216' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
      // '<S157>:1:217'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
      // '<S157>:1:218'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
      // '<S157>:1:219'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
      // '<S157>:1:220'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
      // '<S157>:1:221'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
      // '<S157>:1:222'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
      // '<S157>:1:223'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
      // '<S157>:1:224'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
      // '<S157>:1:225'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
      // '<S157>:1:226'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
      // '<S157>:1:227'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, ... 
      // '<S157>:1:228'             isHyb, Mdis, Ndis, Vdis, numdis, maxdis);
      for (k = 0; k < 10; k++) {
        // MATLAB Function: '<S156>/optimizer'
        ext_mv_idx_2 = 0.0;
        i = 0;
        for (i_0 = 0; i_0 < 6; i_0++) {
          ext_mv_idx_2 += c_a_1[i + k] * y_innov_1[i_0];
          i += 10;
        }

        xk_3[k] = xk_2[k] + ext_mv_idx_2;
      }

      // MATLAB Function: '<S156>/optimizer'
      (void)std::memset(&rtDW.dv[0], 0, 5166U * sizeof(real_T));
      (void)std::memset(&tmp[0], 0, 60U * sizeof(real_T));
      tmp_0[0] = 0.0;
      tmp_0[1] = 0.0;
      tmp_0[2] = 0.0;
      for (k = 0; k < 6; k++) {
        tmp_1[k] = 0.0;
      }

      tmp_2[0] = 0.034121465297356074;
      tmp_2[1] = 0.034121465297356074;
      tmp_2[2] = 0.034121465297356074;
      (void)std::memset(&rtDW.dv1[0], 0, 2520U * sizeof(real_T));
      for (k = 0; k < 6; k++) {
        tmp_3[k] = 1.0;
      }

      tmp_4[0] = 1.0;
      tmp_4[1] = 1.0;
      tmp_4[2] = 1.0;

      // Memory: '<S136>/Memory'
      (void)std::memcpy(&tmp_6[0], &rtDW.Memory_PreviousInput_m[0], 246U *
                        sizeof(boolean_T));

      // MATLAB Function: '<S156>/optimizer'
      (void)std::memcpy(&f_0[0], &f[0], sizeof(real_T) << 4UL);

      // Update for Memory: '<S136>/Memory' incorporates:
      //   Gain: '<S3>/Gain'
      //   Gain: '<S3>/Gain1'
      //   Gain: '<S3>/Gain2'
      //   Gain: '<S3>/Gain3'
      //   Inport: '<Root>/umax'
      //   MATLAB Function: '<S156>/optimizer'
      //   Math: '<S136>/Math Function'
      //   Math: '<S136>/Math Function1'

      mpcblock_optimizer_h(rseq, vseq, rtU.umax, rtb_Gain2, rtb_Gain3,
                           static_cast<int32_T>(std::round(rtb_sig)), xk_3,
                           rtb_Delay2, tmp_6, c, d_1, e_1, rtDW.dv, tmp, tmp_0,
                           tmp_1, 3, f_0, g_1, rtb_ywt, rtb_Gain1, tmp_2, h, l_1,
                           l_1, n_1, rtDW.dv1, o, q, tmp_3, tmp_4, ext_mv,
                           &ext_mv_idx_2, rtb_useq_f, &status,
                           rtDW.Memory_PreviousInput_m);

      // MultiPortSwitch: '<S88>/u0_switch' incorporates:
      //   Constant: '<S114>/constant'
      //   Constant: '<S136>/constant'
      //   Constant: '<S88>/u0_zero'
      //   Constant: '<S92>/constant'
      //   Gain: '<S114>/umin_scale2'
      //   Gain: '<S136>/umin_scale2'
      //   Gain: '<S92>/umin_scale2'

      // '<S157>:1:231' if return_xseq || return_ovseq
      // '<S157>:1:232' [yseq, xseq] = mpc_computeSequence(isLTV, xest, useq, vseq, uoff, yoff, xoff, p, ny, nxQP, nv, A, Bu, Bv, C, Dv); 
      // '<S157>:1:238' if CustomEstimation
      // '<S157>:1:240' else
      //  update x[k+1|k], assuming that above u and v will be applied.
      // '<S157>:1:242' xk1 = A*xk + Bu*(u - uoff) + Bv*v + L*y_innov;
      // '<S157>:1:244' xk1 = xk1 + xoff;
      //  Updated state must include offset
      //  return xest in original value
      // '<S157>:1:247' xest = xest + xoff;
      switch (static_cast<int32_T>(rtb_sig)) {
       case 1:
        rtb_Gain1[0] = rtP.umin_scale2_Gain[0] * rtP.constant_Value[0];
        rtb_Gain1[1] = rtP.umin_scale2_Gain[1] * rtP.constant_Value[1];
        rtb_Gain1[2] = rtP.umin_scale2_Gain[2] * rtP.constant_Value[2];
        break;

       case 2:
        rtb_Gain1[0] = rtP.umin_scale2_Gain_m[0] * rtP.constant_Value_l[0];
        rtb_Gain1[1] = rtP.umin_scale2_Gain_m[1] * rtP.constant_Value_l[1];
        rtb_Gain1[2] = rtP.umin_scale2_Gain_m[2] * rtP.constant_Value_l[2];
        break;

       case 3:
        rtb_Gain1[0] = rtP.umin_scale2_Gain_l[0] * rtP.constant_Value_e[0];
        rtb_Gain1[1] = rtP.umin_scale2_Gain_l[1] * rtP.constant_Value_e[1];
        rtb_Gain1[2] = rtP.umin_scale2_Gain_l[2] * rtP.constant_Value_e[2];
        break;

       default:
        rtb_Gain1[0] = rtP.u0_zero_Value[0];
        rtb_Gain1[1] = rtP.u0_zero_Value[1];
        rtb_Gain1[2] = rtP.u0_zero_Value[2];
        break;
      }

      // End of MultiPortSwitch: '<S88>/u0_switch'

      // Delay: '<S88>/u[k-1]'
      if (rtDW.icLoad) {
        rtDW.uk1_DSTATE[0] = rtb_Gain1[0];
        rtDW.uk1_DSTATE[1] = rtb_Gain1[1];
        rtDW.uk1_DSTATE[2] = rtb_Gain1[2];
      }

      // MultiPortSwitch: '<S88>/mv_switch' incorporates:
      //   Delay: '<S88>/u[k-1]'
      //   Gain: '<S114>/umin_scale1'
      //   Gain: '<S136>/umin_scale1'
      //   Gain: '<S92>/umin_scale1'

      switch (static_cast<int32_T>(rtb_sig)) {
       case 1:
        rtb_Gain1[0] = rtP.umin_scale1_Gain[0] * Sum2_c[0];
        rtb_Gain1[1] = rtP.umin_scale1_Gain[1] * Sum2_c[1];
        rtb_Gain1[2] = rtP.umin_scale1_Gain[2] * Sum2_c[2];
        break;

       case 2:
        rtb_Gain1[0] = rtP.umin_scale1_Gain_l[0] * rtb_u_l[0];
        rtb_Gain1[1] = rtP.umin_scale1_Gain_l[1] * rtb_u_l[1];
        rtb_Gain1[2] = rtP.umin_scale1_Gain_l[2] * rtb_u_l[2];
        break;

       case 3:
        rtb_Gain1[0] = rtP.umin_scale1_Gain_f[0] * ext_mv[0];
        rtb_Gain1[1] = rtP.umin_scale1_Gain_f[1] * ext_mv[1];
        rtb_Gain1[2] = rtP.umin_scale1_Gain_f[2] * ext_mv[2];
        break;

       default:
        rtb_Gain1[0] = rtDW.uk1_DSTATE[0];
        rtb_Gain1[1] = rtDW.uk1_DSTATE[1];
        rtb_Gain1[2] = rtDW.uk1_DSTATE[2];
        break;
      }

      // End of MultiPortSwitch: '<S88>/mv_switch'

      // Saturate: '<S3>/Saturation'
      if (rtb_Gain1[0] > rtP.Saturation_UpperSat_d) {
        rtb_sig = rtP.Saturation_UpperSat_d;
      } else if (rtb_Gain1[0] < rtP.Saturation_LowerSat_k) {
        rtb_sig = rtP.Saturation_LowerSat_k;
      } else {
        rtb_sig = rtb_Gain1[0];
      }

      if (rtb_Gain1[1] > rtP.Saturation_UpperSat_d) {
        status = rtP.Saturation_UpperSat_d;
      } else if (rtb_Gain1[1] < rtP.Saturation_LowerSat_k) {
        status = rtP.Saturation_LowerSat_k;
      } else {
        status = rtb_Gain1[1];
      }

      if (rtb_Gain1[2] > rtP.Saturation_UpperSat_d) {
        rtb_Delay2_tmp = rtP.Saturation_UpperSat_d;
      } else if (rtb_Gain1[2] < rtP.Saturation_LowerSat_k) {
        rtb_Delay2_tmp = rtP.Saturation_LowerSat_k;
      } else {
        rtb_Delay2_tmp = rtb_Gain1[2];
      }

      // Update for Delay: '<S3>/Delay2' incorporates:
      //   Saturate: '<S3>/Saturation'

      rtDW.Delay2_DSTATE[0] = rtb_sig;

      // Update for UnitDelay: '<S92>/last_mv'
      rtDW.last_mv_DSTATE[0] = Sum2_c[0];

      // MATLAB Function: '<S112>/optimizer' incorporates:
      //   UnitDelay: '<S92>/last_mv'

      ext_mv_idx_0 = Sum2_c[0];

      // Update for Delay: '<S3>/Delay2' incorporates:
      //   Saturate: '<S3>/Saturation'

      rtDW.Delay2_DSTATE[1] = status;

      // Update for UnitDelay: '<S92>/last_mv'
      rtDW.last_mv_DSTATE[1] = Sum2_c[1];

      // MATLAB Function: '<S112>/optimizer' incorporates:
      //   UnitDelay: '<S92>/last_mv'

      ext_mv_idx_1 = Sum2_c[1];

      // Update for UnitDelay: '<S92>/last_mv'
      Sum2_c_0 = Sum2_c[2];

      // Update for Delay: '<S3>/Delay2' incorporates:
      //   Saturate: '<S3>/Saturation'

      rtDW.Delay2_DSTATE[2] = rtb_Delay2_tmp;

      // Update for UnitDelay: '<S92>/last_mv'
      rtDW.last_mv_DSTATE[2] = Sum2_c[2];
      for (k = 0; k < 8; k++) {
        // MATLAB Function: '<S112>/optimizer'
        xk_1[k] = 0.0;
        for (i = 0; i < 8; i++) {
          xk_1[k] += b_A[(i << 3UL) + k] * xk[i];
        }

        // Update for Memory: '<S92>/last_x' incorporates:
        //   MATLAB Function: '<S112>/optimizer'

        ext_mv_idx_2 = 0.0;
        for (i = 0; i < 6; i++) {
          ext_mv_idx_2 += d_a[(i << 3UL) + k] * y_innov[i];
        }

        rtDW.last_x_PreviousInput[k] = (((a_9[k + 8] * ext_mv_idx_1 + a_9[k] *
          ext_mv_idx_0) + a_9[k + 16] * Sum2_c_0) + xk_1[k]) + ext_mv_idx_2;

        // End of Update for Memory: '<S92>/last_x'
      }

      // Update for UnitDelay: '<S114>/last_mv'
      rtDW.last_mv_DSTATE_i[0] = rtb_u_l[0];

      // MATLAB Function: '<S134>/optimizer' incorporates:
      //   UnitDelay: '<S114>/last_mv'

      ext_mv_idx_0 = rtb_u_l[0];

      // Update for UnitDelay: '<S114>/last_mv'
      rtDW.last_mv_DSTATE_i[1] = rtb_u_l[1];

      // MATLAB Function: '<S134>/optimizer' incorporates:
      //   UnitDelay: '<S114>/last_mv'

      ext_mv_idx_1 = rtb_u_l[1];

      // Update for UnitDelay: '<S114>/last_mv'
      Sum2_c_0 = rtb_u_l[2];
      rtDW.last_mv_DSTATE_i[2] = rtb_u_l[2];
      for (k = 0; k < 10; k++) {
        // MATLAB Function: '<S134>/optimizer'
        xk_3[k] = 0.0;
        for (i = 0; i < 10; i++) {
          xk_3[k] += b_A_0[10 * i + k] * xk_0[i];
        }

        // Update for Memory: '<S114>/last_x' incorporates:
        //   MATLAB Function: '<S134>/optimizer'

        ext_mv_idx_2 = 0.0;
        for (i = 0; i < 6; i++) {
          ext_mv_idx_2 += d_a_0[10 * i + k] * y_innov_0[i];
        }

        rtDW.last_x_PreviousInput_h[k] = (((a_a[k + 10] * ext_mv_idx_1 + a_a[k] *
          ext_mv_idx_0) + a_a[k + 20] * Sum2_c_0) + xk_3[k]) + ext_mv_idx_2;

        // End of Update for Memory: '<S114>/last_x'
      }

      // Update for UnitDelay: '<S136>/last_mv'
      rtDW.last_mv_DSTATE_g[0] = ext_mv[0];

      // MATLAB Function: '<S156>/optimizer' incorporates:
      //   UnitDelay: '<S136>/last_mv'

      ext_mv_idx_0 = ext_mv[0];

      // Update for UnitDelay: '<S136>/last_mv'
      rtDW.last_mv_DSTATE_g[1] = ext_mv[1];

      // MATLAB Function: '<S156>/optimizer' incorporates:
      //   UnitDelay: '<S136>/last_mv'

      ext_mv_idx_1 = ext_mv[1];

      // Update for UnitDelay: '<S136>/last_mv'
      Sum2_c_0 = ext_mv[2];
      rtDW.last_mv_DSTATE_g[2] = ext_mv[2];
      for (k = 0; k < 10; k++) {
        // MATLAB Function: '<S156>/optimizer'
        xk_3[k] = 0.0;
        for (i = 0; i < 10; i++) {
          xk_3[k] += b_A_1[10 * i + k] * xk_2[i];
        }

        // Update for Memory: '<S136>/last_x' incorporates:
        //   MATLAB Function: '<S156>/optimizer'

        ext_mv_idx_2 = 0.0;
        for (i = 0; i < 6; i++) {
          ext_mv_idx_2 += d_a_1[10 * i + k] * y_innov_1[i];
        }

        rtDW.last_x_PreviousInput_c[k] = (((a_b[k + 10] * ext_mv_idx_1 + a_b[k] *
          ext_mv_idx_0) + a_b[k + 20] * Sum2_c_0) + xk_3[k]) + ext_mv_idx_2;

        // End of Update for Memory: '<S136>/last_x'
      }

      // Update for Delay: '<S88>/u[k-1]'
      rtDW.icLoad = false;
      rtDW.uk1_DSTATE[0] = rtb_Gain1[0];
      rtDW.uk1_DSTATE[1] = rtb_Gain1[1];
      rtDW.uk1_DSTATE[2] = rtb_Gain1[2];

      // Outport: '<Root>/u' incorporates:
      //   Saturate: '<S3>/Saturation'

      rtY.u[0] = rtb_sig;
      rtY.u[1] = status;
      rtY.u[2] = rtb_Delay2_tmp;

      // End of Outputs for SubSystem: '<S1>/gmpc'
      for (i = 0; i < 6; i++) {
        // Outport: '<Root>/ywt' incorporates:
        //   Gain: '<S3>/Gain'

        rtY.ywt[i] = rtb_ywt[i];

        // Outport: '<Root>/currTraj'
        rtY.currTraj[i] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6 +
          i];
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
      rtY.P_e[i] = static_cast<real_T>(tmp_3[i]);
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
    //   SubSystem: '<S1>/gmpc'

    // InitializeConditions for Memory: '<S92>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput[0],
                      &rtP.Memory_InitialCondition_c[0], 246U * sizeof(boolean_T));

    // InitializeConditions for Delay: '<S3>/Delay1'
    rtDW.Delay1_DSTATE[0] = rtP.uwt0[0];
    rtDW.Delay1_DSTATE[1] = rtP.uwt0[1];
    rtDW.Delay1_DSTATE[2] = rtP.uwt0[2];
    for (i = 0; i < 6; i++) {
      // InitializeConditions for Delay: '<S3>/Delay'
      rtDW.Delay_DSTATE[i] = rtP.ywt0[i];
    }

    // InitializeConditions for Delay: '<S3>/Delay2'
    rtDW.Delay2_DSTATE[0] = rtP.Delay2_InitialCondition;

    // InitializeConditions for UnitDelay: '<S92>/last_mv'
    rtDW.last_mv_DSTATE[0] = rtP.last_mv_InitialCondition_n[0];

    // InitializeConditions for Delay: '<S3>/Delay2'
    rtDW.Delay2_DSTATE[1] = rtP.Delay2_InitialCondition;

    // InitializeConditions for UnitDelay: '<S92>/last_mv'
    rtDW.last_mv_DSTATE[1] = rtP.last_mv_InitialCondition_n[1];

    // InitializeConditions for Delay: '<S3>/Delay2'
    rtDW.Delay2_DSTATE[2] = rtP.Delay2_InitialCondition;

    // InitializeConditions for UnitDelay: '<S92>/last_mv'
    rtDW.last_mv_DSTATE[2] = rtP.last_mv_InitialCondition_n[2];

    // InitializeConditions for Memory: '<S92>/last_x'
    (void)std::memcpy(&rtDW.last_x_PreviousInput[0],
                      &rtP.last_x_InitialCondition[0], sizeof(real_T) << 3UL);

    // InitializeConditions for UnitDelay: '<S114>/last_mv'
    rtDW.last_mv_DSTATE_i[0] = rtP.last_mv_InitialCondition_l[0];
    rtDW.last_mv_DSTATE_i[1] = rtP.last_mv_InitialCondition_l[1];
    rtDW.last_mv_DSTATE_i[2] = rtP.last_mv_InitialCondition_l[2];

    // InitializeConditions for Memory: '<S114>/last_x'
    (void)std::memcpy(&rtDW.last_x_PreviousInput_h[0],
                      &rtP.last_x_InitialCondition_k[0], 10U * sizeof(real_T));

    // InitializeConditions for Memory: '<S114>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_g[0],
                      &rtP.Memory_InitialCondition_cu[0], 246U * sizeof
                      (boolean_T));

    // InitializeConditions for UnitDelay: '<S136>/last_mv'
    rtDW.last_mv_DSTATE_g[0] = rtP.last_mv_InitialCondition_g[0];
    rtDW.last_mv_DSTATE_g[1] = rtP.last_mv_InitialCondition_g[1];
    rtDW.last_mv_DSTATE_g[2] = rtP.last_mv_InitialCondition_g[2];

    // InitializeConditions for Memory: '<S136>/last_x'
    (void)std::memcpy(&rtDW.last_x_PreviousInput_c[0],
                      &rtP.last_x_InitialCondition_p[0], 10U * sizeof(real_T));

    // InitializeConditions for Memory: '<S136>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_m[0],
                      &rtP.Memory_InitialCondition_f[0], 246U * sizeof(boolean_T));

    // InitializeConditions for Delay: '<S88>/u[k-1]'
    rtDW.icLoad = true;

    // SystemInitialize for MATLAB Function: '<S3>/MATLAB Function1'
    // 'gainSchSig_:4' sigPrev = 1;
    rtDW.sigPrev = 1.0;
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
