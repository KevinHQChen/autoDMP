//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.2522
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Sun Oct  1 04:42:03 2023
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

// Function for Chart: '<Root>/SupervisoryController'
boolean_T SupervisoryController::any_d(const real_T x[3])
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

void SupervisoryController::binary_expand_op_n(real_T in1[24], int32_T in2,
  const real_T in3[24], int32_T in4, int32_T in5, const real_T in6[24], int32_T
  in7, int32_T in8)
{
  int32_T stride_0_0;
  int32_T stride_1_0;
  int32_T tmp;

  // Chart: '<Root>/SupervisoryController' incorporates:
  //   TriggerPort: '<S1>/measAvail'

  // Chart: '<Root>/SupervisoryController' incorporates:
  //   SubSystem: '<S1>/ampc'

  // MATLAB Function: '<S2>/MATLAB Function2' incorporates:
  //   Outport: '<Root>/theta'

  stride_0_0 = (in5 - in4) + 1 != 1 ? static_cast<int32_T>(1) :
    static_cast<int32_T>(0);
  stride_1_0 = (in8 - in7) + 1 != 1 ? static_cast<int32_T>(1) :
    static_cast<int32_T>(0);
  tmp = in2 << 2UL;
  in1[tmp] = in3[in4] * in6[in7];
  in1[tmp + 1] = in3[in4 + stride_0_0] * in6[in7 + stride_1_0];
  in1[tmp + 2] = in3[(stride_0_0 << 1UL) + in4] * in6[(stride_1_0 << 1UL) + in7];
  in1[tmp + 3] = in3[3 * stride_0_0 + in4] * in6[3 * stride_1_0 + in7];
}

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
void SupervisoryController::mpc_checkhessian(real_T b_H[16], real_T L[16],
  real_T *BadH)
{
  real_T varargin_1[4];
  int32_T Tries;
  int8_T b[16];
  boolean_T guard1{ false };

  *BadH = 0.0;
  (void)std::memcpy(&L[0], &b_H[0], sizeof(real_T) << 4UL);
  Tries = xpotrf(L);
  guard1 = false;
  if (Tries == 0) {
    varargin_1[0] = L[0];
    varargin_1[1] = L[5];
    varargin_1[2] = L[10];
    varargin_1[3] = L[15];
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
          L[j] = b_H[j];
        }

        j = xpotrf(L);
        guard2 = false;
        if (j == 0) {
          varargin_1[0] = L[0];
          varargin_1[1] = L[5];
          varargin_1[2] = L[10];
          varargin_1[3] = L[15];
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
void SupervisoryController::trisolve(const real_T b_A[16], real_T b_B[16])
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
      tmp_0 = b_B[tmp];
      if (tmp_0 != 0.0) {
        b_B[tmp] = tmp_0 / b_A[b_k + kAcol];
        for (int32_T i{b_k + 2}; i < 5; i++) {
          int32_T tmp_1;
          tmp_1 = (i + jBcol) - 1;
          b_B[tmp_1] -= b_A[(i + kAcol) - 1] * b_B[tmp];
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
void SupervisoryController::KWIKfactor(const real_T b_Ac[984], const int32_T iC
  [246], int32_T nA, const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
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
          b_D[b_coltop] = 0.0;
          scalarLB = (k_i + 1) << 2UL;
          for (b_lastv = k_i + 1; b_lastv <= nA; b_lastv++) {
            b_D[b_coltop] += TL[(scalarLB + ii) - 4] * RLinv[(scalarLB + k_i) -
              4];
            scalarLB += 4;
          }
        }

        knt += 4;
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
void SupervisoryController::qpkwik(const real_T b_Linv[16], const real_T b_Hinv
  [16], const real_T f[4], const real_T b_Ac[984], const real_T b[246],
  boolean_T iA[246], int32_T maxiter, real_T FeasTol, real_T x[4], real_T
  lambda[246], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[246];
  real_T RLinv[16];
  real_T U[16];
  real_T b_D[16];
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
        KWIKfactor(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees, RLinv, &Xnorm0);
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
              Opt[i] += b_D[iC_0 + i] * Rhs[iSave + 4];
              iC_0 += 4;
            }
          }

          U_tmp = 0;
          for (i = 0; i < nA; i++) {
            Opt[i + 4] = ((b_D[U_tmp + 1] * Rhs[1] + b_D[U_tmp] * Rhs[0]) +
                          b_D[U_tmp + 2] * Rhs[2]) + b_D[U_tmp + 3] * Rhs[3];
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
              KWIKfactor(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees, RLinv, &cMin);
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
                  r[i] = ((b_D[iSave + 1] * b_Ac[tmp + 246] + b_D[iSave] *
                           b_Ac[tmp]) + b_D[iSave + 2] * b_Ac[tmp + 492]) +
                    b_D[iSave + 3] * b_Ac[tmp + 738];
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

// Function for MATLAB Function: '<S41>/FixedHorizonOptimizer'
void SupervisoryController::mpcblock_optimizer(const real_T rseq[120], const
  real_T vseq[21], const real_T umax[3], const real_T ymin[6], const real_T
  ymax[6], const real_T x[18], const real_T old_u[3], const boolean_T iA[246],
  const real_T b_Mlim[246], real_T b_Mx[4428], real_T b_Mu1[738], real_T b_Mv
  [5166], const real_T b_utarget[60], const real_T b_uoff[3], const real_T
  b_yoff[6], real_T b_H[16], real_T b_Ac[984], const real_T ywt[6], const real_T
  uwt[3], const real_T b_Wdu[3], const real_T b_Jm[180], const real_T b_I1[180],
  const real_T b_A[324], const real_T Bu[1134], const real_T Bv[378], const
  real_T b_C[108], const real_T Dv[126], const int32_T b_Mrows[246], const
  real_T b_RYscale[6], const real_T b_RMVscale[3], real_T u[3], real_T useq[63],
  real_T *status, boolean_T iAout[246])
{
  static const int8_T c_A[400]{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
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

  __m128d tmp;
  real_T Sum_0[360];
  real_T WySuJm[360];
  real_T b_Su1[360];
  real_T Bc[246];
  real_T a__1[246];
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T CA_0[126];
  real_T CA[108];
  real_T CA_1[108];
  real_T b_Kv[63];
  real_T b_Kx[54];
  real_T Sum[18];
  real_T b_Linv[16];
  real_T c_Linv[16];
  real_T b_B[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T rows[6];
  real_T ymin_incr[6];
  real_T f[4];
  real_T zopt[4];
  real_T b_Wu[3];
  real_T Sum_1;
  int32_T CA_tmp;
  int32_T Sum_tmp;
  int32_T a_tmp;
  int32_T i;
  int32_T i1;
  int32_T j2;
  int32_T kidx;
  int16_T ixw;
  int8_T a[3600];
  int8_T c_B[16];
  int8_T b_Su1_tmp[6];
  boolean_T ymax_incr_flag[6];
  boolean_T ymin_incr_flag[6];
  boolean_T umax_incr_flag[3];
  boolean_T exitg1;
  (void)std::memset(&useq[0], 0, 63U * sizeof(real_T));
  (void)std::memset(&iAout[0], 0, 246U * sizeof(boolean_T));
  for (i = 0; i < 6; i++) {
    for (j2 = 0; j2 < 18; j2++) {
      CA_tmp = 6 * j2 + i;
      CA[CA_tmp] = 0.0;
      for (i1 = 0; i1 < 18; i1++) {
        CA[CA_tmp] += b_C[6 * i1 + i] * b_A[18 * j2 + i1];
      }
    }

    for (j2 = 0; j2 < 3; j2++) {
      Sum_tmp = 6 * j2 + i;
      Sum[Sum_tmp] = 0.0;
      for (i1 = 0; i1 < 18; i1++) {
        Sum[Sum_tmp] += b_C[6 * i1 + i] * Bu[18 * j2 + i1];
      }
    }

    rows[i] = 0.0;
    for (j2 = 0; j2 < 18; j2++) {
      rows[i] += b_C[6 * j2 + i] * Bv[j2];
    }

    rtDW.b_Hv[i] = rows[i];
    rtDW.b_Hv[i + 120] = Dv[i];
  }

  i = 0;
  for (j2 = 0; j2 < 19; j2++) {
    for (i1 = 0; i1 < 6; i1++) {
      rtDW.b_Hv[(i1 + i) + 240] = 0.0;
    }

    i += 120;
  }

  i = 0;
  for (j2 = 0; j2 < 21; j2++) {
    (void)std::memset(&rtDW.b_Hv[i + 6], 0, 114U * sizeof(real_T));
    i += 120;
  }

  for (i = 0; i < 18; i++) {
    for (j2 = 0; j2 < 6; j2++) {
      rtDW.b_Sx[j2 + 120 * i] = CA[6 * i + j2];
    }

    (void)std::memset(&rtDW.b_Sx[i * 120 + 6], 0, 114U * sizeof(real_T));
  }

  for (i = 0; i < 6; i++) {
    b_Su1[i] = Sum[i];
    b_Su1[i + 120] = Sum[i + 6];
    b_Su1[i + 240] = Sum[i + 12];
  }

  for (i = 0; i < 114; i++) {
    b_Su1[i + 6] = 0.0;
    b_Su1[i + 126] = 0.0;
    b_Su1[i + 246] = 0.0;
  }

  for (i = 0; i < 6; i++) {
    rtDW.Su[i] = Sum[i];
    rtDW.Su[i + 120] = Sum[i + 6];
    rtDW.Su[i + 240] = Sum[i + 12];
  }

  i = 0;
  for (j2 = 0; j2 < 57; j2++) {
    for (i1 = 0; i1 < 6; i1++) {
      rtDW.Su[(i1 + i) + 360] = 0.0;
    }

    i += 120;
  }

  i = 0;
  for (j2 = 0; j2 < 60; j2++) {
    (void)std::memset(&rtDW.Su[i + 6], 0, 114U * sizeof(real_T));
    i += 120;
  }

  for (kidx = 0; kidx < 19; kidx++) {
    CA_tmp = (kidx + 1) * 6;
    for (i = 0; i < 6; i++) {
      rows[i] = static_cast<real_T>(static_cast<int32_T>(CA_tmp + i)) + 1.0;
      j2 = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 3; a_tmp++) {
        Sum_1 = 0.0;
        Sum_tmp = 0;
        for (int32_T i_0{0}; i_0 < 18; i_0++) {
          Sum_1 += CA[Sum_tmp + i] * Bu[i_0 + i1];
          Sum_tmp += 6;
        }

        Sum_tmp = j2 + i;
        Sum[Sum_tmp] += Sum_1;
        j2 += 6;
        i1 += 18;
      }

      b_Su1_tmp[i] = static_cast<int8_T>(rows[i]);
    }

    for (i = 0; i < 3; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        Sum_tmp = 6 * i + j2;
        Sum_1 = Sum[Sum_tmp];
        b_Su1[(static_cast<int32_T>(b_Su1_tmp[j2]) + 120 * i) - 1] = Sum_1;
        Sum_0[Sum_tmp] = Sum_1;
      }
    }

    for (i = 0; i < 57; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        Sum_0[j2 + 6 * (i + 3)] = rtDW.Su[(120 * i + static_cast<int32_T>
          (b_Su1_tmp[j2])) - 7];
      }
    }

    for (i = 0; i < 60; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        rtDW.Su[(static_cast<int32_T>(rows[j2]) + 120 * i) - 1] = Sum_0[6 * i +
          j2];
      }
    }

    for (i = 0; i < 6; i++) {
      ymin_incr[i] = 0.0;
      j2 = 0;
      for (i1 = 0; i1 < 18; i1++) {
        ymin_incr[i] += CA[j2 + i] * Bv[i1];
        j2 += 6;
      }

      CA_0[i] = ymin_incr[i];
    }

    for (i = 0; i < 20; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        CA_0[j2 + 6 * (i + 1)] = rtDW.b_Hv[(120 * i + static_cast<int32_T>
          (rows[j2])) - 7];
      }
    }

    for (i = 0; i < 21; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        rtDW.b_Hv[(static_cast<int32_T>(rows[j2]) + 120 * i) - 1] = CA_0[6 * i +
          j2];
      }
    }

    for (i = 0; i < 6; i++) {
      j2 = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 18; a_tmp++) {
        CA_tmp = j2 + i;
        CA_1[CA_tmp] = 0.0;
        Sum_tmp = 0;
        for (int32_T i_0{0}; i_0 < 18; i_0++) {
          CA_1[CA_tmp] += CA[Sum_tmp + i] * b_A[i_0 + i1];
          Sum_tmp += 6;
        }

        j2 += 6;
        i1 += 18;
      }
    }

    (void)std::memcpy(&CA[0], &CA_1[0], 108U * sizeof(real_T));
    for (i = 0; i < 18; i++) {
      for (j2 = 0; j2 < 6; j2++) {
        rtDW.b_Sx[(static_cast<int32_T>(rows[j2]) + 120 * i) - 1] = CA[6 * i +
          j2];
      }
    }
  }

  i = 0;
  j2 = 0;
  for (i1 = 0; i1 < 3; i1++) {
    for (a_tmp = 0; a_tmp < 120; a_tmp++) {
      kidx = a_tmp + i;
      Sum_0[kidx] = 0.0;
      Sum_tmp = 0;
      for (int32_T i_0{0}; i_0 < 60; i_0++) {
        Sum_0[kidx] += rtDW.Su[Sum_tmp + a_tmp] * b_Jm[i_0 + j2];
        Sum_tmp += 120;
      }
    }

    i += 120;
    j2 += 60;
  }

  if (b_Mrows[0] > 0) {
    kidx = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (kidx < 246)) {
      if (b_Mrows[kidx] <= 120) {
        j2 = b_Mrows[kidx];
        b_Ac[kidx] = -Sum_0[j2 - 1];
        b_Ac[kidx + 246] = -Sum_0[j2 + 119];
        b_Ac[kidx + 492] = -Sum_0[j2 + 239];
        j2 = b_Mrows[kidx];
        for (i = 0; i < 18; i++) {
          b_Mx[kidx + 246 * i] = -rtDW.b_Sx[(120 * i + j2) - 1];
        }

        j2 = b_Mrows[kidx];
        b_Mu1[kidx] = -b_Su1[j2 - 1];
        b_Mu1[kidx + 246] = -b_Su1[j2 + 119];
        b_Mu1[kidx + 492] = -b_Su1[j2 + 239];
        j2 = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 246 * i] = -rtDW.b_Hv[(120 * i + j2) - 1];
        }

        kidx++;
      } else if (b_Mrows[kidx] <= 240) {
        j2 = b_Mrows[kidx];
        b_Ac[kidx] = Sum_0[j2 - 121];
        b_Ac[kidx + 246] = Sum_0[j2 - 1];
        b_Ac[kidx + 492] = Sum_0[j2 + 119];
        j2 = b_Mrows[kidx];
        for (i = 0; i < 18; i++) {
          b_Mx[kidx + 246 * i] = rtDW.b_Sx[(120 * i + j2) - 121];
        }

        j2 = b_Mrows[kidx];
        b_Mu1[kidx] = b_Su1[j2 - 121];
        b_Mu1[kidx + 246] = b_Su1[j2 - 1];
        b_Mu1[kidx + 492] = b_Su1[j2 + 119];
        j2 = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 246 * i] = rtDW.b_Hv[(120 * i + j2) - 121];
        }

        kidx++;
      } else {
        exitg1 = true;
      }
    }
  }

  for (kidx = 0; kidx < 6; kidx++) {
    Sum_1 = ywt[kidx];
    if (Sum_1 < 0.0) {
      rows[kidx] = 0.0;
    } else {
      rows[kidx] = Sum_1 * Sum_1;
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

  (void)std::memset(&b_B[0], 0, 9U * sizeof(real_T));
  b_B[0] = 1.0;
  b_B[4] = 1.0;
  b_B[8] = 1.0;
  kidx = -1;
  for (i = 0; i < 20; i++) {
    for (j2 = 0; j2 < 3; j2++) {
      for (i1 = 0; i1 < 20; i1++) {
        a_tmp = static_cast<int32_T>(c_A[20 * i + i1]);
        a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * j2]) *
          a_tmp);
        a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * j2 + 1]) *
          a_tmp);
        a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * j2 + 2]) *
          a_tmp);
        kidx += 3;
      }
    }
  }

  i = 0;
  for (j2 = 0; j2 < 3; j2++) {
    for (i1 = 0; i1 < 60; i1++) {
      kidx = i1 + i;
      I2Jm[kidx] = 0.0;
      a_tmp = 0;
      for (Sum_tmp = 0; Sum_tmp < 60; Sum_tmp++) {
        I2Jm[kidx] += static_cast<real_T>(a[a_tmp + i1]) * b_Jm[Sum_tmp + i];
        a_tmp += 60;
      }
    }

    i += 60;
  }

  ixw = 1;
  for (kidx = 0; kidx < 120; kidx++) {
    Sum_1 = rows[ixw - 1];
    WySuJm[kidx] = Sum_1 * Sum_0[kidx];
    WySuJm[kidx + 120] = Sum_0[kidx + 120] * Sum_1;
    WySuJm[kidx + 240] = Sum_0[kidx + 240] * Sum_1;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 6) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    Sum_1 = b_Wu[ixw - 1];
    WuI2Jm[kidx] = Sum_1 * I2Jm[kidx];
    WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * Sum_1;
    WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * Sum_1;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    Sum_1 = b_Wdu[ixw - 1];
    WduJm[kidx] = Sum_1 * b_Jm[kidx];
    WduJm[kidx + 60] = b_Jm[kidx + 60] * Sum_1;
    WduJm[kidx + 120] = b_Jm[kidx + 120] * Sum_1;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j2 = 0; j2 < 3; j2++) {
      kidx = 3 * j2 + i;
      b_B[kidx] = 0.0;
      for (i1 = 0; i1 < 120; i1++) {
        b_B[kidx] += Sum_0[120 * i + i1] * WySuJm[120 * j2 + i1];
      }

      b_Jm_0[kidx] = 0.0;
      Sum_1 = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        a_tmp = 60 * i + i1;
        Sum_tmp = 60 * j2 + i1;
        Sum_1 += I2Jm[a_tmp] * WuI2Jm[Sum_tmp];
        b_Jm_0[kidx] += b_Jm[a_tmp] * WduJm[Sum_tmp];
      }

      b_H[i + (j2 << 2UL)] = (b_B[kidx] + b_Jm_0[kidx]) + Sum_1;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j2 = 0; j2 < 3; j2++) {
      kidx = 3 * j2 + i;
      b_Jm_0[kidx] = 0.0;
      for (i1 = 0; i1 < 120; i1++) {
        b_Jm_0[kidx] += b_Su1[120 * i + i1] * WySuJm[120 * j2 + i1];
      }

      b_I1_0[kidx] = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        b_I1_0[kidx] += b_I1[60 * i + i1] * WuI2Jm[60 * j2 + i1];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
    __m128d tmp_0;
    tmp = _mm_loadu_pd(&b_Jm_0[i]);
    tmp_0 = _mm_loadu_pd(&b_I1_0[i]);
    (void)_mm_storeu_pd(&b_B[i], _mm_add_pd(tmp, tmp_0));
  }

  for (i = 8; i < 9; i++) {
    b_B[i] = b_Jm_0[i] + b_I1_0[i];
  }

  for (i = 0; i <= 178; i += 2) {
    tmp = _mm_loadu_pd(&WuI2Jm[i]);
    (void)_mm_storeu_pd(&WuI2Jm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  i = 0;
  for (j2 = 0; j2 < 18; j2++) {
    i1 = 0;
    a_tmp = 0;
    for (Sum_tmp = 0; Sum_tmp < 3; Sum_tmp++) {
      kidx = i1 + j2;
      b_Kx[kidx] = 0.0;
      for (int32_T i_0{0}; i_0 < 120; i_0++) {
        b_Kx[kidx] += rtDW.b_Sx[i_0 + i] * WySuJm[i_0 + a_tmp];
      }

      i1 += 18;
      a_tmp += 120;
    }

    i += 120;
  }

  i = 0;
  for (j2 = 0; j2 < 21; j2++) {
    i1 = 0;
    a_tmp = 0;
    for (Sum_tmp = 0; Sum_tmp < 3; Sum_tmp++) {
      kidx = i1 + j2;
      b_Kv[kidx] = 0.0;
      for (int32_T i_0{0}; i_0 < 120; i_0++) {
        b_Kv[kidx] += rtDW.b_Hv[i_0 + i] * WySuJm[i_0 + a_tmp];
      }

      i1 += 21;
      a_tmp += 120;
    }

    i += 120;
  }

  for (i = 0; i <= 358; i += 2) {
    tmp = _mm_loadu_pd(&WySuJm[i]);
    (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  (void)std::memcpy(&b_Linv[0], &b_H[0], sizeof(real_T) << 4UL);
  mpc_checkhessian(b_Linv, c_Linv, &Sum_1);
  if (Sum_1 > 1.0) {
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
    for (i = 0; i < 16; i++) {
      c_B[i] = 0;
    }

    c_B[0] = 1;
    c_B[5] = 1;
    c_B[10] = 1;
    c_B[15] = 1;
    CA_tmp = 0;
    for (kidx = 0; kidx < 4; kidx++) {
      b_Linv[CA_tmp] = static_cast<real_T>(c_B[CA_tmp]);
      b_Linv[CA_tmp + 1] = static_cast<real_T>(c_B[CA_tmp + 1]);
      b_Linv[CA_tmp + 2] = static_cast<real_T>(c_B[CA_tmp + 2]);
      b_Linv[CA_tmp + 3] = static_cast<real_T>(c_B[CA_tmp + 3]);
      CA_tmp += 4;
    }

    trisolve(c_Linv, b_Linv);
    for (i = 0; i < 246; i++) {
      Sum_1 = 0.0;
      for (j2 = 0; j2 < 18; j2++) {
        Sum_1 += b_Mx[246 * j2 + i] * x[j2];
      }

      WySuJm_0 = 0.0;
      for (j2 = 0; j2 < 21; j2++) {
        WySuJm_0 += b_Mv[246 * j2 + i] * vseq[j2];
      }

      Bc[i] = -((((b_Mu1[i + 246] * old_u[1] + b_Mu1[i] * old_u[0]) + b_Mu1[i +
                  492] * old_u[2]) + (b_Mlim[i] + Sum_1)) + WySuJm_0);
    }

    for (i = 0; i < 6; i++) {
      ymax_incr_flag[i] = false;
      rows[i] = 0.0;
      ymin_incr_flag[i] = false;
      ymin_incr[i] = 0.0;
    }

    umax_incr_flag[0] = false;
    b_Wu[0] = 0.0;
    umax_incr_flag[1] = false;
    b_Wu[1] = 0.0;
    umax_incr_flag[2] = false;
    b_Wu[2] = 0.0;
    if (b_Mrows[0] > 0) {
      kidx = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (kidx < 246)) {
        if (b_Mrows[kidx] <= 120) {
          boolean_T c_Del_Save_Flag0;
          i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 1, static_cast<
                int32_T>(ny)) * static_cast<int32_T>(ny)) - 1;
          c_Del_Save_Flag0 = ymax_incr_flag[i];
          if (!ymax_incr_flag[i]) {
            Sum_1 = -(b_RYscale[i] * ymax[i] - b_yoff[i]) - (-b_Mlim[kidx]);
            c_Del_Save_Flag0 = true;
          } else {
            Sum_1 = rows[i];
          }

          rows[i] = Sum_1;
          ymax_incr_flag[i] = c_Del_Save_Flag0;
          Bc[kidx] += Sum_1;
          kidx++;
        } else if (b_Mrows[kidx] <= 240) {
          boolean_T c_Del_Save_Flag0;
          i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 121,
                static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 121;
          c_Del_Save_Flag0 = ymin_incr_flag[i];
          if (!ymin_incr_flag[i]) {
            Sum_1 = (b_RYscale[i] * ymin[i] - b_yoff[i]) - (-b_Mlim[kidx]);
            c_Del_Save_Flag0 = true;
          } else {
            Sum_1 = ymin_incr[i];
          }

          ymin_incr[i] = Sum_1;
          ymin_incr_flag[i] = c_Del_Save_Flag0;
          Bc[kidx] += Sum_1;
          kidx++;
        } else if (b_Mrows[kidx] <= 300) {
          boolean_T c_Del_Save_Flag0;
          i = (b_Mrows[kidx] - div_nde_s32_floor(b_Mrows[kidx] - 241,
                static_cast<int32_T>(nu)) * static_cast<int32_T>(nu)) - 241;
          c_Del_Save_Flag0 = umax_incr_flag[i];
          if (!umax_incr_flag[i]) {
            Sum_1 = -(b_RMVscale[i] * umax[i] - b_uoff[i]) - (-b_Mlim[kidx]);
            c_Del_Save_Flag0 = true;
          } else {
            Sum_1 = b_Wu[i];
          }

          b_Wu[i] = Sum_1;
          umax_incr_flag[i] = c_Del_Save_Flag0;
          Bc[kidx] += Sum_1;
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
      real_T b_Kv_0;
      Sum_1 = 0.0;
      for (i = 0; i < 18; i++) {
        Sum_1 += b_Kx[18 * kidx + i] * x[i];
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

      f[kidx] = ((((b_B[3 * kidx + 1] * old_u[1] + b_B[3 * kidx] * old_u[0]) +
                   b_B[3 * kidx + 2] * old_u[2]) + (Sum_1 + WySuJm_0)) + b_Kv_0)
        + WuI2Jm_0;
    }

    (void)std::memcpy(&iAout[0], &iA[0], 246U * sizeof(boolean_T));
    i = 0;
    for (j2 = 0; j2 < 4; j2++) {
      i1 = 0;
      for (a_tmp = 0; a_tmp < 4; a_tmp++) {
        kidx = i1 + j2;
        c_Linv[kidx] = 0.0;
        c_Linv[kidx] += b_Linv[i] * b_Linv[i1];
        c_Linv[kidx] += b_Linv[i + 1] * b_Linv[i1 + 1];
        c_Linv[kidx] += b_Linv[i + 2] * b_Linv[i1 + 2];
        c_Linv[kidx] += b_Linv[i + 3] * b_Linv[i1 + 3];
        i1 += 4;
      }

      i += 4;
    }

    qpkwik(b_Linv, c_Linv, f, b_Ac, Bc, iAout, 1000, 1.0E-6, zopt, a__1, &kidx);
    if ((kidx < 0) || (kidx == 0)) {
      zopt[0] = 0.0;
      zopt[1] = 0.0;
      zopt[2] = 0.0;
    }

    *status = static_cast<real_T>(kidx);
    u[0] = (old_u[0] + zopt[0]) + b_uoff[0];
    u[1] = (old_u[1] + zopt[1]) + b_uoff[1];
    u[2] = (old_u[2] + zopt[2]) + b_uoff[2];
  }
}

// Function for MATLAB Function: '<S45>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv_g(const real_T A[108], const real_T B_1[36],
  real_T Y[108])
{
  __m128d tmp;
  real_T b_A[36];
  real_T smax;
  int32_T b_ix;
  int32_T ijA;
  int32_T ix;
  int32_T iy;
  int32_T jj;
  int32_T kBcol;
  int8_T ipiv[6];
  (void)std::memcpy(&b_A[0], &B_1[0], 36U * sizeof(real_T));
  (void)std::memcpy(&Y[0], &A[0], 108U * sizeof(real_T));
  for (int32_T d_j{0}; d_j < 6; d_j++) {
    ipiv[d_j] = static_cast<int8_T>(d_j + 1);
  }

  for (int32_T d_j{0}; d_j < 5; d_j++) {
    jj = d_j * 7;
    iy = 6 - d_j;
    b_ix = 0;
    smax = std::abs(b_A[jj]);
    for (kBcol = 2; kBcol <= iy; kBcol++) {
      real_T s;
      s = std::abs(b_A[(jj + kBcol) - 1]);
      if (s > smax) {
        b_ix = kBcol - 1;
        smax = s;
      }
    }

    if (b_A[jj + b_ix] != 0.0) {
      if (b_ix != 0) {
        iy = d_j + b_ix;
        ipiv[d_j] = static_cast<int8_T>(iy + 1);
        for (ix = 0; ix < 6; ix++) {
          b_ix = ix * 6 + d_j;
          smax = b_A[b_ix];
          b_A[b_ix] = b_A[iy];
          b_A[iy] = smax;
          iy += 6;
        }
      }

      iy = (jj - d_j) + 6;
      kBcol = (((((iy - jj) - 1) / 2) << 1UL) + jj) + 2;
      ix = kBcol - 2;
      for (b_ix = jj + 2; b_ix <= ix; b_ix += 2) {
        tmp = _mm_loadu_pd(&b_A[b_ix - 1]);
        (void)_mm_storeu_pd(&b_A[b_ix - 1], _mm_div_pd(tmp, _mm_set1_pd(b_A[jj])));
      }

      for (b_ix = kBcol; b_ix <= iy; b_ix++) {
        b_A[b_ix - 1] /= b_A[jj];
      }
    }

    iy = 4 - d_j;
    b_ix = jj + 8;
    for (kBcol = 0; kBcol <= iy; kBcol++) {
      smax = b_A[(kBcol * 6 + jj) + 6];
      if (smax != 0.0) {
        ix = (b_ix - d_j) + 4;
        for (ijA = b_ix; ijA <= ix; ijA++) {
          b_A[ijA - 1] += b_A[((jj + ijA) - b_ix) + 1] * -smax;
        }
      }

      b_ix += 6;
    }
  }

  for (int32_T d_j{0}; d_j < 6; d_j++) {
    jj = 18 * d_j;
    iy = 6 * d_j;
    for (ix = 0; ix < d_j; ix++) {
      kBcol = 18 * ix;
      if (b_A[ix + iy] != 0.0) {
        for (b_ix = 0; b_ix < 18; b_ix++) {
          ijA = b_ix + jj;
          Y[ijA] -= b_A[ix + iy] * Y[b_ix + kBcol];
        }
      }
    }

    smax = 1.0 / b_A[d_j + iy];
    for (iy = 0; iy <= 16; iy += 2) {
      tmp = _mm_loadu_pd(&Y[iy + jj]);
      (void)_mm_storeu_pd(&Y[iy + jj], _mm_mul_pd(tmp, _mm_set1_pd(smax)));
    }
  }

  for (int32_T d_j{5}; d_j >= 0; d_j--) {
    jj = 18 * d_j;
    iy = 6 * d_j - 1;
    for (b_ix = d_j + 2; b_ix < 7; b_ix++) {
      ix = (b_ix - 1) * 18;
      if (b_A[b_ix + iy] != 0.0) {
        for (kBcol = 0; kBcol < 18; kBcol++) {
          ijA = kBcol + jj;
          Y[ijA] -= b_A[b_ix + iy] * Y[kBcol + ix];
        }
      }
    }
  }

  for (int32_T d_j{4}; d_j >= 0; d_j--) {
    int8_T ipiv_0;
    ipiv_0 = ipiv[d_j];
    if (d_j + 1 != static_cast<int32_T>(ipiv_0)) {
      for (iy = 0; iy < 18; iy++) {
        b_ix = 18 * d_j + iy;
        smax = Y[b_ix];
        ijA = (static_cast<int32_T>(ipiv_0) - 1) * 18 + iy;
        Y[b_ix] = Y[ijA];
        Y[ijA] = smax;
      }
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
  static const real_T f[4428]{ -0.975, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.95062499999999994, -0.0, -0.0, -0.0, -0.0, -0.0, -0.92685937499999993,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.90368789062499988, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.88109569335937488, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.85906830102539045, -0.0, -0.0, -0.0, -0.0, -0.0, -0.83759159349975565,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.81665180366226175, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.79623550857070524, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.77632962085643764, -0.0, -0.0, -0.0, -0.0, -0.0, -0.75692138033502665,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.73799834582665091, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.71954838718098457, -0.0, -0.0, -0.0, -0.0, -0.0, -0.70155967750146,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.68402068556392348, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.66692016842482538, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.6502471642142047, -0.0, -0.0, -0.0, -0.0, -0.0, -0.63399098510884955,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.61814121048112824, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.6026876802191, -0.0, -0.0, -0.0, -0.0, -0.0, 0.975, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.95062499999999994, 0.0, 0.0, 0.0, 0.0, 0.0, 0.92685937499999993,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.90368789062499988, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.88109569335937488, 0.0, 0.0, 0.0, 0.0, 0.0, 0.85906830102539045, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.83759159349975565, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.81665180366226175, 0.0, 0.0, 0.0, 0.0, 0.0, 0.79623550857070524, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.77632962085643764, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.75692138033502665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.73799834582665091, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.71954838718098457, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.70155967750146, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68402068556392348, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.66692016842482538, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.6502471642142047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.63399098510884955, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.61814121048112824, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6026876802191,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.975, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.95062499999999994, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.92685937499999993, -0.0, -0.0, -0.0, -0.0, -0.0, -0.90368789062499988,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.88109569335937488, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.85906830102539045, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.83759159349975565, -0.0, -0.0, -0.0, -0.0, -0.0, -0.81665180366226175,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.79623550857070524, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.77632962085643764, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.75692138033502665, -0.0, -0.0, -0.0, -0.0, -0.0, -0.73799834582665091,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.71954838718098457, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.70155967750146, -0.0, -0.0, -0.0, -0.0, -0.0, -0.68402068556392348,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.66692016842482538, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.6502471642142047, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.63399098510884955, -0.0, -0.0, -0.0, -0.0, -0.0, -0.61814121048112824,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.6026876802191, -0.0, -0.0, -0.0, -0.0, 0.0,
    0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95062499999999994, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.92685937499999993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.90368789062499988, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.88109569335937488, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.85906830102539045, 0.0, 0.0, 0.0, 0.0, 0.0, 0.83759159349975565, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.81665180366226175, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.79623550857070524, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77632962085643764, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.75692138033502665, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.73799834582665091, 0.0, 0.0, 0.0, 0.0, 0.0, 0.71954838718098457, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.70155967750146, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.68402068556392348, 0.0, 0.0, 0.0, 0.0, 0.0, 0.66692016842482538, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.6502471642142047, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.63399098510884955, 0.0, 0.0, 0.0, 0.0, 0.0, 0.61814121048112824, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.6026876802191, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.0, -0.0, -0.975, -0.0, -0.0, -0.0, -0.0, -0.0, -0.95062499999999994,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.92685937499999993, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.90368789062499988, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.88109569335937488, -0.0, -0.0, -0.0, -0.0, -0.0, -0.85906830102539045,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.83759159349975565, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.81665180366226175, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.79623550857070524, -0.0, -0.0, -0.0, -0.0, -0.0, -0.77632962085643764,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.75692138033502665, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.73799834582665091, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.71954838718098457, -0.0, -0.0, -0.0, -0.0, -0.0, -0.70155967750146, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.68402068556392348, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.66692016842482538, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6502471642142047,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.63399098510884955, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.61814121048112824, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6026876802191,
    -0.0, -0.0, -0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.95062499999999994, 0.0, 0.0, 0.0, 0.0, 0.0, 0.92685937499999993, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.90368789062499988, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.88109569335937488, 0.0, 0.0, 0.0, 0.0, 0.0, 0.85906830102539045, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.83759159349975565, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.81665180366226175, 0.0, 0.0, 0.0, 0.0, 0.0, 0.79623550857070524, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.77632962085643764, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.75692138033502665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.73799834582665091, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.71954838718098457, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.70155967750146, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68402068556392348, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.66692016842482538, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.6502471642142047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.63399098510884955, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.61814121048112824, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6026876802191,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.975, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.95062499999999994, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.92685937499999993, -0.0, -0.0, -0.0, -0.0, -0.0, -0.90368789062499988,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.88109569335937488, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.85906830102539045, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.83759159349975565, -0.0, -0.0, -0.0, -0.0, -0.0, -0.81665180366226175,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.79623550857070524, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.77632962085643764, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.75692138033502665, -0.0, -0.0, -0.0, -0.0, -0.0, -0.73799834582665091,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.71954838718098457, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.70155967750146, -0.0, -0.0, -0.0, -0.0, -0.0, -0.68402068556392348,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.66692016842482538, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.6502471642142047, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.63399098510884955, -0.0, -0.0, -0.0, -0.0, -0.0, -0.61814121048112824,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.6026876802191, -0.0, -0.0, 0.0, 0.0, 0.0,
    0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95062499999999994, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.92685937499999993, 0.0, 0.0, 0.0, 0.0, 0.0, 0.90368789062499988, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.88109569335937488, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.85906830102539045, 0.0, 0.0, 0.0, 0.0, 0.0, 0.83759159349975565, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.81665180366226175, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.79623550857070524, 0.0, 0.0, 0.0, 0.0, 0.0, 0.77632962085643764, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.75692138033502665, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.73799834582665091, 0.0, 0.0, 0.0, 0.0, 0.0, 0.71954838718098457, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.70155967750146, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.68402068556392348, 0.0, 0.0, 0.0, 0.0, 0.0, 0.66692016842482538, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.6502471642142047, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.63399098510884955, 0.0, 0.0, 0.0, 0.0, 0.0, 0.61814121048112824, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.6026876802191, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0,
    -0.0, -0.0, -0.0, -0.975, -0.0, -0.0, -0.0, -0.0, -0.0, -0.95062499999999994,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.92685937499999993, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.90368789062499988, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.88109569335937488, -0.0, -0.0, -0.0, -0.0, -0.0, -0.85906830102539045,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.83759159349975565, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.81665180366226175, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.79623550857070524, -0.0, -0.0, -0.0, -0.0, -0.0, -0.77632962085643764,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.75692138033502665, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.73799834582665091, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.71954838718098457, -0.0, -0.0, -0.0, -0.0, -0.0, -0.70155967750146, -0.0,
    -0.0, -0.0, -0.0, -0.0, -0.68402068556392348, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.66692016842482538, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6502471642142047,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.63399098510884955, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.61814121048112824, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6026876802191,
    -0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.95062499999999994, 0.0, 0.0, 0.0, 0.0, 0.0, 0.92685937499999993, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.90368789062499988, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.88109569335937488, 0.0, 0.0, 0.0, 0.0, 0.0, 0.85906830102539045, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.83759159349975565, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.81665180366226175, 0.0, 0.0, 0.0, 0.0, 0.0, 0.79623550857070524, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.77632962085643764, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.75692138033502665, 0.0, 0.0, 0.0, 0.0, 0.0, 0.73799834582665091, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.71954838718098457, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.70155967750146, 0.0, 0.0, 0.0, 0.0, 0.0, 0.68402068556392348, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.66692016842482538, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.6502471642142047, 0.0, 0.0, 0.0, 0.0, 0.0, 0.63399098510884955, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.61814121048112824, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6026876802191,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, -0.975,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.95062499999999994, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.92685937499999993, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.90368789062499988, -0.0, -0.0, -0.0, -0.0, -0.0, -0.88109569335937488,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.85906830102539045, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.83759159349975565, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.81665180366226175, -0.0, -0.0, -0.0, -0.0, -0.0, -0.79623550857070524,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.77632962085643764, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.75692138033502665, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.73799834582665091, -0.0, -0.0, -0.0, -0.0, -0.0, -0.71954838718098457,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.70155967750146, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.68402068556392348, -0.0, -0.0, -0.0, -0.0, -0.0,
    -0.66692016842482538, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6502471642142047,
    -0.0, -0.0, -0.0, -0.0, -0.0, -0.63399098510884955, -0.0, -0.0, -0.0, -0.0,
    -0.0, -0.61814121048112824, -0.0, -0.0, -0.0, -0.0, -0.0, -0.6026876802191,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95062499999999994,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.92685937499999993, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.90368789062499988, 0.0, 0.0, 0.0, 0.0, 0.0, 0.88109569335937488, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.85906830102539045, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.83759159349975565, 0.0, 0.0, 0.0, 0.0, 0.0, 0.81665180366226175, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.79623550857070524, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.77632962085643764, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75692138033502665, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.73799834582665091, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.71954838718098457, 0.0, 0.0, 0.0, 0.0, 0.0, 0.70155967750146, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.68402068556392348, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.66692016842482538, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6502471642142047, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.63399098510884955, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.61814121048112824, 0.0, 0.0, 0.0, 0.0, 0.0, 0.6026876802191, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, -0.0125, 0.0125, 0.0125, -0.0, -0.0, -0.0, -0.025, 0.025,
    0.025, -0.0, -0.0, -0.0, -0.037500000000000006, 0.037500000000000006,
    0.037500000000000006, -0.0, -0.0, -0.0, -0.05, 0.05, 0.05, -0.0, -0.0, -0.0,
    -0.0625, 0.0625, 0.0625, -0.0, -0.0, -0.0, -0.075, 0.075, 0.075, -0.0, -0.0,
    -0.0, -0.0875, 0.0875, 0.0875, -0.0, -0.0, -0.0, -0.099999999999999992,
    0.099999999999999992, 0.099999999999999992, -0.0, -0.0, -0.0,
    -0.11249999999999999, 0.11249999999999999, 0.11249999999999999, -0.0, -0.0,
    -0.0, -0.12499999999999999, 0.12499999999999999, 0.12499999999999999, -0.0,
    -0.0, -0.0, -0.13749999999999998, 0.13749999999999998, 0.13749999999999998,
    -0.0, -0.0, -0.0, -0.15, 0.15, 0.15, -0.0, -0.0, -0.0, -0.1625, 0.1625,
    0.1625, -0.0, -0.0, -0.0, -0.17500000000000002, 0.17500000000000002,
    0.17500000000000002, -0.0, -0.0, -0.0, -0.18750000000000003,
    0.18750000000000003, 0.18750000000000003, -0.0, -0.0, -0.0,
    -0.20000000000000004, 0.20000000000000004, 0.20000000000000004, -0.0, -0.0,
    -0.0, -0.21250000000000005, 0.21250000000000005, 0.21250000000000005, -0.0,
    -0.0, -0.0, -0.22500000000000006, 0.22500000000000006, 0.22500000000000006,
    -0.0, -0.0, -0.0, -0.23750000000000007, 0.23750000000000007,
    0.23750000000000007, -0.0, -0.0, -0.0, -0.25000000000000006,
    0.25000000000000006, 0.25000000000000006, -0.0, -0.0, -0.0, 0.0125, -0.0125,
    -0.0125, 0.0, 0.0, 0.0, 0.025, -0.025, -0.025, 0.0, 0.0, 0.0,
    0.037500000000000006, -0.037500000000000006, -0.037500000000000006, 0.0, 0.0,
    0.0, 0.05, -0.05, -0.05, 0.0, 0.0, 0.0, 0.0625, -0.0625, -0.0625, 0.0, 0.0,
    0.0, 0.075, -0.075, -0.075, 0.0, 0.0, 0.0, 0.0875, -0.0875, -0.0875, 0.0,
    0.0, 0.0, 0.099999999999999992, -0.099999999999999992, -0.099999999999999992,
    0.0, 0.0, 0.0, 0.11249999999999999, -0.11249999999999999,
    -0.11249999999999999, 0.0, 0.0, 0.0, 0.12499999999999999,
    -0.12499999999999999, -0.12499999999999999, 0.0, 0.0, 0.0,
    0.13749999999999998, -0.13749999999999998, -0.13749999999999998, 0.0, 0.0,
    0.0, 0.15, -0.15, -0.15, 0.0, 0.0, 0.0, 0.1625, -0.1625, -0.1625, 0.0, 0.0,
    0.0, 0.17500000000000002, -0.17500000000000002, -0.17500000000000002, 0.0,
    0.0, 0.0, 0.18750000000000003, -0.18750000000000003, -0.18750000000000003,
    0.0, 0.0, 0.0, 0.20000000000000004, -0.20000000000000004,
    -0.20000000000000004, 0.0, 0.0, 0.0, 0.21250000000000005,
    -0.21250000000000005, -0.21250000000000005, 0.0, 0.0, 0.0,
    0.22500000000000006, -0.22500000000000006, -0.22500000000000006, 0.0, 0.0,
    0.0, 0.23750000000000007, -0.23750000000000007, -0.23750000000000007, 0.0,
    0.0, 0.0, 0.25000000000000006, -0.25000000000000006, -0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0,
    -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, -0.5, 0.0,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5,
    -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5,
    -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0,
    0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5,
    -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, -0.0125, 0.0125,
    -0.0, -0.0, -0.0, 0.025, -0.025, 0.025, -0.0, -0.0, -0.0,
    0.037500000000000006, -0.037500000000000006, 0.037500000000000006, -0.0,
    -0.0, -0.0, 0.05, -0.05, 0.05, -0.0, -0.0, -0.0, 0.0625, -0.0625, 0.0625,
    -0.0, -0.0, -0.0, 0.075, -0.075, 0.075, -0.0, -0.0, -0.0, 0.0875, -0.0875,
    0.0875, -0.0, -0.0, -0.0, 0.099999999999999992, -0.099999999999999992,
    0.099999999999999992, -0.0, -0.0, -0.0, 0.11249999999999999,
    -0.11249999999999999, 0.11249999999999999, -0.0, -0.0, -0.0,
    0.12499999999999999, -0.12499999999999999, 0.12499999999999999, -0.0, -0.0,
    -0.0, 0.13749999999999998, -0.13749999999999998, 0.13749999999999998, -0.0,
    -0.0, -0.0, 0.15, -0.15, 0.15, -0.0, -0.0, -0.0, 0.1625, -0.1625, 0.1625,
    -0.0, -0.0, -0.0, 0.17500000000000002, -0.17500000000000002,
    0.17500000000000002, -0.0, -0.0, -0.0, 0.18750000000000003,
    -0.18750000000000003, 0.18750000000000003, -0.0, -0.0, -0.0,
    0.20000000000000004, -0.20000000000000004, 0.20000000000000004, -0.0, -0.0,
    -0.0, 0.21250000000000005, -0.21250000000000005, 0.21250000000000005, -0.0,
    -0.0, -0.0, 0.22500000000000006, -0.22500000000000006, 0.22500000000000006,
    -0.0, -0.0, -0.0, 0.23750000000000007, -0.23750000000000007,
    0.23750000000000007, -0.0, -0.0, -0.0, 0.25000000000000006,
    -0.25000000000000006, 0.25000000000000006, -0.0, -0.0, -0.0, -0.0125, 0.0125,
    -0.0125, 0.0, 0.0, 0.0, -0.025, 0.025, -0.025, 0.0, 0.0, 0.0,
    -0.037500000000000006, 0.037500000000000006, -0.037500000000000006, 0.0, 0.0,
    0.0, -0.05, 0.05, -0.05, 0.0, 0.0, 0.0, -0.0625, 0.0625, -0.0625, 0.0, 0.0,
    0.0, -0.075, 0.075, -0.075, 0.0, 0.0, 0.0, -0.0875, 0.0875, -0.0875, 0.0,
    0.0, 0.0, -0.099999999999999992, 0.099999999999999992, -0.099999999999999992,
    0.0, 0.0, 0.0, -0.11249999999999999, 0.11249999999999999,
    -0.11249999999999999, 0.0, 0.0, 0.0, -0.12499999999999999,
    0.12499999999999999, -0.12499999999999999, 0.0, 0.0, 0.0,
    -0.13749999999999998, 0.13749999999999998, -0.13749999999999998, 0.0, 0.0,
    0.0, -0.15, 0.15, -0.15, 0.0, 0.0, 0.0, -0.1625, 0.1625, -0.1625, 0.0, 0.0,
    0.0, -0.17500000000000002, 0.17500000000000002, -0.17500000000000002, 0.0,
    0.0, 0.0, -0.18750000000000003, 0.18750000000000003, -0.18750000000000003,
    0.0, 0.0, 0.0, -0.20000000000000004, 0.20000000000000004,
    -0.20000000000000004, 0.0, 0.0, 0.0, -0.21250000000000005,
    0.21250000000000005, -0.21250000000000005, 0.0, 0.0, 0.0,
    -0.22500000000000006, 0.22500000000000006, -0.22500000000000006, 0.0, 0.0,
    0.0, -0.23750000000000007, 0.23750000000000007, -0.23750000000000007, 0.0,
    0.0, 0.0, -0.25000000000000006, 0.25000000000000006, -0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.5, -0.0, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0,
    0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5,
    0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5,
    -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0,
    -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0,
    0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0,
    0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0,
    0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5,
    -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0,
    0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0,
    -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, 0.0125,
    -0.0125, -0.0, -0.0, -0.0, 0.025, 0.025, -0.025, -0.0, -0.0, -0.0,
    0.037500000000000006, 0.037500000000000006, -0.037500000000000006, -0.0,
    -0.0, -0.0, 0.05, 0.05, -0.05, -0.0, -0.0, -0.0, 0.0625, 0.0625, -0.0625,
    -0.0, -0.0, -0.0, 0.075, 0.075, -0.075, -0.0, -0.0, -0.0, 0.0875, 0.0875,
    -0.0875, -0.0, -0.0, -0.0, 0.099999999999999992, 0.099999999999999992,
    -0.099999999999999992, -0.0, -0.0, -0.0, 0.11249999999999999,
    0.11249999999999999, -0.11249999999999999, -0.0, -0.0, -0.0,
    0.12499999999999999, 0.12499999999999999, -0.12499999999999999, -0.0, -0.0,
    -0.0, 0.13749999999999998, 0.13749999999999998, -0.13749999999999998, -0.0,
    -0.0, -0.0, 0.15, 0.15, -0.15, -0.0, -0.0, -0.0, 0.1625, 0.1625, -0.1625,
    -0.0, -0.0, -0.0, 0.17500000000000002, 0.17500000000000002,
    -0.17500000000000002, -0.0, -0.0, -0.0, 0.18750000000000003,
    0.18750000000000003, -0.18750000000000003, -0.0, -0.0, -0.0,
    0.20000000000000004, 0.20000000000000004, -0.20000000000000004, -0.0, -0.0,
    -0.0, 0.21250000000000005, 0.21250000000000005, -0.21250000000000005, -0.0,
    -0.0, -0.0, 0.22500000000000006, 0.22500000000000006, -0.22500000000000006,
    -0.0, -0.0, -0.0, 0.23750000000000007, 0.23750000000000007,
    -0.23750000000000007, -0.0, -0.0, -0.0, 0.25000000000000006,
    0.25000000000000006, -0.25000000000000006, -0.0, -0.0, -0.0, -0.0125,
    -0.0125, 0.0125, 0.0, 0.0, 0.0, -0.025, -0.025, 0.025, 0.0, 0.0, 0.0,
    -0.037500000000000006, -0.037500000000000006, 0.037500000000000006, 0.0, 0.0,
    0.0, -0.05, -0.05, 0.05, 0.0, 0.0, 0.0, -0.0625, -0.0625, 0.0625, 0.0, 0.0,
    0.0, -0.075, -0.075, 0.075, 0.0, 0.0, 0.0, -0.0875, -0.0875, 0.0875, 0.0,
    0.0, 0.0, -0.099999999999999992, -0.099999999999999992, 0.099999999999999992,
    0.0, 0.0, 0.0, -0.11249999999999999, -0.11249999999999999,
    0.11249999999999999, 0.0, 0.0, 0.0, -0.12499999999999999,
    -0.12499999999999999, 0.12499999999999999, 0.0, 0.0, 0.0,
    -0.13749999999999998, -0.13749999999999998, 0.13749999999999998, 0.0, 0.0,
    0.0, -0.15, -0.15, 0.15, 0.0, 0.0, 0.0, -0.1625, -0.1625, 0.1625, 0.0, 0.0,
    0.0, -0.17500000000000002, -0.17500000000000002, 0.17500000000000002, 0.0,
    0.0, 0.0, -0.18750000000000003, -0.18750000000000003, 0.18750000000000003,
    0.0, 0.0, 0.0, -0.20000000000000004, -0.20000000000000004,
    0.20000000000000004, 0.0, 0.0, 0.0, -0.21250000000000005,
    -0.21250000000000005, 0.21250000000000005, 0.0, 0.0, 0.0,
    -0.22500000000000006, -0.22500000000000006, 0.22500000000000006, 0.0, 0.0,
    0.0, -0.23750000000000007, -0.23750000000000007, 0.23750000000000007, 0.0,
    0.0, 0.0, -0.25000000000000006, -0.25000000000000006, 0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0,
    -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, -0.5, -0.5, 0.5, 0.0,
    0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5,
    0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5,
    0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0,
    -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5,
    0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.0125,
    0.0125, 0.0125, -0.0, -0.0, -0.0, -0.025, 0.025, 0.025, -0.0, -0.0, -0.0,
    -0.037500000000000006, 0.037500000000000006, 0.037500000000000006, -0.0,
    -0.0, -0.0, -0.05, 0.05, 0.05, -0.0, -0.0, -0.0, -0.0625, 0.0625, 0.0625,
    -0.0, -0.0, -0.0, -0.075, 0.075, 0.075, -0.0, -0.0, -0.0, -0.0875, 0.0875,
    0.0875, -0.0, -0.0, -0.0, -0.099999999999999992, 0.099999999999999992,
    0.099999999999999992, -0.0, -0.0, -0.0, -0.11249999999999999,
    0.11249999999999999, 0.11249999999999999, -0.0, -0.0, -0.0,
    -0.12499999999999999, 0.12499999999999999, 0.12499999999999999, -0.0, -0.0,
    -0.0, -0.13749999999999998, 0.13749999999999998, 0.13749999999999998, -0.0,
    -0.0, -0.0, -0.15, 0.15, 0.15, -0.0, -0.0, -0.0, -0.1625, 0.1625, 0.1625,
    -0.0, -0.0, -0.0, -0.17500000000000002, 0.17500000000000002,
    0.17500000000000002, -0.0, -0.0, -0.0, -0.18750000000000003,
    0.18750000000000003, 0.18750000000000003, -0.0, -0.0, -0.0,
    -0.20000000000000004, 0.20000000000000004, 0.20000000000000004, -0.0, -0.0,
    -0.0, -0.21250000000000005, 0.21250000000000005, 0.21250000000000005, -0.0,
    -0.0, -0.0, -0.22500000000000006, 0.22500000000000006, 0.22500000000000006,
    -0.0, -0.0, -0.0, -0.23750000000000007, 0.23750000000000007,
    0.23750000000000007, -0.0, -0.0, -0.0, -0.25000000000000006,
    0.25000000000000006, 0.25000000000000006, 0.0, 0.0, 0.0, 0.0125, -0.0125,
    -0.0125, 0.0, 0.0, 0.0, 0.025, -0.025, -0.025, 0.0, 0.0, 0.0,
    0.037500000000000006, -0.037500000000000006, -0.037500000000000006, 0.0, 0.0,
    0.0, 0.05, -0.05, -0.05, 0.0, 0.0, 0.0, 0.0625, -0.0625, -0.0625, 0.0, 0.0,
    0.0, 0.075, -0.075, -0.075, 0.0, 0.0, 0.0, 0.0875, -0.0875, -0.0875, 0.0,
    0.0, 0.0, 0.099999999999999992, -0.099999999999999992, -0.099999999999999992,
    0.0, 0.0, 0.0, 0.11249999999999999, -0.11249999999999999,
    -0.11249999999999999, 0.0, 0.0, 0.0, 0.12499999999999999,
    -0.12499999999999999, -0.12499999999999999, 0.0, 0.0, 0.0,
    0.13749999999999998, -0.13749999999999998, -0.13749999999999998, 0.0, 0.0,
    0.0, 0.15, -0.15, -0.15, 0.0, 0.0, 0.0, 0.1625, -0.1625, -0.1625, 0.0, 0.0,
    0.0, 0.17500000000000002, -0.17500000000000002, -0.17500000000000002, 0.0,
    0.0, 0.0, 0.18750000000000003, -0.18750000000000003, -0.18750000000000003,
    0.0, 0.0, 0.0, 0.20000000000000004, -0.20000000000000004,
    -0.20000000000000004, 0.0, 0.0, 0.0, 0.21250000000000005,
    -0.21250000000000005, -0.21250000000000005, 0.0, 0.0, 0.0,
    0.22500000000000006, -0.22500000000000006, -0.22500000000000006, 0.0, 0.0,
    0.0, 0.23750000000000007, -0.23750000000000007, -0.23750000000000007, 0.0,
    0.0, 0.0, 0.25000000000000006, -0.25000000000000006, -0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0,
    -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0,
    -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0,
    -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5,
    0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, -0.0, -0.0, -0.0, -0.5, 0.5, 0.5,
    -0.0, -0.0, -0.0, -0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0,
    0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5,
    -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.0125, -0.0125, 0.0125,
    -0.0, -0.0, -0.0, 0.025, -0.025, 0.025, -0.0, -0.0, -0.0,
    0.037500000000000006, -0.037500000000000006, 0.037500000000000006, -0.0,
    -0.0, -0.0, 0.05, -0.05, 0.05, -0.0, -0.0, -0.0, 0.0625, -0.0625, 0.0625,
    -0.0, -0.0, -0.0, 0.075, -0.075, 0.075, -0.0, -0.0, -0.0, 0.0875, -0.0875,
    0.0875, -0.0, -0.0, -0.0, 0.099999999999999992, -0.099999999999999992,
    0.099999999999999992, -0.0, -0.0, -0.0, 0.11249999999999999,
    -0.11249999999999999, 0.11249999999999999, -0.0, -0.0, -0.0,
    0.12499999999999999, -0.12499999999999999, 0.12499999999999999, -0.0, -0.0,
    -0.0, 0.13749999999999998, -0.13749999999999998, 0.13749999999999998, -0.0,
    -0.0, -0.0, 0.15, -0.15, 0.15, -0.0, -0.0, -0.0, 0.1625, -0.1625, 0.1625,
    -0.0, -0.0, -0.0, 0.17500000000000002, -0.17500000000000002,
    0.17500000000000002, -0.0, -0.0, -0.0, 0.18750000000000003,
    -0.18750000000000003, 0.18750000000000003, -0.0, -0.0, -0.0,
    0.20000000000000004, -0.20000000000000004, 0.20000000000000004, -0.0, -0.0,
    -0.0, 0.21250000000000005, -0.21250000000000005, 0.21250000000000005, -0.0,
    -0.0, -0.0, 0.22500000000000006, -0.22500000000000006, 0.22500000000000006,
    -0.0, -0.0, -0.0, 0.23750000000000007, -0.23750000000000007,
    0.23750000000000007, -0.0, -0.0, -0.0, 0.25000000000000006,
    -0.25000000000000006, 0.25000000000000006, 0.0, 0.0, 0.0, -0.0125, 0.0125,
    -0.0125, 0.0, 0.0, 0.0, -0.025, 0.025, -0.025, 0.0, 0.0, 0.0,
    -0.037500000000000006, 0.037500000000000006, -0.037500000000000006, 0.0, 0.0,
    0.0, -0.05, 0.05, -0.05, 0.0, 0.0, 0.0, -0.0625, 0.0625, -0.0625, 0.0, 0.0,
    0.0, -0.075, 0.075, -0.075, 0.0, 0.0, 0.0, -0.0875, 0.0875, -0.0875, 0.0,
    0.0, 0.0, -0.099999999999999992, 0.099999999999999992, -0.099999999999999992,
    0.0, 0.0, 0.0, -0.11249999999999999, 0.11249999999999999,
    -0.11249999999999999, 0.0, 0.0, 0.0, -0.12499999999999999,
    0.12499999999999999, -0.12499999999999999, 0.0, 0.0, 0.0,
    -0.13749999999999998, 0.13749999999999998, -0.13749999999999998, 0.0, 0.0,
    0.0, -0.15, 0.15, -0.15, 0.0, 0.0, 0.0, -0.1625, 0.1625, -0.1625, 0.0, 0.0,
    0.0, -0.17500000000000002, 0.17500000000000002, -0.17500000000000002, 0.0,
    0.0, 0.0, -0.18750000000000003, 0.18750000000000003, -0.18750000000000003,
    0.0, 0.0, 0.0, -0.20000000000000004, 0.20000000000000004,
    -0.20000000000000004, 0.0, 0.0, 0.0, -0.21250000000000005,
    0.21250000000000005, -0.21250000000000005, 0.0, 0.0, 0.0,
    -0.22500000000000006, 0.22500000000000006, -0.22500000000000006, 0.0, 0.0,
    0.0, -0.23750000000000007, 0.23750000000000007, -0.23750000000000007, 0.0,
    0.0, 0.0, -0.25000000000000006, 0.25000000000000006, -0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0,
    0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5,
    0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5,
    -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0,
    -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0,
    -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0,
    0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, -0.0, -0.0, -0.0, 0.5,
    -0.5, 0.5, -0.0, -0.0, -0.0, 0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0,
    0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0,
    0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5,
    0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0,
    0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.0125, 0.0125,
    -0.0125, -0.0, -0.0, -0.0, 0.025, 0.025, -0.025, -0.0, -0.0, -0.0,
    0.037500000000000006, 0.037500000000000006, -0.037500000000000006, -0.0,
    -0.0, -0.0, 0.05, 0.05, -0.05, -0.0, -0.0, -0.0, 0.0625, 0.0625, -0.0625,
    -0.0, -0.0, -0.0, 0.075, 0.075, -0.075, -0.0, -0.0, -0.0, 0.0875, 0.0875,
    -0.0875, -0.0, -0.0, -0.0, 0.099999999999999992, 0.099999999999999992,
    -0.099999999999999992, -0.0, -0.0, -0.0, 0.11249999999999999,
    0.11249999999999999, -0.11249999999999999, -0.0, -0.0, -0.0,
    0.12499999999999999, 0.12499999999999999, -0.12499999999999999, -0.0, -0.0,
    -0.0, 0.13749999999999998, 0.13749999999999998, -0.13749999999999998, -0.0,
    -0.0, -0.0, 0.15, 0.15, -0.15, -0.0, -0.0, -0.0, 0.1625, 0.1625, -0.1625,
    -0.0, -0.0, -0.0, 0.17500000000000002, 0.17500000000000002,
    -0.17500000000000002, -0.0, -0.0, -0.0, 0.18750000000000003,
    0.18750000000000003, -0.18750000000000003, -0.0, -0.0, -0.0,
    0.20000000000000004, 0.20000000000000004, -0.20000000000000004, -0.0, -0.0,
    -0.0, 0.21250000000000005, 0.21250000000000005, -0.21250000000000005, -0.0,
    -0.0, -0.0, 0.22500000000000006, 0.22500000000000006, -0.22500000000000006,
    -0.0, -0.0, -0.0, 0.23750000000000007, 0.23750000000000007,
    -0.23750000000000007, -0.0, -0.0, -0.0, 0.25000000000000006,
    0.25000000000000006, -0.25000000000000006, 0.0, 0.0, 0.0, -0.0125, -0.0125,
    0.0125, 0.0, 0.0, 0.0, -0.025, -0.025, 0.025, 0.0, 0.0, 0.0,
    -0.037500000000000006, -0.037500000000000006, 0.037500000000000006, 0.0, 0.0,
    0.0, -0.05, -0.05, 0.05, 0.0, 0.0, 0.0, -0.0625, -0.0625, 0.0625, 0.0, 0.0,
    0.0, -0.075, -0.075, 0.075, 0.0, 0.0, 0.0, -0.0875, -0.0875, 0.0875, 0.0,
    0.0, 0.0, -0.099999999999999992, -0.099999999999999992, 0.099999999999999992,
    0.0, 0.0, 0.0, -0.11249999999999999, -0.11249999999999999,
    0.11249999999999999, 0.0, 0.0, 0.0, -0.12499999999999999,
    -0.12499999999999999, 0.12499999999999999, 0.0, 0.0, 0.0,
    -0.13749999999999998, -0.13749999999999998, 0.13749999999999998, 0.0, 0.0,
    0.0, -0.15, -0.15, 0.15, 0.0, 0.0, 0.0, -0.1625, -0.1625, 0.1625, 0.0, 0.0,
    0.0, -0.17500000000000002, -0.17500000000000002, 0.17500000000000002, 0.0,
    0.0, 0.0, -0.18750000000000003, -0.18750000000000003, 0.18750000000000003,
    0.0, 0.0, 0.0, -0.20000000000000004, -0.20000000000000004,
    0.20000000000000004, 0.0, 0.0, 0.0, -0.21250000000000005,
    -0.21250000000000005, 0.21250000000000005, 0.0, 0.0, 0.0,
    -0.22500000000000006, -0.22500000000000006, 0.22500000000000006, 0.0, 0.0,
    0.0, -0.23750000000000007, -0.23750000000000007, 0.23750000000000007, 0.0,
    0.0, 0.0, -0.25000000000000006, -0.25000000000000006, 0.25000000000000006,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0,
    -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0,
    -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0,
    0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5,
    -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, -0.0, -0.0, -0.0, 0.5, 0.5, -0.5,
    -0.0, -0.0, -0.0, 0.5, 0.5, -0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0,
    0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5,
    -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5,
    0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0,
    0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5,
    -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5,
    0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0,
    0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5,
    -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T k_0[984]{ -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.220466274921145, -0.220466274921145, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.220466274921145, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.220466274921145, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.220466274921145,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.220466274921145, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.220466274921145,
    0.220466274921145, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.220466274921145, -0.220466274921145,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407, -0.0, -0.0,
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

  static const real_T g[738]{ -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.220466274921145, -0.220466274921145, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.220466274921145, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.220466274921145, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.220466274921145,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.220466274921145, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.014433756729740647,
    0.014433756729740647, -0.028867513459481294, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.028506669541237775,
    0.028506669541237775, -0.05701333908247555, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.042227759532447479,
    0.042227759532447479, -0.084455519064894957, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.055605822273876934,
    0.055605822273876934, -0.11121164454775387, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.068649433446770658,
    0.068649433446770658, -0.13729886689354132, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.081366954340342038,
    0.081366954340342038, -0.16273390868068408, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.09376653721157413,
    0.09376653721157413, -0.18753307442314826, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.10585613051102542,
    0.10585613051102542, -0.21171226102205085, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.11764348397799043,
    0.11764348397799043, -0.23528696795598086, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.12913615360828132,
    0.12913615360828132, -0.25827230721656264, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.14034150649781493,
    0.14034150649781493, -0.28068301299562987, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1512667255651102,
    0.1512667255651102, -0.3025334511302204, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.1619188141557231,
    0.1619188141557231, -0.32383762831144619, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.17230460053157065,
    0.17230460053157065, -0.3446092010631413, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.18243074224802203,
    0.18243074224802203, -0.36486148449604405, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.19230373042156212,
    0.19230373042156212, -0.38460746084312425, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.20192989389076371,
    0.20192989389076371, -0.40385978778152742, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.21131540327323525,
    0.21131540327323525, -0.42263080654647051, 0.220466274921145,
    0.220466274921145, -0.44093254984229, 0.220466274921145, 0.220466274921145,
    -0.44093254984229, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, 0.22938837477785703, 0.22938837477785703,
    -0.45877674955571407, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.014433756729740647, -0.014433756729740647,
    0.028867513459481294, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.028506669541237775, -0.028506669541237775,
    0.05701333908247555, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.042227759532447479, -0.042227759532447479,
    0.084455519064894957, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.055605822273876934, -0.055605822273876934,
    0.11121164454775387, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.068649433446770658, -0.068649433446770658,
    0.13729886689354132, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.081366954340342038, -0.081366954340342038,
    0.16273390868068408, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.09376653721157413, -0.09376653721157413,
    0.18753307442314826, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.10585613051102542, -0.10585613051102542,
    0.21171226102205085, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.11764348397799043, -0.11764348397799043,
    0.23528696795598086, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.12913615360828132, -0.12913615360828132,
    0.25827230721656264, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.14034150649781493, -0.14034150649781493,
    0.28068301299562987, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1512667255651102, -0.1512667255651102,
    0.3025334511302204, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.1619188141557231, -0.1619188141557231,
    0.32383762831144619, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.17230460053157065, -0.17230460053157065,
    0.3446092010631413, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.18243074224802203, -0.18243074224802203,
    0.36486148449604405, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.19230373042156212, -0.19230373042156212,
    0.38460746084312425, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.20192989389076371, -0.20192989389076371,
    0.40385978778152742, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.21131540327323525, -0.21131540327323525,
    0.42263080654647051, -0.220466274921145, -0.220466274921145,
    0.44093254984229, -0.220466274921145, -0.220466274921145, 0.44093254984229,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407,
    -0.22938837477785703, -0.22938837477785703, 0.45877674955571407, -0.0, -0.0,
    -1.0, 0.0, 0.0, 1.0 };

  static const real_T b[324]{ 0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.975, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.025, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T c[288]{ 0.028867513459481294, -0.014433756729740647,
    -0.014433756729740647, 0.028867513459481294, -0.014433756729740647,
    -0.014433756729740647, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, -0.014433756729740647, 0.028867513459481294, -0.014433756729740647,
    -0.014433756729740647, 0.028867513459481294, -0.014433756729740647, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.014433756729740647,
    -0.014433756729740647, 0.028867513459481294, -0.014433756729740647,
    -0.014433756729740647, 0.028867513459481294, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05,
    0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.00062500000000000012, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.05, 0.00062500000000000012, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T l[180]{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

  static const real_T n[180]{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
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

  static const real_T d_0[108]{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.5, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5,
    0.5, -0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, -0.5, -0.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5, 0.5, -0.5, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.5, -0.5, 0.5 };

  static const real_T h[16]{ 0.034121465297356074, 0.0, 0.0, 0.0, 0.0,
    0.034121465297356074, 0.0, 0.0, 0.0, 0.0, 0.034121465297356074, 0.0, 0.0,
    0.0, 0.0, 100000.0 };

  static const int32_T b_Mrows[246]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
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
    197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211,
    212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
    227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241,
    242, 243, 301, 302, 303 };

  static const int16_T e[246]{ 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612, 612,
    612, 612, 612, 612, 612, 80, 80, 80, 0, 0, 0 };

  real_T Bv[378];
  real_T Abar[324];
  real_T rtb_A_b[324];
  real_T rtb_Q_a[324];
  real_T rtb_Transpose2_o[324];
  real_T rtb_Z_d[324];
  real_T rtb_y_n[324];
  real_T rtb_y_n_0[324];
  real_T b_B[288];
  real_T y[270];
  real_T b_Mlim[246];
  real_T b_y[162];
  real_T Product1_j[144];
  real_T Dv[126];
  real_T rseq[120];
  real_T rtb_Add_g[108];
  real_T rtb_C_n[108];
  real_T rtb_N_f[108];
  real_T rtb_Product2_mz[108];
  real_T rtb_Transpose2_l[108];
  real_T D_est[90];
  real_T rtb_R_p_tmp[90];
  real_T rtb_useq_a[63];
  real_T b_utarget[60];
  real_T rtb_B_p[54];
  real_T rtb_A_n[36];
  real_T rtb_N_l[36];
  real_T rtb_R_d[36];
  real_T tmp_6[36];
  real_T prms[24];
  real_T vseq[21];
  real_T b_xoff[18];
  real_T rtb_B_m[18];
  real_T tmp[18];
  real_T tmp_2[18];
  real_T h_0[16];
  real_T Sum_h[12];
  real_T rtb_C_jm[6];
  real_T rtb_ywtT[6];
  real_T tmp_0[6];
  real_T tmp_1[6];
  real_T DiscreteFilter1_tmp[3];
  real_T Sum2_c[3];
  real_T rtb_Sum_gu[3];
  real_T rtb_excitation[3];
  real_T tmp_3[3];
  real_T tmp_4[3];
  real_T tmp_5[3];
  real_T dwt;
  uint16_T waypt;
  int8_T b_I[324];
  int8_T P0_2_tmp[12];
  boolean_T tmp_a[246];
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
      __m128d tmp_7;
      __m128d tmp_8;
      __m128d tmp_9;
      real_T sigmoid_workspace_k_1;
      real_T sigmoid_workspace_x0;
      int32_T b_k;
      int32_T b_k_tmp;
      int32_T d;
      int32_T i;
      int32_T i_0;
      int32_T k;
      int32_T rtb_R_p_tmp_0;
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
      if (any_d(&rtY.currEv.r[0])) {
        rstP2 = false;
        b_k = 3;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (b_k - 3 < 3)) {
          if (rtU.yo[b_k]) {
            rstP2 = true;
            exitg1 = true;
          } else {
            b_k++;
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
              rtDW.P0_2[k + 12 * i] = rtY.P_p[((P0_2_tmp[i] - 1) * 24 +
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
      if (any_d(&rtY.currEv.r[3])) {
        rstP1 = false;
        b_k = 0;
        exitg1 = false;
        while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (b_k < 3)) {
          if (rtU.yo[b_k]) {
            rstP1 = true;
            exitg1 = true;
          } else {
            b_k++;
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
              rtDW.P0_1[k + 12 * i] = rtY.P_p[((P0_2_tmp[i] - 1) * 24 +
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
      b_k = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (b_k < 6)) {
        if (!rtU.enAdapt[b_k]) {
          c_y = false;
          exitg1 = true;
        } else {
          b_k++;
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
      b_k = 0;
      for (i_0 = 0; i_0 < 12; i_0++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S324>/Product1'

        (void)std::memcpy(&rtY.P_p[i], &Product1_j[b_k], 12U * sizeof(real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i_0] = Sum_h[i_0];
        i += 24;
        b_k += 12;
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
      i = 0;
      b_k = 0;
      for (i_0 = 0; i_0 < 12; i_0++) {
        // Outport: '<Root>/P' incorporates:
        //   Product: '<S328>/Product1'

        (void)std::memcpy(&rtY.P_p[i + 300], &Product1_j[b_k], 12U * sizeof
                          (real_T));

        // Outport: '<Root>/theta'
        rtY.theta[i_0 + 12] = Sum_h[i_0];
        i += 24;
        b_k += 12;
      }

      // Outport: '<Root>/prmErr' incorporates:
      //   Sum: '<S328>/Sum2'

      rtY.prmErr[3] = Sum2_c[0];
      rtY.prmErr[4] = Sum2_c[1];
      rtY.prmErr[5] = Sum2_c[2];

      // '<S1>:59:28' currTraj = traj(:,waypt);
      for (k = 0; k < 6; k++) {
        // Outport: '<Root>/currTraj'
        rtY.currTraj[k] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 6 +
          k];
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
      sigmoid_workspace_k_1 = 2.197 / (0.2 * rtU.k_2);

      // 'wtMod_:10' x0 = 0.5*k_2;
      sigmoid_workspace_x0 = 0.5 * rtU.k_2;

      //  midpoint
      // 'wtMod_:11' sigmoid = @(x) 1/(1 + exp(-k_1*(x-x0)));
      // 'wtMod_:13' for i = 1:2*no
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
      for (k = 0; k < 6; k++) {
        real_T rtb_ywtT_j;
        real_T rtb_ywt_d;

        // Delay: '<S9>/Delay'
        rtb_ywtT_j = rtDW.Delay_DSTATE[k];

        // MATLAB Function: '<S9>/MATLAB Function' incorporates:
        //   Inport: '<Root>/k_2'
        //   Outport: '<Root>/currEv'

        // 'wtMod_:14' if yDest(i) ~= 0
        if (rtY.currEv.r[k] != 0.0) {
          //  drive ywt to 1
          // 'wtMod_:16' if (ywtT(i) <= 1)
          if (rtb_ywtT_j <= 1.0) {
            // 'wtMod_:17' ywtT(i) = ywtT(i) + dwt;
            rtb_ywtT_j += dwt;
          }

          // 'wtMod_:19' else
          //  drive ywt to 0
          // 'wtMod_:21' if (ywtT(i) > 0)
        } else if (rtb_ywtT_j > 0.0) {
          // 'wtMod_:22' ywtT(i) = ywtT(i) - dwt;
          rtb_ywtT_j -= dwt;
        } else {
          // no actions
        }

        // 'wtMod_:25' if ywtT(i) <= 0
        if (rtb_ywtT_j <= 0.0) {
          // 'wtMod_:26' ywt(i) = 0;
          rtb_ywt_d = 0.0;
        } else {
          // 'wtMod_:27' else
          // 'wtMod_:28' ywt(i) = sigmoid(ywtT(i)*k_2);
          // 'wtMod_:11' @(x) 1/(1 + exp(-k_1*(x-x0)))
          rtb_ywt_d = 1.0 / (std::exp((rtb_ywtT_j * rtU.k_2 -
            sigmoid_workspace_x0) * -sigmoid_workspace_k_1) + 1.0);
        }

        // Update for Delay: '<S9>/Delay'
        rtDW.Delay_DSTATE[k] = rtb_ywtT_j;

        // Outport: '<Root>/ywt' incorporates:
        //   Gain: '<S9>/Gain'

        rtY.ywt[k] = rtP.beta * rtb_ywt_d;
      }

      // End of Outputs for SubSystem: '<S1>/wtMod'

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // Gain: '<S2>/Gain1' incorporates:
      //   Inport: '<Root>/uwt'
      //   Sum: '<S2>/Sum'

      // '<S1>:59:31' [u, yhat] = ampc(currTraj, y, ymax, ywt, y0, x0, u0, umax, uwt, excitation, theta, thetaSgn); 
      // Simulink Function 'ampc': '<S1>:461'
      rtb_Sum_gu[0] = rtP.beta * rtU.uwt[0];
      rtb_Sum_gu[1] = rtP.beta * rtU.uwt[1];
      rtb_Sum_gu[2] = rtP.beta * rtU.uwt[2];

      // MATLAB Function: '<S2>/MATLAB Function2' incorporates:
      //   DiscreteFilter: '<S2>/Discrete Filter1'
      //   Outport: '<Root>/theta'
      //   RandomNumber: '<S2>/excitation'
      //
      // MATLAB Function 'SupervisoryController/ampc/MATLAB Function2': '<S11>:1' 
      // '<S11>:1:2' [A, B] = theta2ss_(theta, sign, no, ni, np, ns, dt, mdlNum); 
      // 'theta2ss_:3' prms = zeros(np, 2*no);
      //  parameter matrix
      // 'theta2ss_:4' A = zeros(2*ns, 2*ns);
      (void)std::memset(&rtb_A_n[0], 0, 36U * sizeof(real_T));

      // 'theta2ss_:5' B = zeros(2*ns, ni);
      // 'theta2ss_:7' for i=1:2*no
      for (k = 0; k < 6; k++) {
        //  apply correct sign to each np segment of theta and set it as a column vector of prms 
        // 'theta2ss_:8' prms(:,i) = sign( ((i-1)*np + 1):i*np ) .* theta( ((i-1)*np + 1):i*np ); 
        i_0 = k << 2UL;
        i = i_0;
        b_k_tmp = (k + 1) << 2UL;
        if (i_0 + 1 > b_k_tmp) {
          i = 0;
          d = 0;
        } else {
          d = b_k_tmp;
        }

        b_k = i_0;
        if (i_0 + 1 > b_k_tmp) {
          b_k = 0;
          b_k_tmp = 0;
        }

        if (d - i == b_k_tmp - b_k) {
          prms[i_0] = rtDW.thetaSgn[i] * rtY.theta[b_k];
          prms[(k << 2UL) + 1] = rtDW.thetaSgn[i + 1] * rtY.theta[b_k + 1];
          prms[i_0 + 2] = rtDW.thetaSgn[i + 2] * rtY.theta[b_k + 2];
          prms[i_0 + 3] = rtDW.thetaSgn[i + 3] * rtY.theta[b_k + 3];
        } else {
          binary_expand_op_n(prms, k, rtDW.thetaSgn, i, d - 1, rtY.theta, b_k,
                             b_k_tmp - 1);
        }

        // 'theta2ss_:9' A(i, i) = [1 - prms(1,i)];
        rtb_A_n[k + 6 * k] = 1.0 - prms[i_0];

        // 'theta2ss_:10' B(i, :) = prms(2:end, i)';
        rtb_B_m[k] = prms[(k << 2UL) + 1];
        rtb_B_m[k + 6] = prms[(k << 2UL) + 2];
        rtb_B_m[k + 12] = prms[(k << 2UL) + 3];
      }

      // End of MATLAB Function: '<S2>/MATLAB Function2'

      // Delay: '<S43>/MemoryX' incorporates:
      //   Constant: '<S43>/X0'
      //   DataTypeConversion: '<S43>/DataTypeConversionReset'

      rtDW.icLoad_k = ((static_cast<uint32_T>(rtPrevZCX.MemoryX_Reset_ZCE_o) ==
                        POS_ZCSIG) || rtDW.icLoad_k);
      rtPrevZCX.MemoryX_Reset_ZCE_o = 0U;
      if (rtDW.icLoad_k) {
        (void)std::memcpy(&rtDW.MemoryX_DSTATE_d[0], &rtP.X0_Value[0], 18U *
                          sizeof(real_T));
      }

      // MATLAB Function: '<S41>/FixedHorizonOptimizer' incorporates:
      //   BusCreator: '<S2>/Bus Creator1'
      //   Constant: '<S2>/Constant12'
      //   Constant: '<S2>/Constant2'
      //   DiscreteFilter: '<S2>/Discrete Filter1'
      //   Inport: '<Root>/u0'
      //   Inport: '<Root>/x0'
      //   Inport: '<Root>/y0'
      //   Outport: '<Root>/currTraj'
      //   RandomNumber: '<S2>/excitation'
      //
      // MATLAB Function 'Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer': '<S42>:1' 
      // '<S42>:1:18' coder.extrinsic('mpcblock_optimizer_double_mex');
      // '<S42>:1:19' coder.extrinsic('mpcblock_optimizer_single_mex');
      // '<S42>:1:20' coder.extrinsic('mpcblock_refmd_double_mex');
      // '<S42>:1:21' coder.extrinsic('mpcblock_refmd_single_mex');
      //  Inputs (in BlockDataType except iA)
      //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
      // '<S42>:1:25' xk = convertDataType(xk0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    old_u:      last mv (calculated by MPC)
      // '<S42>:1:27' old_u = convertDataType(old_u0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    ym:         current measured output (used only with built-in KF)
      // '<S42>:1:29' ym = convertDataType(ym0,isDouble);
      //    ref:        output reference
      // '<S42>:1:31' ref = convertDataType(ref0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    md:         measured disturbance
      // '<S42>:1:33' md = convertDataType(md0,isDouble);
      //    umin:       run-time MV bound
      // '<S42>:1:35' umin = convertDataType(umin0,isDouble);
      //    umax:       run-time MV bound
      // '<S42>:1:37' umax = convertDataType(umax0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    ymin:       run-time OV bound
      // '<S42>:1:39' ymin = convertDataType(ymin0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    ymax:       run-time OV bound
      // '<S42>:1:41' ymax = convertDataType(ymax0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    E:          run-time mixed constraints
      // '<S42>:1:43' E = convertDataType(E0,isDouble);
      //    F:          run-time mixed constraints
      // '<S42>:1:45' F = convertDataType(F0,isDouble);
      //    G:          run-time mixed constraints
      // '<S42>:1:47' G = convertDataType(G0,isDouble);
      //    S:          run-time mixed constraints
      // '<S42>:1:49' S = convertDataType(S0,isDouble);
      //    switch_in:  if it matches "enable_value", MPC is active in control
      // '<S42>:1:51' switch_in = int32(switch_in0);
      //    ext_mv:     external last mv (actual)
      // '<S42>:1:53' ext_mv = convertDataType(ext_mv0,isDouble);
      //    MVtarget:   MV reference
      // '<S42>:1:55' MVtarget = convertDataType(MVtarget0,isDouble);
      //    ywt:        run-time OV weights
      // '<S42>:1:57' ywt = convertDataType(ywt0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    uwt:        run-time MV weights
      // '<S42>:1:59' uwt = convertDataType(uwt0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    duwt:       run-time DMV weights
      // '<S42>:1:61' duwt = convertDataType(duwt0,isDouble);
      //    rhoeps:     run-time Slack weights
      // '<S42>:1:63' ewt = convertDataType(ewt0,isDouble);
      //    a:          run-time A (must be in DT)
      // '<S42>:1:65' a = convertDataType(a0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    b:          run-time B (must be in DT)
      // '<S42>:1:67' b = convertDataType(b0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    c:          run-time C (must be in DT)
      // '<S42>:1:69' c = convertDataType(c0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    d:          run-time D (must be in DT)
      // '<S42>:1:71' d = convertDataType(d0,isDouble);
      //    U:          run-time nominal value
      // '<S42>:1:73' U = convertDataType(U0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    Y:          run-time nominal value
      // '<S42>:1:75' Y = convertDataType(Y0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    X:          run-time nominal value
      // '<S42>:1:77' X = convertDataType(X0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    DX:         run-time nominal value
      // '<S42>:1:79' DX = convertDataType(DX0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
      //    Pk:         covariance P[k|k-1] (used only with built-in KF)
      // '<S42>:1:81' Pk = convertDataType(Pk0,isDouble);
      // '<S42>:1:317' if isDouble
      //  convert an input signal to double precision when necessary
      // '<S42>:1:319' if isa(u,'double')
      // '<S42>:1:320' y = u;
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
      //    Pk1:        covariance P[k+1|k]
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
      //  Constants
      // '<S42>:1:115' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
      // '<S42>:1:116' isAdaptive = ~isLTV;
      //  isLTV=true forces isAdaptive=false
      // '<S42>:1:117' ZERO = zeros('like',ref);
      // '<S42>:1:118' ONE = ones('like',ref);
      // '<S42>:1:119' hasMD = nv>int32(1);
      //  Pre-allocate all the MEX block outputs for the simulation mode
      // '<S42>:1:123' if isSimulation
      //  Model update
      // '<S42>:1:137' nym = int32(length(myindex));
      // '<S42>:1:138' ai=zeros(nxp,nxp,'like',ref);
      // '<S42>:1:139' bi=zeros(nxp,nup,'like',ref);
      // '<S42>:1:140' ci=zeros(ny,nxp,'like',ref);
      // '<S42>:1:141' di=zeros(ny,nup,'like',ref);
      // '<S42>:1:143' ai(:,:)=a(:,:,1);
      // '<S42>:1:144' bi(:,:)=b(:,:,1);
      // '<S42>:1:145' ci(:,:)=c(:,:,1);
      // '<S42>:1:146' di(:,:)=d(:,:,1);
      //  Allocate matrices. Must allocate 3D matrix also in Adaptive case,
      //  otherwise EML code does not compile.
      // '<S42>:1:150' Bu=zeros(nx,nu,p+1,'like',ref);
      (void)std::memset(&rtDW.Bu[0], 0, 1134U * sizeof(real_T));

      // '<S42>:1:151' Bv=zeros(nx,nv,p+1,'like',ref);
      (void)std::memset(&Bv[0], 0, 378U * sizeof(real_T));

      // '<S42>:1:152' Dv=zeros(ny,nv,p+1,'like',ref);
      (void)std::memset(&Dv[0], 0, 126U * sizeof(real_T));

      // '<S42>:1:153' Dvm=zeros(nym,nv,p+1,'like',ref);
      // '<S42>:1:154' Cm=zeros(nym,nx,p+1,'like',ref);
      // '<S42>:1:155' [A(:,:,1),C(:,:,1),Bu(:,:,1),Bv(:,:,1),Cm(:,:,1),Dv(:,:,1),Dvm(:,:,1),Qk,Rk,Nk] = mpc_plantupdate(... 
      // '<S42>:1:156'     ai,bi,ci,di,A(:,:,1),B(:,:,1),C(:,:,1),D(:,:,1),mvindex,mdindex,unindex,nxp,nup,ny,nu,nv,nxid, ... 
      // '<S42>:1:157'     myindex,Uscale,Yscale,Cid,Did);
      (void)std::memcpy(&rtb_A_b[0], &b[0], 324U * sizeof(real_T));
      (void)std::memcpy(&b_B[0], &c[0], 288U * sizeof(real_T));
      (void)std::memcpy(&rtb_C_n[0], &d_0[0], 108U * sizeof(real_T));
      k = 0;
      i = 0;
      for (b_k = 0; b_k < 6; b_k++) {
        for (i_0 = 0; i_0 < 6; i_0++) {
          b_k_tmp = i_0 + k;
          rtb_C_n[b_k_tmp] = rtP.Constant12_Value[b_k_tmp];
          rtb_A_b[i_0 + i] = rtb_A_n[b_k_tmp];
        }

        k += 6;
        i += 18;
        i_0 = 0;
        b_k_tmp = 0;
        for (d = 0; d < 3; d++) {
          b_B[i_0 + b_k] = rtb_B_m[b_k_tmp + b_k];
          i_0 += 18;
          b_k_tmp += 6;
        }
      }

      (void)std::memcpy(&rtDW.Bu[0], &b_B[0], 54U * sizeof(real_T));
      (void)std::memcpy(&Bv[0], &b_B[54], 18U * sizeof(real_T));
      for (k = 0; k < 6; k++) {
        Dv[k] = 0.0;
      }

      // '<S42>:1:158' if isLTV
      //  Offset update together with Mlim, utarget, Bv and Dv values
      // '<S42>:1:174' [Mlim, utarget, uoff, voff, yoff, myoff, xoff, Bv, Dv] = ... 
      // '<S42>:1:175'     mpc_updateFromNominal(isAdaptive,isQP,Mlim,Mrows,...
      // '<S42>:1:176'        U,Uscale,uoff,mvindex,voff,mdindex,utarget,nu,nv-1,... 
      // '<S42>:1:177'        Y,Yscale,yoff,myoff,myindex,ny,...
      // '<S42>:1:178'        X,xoff,nxp,DX,A,Bu,Bv,C,Dv,nCC);
      for (i = 0; i < 246; i++) {
        b_Mlim[i] = static_cast<real_T>(e[i]);
      }

      (void)std::memset(&b_utarget[0], 0, 60U * sizeof(real_T));
      (void)std::memset(&b_xoff[0], 0, 18U * sizeof(real_T));
      Sum2_c[0] = rtU.u0[0];
      Sum2_c[1] = rtU.u0[1];
      Sum2_c[2] = rtU.u0[2];
      for (i = 0; i < 6; i++) {
        rtb_ywtT[i] = rtU.y0[i];
      }

      DiscreteFilter1_tmp[0] = rtU.u0[0];
      DiscreteFilter1_tmp[1] = rtU.u0[1];
      DiscreteFilter1_tmp[2] = rtU.u0[2];
      for (b_k = 0; b_k < 246; b_k++) {
        dwt = b_Mlim[b_k];
        k = b_Mrows[b_k];
        if (k <= 120) {
          dwt += 0.0 - rtb_ywtT[(k - (k - 1) / static_cast<int32_T>(ny) *
            static_cast<int32_T>(ny)) - 1];
        } else if (k <= 240) {
          dwt -= 0.0 - rtb_ywtT[(k - div_nde_s32_floor(k - 121,
            static_cast<int32_T>(ny)) * static_cast<int32_T>(ny)) - 121];
        } else if (k <= 300) {
          dwt += 0.0 - Sum2_c[(k - div_nde_s32_floor(k - 241, static_cast<
            int32_T>(nu)) * static_cast<int32_T>(nu)) - 241];
        } else {
          dwt -= 0.0 - Sum2_c[(k - div_nde_s32_floor(k - 301,
            static_cast<int32_T>(nu)) * static_cast<int32_T>(nu)) - 301];
        }

        b_Mlim[b_k] = dwt;
      }

      for (b_k = 0; b_k < 3; b_k++) {
        dwt = Sum2_c[b_k];
        k = 0;
        for (i = 0; i < 20; i++) {
          i_0 = k + b_k;
          b_utarget[i_0] -= dwt;
          k += 3;
        }
      }

      for (k = 0; k < 6; k++) {
        b_xoff[k] = rtU.x0[k];
        Bv[k] = rtP.Constant2_Value[k];
      }

      //  Remove last u offset
      // '<S42>:1:181' old_u = old_u - uoff;
      //  Get reference and MD signals -- accounting for previewing
      // '<S42>:1:184' if isSimulation
      // '<S42>:1:190' else
      //  When doing code generation, use M code directly
      // '<S42>:1:192' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
      for (k = 0; k < 21; k++) {
        vseq[k] = 1.0;
      }

      i = 0;
      for (k = 0; k < 20; k++) {
        for (b_k = 0; b_k <= 4; b_k += 2) {
          tmp_8 = _mm_loadu_pd(&rtY.currTraj[b_k]);
          tmp_9 = _mm_loadu_pd(&rtb_ywtT[b_k]);
          (void)_mm_storeu_pd(&rseq[b_k + i], _mm_sub_pd(tmp_8, tmp_9));
        }

        i += 6;
      }

      //  External MV override.
      //  NOTE: old_u and ext_mv input signals are dimensionless and offset-free. 
      // '<S42>:1:197' if no_mv
      //  no external mv: old_u is the optimal u[k-1] from last_mv
      // '<S42>:1:199' delmv = zeros(nu,1,'like',ref);
      //  Obtain x[k|k]
      // '<S42>:1:208' xk = xk - xoff;
      //  Remove offset
      // '<S42>:1:209' if CustomEstimation
      //  Input is x(k|k)
      // '<S42>:1:211' xest = xk;
      //  Real-time MV target override
      //  Note: utargetValue is a vector length p*nu.
      // '<S42>:1:231' if no_uref
      //  no external utarget
      // '<S42>:1:233' utargetValue = utarget;
      //  Real-time custom constraint override (scaled E/F/S)
      // '<S42>:1:241' if ~no_cc
      // '<S42>:1:250' return_sequence = return_mvseq || return_xseq || return_ovseq; 
      // '<S42>:1:251' if isSimulation
      // '<S42>:1:279' else
      //  When doing code generation, use M code directly
      // '<S42>:1:281' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
      // '<S42>:1:282'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
      // '<S42>:1:283'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
      // '<S42>:1:284'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
      // '<S42>:1:285'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
      // '<S42>:1:286'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
      // '<S42>:1:287'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
      // '<S42>:1:288'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
      // '<S42>:1:289'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
      // '<S42>:1:290'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
      // '<S42>:1:291'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
      // '<S42>:1:292'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, isHyb, Mdis, Ndis, Vdis, numdis, maxdis); 
      for (k = 0; k < 6; k++) {
        // Sum: '<S12>/Sum2' incorporates:
        //   Inport: '<Root>/x0'

        tmp[k] = rtU.x0[k];
      }

      // Sum: '<S12>/Sum2' incorporates:
      //   Constant: '<S12>/Constant1'

      (void)std::memcpy(&tmp[6], &rtP.Constant1_Value[0], 12U * sizeof(real_T));

      // End of Outputs for SubSystem: '<S1>/ampc'
      for (k = 0; k <= 4; k += 2) {
        // Inport: '<Root>/ymax'
        tmp_8 = _mm_loadu_pd(&rtU.ymax[k]);

        // Outputs for Function Call SubSystem: '<S1>/ampc'
        // Gain: '<S2>/Gain2' incorporates:
        //   Inport: '<Root>/ymax'

        (void)_mm_storeu_pd(&tmp_0[k], _mm_mul_pd(_mm_set1_pd(rtP.Gain2_Gain),
          tmp_8));

        // Gain: '<S2>/Gain3' incorporates:
        //   Inport: '<Root>/ymax'

        (void)_mm_storeu_pd(&tmp_1[k], _mm_mul_pd(_mm_set1_pd(rtP.Gain3_Gain),
          tmp_8));

        // End of Outputs for SubSystem: '<S1>/ampc'
      }

      for (k = 0; k <= 16; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/ampc'
        // Delay: '<S43>/MemoryX'
        tmp_8 = _mm_loadu_pd(&rtDW.MemoryX_DSTATE_d[k]);

        // Sum: '<S12>/Sum2' incorporates:
        //   Delay: '<S43>/MemoryX'

        tmp_9 = _mm_loadu_pd(&tmp[k]);

        // MATLAB Function: '<S41>/FixedHorizonOptimizer' incorporates:
        //   Delay: '<S43>/MemoryX'

        tmp_7 = _mm_loadu_pd(&b_xoff[k]);
        (void)_mm_storeu_pd(&tmp_2[k], _mm_sub_pd(_mm_add_pd(tmp_8, tmp_9),
          tmp_7));

        // End of Outputs for SubSystem: '<S1>/ampc'
      }

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // MATLAB Function: '<S41>/FixedHorizonOptimizer' incorporates:
      //   Inport: '<Root>/u0'
      //   UnitDelay: '<S13>/last_mv'

      tmp_3[0] = rtDW.last_mv_DSTATE_m[0] - rtU.u0[0];
      tmp_3[1] = rtDW.last_mv_DSTATE_m[1] - rtU.u0[1];
      tmp_3[2] = rtDW.last_mv_DSTATE_m[2] - rtU.u0[2];
      (void)std::memset(&rtDW.dv[0], 0, 5166U * sizeof(real_T));
      tmp_4[0] = 0.034121465297356074;
      tmp_4[1] = 0.034121465297356074;
      tmp_4[2] = 0.034121465297356074;
      for (k = 0; k < 6; k++) {
        rtb_C_jm[k] = 1.0;
      }

      tmp_5[0] = 1.0;
      tmp_5[1] = 1.0;
      tmp_5[2] = 1.0;

      // Memory: '<S13>/Memory'
      (void)std::memcpy(&tmp_a[0], &rtDW.Memory_PreviousInput_j[0], 246U *
                        sizeof(boolean_T));

      // MATLAB Function: '<S41>/FixedHorizonOptimizer'
      (void)std::memcpy(&rtDW.f[0], &f[0], 4428U * sizeof(real_T));
      (void)std::memcpy(&rtDW.g[0], &g[0], 738U * sizeof(real_T));
      (void)std::memcpy(&h_0[0], &h[0], sizeof(real_T) << 4UL);
      (void)std::memcpy(&rtDW.k[0], &k_0[0], 984U * sizeof(real_T));

      // Update for Memory: '<S13>/Memory' incorporates:
      //   Inport: '<Root>/umax'
      //   MATLAB Function: '<S41>/FixedHorizonOptimizer'
      //   Math: '<S13>/Math Function1'
      //   Outport: '<Root>/ywt'
      //   Sum: '<S2>/Sum'

      mpcblock_optimizer(rseq, vseq, rtU.umax, tmp_0, tmp_1, tmp_2, tmp_3, tmp_a,
                         b_Mlim, rtDW.f, rtDW.g, rtDW.dv, b_utarget,
                         DiscreteFilter1_tmp, rtb_ywtT, h_0, rtDW.k, rtY.ywt,
                         rtb_Sum_gu, tmp_4, l, n, rtb_A_b, rtDW.Bu, Bv, rtb_C_n,
                         Dv, b_Mrows, rtb_C_jm, tmp_5, Sum2_c, rtb_useq_a, &dwt,
                         rtDW.Memory_PreviousInput_j);

      // Delay: '<S43>/MemoryP' incorporates:
      //   Constant: '<S43>/P0'
      //   DataTypeConversion: '<S43>/DataTypeConversionReset'

      // '<S42>:1:295' if return_xseq || return_ovseq
      // '<S42>:1:297' else
      // '<S42>:1:298' yseq = zeros(p+1,ny,'like',rseq);
      // '<S42>:1:299' xseq = zeros(p+1,nxQP,'like',rseq);
      // '<S42>:1:302' if CustomEstimation
      // '<S42>:1:303' xk1 = zeros(nx,1,'like',ref);
      // '<S42>:1:304' Pk1 = Pk;
      // '<S42>:1:311' xk1 = xk1 + xoff;
      //  Updated state must include offset
      //  return xest in original value
      // '<S42>:1:314' xest = xest + xoff;
      rtDW.icLoad_j = ((static_cast<uint32_T>(rtPrevZCX.MemoryP_Reset_ZCE_b) ==
                        POS_ZCSIG) || rtDW.icLoad_j);
      rtPrevZCX.MemoryP_Reset_ZCE_b = 0U;
      if (rtDW.icLoad_j) {
        (void)std::memcpy(&rtDW.MemoryP_DSTATE_h4[0], &rtP.P0_Value[0], 324U *
                          sizeof(real_T));
      }

      // MATLAB Function: '<S12>/MATLAB Function' incorporates:
      //   BusCreator: '<S2>/Bus Creator1'
      //   Constant: '<S2>/Constant12'
      //   Constant: '<S2>/Constant13'

      // MATLAB Function 'SupervisoryController/ampc/State Estimator OD (KF)/MATLAB Function': '<S44>:1' 
      // '<S44>:1:2' [A, B, C, D, Q, R, N] = stateEstimator(Ap, Bp, Cp, Dp, Aod, Bod, Cod, Dod, Dmn, no, ni, ns); 
      // 'stateEstimator:3' nsp = 2*ns_;
      //  n_plant_states
      // 'stateEstimator:4' nsod = size(Aod,1);
      //  n_od_states
      // 'stateEstimator:5' ns = nsp + nsod;
      //  n_states = n_plant_states + n_od_states
      // 'stateEstimator:7' A = zeros(ns);
      //  n_states x n_states
      // 'stateEstimator:8' B = zeros(ns,ni);
      (void)std::memset(&rtb_B_p[0], 0, 54U * sizeof(real_T));

      //  n_states  x n_inputs
      // 'stateEstimator:9' C = zeros(2*no,ns);
      //  n_outputs x n_states
      // 'stateEstimator:10' D = zeros(2*no,ni);
      //  n_outputs x n_inputs
      // 'stateEstimator:11' Q = zeros(ns,ns);
      //  n_states  x n_states
      // 'stateEstimator:12' G = eye(ns);
      //  n_states  x n_states
      // 'stateEstimator:13' R = zeros(2*no,2*no);
      //  n_outputs x n_outputs
      // 'stateEstimator:14' N = zeros(ns,2*no);
      //  n_states  x n_outputs
      // 'stateEstimator:15' H = zeros(2*no,ns);
      //  n_outputs x n_states
      //  combine plant and output disturbance model
      //  (force the outputs to fit in preallocated memory)
      // 'stateEstimator:19' A(1:ns, 1:ns) = blkdiag(Ap, Aod);
      (void)std::memset(&rtb_A_b[0], 0, 324U * sizeof(real_T));
      k = 0;
      i = 0;
      for (b_k = 0; b_k < 6; b_k++) {
        for (i_0 = 0; i_0 < 6; i_0++) {
          rtb_A_b[i_0 + k] = rtb_A_n[i_0 + i];
        }

        k += 18;
        i += 6;
      }

      k = 0;
      i = 0;
      for (b_k = 0; b_k < 12; b_k++) {
        (void)std::memcpy(&rtb_A_b[k + 114], &rtP.Aod[i], 12U * sizeof(real_T));
        k += 18;
        i += 12;
      }

      // 'stateEstimator:20' B(1:nsp, 1:ni) = Bp;
      for (k = 0; k < 6; k++) {
        rtb_B_p[k] = rtb_B_m[k];
        rtb_B_p[k + 18] = rtb_B_m[k + 6];
        rtb_B_p[k + 36] = rtb_B_m[k + 12];
      }

      // 'stateEstimator:21' C(1:2*no, 1:ns) = [Cp Cod];
      (void)std::memcpy(&rtb_C_n[0], &rtP.Constant12_Value[0], 36U * sizeof
                        (real_T));
      (void)std::memcpy(&rtb_C_n[36], &rtP.Cod[0], 72U * sizeof(real_T));

      // 'stateEstimator:22' D(1:2*no, 1:ni) = Dp;
      // 'stateEstimator:24' B_est = zeros(ns, ni + 2*no + 2*no);
      (void)std::memset(&y[0], 0, 270U * sizeof(real_T));

      // 'stateEstimator:25' B_est(1:ns, 1:ni+2*no) = blkdiag(Bp, Bod);
      (void)std::memset(&b_y[0], 0, 162U * sizeof(real_T));
      for (k = 0; k < 6; k++) {
        b_y[k] = rtb_B_m[k];
        b_y[k + 18] = rtb_B_m[k + 6];
        b_y[k + 36] = rtb_B_m[k + 12];
      }

      k = 0;
      i = 0;
      for (b_k = 0; b_k < 6; b_k++) {
        (void)std::memcpy(&b_y[k + 60], &rtP.Bod[i], 12U * sizeof(real_T));
        k += 18;
        i += 12;
      }

      (void)std::memcpy(&y[0], &b_y[0], 162U * sizeof(real_T));

      // 'stateEstimator:26' D_est = [Dp Dod Dn];
      (void)std::memcpy(&D_est[0], &rtP.Constant13_Value[0], 18U * sizeof(real_T));
      for (k = 0; k < 36; k++) {
        D_est[k + 18] = rtP.Dod[k];
        D_est[k + 54] = rtP.Dmn[k];
      }

      // 'stateEstimator:27' Q(1:ns, 1:ns) = blkdiag(Bp(1:ns_,1:ni)*Bp(1:ns_,1:ni)', Bp(ns_+1:2*ns_,1:ni)*Bp(ns_+1:2*ns_,1:ni)', Bod*Bod'); 
      (void)std::memset(&rtb_Q_a[0], 0, 324U * sizeof(real_T));
      for (k = 0; k < 3; k++) {
        i = 0;
        for (b_k = 0; b_k < 3; b_k++) {
          d = i + k;
          rtb_Q_a[d] = 0.0;
          rtb_Q_a[d] += rtb_B_m[k] * rtb_B_m[b_k];
          rtb_Q_a[d] += rtb_B_m[k + 6] * rtb_B_m[b_k + 6];
          rtb_Q_a[d] += rtb_B_m[k + 12] * rtb_B_m[b_k + 12];
          rtb_Q_a[d + 57] = 0.0;
          rtb_Q_a[d + 57] += rtb_B_m[k + 3] * rtb_B_m[b_k + 3];
          rtb_Q_a[d + 57] += rtb_B_m[k + 9] * rtb_B_m[b_k + 9];
          rtb_Q_a[d + 57] += rtb_B_m[k + 15] * rtb_B_m[b_k + 15];
          i += 18;
        }
      }

      for (k = 0; k < 12; k++) {
        i = 0;
        for (b_k = 0; b_k < 12; b_k++) {
          d = (i + k) + 114;
          rtb_Q_a[d] = 0.0;
          i_0 = 0;
          for (b_k_tmp = 0; b_k_tmp < 6; b_k_tmp++) {
            rtb_Q_a[d] += rtP.Bod[i_0 + k] * rtP.Bod[i_0 + b_k];
            i_0 += 12;
          }

          i += 18;
        }
      }

      //  Q = B_est * B_est';
      // 'stateEstimator:29' R = D_est * D_est';
      k = 0;
      for (i = 0; i < 6; i++) {
        b_k = 0;
        for (i_0 = 0; i_0 < 15; i_0++) {
          rtb_R_p_tmp[i_0 + k] = D_est[b_k + i];
          b_k += 6;
        }

        k += 15;
      }

      k = 0;
      i = 0;
      for (b_k = 0; b_k < 6; b_k++) {
        for (i_0 = 0; i_0 < 6; i_0++) {
          rtb_R_p_tmp_0 = i_0 + k;
          rtb_R_d[rtb_R_p_tmp_0] = 0.0;
          b_k_tmp = 0;
          for (d = 0; d < 15; d++) {
            rtb_R_d[rtb_R_p_tmp_0] += D_est[b_k_tmp + i_0] * rtb_R_p_tmp[d + i];
            b_k_tmp += 6;
          }
        }

        k += 6;
        i += 15;
      }

      // Outputs for Atomic SubSystem: '<S43>/ScalarExpansionQ'
      // 'stateEstimator:30' N = B_est * D_est';
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
      // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S88>:1'
      //    Copyright 2014-2015 The MathWorks, Inc.
      // '<S88>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
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
      // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S87>:1'
      //    Copyright 2014-2015 The MathWorks, Inc.
      // '<S87>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
      for (k = 0; k < 18; k++) {
        for (i = 0; i < 6; i++) {
          i_0 = 18 * i + k;
          rtb_N_f[i_0] = 0.0;
          for (b_k = 0; b_k < 15; b_k++) {
            rtb_N_f[i_0] += y[18 * b_k + k] * rtb_R_p_tmp[15 * i + b_k];
          }
        }

        for (i = 0; i < 18; i++) {
          // MATLAB Function: '<S65>/ScalarExpansion'
          d = 18 * k + i;
          rtb_y_n[d] = (rtb_Q_a[18 * i + k] + rtb_Q_a[d]) / 2.0;
        }
      }

      // End of MATLAB Function: '<S12>/MATLAB Function'
      // End of Outputs for SubSystem: '<S43>/ScalarExpansionQ'

      // Outputs for Atomic SubSystem: '<S43>/ReducedQRN'
      for (k = 0; k < 18; k++) {
        // Product: '<S63>/Product' incorporates:
        //   Constant: '<S43>/G'
        //   Math: '<S63>/Transpose1'

        i = 0;
        for (b_k = 0; b_k < 18; b_k++) {
          d = i + k;
          rtb_y_n_0[d] = 0.0;
          i_0 = 0;
          for (b_k_tmp = 0; b_k_tmp < 18; b_k_tmp++) {
            rtb_y_n_0[d] += rtb_y_n[i_0 + k] * rtP.G_Value[i_0 + b_k];
            i_0 += 18;
          }

          i += 18;
        }
      }

      // Product: '<S63>/Product' incorporates:
      //   Constant: '<S43>/G'

      k = 0;
      for (i = 0; i < 18; i++) {
        for (b_k = 0; b_k < 18; b_k++) {
          d = b_k + k;
          rtb_Q_a[d] = 0.0;
          i_0 = 0;
          for (b_k_tmp = 0; b_k_tmp < 18; b_k_tmp++) {
            rtb_Q_a[d] += rtP.G_Value[i_0 + b_k] * rtb_y_n_0[b_k_tmp + k];
            i_0 += 18;
          }
        }

        k += 18;
      }

      // Math: '<S63>/Transpose2' incorporates:
      //   Constant: '<S43>/H'

      k = 0;
      for (i = 0; i < 6; i++) {
        b_k = 0;
        for (i_0 = 0; i_0 < 18; i_0++) {
          rtb_Transpose2_l[i_0 + k] = rtP.H_Value[b_k + i];
          b_k += 6;
        }

        k += 18;
      }

      // End of Math: '<S63>/Transpose2'
      for (k = 0; k < 18; k++) {
        // Sum: '<S63>/Add' incorporates:
        //   Math: '<S63>/Transpose2'
        //   Product: '<S63>/Product1'

        i = 0;
        for (b_k = 0; b_k < 6; b_k++) {
          dwt = 0.0;
          i_0 = 0;
          for (b_k_tmp = 0; b_k_tmp < 18; b_k_tmp++) {
            dwt += rtb_y_n[i_0 + k] * rtb_Transpose2_l[b_k_tmp + i];
            i_0 += 18;
          }

          i_0 = i + k;
          rtb_Add_g[i_0] = rtb_N_f[i_0] + dwt;
          i += 18;
        }

        // End of Sum: '<S63>/Add'
      }

      for (k = 0; k < 6; k++) {
        for (i = 0; i < 18; i++) {
          // Product: '<S63>/Product2' incorporates:
          //   Constant: '<S43>/G'
          //   Sum: '<S63>/Add'

          i_0 = 18 * k + i;
          rtb_Product2_mz[i_0] = 0.0;
          for (b_k = 0; b_k < 18; b_k++) {
            rtb_Product2_mz[i_0] += rtP.G_Value[18 * b_k + i] * rtb_Add_g[18 * k
              + b_k];
          }

          // End of Product: '<S63>/Product2'
        }

        for (i = 0; i < 6; i++) {
          // Product: '<S63>/Product3' incorporates:
          //   Constant: '<S43>/H'
          //   Product: '<S63>/Product4'

          b_k = 6 * i + k;
          tmp_6[b_k] = 0.0;

          // Product: '<S63>/Product4'
          rtb_N_l[b_k] = 0.0;
          for (i_0 = 0; i_0 < 18; i_0++) {
            // Product: '<S63>/Product3' incorporates:
            //   Constant: '<S43>/H'
            //   Product: '<S63>/Product4'
            //   Sum: '<S63>/Add'

            b_k_tmp = 18 * i + i_0;
            tmp_6[b_k] += rtP.H_Value[6 * i_0 + k] * rtb_Add_g[b_k_tmp];

            // Product: '<S63>/Product4' incorporates:
            //   Math: '<S63>/Transpose'
            //   Math: '<S63>/Transpose2'

            rtb_N_l[b_k] += rtb_N_f[18 * k + i_0] * rtb_Transpose2_l[b_k_tmp];
          }
        }
      }

      // Sum: '<S63>/Add1' incorporates:
      //   MATLAB Function: '<S66>/ScalarExpansion'

      k = 0;
      for (i = 0; i < 6; i++) {
        b_k = 0;
        for (i_0 = 0; i_0 < 6; i_0++) {
          b_k_tmp = i_0 + k;

          // Outputs for Atomic SubSystem: '<S43>/ScalarExpansionR'
          rtb_A_n[b_k_tmp] = (rtb_R_d[b_k + i] + rtb_R_d[b_k_tmp]) / 2.0 +
            (tmp_6[b_k_tmp] + rtb_N_l[b_k_tmp]);

          // End of Outputs for SubSystem: '<S43>/ScalarExpansionR'
          b_k += 6;
        }

        k += 6;
      }

      // End of Sum: '<S63>/Add1'
      // End of Outputs for SubSystem: '<S43>/ReducedQRN'

      // Outputs for Atomic SubSystem: '<S43>/CalculatePL'
      // MATLAB Function: '<S45>/Discrete-Time KF - Calculate PLMZ' incorporates:
      //   Constant: '<S2>/Constant1'
      //   DataTypeConversion: '<S43>/DataTypeConversionEnable'
      //   Delay: '<S43>/MemoryP'
      //   DiscreteFilter: '<S2>/Discrete Filter1'
      //   Product: '<S63>/Product'
      //   Product: '<S63>/Product2'
      //   Sum: '<S63>/Add1'
      //
      //  See help of ctrlKalmanFilterDTCalculatePL.m
      // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S83>:1' 
      //    Copyright 2014 The MathWorks, Inc.
      // '<S83>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
      if (rtP.Constant1_Value_c != 0.0) {
        k = 0;
        for (i = 0; i < 6; i++) {
          b_k = 0;
          i_0 = 0;
          for (b_k_tmp = 0; b_k_tmp < 18; b_k_tmp++) {
            rtb_R_p_tmp_0 = b_k + i;
            rtb_Transpose2_l[b_k_tmp + k] = rtb_C_n[rtb_R_p_tmp_0];
            rtb_N_f[rtb_R_p_tmp_0] = 0.0;
            d = 0;
            for (int32_T i_1{0}; i_1 < 18; i_1++) {
              rtb_N_f[rtb_R_p_tmp_0] += rtb_C_n[d + i] *
                rtDW.MemoryP_DSTATE_h4[i_1 + i_0];
              d += 6;
            }

            b_k += 6;
            i_0 += 18;
          }

          k += 18;
        }

        for (k = 0; k < 6; k++) {
          i = 0;
          b_k = 0;
          for (i_0 = 0; i_0 < 6; i_0++) {
            dwt = 0.0;
            b_k_tmp = 0;
            for (d = 0; d < 18; d++) {
              dwt += rtb_N_f[b_k_tmp + k] * rtb_Transpose2_l[d + b_k];
              b_k_tmp += 6;
            }

            rtb_R_p_tmp_0 = i + k;
            rtb_R_d[rtb_R_p_tmp_0] = rtb_A_n[rtb_R_p_tmp_0] + dwt;
            i += 6;
            b_k += 18;
          }
        }

        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            i_0 = 18 * i + k;
            rtb_y_n_0[i_0] = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              rtb_y_n_0[i_0] += rtb_A_b[18 * b_k + k] * rtDW.MemoryP_DSTATE_h4
                [18 * i + b_k];
            }
          }

          for (i = 0; i < 6; i++) {
            dwt = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              dwt += rtb_y_n_0[18 * b_k + k] * rtb_Transpose2_l[18 * i + b_k];
            }

            i_0 = 18 * i + k;
            rtb_Add_g[i_0] = rtb_Product2_mz[i_0] + dwt;
          }
        }

        mrdiv_g(rtb_Add_g, rtb_R_d, rtb_N_f);
        k = 0;
        for (i = 0; i < 6; i++) {
          for (b_k = 0; b_k < 18; b_k++) {
            i_0 = b_k + k;
            rtb_Add_g[i_0] = 0.0;
            b_k_tmp = 0;
            for (d = 0; d < 18; d++) {
              rtb_Add_g[i_0] += rtDW.MemoryP_DSTATE_h4[b_k_tmp + b_k] *
                rtb_Transpose2_l[d + k];
              b_k_tmp += 18;
            }
          }

          k += 18;
        }

        mrdiv_g(rtb_Add_g, rtb_R_d, rtb_Transpose2_l);
        (void)std::memset(&b_I[0], 0, 324U * sizeof(int8_T));
        i = 0;
        for (k = 0; k < 18; k++) {
          b_I[i] = 1;
          i += 19;
        }

        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            dwt = 0.0;
            for (b_k = 0; b_k < 6; b_k++) {
              dwt += rtb_Transpose2_l[18 * b_k + k] * rtb_C_n[6 * i + b_k];
            }

            i_0 = 18 * i + k;
            rtb_y_n[i_0] = static_cast<real_T>(b_I[i_0]) - dwt;
          }

          for (i = 0; i < 18; i++) {
            i_0 = 18 * i + k;
            Abar[i_0] = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              Abar[i_0] += rtb_y_n[18 * b_k + k] * rtDW.MemoryP_DSTATE_h4[18 * i
                + b_k];
            }
          }

          for (i = 0; i < 6; i++) {
            rtb_R_p_tmp_0 = 18 * i + k;
            rtb_Add_g[rtb_R_p_tmp_0] = 0.0;
            for (b_k = 0; b_k < 6; b_k++) {
              rtb_Add_g[rtb_R_p_tmp_0] += rtb_Transpose2_l[18 * b_k + k] *
                rtb_A_n[6 * i + b_k];
            }
          }
        }

        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            i_0 = 18 * i + k;
            rtb_y_n_0[i_0] = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              rtb_y_n_0[i_0] += Abar[18 * b_k + k] * rtb_y_n[18 * b_k + i];
            }

            rtb_Transpose2_o[i_0] = 0.0;
            for (b_k = 0; b_k < 6; b_k++) {
              rtb_Transpose2_o[i_0] += rtb_Add_g[18 * b_k + k] *
                rtb_Transpose2_l[18 * b_k + i];
            }
          }
        }

        for (k = 0; k <= 322; k += 2) {
          tmp_8 = _mm_loadu_pd(&rtb_y_n_0[k]);
          tmp_9 = _mm_loadu_pd(&rtb_Transpose2_o[k]);
          (void)_mm_storeu_pd(&rtb_Z_d[k], _mm_add_pd(tmp_8, tmp_9));
        }

        mrdiv_g(rtb_Product2_mz, rtb_A_n, rtb_Transpose2_l);
        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            dwt = 0.0;
            for (b_k = 0; b_k < 6; b_k++) {
              dwt += rtb_Transpose2_l[18 * b_k + k] * rtb_C_n[6 * i + b_k];
            }

            i_0 = 18 * i + k;
            rtb_y_n[i_0] = rtb_A_b[i_0] - dwt;
          }

          for (i = 0; i < 18; i++) {
            i_0 = 18 * i + k;
            Abar[i_0] = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              Abar[i_0] += rtb_y_n[18 * b_k + k] * rtb_Z_d[18 * i + b_k];
            }
          }
        }

        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            dwt = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              dwt += Abar[18 * b_k + k] * rtb_y_n[18 * b_k + i];
            }

            i_0 = 18 * i + k;
            rtb_y_n_0[i_0] = rtb_Q_a[i_0] + dwt;
            rtb_Transpose2_o[i_0] = 0.0;
            for (b_k = 0; b_k < 6; b_k++) {
              rtb_Transpose2_o[i_0] += rtb_Transpose2_l[18 * b_k + k] *
                rtb_Product2_mz[18 * b_k + i];
            }
          }
        }

        for (k = 0; k <= 322; k += 2) {
          tmp_8 = _mm_loadu_pd(&rtb_y_n_0[k]);
          tmp_9 = _mm_loadu_pd(&rtb_Transpose2_o[k]);
          (void)_mm_storeu_pd(&rtb_y_n[k], _mm_sub_pd(tmp_8, tmp_9));
        }
      } else {
        (void)std::memset(&rtb_N_f[0], 0, 108U * sizeof(real_T));
        for (k = 0; k < 18; k++) {
          for (i = 0; i < 18; i++) {
            i_0 = 18 * i + k;
            rtb_y_n_0[i_0] = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              rtb_y_n_0[i_0] += rtb_A_b[18 * b_k + k] * rtDW.MemoryP_DSTATE_h4
                [18 * i + b_k];
            }
          }

          for (i = 0; i < 18; i++) {
            dwt = 0.0;
            for (b_k = 0; b_k < 18; b_k++) {
              dwt += rtb_y_n_0[18 * b_k + k] * rtb_A_b[18 * b_k + i];
            }

            d = 18 * i + k;
            rtb_y_n[d] = rtb_Q_a[d] + dwt;
          }
        }
      }

      // End of MATLAB Function: '<S45>/Discrete-Time KF - Calculate PLMZ'
      // End of Outputs for SubSystem: '<S43>/CalculatePL'

      // RandomNumber: '<S2>/excitation'
      //  Determine if the Square-Root algorithm was used
      // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S85>:1' 
      // '<S85>:1:4' if isSqrtUsed
      rtb_excitation[0] = rtDW.NextOutput[0];
      rtb_excitation[1] = rtDW.NextOutput[1];
      rtb_excitation[2] = rtDW.NextOutput[2];
      for (k = 0; k < 3; k++) {
        // DiscreteFilter: '<S2>/Discrete Filter1'
        i = k * 59;
        sigmoid_workspace_k_1 = rtb_excitation[k] / rtP.lpfDen;
        dwt = rtP.lpfNum[0] * sigmoid_workspace_k_1;
        d = 1;
        for (b_k = 0; b_k < 59; b_k++) {
          dwt += rtDW.DiscreteFilter1_states[i + b_k] * rtP.lpfNum[d];
          d++;
        }

        rtb_Sum_gu[k] = dwt;
        DiscreteFilter1_tmp[k] = sigmoid_workspace_k_1;
      }

      // Sum: '<S2>/Sum' incorporates:
      //   Gain: '<S13>/u_scale'
      //   Inport: '<Root>/excitation'
      //   Product: '<S2>/Product1'

      dwt = rtP.u_scale_Gain[0] * Sum2_c[0] + rtb_Sum_gu[0] * rtU.excitation;

      // Saturate: '<S2>/Saturation'
      if (dwt > rtP.Saturation_UpperSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_UpperSat;
      } else if (dwt < rtP.Saturation_LowerSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_LowerSat;
      } else {
        // no actions
      }

      // Sum: '<S12>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      rtb_Sum_gu[0] = dwt - rtU.u0[0];

      // End of Outputs for SubSystem: '<S1>/ampc'

      // Saturate: '<S2>/Saturation'
      rtb_excitation[0] = dwt;

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // Sum: '<S2>/Sum' incorporates:
      //   Gain: '<S13>/u_scale'
      //   Inport: '<Root>/excitation'
      //   Product: '<S2>/Product1'

      dwt = rtP.u_scale_Gain[1] * Sum2_c[1] + rtb_Sum_gu[1] * rtU.excitation;

      // Saturate: '<S2>/Saturation'
      if (dwt > rtP.Saturation_UpperSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_UpperSat;
      } else if (dwt < rtP.Saturation_LowerSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_LowerSat;
      } else {
        // no actions
      }

      // Sum: '<S12>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      rtb_Sum_gu[1] = dwt - rtU.u0[1];

      // End of Outputs for SubSystem: '<S1>/ampc'

      // Saturate: '<S2>/Saturation'
      rtb_excitation[1] = dwt;

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // Sum: '<S2>/Sum' incorporates:
      //   Gain: '<S13>/u_scale'
      //   Inport: '<Root>/excitation'
      //   Product: '<S2>/Product1'

      dwt = rtP.u_scale_Gain[2] * Sum2_c[2] + rtb_Sum_gu[2] * rtU.excitation;

      // Saturate: '<S2>/Saturation'
      if (dwt > rtP.Saturation_UpperSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_UpperSat;
      } else if (dwt < rtP.Saturation_LowerSat) {
        // Saturate: '<S2>/Saturation'
        dwt = rtP.Saturation_LowerSat;
      } else {
        // no actions
      }

      // Sum: '<S12>/Sum1' incorporates:
      //   Inport: '<Root>/u0'

      rtb_Sum_gu[2] = dwt - rtU.u0[2];

      // End of Outputs for SubSystem: '<S1>/ampc'

      // Saturate: '<S2>/Saturation'
      rtb_excitation[2] = dwt;

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // Outputs for Enabled SubSystem: '<S62>/MeasurementUpdate' incorporates:
      //   EnablePort: '<S86>/Enable'

      // DataTypeConversion: '<S43>/DataTypeConversionEnable' incorporates:
      //   Constant: '<S2>/Constant1'
      //   Constant: '<S2>/Constant13'
      //   Product: '<S86>/C[k]*xhat[k|k-1]'
      //   Product: '<S86>/D[k]*u[k]'

      if (rtP.Constant1_Value_c != 0.0) {
        rtDW.MeasurementUpdate_MODE_b = true;
        for (k = 0; k < 6; k++) {
          // Product: '<S86>/C[k]*xhat[k|k-1]' incorporates:
          //   Delay: '<S43>/MemoryX'

          rtb_C_jm[k] = 0.0;
          i = 0;
          for (b_k = 0; b_k < 18; b_k++) {
            rtb_C_jm[k] += rtb_C_n[i + k] * rtDW.MemoryX_DSTATE_d[b_k];
            i += 6;
          }

          // Product: '<S86>/D[k]*u[k]' incorporates:
          //   Product: '<S86>/C[k]*xhat[k|k-1]'

          tmp_0[k] = 0.0;
          tmp_0[k] += rtP.Constant13_Value[k] * rtb_Sum_gu[0];
          tmp_0[k] += rtP.Constant13_Value[k + 6] * rtb_Sum_gu[1];
          tmp_0[k] += rtP.Constant13_Value[k + 12] * rtb_Sum_gu[2];

          // Sum: '<S86>/Sum' incorporates:
          //   Constant: '<S2>/Constant13'
          //   Inport: '<Root>/y'
          //   Inport: '<Root>/y0'
          //   Product: '<S86>/C[k]*xhat[k|k-1]'
          //   Product: '<S86>/D[k]*u[k]'
          //   Sum: '<S12>/Sum6'
          //   Sum: '<S86>/Add1'

          tmp_1[k] = (rtU.y[k] - rtU.y0[k]) - (rtb_C_jm[k] + tmp_0[k]);
        }

        for (k = 0; k < 18; k++) {
          // Product: '<S86>/Product3'
          rtDW.Product3_n[k] = 0.0;
          i = 0;
          for (b_k = 0; b_k < 6; b_k++) {
            rtDW.Product3_n[k] += rtb_N_f[i + k] * tmp_1[b_k];
            i += 18;
          }
        }
      } else if (rtDW.MeasurementUpdate_MODE_b) {
        for (i = 0; i < 18; i++) {
          // Disable for Product: '<S86>/Product3' incorporates:
          //   Outport: '<S86>/L*(y[k]-yhat[k|k-1])'
          //
          rtDW.Product3_n[i] = rtP.Lykyhatkk1_Y0;
        }

        rtDW.MeasurementUpdate_MODE_b = false;
      } else {
        // no actions
      }

      // End of Outputs for SubSystem: '<S62>/MeasurementUpdate'
      for (k = 0; k < 6; k++) {
        // Product: '<S46>/Product' incorporates:
        //   Delay: '<S43>/MemoryX'

        rtb_C_jm[k] = 0.0;
        i = 0;
        for (b_k = 0; b_k < 18; b_k++) {
          rtb_C_jm[k] += rtb_C_n[i + k] * rtDW.MemoryX_DSTATE_d[b_k];
          i += 6;
        }

        // Product: '<S46>/Product1' incorporates:
        //   Constant: '<S2>/Constant13'
        //   Product: '<S46>/Product'

        tmp_0[k] = 0.0;
        tmp_0[k] += rtP.Constant13_Value[k] * rtb_Sum_gu[0];
        tmp_0[k] += rtP.Constant13_Value[k + 6] * rtb_Sum_gu[1];
        tmp_0[k] += rtP.Constant13_Value[k + 12] * rtb_Sum_gu[2];

        // Sum: '<S46>/Add1' incorporates:
        //   Product: '<S46>/Product'

        rtb_ywtT[k] = rtb_C_jm[k] + tmp_0[k];
      }

      // Update for UnitDelay: '<S13>/last_mv'
      rtDW.last_mv_DSTATE_m[0] = Sum2_c[0];
      rtDW.last_mv_DSTATE_m[1] = Sum2_c[1];
      rtDW.last_mv_DSTATE_m[2] = Sum2_c[2];

      // Update for Delay: '<S43>/MemoryX'
      rtDW.icLoad_k = false;
      for (k = 0; k < 18; k++) {
        // Product: '<S62>/B[k]*u[k]'
        rtb_B_m[k] = 0.0;
        rtb_B_m[k] += rtb_B_p[k] * rtb_Sum_gu[0];
        rtb_B_m[k] += rtb_B_p[k + 18] * rtb_Sum_gu[1];
        rtb_B_m[k] += rtb_B_p[k + 36] * rtb_Sum_gu[2];

        // Product: '<S62>/A[k]*xhat[k|k-1]' incorporates:
        //   Delay: '<S43>/MemoryX'
        //   Product: '<S62>/B[k]*u[k]'

        b_xoff[k] = 0.0;
        i = 0;
        for (b_k = 0; b_k < 18; b_k++) {
          b_xoff[k] += rtb_A_b[i + k] * rtDW.MemoryX_DSTATE_d[b_k];
          i += 18;
        }

        // End of Product: '<S62>/A[k]*xhat[k|k-1]'
      }

      // End of Outputs for SubSystem: '<S1>/ampc'
      for (k = 0; k <= 16; k += 2) {
        // Outputs for Function Call SubSystem: '<S1>/ampc'
        // Sum: '<S62>/Add'
        tmp_8 = _mm_loadu_pd(&rtb_B_m[k]);
        tmp_9 = _mm_loadu_pd(&b_xoff[k]);
        tmp_7 = _mm_loadu_pd(&rtDW.Product3_n[k]);

        // Update for Delay: '<S43>/MemoryX' incorporates:
        //   Sum: '<S62>/Add'

        (void)_mm_storeu_pd(&rtDW.MemoryX_DSTATE_d[k], _mm_add_pd(_mm_add_pd
          (tmp_8, tmp_9), tmp_7));

        // End of Outputs for SubSystem: '<S1>/ampc'
      }

      // Outputs for Function Call SubSystem: '<S1>/ampc'
      // Update for Delay: '<S43>/MemoryP'
      rtDW.icLoad_j = false;
      (void)std::memcpy(&rtDW.MemoryP_DSTATE_h4[0], &rtb_y_n[0], 324U * sizeof
                        (real_T));

      // Update for RandomNumber: '<S2>/excitation'
      rtDW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[0]) *
        rtP.excitation_StdDev[0] + rtP.excitation_Mean[0];
      rtDW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[1]) *
        rtP.excitation_StdDev[1] + rtP.excitation_Mean[1];
      rtDW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[2]) *
        rtP.excitation_StdDev[2] + rtP.excitation_Mean[2];
      for (k = 0; k < 3; k++) {
        // Update for DiscreteFilter: '<S2>/Discrete Filter1'
        i = k * 59;
        for (b_k = 0; b_k < 58; b_k++) {
          i_0 = i - b_k;
          rtDW.DiscreteFilter1_states[i_0 + 58] =
            rtDW.DiscreteFilter1_states[i_0 + 57];
        }

        rtDW.DiscreteFilter1_states[i] = DiscreteFilter1_tmp[k];

        // Outport: '<Root>/u'
        rtY.u[k] = rtb_excitation[k];
      }

      // End of Outputs for SubSystem: '<S1>/ampc'
      for (i = 0; i <= 4; i += 2) {
        // Outputs for Function Call SubSystem: '<S1>/ampc'
        // Sum: '<S12>/Sum3' incorporates:
        //   Outport: '<Root>/yhat'

        tmp_8 = _mm_loadu_pd(&rtb_ywtT[i]);

        // End of Outputs for SubSystem: '<S1>/ampc'

        // Outport: '<Root>/yhat' incorporates:
        //   Inport: '<Root>/y0'

        (void)_mm_storeu_pd(&rtY.yhat[i], _mm_add_pd(tmp_8, _mm_loadu_pd
          (&rtU.y0[i])));
      }

      // [u, ywt, currTraj] = gmpc(traj(:,waypt), currEv.r, y, ymax, umax, uwt, k_2); 
      //  if sig == 1
      //      [u, yhat(1:no)] = mpc1(r_, y__, [0;0;0], 0, u0, umax, uwt, iRST);
      //  elseif sig == 2
      //      [u, yhat(1:no)] = mpc2(r_, y__, [0;0;0], [0;0], u0, umax, uwt, iRST); 
      //  elseif sig == 3
      //      [u, yhat(1:no)] = mpc3(r_, y__, [0;0;0], [0;0], u0, umax, uwt, iRST); 
      //  end
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
    real_T Product3_a[8];
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

    // InitializeConditions for Memory: '<S13>/Memory'
    (void)std::memcpy(&rtDW.Memory_PreviousInput_j[0],
                      &rtP.Memory_InitialCondition[0], 246U * sizeof(boolean_T));

    // InitializeConditions for Delay: '<S43>/MemoryX'
    rtDW.icLoad_k = true;

    // InitializeConditions for Delay: '<S43>/MemoryP'
    rtDW.icLoad_j = true;

    // InitializeConditions for UnitDelay: '<S13>/last_mv'
    rtDW.last_mv_DSTATE_m[0] = rtP.last_mv_InitialCondition[0];

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

    rtDW.RandSeed[0] = tseed;
    rtDW.NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[0]) *
      rtP.excitation_StdDev[0] + rtP.excitation_Mean[0];

    // InitializeConditions for UnitDelay: '<S13>/last_mv'
    rtDW.last_mv_DSTATE_m[1] = rtP.last_mv_InitialCondition[1];

    // InitializeConditions for RandomNumber: '<S2>/excitation'
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

    rtDW.RandSeed[1] = tseed;
    rtDW.NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[1]) *
      rtP.excitation_StdDev[1] + rtP.excitation_Mean[1];

    // InitializeConditions for UnitDelay: '<S13>/last_mv'
    rtDW.last_mv_DSTATE_m[2] = rtP.last_mv_InitialCondition[2];

    // InitializeConditions for RandomNumber: '<S2>/excitation'
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

    rtDW.RandSeed[2] = tseed;
    rtDW.NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&rtDW.RandSeed[2]) *
      rtP.excitation_StdDev[2] + rtP.excitation_Mean[2];
    for (i = 0; i < 177; i++) {
      // InitializeConditions for DiscreteFilter: '<S2>/Discrete Filter1'
      rtDW.DiscreteFilter1_states[i] = rtP.DiscreteFilter1_InitialStates;
    }

    // SystemInitialize for Enabled SubSystem: '<S62>/MeasurementUpdate'
    for (i = 0; i < 18; i++) {
      // SystemInitialize for Product: '<S86>/Product3' incorporates:
      //   Outport: '<S86>/L*(y[k]-yhat[k|k-1])'

      rtDW.Product3_n[i] = rtP.Lykyhatkk1_Y0;
    }

    // End of SystemInitialize for SubSystem: '<S62>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc2'

    // SystemInitialize for Enabled SubSystem: '<S202>/MeasurementUpdate'
    MeasurementUpdate_Init(Product3_a, &rtP.MeasurementUpdate_j);

    // End of SystemInitialize for SubSystem: '<S202>/MeasurementUpdate'

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/mpc3'

    // SystemInitialize for Enabled SubSystem: '<S272>/MeasurementUpdate'
    MeasurementUpdate_Init(Product3_a, &rtP.MeasurementUpdate_c);

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
