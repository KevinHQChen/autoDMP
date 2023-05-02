//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: SupervisoryController.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.720
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 22:26:26 2023
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
#include <cmath>
#include <cstring>
#include <emmintrin.h>
#include "zero_crossing_types.h"
#include <stddef.h>
#include "solver_zc.h"

// Named constants for MATLAB Function: '<S36>/FixedHorizonOptimizer'
const real_T Wy{ 0.018316915599999997 };

const int32_T degrees{ 4 };

const int32_T nu{ 3 };

const int32_T p{ 20 };

const real_T yoff{ 100.0 };

// Named constants for MATLAB Function: '<S141>/FixedHorizonOptimizer'
const int32_T degrees_h{ 4 };

const int32_T nu_g{ 3 };

const int32_T ny{ 2 };

const int32_T p_e{ 20 };

// Named constants for MATLAB Function: '<S285>/FixedHorizonOptimizer'
const int32_T degrees_m{ 4 };

const int32_T nu_h{ 3 };

const int32_T ny_i{ 2 };

const int32_T p_i{ 20 };

// Named constants for Chart: '<Root>/SupervisoryController'
const uint8_T IN_HandleEvent{ 1U };

const uint8_T IN_NO_ACTIVE_CHILD{ 0U };

const uint8_T IN_RequestEvent{ 2U };

const uint8_T IN_State0{ 1U };

const uint8_T IN_State1{ 2U };

const uint8_T IN_State2{ 3U };

const int32_T event_simMeasAvail{ 1 };

const int32_T event_trueMeasAvail{ 0 };

#define NumBitsPerChar                 8U
#include "solver_zc.h"
#ifndef slZcHadEvent
#define slZcHadEvent(ev, zcsDir)       (((ev) & (zcsDir)) != 0x00 )
#endif

#ifndef slZcUnAliasEvents
#define slZcUnAliasEvents(evL, evR)    ((((slZcHadEvent((evL), (SL_ZCS_EVENT_N2Z)) && slZcHadEvent((evR), (SL_ZCS_EVENT_Z2P))) || (slZcHadEvent((evL), (SL_ZCS_EVENT_P2Z)) && slZcHadEvent((evR), (SL_ZCS_EVENT_Z2N)))) ? (SL_ZCS_EVENT_NUL) : (evR)))
#endif

extern real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u);
extern real_T rt_hypotd_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
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

// Function for MATLAB Function: '<S38>/RLS'
real_T SupervisoryController::xnrm2_g(int32_T n, const real_T x[5], int32_T ix0)
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

// Function for MATLAB Function: '<S38>/RLS'
real_T SupervisoryController::qrFactor(const real_T A[4], const real_T S[16],
  real_T Ns)
{
  real_T b_A[5];
  real_T b_S;
  real_T s;
  int32_T aoffset;
  for (int32_T i{0}; i < 4; i++) {
    aoffset = i << 2UL;
    b_A[i] = ((S[aoffset + 1] * A[1] + S[aoffset] * A[0]) + S[aoffset + 2] * A[2])
      + S[aoffset + 3] * A[3];
  }

  b_A[4] = Ns;
  b_S = b_A[0];
  s = xnrm2_g(4, b_A, 2);
  if (s != 0.0) {
    s = rt_hypotd_snf(b_A[0], s);
    if (b_A[0] >= 0.0) {
      s = -s;
    }

    if (std::abs(s) < 1.0020841800044864E-292) {
      aoffset = 0;
      do {
        aoffset++;
        b_A[1] *= 9.9792015476736E+291;
        b_A[2] *= 9.9792015476736E+291;
        b_A[3] *= 9.9792015476736E+291;
        b_A[4] *= 9.9792015476736E+291;
        s *= 9.9792015476736E+291;
        b_S *= 9.9792015476736E+291;
      } while ((std::abs(s) < 1.0020841800044864E-292) && (aoffset < 20));

      s = rt_hypotd_snf(b_S, xnrm2_g(4, b_A, 2));
      if (b_S >= 0.0) {
        s = -s;
      }

      for (int32_T i{0}; i < aoffset; i++) {
        s *= 1.0020841800044864E-292;
      }

      b_S = s;
    } else {
      b_S = s;
    }
  }

  return b_S;
}

// Function for MATLAB Function: '<S38>/RLS'
void SupervisoryController::trisolve_m(real_T A, real_T B_1[4])
{
  if (B_1[0] != 0.0) {
    B_1[0] /= A;
  }

  if (B_1[1] != 0.0) {
    B_1[1] /= A;
  }

  if (B_1[2] != 0.0) {
    B_1[2] /= A;
  }

  if (B_1[3] != 0.0) {
    B_1[3] /= A;
  }
}

// Function for MATLAB Function: '<S38>/RLS'
real_T SupervisoryController::xnrm2_gg(int32_T n, const real_T x[20], int32_T
  ix0)
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

// Function for MATLAB Function: '<S38>/RLS'
void SupervisoryController::xgemv_m(int32_T m, int32_T n, const real_T A[20],
  int32_T ia0, const real_T x[20], int32_T ix0, real_T y[4])
{
  if ((m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      (void)std::memset(&y[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    b = (n - 1) * 5 + ia0;
    for (int32_T b_iy{ia0}; b_iy <= b; b_iy += 5) {
      real_T c;
      int32_T d;
      int32_T iyend;
      c = 0.0;
      d = b_iy + m;
      for (iyend = b_iy; iyend < d; iyend++) {
        c += x[((ix0 + iyend) - b_iy) - 1] * A[iyend - 1];
      }

      iyend = div_nde_s32_floor(b_iy - ia0, 5);
      y[iyend] += c;
    }
  }
}

// Function for MATLAB Function: '<S38>/RLS'
void SupervisoryController::xgerc_n(int32_T m, int32_T n, real_T alpha1, int32_T
  ix0, const real_T y[4], real_T A[20], int32_T ia0)
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
        b = m + jA;
        for (int32_T ijA{jA}; ijA < b; ijA++) {
          A[ijA - 1] += A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 5;
    }
  }
}

// Function for MATLAB Function: '<S38>/RLS'
void SupervisoryController::sqrtMeasurementUpdate(real_T L[16], const real_T H[4],
  real_T a0, real_T K[4])
{
  real_T c_A[20];
  real_T A[16];
  real_T y[16];
  real_T B_0[4];
  real_T Pxy[4];
  real_T K_0;
  real_T K_idx_1;
  real_T K_idx_2;
  real_T Sy;
  real_T prodVal;
  int32_T aoffset;
  int32_T b_i;
  int32_T coffset;
  int32_T ii;
  for (b_i = 0; b_i < 4; b_i++) {
    Pxy[b_i] = 0.0;
    for (ii = 0; ii < 4; ii++) {
      coffset = (ii << 2UL) + b_i;
      A[coffset] = 0.0;
      A[coffset] += L[b_i] * L[ii];
      A[coffset] += L[b_i + 4] * L[ii + 4];
      A[coffset] += L[b_i + 8] * L[ii + 8];
      A[coffset] += L[b_i + 12] * L[ii + 12];
      if (std::isinf(A[coffset])) {
        prodVal = A[coffset];
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        A[coffset] = K_0 * 1.7976931348623157E+308;
      }

      prodVal = A[coffset] * H[ii];
      if (std::isinf(prodVal)) {
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        prodVal = K_0 * 1.7976931348623157E+308;
      }

      prodVal += Pxy[b_i];
      if (std::isinf(prodVal)) {
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        prodVal = K_0 * 1.7976931348623157E+308;
      }

      Pxy[b_i] = prodVal;
    }
  }

  prodVal = std::sqrt(a0);
  Sy = qrFactor(H, L, prodVal);
  B_0[0] = Pxy[0];
  B_0[1] = Pxy[1];
  B_0[2] = Pxy[2];
  B_0[3] = Pxy[3];
  trisolve_m(Sy, B_0);
  K[0] = B_0[0];
  K[1] = B_0[1];
  K[2] = B_0[2];
  K[3] = B_0[3];
  trisolve_m(Sy, K);
  K_0 = K[0];
  if (std::isinf(K[0])) {
    if (std::isnan(K[0])) {
      K_0 = (rtNaN);
    } else if (K[0] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[0] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  Sy = -K_0;
  K[0] = K_0;
  K_0 = K[1];
  if (std::isinf(K[1])) {
    if (std::isnan(K[1])) {
      K_0 = (rtNaN);
    } else if (K[1] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[1] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K_idx_1 = -K_0;
  K[1] = K_0;
  K_0 = K[2];
  if (std::isinf(K[2])) {
    if (std::isnan(K[2])) {
      K_0 = (rtNaN);
    } else if (K[2] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[2] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K_idx_2 = -K_0;
  K[2] = K_0;
  K_0 = K[3];
  if (std::isinf(K[3])) {
    if (std::isnan(K[3])) {
      K_0 = (rtNaN);
    } else if (K[3] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[3] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K[3] = K_0;
  b_i = 0;
  for (ii = 0; ii < 4; ii++) {
    A[b_i] = Sy * H[ii];
    A[b_i + 1] = K_idx_1 * H[ii];
    A[b_i + 2] = K_idx_2 * H[ii];
    A[b_i + 3] = -K_0 * H[ii];
    b_i += 4;
  }

  A[0]++;
  A[5]++;
  A[10]++;
  A[15]++;
  for (b_i = 0; b_i < 4; b_i++) {
    coffset = b_i << 2UL;
    for (ii = 0; ii < 4; ii++) {
      aoffset = ii << 2UL;
      y[coffset + ii] = ((L[aoffset + 1] * A[b_i + 4] + L[aoffset] * A[b_i]) +
                         L[aoffset + 2] * A[b_i + 8]) + L[aoffset + 3] * A[b_i +
        12];
    }
  }

  b_i = 0;
  coffset = 0;
  for (ii = 0; ii < 4; ii++) {
    c_A[b_i] = y[coffset];
    c_A[b_i + 1] = y[coffset + 1];
    c_A[b_i + 2] = y[coffset + 2];
    c_A[b_i + 3] = y[coffset + 3];
    c_A[b_i + 4] = K[ii] * prodVal;
    Pxy[ii] = 0.0;
    b_i += 5;
    coffset += 4;
  }

  for (b_i = 0; b_i < 4; b_i++) {
    int32_T lastv;
    int32_T scalarLB;
    ii = b_i * 5 + b_i;
    Sy = c_A[ii];
    lastv = ii + 2;
    B_0[b_i] = 0.0;
    prodVal = xnrm2_gg(4 - b_i, c_A, ii + 2);
    if (prodVal != 0.0) {
      prodVal = rt_hypotd_snf(c_A[ii], prodVal);
      if (c_A[ii] >= 0.0) {
        prodVal = -prodVal;
      }

      if (std::abs(prodVal) < 1.0020841800044864E-292) {
        __m128d tmp;
        int32_T scalarLB_tmp;
        int32_T vectorUB;
        coffset = 0;
        scalarLB = (ii - b_i) + 5;
        do {
          coffset++;
          scalarLB_tmp = (((((scalarLB - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB_tmp - 2;
          for (aoffset = lastv; aoffset <= vectorUB; aoffset += 2) {
            tmp = _mm_loadu_pd(&c_A[aoffset - 1]);
            (void)_mm_storeu_pd(&c_A[aoffset - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (aoffset = scalarLB_tmp; aoffset <= scalarLB; aoffset++) {
            c_A[aoffset - 1] *= 9.9792015476736E+291;
          }

          prodVal *= 9.9792015476736E+291;
          Sy *= 9.9792015476736E+291;
        } while ((std::abs(prodVal) < 1.0020841800044864E-292) && (coffset < 20));

        prodVal = rt_hypotd_snf(Sy, xnrm2_gg(4 - b_i, c_A, ii + 2));
        if (Sy >= 0.0) {
          prodVal = -prodVal;
        }

        B_0[b_i] = (prodVal - Sy) / prodVal;
        Sy = 1.0 / (Sy - prodVal);
        vectorUB = scalarLB_tmp - 2;
        for (aoffset = lastv; aoffset <= vectorUB; aoffset += 2) {
          tmp = _mm_loadu_pd(&c_A[aoffset - 1]);
          (void)_mm_storeu_pd(&c_A[aoffset - 1], _mm_mul_pd(tmp, _mm_set1_pd(Sy)));
        }

        for (aoffset = scalarLB_tmp; aoffset <= scalarLB; aoffset++) {
          c_A[aoffset - 1] *= Sy;
        }

        for (lastv = 0; lastv < coffset; lastv++) {
          prodVal *= 1.0020841800044864E-292;
        }

        Sy = prodVal;
      } else {
        int32_T vectorUB;
        B_0[b_i] = (prodVal - c_A[ii]) / prodVal;
        Sy = 1.0 / (c_A[ii] - prodVal);
        aoffset = (ii - b_i) + 5;
        scalarLB = (((((aoffset - ii) - 1) / 2) << 1UL) + ii) + 2;
        vectorUB = scalarLB - 2;
        for (coffset = lastv; coffset <= vectorUB; coffset += 2) {
          __m128d tmp;
          tmp = _mm_loadu_pd(&c_A[coffset - 1]);
          (void)_mm_storeu_pd(&c_A[coffset - 1], _mm_mul_pd(tmp, _mm_set1_pd(Sy)));
        }

        for (coffset = scalarLB; coffset <= aoffset; coffset++) {
          c_A[coffset - 1] *= Sy;
        }

        Sy = prodVal;
      }
    }

    c_A[ii] = Sy;
    if (b_i + 1 < 4) {
      prodVal = c_A[ii];
      c_A[ii] = 1.0;
      if (B_0[b_i] != 0.0) {
        boolean_T exitg2;
        lastv = 5 - b_i;
        coffset = (ii - b_i) + 4;
        while ((lastv > 0) && (c_A[coffset] == 0.0)) {
          lastv--;
          coffset--;
        }

        coffset = 3 - b_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (coffset > 0)) {
          int32_T exitg1;
          aoffset = ((coffset - 1) * 5 + ii) + 5;
          scalarLB = aoffset;
          do {
            exitg1 = 0;
            if (scalarLB + 1 <= aoffset + lastv) {
              if (c_A[scalarLB] != 0.0) {
                exitg1 = 1;
              } else {
                scalarLB++;
              }
            } else {
              coffset--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        coffset = 0;
      }

      if (lastv > 0) {
        xgemv_m(lastv, coffset, c_A, ii + 6, c_A, ii + 1, Pxy);
        xgerc_n(lastv, coffset, -B_0[b_i], ii + 1, Pxy, c_A, ii + 6);
      }

      c_A[ii] = prodVal;
    }
  }

  for (b_i = 0; b_i < 4; b_i++) {
    for (ii = 0; ii <= b_i; ii++) {
      A[ii + (b_i << 2UL)] = c_A[5 * b_i + ii];
    }

    for (ii = b_i + 2; ii < 5; ii++) {
      A[(ii + (b_i << 2UL)) - 1] = 0.0;
    }
  }

  for (ii = 0; ii < 4; ii++) {
    b_i = ii << 2UL;
    L[b_i] = A[ii];
    L[b_i + 1] = A[ii + 4];
    L[b_i + 2] = A[ii + 8];
    L[b_i + 3] = A[ii + 12];
    prodVal = L[b_i];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 1];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 1] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 2];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 2] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 3];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 3] = K_0 * 1.7976931348623157E+308;
    }
  }
}

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
void SupervisoryController::KWIKfactor(const real_T b_Ac[184], const int32_T iC
  [46], int32_T nA, const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
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
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 45];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 91];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 137];
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
void SupervisoryController::DropConstraint(int32_T kDrop, boolean_T iA[46],
  int32_T *nA, int32_T iC[46])
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
void SupervisoryController::qpkwik(const real_T b_Linv[16], const real_T b_Hinv
  [16], const real_T f[4], const real_T b_Ac[184], const real_T b[46], boolean_T
  iA[46], int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda[46],
  int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[46];
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
  int32_T iC[46];
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
  for (i = 0; i < 46; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 46; tmp++) {
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
            (void)std::memset(&iC[0], 0, 46U * sizeof(int32_T));
            for (i = 0; i < 46; i++) {
              iA[i] = false;
            }

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
              (void)std::memset(&iC[0], 0, 46U * sizeof(int32_T));
              for (i = 0; i < 46; i++) {
                iA[i] = false;
              }

              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 46U * sizeof(real_T));
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
      for (i = 0; i < 46; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 46] * x[1]);
          z[2] = std::abs(b_Ac[i + 92] * x[2]);
          z[3] = std::abs(b_Ac[i + 138] * x[3]);
          t = std::fmax(t, maximum(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 46] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 92] * x[2])
                   + b_Ac[i + 138] * x[3]) - b[i]) / t;
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
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 138])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 92])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 46])),
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
                    _mm_set1_pd(b_Ac[tmp + 138])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 92])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 46])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((b_D[iSave + 1] * b_Ac[tmp + 46] + b_D[iSave] *
                           b_Ac[tmp]) + b_D[iSave + 2] * b_Ac[tmp + 92]) +
                    b_D[iSave + 3] * b_Ac[tmp + 138];
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

              t = b_Ac[tmp + 46];
              cVal_tmp = b_Ac[tmp + 92];
              cVal_tmp_0 = b_Ac[tmp + 138];
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
                  if ((iC_0 <= 46) && (lambda[iC_0 - 1] < 0.0)) {
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
              for (tmp = 0; tmp < 46; tmp++) {
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

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
void SupervisoryController::mpcblock_optimizer(const real_T rseq[20], const
  real_T vseq[21], const real_T x[2], const real_T old_u[3], const boolean_T iA
  [46], const real_T b_Mlim[46], real_T b_Mx[92], real_T b_Mu1[138], real_T
  b_Mv[966], const real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16],
  real_T b_Ac[184], const real_T b_Wdu[3], const real_T b_Jm[180], const real_T
  b_Wu[3], const real_T b_I1[180], const real_T b_A[4], const real_T Bu[126],
  const real_T Bv[42], const real_T b_C[2], const real_T Dv[21], const int32_T
  b_Mrows[46], real_T u[3], real_T useq[63], real_T *status, boolean_T iAout[46],
  DW_State0controlLawAMPC0 *localDW)
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
  __m128d tmp_0;
  real_T b_Hv[420];
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T b_Kv[63];
  real_T Sum_0[60];
  real_T WySuJm[60];
  real_T b_Su1[60];
  real_T b_Mlim_0[46];
  real_T b_Mlim_1[46];
  real_T b_Mu1_0[46];
  real_T b_Sx[40];
  real_T CA[21];
  real_T L[16];
  real_T b_B[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T b_Kx[6];
  real_T varargin_1[4];
  real_T zopt[4];
  real_T Sum[3];
  real_T CA_idx_0;
  real_T CA_idx_1;
  real_T Sum_1;
  int32_T I2Jm_tmp;
  int32_T Tries;
  int32_T a_tmp;
  int32_T b_SuJm_tmp;
  int32_T i;
  int32_T i1;
  int32_T kidx;
  int16_T ixw;
  int8_T a[3600];
  int8_T b[16];
  boolean_T exitg1;
  boolean_T guard1{ false };

  (void)std::memset(&useq[0], 0, 63U * sizeof(real_T));
  for (i = 0; i < 46; i++) {
    iAout[i] = false;
  }

  CA_idx_0 = b_C[0] * b_A[0] + b_C[1] * b_A[1];
  CA_idx_1 = b_C[0] * b_A[2] + b_C[1] * b_A[3];
  i = 0;
  for (kidx = 0; kidx < 3; kidx++) {
    Sum[kidx] = 0.0;
    Sum[kidx] += Bu[i] * b_C[0];
    Sum[kidx] += Bu[i + 1] * b_C[1];
    i += 2;
  }

  b_Hv[0] = b_C[0] * Bv[0] + b_C[1] * Bv[1];
  b_Hv[20] = Dv[0];
  i = 0;
  for (kidx = 0; kidx < 19; kidx++) {
    b_Hv[i + 40] = 0.0;
    i += 20;
  }

  i = 0;
  for (kidx = 0; kidx < 21; kidx++) {
    (void)std::memset(&b_Hv[i + 1], 0, 19U * sizeof(real_T));
    i += 20;
  }

  b_Sx[0] = CA_idx_0;
  b_Sx[20] = CA_idx_1;
  b_Su1[0] = Sum[0];
  b_Su1[20] = Sum[1];
  b_Su1[40] = Sum[2];
  for (i = 0; i < 19; i++) {
    b_Sx[i + 1] = 0.0;
    b_Sx[i + 21] = 0.0;
    b_Su1[i + 1] = 0.0;
    b_Su1[i + 21] = 0.0;
    b_Su1[i + 41] = 0.0;
  }

  localDW->Su[0] = Sum[0];
  localDW->Su[20] = Sum[1];
  localDW->Su[40] = Sum[2];
  i = 0;
  for (kidx = 0; kidx < 57; kidx++) {
    localDW->Su[i + 60] = 0.0;
    i += 20;
  }

  i = 0;
  for (kidx = 0; kidx < 60; kidx++) {
    (void)std::memset(&localDW->Su[i + 1], 0, 19U * sizeof(real_T));
    i += 20;
  }

  for (kidx = 0; kidx < 19; kidx++) {
    for (i = 0; i < 3; i++) {
      Tries = i << 1UL;
      Sum_1 = (Bu[Tries + 1] * CA_idx_1 + Bu[Tries] * CA_idx_0) + Sum[i];
      b_Su1[(kidx + 20 * i) + 1] = Sum_1;
      Sum_0[i] = Sum_1;
      Sum[i] = Sum_1;
    }

    for (i = 0; i < 57; i++) {
      Sum_0[i + 3] = localDW->Su[20 * i + kidx];
    }

    for (i = 0; i < 60; i++) {
      localDW->Su[(kidx + 20 * i) + 1] = Sum_0[i];
    }

    CA[0] = CA_idx_0 * Bv[0] + CA_idx_1 * Bv[1];
    for (i = 0; i < 20; i++) {
      CA[i + 1] = b_Hv[20 * i + kidx];
    }

    for (i = 0; i < 21; i++) {
      b_Hv[(kidx + 20 * i) + 1] = CA[i];
    }

    Sum_1 = CA_idx_0 * b_A[0] + CA_idx_1 * b_A[1];
    CA_idx_1 = CA_idx_0 * b_A[2] + CA_idx_1 * b_A[3];
    CA_idx_0 = Sum_1;
    b_Sx[kidx + 1] = Sum_1;
    b_Sx[kidx + 21] = CA_idx_1;
  }

  i = 0;
  kidx = 0;
  for (Tries = 0; Tries < 3; Tries++) {
    for (i1 = 0; i1 < 20; i1++) {
      b_SuJm_tmp = i1 + i;
      Sum_0[b_SuJm_tmp] = 0.0;
      a_tmp = 0;
      for (I2Jm_tmp = 0; I2Jm_tmp < 60; I2Jm_tmp++) {
        Sum_0[b_SuJm_tmp] += localDW->Su[a_tmp + i1] * b_Jm[I2Jm_tmp + kidx];
        a_tmp += 20;
      }
    }

    i += 20;
    kidx += 60;
  }

  if (b_Mrows[0] > 0) {
    kidx = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (kidx < 46)) {
      if (b_Mrows[kidx] <= 20) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = -Sum_0[Tries - 1];
        b_Ac[kidx + 46] = -Sum_0[Tries + 19];
        b_Ac[kidx + 92] = -Sum_0[Tries + 39];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = -b_Sx[Tries - 1];
        b_Mx[kidx + 46] = -b_Sx[Tries + 19];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = -b_Su1[Tries - 1];
        b_Mu1[kidx + 46] = -b_Su1[Tries + 19];
        b_Mu1[kidx + 92] = -b_Su1[Tries + 39];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 46 * i] = -b_Hv[(20 * i + Tries) - 1];
        }

        kidx++;
      } else if (b_Mrows[kidx] <= 40) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = Sum_0[Tries - 21];
        b_Ac[kidx + 46] = Sum_0[Tries - 1];
        b_Ac[kidx + 92] = Sum_0[Tries + 19];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = b_Sx[Tries - 21];
        b_Mx[kidx + 46] = b_Sx[Tries - 1];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = b_Su1[Tries - 21];
        b_Mu1[kidx + 46] = b_Su1[Tries - 1];
        b_Mu1[kidx + 92] = b_Su1[Tries + 19];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 46 * i] = b_Hv[(20 * i + Tries) - 21];
        }

        kidx++;
      } else {
        exitg1 = true;
      }
    }
  }

  (void)std::memset(&b_B[0], 0, 9U * sizeof(real_T));
  b_B[0] = 1.0;
  b_B[4] = 1.0;
  b_B[8] = 1.0;
  kidx = -1;
  for (Tries = 0; Tries < 20; Tries++) {
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 20; i1++) {
        a_tmp = static_cast<int32_T>(c_A[20 * Tries + i1]);
        a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i]) *
          a_tmp);
        a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 1]) *
          a_tmp);
        a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 2]) *
          a_tmp);
        kidx += 3;
      }
    }
  }

  i = 0;
  for (kidx = 0; kidx < 3; kidx++) {
    for (Tries = 0; Tries < 60; Tries++) {
      I2Jm_tmp = Tries + i;
      I2Jm[I2Jm_tmp] = 0.0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 60; a_tmp++) {
        I2Jm[I2Jm_tmp] += static_cast<real_T>(a[i1 + Tries]) * b_Jm[a_tmp + i];
        i1 += 60;
      }
    }

    i += 60;
  }

  for (kidx = 0; kidx <= 18; kidx += 2) {
    tmp = _mm_loadu_pd(&Sum_0[kidx]);
    tmp_0 = _mm_set1_pd(Wy);
    (void)_mm_storeu_pd(&WySuJm[kidx], _mm_mul_pd(tmp_0, tmp));
    tmp = _mm_loadu_pd(&Sum_0[kidx + 20]);
    (void)_mm_storeu_pd(&WySuJm[kidx + 20], _mm_mul_pd(tmp, tmp_0));
    tmp = _mm_loadu_pd(&Sum_0[kidx + 40]);
    (void)_mm_storeu_pd(&WySuJm[kidx + 40], _mm_mul_pd(tmp, tmp_0));
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    CA_idx_0 = b_Wu[ixw - 1];
    WuI2Jm[kidx] = CA_idx_0 * I2Jm[kidx];
    WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * CA_idx_0;
    WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * CA_idx_0;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    CA_idx_0 = b_Wdu[ixw - 1];
    WduJm[kidx] = CA_idx_0 * b_Jm[kidx];
    WduJm[kidx + 60] = b_Jm[kidx + 60] * CA_idx_0;
    WduJm[kidx + 120] = b_Jm[kidx + 120] * CA_idx_0;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  for (i = 0; i < 3; i++) {
    for (kidx = 0; kidx < 3; kidx++) {
      b_SuJm_tmp = 3 * kidx + i;
      b_B[b_SuJm_tmp] = 0.0;
      for (Tries = 0; Tries < 20; Tries++) {
        b_B[b_SuJm_tmp] += Sum_0[20 * i + Tries] * WySuJm[20 * kidx + Tries];
      }

      b_Jm_0[b_SuJm_tmp] = 0.0;
      CA_idx_0 = 0.0;
      for (Tries = 0; Tries < 60; Tries++) {
        i1 = 60 * i + Tries;
        a_tmp = 60 * kidx + Tries;
        CA_idx_0 += I2Jm[i1] * WuI2Jm[a_tmp];
        b_Jm_0[b_SuJm_tmp] += b_Jm[i1] * WduJm[a_tmp];
      }

      b_H[i + (kidx << 2UL)] = (b_B[b_SuJm_tmp] + b_Jm_0[b_SuJm_tmp]) + CA_idx_0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (kidx = 0; kidx < 3; kidx++) {
      i1 = 3 * kidx + i;
      b_Jm_0[i1] = 0.0;
      for (Tries = 0; Tries < 20; Tries++) {
        b_Jm_0[i1] += b_Su1[20 * i + Tries] * WySuJm[20 * kidx + Tries];
      }

      b_I1_0[i1] = 0.0;
      for (Tries = 0; Tries < 60; Tries++) {
        b_I1_0[i1] += b_I1[60 * i + Tries] * WuI2Jm[60 * kidx + Tries];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
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
  for (kidx = 0; kidx < 2; kidx++) {
    Tries = 0;
    i1 = 0;
    for (a_tmp = 0; a_tmp < 3; a_tmp++) {
      b_SuJm_tmp = Tries + kidx;
      b_Kx[b_SuJm_tmp] = 0.0;
      for (I2Jm_tmp = 0; I2Jm_tmp < 20; I2Jm_tmp++) {
        b_Kx[b_SuJm_tmp] += b_Sx[I2Jm_tmp + i] * WySuJm[I2Jm_tmp + i1];
      }

      Tries += 2;
      i1 += 20;
    }

    i += 20;
  }

  i = 0;
  for (kidx = 0; kidx < 21; kidx++) {
    Tries = 0;
    i1 = 0;
    for (a_tmp = 0; a_tmp < 3; a_tmp++) {
      b_SuJm_tmp = Tries + kidx;
      b_Kv[b_SuJm_tmp] = 0.0;
      for (I2Jm_tmp = 0; I2Jm_tmp < 20; I2Jm_tmp++) {
        b_Kv[b_SuJm_tmp] += b_Hv[I2Jm_tmp + i] * WySuJm[I2Jm_tmp + i1];
      }

      Tries += 21;
      i1 += 20;
    }

    i += 20;
  }

  for (i = 0; i <= 58; i += 2) {
    tmp = _mm_loadu_pd(&WySuJm[i]);
    (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  kidx = 0;
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
    boolean_T exitg2;
    CA_idx_0 = 0.0;
    Tries = 0;
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (Tries < 4)) {
      CA_idx_1 = ((std::abs(b_H[Tries + 4]) + std::abs(b_H[Tries])) + std::abs
                  (b_H[Tries + 8])) + std::abs(b_H[Tries + 12]);
      if (std::isnan(CA_idx_1)) {
        CA_idx_0 = (rtNaN);
        exitg2 = true;
      } else {
        if (CA_idx_1 > CA_idx_0) {
          CA_idx_0 = CA_idx_1;
        }

        Tries++;
      }
    }

    if (CA_idx_0 >= 1.0E+10) {
      kidx = 2;
    } else {
      Tries = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (Tries <= 4)) {
        boolean_T guard2{ false };

        CA_idx_0 = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
          1.4901161193847656E-7;
        for (i = 0; i < 16; i++) {
          b[i] = 0;
        }

        b[0] = 1;
        b[5] = 1;
        b[10] = 1;
        b[15] = 1;
        for (i = 0; i < 16; i++) {
          b_H[i] += CA_idx_0 * static_cast<real_T>(b[i]);
          L[i] = b_H[i];
        }

        kidx = xpotrf(L);
        guard2 = false;
        if (kidx == 0) {
          varargin_1[0] = L[0];
          varargin_1[1] = L[5];
          varargin_1[2] = L[10];
          varargin_1[3] = L[15];
          if (minimum(varargin_1) > 1.4901161193847656E-7) {
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
    for (i = 0; i < 16; i++) {
      b[i] = 0;
    }

    b[0] = 1;
    b[5] = 1;
    b[10] = 1;
    b[15] = 1;
    i = 0;
    for (kidx = 0; kidx < 4; kidx++) {
      b_H[i] = static_cast<real_T>(b[i]);
      b_H[i + 1] = static_cast<real_T>(b[i + 1]);
      b_H[i + 2] = static_cast<real_T>(b[i + 2]);
      b_H[i + 3] = static_cast<real_T>(b[i + 3]);
      i += 4;
    }

    trisolve(L, b_H);
    varargin_1[0] = 0.0;
    varargin_1[1] = 0.0;
    varargin_1[2] = 0.0;
    varargin_1[3] = 0.0;
    for (kidx = 0; kidx < 3; kidx++) {
      CA_idx_0 = 0.0;
      for (i = 0; i < 20; i++) {
        CA_idx_0 += WySuJm[20 * kidx + i] * rseq[i];
      }

      CA_idx_1 = 0.0;
      for (i = 0; i < 21; i++) {
        CA_idx_1 += b_Kv[21 * kidx + i] * vseq[i];
      }

      Sum_1 = 0.0;
      for (i = 0; i < 60; i++) {
        Sum_1 += WuI2Jm[60 * kidx + i] * b_utarget[i];
      }

      i = kidx << 1UL;
      varargin_1[kidx] = ((((b_B[3 * kidx + 1] * old_u[1] + b_B[3 * kidx] *
        old_u[0]) + b_B[3 * kidx + 2] * old_u[2]) + ((b_Kx[i + 1] * x[1] +
        b_Kx[i] * x[0]) + CA_idx_0)) + CA_idx_1) + Sum_1;
    }

    for (i = 0; i < 46; i++) {
      iAout[i] = iA[i];
      b_Mlim_0[i] = (b_Mx[i + 46] * x[1] + b_Mx[i] * x[0]) + b_Mlim[i];
      b_Mu1_0[i] = 0.0;
      b_Mu1_0[i] += b_Mu1[i] * old_u[0];
      b_Mu1_0[i] += b_Mu1[i + 46] * old_u[1];
      b_Mu1_0[i] += b_Mu1[i + 92] * old_u[2];
    }

    i = 0;
    for (kidx = 0; kidx < 4; kidx++) {
      Tries = 0;
      for (i1 = 0; i1 < 4; i1++) {
        a_tmp = Tries + kidx;
        L[a_tmp] = 0.0;
        L[a_tmp] += b_H[i] * b_H[Tries];
        L[a_tmp] += b_H[i + 1] * b_H[Tries + 1];
        L[a_tmp] += b_H[i + 2] * b_H[Tries + 2];
        L[a_tmp] += b_H[i + 3] * b_H[Tries + 3];
        Tries += 4;
      }

      i += 4;
    }

    for (i = 0; i < 46; i++) {
      CA_idx_0 = 0.0;
      kidx = 0;
      for (Tries = 0; Tries < 21; Tries++) {
        CA_idx_0 += b_Mv[kidx + i] * vseq[Tries];
        kidx += 46;
      }

      b_Mlim_1[i] = -((b_Mlim_0[i] + b_Mu1_0[i]) + CA_idx_0);
    }

    qpkwik(b_H, L, varargin_1, b_Ac, b_Mlim_1, iAout, 200, 1.0E-6, zopt,
           b_Mlim_0, &kidx);
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

// System initialize for function-call system: '<S1>/State0.controlLaw.AMPC0'
void SupervisoryController::State0controlLawAMPC0_Init(real_T rty_u[3], real_T
  *rty_yhat, DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0 *localP)
{
  int32_T i;
  int32_T t_tmp;
  uint32_T tseed;

  // InitializeConditions for Memory: '<S8>/Memory'
  for (i = 0; i < 46; i++) {
    localDW->Memory_PreviousInput[i] = localP->Memory_InitialCondition[i];
  }

  // End of InitializeConditions for Memory: '<S8>/Memory'

  // InitializeConditions for UnitDelay: '<S6>/Unit Delay3'
  localDW->UnitDelay3_DSTATE = localP->UnitDelay3_InitialCondition;

  // InitializeConditions for Delay: '<S38>/delayBuffery'
  localDW->icLoad = true;

  // InitializeConditions for Delay: '<S38>/delayBufferH'
  localDW->icLoad_n = true;

  // InitializeConditions for Delay: '<S38>/delayTheta'
  localDW->icLoad_m = true;

  // InitializeConditions for Delay: '<S38>/delayL'
  localDW->icLoad_e = true;

  // InitializeConditions for Delay: '<S64>/MemoryX'
  localDW->icLoad_k = true;

  // InitializeConditions for Delay: '<S64>/MemoryP'
  localDW->icLoad_f = true;

  // InitializeConditions for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[0] = localP->last_mv_InitialCondition[0];

  // InitializeConditions for UnitDelay: '<S6>/Unit Delay2'
  localDW->UnitDelay2_DSTATE[0] = localP->UnitDelay2_InitialCondition;

  // InitializeConditions for RandomNumber: '<S2>/Measurement Noise'
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  t_tmp = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) & 32768U));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[0] = tseed;
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];

  // InitializeConditions for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[1] = localP->last_mv_InitialCondition[1];

  // InitializeConditions for UnitDelay: '<S6>/Unit Delay2'
  localDW->UnitDelay2_DSTATE[1] = localP->UnitDelay2_InitialCondition;

  // InitializeConditions for RandomNumber: '<S2>/Measurement Noise'
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[1] = tseed;
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];

  // InitializeConditions for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[2] = localP->last_mv_InitialCondition[2];

  // InitializeConditions for UnitDelay: '<S6>/Unit Delay2'
  localDW->UnitDelay2_DSTATE[2] = localP->UnitDelay2_InitialCondition;

  // InitializeConditions for RandomNumber: '<S2>/Measurement Noise'
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[2] = tseed;
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];

  // SystemInitialize for Enabled SubSystem: '<S83>/MeasurementUpdate'
  // SystemInitialize for Product: '<S107>/Product3' incorporates:
  //   Outport: '<S107>/L*(y[k]-yhat[k|k-1])'

  localDW->Product3[0] = localP->Lykyhatkk1_Y0;
  localDW->Product3[1] = localP->Lykyhatkk1_Y0;

  // End of SystemInitialize for SubSystem: '<S83>/MeasurementUpdate'

  // SystemInitialize for Outport: '<S2>/u'
  rty_u[0] = localP->u_Y0;
  rty_u[1] = localP->u_Y0;
  rty_u[2] = localP->u_Y0;

  // SystemInitialize for Outport: '<S2>/yhat'
  *rty_yhat = localP->yhat_Y0;
}

// Disable for function-call system: '<S1>/State0.controlLaw.AMPC0'
void SupervisoryController::State0controlLawAMPC0_Disable
  (DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0 *localP)
{
  // Disable for Enabled SubSystem: '<S83>/MeasurementUpdate'
  if (localDW->MeasurementUpdate_MODE) {
    // Disable for Product: '<S107>/Product3' incorporates:
    //   Outport: '<S107>/L*(y[k]-yhat[k|k-1])'

    localDW->Product3[0] = localP->Lykyhatkk1_Y0;
    localDW->Product3[1] = localP->Lykyhatkk1_Y0;
    localDW->MeasurementUpdate_MODE = false;
  }

  // End of Disable for SubSystem: '<S83>/MeasurementUpdate'
}

// Output and update for function-call system: '<S1>/State0.controlLaw.AMPC0'
void SupervisoryController::State0controlLawAMPC0(real_T rtu_r, real_T rtu_y,
  real_T rtu_y0, const real_T rtu_u0[3], real_T rtu_paramEst, real_T rty_u[3],
  real_T *rty_yhat, DW_State0controlLawAMPC0 *localDW, P_State0controlLawAMPC0
  *localP, P *rtP, ZCE_State0controlLawAMPC0 *localZCE)
{
  static const real_T h[184]{ -0.7084943675087555, -1.4170115739006484,
    -2.1255516199119082, -2.8341145062787887, -3.5427002337375666,
    -4.2513088030245427, -4.9599402148760419, -5.6685944700284132,
    -6.3772715692180277, -7.0859715131812822, -7.7946943026545954,
    -8.5034399383744113, -9.2122084210771966, -9.92099975149944,
    -10.62981393037766, -11.338650958448392, -12.047510836448199,
    -12.756393565113667, -13.465299145181405, -14.174227577388047,
    0.7084943675087555, 1.4170115739006484, 2.1255516199119082,
    2.8341145062787887, 3.5427002337375666, 4.2513088030245427,
    4.9599402148760419, 5.6685944700284132, 6.3772715692180277,
    7.0859715131812822, 7.7946943026545954, 8.5034399383744113,
    9.2122084210771966, 9.92099975149944, 10.62981393037766, 11.338650958448392,
    12.047510836448199, 12.756393565113667, 13.465299145181405,
    14.174227577388047, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.34693844974269455,
    0.69388808332406571, 1.0408489011046336, 1.3878209034449296,
    1.7348040907054971, 2.0817984632468911, 2.428804021429678,
    2.7758207656144358, 3.1228486961617543, 3.4698878134322348,
    3.8169381177864903, 4.1639996095851455, 4.511072289188836, 4.858156156958211,
    5.2052512132539288, 5.5523574584366608, 5.89947489286709, 6.24660351690591,
    6.5937433309138278, 6.9408943352515609, -0.34693844974269455,
    -0.69388808332406571, -1.0408489011046336, -1.3878209034449296,
    -1.7348040907054971, -2.0817984632468911, -2.428804021429678,
    -2.7758207656144358, -3.1228486961617543, -3.4698878134322348,
    -3.8169381177864903, -4.1639996095851455, -4.511072289188836,
    -4.858156156958211, -5.2052512132539288, -5.5523574584366608,
    -5.89947489286709, -6.24660351690591, -6.5937433309138278,
    -6.9408943352515609, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.3458694000362138,
    0.69174994944943113, 1.0376416485990609, 1.3835444978445242,
    1.7294584975452529, 2.0753836480606909, 2.421319949750294,
    2.7672674029735287, 3.1132260080898742, 3.45919576545882, 3.8051766754398688,
    4.1511687383925331, 4.4971719546763387, 4.8431863246508211,
    5.1892118486755292, 5.5352485271100225, 5.8812963603138728,
    6.227355348646662, 6.5734254924679858, 6.9195067921374491,
    -0.3458694000362138, -0.69174994944943113, -1.0376416485990609,
    -1.3835444978445242, -1.7294584975452529, -2.0753836480606909,
    -2.421319949750294, -2.7672674029735287, -3.1132260080898742,
    -3.45919576545882, -3.8051766754398688, -4.1511687383925331,
    -4.4971719546763387, -4.8431863246508211, -5.1892118486755292,
    -5.5352485271100225, -5.8812963603138728, -6.227355348646662,
    -6.5734254924679858, -6.9195067921374491, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0 };

  static const real_T j[180]{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
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

  static const real_T k[180]{ 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
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

  static const real_T f[138]{ -0.7084943675087555, -1.4170115739006484,
    -2.1255516199119082, -2.8341145062787887, -3.5427002337375666,
    -4.2513088030245427, -4.9599402148760419, -5.6685944700284132,
    -6.3772715692180277, -7.0859715131812822, -7.7946943026545954,
    -8.5034399383744113, -9.2122084210771966, -9.92099975149944,
    -10.62981393037766, -11.338650958448392, -12.047510836448199,
    -12.756393565113667, -13.465299145181405, -14.174227577388047,
    0.7084943675087555, 1.4170115739006484, 2.1255516199119082,
    2.8341145062787887, 3.5427002337375666, 4.2513088030245427,
    4.9599402148760419, 5.6685944700284132, 6.3772715692180277,
    7.0859715131812822, 7.7946943026545954, 8.5034399383744113,
    9.2122084210771966, 9.92099975149944, 10.62981393037766, 11.338650958448392,
    12.047510836448199, 12.756393565113667, 13.465299145181405,
    14.174227577388047, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, 0.34693844974269455,
    0.69388808332406571, 1.0408489011046336, 1.3878209034449296,
    1.7348040907054971, 2.0817984632468911, 2.428804021429678,
    2.7758207656144358, 3.1228486961617543, 3.4698878134322348,
    3.8169381177864903, 4.1639996095851455, 4.511072289188836, 4.858156156958211,
    5.2052512132539288, 5.5523574584366608, 5.89947489286709, 6.24660351690591,
    6.5937433309138278, 6.9408943352515609, -0.34693844974269455,
    -0.69388808332406571, -1.0408489011046336, -1.3878209034449296,
    -1.7348040907054971, -2.0817984632468911, -2.428804021429678,
    -2.7758207656144358, -3.1228486961617543, -3.4698878134322348,
    -3.8169381177864903, -4.1639996095851455, -4.511072289188836,
    -4.858156156958211, -5.2052512132539288, -5.5523574584366608,
    -5.89947489286709, -6.24660351690591, -6.5937433309138278,
    -6.9408943352515609, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.3458694000362138,
    0.69174994944943113, 1.0376416485990609, 1.3835444978445242,
    1.7294584975452529, 2.0753836480606909, 2.421319949750294,
    2.7672674029735287, 3.1132260080898742, 3.45919576545882, 3.8051766754398688,
    4.1511687383925331, 4.4971719546763387, 4.8431863246508211,
    5.1892118486755292, 5.5352485271100225, 5.8812963603138728,
    6.227355348646662, 6.5734254924679858, 6.9195067921374491,
    -0.3458694000362138, -0.69174994944943113, -1.0376416485990609,
    -1.3835444978445242, -1.7294584975452529, -2.0753836480606909,
    -2.421319949750294, -2.7672674029735287, -3.1132260080898742,
    -3.45919576545882, -3.8051766754398688, -4.1511687383925331,
    -4.4971719546763387, -4.8431863246508211, -5.1892118486755292,
    -5.5352485271100225, -5.8812963603138728, -6.227355348646662,
    -6.5734254924679858, -6.9195067921374491, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

  static const real_T e[92]{ -1.000032235800572, -1.0000644726402907,
    -1.0000967105191898, -1.0001289494373029, -1.0001611893946631,
    -1.0001934303913043, -1.0002256724272598, -1.000257915502563,
    -1.0002901596172478, -1.0003224047713473, -1.0003546509648951,
    -1.0003868981979249, -1.0004191464704701, -1.0004513957825643,
    -1.0004836461342406, -1.0005158975255328, -1.0005481499564746,
    -1.0005804034270993, -1.0006126579374404, -1.0006449134875315,
    1.000032235800572, 1.0000644726402907, 1.0000967105191898,
    1.0001289494373029, 1.0001611893946631, 1.0001934303913043,
    1.0002256724272598, 1.000257915502563, 1.0002901596172478,
    1.0003224047713473, 1.0003546509648951, 1.0003868981979249,
    1.0004191464704701, 1.0004513957825643, 1.0004836461342406,
    1.0005158975255328, 1.0005481499564746, 1.0005804034270993,
    1.0006126579374404, 1.0006449134875315, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0 };

  static const real_T g[16]{ 26.946201635797934, -12.927787530354905,
    -12.887952085551882, 0.0, -12.927787530354905, 6.8764758668504413,
    6.3110256368611726, 0.0, -12.887952085551882, 6.3110256368611726,
    6.8375224042336367, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const real_T c[12]{ 0.7084943675087555, 0.0, -0.34693844974269455, 0.0,
    -0.3458694000362138, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0 };

  static const int32_T b_Mrows[46]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 101, 102, 103 };

  static const int16_T d[46]{ 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    40, 40, 40, 40, 0, 0 };

  static const int8_T old_mvoff_0[3]{ 40, 0, 0 };

  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_6;
  real_T h_0[184];
  real_T f_0[138];
  real_T Bu[126];
  real_T e_0[92];
  real_T rtb_useq[63];
  real_T b_utarget[60];
  real_T b_Mlim[46];
  real_T Bv[42];
  real_T Dv[21];
  real_T vseq[21];
  real_T rseq[20];
  real_T g_0[16];
  real_T rtb_L_i[16];
  real_T b_B[12];
  real_T B_est[10];
  real_T rtb_x_0[6];
  real_T D_est[5];
  real_T Abar[4];
  real_T rtb_Q[4];
  real_T rtb_TmpSignalConversionAtSFunct[4];
  real_T rtb_x[4];
  real_T rtb_y_e[4];
  real_T rtb_y_k[4];
  real_T U[3];
  real_T rtb_Sum1[3];
  real_T tmp[3];
  real_T tmp_0[3];
  real_T tmp_1[3];
  real_T rtb_Bkuk[2];
  real_T rtb_C[2];
  real_T rtb_N[2];
  real_T rtb_Add1_tmp;
  real_T rtb_Product1;
  real_T rtb_Product2_gf;
  real_T sqrtFF;
  real_T u;
  int32_T b_utarget_tmp;
  int32_T i;
  int32_T kk;
  int32_T old_mvoff;
  int32_T rtb_L_k_tmp;
  boolean_T tmp_5[46];

  // Sum: '<S6>/Add1' incorporates:
  //   Sum: '<S7>/Sum6'

  rtb_Add1_tmp = rtu_y - rtu_y0;

  // Delay: '<S38>/delayTheta' incorporates:
  //   Constant: '<S38>/InitialParameters'

  if (localDW->icLoad_m) {
    localDW->delayTheta_DSTATE[0] = localP->InitialParameters_Value[0];
    localDW->delayTheta_DSTATE[1] = localP->InitialParameters_Value[1];
    localDW->delayTheta_DSTATE[2] = localP->InitialParameters_Value[2];
    localDW->delayTheta_DSTATE[3] = localP->InitialParameters_Value[3];
  }

  rtb_x[0] = localDW->delayTheta_DSTATE[0];
  rtb_x[1] = localDW->delayTheta_DSTATE[1];
  rtb_x[2] = localDW->delayTheta_DSTATE[2];
  rtb_x[3] = localDW->delayTheta_DSTATE[3];

  // Delay: '<S38>/delayL' incorporates:
  //   Constant: '<S38>/InitialCovariance'

  for (i = 0; i < 16; i++) {
    if (localDW->icLoad_e) {
      localDW->delayL_DSTATE[i] = localP->InitialCovariance_Value[i];
    }

    rtb_L_i[i] = localDW->delayL_DSTATE[i];
  }

  // End of Delay: '<S38>/delayL'

  // SignalConversion generated from: '<S61>/ SFunction ' incorporates:
  //   MATLAB Function: '<S38>/RLS'
  //   UnitDelay: '<S6>/Unit Delay2'
  //   UnitDelay: '<S6>/Unit Delay3'

  rtb_TmpSignalConversionAtSFunct[0] = localDW->UnitDelay3_DSTATE;
  rtb_TmpSignalConversionAtSFunct[1] = localDW->UnitDelay2_DSTATE[0];
  rtb_TmpSignalConversionAtSFunct[2] = localDW->UnitDelay2_DSTATE[1];
  rtb_TmpSignalConversionAtSFunct[3] = localDW->UnitDelay2_DSTATE[2];

  // MATLAB Function: '<S38>/RLS' incorporates:
  //   Constant: '<S38>/Forgetting Factor'
  //   DataTypeConversion: '<S38>/DataTypeConversionEnable'
  //   Delay: '<S38>/delayTheta'
  //   SignalConversion generated from: '<S61>/ SFunction '
  //   Sum: '<S6>/Add1'
  //   UnitDelay: '<S6>/Unit Delay2'
  //   UnitDelay: '<S6>/Unit Delay3'

  //  % I/O data
  //  % Helper variables
  //  % States
  //    % Was a reset trigered?
  // MATLAB Function 'Recursive Least Squares Estimator/RLS': '<S61>:1'
  //    Copyright 2013-2018 The MathWorks, Inc.
  //  Get the estimator object
  // '<S61>:1:14' if isempty(rlsEstimator)
  if (!localDW->rlsEstimator_not_empty) {
    // '<S61>:1:15' rlsEstimator = controllib.internal.blocks.rls.getEstimator(algorithmParams); 
    localDW->rlsEstimator_not_empty = true;
  }

  //  If user triggered a reset, communicate this to the estimator
  // '<S61>:1:19' if hasTriggeredReset
  //  Estimate
  // '<S61>:1:24' [yBuffer,HBuffer,x,L,e] = rlsEstimator.estimate(...
  // '<S61>:1:25'     yBuffer, HBuffer, x, L,... % states
  // '<S61>:1:26'     y, H,... % new IO data
  // '<S61>:1:27'     isEnabled, adg1, adg2, algorithmParams);
  localDW->rlsEstimator.DataIterator.IteratorPosition = 1;
  if (localDW->rlsEstimator.DataIterator.IteratorPosition + 1 <= 2) {
    localDW->rlsEstimator.DataIterator.IteratorPosition++;
  } else {
    localDW->rlsEstimator.DataIterator.IteratorPosition = 2;
  }

  rtb_Product1 = rtb_Add1_tmp - (((localDW->UnitDelay3_DSTATE *
    localDW->delayTheta_DSTATE[0] + localDW->UnitDelay2_DSTATE[0] *
    localDW->delayTheta_DSTATE[1]) + localDW->UnitDelay2_DSTATE[1] *
    localDW->delayTheta_DSTATE[2]) + localDW->UnitDelay2_DSTATE[2] *
    localDW->delayTheta_DSTATE[3]);
  if (rtu_paramEst != 0.0) {
    sqrtMeasurementUpdate(rtb_L_i, rtb_TmpSignalConversionAtSFunct,
                          localP->ForgettingFactor_Value, rtb_x);
    sqrtFF = std::sqrt(localP->ForgettingFactor_Value);
    for (kk = 0; kk < 4; kk++) {
      for (i = 0; i <= kk; i++) {
        rtb_L_k_tmp = ((((3 - kk) << 2UL) + i) - kk) + 3;
        rtb_L_i[rtb_L_k_tmp] /= sqrtFF;
      }
    }

    for (kk = 0; kk < 4; kk++) {
      if (std::isinf(rtb_L_i[kk])) {
        u = rtb_L_i[kk];
        if (std::isnan(u)) {
          u = (rtNaN);
        } else if (u < 0.0) {
          u = -1.0;
        } else {
          u = static_cast<real_T>(u > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rtb_L_i[kk] = u * 1.7976931348623157E+308;
      }

      u = rtb_L_i[kk + 4];
      if (std::isinf(u)) {
        if (std::isnan(u)) {
          u = (rtNaN);
        } else if (u < 0.0) {
          u = -1.0;
        } else {
          u = static_cast<real_T>(u > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rtb_L_i[kk + 4] = u * 1.7976931348623157E+308;
      }

      u = rtb_L_i[kk + 8];
      if (std::isinf(u)) {
        if (std::isnan(u)) {
          u = (rtNaN);
        } else if (u < 0.0) {
          u = -1.0;
        } else {
          u = static_cast<real_T>(u > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rtb_L_i[kk + 8] = u * 1.7976931348623157E+308;
      }

      u = rtb_L_i[kk + 12];
      if (std::isinf(u)) {
        if (std::isnan(u)) {
          u = (rtNaN);
        } else if (u < 0.0) {
          u = -1.0;
        } else {
          u = static_cast<real_T>(u > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rtb_L_i[kk + 12] = u * 1.7976931348623157E+308;
      }

      rtb_x[kk] = rtb_x[kk] * rtb_Product1 + localDW->delayTheta_DSTATE[kk];
    }
  }

  // Product: '<S6>/Product' incorporates:
  //   Constant: '<S6>/Constant11'

  //  % states
  //  % new IO data
  //  helper variables
  rtb_Product1 = localP->Constant11_Value * rtu_y0;

  // Delay: '<S64>/MemoryX' incorporates:
  //   Constant: '<S64>/X0'
  //   DataTypeConversion: '<S64>/DataTypeConversionReset'

  localDW->icLoad_k = ((static_cast<uint32_T>(localZCE->MemoryX_Reset_ZCE_j) ==
                        POS_ZCSIG) || localDW->icLoad_k);
  localZCE->MemoryX_Reset_ZCE_j = 0U;
  if (localDW->icLoad_k) {
    localDW->MemoryX_DSTATE[0] = localP->X0_Value[0];
    localDW->MemoryX_DSTATE[1] = localP->X0_Value[1];
  }

  // MATLAB Function: '<S36>/FixedHorizonOptimizer' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Constant: '<S6>/Constant10'
  //   Constant: '<S6>/Constant12'
  //   RandomNumber: '<S2>/Measurement Noise'
  //
  // MATLAB Function 'Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer': '<S37>:1' 
  // '<S37>:1:18' coder.extrinsic('mpcblock_optimizer_double_mex');
  // '<S37>:1:19' coder.extrinsic('mpcblock_optimizer_single_mex');
  // '<S37>:1:20' coder.extrinsic('mpcblock_refmd_double_mex');
  // '<S37>:1:21' coder.extrinsic('mpcblock_refmd_single_mex');
  //  Inputs (in BlockDataType except iA)
  //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
  // '<S37>:1:25' xk = convertDataType(xk0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    old_u:      last mv (calculated by MPC)
  // '<S37>:1:27' old_u = convertDataType(old_u0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    ym:         current measured output (used only with built-in KF)
  // '<S37>:1:29' ym = convertDataType(ym0,isDouble);
  //    ref:        output reference
  // '<S37>:1:31' ref = convertDataType(ref0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    md:         measured disturbance
  // '<S37>:1:33' md = convertDataType(md0,isDouble);
  //    umin:       run-time MV bound
  // '<S37>:1:35' umin = convertDataType(umin0,isDouble);
  //    umax:       run-time MV bound
  // '<S37>:1:37' umax = convertDataType(umax0,isDouble);
  //    ymin:       run-time OV bound
  // '<S37>:1:39' ymin = convertDataType(ymin0,isDouble);
  //    ymax:       run-time OV bound
  // '<S37>:1:41' ymax = convertDataType(ymax0,isDouble);
  //    E:          run-time mixed constraints
  // '<S37>:1:43' E = convertDataType(E0,isDouble);
  //    F:          run-time mixed constraints
  // '<S37>:1:45' F = convertDataType(F0,isDouble);
  //    G:          run-time mixed constraints
  // '<S37>:1:47' G = convertDataType(G0,isDouble);
  //    S:          run-time mixed constraints
  // '<S37>:1:49' S = convertDataType(S0,isDouble);
  //    switch_in:  if it matches "enable_value", MPC is active in control
  // '<S37>:1:51' switch_in = int32(switch_in0);
  //    ext_mv:     external last mv (actual)
  // '<S37>:1:53' ext_mv = convertDataType(ext_mv0,isDouble);
  //    MVtarget:   MV reference
  // '<S37>:1:55' MVtarget = convertDataType(MVtarget0,isDouble);
  //    ywt:        run-time OV weights
  // '<S37>:1:57' ywt = convertDataType(ywt0,isDouble);
  //    uwt:        run-time MV weights
  // '<S37>:1:59' uwt = convertDataType(uwt0,isDouble);
  //    duwt:       run-time DMV weights
  // '<S37>:1:61' duwt = convertDataType(duwt0,isDouble);
  //    rhoeps:     run-time Slack weights
  // '<S37>:1:63' ewt = convertDataType(ewt0,isDouble);
  //    a:          run-time A (must be in DT)
  // '<S37>:1:65' a = convertDataType(a0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    b:          run-time B (must be in DT)
  // '<S37>:1:67' b = convertDataType(b0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    c:          run-time C (must be in DT)
  // '<S37>:1:69' c = convertDataType(c0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    d:          run-time D (must be in DT)
  // '<S37>:1:71' d = convertDataType(d0,isDouble);
  //    U:          run-time nominal value
  // '<S37>:1:73' U = convertDataType(U0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    Y:          run-time nominal value
  // '<S37>:1:75' Y = convertDataType(Y0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    X:          run-time nominal value
  // '<S37>:1:77' X = convertDataType(X0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    DX:         run-time nominal value
  // '<S37>:1:79' DX = convertDataType(DX0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
  //    Pk:         covariance P[k|k-1] (used only with built-in KF)
  // '<S37>:1:81' Pk = convertDataType(Pk0,isDouble);
  // '<S37>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S37>:1:319' if isa(u,'double')
  // '<S37>:1:320' y = u;
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
  // '<S37>:1:115' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
  // '<S37>:1:116' isAdaptive = ~isLTV;
  //  isLTV=true forces isAdaptive=false
  // '<S37>:1:117' ZERO = zeros('like',ref);
  // '<S37>:1:118' ONE = ones('like',ref);
  // '<S37>:1:119' hasMD = nv>int32(1);
  //  Pre-allocate all the MEX block outputs for the simulation mode
  // '<S37>:1:123' if isSimulation
  //  Model update
  // '<S37>:1:137' nym = int32(length(myindex));
  // '<S37>:1:138' ai=zeros(nxp,nxp,'like',ref);
  // '<S37>:1:139' bi=zeros(nxp,nup,'like',ref);
  // '<S37>:1:140' ci=zeros(ny,nxp,'like',ref);
  // '<S37>:1:141' di=zeros(ny,nup,'like',ref);
  // '<S37>:1:143' ai(:,:)=a(:,:,1);
  // '<S37>:1:144' bi(:,:)=b(:,:,1);
  // '<S37>:1:145' ci(:,:)=c(:,:,1);
  // '<S37>:1:146' di(:,:)=d(:,:,1);
  //  Allocate matrices. Must allocate 3D matrix also in Adaptive case,
  //  otherwise EML code does not compile.
  // '<S37>:1:150' Bu=zeros(nx,nu,p+1,'like',ref);
  (void)std::memset(&Bu[0], 0, 126U * sizeof(real_T));

  // '<S37>:1:151' Bv=zeros(nx,nv,p+1,'like',ref);
  (void)std::memset(&Bv[0], 0, 42U * sizeof(real_T));

  // '<S37>:1:152' Dv=zeros(ny,nv,p+1,'like',ref);
  (void)std::memset(&Dv[0], 0, 21U * sizeof(real_T));

  // '<S37>:1:153' Dvm=zeros(nym,nv,p+1,'like',ref);
  // '<S37>:1:154' Cm=zeros(nym,nx,p+1,'like',ref);
  // '<S37>:1:155' [A(:,:,1),C(:,:,1),Bu(:,:,1),Bv(:,:,1),Cm(:,:,1),Dv(:,:,1),Dvm(:,:,1),Qk,Rk,Nk] = mpc_plantupdate(... 
  // '<S37>:1:156'     ai,bi,ci,di,A(:,:,1),B(:,:,1),C(:,:,1),D(:,:,1),mvindex,mdindex,unindex,nxp,nup,ny,nu,nv,nxid, ... 
  // '<S37>:1:157'     myindex,Uscale,Yscale,Cid,Did);
  rtb_TmpSignalConversionAtSFunct[1] = 0.0;
  rtb_TmpSignalConversionAtSFunct[2] = 0.0;
  rtb_TmpSignalConversionAtSFunct[3] = 1.0;
  (void)std::memcpy(&b_B[0], &c[0], 12U * sizeof(real_T));
  rtb_C[1] = 1.0;
  rtb_C[0] = localP->Constant12_Value;
  rtb_TmpSignalConversionAtSFunct[0] = rtb_x[0];
  i = 0;
  for (rtb_L_k_tmp = 0; rtb_L_k_tmp < 3; rtb_L_k_tmp++) {
    b_B[i] = rtb_x[rtb_L_k_tmp + 1];
    Bu[i] = b_B[i];
    Bu[i + 1] = b_B[i + 1];
    i += 2;
  }

  Bv[1] = b_B[7];
  Dv[0] = 0.0;

  // '<S37>:1:158' if isLTV
  //  Offset update together with Mlim, utarget, Bv and Dv values
  // '<S37>:1:174' [Mlim, utarget, uoff, voff, yoff, myoff, xoff, Bv, Dv] = ...
  // '<S37>:1:175'     mpc_updateFromNominal(isAdaptive,isQP,Mlim,Mrows,...
  // '<S37>:1:176'        U,Uscale,uoff,mvindex,voff,mdindex,utarget,nu,nv-1,... 
  // '<S37>:1:177'        Y,Yscale,yoff,myoff,myindex,ny,...
  // '<S37>:1:178'        X,xoff,nxp,DX,A,Bu,Bv,C,Dv,nCC);
  for (i = 0; i < 46; i++) {
    b_Mlim[i] = static_cast<real_T>(d[i]);
  }

  (void)std::memset(&b_utarget[0], 0, 60U * sizeof(real_T));
  U[0] = rtu_u0[0];
  rtb_Sum1[0] = rtu_u0[0];
  U[1] = rtu_u0[1];
  rtb_Sum1[1] = rtu_u0[1];
  U[2] = rtu_u0[2];
  rtb_Sum1[2] = rtu_u0[2];
  for (kk = 0; kk < 46; kk++) {
    u = b_Mlim[kk];
    i = b_Mrows[kk];
    if (i <= 20) {
      u += yoff - rtu_y0;
    } else if (i <= 40) {
      u -= yoff - rtu_y0;
    } else if (i <= 100) {
      i = (i - div_nde_s32_floor(i - 41, static_cast<int32_T>(nu)) *
           static_cast<int32_T>(nu)) - 41;
      u += static_cast<real_T>(old_mvoff_0[i]) - U[i];
    } else {
      i = (i - div_nde_s32_floor(i - 101, static_cast<int32_T>(nu)) *
           static_cast<int32_T>(nu)) - 101;
      u -= static_cast<real_T>(old_mvoff_0[i]) - U[i];
    }

    b_Mlim[kk] = u;
  }

  for (kk = 0; kk < 3; kk++) {
    old_mvoff = static_cast<int32_T>(old_mvoff_0[kk]);
    u = U[kk];
    i = 0;
    for (rtb_L_k_tmp = 0; rtb_L_k_tmp < 20; rtb_L_k_tmp++) {
      b_utarget_tmp = i + kk;
      b_utarget[b_utarget_tmp] = (b_utarget[b_utarget_tmp] + static_cast<real_T>
        (old_mvoff)) - u;
      i += 3;
    }
  }

  Bv[0] = localP->Constant10_Value;

  //  Remove last u offset
  // '<S37>:1:181' old_u = old_u - uoff;
  //  Get reference and MD signals -- accounting for previewing
  // '<S37>:1:184' if isSimulation
  // '<S37>:1:190' else
  //  When doing code generation, use M code directly
  // '<S37>:1:192' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
  for (kk = 0; kk < 21; kk++) {
    vseq[kk] = 1.0;
  }

  for (kk = 0; kk < 20; kk++) {
    rseq[kk] = rtu_r - rtu_y0;
  }

  // Sum: '<S7>/Sum2' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Constant: '<S7>/Constant'
  //   Delay: '<S64>/MemoryX'
  //   MATLAB Function: '<S36>/FixedHorizonOptimizer'

  //  External MV override.
  //  NOTE: old_u and ext_mv input signals are dimensionless and offset-free.
  // '<S37>:1:197' if no_mv
  //  no external mv: old_u is the optimal u[k-1] from last_mv
  // '<S37>:1:199' delmv = zeros(nu,1,'like',ref);
  //  Obtain x[k|k]
  // '<S37>:1:208' xk = xk - xoff;
  //  Remove offset
  // '<S37>:1:209' if CustomEstimation
  //  Input is x(k|k)
  // '<S37>:1:211' xest = xk;
  //  Real-time MV target override
  //  Note: utargetValue is a vector length p*nu.
  // '<S37>:1:231' if no_uref
  //  no external utarget
  // '<S37>:1:233' utargetValue = utarget;
  //  Real-time custom constraint override (scaled E/F/S)
  // '<S37>:1:241' if ~no_cc
  // '<S37>:1:250' return_sequence = return_mvseq || return_xseq || return_ovseq; 
  // '<S37>:1:251' if isSimulation
  // '<S37>:1:279' else
  //  When doing code generation, use M code directly
  // '<S37>:1:281' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
  // '<S37>:1:282'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
  // '<S37>:1:283'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
  // '<S37>:1:284'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
  // '<S37>:1:285'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
  // '<S37>:1:286'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
  // '<S37>:1:287'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
  // '<S37>:1:288'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
  // '<S37>:1:289'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
  // '<S37>:1:290'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
  // '<S37>:1:291'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
  // '<S37>:1:292'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, isHyb, Mdis, Ndis, Vdis, numdis, maxdis); 
  rtb_N[0] = (localDW->MemoryX_DSTATE[0] + rtb_Product1) - rtb_Product1;
  rtb_N[1] = localDW->MemoryX_DSTATE[1] + localP->Constant_Value;

  // MATLAB Function: '<S36>/FixedHorizonOptimizer' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   UnitDelay: '<S8>/last_mv'

  tmp[0] = localDW->last_mv_DSTATE[0] - rtu_u0[0];
  tmp[1] = localDW->last_mv_DSTATE[1] - rtu_u0[1];
  tmp[2] = localDW->last_mv_DSTATE[2] - rtu_u0[2];
  (void)std::memset(&localDW->dv[0], 0, 966U * sizeof(real_T));
  tmp_0[0] = 0.54594344475769718;
  tmp_1[0] = 0.0;
  tmp_0[1] = 0.54594344475769718;
  tmp_1[1] = 0.0;
  tmp_0[2] = 0.54594344475769718;
  tmp_1[2] = 0.0;

  // Memory: '<S8>/Memory'
  for (i = 0; i < 46; i++) {
    tmp_5[i] = localDW->Memory_PreviousInput[i];
  }

  // End of Memory: '<S8>/Memory'

  // MATLAB Function: '<S36>/FixedHorizonOptimizer'
  (void)std::memcpy(&e_0[0], &e[0], 92U * sizeof(real_T));
  (void)std::memcpy(&f_0[0], &f[0], 138U * sizeof(real_T));
  (void)std::memcpy(&g_0[0], &g[0], sizeof(real_T) << 4UL);
  (void)std::memcpy(&h_0[0], &h[0], 184U * sizeof(real_T));

  // Update for Memory: '<S8>/Memory' incorporates:
  //   MATLAB Function: '<S36>/FixedHorizonOptimizer'

  mpcblock_optimizer(rseq, vseq, rtb_N, tmp, tmp_5, b_Mlim, e_0, f_0,
                     localDW->dv, b_utarget, rtb_Sum1, g_0, h_0, tmp_0, j, tmp_1,
                     k, rtb_TmpSignalConversionAtSFunct, Bu, Bv, rtb_C, Dv,
                     b_Mrows, U, rtb_useq, &sqrtFF,
                     localDW->Memory_PreviousInput, localDW);

  // MATLAB Function: '<S7>/MATLAB Function' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Constant: '<S6>/Constant12'
  //   Constant: '<S6>/Constant13'
  //   MATLAB Function: '<S86>/ScalarExpansion'

  // '<S37>:1:295' if return_xseq || return_ovseq
  // '<S37>:1:297' else
  // '<S37>:1:298' yseq = zeros(p+1,ny,'like',rseq);
  // '<S37>:1:299' xseq = zeros(p+1,nxQP,'like',rseq);
  // '<S37>:1:302' if CustomEstimation
  // '<S37>:1:303' xk1 = zeros(nx,1,'like',ref);
  // '<S37>:1:304' Pk1 = Pk;
  // '<S37>:1:311' xk1 = xk1 + xoff;
  //  Updated state must include offset
  //  return xest in original value
  // '<S37>:1:314' xest = xest + xoff;
  // MATLAB Function 'SupervisoryController/State0.controlLaw.AMPC0/State Estimator (KF)/MATLAB Function': '<S65>:1' 
  // '<S65>:1:2' [A, B, C, D, Q, R, N] = stateEstimator(Ap, Bp, Cp, Dp, Aod0, Bod0, Cod0, Dod0, Dmn0); 
  // 'stateEstimator:3' no = size(Cp,1);
  //  n_outputs
  // 'stateEstimator:4' ni = 3;
  //  n_inputs
  // 'stateEstimator:5' nsp = size(Ap,1);
  //  n_plant_states
  // 'stateEstimator:6' ns = nsp + no;
  //  n_states = n_plant_states + n_outputs
  // 'stateEstimator:8' A = zeros(ns);
  //  n_states x n_states
  // 'stateEstimator:9' B = zeros(ns,ni);
  //  n_states  x n_inputs
  // 'stateEstimator:10' C = zeros(no,ns);
  //  n_outputs x n_states
  // 'stateEstimator:11' D = zeros(no,ni);
  //  n_outputs x n_inputs
  // 'stateEstimator:12' Q = zeros(ns,ns);
  //  n_states  x n_states
  // 'stateEstimator:13' G = eye(ns);
  //  n_states  x n_states
  // 'stateEstimator:14' R = zeros(no,no);
  //  n_outputs x n_outputs
  // 'stateEstimator:15' N = zeros(ns,no);
  //  n_states  x n_outputs
  // 'stateEstimator:16' H = zeros(no,ns);
  //  n_outputs x n_states
  //  combine plant and output disturbance model
  // 'stateEstimator:19' A = [Ap, zeros(nsp,no);
  // 'stateEstimator:20'      zeros(no,nsp), Aod];
  rtb_TmpSignalConversionAtSFunct[0] = rtb_x[0];
  rtb_TmpSignalConversionAtSFunct[2] = 0.0;
  rtb_TmpSignalConversionAtSFunct[1] = 0.0;
  rtb_TmpSignalConversionAtSFunct[3] = rtP->Aod0;

  // 'stateEstimator:21' B = [Bp;
  // 'stateEstimator:22'      zeros(no,ni)];
  // 'stateEstimator:23' C = [Cp Cod];
  rtb_C[0] = localP->Constant12_Value;
  rtb_C[1] = rtP->Cod0;

  // 'stateEstimator:24' D = Dp;
  // 'stateEstimator:26' B_est = [ [Bp; zeros(no,ni)] [zeros(size(Bp,1), size(Cp,1)); Bod] [zeros(size(Bp,1) + size(Bod,1), size(Cp,1))] ]; 
  B_est[0] = rtb_x[1];
  B_est[1] = 0.0;
  B_est[2] = rtb_x[2];
  B_est[3] = 0.0;
  B_est[4] = rtb_x[3];
  B_est[5] = 0.0;
  B_est[6] = 0.0;
  B_est[7] = rtP->Bod0;
  B_est[8] = 0.0;
  B_est[9] = 0.0;

  // 'stateEstimator:27' D_est = [Dp Dod Dn];
  D_est[0] = localP->Constant13_Value[0];
  D_est[1] = localP->Constant13_Value[1];
  D_est[2] = localP->Constant13_Value[2];
  D_est[3] = rtP->Dod0;
  D_est[4] = rtP->Dmn0;

  // 'stateEstimator:28' Q = B_est * B_est';
  for (i = 0; i < 2; i++) {
    rtb_L_k_tmp = 0;
    for (kk = 0; kk < 2; kk++) {
      int32_T rtb_Q_tmp;
      rtb_Q_tmp = rtb_L_k_tmp + i;
      rtb_Q[rtb_Q_tmp] = 0.0;
      old_mvoff = 0;
      for (b_utarget_tmp = 0; b_utarget_tmp < 5; b_utarget_tmp++) {
        rtb_Q[rtb_Q_tmp] += B_est[old_mvoff + i] * B_est[old_mvoff + kk];
        old_mvoff += 2;
      }

      rtb_L_k_tmp += 2;
    }
  }

  // 'stateEstimator:29' R = D_est * D_est';
  rtb_Product1 = 0.0;
  for (i = 0; i < 5; i++) {
    u = D_est[i];
    rtb_Product1 += u * u;
  }

  // Outputs for Atomic SubSystem: '<S64>/ScalarExpansionQ'
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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S108>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S108>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  i = 0;
  for (rtb_L_k_tmp = 0; rtb_L_k_tmp < 2; rtb_L_k_tmp++) {
    rtb_Bkuk[rtb_L_k_tmp] = 0.0;
    kk = 0;
    for (old_mvoff = 0; old_mvoff < 5; old_mvoff++) {
      rtb_Bkuk[rtb_L_k_tmp] += B_est[kk + rtb_L_k_tmp] * D_est[old_mvoff];
      kk += 2;
    }

    rtb_N[rtb_L_k_tmp] = rtb_Bkuk[rtb_L_k_tmp];
    rtb_y_k[i] = (rtb_Q[i] + rtb_Q[rtb_L_k_tmp]) / 2.0;
    rtb_y_k[i + 1] = (rtb_Q[i + 1] + rtb_Q[rtb_L_k_tmp + 2]) / 2.0;
    i += 2;
  }

  // End of Outputs for SubSystem: '<S64>/ScalarExpansionQ'

  // Outputs for Atomic SubSystem: '<S64>/ScalarExpansionR'
  // MATLAB Function: '<S87>/ScalarExpansion' incorporates:
  //   MATLAB Function: '<S7>/MATLAB Function'

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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S109>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S109>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  if (rtb_Product1 <= 0.0) {
    rtb_Product1 = 2.2204460492503131E-16;
  }

  // End of Outputs for SubSystem: '<S64>/ScalarExpansionR'
  for (i = 0; i <= 0; i += 2) {
    // Outputs for Atomic SubSystem: '<S64>/ReducedQRN'
    // Product: '<S84>/Product' incorporates:
    //   Constant: '<S64>/G'
    //   Math: '<S84>/Transpose1'

    tmp_4 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_y_e[i], tmp_4);
    tmp_2 = _mm_loadu_pd(&rtb_y_k[i]);
    tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
    (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
      (localP->G_Value[0])), tmp_3));
    tmp_2 = _mm_loadu_pd(&rtb_y_k[i + 2]);
    tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
    (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
      (localP->G_Value[2])), tmp_3));
    (void)_mm_storeu_pd(&rtb_y_e[i + 2], tmp_4);
    tmp_4 = _mm_loadu_pd(&rtb_y_k[i]);
    tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
    (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(tmp_2, _mm_mul_pd(tmp_4,
      _mm_set1_pd(localP->G_Value[1]))));
    tmp_4 = _mm_loadu_pd(&rtb_y_k[i + 2]);
    tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
    (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(_mm_mul_pd(tmp_4,
      _mm_set1_pd(localP->G_Value[3])), tmp_2));

    // End of Outputs for SubSystem: '<S64>/ReducedQRN'
  }

  for (i = 0; i <= 0; i += 2) {
    // Outputs for Atomic SubSystem: '<S64>/ReducedQRN'
    // Product: '<S84>/Product' incorporates:
    //   Constant: '<S64>/G'

    tmp_4 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_Q[i], tmp_4);
    tmp_2 = _mm_loadu_pd(&rtb_Q[i]);

    // Constant: '<S64>/G' incorporates:
    //   Product: '<S84>/Product'

    tmp_3 = _mm_loadu_pd(&localP->G_Value[i]);

    // Product: '<S84>/Product' incorporates:
    //   Constant: '<S64>/G'

    (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_y_e[0])), tmp_2));
    tmp_2 = _mm_loadu_pd(&rtb_Q[i]);

    // Constant: '<S64>/G' incorporates:
    //   Product: '<S84>/Product'

    tmp_6 = _mm_loadu_pd(&localP->G_Value[i + 2]);

    // Product: '<S84>/Product' incorporates:
    //   Constant: '<S64>/G'

    (void)_mm_storeu_pd(&rtb_Q[i], _mm_add_pd(_mm_mul_pd(tmp_6, _mm_set1_pd
      (rtb_y_e[1])), tmp_2));
    (void)_mm_storeu_pd(&rtb_Q[i + 2], tmp_4);
    tmp_4 = _mm_loadu_pd(&rtb_Q[i + 2]);
    (void)_mm_storeu_pd(&rtb_Q[i + 2], _mm_add_pd(tmp_4, _mm_mul_pd(tmp_3,
      _mm_set1_pd(rtb_y_e[2]))));
    tmp_4 = _mm_loadu_pd(&rtb_Q[i + 2]);
    (void)_mm_storeu_pd(&rtb_Q[i + 2], _mm_add_pd(_mm_mul_pd(tmp_6, _mm_set1_pd
      (rtb_y_e[3])), tmp_4));

    // Sum: '<S84>/Add' incorporates:
    //   Constant: '<S64>/G'
    //   Constant: '<S64>/H'
    //   Math: '<S84>/Transpose2'
    //   Product: '<S83>/B[k]*u[k]'
    //   Product: '<S84>/Product'
    //   Product: '<S84>/Product1'

    tmp_4 = _mm_loadu_pd(&rtb_y_k[i + 2]);
    tmp_2 = _mm_loadu_pd(&rtb_y_k[i]);
    tmp_3 = _mm_loadu_pd(&rtb_N[i]);
    (void)_mm_storeu_pd(&rtb_Bkuk[i], _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_4,
      _mm_set1_pd(localP->H_Value[1])), _mm_mul_pd(tmp_2, _mm_set1_pd
      (localP->H_Value[0]))), tmp_3));

    // End of Outputs for SubSystem: '<S64>/ReducedQRN'
  }

  // Outputs for Atomic SubSystem: '<S64>/ReducedQRN'
  // Product: '<S84>/Product2' incorporates:
  //   Constant: '<S64>/G'
  //   Product: '<S83>/B[k]*u[k]'

  u = localP->G_Value[0] * rtb_Bkuk[0] + rtb_Bkuk[1] * localP->G_Value[2];
  rtb_Product2_gf = rtb_Bkuk[0] * localP->G_Value[1] + rtb_Bkuk[1] *
    localP->G_Value[3];

  // Outputs for Atomic SubSystem: '<S64>/ScalarExpansionR'
  // Sum: '<S84>/Add1' incorporates:
  //   Constant: '<S64>/H'
  //   MATLAB Function: '<S87>/ScalarExpansion'
  //   Math: '<S84>/Transpose'
  //   Math: '<S84>/Transpose2'
  //   Product: '<S83>/B[k]*u[k]'
  //   Product: '<S84>/Product3'
  //   Product: '<S84>/Product4'

  rtb_Product1 += (localP->H_Value[0] * rtb_Bkuk[0] + localP->H_Value[1] *
                   rtb_Bkuk[1]) + (rtb_N[0] * localP->H_Value[0] + rtb_N[1] *
    localP->H_Value[1]);

  // End of Outputs for SubSystem: '<S64>/ScalarExpansionR'
  // End of Outputs for SubSystem: '<S64>/ReducedQRN'

  // Delay: '<S64>/MemoryP' incorporates:
  //   Constant: '<S64>/P0'
  //   DataTypeConversion: '<S64>/DataTypeConversionReset'

  localDW->icLoad_f = ((static_cast<uint32_T>(localZCE->MemoryP_Reset_ZCE_a) ==
                        POS_ZCSIG) || localDW->icLoad_f);
  localZCE->MemoryP_Reset_ZCE_a = 0U;
  if (localDW->icLoad_f) {
    localDW->MemoryP_DSTATE[0] = localP->P0_Value[0];
    localDW->MemoryP_DSTATE[1] = localP->P0_Value[1];
    localDW->MemoryP_DSTATE[2] = localP->P0_Value[2];
    localDW->MemoryP_DSTATE[3] = localP->P0_Value[3];
  }

  // Outputs for Atomic SubSystem: '<S64>/CalculatePL'
  // MATLAB Function: '<S66>/Discrete-Time KF - Calculate PLMZ' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Constant: '<S2>/Constant'
  //   Constant: '<S6>/Constant12'
  //   DataTypeConversion: '<S64>/DataTypeConversionEnable'
  //   Delay: '<S64>/MemoryP'
  //   MATLAB Function: '<S7>/MATLAB Function'
  //   Product: '<S84>/Product'

  //  See help of ctrlKalmanFilterDTCalculatePL.m
  // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S104>:1' 
  //    Copyright 2014 The MathWorks, Inc.
  // '<S104>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
  if (localP->Constant_Value_b != 0.0) {
    real_T rtb_Z_idx_1;
    real_T rtb_Z_idx_2;
    real_T rtb_Z_idx_3;
    sqrtFF = 0.0;
    i = 0;
    for (rtb_L_k_tmp = 0; rtb_L_k_tmp < 2; rtb_L_k_tmp++) {
      rtb_N[rtb_L_k_tmp] = 0.0;
      rtb_N[rtb_L_k_tmp] += localDW->MemoryP_DSTATE[i] *
        localP->Constant12_Value;
      rtb_y_e[rtb_L_k_tmp] = 0.0;
      rtb_y_e[rtb_L_k_tmp] += rtb_TmpSignalConversionAtSFunct[rtb_L_k_tmp] *
        localDW->MemoryP_DSTATE[0];
      rtb_Z_idx_1 = rtb_TmpSignalConversionAtSFunct[rtb_L_k_tmp + 2];
      rtb_y_e[rtb_L_k_tmp] += rtb_Z_idx_1 * localDW->MemoryP_DSTATE[1];
      rtb_N[rtb_L_k_tmp] += localDW->MemoryP_DSTATE[i + 1] * rtP->Cod0;
      rtb_y_e[rtb_L_k_tmp + 2] = 0.0;
      rtb_y_e[rtb_L_k_tmp + 2] += rtb_TmpSignalConversionAtSFunct[rtb_L_k_tmp] *
        localDW->MemoryP_DSTATE[2];
      rtb_y_e[rtb_L_k_tmp + 2] += rtb_Z_idx_1 * localDW->MemoryP_DSTATE[3];
      sqrtFF += rtb_N[rtb_L_k_tmp] * rtb_C[rtb_L_k_tmp];
      i += 2;
    }

    sqrtFF += rtb_Product1;
    rtb_N[0] = ((rtb_y_e[0] * localP->Constant12_Value + rtb_y_e[2] * rtP->Cod0)
                + u) / sqrtFF;
    rtb_Bkuk[0] = (localDW->MemoryP_DSTATE[0] * localP->Constant12_Value +
                   localDW->MemoryP_DSTATE[2] * rtP->Cod0) / sqrtFF;
    rtb_N[1] = ((rtb_y_e[1] * localP->Constant12_Value + rtb_y_e[3] * rtP->Cod0)
                + rtb_Product2_gf) / sqrtFF;
    rtb_Bkuk[1] = (localDW->MemoryP_DSTATE[1] * localP->Constant12_Value +
                   localDW->MemoryP_DSTATE[3] * rtP->Cod0) / sqrtFF;
    rtb_y_k[0] = 1.0 - rtb_Bkuk[0] * localP->Constant12_Value;
    rtb_y_k[1] = 0.0 - localP->Constant12_Value * rtb_Bkuk[1];
    rtb_y_k[2] = 0.0 - rtb_Bkuk[0] * rtP->Cod0;
    rtb_y_k[3] = 1.0 - rtb_Bkuk[1] * rtP->Cod0;
    for (i = 0; i <= 0; i += 2) {
      tmp_4 = _mm_set1_pd(0.0);
      (void)_mm_storeu_pd(&Abar[i], tmp_4);
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i]);
      tmp_3 = _mm_loadu_pd(&Abar[i]);
      (void)_mm_storeu_pd(&Abar[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (localDW->MemoryP_DSTATE[0])), tmp_3));
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i + 2]);
      tmp_3 = _mm_loadu_pd(&Abar[i]);
      (void)_mm_storeu_pd(&Abar[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (localDW->MemoryP_DSTATE[1])), tmp_3));
      (void)_mm_storeu_pd(&Abar[i + 2], tmp_4);
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i]);
      tmp_3 = _mm_loadu_pd(&Abar[i + 2]);
      (void)_mm_storeu_pd(&Abar[i + 2], _mm_add_pd(tmp_3, _mm_mul_pd(tmp_2,
        _mm_set1_pd(localDW->MemoryP_DSTATE[2]))));
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i + 2]);
      tmp_3 = _mm_loadu_pd(&Abar[i + 2]);
      (void)_mm_storeu_pd(&Abar[i + 2], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (localDW->MemoryP_DSTATE[3])), tmp_3));
      (void)_mm_storeu_pd(&rtb_y_e[i], tmp_4);
      tmp_2 = _mm_loadu_pd(&Abar[i]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (rtb_y_k[0])), tmp_3));
      tmp_2 = _mm_loadu_pd(&Abar[i + 2]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (rtb_y_k[2])), tmp_3));
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], tmp_4);
      tmp_4 = _mm_loadu_pd(&Abar[i]);
      tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(tmp_2, _mm_mul_pd(tmp_4,
        _mm_set1_pd(rtb_y_k[1]))));
      tmp_4 = _mm_loadu_pd(&Abar[i + 2]);
      tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(_mm_mul_pd(tmp_4,
        _mm_set1_pd(rtb_y_k[3])), tmp_2));
    }

    rtb_Z_idx_2 = rtb_Bkuk[0] * rtb_Product1;
    sqrtFF = rtb_Z_idx_2 * rtb_Bkuk[0] + rtb_y_e[0];
    rtb_Z_idx_3 = rtb_Bkuk[1] * rtb_Product1;
    rtb_Z_idx_1 = rtb_Z_idx_3 * rtb_Bkuk[0] + rtb_y_e[1];
    rtb_Z_idx_2 = rtb_Z_idx_2 * rtb_Bkuk[1] + rtb_y_e[2];
    rtb_Z_idx_3 = rtb_Z_idx_3 * rtb_Bkuk[1] + rtb_y_e[3];
    rtb_Bkuk[0] = u / rtb_Product1;
    rtb_Bkuk[1] = rtb_Product2_gf / rtb_Product1;
    rtb_y_k[0] = rtb_x[0] - rtb_Bkuk[0] * localP->Constant12_Value;
    rtb_y_k[1] = 0.0 - localP->Constant12_Value * rtb_Bkuk[1];
    rtb_y_k[2] = 0.0 - rtb_Bkuk[0] * rtP->Cod0;
    rtb_y_k[3] = rtP->Aod0 - rtb_Bkuk[1] * rtP->Cod0;
    for (i = 0; i <= 0; i += 2) {
      tmp_4 = _mm_set1_pd(0.0);
      (void)_mm_storeu_pd(&Abar[i], tmp_4);
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i]);
      tmp_3 = _mm_loadu_pd(&Abar[i]);
      (void)_mm_storeu_pd(&Abar[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (sqrtFF)), tmp_3));
      tmp_2 = _mm_loadu_pd(&rtb_y_k[i + 2]);
      tmp_3 = _mm_loadu_pd(&Abar[i]);
      (void)_mm_storeu_pd(&Abar[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (rtb_Z_idx_1)), tmp_3));
      (void)_mm_storeu_pd(&Abar[i + 2], tmp_4);
      tmp_4 = _mm_loadu_pd(&rtb_y_k[i]);
      tmp_2 = _mm_loadu_pd(&Abar[i + 2]);
      (void)_mm_storeu_pd(&Abar[i + 2], _mm_add_pd(tmp_2, _mm_mul_pd(tmp_4,
        _mm_set1_pd(rtb_Z_idx_2))));
      tmp_4 = _mm_loadu_pd(&rtb_y_k[i + 2]);
      tmp_2 = _mm_loadu_pd(&Abar[i + 2]);
      (void)_mm_storeu_pd(&Abar[i + 2], _mm_add_pd(_mm_mul_pd(tmp_4, _mm_set1_pd
        (rtb_Z_idx_3)), tmp_2));
      tmp_4 = _mm_loadu_pd(&Abar[i + 2]);
      tmp_2 = _mm_loadu_pd(&Abar[i]);
      tmp_3 = _mm_loadu_pd(&rtb_Q[i]);
      (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_4,
        _mm_set1_pd(rtb_y_k[2])), _mm_mul_pd(tmp_2, _mm_set1_pd(rtb_y_k[0]))),
        tmp_3));
      tmp_4 = _mm_loadu_pd(&Abar[i + 2]);
      tmp_2 = _mm_loadu_pd(&Abar[i]);
      tmp_3 = _mm_loadu_pd(&rtb_Q[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(_mm_add_pd(_mm_mul_pd
        (tmp_4, _mm_set1_pd(rtb_y_k[3])), _mm_mul_pd(tmp_2, _mm_set1_pd(rtb_y_k
        [1]))), tmp_3));
    }

    rtb_y_k[0] = rtb_y_e[0] - rtb_Bkuk[0] * u;
    rtb_y_k[1] = rtb_y_e[1] - rtb_Bkuk[1] * u;
    rtb_y_k[2] = rtb_y_e[2] - rtb_Bkuk[0] * rtb_Product2_gf;
    rtb_y_k[3] = rtb_y_e[3] - rtb_Bkuk[1] * rtb_Product2_gf;
  } else {
    for (i = 0; i <= 0; i += 2) {
      tmp_4 = _mm_set1_pd(0.0);
      (void)_mm_storeu_pd(&rtb_N[i], tmp_4);
      (void)_mm_storeu_pd(&rtb_y_e[i], tmp_4);
      tmp_2 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (localDW->MemoryP_DSTATE[0])), tmp_3));
      tmp_2 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i + 2]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      (void)_mm_storeu_pd(&rtb_y_e[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (localDW->MemoryP_DSTATE[1])), tmp_3));
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], tmp_4);
      tmp_2 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(tmp_3, _mm_mul_pd(tmp_2,
        _mm_set1_pd(localDW->MemoryP_DSTATE[2]))));
      tmp_2 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i + 2]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_e[i + 2], _mm_add_pd(_mm_mul_pd(tmp_2,
        _mm_set1_pd(localDW->MemoryP_DSTATE[3])), tmp_3));
      tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      tmp_6 = _mm_loadu_pd(&rtb_Q[i]);
      (void)_mm_storeu_pd(&rtb_y_k[i], _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_2,
        tmp_4), _mm_mul_pd(tmp_3, _mm_set1_pd(rtb_TmpSignalConversionAtSFunct[0]))),
        tmp_6));
      tmp_2 = _mm_loadu_pd(&rtb_y_e[i + 2]);
      tmp_3 = _mm_loadu_pd(&rtb_y_e[i]);
      tmp_6 = _mm_loadu_pd(&rtb_Q[i + 2]);
      (void)_mm_storeu_pd(&rtb_y_k[i + 2], _mm_add_pd(_mm_add_pd(_mm_mul_pd
        (tmp_2, _mm_set1_pd(rtP->Aod0)), _mm_mul_pd(tmp_3, tmp_4)), tmp_6));
    }
  }

  // End of MATLAB Function: '<S66>/Discrete-Time KF - Calculate PLMZ'
  // End of Outputs for SubSystem: '<S64>/CalculatePL'

  // RandomNumber: '<S2>/Measurement Noise'
  //  Determine if the Square-Root algorithm was used
  // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S106>:1' 
  // '<S106>:1:4' if isSqrtUsed
  rty_u[0] = localDW->NextOutput[0];

  // Sum: '<S2>/Sum' incorporates:
  //   Gain: '<S8>/u_scale'
  //   Product: '<S2>/Product'

  rty_u[0] = localP->u_scale_Gain[0] * U[0] + rtu_paramEst * rty_u[0];

  // Sum: '<S7>/Sum1' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'

  rtb_Sum1[0] = rty_u[0] - rtu_u0[0];

  // RandomNumber: '<S2>/Measurement Noise'
  rty_u[1] = localDW->NextOutput[1];

  // Sum: '<S2>/Sum' incorporates:
  //   Gain: '<S8>/u_scale'
  //   Product: '<S2>/Product'

  rty_u[1] = localP->u_scale_Gain[1] * U[1] + rtu_paramEst * rty_u[1];

  // Sum: '<S7>/Sum1' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'

  rtb_Sum1[1] = rty_u[1] - rtu_u0[1];

  // RandomNumber: '<S2>/Measurement Noise'
  rty_u[2] = localDW->NextOutput[2];

  // Sum: '<S2>/Sum' incorporates:
  //   Gain: '<S8>/u_scale'
  //   Product: '<S2>/Product'

  rty_u[2] = localP->u_scale_Gain[2] * U[2] + rtu_paramEst * rty_u[2];

  // Sum: '<S7>/Sum1' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Sum: '<S6>/Add3'

  u = rty_u[2] - rtu_u0[2];

  // Outputs for Enabled SubSystem: '<S83>/MeasurementUpdate' incorporates:
  //   EnablePort: '<S107>/Enable'

  // DataTypeConversion: '<S64>/DataTypeConversionEnable' incorporates:
  //   Constant: '<S2>/Constant'

  if (localP->Constant_Value_b != 0.0) {
    localDW->MeasurementUpdate_MODE = true;

    // Sum: '<S107>/Sum' incorporates:
    //   BusCreator: '<S6>/Bus Creator1'
    //   Constant: '<S6>/Constant13'
    //   Delay: '<S64>/MemoryX'
    //   Product: '<S107>/C[k]*xhat[k|k-1]'
    //   Product: '<S107>/D[k]*u[k]'
    //   Sum: '<S107>/Add1'
    //   Sum: '<S7>/Sum1'

    rtb_Product1 = rtb_Add1_tmp - (((localP->Constant13_Value[0] * rtb_Sum1[0] +
      localP->Constant13_Value[1] * rtb_Sum1[1]) + localP->Constant13_Value[2] *
      u) + (rtb_C[0] * localDW->MemoryX_DSTATE[0] + rtb_C[1] *
            localDW->MemoryX_DSTATE[1]));

    // Product: '<S107>/Product3'
    localDW->Product3[0] = rtb_N[0] * rtb_Product1;
    localDW->Product3[1] = rtb_N[1] * rtb_Product1;
  } else if (localDW->MeasurementUpdate_MODE) {
    // Disable for Product: '<S107>/Product3' incorporates:
    //   Outport: '<S107>/L*(y[k]-yhat[k|k-1])'

    localDW->Product3[0] = localP->Lykyhatkk1_Y0;
    localDW->Product3[1] = localP->Lykyhatkk1_Y0;
    localDW->MeasurementUpdate_MODE = false;
  } else {
    // no actions
  }

  // End of Outputs for SubSystem: '<S83>/MeasurementUpdate'

  // Update for UnitDelay: '<S6>/Unit Delay3' incorporates:
  //   Sum: '<S6>/Add1'

  localDW->UnitDelay3_DSTATE = rtb_Add1_tmp;

  // Update for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[0] = U[0];

  // Update for UnitDelay: '<S6>/Unit Delay2' incorporates:
  //   Sum: '<S6>/Add3'

  localDW->UnitDelay2_DSTATE[0] = rty_u[0] - rtu_u0[0];

  // Update for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[1] = U[1];

  // Update for UnitDelay: '<S6>/Unit Delay2' incorporates:
  //   Sum: '<S6>/Add3'

  localDW->UnitDelay2_DSTATE[1] = rty_u[1] - rtu_u0[1];

  // Update for UnitDelay: '<S8>/last_mv'
  localDW->last_mv_DSTATE[2] = U[2];

  // Update for UnitDelay: '<S6>/Unit Delay2'
  localDW->UnitDelay2_DSTATE[2] = u;

  // Sum: '<S7>/Sum3' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'
  //   Constant: '<S6>/Constant13'
  //   Delay: '<S64>/MemoryX'
  //   Product: '<S67>/Product'
  //   Product: '<S67>/Product1'
  //   Sum: '<S67>/Add1'
  //   Sum: '<S7>/Sum1'

  *rty_yhat = (((localP->Constant13_Value[0] * rtb_Sum1[0] +
                 localP->Constant13_Value[1] * rtb_Sum1[1]) +
                localP->Constant13_Value[2] * u) + (rtb_C[0] *
    localDW->MemoryX_DSTATE[0] + rtb_C[1] * localDW->MemoryX_DSTATE[1])) +
    rtu_y0;

  // Update for Delay: '<S38>/delayBuffery'
  localDW->icLoad = false;

  // Update for Delay: '<S38>/delayBufferH'
  localDW->icLoad_n = false;

  // Update for Delay: '<S38>/delayTheta'
  localDW->icLoad_m = false;
  localDW->delayTheta_DSTATE[0] = rtb_x[0];
  localDW->delayTheta_DSTATE[1] = rtb_x[1];
  localDW->delayTheta_DSTATE[2] = rtb_x[2];
  localDW->delayTheta_DSTATE[3] = rtb_x[3];

  // Update for Delay: '<S38>/delayL'
  localDW->icLoad_e = false;
  (void)std::memcpy(&localDW->delayL_DSTATE[0], &rtb_L_i[0], sizeof(real_T) <<
                    4UL);

  // Update for Delay: '<S64>/MemoryX'
  localDW->icLoad_k = false;

  // MATLAB Function: '<S7>/MATLAB Function' incorporates:
  //   BusCreator: '<S6>/Bus Creator1'

  rtb_x_0[0] = rtb_x[1];
  rtb_x_0[1] = 0.0;
  rtb_x_0[2] = rtb_x[2];
  rtb_x_0[3] = 0.0;
  rtb_x_0[4] = rtb_x[3];
  rtb_x_0[5] = 0.0;
  for (i = 0; i <= 0; i += 2) {
    // Product: '<S83>/B[k]*u[k]' incorporates:
    //   Sum: '<S7>/Sum1'

    tmp_4 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_C[i], tmp_4);
    tmp_2 = _mm_loadu_pd(&rtb_x_0[i]);
    tmp_3 = _mm_loadu_pd(&rtb_C[i]);
    (void)_mm_storeu_pd(&rtb_C[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
      (rtb_Sum1[0])), tmp_3));
    tmp_2 = _mm_loadu_pd(&rtb_x_0[i + 2]);
    tmp_3 = _mm_loadu_pd(&rtb_C[i]);
    (void)_mm_storeu_pd(&rtb_C[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
      (rtb_Sum1[1])), tmp_3));
    tmp_2 = _mm_loadu_pd(&rtb_x_0[i + 4]);
    tmp_3 = _mm_loadu_pd(&rtb_C[i]);
    (void)_mm_storeu_pd(&rtb_C[i], _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd(u)),
      tmp_3));

    // Product: '<S83>/A[k]*xhat[k|k-1]' incorporates:
    //   Delay: '<S64>/MemoryX'
    //   Product: '<S83>/B[k]*u[k]'

    (void)_mm_storeu_pd(&rtb_N[i], tmp_4);
    tmp_4 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i]);
    tmp_2 = _mm_loadu_pd(&rtb_N[i]);
    (void)_mm_storeu_pd(&rtb_N[i], _mm_add_pd(_mm_mul_pd(tmp_4, _mm_set1_pd
      (localDW->MemoryX_DSTATE[0])), tmp_2));
    tmp_4 = _mm_loadu_pd(&rtb_TmpSignalConversionAtSFunct[i + 2]);
    tmp_2 = _mm_loadu_pd(&rtb_N[i]);
    (void)_mm_storeu_pd(&rtb_N[i], _mm_add_pd(_mm_mul_pd(tmp_4, _mm_set1_pd
      (localDW->MemoryX_DSTATE[1])), tmp_2));
  }

  // Update for Delay: '<S64>/MemoryX' incorporates:
  //   Sum: '<S83>/Add'

  localDW->MemoryX_DSTATE[0] = (rtb_C[0] + rtb_N[0]) + localDW->Product3[0];
  localDW->MemoryX_DSTATE[1] = (rtb_C[1] + rtb_N[1]) + localDW->Product3[1];

  // Update for Delay: '<S64>/MemoryP'
  localDW->icLoad_f = false;
  localDW->MemoryP_DSTATE[0] = rtb_y_k[0];
  localDW->MemoryP_DSTATE[1] = rtb_y_k[1];
  localDW->MemoryP_DSTATE[2] = rtb_y_k[2];
  localDW->MemoryP_DSTATE[3] = rtb_y_k[3];

  // Update for RandomNumber: '<S2>/Measurement Noise'
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];
}

//
// Output and update for atomic system:
//    '<S111>/MATLAB Function'
//    '<S255>/MATLAB Function'
//
void SupervisoryController::MATLABFunction(real_T rtu_A1, const real_T rtu_B1[3],
  real_T rtu_C1, const real_T rtu_D1[3], real_T rtu_A2, const real_T rtu_B2[3],
  real_T rtu_C2, const real_T rtu_D2[3], real_T rty_A[4], real_T rty_B[6],
  real_T rty_C[4], real_T rty_D[6])
{
  // MATLAB Function 'SupervisoryController/State1.controlLaw.AMPC1/Model Estimator/MATLAB Function': '<S143>:1' 
  // '<S143>:1:3' A = zeros(2,2);
  // '<S143>:1:4' B = zeros(2,3);
  // '<S143>:1:5' C = zeros(2,2);
  // '<S143>:1:6' D = zeros(2,3);
  // '<S143>:1:8' A = [A1, 0;
  // '<S143>:1:9'      0, A2];
  rty_A[0] = rtu_A1;
  rty_A[2] = 0.0;
  rty_A[1] = 0.0;
  rty_A[3] = rtu_A2;

  // '<S143>:1:10' B = [B1; B2];
  // '<S143>:1:11' C = [C1, 0;
  // '<S143>:1:12'      0, C2];
  rty_C[0] = rtu_C1;
  rty_C[2] = 0.0;
  rty_C[1] = 0.0;
  rty_C[3] = rtu_C2;

  // '<S143>:1:13' D = [D1; D2];
  rty_B[0] = rtu_B1[0];
  rty_B[1] = rtu_B2[0];
  rty_D[0] = rtu_D1[0];
  rty_D[1] = rtu_D2[0];
  rty_B[2] = rtu_B1[1];
  rty_B[3] = rtu_B2[1];
  rty_D[2] = rtu_D1[1];
  rty_D[3] = rtu_D2[1];
  rty_B[4] = rtu_B1[2];
  rty_B[5] = rtu_B2[2];
  rty_D[4] = rtu_D1[2];
  rty_D[5] = rtu_D2[2];

  //  disp(A);
  //  disp(B);
  //  disp(C);
  //  disp(D);
}

//
// Output and update for atomic system:
//    '<S144>/Model Type Converter'
//    '<S145>/Model Type Converter'
//    '<S288>/Model Type Converter'
//    '<S289>/Model Type Converter'
//
void SupervisoryController::ModelTypeConverter(const real_T rtu_num[6], const
  real_T rtu_den[2], real_T *rty_A, real_T rty_B[3], real_T *rty_C, real_T
  rty_D[3])
{
  real_T num[6];
  real_T den_idx_0;
  real_T den_idx_1;
  for (int32_T i{0}; i < 6; i++) {
    num[i] = rtu_num[i];
  }

  den_idx_0 = rtu_den[0];
  den_idx_1 = rtu_den[1];

  //  idModelTypeConverterTF2SS Convert a MISO transfer function described with
  //  variable z^-1 to state-space representation
  // MATLAB Function 'Estimators/Model Type Converter/Model Type Converter': '<S148>:1' 
  //    Copyright 2013-2016 The MathWorks, Inc.
  //  #codegen
  // '<S148>:1:8' [A,B,C,D] = idModelTypeConverterTF2SS(num,den);
  if (rtu_den[0] == 0.0) {
    den_idx_0 = 2.2204460492503131E-16;
  }

  if (den_idx_0 != 1.0) {
    for (int32_T i{0}; i <= 4; i += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&rtu_num[i]);
      (void)_mm_storeu_pd(&num[i], _mm_div_pd(tmp, _mm_set1_pd(den_idx_0)));
    }

    den_idx_1 = rtu_den[1] / den_idx_0;
  }

  *rty_C = 1.0;
  *rty_A = -den_idx_1;
  rty_B[0] = -den_idx_1 * num[0] + num[3];
  rty_D[0] = num[0];
  rty_B[1] = -den_idx_1 * num[1] + num[4];
  rty_D[1] = num[1];
  rty_B[2] = -den_idx_1 * num[2] + num[5];
  rty_D[2] = num[2];
}

// Function for MATLAB Function: '<S169>/rarx'
real_T SupervisoryController::xnrm2_j(int32_T n, const real_T x[5], int32_T ix0)
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

// Function for MATLAB Function: '<S169>/rarx'
real_T SupervisoryController::qrFactor_b(const real_T A[4], const real_T S[16],
  real_T Ns)
{
  real_T b_A[5];
  real_T b_S;
  real_T s;
  int32_T aoffset;
  for (int32_T i{0}; i < 4; i++) {
    aoffset = i << 2UL;
    b_A[i] = ((S[aoffset + 1] * A[1] + S[aoffset] * A[0]) + S[aoffset + 2] * A[2])
      + S[aoffset + 3] * A[3];
  }

  b_A[4] = Ns;
  b_S = b_A[0];
  s = xnrm2_j(4, b_A, 2);
  if (s != 0.0) {
    s = rt_hypotd_snf(b_A[0], s);
    if (b_A[0] >= 0.0) {
      s = -s;
    }

    if (std::abs(s) < 1.0020841800044864E-292) {
      aoffset = 0;
      do {
        aoffset++;
        b_A[1] *= 9.9792015476736E+291;
        b_A[2] *= 9.9792015476736E+291;
        b_A[3] *= 9.9792015476736E+291;
        b_A[4] *= 9.9792015476736E+291;
        s *= 9.9792015476736E+291;
        b_S *= 9.9792015476736E+291;
      } while ((std::abs(s) < 1.0020841800044864E-292) && (aoffset < 20));

      s = rt_hypotd_snf(b_S, xnrm2_j(4, b_A, 2));
      if (b_S >= 0.0) {
        s = -s;
      }

      for (int32_T i{0}; i < aoffset; i++) {
        s *= 1.0020841800044864E-292;
      }

      b_S = s;
    } else {
      b_S = s;
    }
  }

  return b_S;
}

// Function for MATLAB Function: '<S169>/rarx'
void SupervisoryController::trisolve_b(real_T A, real_T B_3[4])
{
  if (B_3[0] != 0.0) {
    B_3[0] /= A;
  }

  if (B_3[1] != 0.0) {
    B_3[1] /= A;
  }

  if (B_3[2] != 0.0) {
    B_3[2] /= A;
  }

  if (B_3[3] != 0.0) {
    B_3[3] /= A;
  }
}

// Function for MATLAB Function: '<S169>/rarx'
real_T SupervisoryController::xnrm2_ji(int32_T n, const real_T x[20], int32_T
  ix0)
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

// Function for MATLAB Function: '<S169>/rarx'
void SupervisoryController::xgemv_mj(int32_T m, int32_T n, const real_T A[20],
  int32_T ia0, const real_T x[20], int32_T ix0, real_T y[4])
{
  if ((m != 0) && (n != 0)) {
    int32_T b;
    if (n - 1 >= 0) {
      (void)std::memset(&y[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    b = (n - 1) * 5 + ia0;
    for (int32_T b_iy{ia0}; b_iy <= b; b_iy += 5) {
      real_T c;
      int32_T d;
      int32_T iyend;
      c = 0.0;
      d = b_iy + m;
      for (iyend = b_iy; iyend < d; iyend++) {
        c += x[((ix0 + iyend) - b_iy) - 1] * A[iyend - 1];
      }

      iyend = div_nde_s32_floor(b_iy - ia0, 5);
      y[iyend] += c;
    }
  }
}

// Function for MATLAB Function: '<S169>/rarx'
void SupervisoryController::xgerc_c(int32_T m, int32_T n, real_T alpha1, int32_T
  ix0, const real_T y[4], real_T A[20], int32_T ia0)
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
        b = m + jA;
        for (int32_T ijA{jA}; ijA < b; ijA++) {
          A[ijA - 1] += A[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += 5;
    }
  }
}

// Function for MATLAB Function: '<S169>/rarx'
void SupervisoryController::sqrtMeasurementUpdate_g(real_T L[16], const real_T
  H[4], real_T a0, real_T K[4])
{
  real_T c_A[20];
  real_T A[16];
  real_T y[16];
  real_T B_2[4];
  real_T Pxy[4];
  real_T K_0;
  real_T K_idx_1;
  real_T K_idx_2;
  real_T Sy;
  real_T prodVal;
  int32_T aoffset;
  int32_T b_i;
  int32_T coffset;
  int32_T ii;
  for (b_i = 0; b_i < 4; b_i++) {
    Pxy[b_i] = 0.0;
    for (ii = 0; ii < 4; ii++) {
      coffset = (ii << 2UL) + b_i;
      A[coffset] = 0.0;
      A[coffset] += L[b_i] * L[ii];
      A[coffset] += L[b_i + 4] * L[ii + 4];
      A[coffset] += L[b_i + 8] * L[ii + 8];
      A[coffset] += L[b_i + 12] * L[ii + 12];
      if (std::isinf(A[coffset])) {
        prodVal = A[coffset];
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        A[coffset] = K_0 * 1.7976931348623157E+308;
      }

      prodVal = A[coffset] * H[ii];
      if (std::isinf(prodVal)) {
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        prodVal = K_0 * 1.7976931348623157E+308;
      }

      prodVal += Pxy[b_i];
      if (std::isinf(prodVal)) {
        if (std::isnan(prodVal)) {
          K_0 = (rtNaN);
        } else if (prodVal < 0.0) {
          K_0 = -1.0;
        } else {
          K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        prodVal = K_0 * 1.7976931348623157E+308;
      }

      Pxy[b_i] = prodVal;
    }
  }

  prodVal = std::sqrt(a0);
  Sy = qrFactor_b(H, L, prodVal);
  B_2[0] = Pxy[0];
  B_2[1] = Pxy[1];
  B_2[2] = Pxy[2];
  B_2[3] = Pxy[3];
  trisolve_b(Sy, B_2);
  K[0] = B_2[0];
  K[1] = B_2[1];
  K[2] = B_2[2];
  K[3] = B_2[3];
  trisolve_b(Sy, K);
  K_0 = K[0];
  if (std::isinf(K[0])) {
    if (std::isnan(K[0])) {
      K_0 = (rtNaN);
    } else if (K[0] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[0] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  Sy = -K_0;
  K[0] = K_0;
  K_0 = K[1];
  if (std::isinf(K[1])) {
    if (std::isnan(K[1])) {
      K_0 = (rtNaN);
    } else if (K[1] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[1] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K_idx_1 = -K_0;
  K[1] = K_0;
  K_0 = K[2];
  if (std::isinf(K[2])) {
    if (std::isnan(K[2])) {
      K_0 = (rtNaN);
    } else if (K[2] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[2] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K_idx_2 = -K_0;
  K[2] = K_0;
  K_0 = K[3];
  if (std::isinf(K[3])) {
    if (std::isnan(K[3])) {
      K_0 = (rtNaN);
    } else if (K[3] < 0.0) {
      K_0 = -1.0;
    } else {
      K_0 = static_cast<real_T>(K[3] > 0.0 ? static_cast<int32_T>(1) :
        static_cast<int32_T>(0));
    }

    K_0 *= 1.7976931348623157E+308;
  }

  K[3] = K_0;
  b_i = 0;
  for (ii = 0; ii < 4; ii++) {
    A[b_i] = Sy * H[ii];
    A[b_i + 1] = K_idx_1 * H[ii];
    A[b_i + 2] = K_idx_2 * H[ii];
    A[b_i + 3] = -K_0 * H[ii];
    b_i += 4;
  }

  A[0]++;
  A[5]++;
  A[10]++;
  A[15]++;
  for (b_i = 0; b_i < 4; b_i++) {
    coffset = b_i << 2UL;
    for (ii = 0; ii < 4; ii++) {
      aoffset = ii << 2UL;
      y[coffset + ii] = ((L[aoffset + 1] * A[b_i + 4] + L[aoffset] * A[b_i]) +
                         L[aoffset + 2] * A[b_i + 8]) + L[aoffset + 3] * A[b_i +
        12];
    }
  }

  b_i = 0;
  coffset = 0;
  for (ii = 0; ii < 4; ii++) {
    c_A[b_i] = y[coffset];
    c_A[b_i + 1] = y[coffset + 1];
    c_A[b_i + 2] = y[coffset + 2];
    c_A[b_i + 3] = y[coffset + 3];
    c_A[b_i + 4] = K[ii] * prodVal;
    Pxy[ii] = 0.0;
    b_i += 5;
    coffset += 4;
  }

  for (b_i = 0; b_i < 4; b_i++) {
    int32_T lastv;
    int32_T scalarLB;
    ii = b_i * 5 + b_i;
    Sy = c_A[ii];
    lastv = ii + 2;
    B_2[b_i] = 0.0;
    prodVal = xnrm2_ji(4 - b_i, c_A, ii + 2);
    if (prodVal != 0.0) {
      prodVal = rt_hypotd_snf(c_A[ii], prodVal);
      if (c_A[ii] >= 0.0) {
        prodVal = -prodVal;
      }

      if (std::abs(prodVal) < 1.0020841800044864E-292) {
        __m128d tmp;
        int32_T scalarLB_tmp;
        int32_T vectorUB;
        coffset = 0;
        scalarLB = (ii - b_i) + 5;
        do {
          coffset++;
          scalarLB_tmp = (((((scalarLB - ii) - 1) / 2) << 1UL) + ii) + 2;
          vectorUB = scalarLB_tmp - 2;
          for (aoffset = lastv; aoffset <= vectorUB; aoffset += 2) {
            tmp = _mm_loadu_pd(&c_A[aoffset - 1]);
            (void)_mm_storeu_pd(&c_A[aoffset - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (aoffset = scalarLB_tmp; aoffset <= scalarLB; aoffset++) {
            c_A[aoffset - 1] *= 9.9792015476736E+291;
          }

          prodVal *= 9.9792015476736E+291;
          Sy *= 9.9792015476736E+291;
        } while ((std::abs(prodVal) < 1.0020841800044864E-292) && (coffset < 20));

        prodVal = rt_hypotd_snf(Sy, xnrm2_ji(4 - b_i, c_A, ii + 2));
        if (Sy >= 0.0) {
          prodVal = -prodVal;
        }

        B_2[b_i] = (prodVal - Sy) / prodVal;
        Sy = 1.0 / (Sy - prodVal);
        vectorUB = scalarLB_tmp - 2;
        for (aoffset = lastv; aoffset <= vectorUB; aoffset += 2) {
          tmp = _mm_loadu_pd(&c_A[aoffset - 1]);
          (void)_mm_storeu_pd(&c_A[aoffset - 1], _mm_mul_pd(tmp, _mm_set1_pd(Sy)));
        }

        for (aoffset = scalarLB_tmp; aoffset <= scalarLB; aoffset++) {
          c_A[aoffset - 1] *= Sy;
        }

        for (lastv = 0; lastv < coffset; lastv++) {
          prodVal *= 1.0020841800044864E-292;
        }

        Sy = prodVal;
      } else {
        int32_T vectorUB;
        B_2[b_i] = (prodVal - c_A[ii]) / prodVal;
        Sy = 1.0 / (c_A[ii] - prodVal);
        aoffset = (ii - b_i) + 5;
        scalarLB = (((((aoffset - ii) - 1) / 2) << 1UL) + ii) + 2;
        vectorUB = scalarLB - 2;
        for (coffset = lastv; coffset <= vectorUB; coffset += 2) {
          __m128d tmp;
          tmp = _mm_loadu_pd(&c_A[coffset - 1]);
          (void)_mm_storeu_pd(&c_A[coffset - 1], _mm_mul_pd(tmp, _mm_set1_pd(Sy)));
        }

        for (coffset = scalarLB; coffset <= aoffset; coffset++) {
          c_A[coffset - 1] *= Sy;
        }

        Sy = prodVal;
      }
    }

    c_A[ii] = Sy;
    if (b_i + 1 < 4) {
      prodVal = c_A[ii];
      c_A[ii] = 1.0;
      if (B_2[b_i] != 0.0) {
        boolean_T exitg2;
        lastv = 5 - b_i;
        coffset = (ii - b_i) + 4;
        while ((lastv > 0) && (c_A[coffset] == 0.0)) {
          lastv--;
          coffset--;
        }

        coffset = 3 - b_i;
        exitg2 = false;
        while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                == false) && (coffset > 0)) {
          int32_T exitg1;
          aoffset = ((coffset - 1) * 5 + ii) + 5;
          scalarLB = aoffset;
          do {
            exitg1 = 0;
            if (scalarLB + 1 <= aoffset + lastv) {
              if (c_A[scalarLB] != 0.0) {
                exitg1 = 1;
              } else {
                scalarLB++;
              }
            } else {
              coffset--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        lastv = 0;
        coffset = 0;
      }

      if (lastv > 0) {
        xgemv_mj(lastv, coffset, c_A, ii + 6, c_A, ii + 1, Pxy);
        xgerc_c(lastv, coffset, -B_2[b_i], ii + 1, Pxy, c_A, ii + 6);
      }

      c_A[ii] = prodVal;
    }
  }

  for (b_i = 0; b_i < 4; b_i++) {
    for (ii = 0; ii <= b_i; ii++) {
      A[ii + (b_i << 2UL)] = c_A[5 * b_i + ii];
    }

    for (ii = b_i + 2; ii < 5; ii++) {
      A[(ii + (b_i << 2UL)) - 1] = 0.0;
    }
  }

  for (ii = 0; ii < 4; ii++) {
    b_i = ii << 2UL;
    L[b_i] = A[ii];
    L[b_i + 1] = A[ii + 4];
    L[b_i + 2] = A[ii + 8];
    L[b_i + 3] = A[ii + 12];
    prodVal = L[b_i];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 1];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 1] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 2];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 2] = K_0 * 1.7976931348623157E+308;
    }

    prodVal = L[b_i + 3];
    if (std::isinf(prodVal)) {
      if (std::isnan(prodVal)) {
        K_0 = (rtNaN);
      } else if (prodVal < 0.0) {
        K_0 = -1.0;
      } else {
        K_0 = static_cast<real_T>(prodVal > 0.0 ? static_cast<int32_T>(1) :
          static_cast<int32_T>(0));
      }

      L[b_i + 3] = K_0 * 1.7976931348623157E+308;
    }
  }
}

//
// Output and update for atomic system:
//    '<S169>/rarx'
//    '<S198>/rarx'
//    '<S313>/rarx'
//    '<S342>/rarx'
//
void SupervisoryController::rarx(const real_T rtu_uMeas[3], real_T rtu_yMeas,
  boolean_T rtu_isEnabled, real_T rtu_adg1, real_T rty_A[2], real_T rty_B[6],
  real_T *rty_e, real_T rty_theta[4], real_T rty_L[16], real_T rty_phiMemory[7],
  DW_rarx *localDW)
{
  real_T a[4];
  real_T rty_phiMemory_0[4];
  real_T b;

  //    % Was a reset trigered?
  // MATLAB Function 'Recursive Estimation/rarx/rarx': '<S177>:1'
  //    Copyright 2015-2018 The MathWorks, Inc.
  //  Initialization and reset
  // '<S177>:1:17' if isempty(arxEstimator)
  if (!localDW->arxEstimator_not_empty) {
    //  Get the estimator object
    // '<S177>:1:19' arxEstimator = recursiveEstimation.internal.rpoly.getEstimator('arx',algorithmParams); 
    localDW->arxEstimator_not_empty = true;

    // '<S177>:1:20' [theta, L] = arxEstimator.reset(phiMemory, uMemory, yMemory, theta, L); 
  }

  //  Update algorithm states, perform estimation
  // '<S177>:1:27' [phiMemory, uMemory, yMemory, ...
  // '<S177>:1:28'     theta, L, ...
  // '<S177>:1:29'     A, B, e] = ...
  // '<S177>:1:30'     arxEstimator.estimate(...
  // '<S177>:1:31'     phiMemory, uMemory, yMemory, theta, L, ... % States
  // '<S177>:1:32'     yMeas, uMeas,... % new IO data
  // '<S177>:1:33'     isEnabled, adg1, adg2, algorithmParams);
  localDW->arxEstimator.OutputDataIterator.IteratorPosition = 1;
  localDW->arxEstimator.InputDataIterator.IteratorPosition = 1;
  if (localDW->arxEstimator.OutputDataIterator.IteratorPosition + 1 <= 2) {
    localDW->arxEstimator.OutputDataIterator.IteratorPosition++;
  } else {
    localDW->arxEstimator.OutputDataIterator.IteratorPosition = 2;
  }

  if (localDW->arxEstimator.InputDataIterator.IteratorPosition + 1 <= 2) {
    localDW->arxEstimator.InputDataIterator.IteratorPosition++;
  } else {
    localDW->arxEstimator.InputDataIterator.IteratorPosition = 2;
  }

  rty_phiMemory[1] = rtu_uMeas[0];
  rty_phiMemory[3] = rtu_uMeas[1];
  rty_phiMemory[5] = rtu_uMeas[2];
  b = rtu_yMeas - (((rty_phiMemory[0] * rty_theta[0] + rty_theta[1] *
                     rty_phiMemory[2]) + rty_theta[2] * rty_phiMemory[4]) +
                   rty_theta[3] * rty_phiMemory[6]);
  if (rtu_isEnabled) {
    real_T sqrtFF;
    rty_phiMemory_0[0] = rty_phiMemory[0];
    rty_phiMemory_0[1] = rty_phiMemory[2];
    rty_phiMemory_0[2] = rty_phiMemory[4];
    rty_phiMemory_0[3] = rty_phiMemory[6];
    sqrtMeasurementUpdate_g(rty_L, rty_phiMemory_0, rtu_adg1, a);
    sqrtFF = std::sqrt(rtu_adg1);
    for (int32_T kk{0}; kk < 4; kk++) {
      for (int32_T b_kk2{0}; b_kk2 <= kk; b_kk2++) {
        int32_T tmp;
        tmp = ((((3 - kk) << 2UL) + b_kk2) - kk) + 3;
        rty_L[tmp] /= sqrtFF;
      }
    }

    for (int32_T kk{0}; kk < 4; kk++) {
      if (std::isinf(rty_L[kk])) {
        sqrtFF = rty_L[kk];
        if (std::isnan(sqrtFF)) {
          sqrtFF = (rtNaN);
        } else if (sqrtFF < 0.0) {
          sqrtFF = -1.0;
        } else {
          sqrtFF = static_cast<real_T>(sqrtFF > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rty_L[kk] = sqrtFF * 1.7976931348623157E+308;
      }

      sqrtFF = rty_L[kk + 4];
      if (std::isinf(sqrtFF)) {
        if (std::isnan(sqrtFF)) {
          sqrtFF = (rtNaN);
        } else if (sqrtFF < 0.0) {
          sqrtFF = -1.0;
        } else {
          sqrtFF = static_cast<real_T>(sqrtFF > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rty_L[kk + 4] = sqrtFF * 1.7976931348623157E+308;
      }

      sqrtFF = rty_L[kk + 8];
      if (std::isinf(sqrtFF)) {
        if (std::isnan(sqrtFF)) {
          sqrtFF = (rtNaN);
        } else if (sqrtFF < 0.0) {
          sqrtFF = -1.0;
        } else {
          sqrtFF = static_cast<real_T>(sqrtFF > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rty_L[kk + 8] = sqrtFF * 1.7976931348623157E+308;
      }

      sqrtFF = rty_L[kk + 12];
      if (std::isinf(sqrtFF)) {
        if (std::isnan(sqrtFF)) {
          sqrtFF = (rtNaN);
        } else if (sqrtFF < 0.0) {
          sqrtFF = -1.0;
        } else {
          sqrtFF = static_cast<real_T>(sqrtFF > 0.0 ? static_cast<int32_T>(1) :
            static_cast<int32_T>(0));
        }

        rty_L[kk + 12] = sqrtFF * 1.7976931348623157E+308;
      }

      rty_theta[kk] += a[kk] * b;
    }
  }

  rty_A[0] = 1.0;
  rty_A[1] = rty_theta[0];
  for (int32_T kk{0}; kk < 6; kk++) {
    rty_phiMemory[6 - kk] = rty_phiMemory[5 - kk];
    rty_B[kk] = 0.0;
  }

  rty_phiMemory[0] = -rtu_yMeas;
  rty_B[3] = rty_theta[1];
  rty_B[4] = rty_theta[2];
  rty_B[5] = rty_theta[3];
  *rty_e = b;

  //  % States
  //  % new IO data
  //  helper variables
}

// Function for MATLAB Function: '<S210>/Discrete-Time KF - Calculate PLMZ'
void SupervisoryController::mrdiv(const real_T A[8], const real_T B_4[4], real_T
  Y[8])
{
  real_T a21;
  real_T a22;
  real_T a22_tmp;
  int32_T Y_tmp;
  int32_T r1;
  int32_T r2;
  if (std::abs(B_4[1]) > std::abs(B_4[0])) {
    r1 = 1;
    r2 = 0;
  } else {
    r1 = 0;
    r2 = 1;
  }

  a21 = B_4[r2] / B_4[r1];
  a22_tmp = B_4[r1 + 2];
  a22 = B_4[r2 + 2] - a22_tmp * a21;
  Y_tmp = r1 << 2UL;
  Y[Y_tmp] = A[0] / B_4[r1];
  r2 <<= 2UL;
  Y[r2] = (A[4] - Y[Y_tmp] * a22_tmp) / a22;
  Y[Y_tmp] -= Y[r2] * a21;
  Y[Y_tmp + 1] = A[1] / B_4[r1];
  Y[r2 + 1] = (A[5] - Y[Y_tmp + 1] * a22_tmp) / a22;
  Y[Y_tmp + 1] -= Y[r2 + 1] * a21;
  Y[Y_tmp + 2] = A[2] / B_4[r1];
  Y[r2 + 2] = (A[6] - Y[Y_tmp + 2] * a22_tmp) / a22;
  Y[Y_tmp + 2] -= Y[r2 + 2] * a21;
  Y[Y_tmp + 3] = A[3] / B_4[r1];
  Y[r2 + 3] = (A[7] - Y[Y_tmp + 3] * a22_tmp) / a22;
  Y[Y_tmp + 3] -= Y[r2 + 3] * a21;
}

//
// Output and update for atomic system:
//    '<S208>/CalculatePL'
//    '<S352>/CalculatePL'
//
void SupervisoryController::CalculatePL(const real_T rtu_Ak[16], const real_T
  rtu_Ck[8], const real_T rtu_Qbark[16], const real_T rtu_Rbark[4], const real_T
  rtu_Nbark[8], boolean_T rtu_Enablek, const real_T rtu_Pk[16], real_T rty_Mk[8],
  real_T rty_Lk[8], real_T rty_Zk[16], real_T rty_Pk1[16])
{
  real_T Abar[16];
  real_T Abar_0[16];
  real_T Abar_1[16];
  real_T rty_Mk_0[16];
  real_T NRinv[8];
  real_T rtu_Ck_0[8];
  real_T yCov[4];
  int8_T b_I[16];

  // MATLAB Function: '<S210>/Discrete-Time KF - Calculate PLMZ'
  //  See help of ctrlKalmanFilterDTCalculatePL.m
  // MATLAB Function 'KalmanFilterUtilities/DTCalculatePL/Discrete-Time KF - Calculate PLMZ': '<S248>:1' 
  //    Copyright 2014 The MathWorks, Inc.
  // '<S248>:1:7' [L,M,Z,PNew] = ctrlKalmanFilterDTCalculatePL(A,C,Q,R,N,P,isEnabled); 
  if (rtu_Enablek) {
    __m128d tmp;
    __m128d tmp_0;
    int32_T i;
    int32_T i_0;
    int32_T rtu_Ak_tmp;
    int32_T rtu_Ak_tmp_tmp;
    int32_T yCov_tmp;
    i = 0;
    for (i_0 = 0; i_0 < 2; i_0++) {
      rtu_Ak_tmp = 0;
      rtu_Ak_tmp_tmp = 0;
      for (yCov_tmp = 0; yCov_tmp < 4; yCov_tmp++) {
        int32_T NRinv_tmp;
        NRinv_tmp = rtu_Ak_tmp + i_0;
        NRinv[yCov_tmp + i] = rtu_Ck[NRinv_tmp];
        rtu_Ck_0[NRinv_tmp] = 0.0;
        rtu_Ck_0[NRinv_tmp] += rtu_Pk[rtu_Ak_tmp_tmp] * rtu_Ck[i_0];
        rtu_Ck_0[NRinv_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 1] * rtu_Ck[i_0 + 2];
        rtu_Ck_0[NRinv_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 2] * rtu_Ck[i_0 + 4];
        rtu_Ck_0[NRinv_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 3] * rtu_Ck[i_0 + 6];
        rtu_Ak_tmp += 2;
        rtu_Ak_tmp_tmp += 4;
      }

      i += 4;
    }

    for (i = 0; i < 2; i++) {
      i_0 = 0;
      rtu_Ak_tmp = 0;
      for (rtu_Ak_tmp_tmp = 0; rtu_Ak_tmp_tmp < 2; rtu_Ak_tmp_tmp++) {
        yCov_tmp = i_0 + i;
        yCov[yCov_tmp] = (((NRinv[rtu_Ak_tmp + 1] * rtu_Ck_0[i + 2] +
                            NRinv[rtu_Ak_tmp] * rtu_Ck_0[i]) + NRinv[rtu_Ak_tmp
                           + 2] * rtu_Ck_0[i + 4]) + NRinv[rtu_Ak_tmp + 3] *
                          rtu_Ck_0[i + 6]) + rtu_Rbark[yCov_tmp];
        i_0 += 2;
        rtu_Ak_tmp += 4;
      }
    }

    for (i = 0; i < 4; i++) {
      for (i_0 = 0; i_0 < 4; i_0++) {
        rtu_Ak_tmp_tmp = i_0 << 2UL;
        rtu_Ak_tmp = i + rtu_Ak_tmp_tmp;
        Abar[rtu_Ak_tmp] = 0.0;
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp] * rtu_Ak[i];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 1] * rtu_Ak[i + 4];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 2] * rtu_Ak[i + 8];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 3] * rtu_Ak[i + 12];
      }

      for (i_0 = 0; i_0 < 2; i_0++) {
        rtu_Ak_tmp = i_0 << 2UL;
        rtu_Ak_tmp_tmp = rtu_Ak_tmp + i;
        rtu_Ck_0[rtu_Ak_tmp_tmp] = (((NRinv[rtu_Ak_tmp + 1] * Abar[i + 4] +
          NRinv[rtu_Ak_tmp] * Abar[i]) + NRinv[rtu_Ak_tmp + 2] * Abar[i + 8]) +
          NRinv[rtu_Ak_tmp + 3] * Abar[i + 12]) + rtu_Nbark[rtu_Ak_tmp_tmp];
      }
    }

    mrdiv(rtu_Ck_0, yCov, rty_Lk);
    i = 0;
    for (i_0 = 0; i_0 < 2; i_0++) {
      for (rtu_Ak_tmp = 0; rtu_Ak_tmp <= 2; rtu_Ak_tmp += 2) {
        rtu_Ak_tmp_tmp = rtu_Ak_tmp + i;
        (void)_mm_storeu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&rtu_Pk[rtu_Ak_tmp]);
        tmp_0 = _mm_loadu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp]);
        (void)_mm_storeu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp], _mm_add_pd(tmp_0,
          _mm_mul_pd(_mm_set1_pd(NRinv[i]), tmp)));
        tmp = _mm_loadu_pd(&rtu_Pk[rtu_Ak_tmp + 4]);
        tmp_0 = _mm_loadu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp]);
        (void)_mm_storeu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp], _mm_add_pd(_mm_mul_pd
          (_mm_set1_pd(NRinv[i + 1]), tmp), tmp_0));
        tmp = _mm_loadu_pd(&rtu_Pk[rtu_Ak_tmp + 8]);
        tmp_0 = _mm_loadu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp]);
        (void)_mm_storeu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp], _mm_add_pd(_mm_mul_pd
          (_mm_set1_pd(NRinv[i + 2]), tmp), tmp_0));
        tmp = _mm_loadu_pd(&rtu_Pk[rtu_Ak_tmp + 12]);
        tmp_0 = _mm_loadu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp]);
        (void)_mm_storeu_pd(&rtu_Ck_0[rtu_Ak_tmp_tmp], _mm_add_pd(_mm_mul_pd
          (_mm_set1_pd(NRinv[i + 3]), tmp), tmp_0));
      }

      i += 4;
    }

    mrdiv(rtu_Ck_0, yCov, rty_Mk);
    for (i = 0; i < 16; i++) {
      b_I[i] = 0;
    }

    b_I[0] = 1;
    b_I[5] = 1;
    b_I[10] = 1;
    b_I[15] = 1;
    for (i = 0; i < 4; i++) {
      for (i_0 = 0; i_0 < 4; i_0++) {
        rtu_Ak_tmp_tmp = i_0 << 1UL;
        rtu_Ak_tmp = (i_0 << 2UL) + i;
        Abar[rtu_Ak_tmp] = static_cast<real_T>(b_I[rtu_Ak_tmp]) -
          (rtu_Ck[rtu_Ak_tmp_tmp + 1] * rty_Mk[i + 4] + rtu_Ck[rtu_Ak_tmp_tmp] *
           rty_Mk[i]);
      }

      for (i_0 = 0; i_0 < 4; i_0++) {
        rtu_Ak_tmp = i_0 << 2UL;
        rtu_Ak_tmp_tmp = i + rtu_Ak_tmp;
        Abar_0[rtu_Ak_tmp_tmp] = 0.0;
        Abar_0[rtu_Ak_tmp_tmp] += rtu_Pk[rtu_Ak_tmp] * Abar[i];
        Abar_0[rtu_Ak_tmp_tmp] += rtu_Pk[rtu_Ak_tmp + 1] * Abar[i + 4];
        Abar_0[rtu_Ak_tmp_tmp] += rtu_Pk[rtu_Ak_tmp + 2] * Abar[i + 8];
        Abar_0[rtu_Ak_tmp_tmp] += rtu_Pk[rtu_Ak_tmp + 3] * Abar[i + 12];
      }

      real_T rty_Mk_tmp;
      NRinv[i] = 0.0;
      NRinv[i] += rty_Mk[i] * rtu_Rbark[0];
      rty_Mk_tmp = rty_Mk[i + 4];
      NRinv[i] += rty_Mk_tmp * rtu_Rbark[1];
      NRinv[i + 4] = 0.0;
      NRinv[i + 4] += rty_Mk[i] * rtu_Rbark[2];
      NRinv[i + 4] += rty_Mk_tmp * rtu_Rbark[3];
    }

    for (i = 0; i < 4; i++) {
      i_0 = 0;
      for (rtu_Ak_tmp = 0; rtu_Ak_tmp < 4; rtu_Ak_tmp++) {
        rtu_Ak_tmp_tmp = i_0 + i;
        Abar_1[rtu_Ak_tmp_tmp] = 0.0;
        Abar_1[rtu_Ak_tmp_tmp] += Abar_0[i] * Abar[rtu_Ak_tmp];
        Abar_1[rtu_Ak_tmp_tmp] += Abar_0[i + 4] * Abar[rtu_Ak_tmp + 4];
        Abar_1[rtu_Ak_tmp_tmp] += Abar_0[i + 8] * Abar[rtu_Ak_tmp + 8];
        Abar_1[rtu_Ak_tmp_tmp] += Abar_0[i + 12] * Abar[rtu_Ak_tmp + 12];
        rty_Mk_0[rtu_Ak_tmp_tmp] = 0.0;
        rty_Mk_0[rtu_Ak_tmp_tmp] += NRinv[i] * rty_Mk[rtu_Ak_tmp];
        rty_Mk_0[rtu_Ak_tmp_tmp] += NRinv[i + 4] * rty_Mk[rtu_Ak_tmp + 4];
        i_0 += 4;
      }
    }

    for (i = 0; i <= 14; i += 2) {
      tmp = _mm_loadu_pd(&Abar_1[i]);
      tmp_0 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Zk[i], _mm_add_pd(tmp, tmp_0));
    }

    mrdiv(rtu_Nbark, rtu_Rbark, NRinv);
    for (i = 0; i < 4; i++) {
      for (i_0 = 0; i_0 < 4; i_0++) {
        rtu_Ak_tmp_tmp = i_0 << 1UL;
        rtu_Ak_tmp = (i_0 << 2UL) + i;
        Abar[rtu_Ak_tmp] = rtu_Ak[rtu_Ak_tmp] - (rtu_Ck[rtu_Ak_tmp_tmp + 1] *
          NRinv[i + 4] + rtu_Ck[rtu_Ak_tmp_tmp] * NRinv[i]);
      }

      for (i_0 = 0; i_0 < 4; i_0++) {
        rtu_Ak_tmp = i_0 << 2UL;
        rtu_Ak_tmp_tmp = i + rtu_Ak_tmp;
        Abar_0[rtu_Ak_tmp_tmp] = 0.0;
        Abar_0[rtu_Ak_tmp_tmp] += rty_Zk[rtu_Ak_tmp] * Abar[i];
        Abar_0[rtu_Ak_tmp_tmp] += rty_Zk[rtu_Ak_tmp + 1] * Abar[i + 4];
        Abar_0[rtu_Ak_tmp_tmp] += rty_Zk[rtu_Ak_tmp + 2] * Abar[i + 8];
        Abar_0[rtu_Ak_tmp_tmp] += rty_Zk[rtu_Ak_tmp + 3] * Abar[i + 12];
      }
    }

    for (i = 0; i < 4; i++) {
      i_0 = 0;
      for (rtu_Ak_tmp = 0; rtu_Ak_tmp < 4; rtu_Ak_tmp++) {
        rtu_Ak_tmp_tmp = i_0 + i;
        Abar_1[rtu_Ak_tmp_tmp] = (((Abar_0[i + 4] * Abar[rtu_Ak_tmp + 4] +
          Abar_0[i] * Abar[rtu_Ak_tmp]) + Abar_0[i + 8] * Abar[rtu_Ak_tmp + 8])
          + Abar_0[i + 12] * Abar[rtu_Ak_tmp + 12]) + rtu_Qbark[rtu_Ak_tmp_tmp];
        rty_Mk_0[rtu_Ak_tmp_tmp] = 0.0;
        rty_Mk_0[rtu_Ak_tmp_tmp] += NRinv[i] * rtu_Nbark[rtu_Ak_tmp];
        rty_Mk_0[rtu_Ak_tmp_tmp] += NRinv[i + 4] * rtu_Nbark[rtu_Ak_tmp + 4];
        i_0 += 4;
      }
    }

    for (i = 0; i <= 14; i += 2) {
      tmp = _mm_loadu_pd(&Abar_1[i]);
      tmp_0 = _mm_loadu_pd(&rty_Mk_0[i]);
      (void)_mm_storeu_pd(&rty_Pk1[i], _mm_sub_pd(tmp, tmp_0));
    }
  } else {
    (void)std::memset(&rty_Lk[0], 0, sizeof(real_T) << 3UL);
    (void)std::memset(&rty_Mk[0], 0, sizeof(real_T) << 3UL);
    (void)std::memcpy(&rty_Zk[0], &rtu_Pk[0], sizeof(real_T) << 4UL);
    for (int32_T i{0}; i < 4; i++) {
      int32_T rtu_Ak_tmp;
      for (int32_T i_0{0}; i_0 < 4; i_0++) {
        int32_T rtu_Ak_tmp_tmp;
        rtu_Ak_tmp_tmp = i_0 << 2UL;
        rtu_Ak_tmp = i + rtu_Ak_tmp_tmp;
        Abar[rtu_Ak_tmp] = 0.0;
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp] * rtu_Ak[i];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 1] * rtu_Ak[i + 4];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 2] * rtu_Ak[i + 8];
        Abar[rtu_Ak_tmp] += rtu_Pk[rtu_Ak_tmp_tmp + 3] * rtu_Ak[i + 12];
      }

      for (int32_T i_0{0}; i_0 < 4; i_0++) {
        rtu_Ak_tmp = (i_0 << 2UL) + i;
        rty_Pk1[rtu_Ak_tmp] = (((Abar[i + 4] * rtu_Ak[i_0 + 4] + Abar[i] *
          rtu_Ak[i_0]) + Abar[i + 8] * rtu_Ak[i_0 + 8]) + Abar[i + 12] *
          rtu_Ak[i_0 + 12]) + rtu_Qbark[rtu_Ak_tmp];
      }
    }
  }

  // End of MATLAB Function: '<S210>/Discrete-Time KF - Calculate PLMZ'
}

//
// Output and update for atomic system:
//    '<S249>/SqrtUsedFcn'
//    '<S393>/SqrtUsedFcn'
//
void SupervisoryController::SqrtUsedFcn(const real_T rtu_u[16], boolean_T
  rtu_isSqrtUsed, real_T rty_P[16])
{
  //  Determine if the Square-Root algorithm was used
  // MATLAB Function 'Kalman Filter/CovarianceOutputConfigurator/decideOutput/SqrtUsedFcn': '<S250>:1' 
  // '<S250>:1:4' if isSqrtUsed
  if (rtu_isSqrtUsed) {
    // '<S250>:1:5' P = u*u.';
    for (int32_T i{0}; i < 4; i++) {
      int32_T tmp;
      tmp = 0;
      for (int32_T i_0{0}; i_0 < 4; i_0++) {
        int32_T tmp_0;
        tmp_0 = tmp + i;
        rty_P[tmp_0] = 0.0;
        rty_P[tmp_0] += rtu_u[i] * rtu_u[i_0];
        rty_P[tmp_0] += rtu_u[i + 4] * rtu_u[i_0 + 4];
        rty_P[tmp_0] += rtu_u[i + 8] * rtu_u[i_0 + 8];
        rty_P[tmp_0] += rtu_u[i + 12] * rtu_u[i_0 + 12];
        tmp += 4;
      }
    }
  } else {
    // '<S250>:1:6' else
    // '<S250>:1:7' P = u;
    (void)std::memcpy(&rty_P[0], &rtu_u[0], sizeof(real_T) << 4UL);
  }
}

//
// System initialize for enable system:
//    '<S227>/MeasurementUpdate'
//    '<S371>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Init(real_T rty_Lykyhatkk1[4],
  P_MeasurementUpdate *localP)
{
  // SystemInitialize for Outport: '<S251>/L*(y[k]-yhat[k|k-1])'
  rty_Lykyhatkk1[0] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[1] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[2] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[3] = localP->Lykyhatkk1_Y0;
}

//
// Disable for enable system:
//    '<S227>/MeasurementUpdate'
//    '<S371>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate_Disable(real_T rty_Lykyhatkk1[4],
  DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  // Outputs for Enabled SubSystem: '<S227>/MeasurementUpdate' incorporates:
  //   EnablePort: '<S251>/Enable'

  // Disable for Outport: '<S251>/L*(y[k]-yhat[k|k-1])'
  rty_Lykyhatkk1[0] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[1] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[2] = localP->Lykyhatkk1_Y0;
  rty_Lykyhatkk1[3] = localP->Lykyhatkk1_Y0;

  // End of Outputs for SubSystem: '<S227>/MeasurementUpdate'
  localDW->MeasurementUpdate_MODE = false;
}

//
// Output and update for enable system:
//    '<S227>/MeasurementUpdate'
//    '<S371>/MeasurementUpdate'
//
void SupervisoryController::MeasurementUpdate(boolean_T rtu_Enable, const real_T
  rtu_Lk[8], const real_T rtu_yk[2], const real_T rtu_Ck[8], const real_T
  rtu_xhatkk1[4], const real_T rtu_Dk[6], const real_T rtu_uk[3], real_T
  rty_Lykyhatkk1[4], DW_MeasurementUpdate *localDW, P_MeasurementUpdate *localP)
{
  real_T rtu_Ck_0[2];
  real_T rtu_Dk_0[2];
  real_T rtu_yk_0[2];

  // Outputs for Enabled SubSystem: '<S227>/MeasurementUpdate' incorporates:
  //   EnablePort: '<S251>/Enable'

  if (rtu_Enable) {
    __m128d tmp;
    __m128d tmp_1;
    localDW->MeasurementUpdate_MODE = true;
    for (int32_T i{0}; i <= 0; i += 2) {
      __m128d tmp_0;

      // Product: '<S251>/C[k]*xhat[k|k-1]'
      tmp_1 = _mm_set1_pd(0.0);
      (void)_mm_storeu_pd(&rtu_Ck_0[i], tmp_1);
      tmp = _mm_loadu_pd(&rtu_Ck[i]);
      tmp_0 = _mm_loadu_pd(&rtu_Ck_0[i]);
      (void)_mm_storeu_pd(&rtu_Ck_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
        (rtu_xhatkk1[0])), tmp_0));
      tmp = _mm_loadu_pd(&rtu_Ck[i + 2]);
      tmp_0 = _mm_loadu_pd(&rtu_Ck_0[i]);
      (void)_mm_storeu_pd(&rtu_Ck_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
        (rtu_xhatkk1[1])), tmp_0));
      tmp = _mm_loadu_pd(&rtu_Ck[i + 4]);
      tmp_0 = _mm_loadu_pd(&rtu_Ck_0[i]);
      (void)_mm_storeu_pd(&rtu_Ck_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
        (rtu_xhatkk1[2])), tmp_0));
      tmp = _mm_loadu_pd(&rtu_Ck[i + 6]);
      tmp_0 = _mm_loadu_pd(&rtu_Ck_0[i]);
      (void)_mm_storeu_pd(&rtu_Ck_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
        (rtu_xhatkk1[3])), tmp_0));

      // Product: '<S251>/D[k]*u[k]' incorporates:
      //   Product: '<S251>/C[k]*xhat[k|k-1]'

      (void)_mm_storeu_pd(&rtu_Dk_0[i], tmp_1);
      tmp_1 = _mm_loadu_pd(&rtu_Dk[i]);
      tmp = _mm_loadu_pd(&rtu_Dk_0[i]);
      (void)_mm_storeu_pd(&rtu_Dk_0[i], _mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd
        (rtu_uk[0])), tmp));
      tmp_1 = _mm_loadu_pd(&rtu_Dk[i + 2]);
      tmp = _mm_loadu_pd(&rtu_Dk_0[i]);
      (void)_mm_storeu_pd(&rtu_Dk_0[i], _mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd
        (rtu_uk[1])), tmp));
      tmp_1 = _mm_loadu_pd(&rtu_Dk[i + 4]);
      tmp = _mm_loadu_pd(&rtu_Dk_0[i]);
      (void)_mm_storeu_pd(&rtu_Dk_0[i], _mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd
        (rtu_uk[2])), tmp));

      // Product: '<S251>/C[k]*xhat[k|k-1]'
      tmp_1 = _mm_loadu_pd(&rtu_Ck_0[i]);

      // Product: '<S251>/D[k]*u[k]' incorporates:
      //   Product: '<S251>/C[k]*xhat[k|k-1]'

      tmp = _mm_loadu_pd(&rtu_Dk_0[i]);

      // Sum: '<S251>/Sum' incorporates:
      //   Product: '<S251>/C[k]*xhat[k|k-1]'

      tmp_0 = _mm_loadu_pd(&rtu_yk[i]);
      (void)_mm_storeu_pd(&rtu_yk_0[i], _mm_sub_pd(tmp_0, _mm_add_pd(tmp_1, tmp)));
    }

    for (int32_T i{0}; i <= 2; i += 2) {
      // Product: '<S251>/Product3'
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_set1_pd(0.0));
      tmp_1 = _mm_loadu_pd(&rtu_Lk[i]);
      tmp = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_1,
        _mm_set1_pd(rtu_yk_0[0])), tmp));
      tmp_1 = _mm_loadu_pd(&rtu_Lk[i + 4]);
      tmp = _mm_loadu_pd(&rty_Lykyhatkk1[i]);
      (void)_mm_storeu_pd(&rty_Lykyhatkk1[i], _mm_add_pd(_mm_mul_pd(tmp_1,
        _mm_set1_pd(rtu_yk_0[1])), tmp));
    }
  } else if (localDW->MeasurementUpdate_MODE) {
    MeasurementUpdate_Disable(rty_Lykyhatkk1, localDW, localP);
  } else {
    // no actions
  }

  // End of Outputs for SubSystem: '<S227>/MeasurementUpdate'
}

//
// Output and update for atomic system:
//    '<S208>/ReducedQRN'
//    '<S352>/ReducedQRN'
//
void SupervisoryController::ReducedQRN(const real_T rtu_G[16], const real_T
  rtu_H[8], const real_T rtu_Q[16], const real_T rtu_R[4], const real_T rtu_N[8],
  real_T rty_Qbar[16], real_T rty_Rbar[4], real_T rty_Nbar[8])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T rtu_Q_0[16];
  real_T rtb_Add_a[8];
  real_T rtb_Transpose2[8];
  real_T rtu_H_0[4];
  real_T rtu_N_0[4];
  int32_T i;
  int32_T i_0;
  int32_T rtu_H_tmp;
  int32_T rtu_Q_tmp;

  // Product: '<S228>/Product' incorporates:
  //   Math: '<S228>/Transpose1'

  for (i = 0; i < 4; i++) {
    i_0 = 0;
    for (rtu_H_tmp = 0; rtu_H_tmp < 4; rtu_H_tmp++) {
      rtu_Q_tmp = i_0 + i;
      rtu_Q_0[rtu_Q_tmp] = 0.0;
      rtu_Q_0[rtu_Q_tmp] += rtu_Q[i] * rtu_G[rtu_H_tmp];
      rtu_Q_0[rtu_Q_tmp] += rtu_Q[i + 4] * rtu_G[rtu_H_tmp + 4];
      rtu_Q_0[rtu_Q_tmp] += rtu_Q[i + 8] * rtu_G[rtu_H_tmp + 8];
      rtu_Q_0[rtu_Q_tmp] += rtu_Q[i + 12] * rtu_G[rtu_H_tmp + 12];
      i_0 += 4;
    }
  }

  i = 0;
  for (i_0 = 0; i_0 < 4; i_0++) {
    for (rtu_H_tmp = 0; rtu_H_tmp <= 2; rtu_H_tmp += 2) {
      rtu_Q_tmp = rtu_H_tmp + i;
      (void)_mm_storeu_pd(&rty_Qbar[rtu_Q_tmp], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&rtu_G[rtu_H_tmp]);
      tmp_0 = _mm_loadu_pd(&rty_Qbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Qbar[rtu_Q_tmp], _mm_add_pd(tmp_0, _mm_mul_pd
        (_mm_set1_pd(rtu_Q_0[i]), tmp)));
      tmp = _mm_loadu_pd(&rtu_G[rtu_H_tmp + 4]);
      tmp_0 = _mm_loadu_pd(&rty_Qbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Qbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtu_Q_0[i + 1]), tmp), tmp_0));
      tmp = _mm_loadu_pd(&rtu_G[rtu_H_tmp + 8]);
      tmp_0 = _mm_loadu_pd(&rty_Qbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Qbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtu_Q_0[i + 2]), tmp), tmp_0));
      tmp = _mm_loadu_pd(&rtu_G[rtu_H_tmp + 12]);
      tmp_0 = _mm_loadu_pd(&rty_Qbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Qbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtu_Q_0[i + 3]), tmp), tmp_0));
    }

    i += 4;
  }

  // End of Product: '<S228>/Product'

  // Math: '<S228>/Transpose2'
  i = 0;
  for (i_0 = 0; i_0 < 2; i_0++) {
    rtb_Transpose2[i] = rtu_H[i_0];
    rtb_Transpose2[i + 1] = rtu_H[i_0 + 2];
    rtb_Transpose2[i + 2] = rtu_H[i_0 + 4];
    rtb_Transpose2[i + 3] = rtu_H[i_0 + 6];
    i += 4;
  }

  // End of Math: '<S228>/Transpose2'

  // Sum: '<S228>/Add' incorporates:
  //   Math: '<S228>/Transpose2'
  //   Product: '<S228>/Product1'

  for (i = 0; i < 4; i++) {
    i_0 = 0;
    for (rtu_H_tmp = 0; rtu_H_tmp < 2; rtu_H_tmp++) {
      rtu_Q_tmp = i_0 + i;
      rtb_Add_a[rtu_Q_tmp] = (((rtb_Transpose2[i_0 + 1] * rtu_Q[i + 4] +
        rtb_Transpose2[i_0] * rtu_Q[i]) + rtb_Transpose2[i_0 + 2] * rtu_Q[i + 8])
        + rtb_Transpose2[i_0 + 3] * rtu_Q[i + 12]) + rtu_N[rtu_Q_tmp];
      i_0 += 4;
    }
  }

  // End of Sum: '<S228>/Add'
  for (i = 0; i < 2; i++) {
    for (i_0 = 0; i_0 <= 2; i_0 += 2) {
      // Product: '<S228>/Product2' incorporates:
      //   Sum: '<S228>/Add'

      rtu_H_tmp = i << 2UL;
      rtu_Q_tmp = i_0 + rtu_H_tmp;
      (void)_mm_storeu_pd(&rty_Nbar[rtu_Q_tmp], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&rtu_G[i_0]);
      tmp_0 = _mm_loadu_pd(&rty_Nbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Nbar[rtu_Q_tmp], _mm_add_pd(tmp_0, _mm_mul_pd
        (_mm_set1_pd(rtb_Add_a[rtu_H_tmp]), tmp)));
      tmp = _mm_loadu_pd(&rtu_G[i_0 + 4]);
      tmp_0 = _mm_loadu_pd(&rty_Nbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Nbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtb_Add_a[rtu_H_tmp + 1]), tmp), tmp_0));
      tmp = _mm_loadu_pd(&rtu_G[i_0 + 8]);
      tmp_0 = _mm_loadu_pd(&rty_Nbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Nbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtb_Add_a[rtu_H_tmp + 2]), tmp), tmp_0));
      tmp = _mm_loadu_pd(&rtu_G[i_0 + 12]);
      tmp_0 = _mm_loadu_pd(&rty_Nbar[rtu_Q_tmp]);
      (void)_mm_storeu_pd(&rty_Nbar[rtu_Q_tmp], _mm_add_pd(_mm_mul_pd
        (_mm_set1_pd(rtb_Add_a[rtu_H_tmp + 3]), tmp), tmp_0));
    }

    for (i_0 = 0; i_0 < 2; i_0++) {
      int32_T rtu_N_tmp;

      // Product: '<S228>/Product3' incorporates:
      //   Product: '<S228>/Product4'

      rtu_H_tmp = (i_0 << 1UL) + i;
      rtu_H_0[rtu_H_tmp] = 0.0;

      // Product: '<S228>/Product4'
      rtu_N_0[rtu_H_tmp] = 0.0;

      // Product: '<S228>/Product3' incorporates:
      //   Product: '<S228>/Product4'
      //   Sum: '<S228>/Add'

      rtu_Q_tmp = i_0 << 2UL;
      rtu_H_0[rtu_H_tmp] += rtb_Add_a[rtu_Q_tmp] * rtu_H[i];

      // Product: '<S228>/Product4' incorporates:
      //   Math: '<S228>/Transpose'
      //   Math: '<S228>/Transpose2'

      rtu_N_tmp = i << 2UL;
      rtu_N_0[rtu_H_tmp] += rtu_N[rtu_N_tmp] * rtb_Transpose2[rtu_Q_tmp];

      // Product: '<S228>/Product3' incorporates:
      //   Sum: '<S228>/Add'

      rtu_H_0[rtu_H_tmp] += rtb_Add_a[rtu_Q_tmp + 1] * rtu_H[i + 2];

      // Product: '<S228>/Product4' incorporates:
      //   Math: '<S228>/Transpose'
      //   Math: '<S228>/Transpose2'
      //   Product: '<S228>/Product3'

      rtu_N_0[rtu_H_tmp] += rtu_N[rtu_N_tmp + 1] * rtb_Transpose2[rtu_Q_tmp + 1];

      // Product: '<S228>/Product3' incorporates:
      //   Sum: '<S228>/Add'

      rtu_H_0[rtu_H_tmp] += rtb_Add_a[rtu_Q_tmp + 2] * rtu_H[i + 4];

      // Product: '<S228>/Product4' incorporates:
      //   Math: '<S228>/Transpose'
      //   Math: '<S228>/Transpose2'
      //   Product: '<S228>/Product3'

      rtu_N_0[rtu_H_tmp] += rtu_N[rtu_N_tmp + 2] * rtb_Transpose2[rtu_Q_tmp + 2];

      // Product: '<S228>/Product3' incorporates:
      //   Sum: '<S228>/Add'

      rtu_H_0[rtu_H_tmp] += rtb_Add_a[rtu_Q_tmp + 3] * rtu_H[i + 6];

      // Product: '<S228>/Product4' incorporates:
      //   Math: '<S228>/Transpose'
      //   Math: '<S228>/Transpose2'
      //   Product: '<S228>/Product3'

      rtu_N_0[rtu_H_tmp] += rtu_N[rtu_N_tmp + 3] * rtb_Transpose2[rtu_Q_tmp + 3];
    }
  }

  // Sum: '<S228>/Add1'
  rty_Rbar[0] = (rtu_H_0[0] + rtu_N_0[0]) + rtu_R[0];
  rty_Rbar[1] = (rtu_H_0[1] + rtu_N_0[1]) + rtu_R[1];
  rty_Rbar[2] = (rtu_H_0[2] + rtu_N_0[2]) + rtu_R[2];
  rty_Rbar[3] = (rtu_H_0[3] + rtu_N_0[3]) + rtu_R[3];
}

//
// Output and update for atomic system:
//    '<S208>/ScalarExpansionQ'
//    '<S352>/ScalarExpansionQ'
//
void SupervisoryController::ScalarExpansionQ(const real_T rtu_u[16], real_T
  rty_y[16])
{
  int32_T tmp;

  // MATLAB Function: '<S230>/ScalarExpansion'
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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S252>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S252>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  tmp = 0;
  for (int32_T i{0}; i < 4; i++) {
    rty_y[tmp] = (rtu_u[tmp] + rtu_u[i]) / 2.0;
    rty_y[tmp + 1] = (rtu_u[tmp + 1] + rtu_u[i + 4]) / 2.0;
    rty_y[tmp + 2] = (rtu_u[tmp + 2] + rtu_u[i + 8]) / 2.0;
    rty_y[tmp + 3] = (rtu_u[tmp + 3] + rtu_u[i + 12]) / 2.0;
    tmp += 4;
  }

  // End of MATLAB Function: '<S230>/ScalarExpansion'
}

//
// Output and update for atomic system:
//    '<S208>/ScalarExpansionR'
//    '<S352>/ScalarExpansionR'
//
void SupervisoryController::ScalarExpansionR(const real_T rtu_u[4], real_T
  rty_y[4])
{
  real_T tmp;

  // MATLAB Function: '<S231>/ScalarExpansion'
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
  // MATLAB Function 'Utilities/ScalarExpansion/ScalarExpansion': '<S253>:1'
  //    Copyright 2014-2015 The MathWorks, Inc.
  // '<S253>:1:16' y = ctrlScalarExpansion(u,n,IsStrictPositiveDefinite,OutputSquareRootY); 
  rty_y[0] = (rtu_u[0] + rtu_u[0]) / 2.0;
  tmp = (rtu_u[1] + rtu_u[2]) / 2.0;
  rty_y[1] = tmp;
  rty_y[2] = tmp;
  rty_y[3] = (rtu_u[3] + rtu_u[3]) / 2.0;
}

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
int32_T SupervisoryController::xpotrf_n(real_T b_A[16])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
real_T SupervisoryController::minimum_i(const real_T x[4])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::trisolve_i(const real_T b_A[16], real_T b_B[16])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
real_T SupervisoryController::norm_h(const real_T x[4])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
real_T SupervisoryController::xnrm2_n(int32_T n, const real_T x[16], int32_T ix0)
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::xgemv_o(int32_T b_m, int32_T n, const real_T b_A[16],
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::xgerc_j(int32_T b_m, int32_T n, real_T alpha1,
  int32_T ix0, const real_T y[4], real_T b_A[16], int32_T ia0)
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::KWIKfactor_h(const real_T b_Ac[344], const int32_T
  iC[86], int32_T nA, const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
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
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 85];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 171];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 257];
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
      beta1 = xnrm2_n(3 - k_i, b_A, ii + 2);
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

          beta1 = rt_hypotd_snf(atmp, xnrm2_n(3 - k_i, b_A, ii + 2));
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
        xgemv_o(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        xgerc_j(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
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
        xgemv_o(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 4, work);
        xgerc_j(knt, b_coltop, -tau[k_i], b_lastv - 4, work, b_A, b_lastv);
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
real_T SupervisoryController::maximum_l(const real_T x[4])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::DropConstraint_j(int32_T kDrop, boolean_T iA[86],
  int32_T *nA, int32_T iC[86])
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::qpkwik_i(const real_T b_Linv[16], const real_T
  b_Hinv[16], const real_T f[4], const real_T b_Ac[344], const real_T b[86],
  boolean_T iA[86], int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda
  [86], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[86];
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
  int32_T iC[86];
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
  for (i = 0; i < 86; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 86; tmp++) {
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
        KWIKfactor_h(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees_h, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            (void)std::memset(&iA[0], 0, 86U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 86U * sizeof(int32_T));
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
              (void)std::memset(&iA[0], 0, 86U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 86U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint_j(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 86U * sizeof(real_T));
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
    Xnorm0 = norm_h(x);
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (*status <= maxiter)) {
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 86; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 86] * x[1]);
          z[2] = std::abs(b_Ac[i + 172] * x[2]);
          z[3] = std::abs(b_Ac[i + 258] * x[3]);
          t = std::fmax(t, maximum_l(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 86] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 172] * x[2])
                   + b_Ac[i + 258] * x[3]) - b[i]) / t;
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
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 258])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 172])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 86])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0]), _mm_set1_pd
                  (b_Ac[tmp])), _mm_set1_pd(0.0))))));
              }

              guard2 = true;
            } else {
              KWIKfactor_h(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees_h, RLinv,
                           &cMin);
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
                    _mm_set1_pd(b_Ac[tmp + 258])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 172])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 86])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((b_D[iSave + 1] * b_Ac[tmp + 86] + b_D[iSave] *
                           b_Ac[tmp]) + b_D[iSave + 2] * b_Ac[tmp + 172]) +
                    b_D[iSave + 3] * b_Ac[tmp + 258];
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

              t = b_Ac[tmp + 86];
              cVal_tmp = b_Ac[tmp + 172];
              cVal_tmp_0 = b_Ac[tmp + 258];
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
                  if ((iC_0 <= 86) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                (void)std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint_j(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  (void)std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == static_cast<int32_T>(degrees_h)) {
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
            cMin = norm_h(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 86; tmp++) {
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

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void SupervisoryController::mpcblock_optimizer_d(const real_T rseq[40], const
  real_T vseq[21], const real_T x[4], const real_T old_u[3], const boolean_T iA
  [86], const real_T b_Mlim[86], real_T b_Mx[344], real_T b_Mu1[258], real_T
  b_Mv[1806], const real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16],
  real_T b_Ac[344], const real_T b_Wy[2], const real_T b_Wdu[3], const real_T
  b_Jm[180], const real_T b_Wu[3], const real_T b_I1[180], const real_T b_A[16],
  const real_T Bu[252], const real_T Bv[84], const real_T b_C[8], const real_T
  Dv[42], const int32_T b_Mrows[86], real_T u[3], real_T useq[63], real_T
  *status, boolean_T iAout[86], DW_State1controlLawAMPC1 *localDW)
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
  __m128d tmp_0;
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T b_Sx[160];
  real_T Sum_0[120];
  real_T WySuJm[120];
  real_T b_Su1[120];
  real_T b_Mlim_0[86];
  real_T b_Mlim_1[86];
  real_T b_Mu1_0[86];
  real_T b_Kv[63];
  real_T CA_0[42];
  real_T L[16];
  real_T b_Kx[12];
  real_T b_B[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T CA[8];
  real_T CA_1[8];
  real_T Sum[6];
  real_T varargin_1[4];
  real_T zopt[4];
  real_T b_C_0[2];
  real_T normH;
  real_T s;
  int32_T CA_tmp;
  int32_T Tries;
  int32_T a_tmp;
  int32_T b_SuJm_tmp;
  int32_T i;
  int32_T i1;
  int32_T kidx;
  int16_T ixw;
  int8_T a[3600];
  int8_T b[16];
  int8_T rows[2];
  boolean_T exitg1;
  boolean_T guard1{ false };

  (void)std::memset(&useq[0], 0, 63U * sizeof(real_T));
  (void)std::memset(&iAout[0], 0, 86U * sizeof(boolean_T));
  for (i = 0; i < 2; i++) {
    for (Tries = 0; Tries < 4; Tries++) {
      CA_tmp = (Tries << 1UL) + i;
      CA[CA_tmp] = 0.0;
      kidx = Tries << 2UL;
      CA[CA_tmp] += b_A[kidx] * b_C[i];
      CA[CA_tmp] += b_A[kidx + 1] * b_C[i + 2];
      CA[CA_tmp] += b_A[kidx + 2] * b_C[i + 4];
      CA[CA_tmp] += b_A[kidx + 3] * b_C[i + 6];
    }

    for (Tries = 0; Tries < 3; Tries++) {
      CA_tmp = (Tries << 1UL) + i;
      Sum[CA_tmp] = 0.0;
      kidx = Tries << 2UL;
      Sum[CA_tmp] += Bu[kidx] * b_C[i];
      Sum[CA_tmp] += Bu[kidx + 1] * b_C[i + 2];
      Sum[CA_tmp] += Bu[kidx + 2] * b_C[i + 4];
      Sum[CA_tmp] += Bu[kidx + 3] * b_C[i + 6];
    }

    b_C_0[i] = 0.0;
    b_C_0[i] += b_C[i] * Bv[0];
    b_C_0[i] += b_C[i + 2] * Bv[1];
    b_C_0[i] += b_C[i + 4] * Bv[2];
    b_C_0[i] += b_C[i + 6] * Bv[3];
    localDW->b_Hv[i] = b_C_0[i];
    localDW->b_Hv[i + 40] = Dv[i];
  }

  i = 0;
  for (Tries = 0; Tries < 19; Tries++) {
    localDW->b_Hv[i + 80] = 0.0;
    localDW->b_Hv[i + 81] = 0.0;
    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 21; Tries++) {
    (void)std::memset(&localDW->b_Hv[i + 2], 0, 38U * sizeof(real_T));
    i += 40;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 4; i1++) {
    b_Sx[i] = CA[Tries];
    b_Sx[i + 1] = CA[Tries + 1];
    i += 40;
    Tries += 2;
  }

  for (i = 0; i < 38; i++) {
    b_Sx[i + 2] = 0.0;
    b_Sx[i + 42] = 0.0;
    b_Sx[i + 82] = 0.0;
    b_Sx[i + 122] = 0.0;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    b_Su1[i] = Sum[Tries];
    b_Su1[i + 1] = Sum[Tries + 1];
    i += 40;
    Tries += 2;
  }

  for (i = 0; i < 38; i++) {
    b_Su1[i + 2] = 0.0;
    b_Su1[i + 42] = 0.0;
    b_Su1[i + 82] = 0.0;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    localDW->Su[i] = Sum[Tries];
    localDW->Su[i + 1] = Sum[Tries + 1];
    i += 40;
    Tries += 2;
  }

  i = 0;
  for (Tries = 0; Tries < 57; Tries++) {
    localDW->Su[i + 120] = 0.0;
    localDW->Su[i + 121] = 0.0;
    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 60; Tries++) {
    (void)std::memset(&localDW->Su[i + 2], 0, 38U * sizeof(real_T));
    i += 40;
  }

  for (kidx = 0; kidx < 19; kidx++) {
    int8_T kidx_0;
    kidx_0 = static_cast<int8_T>(((kidx + 1) << 1UL) + 1);
    for (i = 0; i < 2; i++) {
      rows[i] = static_cast<int8_T>(i + static_cast<int32_T>(kidx_0));
      Tries = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 3; a_tmp++) {
        CA_tmp = Tries + i;
        Sum[CA_tmp] += ((Bu[i1 + 1] * CA[i + 2] + Bu[i1] * CA[i]) + Bu[i1 + 2] *
                        CA[i + 4]) + Bu[i1 + 3] * CA[i + 6];
        Tries += 2;
        i1 += 4;
      }
    }

    for (i = 0; i < 3; i++) {
      CA_tmp = i << 1UL;
      normH = Sum[CA_tmp];
      b_Su1[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = normH;
      Sum_0[CA_tmp] = normH;
      normH = Sum[CA_tmp + 1];
      b_Su1[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = normH;
      Sum_0[CA_tmp + 1] = normH;
    }

    for (i = 0; i < 57; i++) {
      CA_tmp = (i + 3) << 1UL;
      Sum_0[CA_tmp] = localDW->Su[(40 * i + static_cast<int32_T>(rows[0])) - 3];
      Sum_0[CA_tmp + 1] = localDW->Su[(40 * i + static_cast<int32_T>(rows[1])) -
        3];
    }

    for (i = 0; i < 60; i++) {
      Tries = i << 1UL;
      localDW->Su[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = Sum_0[Tries];
      localDW->Su[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = Sum_0[Tries +
        1];
    }

    for (i = 0; i <= 0; i += 2) {
      (void)_mm_storeu_pd(&b_C_0[i], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&CA[i]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[0])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 2]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[1])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 4]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[2])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 6]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[3])),
        tmp_0));
      tmp = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&CA_0[i], tmp);
    }

    for (i = 0; i < 20; i++) {
      CA_tmp = (i + 1) << 1UL;
      CA_0[CA_tmp] = localDW->b_Hv[(40 * i + static_cast<int32_T>(rows[0])) - 3];
      CA_0[CA_tmp + 1] = localDW->b_Hv[(40 * i + static_cast<int32_T>(rows[1]))
        - 3];
    }

    for (i = 0; i < 21; i++) {
      Tries = i << 1UL;
      localDW->b_Hv[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = CA_0[Tries];
      localDW->b_Hv[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = CA_0[Tries +
        1];
    }

    for (i = 0; i < 2; i++) {
      Tries = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 4; a_tmp++) {
        CA_tmp = Tries + i;
        CA_1[CA_tmp] = 0.0;
        CA_1[CA_tmp] += b_A[i1] * CA[i];
        CA_1[CA_tmp] += b_A[i1 + 1] * CA[i + 2];
        CA_1[CA_tmp] += b_A[i1 + 2] * CA[i + 4];
        CA_1[CA_tmp] += b_A[i1 + 3] * CA[i + 6];
        Tries += 2;
        i1 += 4;
      }
    }

    (void)std::memcpy(&CA[0], &CA_1[0], sizeof(real_T) << 3UL);
    for (i = 0; i < 4; i++) {
      Tries = i << 1UL;
      b_Sx[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = CA[Tries];
      b_Sx[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = CA[Tries + 1];
    }
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    for (a_tmp = 0; a_tmp < 40; a_tmp++) {
      b_SuJm_tmp = a_tmp + i;
      Sum_0[b_SuJm_tmp] = 0.0;
      kidx = 0;
      for (CA_tmp = 0; CA_tmp < 60; CA_tmp++) {
        Sum_0[b_SuJm_tmp] += localDW->Su[kidx + a_tmp] * b_Jm[CA_tmp + Tries];
        kidx += 40;
      }
    }

    i += 40;
    Tries += 60;
  }

  if (b_Mrows[0] > 0) {
    kidx = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (kidx < 86)) {
      if (b_Mrows[kidx] <= 40) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = -Sum_0[Tries - 1];
        b_Ac[kidx + 86] = -Sum_0[Tries + 39];
        b_Ac[kidx + 172] = -Sum_0[Tries + 79];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = -b_Sx[Tries - 1];
        b_Mx[kidx + 86] = -b_Sx[Tries + 39];
        b_Mx[kidx + 172] = -b_Sx[Tries + 79];
        b_Mx[kidx + 258] = -b_Sx[Tries + 119];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = -b_Su1[Tries - 1];
        b_Mu1[kidx + 86] = -b_Su1[Tries + 39];
        b_Mu1[kidx + 172] = -b_Su1[Tries + 79];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 86 * i] = -localDW->b_Hv[(40 * i + Tries) - 1];
        }

        kidx++;
      } else if (b_Mrows[kidx] <= 80) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = Sum_0[Tries - 41];
        b_Ac[kidx + 86] = Sum_0[Tries - 1];
        b_Ac[kidx + 172] = Sum_0[Tries + 39];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = b_Sx[Tries - 41];
        b_Mx[kidx + 86] = b_Sx[Tries - 1];
        b_Mx[kidx + 172] = b_Sx[Tries + 39];
        b_Mx[kidx + 258] = b_Sx[Tries + 79];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = b_Su1[Tries - 41];
        b_Mu1[kidx + 86] = b_Su1[Tries - 1];
        b_Mu1[kidx + 172] = b_Su1[Tries + 39];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 86 * i] = localDW->b_Hv[(40 * i + Tries) - 41];
        }

        kidx++;
      } else {
        exitg1 = true;
      }
    }
  }

  (void)std::memset(&b_B[0], 0, 9U * sizeof(real_T));
  b_B[0] = 1.0;
  b_B[4] = 1.0;
  b_B[8] = 1.0;
  kidx = -1;
  for (Tries = 0; Tries < 20; Tries++) {
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 20; i1++) {
        a_tmp = static_cast<int32_T>(c_A[20 * Tries + i1]);
        a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i]) *
          a_tmp);
        a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 1]) *
          a_tmp);
        a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 2]) *
          a_tmp);
        kidx += 3;
      }
    }
  }

  i = 0;
  for (Tries = 0; Tries < 3; Tries++) {
    for (i1 = 0; i1 < 60; i1++) {
      CA_tmp = i1 + i;
      I2Jm[CA_tmp] = 0.0;
      a_tmp = 0;
      for (kidx = 0; kidx < 60; kidx++) {
        I2Jm[CA_tmp] += static_cast<real_T>(a[a_tmp + i1]) * b_Jm[kidx + i];
        a_tmp += 60;
      }
    }

    i += 60;
  }

  ixw = 1;
  for (kidx = 0; kidx < 40; kidx++) {
    normH = b_Wy[ixw - 1];
    WySuJm[kidx] = normH * Sum_0[kidx];
    WySuJm[kidx + 40] = Sum_0[kidx + 40] * normH;
    WySuJm[kidx + 80] = Sum_0[kidx + 80] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 2) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    normH = b_Wu[ixw - 1];
    WuI2Jm[kidx] = normH * I2Jm[kidx];
    WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * normH;
    WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    normH = b_Wdu[ixw - 1];
    WduJm[kidx] = normH * b_Jm[kidx];
    WduJm[kidx + 60] = b_Jm[kidx + 60] * normH;
    WduJm[kidx + 120] = b_Jm[kidx + 120] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  for (i = 0; i < 3; i++) {
    for (Tries = 0; Tries < 3; Tries++) {
      b_SuJm_tmp = 3 * Tries + i;
      b_B[b_SuJm_tmp] = 0.0;
      for (i1 = 0; i1 < 40; i1++) {
        b_B[b_SuJm_tmp] += Sum_0[40 * i + i1] * WySuJm[40 * Tries + i1];
      }

      b_Jm_0[b_SuJm_tmp] = 0.0;
      normH = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        a_tmp = 60 * i + i1;
        kidx = 60 * Tries + i1;
        normH += I2Jm[a_tmp] * WuI2Jm[kidx];
        b_Jm_0[b_SuJm_tmp] += b_Jm[a_tmp] * WduJm[kidx];
      }

      b_H[i + (Tries << 2UL)] = (b_B[b_SuJm_tmp] + b_Jm_0[b_SuJm_tmp]) + normH;
    }
  }

  for (i = 0; i < 3; i++) {
    for (Tries = 0; Tries < 3; Tries++) {
      kidx = 3 * Tries + i;
      b_Jm_0[kidx] = 0.0;
      for (i1 = 0; i1 < 40; i1++) {
        b_Jm_0[kidx] += b_Su1[40 * i + i1] * WySuJm[40 * Tries + i1];
      }

      b_I1_0[kidx] = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        b_I1_0[kidx] += b_I1[60 * i + i1] * WuI2Jm[60 * Tries + i1];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
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
  for (Tries = 0; Tries < 4; Tries++) {
    i1 = 0;
    a_tmp = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      b_SuJm_tmp = i1 + Tries;
      b_Kx[b_SuJm_tmp] = 0.0;
      for (CA_tmp = 0; CA_tmp < 40; CA_tmp++) {
        b_Kx[b_SuJm_tmp] += b_Sx[CA_tmp + i] * WySuJm[CA_tmp + a_tmp];
      }

      i1 += 4;
      a_tmp += 40;
    }

    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 21; Tries++) {
    i1 = 0;
    a_tmp = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      b_SuJm_tmp = i1 + Tries;
      b_Kv[b_SuJm_tmp] = 0.0;
      for (CA_tmp = 0; CA_tmp < 40; CA_tmp++) {
        b_Kv[b_SuJm_tmp] += localDW->b_Hv[CA_tmp + i] * WySuJm[CA_tmp + a_tmp];
      }

      i1 += 21;
      a_tmp += 40;
    }

    i += 40;
  }

  for (i = 0; i <= 118; i += 2) {
    tmp = _mm_loadu_pd(&WySuJm[i]);
    (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  kidx = 0;
  (void)std::memcpy(&L[0], &b_H[0], sizeof(real_T) << 4UL);
  Tries = xpotrf_n(L);
  guard1 = false;
  if (Tries == 0) {
    varargin_1[0] = L[0];
    varargin_1[1] = L[5];
    varargin_1[2] = L[10];
    varargin_1[3] = L[15];
    if (minimum_i(varargin_1) > 1.4901161193847656E-7) {
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    normH = 0.0;
    Tries = 0;
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (Tries < 4)) {
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
      kidx = 2;
    } else {
      Tries = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (Tries <= 4)) {
        boolean_T guard2{ false };

        normH = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
          1.4901161193847656E-7;
        for (i = 0; i < 16; i++) {
          b[i] = 0;
        }

        b[0] = 1;
        b[5] = 1;
        b[10] = 1;
        b[15] = 1;
        for (i = 0; i < 16; i++) {
          b_H[i] += normH * static_cast<real_T>(b[i]);
          L[i] = b_H[i];
        }

        kidx = xpotrf_n(L);
        guard2 = false;
        if (kidx == 0) {
          varargin_1[0] = L[0];
          varargin_1[1] = L[5];
          varargin_1[2] = L[10];
          varargin_1[3] = L[15];
          if (minimum_i(varargin_1) > 1.4901161193847656E-7) {
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
    for (i = 0; i < 16; i++) {
      b[i] = 0;
    }

    b[0] = 1;
    b[5] = 1;
    b[10] = 1;
    b[15] = 1;
    i = 0;
    for (kidx = 0; kidx < 4; kidx++) {
      b_H[i] = static_cast<real_T>(b[i]);
      b_H[i + 1] = static_cast<real_T>(b[i + 1]);
      b_H[i + 2] = static_cast<real_T>(b[i + 2]);
      b_H[i + 3] = static_cast<real_T>(b[i + 3]);
      i += 4;
    }

    trisolve_i(L, b_H);
    varargin_1[0] = 0.0;
    varargin_1[1] = 0.0;
    varargin_1[2] = 0.0;
    varargin_1[3] = 0.0;
    for (kidx = 0; kidx < 3; kidx++) {
      real_T WuI2Jm_0;
      normH = 0.0;
      for (i = 0; i < 40; i++) {
        normH += WySuJm[40 * kidx + i] * rseq[i];
      }

      s = 0.0;
      for (i = 0; i < 21; i++) {
        s += b_Kv[21 * kidx + i] * vseq[i];
      }

      WuI2Jm_0 = 0.0;
      for (i = 0; i < 60; i++) {
        WuI2Jm_0 += WuI2Jm[60 * kidx + i] * b_utarget[i];
      }

      i = kidx << 2UL;
      varargin_1[kidx] = ((((((b_Kx[i + 1] * x[1] + b_Kx[i] * x[0]) + b_Kx[i + 2]
        * x[2]) + b_Kx[i + 3] * x[3]) + normH) + ((b_B[3 * kidx + 1] * old_u[1]
        + b_B[3 * kidx] * old_u[0]) + b_B[3 * kidx + 2] * old_u[2])) + s) +
        WuI2Jm_0;
    }

    for (i = 0; i < 86; i++) {
      iAout[i] = iA[i];
      b_Mlim_0[i] = (((b_Mx[i + 86] * x[1] + b_Mx[i] * x[0]) + b_Mx[i + 172] *
                      x[2]) + b_Mx[i + 258] * x[3]) + b_Mlim[i];
      b_Mu1_0[i] = 0.0;
      b_Mu1_0[i] += b_Mu1[i] * old_u[0];
      b_Mu1_0[i] += b_Mu1[i + 86] * old_u[1];
      b_Mu1_0[i] += b_Mu1[i + 172] * old_u[2];
    }

    i = 0;
    for (Tries = 0; Tries < 4; Tries++) {
      i1 = 0;
      for (a_tmp = 0; a_tmp < 4; a_tmp++) {
        kidx = i1 + Tries;
        L[kidx] = 0.0;
        L[kidx] += b_H[i] * b_H[i1];
        L[kidx] += b_H[i + 1] * b_H[i1 + 1];
        L[kidx] += b_H[i + 2] * b_H[i1 + 2];
        L[kidx] += b_H[i + 3] * b_H[i1 + 3];
        i1 += 4;
      }

      i += 4;
    }

    for (i = 0; i < 86; i++) {
      normH = 0.0;
      Tries = 0;
      for (i1 = 0; i1 < 21; i1++) {
        normH += b_Mv[Tries + i] * vseq[i1];
        Tries += 86;
      }

      b_Mlim_1[i] = -((b_Mlim_0[i] + b_Mu1_0[i]) + normH);
    }

    qpkwik_i(b_H, L, varargin_1, b_Ac, b_Mlim_1, iAout, 360, 1.0E-6, zopt,
             b_Mlim_0, &kidx);
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

// System initialize for function-call system: '<S1>/State1.controlLaw.AMPC1'
void SupervisoryController::State1controlLawAMPC1_Init(real_T rty_u[3], real_T
  rty_yhat[2], DW_State1controlLawAMPC1 *localDW, P_State1controlLawAMPC1
  *localP)
{
  int32_T i;
  int32_T t_tmp;
  uint32_T tseed;

  // InitializeConditions for Memory: '<S113>/Memory'
  (void)std::memcpy(&localDW->Memory_PreviousInput[0],
                    &localP->Memory_InitialCondition[0], 86U * sizeof(boolean_T));

  // InitializeConditions for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[0] = localP->last_mv_InitialCondition[0];

  // InitializeConditions for UnitDelay: '<S111>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[0] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[1] = localP->last_mv_InitialCondition[1];

  // InitializeConditions for UnitDelay: '<S111>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[1] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[2] = localP->last_mv_InitialCondition[2];

  // InitializeConditions for UnitDelay: '<S111>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[2] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S111>/Unit Delay4'
  localDW->UnitDelay4_DSTATE[0] = localP->UnitDelay4_InitialCondition[0];
  localDW->UnitDelay4_DSTATE[1] = localP->UnitDelay4_InitialCondition[1];

  // InitializeConditions for Delay: '<S146>/delayTheta'
  localDW->icLoad = true;

  // InitializeConditions for Delay: '<S146>/delayPhi'
  localDW->icLoad_o = true;

  // InitializeConditions for Delay: '<S146>/delayBufferU'
  localDW->icLoad_a = true;

  // InitializeConditions for Delay: '<S146>/delayBufferY'
  localDW->icLoad_p = true;

  // InitializeConditions for Delay: '<S146>/delayL'
  localDW->icLoad_n = true;

  // InitializeConditions for Delay: '<S147>/delayTheta'
  localDW->icLoad_pm = true;

  // InitializeConditions for Delay: '<S147>/delayPhi'
  localDW->icLoad_d = true;

  // InitializeConditions for Delay: '<S147>/delayBufferU'
  localDW->icLoad_nq = true;

  // InitializeConditions for Delay: '<S147>/delayBufferY'
  localDW->icLoad_k = true;

  // InitializeConditions for Delay: '<S147>/delayL'
  localDW->icLoad_g = true;

  // InitializeConditions for Delay: '<S208>/MemoryX'
  localDW->icLoad_f = true;

  // InitializeConditions for Delay: '<S208>/MemoryP'
  localDW->icLoad_a3 = true;

  // InitializeConditions for RandomNumber: '<S3>/Measurement Noise'
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  t_tmp = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) & 32768U));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[0] = tseed;
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[1] = tseed;
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[2] = tseed;
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];

  // End of InitializeConditions for RandomNumber: '<S3>/Measurement Noise'

  // SystemInitialize for Enabled SubSystem: '<S227>/MeasurementUpdate'
  MeasurementUpdate_Init(localDW->Product3, &localP->MeasurementUpdate_n);

  // End of SystemInitialize for SubSystem: '<S227>/MeasurementUpdate'

  // SystemInitialize for Outport: '<S3>/u'
  rty_u[0] = localP->u_Y0;
  rty_u[1] = localP->u_Y0;
  rty_u[2] = localP->u_Y0;

  // SystemInitialize for Outport: '<S3>/yhat'
  rty_yhat[0] = localP->yhat_Y0;
  rty_yhat[1] = localP->yhat_Y0;
}

// Disable for function-call system: '<S1>/State1.controlLaw.AMPC1'
void SupervisoryController::State1controlLawAMPC1_Disable
  (DW_State1controlLawAMPC1 *localDW, P_State1controlLawAMPC1 *localP)
{
  // Disable for Enabled SubSystem: '<S227>/MeasurementUpdate'
  if (localDW->MeasurementUpdate_n.MeasurementUpdate_MODE) {
    MeasurementUpdate_Disable(localDW->Product3, &localDW->MeasurementUpdate_n,
      &localP->MeasurementUpdate_n);
  }

  // End of Disable for SubSystem: '<S227>/MeasurementUpdate'
}

// Output and update for function-call system: '<S1>/State1.controlLaw.AMPC1'
void SupervisoryController::State1controlLawAMPC1(const real_T rtu_r[2], const
  real_T rtu_y[2], const real_T rtu_y0[2], const real_T rtu_u0[3], real_T
  rtu_paramEst, real_T rty_u[3], real_T rty_yhat[2], DW_State1controlLawAMPC1
  *localDW, P_State1controlLawAMPC1 *localP, P *rtP, ZCE_State1controlLawAMPC1
  *localZCE)
{
  static const real_T f[344]{ -1.0000306206339091, -1.9691324621873285E-5,
    -1.000061241772034, -3.9382763869063235E-5, -1.0000918634143874,
    -5.9074317748655875E-5, -1.000122485560982, -7.87659862677373E-5,
    -1.0001531082118313, -9.8457769433393677E-5, -1.0001837313669477,
    -0.00011814966725271126, -1.0002143550263445, -0.00013784167973277642,
    -1.0002449791900341, -0.00015753380688067554, -1.0002756038580298,
    -0.00017722604870349512, -1.0003062290303446, -0.0001969184052083217,
    -1.0003368547069913, -0.00021661087640224196, -1.0003674808879828,
    -0.00023630346229234261, -1.0003981075733319, -0.00025599616288571049,
    -1.0004287347630518, -0.00027568897818943246, -1.0004593624571549,
    -0.00029538190821059548, -1.0004899906556548, -0.00031507495295628662,
    -1.0005206193585638, -0.00033476811243359292, -1.0005512485658954,
    -0.00035446138664960169, -1.0005818782776621, -0.00037415477561140016,
    -1.0006125084938771, -0.00039384827932607564, 1.0000306206339091,
    1.9691324621873285E-5, 1.000061241772034, 3.9382763869063235E-5,
    1.0000918634143874, 5.9074317748655875E-5, 1.000122485560982,
    7.87659862677373E-5, 1.0001531082118313, 9.8457769433393677E-5,
    1.0001837313669477, 0.00011814966725271126, 1.0002143550263445,
    0.00013784167973277642, 1.0002449791900341, 0.00015753380688067554,
    1.0002756038580298, 0.00017722604870349512, 1.0003062290303446,
    0.0001969184052083217, 1.0003368547069913, 0.00021661087640224196,
    1.0003674808879828, 0.00023630346229234261, 1.0003981075733319,
    0.00025599616288571049, 1.0004287347630518, 0.00027568897818943246,
    1.0004593624571549, 0.00029538190821059548, 1.0004899906556548,
    0.00031507495295628662, 1.0005206193585638, 0.00033476811243359292,
    1.0005512485658954, 0.00035446138664960169, 1.0005818782776621,
    0.00037415477561140016, 1.0006125084938771, 0.00039384827932607564, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 2.2010078492673557E-5, -0.99997520047355148,
    4.4020285108379236E-5, -0.99995040112871181, 6.6030619855037473E-5,
    -0.99992560196547409, 8.8041082740568784E-5, -0.99990080298383122,
    0.00011005167377289379, -0.9998760041837762, 0.00013206239295993321,
    -0.999851205565302, 0.00015407324030960778, -0.99982640712840154,
    0.00017608421582983839, -0.99980160887306779, 0.00019809531952854606,
    -0.99977681079929381, 0.00022010655141365177, -0.9997520129070725,
    0.00024211791149307668, -0.99972721519639685, 0.000264129399774742,
    -0.99970241766725987, 0.00028614101626656905, -0.99967762031965446,
    0.00030815276097647927, -0.99965282315357373, 0.00033016463391239413,
    -0.99962802616901059, 0.00035217663508223517, -0.99960322936595791,
    0.00037418876449392413, -0.99957843274440872, 0.00039620102215538267,
    -0.99955363630435612, 0.00041821340807453268, -0.9995288400457929,
    0.00044022592225929609, -0.99950404396871217, -2.2010078492673557E-5,
    0.99997520047355148, -4.4020285108379236E-5, 0.99995040112871181,
    -6.6030619855037473E-5, 0.99992560196547409, -8.8041082740568784E-5,
    0.99990080298383122, -0.00011005167377289379, 0.9998760041837762,
    -0.00013206239295993321, 0.999851205565302, -0.00015407324030960778,
    0.99982640712840154, -0.00017608421582983839, 0.99980160887306779,
    -0.00019809531952854606, 0.99977681079929381, -0.00022010655141365177,
    0.9997520129070725, -0.00024211791149307668, 0.99972721519639685,
    -0.000264129399774742, 0.99970241766725987, -0.00028614101626656905,
    0.99967762031965446, -0.00030815276097647927, 0.99965282315357373,
    -0.00033016463391239413, 0.99962802616901059, -0.00035217663508223517,
    0.99960322936595791, -0.00037418876449392413, 0.99957843274440872,
    -0.00039620102215538267, 0.99955363630435612, -0.00041821340807453268,
    0.9995288400457929, -0.00044022592225929609, 0.99950404396871217, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T k[344]{ 0.34290108520429963, 0.3552101510579076,
    0.68580485205389019, 0.71041824524886121, 1.0287113006761572,
    1.0656242826766753, 1.3716204311984876, 1.4208282634451641,
    1.7145322437482702, 1.7760301876581424, 2.0574467384528954,
    2.1312300554194241, 2.4003639154397551, 2.4864278668328237,
    2.7432837748362431, 2.8416236220021549, 3.0862063167697542,
    3.1968173210312321, 3.4291315413676853, 3.5520089640238695,
    3.7720594487574348, 3.9071985510838805, 4.1149900390664023,
    4.2623860823150794, 4.4579233124219906, 4.61757155782128, 4.8008592689516014,
    4.9727549777062956, 5.1437979087826413, 5.32793634207394, 5.486739232042515,
    5.6831156510280278, 5.8296832388586317, 6.038292904672371,
    6.1726299293584015, 6.3934681031107843, 6.5155793036692344,
    6.7486412464470806, 6.8585313619185442, 7.103812334785073,
    -0.34290108520429963, -0.3552101510579076, -0.68580485205389019,
    -0.71041824524886121, -1.0287113006761572, -1.0656242826766753,
    -1.3716204311984876, -1.4208282634451641, -1.7145322437482702,
    -1.7760301876581424, -2.0574467384528954, -2.1312300554194241,
    -2.4003639154397551, -2.4864278668328237, -2.7432837748362431,
    -2.8416236220021549, -3.0862063167697542, -3.1968173210312321,
    -3.4291315413676853, -3.5520089640238695, -3.7720594487574348,
    -3.9071985510838805, -4.1149900390664023, -4.2623860823150794,
    -4.4579233124219906, -4.61757155782128, -4.8008592689516014,
    -4.9727549777062956, -5.1437979087826413, -5.32793634207394,
    -5.486739232042515, -5.6831156510280278, -5.8296832388586317,
    -6.038292904672371, -6.1726299293584015, -6.3934681031107843,
    -6.5155793036692344, -6.7486412464470806, -6.8585313619185442,
    -7.103812334785073, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, -0.43337106536580761,
    0.0749598690059731, -0.86675705070095521, 0.14990934539236234,
    -1.3001579562335595, 0.22484842912310582, -1.733573782191743,
    0.29977712016213809, -2.1670045288036346, 0.37469541847338994,
    -2.60045019629737, 0.44960332402078873, -3.0339107849010896,
    0.52450083676825809, -3.4673862948429415, 0.59938795667971823,
    -3.9008767263510791, 0.67426468371908554, -4.334382079653663,
    0.74913101785027292, -4.7679023549788582, 0.82398695903718977,
    -5.2014375525548378, 0.89883250724374175, -5.63498767260978,
    0.97366766243383085, -6.06855271537187, 1.0484924245713558,
    -6.5021326810692974, 1.1233067936202112, -6.9357275699302594,
    1.1981107695442887, -7.3693373821829606, 1.2729043523074759,
    -7.8029621180556088, 1.3476875418736569, -8.2366017777764213,
    1.4224603382067122, -8.67025636157362, 1.4972227412705188,
    0.43337106536580761, -0.0749598690059731, 0.86675705070095521,
    -0.14990934539236234, 1.3001579562335595, -0.22484842912310582,
    1.733573782191743, -0.29977712016213809, 2.1670045288036346,
    -0.37469541847338994, 2.60045019629737, -0.44960332402078873,
    3.0339107849010896, -0.52450083676825809, 3.4673862948429415,
    -0.59938795667971823, 3.9008767263510791, -0.67426468371908554,
    4.334382079653663, -0.74913101785027292, 4.7679023549788582,
    -0.82398695903718977, 5.2014375525548378, -0.89883250724374175,
    5.63498767260978, -0.97366766243383085, 6.06855271537187,
    -1.0484924245713558, 6.5021326810692974, -1.1233067936202112,
    6.9357275699302594, -1.1981107695442887, 7.3693373821829606,
    -1.2729043523074759, 7.8029621180556088, -1.3476875418736569,
    8.2366017777764213, -1.4224603382067122, 8.67025636157362,
    -1.4972227412705188, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.079433657162736954,
    -0.43031643660853247, 0.1588792179329564, -0.860620637419285,
    0.23833668240584346, -1.2909126025013016, 0.31780605067658774,
    -1.7211923319236231, 0.39728732284038321, -2.1514598257552864,
    0.47678049899242836, -2.5817150840653249, 0.556285579227926,
    -3.0119581069227683, 0.63580256364208343, -3.4421888943966432,
    0.71533145233011242, -3.8724074465559721, 0.79487224538722923,
    -4.3026137634697745, 0.87442494290865436, -4.7328078452070654,
    0.953989544989613, -5.1629896918368567, 1.0335660517253344,
    -5.5931593034281564, 1.1131544632110526, -6.02331668004997,
    1.1927547795420059, -6.453461821771298, 1.2723670008134373,
    -6.883594728661139, 1.351991127120594, -7.3137154007884861,
    1.4316271585587275, -7.74382383822233, 1.5112750952230938,
    -8.173920041031657, 1.5909349372089536, -8.60400400928545,
    -0.079433657162736954, 0.43031643660853247, -0.1588792179329564,
    0.860620637419285, -0.23833668240584346, 1.2909126025013016,
    -0.31780605067658774, 1.7211923319236231, -0.39728732284038321,
    2.1514598257552864, -0.47678049899242836, 2.5817150840653249,
    -0.556285579227926, 3.0119581069227683, -0.63580256364208343,
    3.4421888943966432, -0.71533145233011242, 3.8724074465559721,
    -0.79487224538722923, 4.3026137634697745, -0.87442494290865436,
    4.7328078452070654, -0.953989544989613, 5.1629896918368567,
    -1.0335660517253344, 5.5931593034281564, -1.1131544632110526,
    6.02331668004997, -1.1927547795420059, 6.453461821771298,
    -1.2723670008134373, 6.883594728661139, -1.351991127120594,
    7.3137154007884861, -1.4316271585587275, 7.74382383822233,
    -1.5112750952230938, 8.173920041031657, -1.5909349372089536,
    8.60400400928545, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T g[258]{ 0.34290108520429963, 0.3552101510579076,
    0.68580485205389019, 0.71041824524886121, 1.0287113006761572,
    1.0656242826766753, 1.3716204311984876, 1.4208282634451641,
    1.7145322437482702, 1.7760301876581424, 2.0574467384528954,
    2.1312300554194241, 2.4003639154397551, 2.4864278668328237,
    2.7432837748362431, 2.8416236220021549, 3.0862063167697542,
    3.1968173210312321, 3.4291315413676853, 3.5520089640238695,
    3.7720594487574348, 3.9071985510838805, 4.1149900390664023,
    4.2623860823150794, 4.4579233124219906, 4.61757155782128, 4.8008592689516014,
    4.9727549777062956, 5.1437979087826413, 5.32793634207394, 5.486739232042515,
    5.6831156510280278, 5.8296832388586317, 6.038292904672371,
    6.1726299293584015, 6.3934681031107843, 6.5155793036692344,
    6.7486412464470806, 6.8585313619185442, 7.103812334785073,
    -0.34290108520429963, -0.3552101510579076, -0.68580485205389019,
    -0.71041824524886121, -1.0287113006761572, -1.0656242826766753,
    -1.3716204311984876, -1.4208282634451641, -1.7145322437482702,
    -1.7760301876581424, -2.0574467384528954, -2.1312300554194241,
    -2.4003639154397551, -2.4864278668328237, -2.7432837748362431,
    -2.8416236220021549, -3.0862063167697542, -3.1968173210312321,
    -3.4291315413676853, -3.5520089640238695, -3.7720594487574348,
    -3.9071985510838805, -4.1149900390664023, -4.2623860823150794,
    -4.4579233124219906, -4.61757155782128, -4.8008592689516014,
    -4.9727549777062956, -5.1437979087826413, -5.32793634207394,
    -5.486739232042515, -5.6831156510280278, -5.8296832388586317,
    -6.038292904672371, -6.1726299293584015, -6.3934681031107843,
    -6.5155793036692344, -6.7486412464470806, -6.8585313619185442,
    -7.103812334785073, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0, -0.43337106536580761,
    0.0749598690059731, -0.86675705070095521, 0.14990934539236234,
    -1.3001579562335595, 0.22484842912310582, -1.733573782191743,
    0.29977712016213809, -2.1670045288036346, 0.37469541847338994,
    -2.60045019629737, 0.44960332402078873, -3.0339107849010896,
    0.52450083676825809, -3.4673862948429415, 0.59938795667971823,
    -3.9008767263510791, 0.67426468371908554, -4.334382079653663,
    0.74913101785027292, -4.7679023549788582, 0.82398695903718977,
    -5.2014375525548378, 0.89883250724374175, -5.63498767260978,
    0.97366766243383085, -6.06855271537187, 1.0484924245713558,
    -6.5021326810692974, 1.1233067936202112, -6.9357275699302594,
    1.1981107695442887, -7.3693373821829606, 1.2729043523074759,
    -7.8029621180556088, 1.3476875418736569, -8.2366017777764213,
    1.4224603382067122, -8.67025636157362, 1.4972227412705188,
    0.43337106536580761, -0.0749598690059731, 0.86675705070095521,
    -0.14990934539236234, 1.3001579562335595, -0.22484842912310582,
    1.733573782191743, -0.29977712016213809, 2.1670045288036346,
    -0.37469541847338994, 2.60045019629737, -0.44960332402078873,
    3.0339107849010896, -0.52450083676825809, 3.4673862948429415,
    -0.59938795667971823, 3.9008767263510791, -0.67426468371908554,
    4.334382079653663, -0.74913101785027292, 4.7679023549788582,
    -0.82398695903718977, 5.2014375525548378, -0.89883250724374175,
    5.63498767260978, -0.97366766243383085, 6.06855271537187,
    -1.0484924245713558, 6.5021326810692974, -1.1233067936202112,
    6.9357275699302594, -1.1981107695442887, 7.3693373821829606,
    -1.2729043523074759, 7.8029621180556088, -1.3476875418736569,
    8.2366017777764213, -1.4224603382067122, 8.67025636157362,
    -1.4972227412705188, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.079433657162736954,
    -0.43031643660853247, 0.1588792179329564, -0.860620637419285,
    0.23833668240584346, -1.2909126025013016, 0.31780605067658774,
    -1.7211923319236231, 0.39728732284038321, -2.1514598257552864,
    0.47678049899242836, -2.5817150840653249, 0.556285579227926,
    -3.0119581069227683, 0.63580256364208343, -3.4421888943966432,
    0.71533145233011242, -3.8724074465559721, 0.79487224538722923,
    -4.3026137634697745, 0.87442494290865436, -4.7328078452070654,
    0.953989544989613, -5.1629896918368567, 1.0335660517253344,
    -5.5931593034281564, 1.1131544632110526, -6.02331668004997,
    1.1927547795420059, -6.453461821771298, 1.2723670008134373,
    -6.883594728661139, 1.351991127120594, -7.3137154007884861,
    1.4316271585587275, -7.74382383822233, 1.5112750952230938,
    -8.173920041031657, 1.5909349372089536, -8.60400400928545,
    -0.079433657162736954, 0.43031643660853247, -0.1588792179329564,
    0.860620637419285, -0.23833668240584346, 1.2909126025013016,
    -0.31780605067658774, 1.7211923319236231, -0.39728732284038321,
    2.1514598257552864, -0.47678049899242836, 2.5817150840653249,
    -0.556285579227926, 3.0119581069227683, -0.63580256364208343,
    3.4421888943966432, -0.71533145233011242, 3.8724074465559721,
    -0.79487224538722923, 4.3026137634697745, -0.87442494290865436,
    4.7328078452070654, -0.953989544989613, 5.1629896918368567,
    -1.0335660517253344, 5.5931593034281564, -1.1131544632110526,
    6.02331668004997, -1.1927547795420059, 6.453461821771298,
    -1.2723670008134373, 6.883594728661139, -1.351991127120594,
    7.3137154007884861, -1.4316271585587275, 7.74382383822233,
    -1.5112750952230938, 8.173920041031657, -1.5909349372089536,
    8.60400400928545, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

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

  static const real_T c[32]{ -0.34290108520429963, -0.3552101510579076, 0.0, 0.0,
    0.43337106536580761, -0.0749598690059731, 0.0, 0.0, -0.079433657162736954,
    0.43031643660853247, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0,
    0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b[16]{ 1.0000306206339091, 1.9691324621873285E-5, 0.0, 0.0,
    -2.2010078492673557E-5, 0.99997520047355148, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 1.0 };

  static const real_T h[16]{ 13.360197886127139, -6.41609373732472,
    -6.5999203145628318, 0.0, -6.41609373732472, 10.718738864841944,
    -3.5057363011318934, 0.0, -6.5999203145628318, -3.5057363011318934,
    10.608801223672017, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const int32_T b_Mrows[86]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 141, 142, 143 };

  static const int16_T e[86]{ 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 40, 40, 40, 40, 0, 0 };

  static const int8_T d[8]{ 1, 0, 0, 1, 1, 0, 0, 1 };

  static const int8_T old_mvoff_0[3]{ 40, 0, 0 };

  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  real_T f_0[344];
  real_T k_0[344];
  real_T g_0[258];
  real_T Bu[252];
  real_T b_Mlim[86];
  real_T Bv[84];
  real_T rtb_useq[63];
  real_T b_utarget[60];
  real_T Dv[42];
  real_T rseq[40];
  real_T b_B[32];
  real_T B_est[28];
  real_T vseq[21];
  real_T rtb_A_m[16];
  real_T rtb_Product_on[16];
  real_T rtb_Q[16];
  real_T rtb_Z_c[16];
  real_T D_est[14];
  real_T rtb_R_h_tmp[14];
  real_T rtb_B_k[12];
  real_T rtb_C_d1[8];
  real_T rtb_L_i[8];
  real_T rtb_N[8];
  real_T rtb_Product2_o[8];
  real_T rtb_B_fo[6];
  real_T rtb_D_oi[6];
  real_T rtb_A_kh[4];
  real_T rtb_C_j[4];
  real_T rtb_MemoryX_c[4];
  real_T tmp[4];
  real_T rtb_B_n[3];
  real_T rtb_D_c[3];
  real_T rtb_D_i[3];
  real_T tmp_0[3];
  real_T tmp_1[3];
  real_T Y[2];
  real_T rtb_Add1_gs[2];
  real_T rtb_C;
  real_T rtb_C_ob;
  real_T rtb_e_a;
  real_T rtb_e_h;
  int32_T b_utarget_tmp;
  int32_T i;
  int32_T i_0;
  int32_T old_mvoff;
  boolean_T tmp_6[86];
  boolean_T rtb_DataTypeConversionEnable_k;

  // Delay: '<S146>/delayTheta' incorporates:
  //   Constant: '<S146>/InitialParameters'

  if (localDW->icLoad) {
    localDW->delayTheta_DSTATE[0] = localP->InitialParameters_Value[0];
    localDW->delayTheta_DSTATE[1] = localP->InitialParameters_Value[1];
    localDW->delayTheta_DSTATE[2] = localP->InitialParameters_Value[2];
    localDW->delayTheta_DSTATE[3] = localP->InitialParameters_Value[3];
  }

  // End of Delay: '<S146>/delayTheta'

  // Delay: '<S146>/delayPhi' incorporates:
  //   Constant: '<S146>/InitialPhiMemory'

  if (localDW->icLoad_o) {
    for (i = 0; i < 7; i++) {
      localDW->delayPhi_DSTATE[i] = localP->InitialPhiMemory_Value[i];
    }
  }

  // End of Delay: '<S146>/delayPhi'

  // Delay: '<S146>/delayBufferU' incorporates:
  //   Constant: '<S146>/InitialInputs'

  if (localDW->icLoad_a) {
    localDW->delayBufferU_DSTATE = localP->InitialInputs_Value;
  }

  // End of Delay: '<S146>/delayBufferU'

  // Delay: '<S146>/delayBufferY' incorporates:
  //   Constant: '<S146>/InitialOutputs'

  if (localDW->icLoad_p) {
    localDW->delayBufferY_DSTATE = localP->InitialOutputs_Value;
  }

  // End of Delay: '<S146>/delayBufferY'

  // Delay: '<S146>/delayL' incorporates:
  //   Constant: '<S146>/InitialCovariance'

  if (localDW->icLoad_n) {
    (void)std::memcpy(&localDW->delayL_DSTATE[0],
                      &localP->InitialCovariance_Value[0], sizeof(real_T) << 4UL);
  }

  // End of Delay: '<S146>/delayL'

  // Update for Delay: '<S146>/delayL' incorporates:
  //   Constant: '<S146>/Forgetting Factor'
  //   DataTypeConversion: '<S146>/DataTypeConversionEnable'
  //   Delay: '<S146>/delayBufferU'
  //   Delay: '<S146>/delayBufferY'
  //   Delay: '<S146>/delayPhi'
  //   Delay: '<S146>/delayTheta'
  //   MATLAB Function: '<S169>/rarx'
  //   UnitDelay: '<S111>/Unit Delay4'
  //   UnitDelay: '<S111>/Unit Delay5'

  rarx(localDW->UnitDelay5_DSTATE, localDW->UnitDelay4_DSTATE[0], rtu_paramEst
       != 0.0, localP->ForgettingFactor_Value, rtb_Add1_gs, rtb_B_fo, &rtb_e_a,
       localDW->delayTheta_DSTATE, localDW->delayL_DSTATE,
       localDW->delayPhi_DSTATE, &localDW->sf_rarx);

  // MATLAB Function: '<S144>/Model Type Converter'
  ModelTypeConverter(rtb_B_fo, rtb_Add1_gs, &rtb_e_a, rtb_B_n, &rtb_C, rtb_D_c);

  // Delay: '<S147>/delayTheta' incorporates:
  //   Constant: '<S147>/InitialParameters'

  if (localDW->icLoad_pm) {
    localDW->delayTheta_DSTATE_l[0] = localP->InitialParameters_Value_n[0];
    localDW->delayTheta_DSTATE_l[1] = localP->InitialParameters_Value_n[1];
    localDW->delayTheta_DSTATE_l[2] = localP->InitialParameters_Value_n[2];
    localDW->delayTheta_DSTATE_l[3] = localP->InitialParameters_Value_n[3];
  }

  // End of Delay: '<S147>/delayTheta'

  // Delay: '<S147>/delayPhi' incorporates:
  //   Constant: '<S147>/InitialPhiMemory'

  if (localDW->icLoad_d) {
    for (i = 0; i < 7; i++) {
      localDW->delayPhi_DSTATE_j[i] = localP->InitialPhiMemory_Value_o[i];
    }
  }

  // End of Delay: '<S147>/delayPhi'

  // Delay: '<S147>/delayBufferU' incorporates:
  //   Constant: '<S147>/InitialInputs'

  if (localDW->icLoad_nq) {
    localDW->delayBufferU_DSTATE_c = localP->InitialInputs_Value_n;
  }

  // End of Delay: '<S147>/delayBufferU'

  // Delay: '<S147>/delayBufferY' incorporates:
  //   Constant: '<S147>/InitialOutputs'

  if (localDW->icLoad_k) {
    localDW->delayBufferY_DSTATE_j = localP->InitialOutputs_Value_n;
  }

  // End of Delay: '<S147>/delayBufferY'

  // Delay: '<S147>/delayL' incorporates:
  //   Constant: '<S147>/InitialCovariance'

  if (localDW->icLoad_g) {
    (void)std::memcpy(&localDW->delayL_DSTATE_j[0],
                      &localP->InitialCovariance_Value_a[0], sizeof(real_T) <<
                      4UL);
  }

  // End of Delay: '<S147>/delayL'

  // Update for Delay: '<S147>/delayL' incorporates:
  //   Constant: '<S147>/Forgetting Factor'
  //   DataTypeConversion: '<S147>/DataTypeConversionEnable'
  //   Delay: '<S147>/delayBufferU'
  //   Delay: '<S147>/delayBufferY'
  //   Delay: '<S147>/delayPhi'
  //   Delay: '<S147>/delayTheta'
  //   MATLAB Function: '<S198>/rarx'
  //   UnitDelay: '<S111>/Unit Delay4'
  //   UnitDelay: '<S111>/Unit Delay5'

  rarx(localDW->UnitDelay5_DSTATE, localDW->UnitDelay4_DSTATE[1], rtu_paramEst
       != 0.0, localP->ForgettingFactor_Value_o, rtb_Add1_gs, rtb_B_fo, &rtb_e_h,
       localDW->delayTheta_DSTATE_l, localDW->delayL_DSTATE_j,
       localDW->delayPhi_DSTATE_j, &localDW->sf_rarx_l);

  // MATLAB Function: '<S145>/Model Type Converter'
  ModelTypeConverter(rtb_B_fo, rtb_Add1_gs, &rtb_e_h, rty_u, &rtb_C_ob, rtb_D_i);

  // MATLAB Function: '<S111>/MATLAB Function'
  MATLABFunction(rtb_e_a, rtb_B_n, rtb_C, rtb_D_c, rtb_e_h, rty_u, rtb_C_ob,
                 rtb_D_i, rtb_A_kh, rtb_B_fo, rtb_C_j, rtb_D_oi);

  // Product: '<S111>/Product1' incorporates:
  //   Constant: '<S111>/Constant1'
  //   Product: '<S211>/Product1'

  rtb_Add1_gs[0] = localP->Constant1_Value * rtu_y0[0];
  rtb_Add1_gs[1] = localP->Constant1_Value * rtu_y0[1];

  // Delay: '<S208>/MemoryX' incorporates:
  //   Constant: '<S208>/X0'
  //   DataTypeConversion: '<S208>/DataTypeConversionReset'

  localDW->icLoad_f = ((static_cast<uint32_T>(localZCE->MemoryX_Reset_ZCE_b) ==
                        POS_ZCSIG) || localDW->icLoad_f);
  localZCE->MemoryX_Reset_ZCE_b = 0U;
  if (localDW->icLoad_f) {
    localDW->MemoryX_DSTATE[0] = localP->X0_Value[0];
    localDW->MemoryX_DSTATE[1] = localP->X0_Value[1];
    localDW->MemoryX_DSTATE[2] = localP->X0_Value[2];
    localDW->MemoryX_DSTATE[3] = localP->X0_Value[3];
  }

  rtb_MemoryX_c[0] = localDW->MemoryX_DSTATE[0];
  rtb_MemoryX_c[1] = localDW->MemoryX_DSTATE[1];
  rtb_MemoryX_c[2] = localDW->MemoryX_DSTATE[2];
  rtb_MemoryX_c[3] = localDW->MemoryX_DSTATE[3];

  // MATLAB Function: '<S141>/FixedHorizonOptimizer' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'
  //   Constant: '<S111>/Constant9'
  //   Constant: '<S112>/Constant'
  //   Delay: '<S208>/MemoryX'
  //   Product: '<S211>/Product1'
  //   RandomNumber: '<S3>/Measurement Noise'
  //   Sum: '<S112>/Sum2'
  //   UnitDelay: '<S113>/last_mv'
  //
  // MATLAB Function 'Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer': '<S142>:1' 
  // '<S142>:1:18' coder.extrinsic('mpcblock_optimizer_double_mex');
  // '<S142>:1:19' coder.extrinsic('mpcblock_optimizer_single_mex');
  // '<S142>:1:20' coder.extrinsic('mpcblock_refmd_double_mex');
  // '<S142>:1:21' coder.extrinsic('mpcblock_refmd_single_mex');
  //  Inputs (in BlockDataType except iA)
  //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
  // '<S142>:1:25' xk = convertDataType(xk0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    old_u:      last mv (calculated by MPC)
  // '<S142>:1:27' old_u = convertDataType(old_u0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    ym:         current measured output (used only with built-in KF)
  // '<S142>:1:29' ym = convertDataType(ym0,isDouble);
  //    ref:        output reference
  // '<S142>:1:31' ref = convertDataType(ref0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    md:         measured disturbance
  // '<S142>:1:33' md = convertDataType(md0,isDouble);
  //    umin:       run-time MV bound
  // '<S142>:1:35' umin = convertDataType(umin0,isDouble);
  //    umax:       run-time MV bound
  // '<S142>:1:37' umax = convertDataType(umax0,isDouble);
  //    ymin:       run-time OV bound
  // '<S142>:1:39' ymin = convertDataType(ymin0,isDouble);
  //    ymax:       run-time OV bound
  // '<S142>:1:41' ymax = convertDataType(ymax0,isDouble);
  //    E:          run-time mixed constraints
  // '<S142>:1:43' E = convertDataType(E0,isDouble);
  //    F:          run-time mixed constraints
  // '<S142>:1:45' F = convertDataType(F0,isDouble);
  //    G:          run-time mixed constraints
  // '<S142>:1:47' G = convertDataType(G0,isDouble);
  //    S:          run-time mixed constraints
  // '<S142>:1:49' S = convertDataType(S0,isDouble);
  //    switch_in:  if it matches "enable_value", MPC is active in control
  // '<S142>:1:51' switch_in = int32(switch_in0);
  //    ext_mv:     external last mv (actual)
  // '<S142>:1:53' ext_mv = convertDataType(ext_mv0,isDouble);
  //    MVtarget:   MV reference
  // '<S142>:1:55' MVtarget = convertDataType(MVtarget0,isDouble);
  //    ywt:        run-time OV weights
  // '<S142>:1:57' ywt = convertDataType(ywt0,isDouble);
  //    uwt:        run-time MV weights
  // '<S142>:1:59' uwt = convertDataType(uwt0,isDouble);
  //    duwt:       run-time DMV weights
  // '<S142>:1:61' duwt = convertDataType(duwt0,isDouble);
  //    rhoeps:     run-time Slack weights
  // '<S142>:1:63' ewt = convertDataType(ewt0,isDouble);
  //    a:          run-time A (must be in DT)
  // '<S142>:1:65' a = convertDataType(a0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    b:          run-time B (must be in DT)
  // '<S142>:1:67' b = convertDataType(b0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    c:          run-time C (must be in DT)
  // '<S142>:1:69' c = convertDataType(c0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    d:          run-time D (must be in DT)
  // '<S142>:1:71' d = convertDataType(d0,isDouble);
  //    U:          run-time nominal value
  // '<S142>:1:73' U = convertDataType(U0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    Y:          run-time nominal value
  // '<S142>:1:75' Y = convertDataType(Y0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    X:          run-time nominal value
  // '<S142>:1:77' X = convertDataType(X0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    DX:         run-time nominal value
  // '<S142>:1:79' DX = convertDataType(DX0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
  //    Pk:         covariance P[k|k-1] (used only with built-in KF)
  // '<S142>:1:81' Pk = convertDataType(Pk0,isDouble);
  // '<S142>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S142>:1:319' if isa(u,'double')
  // '<S142>:1:320' y = u;
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
  // '<S142>:1:115' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
  // '<S142>:1:116' isAdaptive = ~isLTV;
  //  isLTV=true forces isAdaptive=false
  // '<S142>:1:117' ZERO = zeros('like',ref);
  // '<S142>:1:118' ONE = ones('like',ref);
  // '<S142>:1:119' hasMD = nv>int32(1);
  //  Pre-allocate all the MEX block outputs for the simulation mode
  // '<S142>:1:123' if isSimulation
  //  Model update
  // '<S142>:1:137' nym = int32(length(myindex));
  // '<S142>:1:138' ai=zeros(nxp,nxp,'like',ref);
  // '<S142>:1:139' bi=zeros(nxp,nup,'like',ref);
  // '<S142>:1:140' ci=zeros(ny,nxp,'like',ref);
  // '<S142>:1:141' di=zeros(ny,nup,'like',ref);
  // '<S142>:1:143' ai(:,:)=a(:,:,1);
  // '<S142>:1:144' bi(:,:)=b(:,:,1);
  // '<S142>:1:145' ci(:,:)=c(:,:,1);
  // '<S142>:1:146' di(:,:)=d(:,:,1);
  //  Allocate matrices. Must allocate 3D matrix also in Adaptive case,
  //  otherwise EML code does not compile.
  // '<S142>:1:150' Bu=zeros(nx,nu,p+1,'like',ref);
  (void)std::memset(&Bu[0], 0, 252U * sizeof(real_T));

  // '<S142>:1:151' Bv=zeros(nx,nv,p+1,'like',ref);
  (void)std::memset(&Bv[0], 0, 84U * sizeof(real_T));

  // '<S142>:1:152' Dv=zeros(ny,nv,p+1,'like',ref);
  (void)std::memset(&Dv[0], 0, 42U * sizeof(real_T));

  // '<S142>:1:153' Dvm=zeros(nym,nv,p+1,'like',ref);
  // '<S142>:1:154' Cm=zeros(nym,nx,p+1,'like',ref);
  // '<S142>:1:155' [A(:,:,1),C(:,:,1),Bu(:,:,1),Bv(:,:,1),Cm(:,:,1),Dv(:,:,1),Dvm(:,:,1),Qk,Rk,Nk] = mpc_plantupdate(... 
  // '<S142>:1:156'     ai,bi,ci,di,A(:,:,1),B(:,:,1),C(:,:,1),D(:,:,1),mvindex,mdindex,unindex,nxp,nup,ny,nu,nv,nxid, ... 
  // '<S142>:1:157'     myindex,Uscale,Yscale,Cid,Did);
  (void)std::memcpy(&rtb_A_m[0], &b[0], sizeof(real_T) << 4UL);
  (void)std::memcpy(&b_B[0], &c[0], sizeof(real_T) << 5UL);
  for (i = 0; i < 8; i++) {
    rtb_C_d1[i] = static_cast<real_T>(d[i]);
  }

  rtb_C_d1[0] = rtb_C_j[0];
  rtb_A_m[0] = rtb_A_kh[0];
  rtb_C_d1[1] = rtb_C_j[1];
  rtb_A_m[1] = rtb_A_kh[1];
  rtb_C_d1[2] = rtb_C_j[2];
  rtb_A_m[4] = rtb_A_kh[2];
  rtb_C_d1[3] = rtb_C_j[3];
  rtb_A_m[5] = rtb_A_kh[3];
  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    b_B[i] = rtb_B_fo[i_0];
    b_B[i + 1] = rtb_B_fo[i_0 + 1];
    Bu[i] = b_B[i];
    Bu[i + 1] = b_B[i + 1];
    Bu[i + 2] = b_B[i + 2];
    Bu[i + 3] = b_B[i + 3];
    i += 4;
    i_0 += 2;
  }

  Bv[2] = b_B[14];
  Bv[3] = b_B[15];
  Dv[0] = 0.0;
  Dv[1] = 0.0;

  // '<S142>:1:158' if isLTV
  //  Offset update together with Mlim, utarget, Bv and Dv values
  // '<S142>:1:174' [Mlim, utarget, uoff, voff, yoff, myoff, xoff, Bv, Dv] = ... 
  // '<S142>:1:175'     mpc_updateFromNominal(isAdaptive,isQP,Mlim,Mrows,...
  // '<S142>:1:176'        U,Uscale,uoff,mvindex,voff,mdindex,utarget,nu,nv-1,... 
  // '<S142>:1:177'        Y,Yscale,yoff,myoff,myindex,ny,...
  // '<S142>:1:178'        X,xoff,nxp,DX,A,Bu,Bv,C,Dv,nCC);
  for (i = 0; i < 86; i++) {
    b_Mlim[i] = static_cast<real_T>(e[i]);
  }

  (void)std::memset(&b_utarget[0], 0, 60U * sizeof(real_T));
  rtb_D_c[0] = rtu_u0[0];
  rtb_D_c[1] = rtu_u0[1];
  rtb_D_c[2] = rtu_u0[2];
  Y[0] = rtu_y0[0];
  Y[1] = rtu_y0[1];
  rtb_B_n[0] = rtu_u0[0];
  rtb_B_n[1] = rtu_u0[1];
  rtb_B_n[2] = rtu_u0[2];
  for (int32_T b_ct{0}; b_ct < 86; b_ct++) {
    rtb_e_a = b_Mlim[b_ct];
    i = b_Mrows[b_ct];
    if (i <= 40) {
      rtb_e_a += 100.0 - Y[(i - (((i - 1) / static_cast<int32_T>(ny)) << 1UL)) -
        1];
    } else if (i <= 80) {
      rtb_e_a -= 100.0 - Y[(i - (((i - 41) >> 1UL) << 1UL)) - 41];
    } else if (i <= 140) {
      i = (i - div_nde_s32_floor(i - 81, static_cast<int32_T>(nu_g)) *
           static_cast<int32_T>(nu_g)) - 81;
      rtb_e_a += static_cast<real_T>(old_mvoff_0[i]) - rtb_D_c[i];
    } else {
      i = (i - div_nde_s32_floor(i - 141, static_cast<int32_T>(nu_g)) *
           static_cast<int32_T>(nu_g)) - 141;
      rtb_e_a -= static_cast<real_T>(old_mvoff_0[i]) - rtb_D_c[i];
    }

    b_Mlim[b_ct] = rtb_e_a;
  }

  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    old_mvoff = static_cast<int32_T>(old_mvoff_0[b_ct]);
    rtb_e_a = rtb_D_c[b_ct];
    i = 0;
    for (i_0 = 0; i_0 < 20; i_0++) {
      b_utarget_tmp = i + b_ct;
      b_utarget[b_utarget_tmp] = (b_utarget[b_utarget_tmp] + static_cast<real_T>
        (old_mvoff)) - rtb_e_a;
      i += 3;
    }
  }

  Bv[0] = localP->Constant9_Value[0];
  Bv[1] = localP->Constant9_Value[1];

  //  Remove last u offset
  // '<S142>:1:181' old_u = old_u - uoff;
  //  Get reference and MD signals -- accounting for previewing
  // '<S142>:1:184' if isSimulation
  // '<S142>:1:190' else
  //  When doing code generation, use M code directly
  // '<S142>:1:192' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
  for (int32_T b_ct{0}; b_ct < 21; b_ct++) {
    vseq[b_ct] = 1.0;
  }

  i = 0;
  for (int32_T b_ct{0}; b_ct < 20; b_ct++) {
    rseq[i] = rtu_r[0] - Y[0];
    rseq[i + 1] = rtu_r[1] - Y[1];
    i += 2;
  }

  //  External MV override.
  //  NOTE: old_u and ext_mv input signals are dimensionless and offset-free.
  // '<S142>:1:197' if no_mv
  //  no external mv: old_u is the optimal u[k-1] from last_mv
  // '<S142>:1:199' delmv = zeros(nu,1,'like',ref);
  //  Obtain x[k|k]
  // '<S142>:1:208' xk = xk - xoff;
  //  Remove offset
  // '<S142>:1:209' if CustomEstimation
  //  Input is x(k|k)
  // '<S142>:1:211' xest = xk;
  //  Real-time MV target override
  //  Note: utargetValue is a vector length p*nu.
  // '<S142>:1:231' if no_uref
  //  no external utarget
  // '<S142>:1:233' utargetValue = utarget;
  //  Real-time custom constraint override (scaled E/F/S)
  // '<S142>:1:241' if ~no_cc
  // '<S142>:1:250' return_sequence = return_mvseq || return_xseq || return_ovseq; 
  // '<S142>:1:251' if isSimulation
  // '<S142>:1:279' else
  //  When doing code generation, use M code directly
  // '<S142>:1:281' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
  // '<S142>:1:282'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
  // '<S142>:1:283'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
  // '<S142>:1:284'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
  // '<S142>:1:285'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
  // '<S142>:1:286'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
  // '<S142>:1:287'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
  // '<S142>:1:288'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
  // '<S142>:1:289'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
  // '<S142>:1:290'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
  // '<S142>:1:291'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
  // '<S142>:1:292'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, isHyb, Mdis, Ndis, Vdis, numdis, maxdis); 
  tmp[0] = (localDW->MemoryX_DSTATE[0] + rtb_Add1_gs[0]) - rtb_Add1_gs[0];
  tmp[1] = (localDW->MemoryX_DSTATE[1] + rtb_Add1_gs[1]) - rtb_Add1_gs[1];
  tmp[2] = localP->Constant_Value[0] + localDW->MemoryX_DSTATE[2];
  tmp[3] = localP->Constant_Value[1] + localDW->MemoryX_DSTATE[3];
  rtb_D_i[0] = localDW->last_mv_DSTATE[0] - rtu_u0[0];
  rtb_D_i[1] = localDW->last_mv_DSTATE[1] - rtu_u0[1];
  rtb_D_i[2] = localDW->last_mv_DSTATE[2] - rtu_u0[2];
  (void)std::memset(&localDW->dv[0], 0, 1806U * sizeof(real_T));
  Y[0] = 0.018316915599999997;
  Y[1] = 0.018316915599999997;
  tmp_0[0] = 0.54594344475769718;
  tmp_1[0] = 0.0;
  tmp_0[1] = 0.54594344475769718;
  tmp_1[1] = 0.0;
  tmp_0[2] = 0.54594344475769718;
  tmp_1[2] = 0.0;

  // Memory: '<S113>/Memory'
  (void)std::memcpy(&tmp_6[0], &localDW->Memory_PreviousInput[0], 86U * sizeof
                    (boolean_T));

  // MATLAB Function: '<S141>/FixedHorizonOptimizer'
  (void)std::memcpy(&f_0[0], &f[0], 344U * sizeof(real_T));
  (void)std::memcpy(&g_0[0], &g[0], 258U * sizeof(real_T));
  (void)std::memcpy(&rtb_Q[0], &h[0], sizeof(real_T) << 4UL);
  (void)std::memcpy(&k_0[0], &k[0], 344U * sizeof(real_T));

  // Update for Memory: '<S113>/Memory' incorporates:
  //   MATLAB Function: '<S141>/FixedHorizonOptimizer'

  mpcblock_optimizer_d(rseq, vseq, tmp, rtb_D_i, tmp_6, b_Mlim, f_0, g_0,
                       localDW->dv, b_utarget, rtb_B_n, rtb_Q, k_0, Y, tmp_0, l,
                       tmp_1, n, rtb_A_m, Bu, Bv, rtb_C_d1, Dv, b_Mrows, rtb_D_c,
                       rtb_useq, &rtb_e_a, localDW->Memory_PreviousInput,
                       localDW);

  // Delay: '<S208>/MemoryP' incorporates:
  //   Constant: '<S208>/P0'
  //   DataTypeConversion: '<S208>/DataTypeConversionReset'

  // '<S142>:1:295' if return_xseq || return_ovseq
  // '<S142>:1:297' else
  // '<S142>:1:298' yseq = zeros(p+1,ny,'like',rseq);
  // '<S142>:1:299' xseq = zeros(p+1,nxQP,'like',rseq);
  // '<S142>:1:302' if CustomEstimation
  // '<S142>:1:303' xk1 = zeros(nx,1,'like',ref);
  // '<S142>:1:304' Pk1 = Pk;
  // '<S142>:1:311' xk1 = xk1 + xoff;
  //  Updated state must include offset
  //  return xest in original value
  // '<S142>:1:314' xest = xest + xoff;
  localDW->icLoad_a3 = ((static_cast<uint32_T>(localZCE->MemoryP_Reset_ZCE_k) ==
    POS_ZCSIG) || localDW->icLoad_a3);
  localZCE->MemoryP_Reset_ZCE_k = 0U;
  if (localDW->icLoad_a3) {
    (void)std::memcpy(&localDW->MemoryP_DSTATE[0], &localP->P0_Value[0], sizeof
                      (real_T) << 4UL);
  }

  // DataTypeConversion: '<S208>/DataTypeConversionEnable' incorporates:
  //   Constant: '<S3>/Constant'

  rtb_DataTypeConversionEnable_k = (localP->Constant_Value_f != 0.0);

  // MATLAB Function: '<S112>/MATLAB Function' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'

  // MATLAB Function 'SupervisoryController/State1.controlLaw.AMPC1/State Estimator (KF)/MATLAB Function': '<S209>:1' 
  // '<S209>:1:2' [A, B, C, D, Q, R, N] = stateEstimator(Ap, Bp, Cp, Dp, Aod1, Bod1, Cod1, Dod1, Dmn1); 
  // 'stateEstimator:3' no = size(Cp,1);
  //  n_outputs
  // 'stateEstimator:4' ni = 3;
  //  n_inputs
  // 'stateEstimator:5' nsp = size(Ap,1);
  //  n_plant_states
  // 'stateEstimator:6' ns = nsp + no;
  //  n_states = n_plant_states + n_outputs
  // 'stateEstimator:8' A = zeros(ns);
  //  n_states x n_states
  // 'stateEstimator:9' B = zeros(ns,ni);
  //  n_states  x n_inputs
  // 'stateEstimator:10' C = zeros(no,ns);
  //  n_outputs x n_states
  // 'stateEstimator:11' D = zeros(no,ni);
  //  n_outputs x n_inputs
  // 'stateEstimator:12' Q = zeros(ns,ns);
  //  n_states  x n_states
  // 'stateEstimator:13' G = eye(ns);
  //  n_states  x n_states
  // 'stateEstimator:14' R = zeros(no,no);
  //  n_outputs x n_outputs
  // 'stateEstimator:15' N = zeros(ns,no);
  //  n_states  x n_outputs
  // 'stateEstimator:16' H = zeros(no,ns);
  //  n_outputs x n_states
  //  combine plant and output disturbance model
  // 'stateEstimator:19' A = [Ap, zeros(nsp,no);
  // 'stateEstimator:20'      zeros(no,nsp), Aod];
  rtb_A_m[0] = rtb_A_kh[0];
  rtb_A_m[8] = 0.0;
  rtb_A_m[2] = 0.0;
  rtb_A_m[10] = rtP->Aod1[0];
  rtb_A_m[1] = rtb_A_kh[1];
  rtb_A_m[9] = 0.0;
  rtb_A_m[3] = 0.0;
  rtb_A_m[11] = rtP->Aod1[1];
  rtb_A_m[4] = rtb_A_kh[2];
  rtb_A_m[12] = 0.0;
  rtb_A_m[6] = 0.0;
  rtb_A_m[14] = rtP->Aod1[2];
  rtb_A_m[5] = rtb_A_kh[3];
  rtb_A_m[13] = 0.0;
  rtb_A_m[7] = 0.0;
  rtb_A_m[15] = rtP->Aod1[3];

  // 'stateEstimator:21' B = [Bp;
  // 'stateEstimator:22'      zeros(no,ni)];
  // 'stateEstimator:23' C = [Cp Cod];
  rtb_C_d1[0] = rtb_C_j[0];
  rtb_C_d1[4] = rtP->Cod1[0];
  rtb_C_d1[1] = rtb_C_j[1];
  rtb_C_d1[5] = rtP->Cod1[1];
  rtb_C_d1[2] = rtb_C_j[2];
  rtb_C_d1[6] = rtP->Cod1[2];
  rtb_C_d1[3] = rtb_C_j[3];
  rtb_C_d1[7] = rtP->Cod1[3];

  // 'stateEstimator:24' D = Dp;
  // 'stateEstimator:26' B_est = [ [Bp; zeros(no,ni)] [zeros(size(Bp,1), size(Cp,1)); Bod] [zeros(size(Bp,1) + size(Bod,1), size(Cp,1))] ]; 
  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    B_est[i] = rtb_B_fo[i_0];
    B_est[i + 2] = 0.0;
    B_est[i + 1] = rtb_B_fo[i_0 + 1];
    B_est[i + 3] = 0.0;
    i += 4;
    i_0 += 2;
  }

  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 2; b_ct++) {
    B_est[i + 12] = 0.0;
    B_est[i + 14] = rtP->Bod1[i_0];
    B_est[i + 13] = 0.0;
    B_est[i + 15] = rtP->Bod1[i_0 + 1];
    B_est[i + 20] = 0.0;
    B_est[i + 21] = 0.0;
    B_est[i + 22] = 0.0;
    B_est[i + 23] = 0.0;
    i += 4;
    i_0 += 2;
  }

  // 'stateEstimator:27' D_est = [Dp Dod Dn];
  for (i = 0; i < 6; i++) {
    D_est[i] = rtb_D_oi[i];
  }

  // 'stateEstimator:28' Q = B_est * B_est';
  for (i = 0; i < 4; i++) {
    D_est[i + 6] = rtP->Dod1[i];
    D_est[i + 10] = rtP->Dmn1[i];
    i_0 = 0;
    for (int32_T b_ct{0}; b_ct < 4; b_ct++) {
      int32_T tmp_2;
      b_utarget_tmp = i_0 + i;
      rtb_Q[b_utarget_tmp] = 0.0;
      tmp_2 = 0;
      for (old_mvoff = 0; old_mvoff < 7; old_mvoff++) {
        rtb_Q[b_utarget_tmp] += B_est[tmp_2 + i] * B_est[tmp_2 + b_ct];
        tmp_2 += 4;
      }

      i_0 += 4;
    }
  }

  // 'stateEstimator:29' R = D_est * D_est';
  i = 0;
  for (i_0 = 0; i_0 < 7; i_0++) {
    rtb_R_h_tmp[i_0] = D_est[i];
    rtb_R_h_tmp[i_0 + 7] = D_est[i + 1];
    i += 2;
  }

  // 'stateEstimator:30' N = B_est * D_est';
  for (i = 0; i < 2; i++) {
    for (i_0 = 0; i_0 < 2; i_0++) {
      old_mvoff = (i << 1UL) + i_0;
      rtb_C_j[old_mvoff] = 0.0;
      for (int32_T b_ct{0}; b_ct < 7; b_ct++) {
        rtb_C_j[old_mvoff] += D_est[(b_ct << 1UL) + i_0] * rtb_R_h_tmp[7 * i +
          b_ct];
      }
    }

    for (i_0 = 0; i_0 < 4; i_0++) {
      old_mvoff = (i << 2UL) + i_0;
      rtb_N[old_mvoff] = 0.0;
      for (int32_T b_ct{0}; b_ct < 7; b_ct++) {
        rtb_N[old_mvoff] += B_est[(b_ct << 2UL) + i_0] * rtb_R_h_tmp[7 * i +
          b_ct];
      }
    }
  }

  // Outputs for Atomic SubSystem: '<S208>/ScalarExpansionR'
  //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
  //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
  //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
  //  yhat = C*xhat + D*u;
  ScalarExpansionR(rtb_C_j, rtb_A_kh);

  // End of Outputs for SubSystem: '<S208>/ScalarExpansionR'

  // Outputs for Atomic SubSystem: '<S208>/ScalarExpansionQ'
  ScalarExpansionQ(rtb_Q, rtb_Z_c);

  // End of Outputs for SubSystem: '<S208>/ScalarExpansionQ'

  // Outputs for Atomic SubSystem: '<S208>/ReducedQRN'
  // Constant: '<S208>/G' incorporates:
  //   Constant: '<S208>/H'

  ReducedQRN(localP->G_Value, localP->H_Value, rtb_Z_c, rtb_A_kh, rtb_N,
             rtb_Product_on, rtb_C_j, rtb_Product2_o);

  // End of Outputs for SubSystem: '<S208>/ReducedQRN'

  // Outputs for Atomic SubSystem: '<S208>/CalculatePL'
  // DataTypeConversion: '<S208>/DataTypeConversionEnable' incorporates:
  //   Constant: '<S3>/Constant'
  //   Delay: '<S208>/MemoryP'

  CalculatePL(rtb_A_m, rtb_C_d1, rtb_Product_on, rtb_C_j, rtb_Product2_o,
              localP->Constant_Value_f != 0.0, localDW->MemoryP_DSTATE, rtb_N,
              rtb_L_i, rtb_Z_c, rtb_Q);

  // End of Outputs for SubSystem: '<S208>/CalculatePL'

  // MATLAB Function: '<S249>/SqrtUsedFcn' incorporates:
  //   Constant: '<S249>/isSqrtUsed'

  SqrtUsedFcn(rtb_Z_c, localP->isSqrtUsed_Value, rtb_Product_on);

  // RandomNumber: '<S3>/Measurement Noise'
  rty_u[0] = localDW->NextOutput[0];

  // Sum: '<S3>/Sum' incorporates:
  //   Gain: '<S113>/u_scale'
  //   Product: '<S3>/Product'

  rty_u[0] = localP->u_scale_Gain[0] * rtb_D_c[0] + rtu_paramEst * rty_u[0];

  // Sum: '<S112>/Sum1' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'

  rtb_B_n[0] = rty_u[0] - rtu_u0[0];

  // RandomNumber: '<S3>/Measurement Noise'
  rty_u[1] = localDW->NextOutput[1];

  // Sum: '<S3>/Sum' incorporates:
  //   Gain: '<S113>/u_scale'
  //   Product: '<S3>/Product'

  rty_u[1] = localP->u_scale_Gain[1] * rtb_D_c[1] + rtu_paramEst * rty_u[1];

  // Sum: '<S112>/Sum1' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'

  rtb_B_n[1] = rty_u[1] - rtu_u0[1];

  // RandomNumber: '<S3>/Measurement Noise'
  rty_u[2] = localDW->NextOutput[2];

  // Sum: '<S3>/Sum' incorporates:
  //   Gain: '<S113>/u_scale'
  //   Product: '<S3>/Product'

  rty_u[2] = localP->u_scale_Gain[2] * rtb_D_c[2] + rtu_paramEst * rty_u[2];

  // Sum: '<S112>/Sum1' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'
  //   Sum: '<S111>/Add4'

  rtb_e_a = rty_u[2] - rtu_u0[2];
  rtb_B_n[2] = rtb_e_a;

  // Sum: '<S112>/Sum6' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'
  //   Sum: '<S111>/Add2'

  rtb_C = rtu_y[0] - rtu_y0[0];
  rtb_Add1_gs[0] = rtb_C;
  rtb_e_h = rtu_y[1] - rtu_y0[1];
  rtb_Add1_gs[1] = rtb_e_h;

  // Outputs for Enabled SubSystem: '<S227>/MeasurementUpdate'
  // BusCreator: '<S111>/Bus Creator2'
  MeasurementUpdate(rtb_DataTypeConversionEnable_k, rtb_L_i, rtb_Add1_gs,
                    rtb_C_d1, rtb_MemoryX_c, rtb_D_oi, rtb_B_n,
                    localDW->Product3, &localDW->MeasurementUpdate_n,
                    &localP->MeasurementUpdate_n);

  // End of Outputs for SubSystem: '<S227>/MeasurementUpdate'

  // MATLAB Function: '<S112>/MATLAB Function' incorporates:
  //   BusCreator: '<S111>/Bus Creator2'

  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    rtb_B_k[i] = rtb_B_fo[i_0];
    rtb_B_k[i + 2] = 0.0;
    rtb_B_k[i + 1] = rtb_B_fo[i_0 + 1];
    rtb_B_k[i + 3] = 0.0;
    i += 4;
    i_0 += 2;
  }

  for (i = 0; i <= 2; i += 2) {
    // Product: '<S227>/B[k]*u[k]' incorporates:
    //   Sum: '<S112>/Sum1'

    tmp_5 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_A_kh[i], tmp_5);
    tmp_3 = _mm_loadu_pd(&rtb_B_k[i]);
    tmp_4 = _mm_loadu_pd(&rtb_A_kh[i]);
    (void)_mm_storeu_pd(&rtb_A_kh[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_B_n[0])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_B_k[i + 4]);
    tmp_4 = _mm_loadu_pd(&rtb_A_kh[i]);
    (void)_mm_storeu_pd(&rtb_A_kh[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_B_n[1])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_B_k[i + 8]);
    tmp_4 = _mm_loadu_pd(&rtb_A_kh[i]);
    (void)_mm_storeu_pd(&rtb_A_kh[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_e_a)), tmp_4));

    // Product: '<S227>/A[k]*xhat[k|k-1]' incorporates:
    //   Product: '<S227>/B[k]*u[k]'

    (void)_mm_storeu_pd(&rtb_C_j[i], tmp_5);
    tmp_5 = _mm_loadu_pd(&rtb_A_m[i]);
    tmp_3 = _mm_loadu_pd(&rtb_C_j[i]);
    (void)_mm_storeu_pd(&rtb_C_j[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_c[0])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_m[i + 4]);
    tmp_3 = _mm_loadu_pd(&rtb_C_j[i]);
    (void)_mm_storeu_pd(&rtb_C_j[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_c[1])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_m[i + 8]);
    tmp_3 = _mm_loadu_pd(&rtb_C_j[i]);
    (void)_mm_storeu_pd(&rtb_C_j[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_c[2])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_m[i + 12]);
    tmp_3 = _mm_loadu_pd(&rtb_C_j[i]);
    (void)_mm_storeu_pd(&rtb_C_j[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_c[3])), tmp_3));

    // Product: '<S227>/B[k]*u[k]'
    tmp_5 = _mm_loadu_pd(&rtb_A_kh[i]);

    // Product: '<S227>/A[k]*xhat[k|k-1]' incorporates:
    //   Product: '<S227>/B[k]*u[k]'

    tmp_3 = _mm_loadu_pd(&rtb_C_j[i]);

    // Sum: '<S227>/Add' incorporates:
    //   Product: '<S227>/B[k]*u[k]'

    tmp_4 = _mm_loadu_pd(&localDW->Product3[i]);

    // Update for Delay: '<S208>/MemoryX' incorporates:
    //   Product: '<S227>/B[k]*u[k]'

    (void)_mm_storeu_pd(&localDW->MemoryX_DSTATE[i], _mm_add_pd(_mm_add_pd(tmp_5,
      tmp_3), tmp_4));
  }

  for (i = 0; i <= 0; i += 2) {
    // Product: '<S211>/Product'
    tmp_5 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_Add1_gs[i], tmp_5);
    tmp_3 = _mm_loadu_pd(&rtb_C_d1[i]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_gs[i]);
    (void)_mm_storeu_pd(&rtb_Add1_gs[i], _mm_add_pd(_mm_mul_pd(tmp_3,
      _mm_set1_pd(rtb_MemoryX_c[0])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_d1[i + 2]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_gs[i]);
    (void)_mm_storeu_pd(&rtb_Add1_gs[i], _mm_add_pd(_mm_mul_pd(tmp_3,
      _mm_set1_pd(rtb_MemoryX_c[1])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_d1[i + 4]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_gs[i]);
    (void)_mm_storeu_pd(&rtb_Add1_gs[i], _mm_add_pd(_mm_mul_pd(tmp_3,
      _mm_set1_pd(rtb_MemoryX_c[2])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_d1[i + 6]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_gs[i]);
    (void)_mm_storeu_pd(&rtb_Add1_gs[i], _mm_add_pd(_mm_mul_pd(tmp_3,
      _mm_set1_pd(rtb_MemoryX_c[3])), tmp_4));

    // Product: '<S211>/Product1' incorporates:
    //   BusCreator: '<S111>/Bus Creator2'
    //   Product: '<S211>/Product'
    //   Sum: '<S112>/Sum1'

    (void)_mm_storeu_pd(&Y[i], tmp_5);
    tmp_5 = _mm_loadu_pd(&rtb_D_oi[i]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd(rtb_B_n
      [0])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_D_oi[i + 2]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd(rtb_B_n
      [1])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_D_oi[i + 4]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd(rtb_e_a)),
      tmp_3));

    // Product: '<S211>/Product'
    tmp_5 = _mm_loadu_pd(&rtb_Add1_gs[i]);

    // Product: '<S211>/Product1' incorporates:
    //   Product: '<S211>/Product'

    tmp_3 = _mm_loadu_pd(&Y[i]);

    // BusCreator: '<S111>/Bus Creator2' incorporates:
    //   Product: '<S211>/Product'

    tmp_4 = _mm_loadu_pd(&rtu_y0[i]);

    // Sum: '<S112>/Sum3' incorporates:
    //   Product: '<S211>/Product'

    (void)_mm_storeu_pd(&rty_yhat[i], _mm_add_pd(_mm_add_pd(tmp_5, tmp_3), tmp_4));
  }

  // Update for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[0] = rtb_D_c[0];

  // Update for UnitDelay: '<S111>/Unit Delay5' incorporates:
  //   Sum: '<S111>/Add4'

  localDW->UnitDelay5_DSTATE[0] = rty_u[0] - rtu_u0[0];

  // Update for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[1] = rtb_D_c[1];

  // Update for UnitDelay: '<S111>/Unit Delay5' incorporates:
  //   Sum: '<S111>/Add4'

  localDW->UnitDelay5_DSTATE[1] = rty_u[1] - rtu_u0[1];

  // Update for UnitDelay: '<S113>/last_mv'
  localDW->last_mv_DSTATE[2] = rtb_D_c[2];

  // Update for UnitDelay: '<S111>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[2] = rtb_e_a;

  // Update for UnitDelay: '<S111>/Unit Delay4'
  localDW->UnitDelay4_DSTATE[0] = rtb_C;
  localDW->UnitDelay4_DSTATE[1] = rtb_e_h;

  // Update for Delay: '<S146>/delayTheta'
  localDW->icLoad = false;

  // Update for Delay: '<S146>/delayPhi'
  localDW->icLoad_o = false;

  // Update for Delay: '<S146>/delayBufferU'
  localDW->icLoad_a = false;

  // Update for Delay: '<S146>/delayBufferY'
  localDW->icLoad_p = false;

  // Update for Delay: '<S146>/delayL'
  localDW->icLoad_n = false;

  // Update for Delay: '<S147>/delayTheta'
  localDW->icLoad_pm = false;

  // Update for Delay: '<S147>/delayPhi'
  localDW->icLoad_d = false;

  // Update for Delay: '<S147>/delayBufferU'
  localDW->icLoad_nq = false;

  // Update for Delay: '<S147>/delayBufferY'
  localDW->icLoad_k = false;

  // Update for Delay: '<S147>/delayL'
  localDW->icLoad_g = false;

  // Update for Delay: '<S208>/MemoryX'
  localDW->icLoad_f = false;

  // Update for Delay: '<S208>/MemoryP'
  localDW->icLoad_a3 = false;
  (void)std::memcpy(&localDW->MemoryP_DSTATE[0], &rtb_Q[0], sizeof(real_T) <<
                    4UL);

  // Update for RandomNumber: '<S3>/Measurement Noise'
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];
}

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
int32_T SupervisoryController::xpotrf_b(real_T b_A[16])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
real_T SupervisoryController::minimum_k(const real_T x[4])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::trisolve_a(const real_T b_A[16], real_T b_B[16])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
real_T SupervisoryController::norm_i(const real_T x[4])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
real_T SupervisoryController::xnrm2_i(int32_T n, const real_T x[16], int32_T ix0)
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::xgemv_b(int32_T b_m, int32_T n, const real_T b_A[16],
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::xgerc_h(int32_T b_m, int32_T n, real_T alpha1,
  int32_T ix0, const real_T y[4], real_T b_A[16], int32_T ia0)
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::KWIKfactor_g(const real_T b_Ac[344], const int32_T
  iC[86], int32_T nA, const real_T b_Linv[16], real_T b_D[16], real_T b_H[16],
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
      RLinv[knt] += b_Linv[b_coltop + 4] * b_Ac[b_lastv + 85];
      RLinv[knt] += b_Linv[b_coltop + 8] * b_Ac[b_lastv + 171];
      RLinv[knt] += b_Linv[b_coltop + 12] * b_Ac[b_lastv + 257];
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
      beta1 = xnrm2_i(3 - k_i, b_A, ii + 2);
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

          beta1 = rt_hypotd_snf(atmp, xnrm2_i(3 - k_i, b_A, ii + 2));
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
        xgemv_b(b_lastv, knt, b_A, ii + 5, b_A, ii + 1, work);
        xgerc_h(b_lastv, knt, -tau[k_i], ii + 1, work, b_A, ii + 5);
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
        xgemv_b(knt, b_coltop, b_A, b_lastv, b_A, b_lastv - 4, work);
        xgerc_h(knt, b_coltop, -tau[k_i], b_lastv - 4, work, b_A, b_lastv);
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
real_T SupervisoryController::maximum_m(const real_T x[4])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::DropConstraint_a(int32_T kDrop, boolean_T iA[86],
  int32_T *nA, int32_T iC[86])
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::qpkwik_l(const real_T b_Linv[16], const real_T
  b_Hinv[16], const real_T f[4], const real_T b_Ac[344], const real_T b[86],
  boolean_T iA[86], int32_T maxiter, real_T FeasTol, real_T x[4], real_T lambda
  [86], int32_T *status)
{
  __m128d tmp_3;
  real_T cTol[86];
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
  int32_T iC[86];
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
  for (i = 0; i < 86; i++) {
    lambda[i] = 0.0;
    cTol[i] = 1.0;
    iC[i] = 0;
  }

  nA = 0;
  for (tmp = 0; tmp < 86; tmp++) {
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
        KWIKfactor_g(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees_m, RLinv, &Xnorm0);
        if (Xnorm0 < 0.0) {
          if (ColdReset) {
            *status = -2;
            exitg3 = 2;
          } else {
            nA = 0;
            (void)std::memset(&iA[0], 0, 86U * sizeof(boolean_T));
            (void)std::memset(&iC[0], 0, 86U * sizeof(int32_T));
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
              (void)std::memset(&iA[0], 0, 86U * sizeof(boolean_T));
              (void)std::memset(&iC[0], 0, 86U * sizeof(int32_T));
              ColdReset = true;
            } else {
              lambda[iC[i] - 1] = 0.0;
              DropConstraint_a(i + 1, iA, &nA, iC);
            }
          }
        }
      } else {
        if (nA <= 0) {
          (void)std::memset(&lambda[0], 0, 86U * sizeof(real_T));
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
    Xnorm0 = norm_i(x);
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (*status <= maxiter)) {
      real_T cVal;
      real_T t;
      cMin = -FeasTol;
      tmp = -1;
      for (i = 0; i < 86; i++) {
        t = cTol[i];
        if (!cTolComputed) {
          z[0] = std::abs(b_Ac[i] * x[0]);
          z[1] = std::abs(b_Ac[i + 86] * x[1]);
          z[2] = std::abs(b_Ac[i + 172] * x[2]);
          z[3] = std::abs(b_Ac[i + 258] * x[3]);
          t = std::fmax(t, maximum_m(z));
        }

        if (!iA[i]) {
          cVal = ((((b_Ac[i + 86] * x[1] + b_Ac[i] * x[0]) + b_Ac[i + 172] * x[2])
                   + b_Ac[i + 258] * x[3]) - b[i]) / t;
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
                  (&b_Hinv[iC_0 + 12]), _mm_set1_pd(b_Ac[tmp + 258])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0 + 8]),
                  _mm_set1_pd(b_Ac[tmp + 172])), _mm_add_pd(_mm_mul_pd
                  (_mm_loadu_pd(&b_Hinv[iC_0 + 4]), _mm_set1_pd(b_Ac[tmp + 86])),
                  _mm_add_pd(_mm_mul_pd(_mm_loadu_pd(&b_Hinv[iC_0]), _mm_set1_pd
                  (b_Ac[tmp])), _mm_set1_pd(0.0))))));
              }

              guard2 = true;
            } else {
              KWIKfactor_g(b_Ac, iC, nA, b_Linv, b_D, b_H, degrees_m, RLinv,
                           &cMin);
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
                    _mm_set1_pd(b_Ac[tmp + 258])), _mm_add_pd(_mm_mul_pd(tmp_1,
                    _mm_set1_pd(b_Ac[tmp + 172])), _mm_add_pd(_mm_mul_pd(tmp_0,
                    _mm_set1_pd(b_Ac[tmp + 86])), _mm_add_pd(_mm_mul_pd(tmp_3,
                    _mm_set1_pd(b_Ac[tmp])), _mm_set1_pd(0.0))))));
                }

                for (i = 0; i < nA; i++) {
                  iSave = i << 2UL;
                  r[i] = ((b_D[iSave + 1] * b_Ac[tmp + 86] + b_D[iSave] *
                           b_Ac[tmp]) + b_D[iSave + 2] * b_Ac[tmp + 172]) +
                    b_D[iSave + 3] * b_Ac[tmp + 258];
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

              t = b_Ac[tmp + 86];
              cVal_tmp = b_Ac[tmp + 172];
              cVal_tmp_0 = b_Ac[tmp + 258];
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
                  if ((iC_0 <= 86) && (lambda[iC_0 - 1] < 0.0)) {
                    lambda[iC_0 - 1] = 0.0;
                  }
                }

                lambda[tmp] += t;
                (void)std::frexp(1.0, &exponent);
                if (std::abs(t - cMin) < 2.2204460492503131E-16) {
                  DropConstraint_a(i, iA, &nA, iC);
                }

                if (!ColdReset) {
                  x[0] += t * z[0];
                  x[1] += t * z[1];
                  x[2] += t * z[2];
                  x[3] += t * z[3];
                  (void)std::frexp(1.0, &b_exponent);
                  if (std::abs(t - cVal) < 2.2204460492503131E-16) {
                    if (nA == static_cast<int32_T>(degrees_m)) {
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
            cMin = norm_i(x);
            if (std::abs(cMin - Xnorm0) > 0.001) {
              Xnorm0 = cMin;
              for (tmp = 0; tmp < 86; tmp++) {
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

// Function for MATLAB Function: '<S285>/FixedHorizonOptimizer'
void SupervisoryController::mpcblock_optimizer_l(const real_T rseq[40], const
  real_T vseq[21], const real_T x[4], const real_T old_u[3], const boolean_T iA
  [86], const real_T b_Mlim[86], real_T b_Mx[344], real_T b_Mu1[258], real_T
  b_Mv[1806], const real_T b_utarget[60], const real_T b_uoff[3], real_T b_H[16],
  real_T b_Ac[344], const real_T b_Wy[2], const real_T b_Wdu[3], const real_T
  b_Jm[180], const real_T b_Wu[3], const real_T b_I1[180], const real_T b_A[16],
  const real_T Bu[252], const real_T Bv[84], const real_T b_C[8], const real_T
  Dv[42], const int32_T b_Mrows[86], real_T u[3], real_T useq[63], real_T
  *status, boolean_T iAout[86], DW_State2controlLawAMPC2 *localDW)
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
  __m128d tmp_0;
  real_T I2Jm[180];
  real_T WduJm[180];
  real_T WuI2Jm[180];
  real_T b_Sx[160];
  real_T Sum_0[120];
  real_T WySuJm[120];
  real_T b_Su1[120];
  real_T b_Mlim_0[86];
  real_T b_Mlim_1[86];
  real_T b_Mu1_0[86];
  real_T b_Kv[63];
  real_T CA_0[42];
  real_T L[16];
  real_T b_Kx[12];
  real_T b_B[9];
  real_T b_I1_0[9];
  real_T b_Jm_0[9];
  real_T CA[8];
  real_T CA_1[8];
  real_T Sum[6];
  real_T varargin_1[4];
  real_T zopt[4];
  real_T b_C_0[2];
  real_T normH;
  real_T s;
  int32_T CA_tmp;
  int32_T Tries;
  int32_T a_tmp;
  int32_T b_SuJm_tmp;
  int32_T i;
  int32_T i1;
  int32_T kidx;
  int16_T ixw;
  int8_T a[3600];
  int8_T b[16];
  int8_T rows[2];
  boolean_T exitg1;
  boolean_T guard1{ false };

  (void)std::memset(&useq[0], 0, 63U * sizeof(real_T));
  (void)std::memset(&iAout[0], 0, 86U * sizeof(boolean_T));
  for (i = 0; i < 2; i++) {
    for (Tries = 0; Tries < 4; Tries++) {
      CA_tmp = (Tries << 1UL) + i;
      CA[CA_tmp] = 0.0;
      kidx = Tries << 2UL;
      CA[CA_tmp] += b_A[kidx] * b_C[i];
      CA[CA_tmp] += b_A[kidx + 1] * b_C[i + 2];
      CA[CA_tmp] += b_A[kidx + 2] * b_C[i + 4];
      CA[CA_tmp] += b_A[kidx + 3] * b_C[i + 6];
    }

    for (Tries = 0; Tries < 3; Tries++) {
      CA_tmp = (Tries << 1UL) + i;
      Sum[CA_tmp] = 0.0;
      kidx = Tries << 2UL;
      Sum[CA_tmp] += Bu[kidx] * b_C[i];
      Sum[CA_tmp] += Bu[kidx + 1] * b_C[i + 2];
      Sum[CA_tmp] += Bu[kidx + 2] * b_C[i + 4];
      Sum[CA_tmp] += Bu[kidx + 3] * b_C[i + 6];
    }

    b_C_0[i] = 0.0;
    b_C_0[i] += b_C[i] * Bv[0];
    b_C_0[i] += b_C[i + 2] * Bv[1];
    b_C_0[i] += b_C[i + 4] * Bv[2];
    b_C_0[i] += b_C[i + 6] * Bv[3];
    localDW->b_Hv[i] = b_C_0[i];
    localDW->b_Hv[i + 40] = Dv[i];
  }

  i = 0;
  for (Tries = 0; Tries < 19; Tries++) {
    localDW->b_Hv[i + 80] = 0.0;
    localDW->b_Hv[i + 81] = 0.0;
    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 21; Tries++) {
    (void)std::memset(&localDW->b_Hv[i + 2], 0, 38U * sizeof(real_T));
    i += 40;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 4; i1++) {
    b_Sx[i] = CA[Tries];
    b_Sx[i + 1] = CA[Tries + 1];
    i += 40;
    Tries += 2;
  }

  for (i = 0; i < 38; i++) {
    b_Sx[i + 2] = 0.0;
    b_Sx[i + 42] = 0.0;
    b_Sx[i + 82] = 0.0;
    b_Sx[i + 122] = 0.0;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    b_Su1[i] = Sum[Tries];
    b_Su1[i + 1] = Sum[Tries + 1];
    i += 40;
    Tries += 2;
  }

  for (i = 0; i < 38; i++) {
    b_Su1[i + 2] = 0.0;
    b_Su1[i + 42] = 0.0;
    b_Su1[i + 82] = 0.0;
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    localDW->Su[i] = Sum[Tries];
    localDW->Su[i + 1] = Sum[Tries + 1];
    i += 40;
    Tries += 2;
  }

  i = 0;
  for (Tries = 0; Tries < 57; Tries++) {
    localDW->Su[i + 120] = 0.0;
    localDW->Su[i + 121] = 0.0;
    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 60; Tries++) {
    (void)std::memset(&localDW->Su[i + 2], 0, 38U * sizeof(real_T));
    i += 40;
  }

  for (kidx = 0; kidx < 19; kidx++) {
    int8_T kidx_0;
    kidx_0 = static_cast<int8_T>(((kidx + 1) << 1UL) + 1);
    for (i = 0; i < 2; i++) {
      rows[i] = static_cast<int8_T>(i + static_cast<int32_T>(kidx_0));
      Tries = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 3; a_tmp++) {
        CA_tmp = Tries + i;
        Sum[CA_tmp] += ((Bu[i1 + 1] * CA[i + 2] + Bu[i1] * CA[i]) + Bu[i1 + 2] *
                        CA[i + 4]) + Bu[i1 + 3] * CA[i + 6];
        Tries += 2;
        i1 += 4;
      }
    }

    for (i = 0; i < 3; i++) {
      CA_tmp = i << 1UL;
      normH = Sum[CA_tmp];
      b_Su1[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = normH;
      Sum_0[CA_tmp] = normH;
      normH = Sum[CA_tmp + 1];
      b_Su1[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = normH;
      Sum_0[CA_tmp + 1] = normH;
    }

    for (i = 0; i < 57; i++) {
      CA_tmp = (i + 3) << 1UL;
      Sum_0[CA_tmp] = localDW->Su[(40 * i + static_cast<int32_T>(rows[0])) - 3];
      Sum_0[CA_tmp + 1] = localDW->Su[(40 * i + static_cast<int32_T>(rows[1])) -
        3];
    }

    for (i = 0; i < 60; i++) {
      Tries = i << 1UL;
      localDW->Su[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = Sum_0[Tries];
      localDW->Su[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = Sum_0[Tries +
        1];
    }

    for (i = 0; i <= 0; i += 2) {
      (void)_mm_storeu_pd(&b_C_0[i], _mm_set1_pd(0.0));
      tmp = _mm_loadu_pd(&CA[i]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[0])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 2]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[1])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 4]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[2])),
        tmp_0));
      tmp = _mm_loadu_pd(&CA[i + 6]);
      tmp_0 = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&b_C_0[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(Bv[3])),
        tmp_0));
      tmp = _mm_loadu_pd(&b_C_0[i]);
      (void)_mm_storeu_pd(&CA_0[i], tmp);
    }

    for (i = 0; i < 20; i++) {
      CA_tmp = (i + 1) << 1UL;
      CA_0[CA_tmp] = localDW->b_Hv[(40 * i + static_cast<int32_T>(rows[0])) - 3];
      CA_0[CA_tmp + 1] = localDW->b_Hv[(40 * i + static_cast<int32_T>(rows[1]))
        - 3];
    }

    for (i = 0; i < 21; i++) {
      Tries = i << 1UL;
      localDW->b_Hv[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = CA_0[Tries];
      localDW->b_Hv[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = CA_0[Tries +
        1];
    }

    for (i = 0; i < 2; i++) {
      Tries = 0;
      i1 = 0;
      for (a_tmp = 0; a_tmp < 4; a_tmp++) {
        CA_tmp = Tries + i;
        CA_1[CA_tmp] = 0.0;
        CA_1[CA_tmp] += b_A[i1] * CA[i];
        CA_1[CA_tmp] += b_A[i1 + 1] * CA[i + 2];
        CA_1[CA_tmp] += b_A[i1 + 2] * CA[i + 4];
        CA_1[CA_tmp] += b_A[i1 + 3] * CA[i + 6];
        Tries += 2;
        i1 += 4;
      }
    }

    (void)std::memcpy(&CA[0], &CA_1[0], sizeof(real_T) << 3UL);
    for (i = 0; i < 4; i++) {
      Tries = i << 1UL;
      b_Sx[(static_cast<int32_T>(rows[0]) + 40 * i) - 1] = CA[Tries];
      b_Sx[(static_cast<int32_T>(rows[1]) + 40 * i) - 1] = CA[Tries + 1];
    }
  }

  i = 0;
  Tries = 0;
  for (i1 = 0; i1 < 3; i1++) {
    for (a_tmp = 0; a_tmp < 40; a_tmp++) {
      b_SuJm_tmp = a_tmp + i;
      Sum_0[b_SuJm_tmp] = 0.0;
      kidx = 0;
      for (CA_tmp = 0; CA_tmp < 60; CA_tmp++) {
        Sum_0[b_SuJm_tmp] += localDW->Su[kidx + a_tmp] * b_Jm[CA_tmp + Tries];
        kidx += 40;
      }
    }

    i += 40;
    Tries += 60;
  }

  if (b_Mrows[0] > 0) {
    kidx = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (kidx < 86)) {
      if (b_Mrows[kidx] <= 40) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = -Sum_0[Tries - 1];
        b_Ac[kidx + 86] = -Sum_0[Tries + 39];
        b_Ac[kidx + 172] = -Sum_0[Tries + 79];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = -b_Sx[Tries - 1];
        b_Mx[kidx + 86] = -b_Sx[Tries + 39];
        b_Mx[kidx + 172] = -b_Sx[Tries + 79];
        b_Mx[kidx + 258] = -b_Sx[Tries + 119];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = -b_Su1[Tries - 1];
        b_Mu1[kidx + 86] = -b_Su1[Tries + 39];
        b_Mu1[kidx + 172] = -b_Su1[Tries + 79];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 86 * i] = -localDW->b_Hv[(40 * i + Tries) - 1];
        }

        kidx++;
      } else if (b_Mrows[kidx] <= 80) {
        Tries = b_Mrows[kidx];
        b_Ac[kidx] = Sum_0[Tries - 41];
        b_Ac[kidx + 86] = Sum_0[Tries - 1];
        b_Ac[kidx + 172] = Sum_0[Tries + 39];
        Tries = b_Mrows[kidx];
        b_Mx[kidx] = b_Sx[Tries - 41];
        b_Mx[kidx + 86] = b_Sx[Tries - 1];
        b_Mx[kidx + 172] = b_Sx[Tries + 39];
        b_Mx[kidx + 258] = b_Sx[Tries + 79];
        Tries = b_Mrows[kidx];
        b_Mu1[kidx] = b_Su1[Tries - 41];
        b_Mu1[kidx + 86] = b_Su1[Tries - 1];
        b_Mu1[kidx + 172] = b_Su1[Tries + 39];
        Tries = b_Mrows[kidx];
        for (i = 0; i < 21; i++) {
          b_Mv[kidx + 86 * i] = localDW->b_Hv[(40 * i + Tries) - 41];
        }

        kidx++;
      } else {
        exitg1 = true;
      }
    }
  }

  (void)std::memset(&b_B[0], 0, 9U * sizeof(real_T));
  b_B[0] = 1.0;
  b_B[4] = 1.0;
  b_B[8] = 1.0;
  kidx = -1;
  for (Tries = 0; Tries < 20; Tries++) {
    for (i = 0; i < 3; i++) {
      for (i1 = 0; i1 < 20; i1++) {
        a_tmp = static_cast<int32_T>(c_A[20 * Tries + i1]);
        a[kidx + 1] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i]) *
          a_tmp);
        a[kidx + 2] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 1]) *
          a_tmp);
        a[kidx + 3] = static_cast<int8_T>(static_cast<int32_T>(b_B[3 * i + 2]) *
          a_tmp);
        kidx += 3;
      }
    }
  }

  i = 0;
  for (Tries = 0; Tries < 3; Tries++) {
    for (i1 = 0; i1 < 60; i1++) {
      CA_tmp = i1 + i;
      I2Jm[CA_tmp] = 0.0;
      a_tmp = 0;
      for (kidx = 0; kidx < 60; kidx++) {
        I2Jm[CA_tmp] += static_cast<real_T>(a[a_tmp + i1]) * b_Jm[kidx + i];
        a_tmp += 60;
      }
    }

    i += 60;
  }

  ixw = 1;
  for (kidx = 0; kidx < 40; kidx++) {
    normH = b_Wy[ixw - 1];
    WySuJm[kidx] = normH * Sum_0[kidx];
    WySuJm[kidx + 40] = Sum_0[kidx + 40] * normH;
    WySuJm[kidx + 80] = Sum_0[kidx + 80] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 2) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    normH = b_Wu[ixw - 1];
    WuI2Jm[kidx] = normH * I2Jm[kidx];
    WuI2Jm[kidx + 60] = I2Jm[kidx + 60] * normH;
    WuI2Jm[kidx + 120] = I2Jm[kidx + 120] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  ixw = 1;
  for (kidx = 0; kidx < 60; kidx++) {
    normH = b_Wdu[ixw - 1];
    WduJm[kidx] = normH * b_Jm[kidx];
    WduJm[kidx + 60] = b_Jm[kidx + 60] * normH;
    WduJm[kidx + 120] = b_Jm[kidx + 120] * normH;
    ixw = static_cast<int16_T>(ixw + 1);
    if (ixw > 3) {
      ixw = 1;
    }
  }

  for (i = 0; i < 3; i++) {
    for (Tries = 0; Tries < 3; Tries++) {
      b_SuJm_tmp = 3 * Tries + i;
      b_B[b_SuJm_tmp] = 0.0;
      for (i1 = 0; i1 < 40; i1++) {
        b_B[b_SuJm_tmp] += Sum_0[40 * i + i1] * WySuJm[40 * Tries + i1];
      }

      b_Jm_0[b_SuJm_tmp] = 0.0;
      normH = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        a_tmp = 60 * i + i1;
        kidx = 60 * Tries + i1;
        normH += I2Jm[a_tmp] * WuI2Jm[kidx];
        b_Jm_0[b_SuJm_tmp] += b_Jm[a_tmp] * WduJm[kidx];
      }

      b_H[i + (Tries << 2UL)] = (b_B[b_SuJm_tmp] + b_Jm_0[b_SuJm_tmp]) + normH;
    }
  }

  for (i = 0; i < 3; i++) {
    for (Tries = 0; Tries < 3; Tries++) {
      kidx = 3 * Tries + i;
      b_Jm_0[kidx] = 0.0;
      for (i1 = 0; i1 < 40; i1++) {
        b_Jm_0[kidx] += b_Su1[40 * i + i1] * WySuJm[40 * Tries + i1];
      }

      b_I1_0[kidx] = 0.0;
      for (i1 = 0; i1 < 60; i1++) {
        b_I1_0[kidx] += b_I1[60 * i + i1] * WuI2Jm[60 * Tries + i1];
      }
    }
  }

  for (i = 0; i <= 6; i += 2) {
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
  for (Tries = 0; Tries < 4; Tries++) {
    i1 = 0;
    a_tmp = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      b_SuJm_tmp = i1 + Tries;
      b_Kx[b_SuJm_tmp] = 0.0;
      for (CA_tmp = 0; CA_tmp < 40; CA_tmp++) {
        b_Kx[b_SuJm_tmp] += b_Sx[CA_tmp + i] * WySuJm[CA_tmp + a_tmp];
      }

      i1 += 4;
      a_tmp += 40;
    }

    i += 40;
  }

  i = 0;
  for (Tries = 0; Tries < 21; Tries++) {
    i1 = 0;
    a_tmp = 0;
    for (kidx = 0; kidx < 3; kidx++) {
      b_SuJm_tmp = i1 + Tries;
      b_Kv[b_SuJm_tmp] = 0.0;
      for (CA_tmp = 0; CA_tmp < 40; CA_tmp++) {
        b_Kv[b_SuJm_tmp] += localDW->b_Hv[CA_tmp + i] * WySuJm[CA_tmp + a_tmp];
      }

      i1 += 21;
      a_tmp += 40;
    }

    i += 40;
  }

  for (i = 0; i <= 118; i += 2) {
    tmp = _mm_loadu_pd(&WySuJm[i]);
    (void)_mm_storeu_pd(&WySuJm[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  kidx = 0;
  (void)std::memcpy(&L[0], &b_H[0], sizeof(real_T) << 4UL);
  Tries = xpotrf_b(L);
  guard1 = false;
  if (Tries == 0) {
    varargin_1[0] = L[0];
    varargin_1[1] = L[5];
    varargin_1[2] = L[10];
    varargin_1[3] = L[15];
    if (minimum_k(varargin_1) > 1.4901161193847656E-7) {
    } else {
      guard1 = true;
    }
  } else {
    guard1 = true;
  }

  if (guard1) {
    boolean_T exitg2;
    normH = 0.0;
    Tries = 0;
    exitg2 = false;
    while (((exitg2 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (Tries < 4)) {
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
      kidx = 2;
    } else {
      Tries = 0;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (Tries <= 4)) {
        boolean_T guard2{ false };

        normH = rt_powd_snf(10.0, static_cast<real_T>(Tries)) *
          1.4901161193847656E-7;
        for (i = 0; i < 16; i++) {
          b[i] = 0;
        }

        b[0] = 1;
        b[5] = 1;
        b[10] = 1;
        b[15] = 1;
        for (i = 0; i < 16; i++) {
          b_H[i] += normH * static_cast<real_T>(b[i]);
          L[i] = b_H[i];
        }

        kidx = xpotrf_b(L);
        guard2 = false;
        if (kidx == 0) {
          varargin_1[0] = L[0];
          varargin_1[1] = L[5];
          varargin_1[2] = L[10];
          varargin_1[3] = L[15];
          if (minimum_k(varargin_1) > 1.4901161193847656E-7) {
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
    for (i = 0; i < 16; i++) {
      b[i] = 0;
    }

    b[0] = 1;
    b[5] = 1;
    b[10] = 1;
    b[15] = 1;
    i = 0;
    for (kidx = 0; kidx < 4; kidx++) {
      b_H[i] = static_cast<real_T>(b[i]);
      b_H[i + 1] = static_cast<real_T>(b[i + 1]);
      b_H[i + 2] = static_cast<real_T>(b[i + 2]);
      b_H[i + 3] = static_cast<real_T>(b[i + 3]);
      i += 4;
    }

    trisolve_a(L, b_H);
    varargin_1[0] = 0.0;
    varargin_1[1] = 0.0;
    varargin_1[2] = 0.0;
    varargin_1[3] = 0.0;
    for (kidx = 0; kidx < 3; kidx++) {
      real_T WuI2Jm_0;
      normH = 0.0;
      for (i = 0; i < 40; i++) {
        normH += WySuJm[40 * kidx + i] * rseq[i];
      }

      s = 0.0;
      for (i = 0; i < 21; i++) {
        s += b_Kv[21 * kidx + i] * vseq[i];
      }

      WuI2Jm_0 = 0.0;
      for (i = 0; i < 60; i++) {
        WuI2Jm_0 += WuI2Jm[60 * kidx + i] * b_utarget[i];
      }

      i = kidx << 2UL;
      varargin_1[kidx] = ((((((b_Kx[i + 1] * x[1] + b_Kx[i] * x[0]) + b_Kx[i + 2]
        * x[2]) + b_Kx[i + 3] * x[3]) + normH) + ((b_B[3 * kidx + 1] * old_u[1]
        + b_B[3 * kidx] * old_u[0]) + b_B[3 * kidx + 2] * old_u[2])) + s) +
        WuI2Jm_0;
    }

    for (i = 0; i < 86; i++) {
      iAout[i] = iA[i];
      b_Mlim_0[i] = (((b_Mx[i + 86] * x[1] + b_Mx[i] * x[0]) + b_Mx[i + 172] *
                      x[2]) + b_Mx[i + 258] * x[3]) + b_Mlim[i];
      b_Mu1_0[i] = 0.0;
      b_Mu1_0[i] += b_Mu1[i] * old_u[0];
      b_Mu1_0[i] += b_Mu1[i + 86] * old_u[1];
      b_Mu1_0[i] += b_Mu1[i + 172] * old_u[2];
    }

    i = 0;
    for (Tries = 0; Tries < 4; Tries++) {
      i1 = 0;
      for (a_tmp = 0; a_tmp < 4; a_tmp++) {
        kidx = i1 + Tries;
        L[kidx] = 0.0;
        L[kidx] += b_H[i] * b_H[i1];
        L[kidx] += b_H[i + 1] * b_H[i1 + 1];
        L[kidx] += b_H[i + 2] * b_H[i1 + 2];
        L[kidx] += b_H[i + 3] * b_H[i1 + 3];
        i1 += 4;
      }

      i += 4;
    }

    for (i = 0; i < 86; i++) {
      normH = 0.0;
      Tries = 0;
      for (i1 = 0; i1 < 21; i1++) {
        normH += b_Mv[Tries + i] * vseq[i1];
        Tries += 86;
      }

      b_Mlim_1[i] = -((b_Mlim_0[i] + b_Mu1_0[i]) + normH);
    }

    qpkwik_l(b_H, L, varargin_1, b_Ac, b_Mlim_1, iAout, 360, 1.0E-6, zopt,
             b_Mlim_0, &kidx);
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

// System initialize for function-call system: '<S1>/State2.controlLaw.AMPC2'
void SupervisoryController::State2controlLawAMPC2_Init(real_T rty_u[3], real_T
  rty_yhat[2], DW_State2controlLawAMPC2 *localDW, P_State2controlLawAMPC2
  *localP)
{
  int32_T i;
  int32_T t_tmp;
  uint32_T tseed;

  // InitializeConditions for Memory: '<S257>/Memory'
  (void)std::memcpy(&localDW->Memory_PreviousInput[0],
                    &localP->Memory_InitialCondition[0], 86U * sizeof(boolean_T));

  // InitializeConditions for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[0] = localP->last_mv_InitialCondition[0];

  // InitializeConditions for UnitDelay: '<S255>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[0] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[1] = localP->last_mv_InitialCondition[1];

  // InitializeConditions for UnitDelay: '<S255>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[1] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[2] = localP->last_mv_InitialCondition[2];

  // InitializeConditions for UnitDelay: '<S255>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[2] = localP->UnitDelay5_InitialCondition;

  // InitializeConditions for UnitDelay: '<S255>/Unit Delay4'
  localDW->UnitDelay4_DSTATE[0] = localP->UnitDelay4_InitialCondition[0];
  localDW->UnitDelay4_DSTATE[1] = localP->UnitDelay4_InitialCondition[1];

  // InitializeConditions for Delay: '<S290>/delayTheta'
  localDW->icLoad = true;

  // InitializeConditions for Delay: '<S290>/delayPhi'
  localDW->icLoad_o = true;

  // InitializeConditions for Delay: '<S290>/delayBufferU'
  localDW->icLoad_e = true;

  // InitializeConditions for Delay: '<S290>/delayBufferY'
  localDW->icLoad_d = true;

  // InitializeConditions for Delay: '<S290>/delayL'
  localDW->icLoad_e2 = true;

  // InitializeConditions for Delay: '<S291>/delayTheta'
  localDW->icLoad_f = true;

  // InitializeConditions for Delay: '<S291>/delayPhi'
  localDW->icLoad_ej = true;

  // InitializeConditions for Delay: '<S291>/delayBufferU'
  localDW->icLoad_j = true;

  // InitializeConditions for Delay: '<S291>/delayBufferY'
  localDW->icLoad_i = true;

  // InitializeConditions for Delay: '<S291>/delayL'
  localDW->icLoad_jb = true;

  // InitializeConditions for Delay: '<S352>/MemoryX'
  localDW->icLoad_l = true;

  // InitializeConditions for Delay: '<S352>/MemoryP'
  localDW->icLoad_h = true;

  // InitializeConditions for RandomNumber: '<S4>/Measurement Noise'
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  t_tmp = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) & 32768U));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[0] = tseed;
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[1] = tseed;
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];
  i = static_cast<int32_T>(static_cast<uint32_T>(static_cast<uint32_T>
    (localP->MeasurementNoise_Seed) >> 16UL));
  tseed = ((((static_cast<uint32_T>(localP->MeasurementNoise_Seed) - (
    static_cast<uint32_T>(i) << 16UL)) + static_cast<uint32_T>(t_tmp)) << 16UL)
           + static_cast<uint32_T>(t_tmp)) + static_cast<uint32_T>(i);
  if (tseed < 1U) {
    tseed = 1144108930U;
  } else if (tseed > 2147483646U) {
    tseed = 2147483646U;
  } else {
    // no actions
  }

  localDW->RandSeed[2] = tseed;
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];

  // End of InitializeConditions for RandomNumber: '<S4>/Measurement Noise'

  // SystemInitialize for Enabled SubSystem: '<S371>/MeasurementUpdate'
  MeasurementUpdate_Init(localDW->Product3, &localP->MeasurementUpdate_a);

  // End of SystemInitialize for SubSystem: '<S371>/MeasurementUpdate'

  // SystemInitialize for Outport: '<S4>/u'
  rty_u[0] = localP->u_Y0;
  rty_u[1] = localP->u_Y0;
  rty_u[2] = localP->u_Y0;

  // SystemInitialize for Outport: '<S4>/yhat'
  rty_yhat[0] = localP->yhat_Y0;
  rty_yhat[1] = localP->yhat_Y0;
}

// Disable for function-call system: '<S1>/State2.controlLaw.AMPC2'
void SupervisoryController::State2controlLawAMPC2_Disable
  (DW_State2controlLawAMPC2 *localDW, P_State2controlLawAMPC2 *localP)
{
  // Disable for Enabled SubSystem: '<S371>/MeasurementUpdate'
  if (localDW->MeasurementUpdate_a.MeasurementUpdate_MODE) {
    MeasurementUpdate_Disable(localDW->Product3, &localDW->MeasurementUpdate_a,
      &localP->MeasurementUpdate_a);
  }

  // End of Disable for SubSystem: '<S371>/MeasurementUpdate'
}

// Output and update for function-call system: '<S1>/State2.controlLaw.AMPC2'
void SupervisoryController::State2controlLawAMPC2(const real_T rtu_r[2], const
  real_T rtu_y[2], const real_T rtu_y0[2], const real_T rtu_u0[3], real_T
  rtu_paramEst, real_T rty_u[3], real_T rty_yhat[2], DW_State2controlLawAMPC2
  *localDW, P_State2controlLawAMPC2 *localP, P *rtP, ZCE_State2controlLawAMPC2
  *localZCE)
{
  static const real_T f[344]{ -0.99996947351391963, 1.9978365398021225E-5,
    -0.99993894863151356, 3.9955226117636394E-5, -0.99990842535268221,
    5.99305822582772E-5, -0.999877903677326, 7.9904433919368887E-5,
    -0.99984738360534542, 9.9876781200330258E-5, -0.99981686513664081,
    0.00011984762420057368, -0.9997863482711129, 0.00013981696301950509,
    -0.99975583300866211, 0.00015978479775652393, -0.99972531934918873,
    0.00017975112851102327, -0.99969480729259341, 0.0001997159553823897,
    -0.99966429683877656, 0.0002196792784700034, -0.9996337879876388,
    0.00023964109787323807, -0.99960328073908056, 0.000259601413691461,
    -0.99957277509300235, 0.000279560226024033, -0.99954227104930482,
    0.00029951753497030862, -0.99951176860788837, 0.00031947334062963569,
    -0.99948126776865354, 0.00033942764310135575, -0.99945076853150094,
    0.00035938044248480397, -0.999420270896331, 0.000379331738879309,
    -0.99938977486304448, 0.00039928153238419297, 0.99996947351391963,
    -1.9978365398021225E-5, 0.99993894863151356, -3.9955226117636394E-5,
    0.99990842535268221, -5.99305822582772E-5, 0.999877903677326,
    -7.9904433919368887E-5, 0.99984738360534542, -9.9876781200330258E-5,
    0.99981686513664081, -0.00011984762420057368, 0.9997863482711129,
    -0.00013981696301950509, 0.99975583300866211, -0.00015978479775652393,
    0.99972531934918873, -0.00017975112851102327, 0.99969480729259341,
    -0.0001997159553823897, 0.99966429683877656, -0.0002196792784700034,
    0.9996337879876388, -0.00023964109787323807, 0.99960328073908056,
    -0.000259601413691461, 0.99957277509300235, -0.000279560226024033,
    0.99954227104930482, -0.00029951753497030862, 0.99951176860788837,
    -0.00031947334062963569, 0.99948126776865354, -0.00033942764310135575,
    0.99945076853150094, -0.00035938044248480397, 0.999420270896331,
    -0.000379331738879309, 0.99938977486304448, -0.00039928153238419297, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 3.3626770675199147E-5, -0.99995521109485208,
    6.7251008737009146E-5, -0.9999104248675581, 0.00010087271435278937,
    -0.99986564131794753, 0.00013449188768988832, -0.99982086044584983,
    0.00016810852891564372, -0.99977608225109449, 0.00020172263819738238,
    -0.999731306733511, 0.0002353342157024203, -0.99968653389292883,
    0.00026894326159806266, -0.99964176372917768, 0.00030254977605160372,
    -0.99959699624208687, 0.00033615375923032705, -0.999552231431486,
    0.0003697552113015052, -0.99950746929720458, 0.0004033541324324,
    -0.99946270983907226, 0.00043695052279026244, -0.99941795305691861,
    0.00047054438254233265, -0.99937319895057319, 0.00050413571185583986,
    -0.99932844751986571, 0.00053772451089800256, -0.99928369876462564,
    0.00057131077983602844, -0.99923895268468255, 0.00060489451883711414,
    -0.99919420927986613, 0.00063847572806844581, -0.99914946855000608,
    0.00067205440769719843, -0.999104730494932, -3.3626770675199147E-5,
    0.99995521109485208, -6.7251008737009146E-5, 0.9999104248675581,
    -0.00010087271435278937, 0.99986564131794753, -0.00013449188768988832,
    0.99982086044584983, -0.00016810852891564372, 0.99977608225109449,
    -0.00020172263819738238, 0.999731306733511, -0.0002353342157024203,
    0.99968653389292883, -0.00026894326159806266, 0.99964176372917768,
    -0.00030254977605160372, 0.99959699624208687, -0.00033615375923032705,
    0.999552231431486, -0.0003697552113015052, 0.99950746929720458,
    -0.0004033541324324, 0.99946270983907226, -0.00043695052279026244,
    0.99941795305691861, -0.00047054438254233265, 0.99937319895057319,
    -0.00050413571185583986, 0.99932844751986571, -0.00053772451089800256,
    0.99928369876462564, -0.00057131077983602844, 0.99923895268468255,
    -0.00060489451883711414, 0.99919420927986613, -0.00063847572806844581,
    0.99914946855000608, -0.00067205440769719843, 0.999104730494932, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0,
    -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, -0.0,
    -1.0, -0.0, -1.0, -0.0, -1.0, -0.0, -1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
    0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T k[344]{ -0.70746101068200651, 0.35521300513038595,
    -1.4149123697315824, 0.71042423457375414, -2.1223540773836476,
    1.0656336882168116, -2.8297861338731107, 1.4208413659462755,
    -3.53720853943487, 1.7760472676488726, -4.244621294303812,
    2.1312513932113393, -4.952024398714812, 2.4864537425204221,
    -5.6594178529027355, 2.8416543154628768, -6.3668016571024371,
    3.1968531119254697, -7.0741758115487583, 3.5520501317949762,
    -7.7815403164765327, 3.9072453749581815, -8.48889517212058,
    4.262438841301881, -9.19624037871571, 4.6176305307128791,
    -9.9035759364967237, 4.9728204430779908, -10.610901845698407,
    5.3280085782840407, -11.318218106555539, 5.6831949362178626,
    -12.025524719302885, 6.0383795167663008, -12.732821684175198,
    6.3935623198162084, -13.440109001407224, 6.7487433452544492,
    -14.147386671233697, 7.103922592967896, 0.70746101068200651,
    -0.35521300513038595, 1.4149123697315824, -0.71042423457375414,
    2.1223540773836476, -1.0656336882168116, 2.8297861338731107,
    -1.4208413659462755, 3.53720853943487, -1.7760472676488726,
    4.244621294303812, -2.1312513932113393, 4.952024398714812,
    -2.4864537425204221, 5.6594178529027355, -2.8416543154628768,
    6.3668016571024371, -3.1968531119254697, 7.0741758115487583,
    -3.5520501317949762, 7.7815403164765327, -3.9072453749581815,
    8.48889517212058, -4.262438841301881, 9.19624037871571, -4.6176305307128791,
    9.9035759364967237, -4.9728204430779908, 10.610901845698407,
    -5.3280085782840407, 11.318218106555539, -5.6831949362178626,
    12.025524719302885, -6.0383795167663008, 12.732821684175198,
    -6.3935623198162084, 13.440109001407224, -6.7487433452544492,
    14.147386671233697, -7.103922592967896, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.34920656360027097, 0.074930066385699018, 0.69839994749507994,
    0.14984980015943497, 1.0475801524342094, 0.22475920204730321,
    1.3967471791673947, 0.29965827277535156, 1.7459010284443235,
    0.37454701306958038, 2.0950417010146363, 0.44942542365594257,
    2.4441691976279265, 0.52429350526034346, 2.79328351903374,
    0.59915125860864094, 3.142384665981576, 0.67399868442664534,
    3.491472639220885, 0.74883578344011958, 3.8405474395010719,
    0.82366255637477914, 4.1896090675714932, 0.898479003956292,
    4.538657524181458, 0.97328512691027858, 4.887692810080229,
    1.0480809259623118, 5.2367149260170214, 1.1228664018379175,
    5.5857238727410028, 1.1976415552625737, 5.9347196510012941,
    1.2724063869617108, 6.283702261546968, 1.3471608976607121, 6.63267170512705,
    1.4219050880849133, 6.9816279824905205, 1.4966389589596025,
    -0.34920656360027097, -0.074930066385699018, -0.69839994749507994,
    -0.14984980015943497, -1.0475801524342094, -0.22475920204730321,
    -1.3967471791673947, -0.29965827277535156, -1.7459010284443235,
    -0.37454701306958038, -2.0950417010146363, -0.44942542365594257,
    -2.4441691976279265, -0.52429350526034346, -2.79328351903374,
    -0.59915125860864094, -3.142384665981576, -0.67399868442664534,
    -3.491472639220885, -0.74883578344011958, -3.8405474395010719,
    -0.82366255637477914, -4.1896090675714932, -0.898479003956292,
    -4.538657524181458, -0.97328512691027858, -4.887692810080229,
    -1.0480809259623118, -5.2367149260170214, -1.1228664018379175,
    -5.5857238727410028, -1.1976415552625737, -5.9347196510012941,
    -1.2724063869617108, -6.283702261546968, -1.3471608976607121,
    -6.63267170512705, -1.4219050880849133, -6.9816279824905205,
    -1.4966389589596025, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.34258808392522866,
    -0.43035445788055782, 0.68518018127074631, -0.86069648500604323,
    1.0277762914960311, -1.2910260820133979, 1.3703764140605992,
    -1.7213432495395238, 1.7129805484240048, -2.151647988221284,
    2.0555886940458397, -2.5819402986955016, 2.3982008503857339,
    -3.0122201815989613, 2.740817016903355, -3.4424876375684077,
    3.0834371930584092, -3.8727426672405465, 3.4260613783106395,
    -4.3029852712520436, 3.7686895721198277, -4.7332154502395261,
    4.111321773945793, -5.1634332048395821, 4.4539579832483938,
    -5.5936385356887595, 4.7965981994875238, -6.0238314434235676,
    5.1392424221231163, -6.454011928680476, 5.481890650615143,
    -6.8841799920959144, 5.8245428844236127, -7.3143356343062749,
    6.1671991230085714, -7.744478855947909, 6.5098593658301045,
    -8.174609657657129, 6.8525236123483344, -8.6047280400702082,
    -0.34258808392522866, 0.43035445788055782, -0.68518018127074631,
    0.86069648500604323, -1.0277762914960311, 1.2910260820133979,
    -1.3703764140605992, 1.7213432495395238, -1.7129805484240048,
    2.151647988221284, -2.0555886940458397, 2.5819402986955016,
    -2.3982008503857339, 3.0122201815989613, -2.740817016903355,
    3.4424876375684077, -3.0834371930584092, 3.8727426672405465,
    -3.4260613783106395, 4.3029852712520436, -3.7686895721198277,
    4.7332154502395261, -4.111321773945793, 5.1634332048395821,
    -4.4539579832483938, 5.5936385356887595, -4.7965981994875238,
    6.0238314434235676, -5.1392424221231163, 6.454011928680476,
    -5.481890650615143, 6.8841799920959144, -5.8245428844236127,
    7.3143356343062749, -6.1671991230085714, 7.744478855947909,
    -6.5098593658301045, 8.174609657657129, -6.8525236123483344,
    8.6047280400702082, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T g[258]{ -0.70746101068200651, 0.35521300513038595,
    -1.4149123697315824, 0.71042423457375414, -2.1223540773836476,
    1.0656336882168116, -2.8297861338731107, 1.4208413659462755,
    -3.53720853943487, 1.7760472676488726, -4.244621294303812,
    2.1312513932113393, -4.952024398714812, 2.4864537425204221,
    -5.6594178529027355, 2.8416543154628768, -6.3668016571024371,
    3.1968531119254697, -7.0741758115487583, 3.5520501317949762,
    -7.7815403164765327, 3.9072453749581815, -8.48889517212058,
    4.262438841301881, -9.19624037871571, 4.6176305307128791,
    -9.9035759364967237, 4.9728204430779908, -10.610901845698407,
    5.3280085782840407, -11.318218106555539, 5.6831949362178626,
    -12.025524719302885, 6.0383795167663008, -12.732821684175198,
    6.3935623198162084, -13.440109001407224, 6.7487433452544492,
    -14.147386671233697, 7.103922592967896, 0.70746101068200651,
    -0.35521300513038595, 1.4149123697315824, -0.71042423457375414,
    2.1223540773836476, -1.0656336882168116, 2.8297861338731107,
    -1.4208413659462755, 3.53720853943487, -1.7760472676488726,
    4.244621294303812, -2.1312513932113393, 4.952024398714812,
    -2.4864537425204221, 5.6594178529027355, -2.8416543154628768,
    6.3668016571024371, -3.1968531119254697, 7.0741758115487583,
    -3.5520501317949762, 7.7815403164765327, -3.9072453749581815,
    8.48889517212058, -4.262438841301881, 9.19624037871571, -4.6176305307128791,
    9.9035759364967237, -4.9728204430779908, 10.610901845698407,
    -5.3280085782840407, 11.318218106555539, -5.6831949362178626,
    12.025524719302885, -6.0383795167663008, 12.732821684175198,
    -6.3935623198162084, 13.440109001407224, -6.7487433452544492,
    14.147386671233697, -7.103922592967896, -1.0, -0.0, -0.0, 1.0, 0.0, 0.0,
    0.34920656360027097, 0.074930066385699018, 0.69839994749507994,
    0.14984980015943497, 1.0475801524342094, 0.22475920204730321,
    1.3967471791673947, 0.29965827277535156, 1.7459010284443235,
    0.37454701306958038, 2.0950417010146363, 0.44942542365594257,
    2.4441691976279265, 0.52429350526034346, 2.79328351903374,
    0.59915125860864094, 3.142384665981576, 0.67399868442664534,
    3.491472639220885, 0.74883578344011958, 3.8405474395010719,
    0.82366255637477914, 4.1896090675714932, 0.898479003956292,
    4.538657524181458, 0.97328512691027858, 4.887692810080229,
    1.0480809259623118, 5.2367149260170214, 1.1228664018379175,
    5.5857238727410028, 1.1976415552625737, 5.9347196510012941,
    1.2724063869617108, 6.283702261546968, 1.3471608976607121, 6.63267170512705,
    1.4219050880849133, 6.9816279824905205, 1.4966389589596025,
    -0.34920656360027097, -0.074930066385699018, -0.69839994749507994,
    -0.14984980015943497, -1.0475801524342094, -0.22475920204730321,
    -1.3967471791673947, -0.29965827277535156, -1.7459010284443235,
    -0.37454701306958038, -2.0950417010146363, -0.44942542365594257,
    -2.4441691976279265, -0.52429350526034346, -2.79328351903374,
    -0.59915125860864094, -3.142384665981576, -0.67399868442664534,
    -3.491472639220885, -0.74883578344011958, -3.8405474395010719,
    -0.82366255637477914, -4.1896090675714932, -0.898479003956292,
    -4.538657524181458, -0.97328512691027858, -4.887692810080229,
    -1.0480809259623118, -5.2367149260170214, -1.1228664018379175,
    -5.5857238727410028, -1.1976415552625737, -5.9347196510012941,
    -1.2724063869617108, -6.283702261546968, -1.3471608976607121,
    -6.63267170512705, -1.4219050880849133, -6.9816279824905205,
    -1.4966389589596025, -0.0, -1.0, -0.0, 0.0, 1.0, 0.0, 0.34258808392522866,
    -0.43035445788055782, 0.68518018127074631, -0.86069648500604323,
    1.0277762914960311, -1.2910260820133979, 1.3703764140605992,
    -1.7213432495395238, 1.7129805484240048, -2.151647988221284,
    2.0555886940458397, -2.5819402986955016, 2.3982008503857339,
    -3.0122201815989613, 2.740817016903355, -3.4424876375684077,
    3.0834371930584092, -3.8727426672405465, 3.4260613783106395,
    -4.3029852712520436, 3.7686895721198277, -4.7332154502395261,
    4.111321773945793, -5.1634332048395821, 4.4539579832483938,
    -5.5936385356887595, 4.7965981994875238, -6.0238314434235676,
    5.1392424221231163, -6.454011928680476, 5.481890650615143,
    -6.8841799920959144, 5.8245428844236127, -7.3143356343062749,
    6.1671991230085714, -7.744478855947909, 6.5098593658301045,
    -8.174609657657129, 6.8525236123483344, -8.6047280400702082,
    -0.34258808392522866, 0.43035445788055782, -0.68518018127074631,
    0.86069648500604323, -1.0277762914960311, 1.2910260820133979,
    -1.3703764140605992, 1.7213432495395238, -1.7129805484240048,
    2.151647988221284, -2.0555886940458397, 2.5819402986955016,
    -2.3982008503857339, 3.0122201815989613, -2.740817016903355,
    3.4424876375684077, -3.0834371930584092, 3.8727426672405465,
    -3.4260613783106395, 4.3029852712520436, -3.7686895721198277,
    4.7332154502395261, -4.111321773945793, 5.1634332048395821,
    -4.4539579832483938, 5.5936385356887595, -4.7965981994875238,
    6.0238314434235676, -5.1392424221231163, 6.454011928680476,
    -5.481890650615143, 6.8841799920959144, -5.8245428844236127,
    7.3143356343062749, -6.1671991230085714, 7.744478855947909,
    -6.5098593658301045, 8.174609657657129, -6.8525236123483344,
    8.6047280400702082, -0.0, -0.0, -1.0, 0.0, 0.0, 1.0 };

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

  static const real_T c[32]{ 0.70746101068200651, -0.35521300513038595, 0.0, 0.0,
    -0.34920656360027097, -0.074930066385699018, 0.0, 0.0, -0.34258808392522866,
    0.43035445788055782, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.025, 0.0, 0.0,
    0.0, 0.0, 0.025, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  static const real_T b[16]{ 0.99996947351391963, -1.9978365398021225E-5, 0.0,
    0.0, -3.3626770675199147E-5, 0.99995521109485208, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 1.0 };

  static const real_T h[16]{ 33.484456622589349, -11.58475376880437,
    -20.77520380086132, 0.0, -11.58475376880437, 7.247642711486316,
    4.5947763709116511, 0.0, -20.77520380086132, 4.5947763709116511,
    16.448992358937382, 0.0, 0.0, 0.0, 0.0, 100000.0 };

  static const int32_T b_Mrows[86]{ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
    52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
    71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 141, 142, 143 };

  static const int16_T e[86]{ 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 40, 40, 40, 40, 0, 0 };

  static const int8_T d[8]{ 1, 0, 0, 1, 1, 0, 0, 1 };

  static const int8_T old_mvoff_0[3]{ 40, 0, 0 };

  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  real_T f_0[344];
  real_T k_0[344];
  real_T g_0[258];
  real_T Bu[252];
  real_T b_Mlim[86];
  real_T Bv[84];
  real_T rtb_useq[63];
  real_T b_utarget[60];
  real_T Dv[42];
  real_T rseq[40];
  real_T b_B[32];
  real_T B_est[28];
  real_T vseq[21];
  real_T rtb_A_i[16];
  real_T rtb_Product_d[16];
  real_T rtb_Q[16];
  real_T rtb_Z_b[16];
  real_T D_est[14];
  real_T rtb_R_f_tmp[14];
  real_T rtb_B_e[12];
  real_T rtb_C_k[8];
  real_T rtb_L_g[8];
  real_T rtb_N[8];
  real_T rtb_Product2_k[8];
  real_T rtb_B_em[6];
  real_T rtb_D_m[6];
  real_T rtb_A_p[4];
  real_T rtb_C_fa[4];
  real_T rtb_MemoryX_n[4];
  real_T tmp[4];
  real_T rtb_B_mk[3];
  real_T rtb_D_a[3];
  real_T rtb_D_p[3];
  real_T tmp_0[3];
  real_T tmp_1[3];
  real_T Y[2];
  real_T rtb_Add1_i[2];
  real_T rtb_C;
  real_T rtb_C_o;
  real_T rtb_e_h;
  real_T rtb_e_m;
  int32_T b_utarget_tmp;
  int32_T i;
  int32_T i_0;
  int32_T old_mvoff;
  boolean_T tmp_6[86];
  boolean_T rtb_DataTypeConversionEnable_a1;

  // Delay: '<S290>/delayTheta' incorporates:
  //   Constant: '<S290>/InitialParameters'

  if (localDW->icLoad) {
    localDW->delayTheta_DSTATE[0] = localP->InitialParameters_Value[0];
    localDW->delayTheta_DSTATE[1] = localP->InitialParameters_Value[1];
    localDW->delayTheta_DSTATE[2] = localP->InitialParameters_Value[2];
    localDW->delayTheta_DSTATE[3] = localP->InitialParameters_Value[3];
  }

  // End of Delay: '<S290>/delayTheta'

  // Delay: '<S290>/delayPhi' incorporates:
  //   Constant: '<S290>/InitialPhiMemory'

  if (localDW->icLoad_o) {
    for (i = 0; i < 7; i++) {
      localDW->delayPhi_DSTATE[i] = localP->InitialPhiMemory_Value[i];
    }
  }

  // End of Delay: '<S290>/delayPhi'

  // Delay: '<S290>/delayBufferU' incorporates:
  //   Constant: '<S290>/InitialInputs'

  if (localDW->icLoad_e) {
    localDW->delayBufferU_DSTATE = localP->InitialInputs_Value;
  }

  // End of Delay: '<S290>/delayBufferU'

  // Delay: '<S290>/delayBufferY' incorporates:
  //   Constant: '<S290>/InitialOutputs'

  if (localDW->icLoad_d) {
    localDW->delayBufferY_DSTATE = localP->InitialOutputs_Value;
  }

  // End of Delay: '<S290>/delayBufferY'

  // Delay: '<S290>/delayL' incorporates:
  //   Constant: '<S290>/InitialCovariance'

  if (localDW->icLoad_e2) {
    (void)std::memcpy(&localDW->delayL_DSTATE[0],
                      &localP->InitialCovariance_Value[0], sizeof(real_T) << 4UL);
  }

  // End of Delay: '<S290>/delayL'

  // Update for Delay: '<S290>/delayL' incorporates:
  //   Constant: '<S290>/Forgetting Factor'
  //   DataTypeConversion: '<S290>/DataTypeConversionEnable'
  //   Delay: '<S290>/delayBufferU'
  //   Delay: '<S290>/delayBufferY'
  //   Delay: '<S290>/delayPhi'
  //   Delay: '<S290>/delayTheta'
  //   MATLAB Function: '<S313>/rarx'
  //   UnitDelay: '<S255>/Unit Delay4'
  //   UnitDelay: '<S255>/Unit Delay5'

  rarx(localDW->UnitDelay5_DSTATE, localDW->UnitDelay4_DSTATE[0], rtu_paramEst
       != 0.0, localP->ForgettingFactor_Value, rtb_Add1_i, rtb_B_em, &rtb_e_m,
       localDW->delayTheta_DSTATE, localDW->delayL_DSTATE,
       localDW->delayPhi_DSTATE, &localDW->sf_rarx);

  // MATLAB Function: '<S288>/Model Type Converter'
  ModelTypeConverter(rtb_B_em, rtb_Add1_i, &rtb_e_m, rtb_B_mk, &rtb_C, rtb_D_p);

  // Delay: '<S291>/delayTheta' incorporates:
  //   Constant: '<S291>/InitialParameters'

  if (localDW->icLoad_f) {
    localDW->delayTheta_DSTATE_h[0] = localP->InitialParameters_Value_j[0];
    localDW->delayTheta_DSTATE_h[1] = localP->InitialParameters_Value_j[1];
    localDW->delayTheta_DSTATE_h[2] = localP->InitialParameters_Value_j[2];
    localDW->delayTheta_DSTATE_h[3] = localP->InitialParameters_Value_j[3];
  }

  // End of Delay: '<S291>/delayTheta'

  // Delay: '<S291>/delayPhi' incorporates:
  //   Constant: '<S291>/InitialPhiMemory'

  if (localDW->icLoad_ej) {
    for (i = 0; i < 7; i++) {
      localDW->delayPhi_DSTATE_c[i] = localP->InitialPhiMemory_Value_n[i];
    }
  }

  // End of Delay: '<S291>/delayPhi'

  // Delay: '<S291>/delayBufferU' incorporates:
  //   Constant: '<S291>/InitialInputs'

  if (localDW->icLoad_j) {
    localDW->delayBufferU_DSTATE_l = localP->InitialInputs_Value_f;
  }

  // End of Delay: '<S291>/delayBufferU'

  // Delay: '<S291>/delayBufferY' incorporates:
  //   Constant: '<S291>/InitialOutputs'

  if (localDW->icLoad_i) {
    localDW->delayBufferY_DSTATE_f = localP->InitialOutputs_Value_j;
  }

  // End of Delay: '<S291>/delayBufferY'

  // Delay: '<S291>/delayL' incorporates:
  //   Constant: '<S291>/InitialCovariance'

  if (localDW->icLoad_jb) {
    (void)std::memcpy(&localDW->delayL_DSTATE_i[0],
                      &localP->InitialCovariance_Value_m[0], sizeof(real_T) <<
                      4UL);
  }

  // End of Delay: '<S291>/delayL'

  // Update for Delay: '<S291>/delayL' incorporates:
  //   Constant: '<S291>/Forgetting Factor'
  //   DataTypeConversion: '<S291>/DataTypeConversionEnable'
  //   Delay: '<S291>/delayBufferU'
  //   Delay: '<S291>/delayBufferY'
  //   Delay: '<S291>/delayPhi'
  //   Delay: '<S291>/delayTheta'
  //   MATLAB Function: '<S342>/rarx'
  //   UnitDelay: '<S255>/Unit Delay4'
  //   UnitDelay: '<S255>/Unit Delay5'

  rarx(localDW->UnitDelay5_DSTATE, localDW->UnitDelay4_DSTATE[1], rtu_paramEst
       != 0.0, localP->ForgettingFactor_Value_p, rtb_Add1_i, rtb_B_em, &rtb_e_h,
       localDW->delayTheta_DSTATE_h, localDW->delayL_DSTATE_i,
       localDW->delayPhi_DSTATE_c, &localDW->sf_rarx_o);

  // MATLAB Function: '<S289>/Model Type Converter'
  ModelTypeConverter(rtb_B_em, rtb_Add1_i, &rtb_e_h, rty_u, &rtb_C_o, rtb_D_a);

  // MATLAB Function: '<S255>/MATLAB Function'
  MATLABFunction(rtb_e_m, rtb_B_mk, rtb_C, rtb_D_p, rtb_e_h, rty_u, rtb_C_o,
                 rtb_D_a, rtb_A_p, rtb_B_em, rtb_C_fa, rtb_D_m);

  // Product: '<S255>/Product1' incorporates:
  //   Constant: '<S255>/Constant1'
  //   Product: '<S355>/Product1'

  rtb_Add1_i[0] = localP->Constant1_Value * rtu_y0[0];
  rtb_Add1_i[1] = localP->Constant1_Value * rtu_y0[1];

  // Delay: '<S352>/MemoryX' incorporates:
  //   Constant: '<S352>/X0'
  //   DataTypeConversion: '<S352>/DataTypeConversionReset'

  localDW->icLoad_l = ((static_cast<uint32_T>(localZCE->MemoryX_Reset_ZCE) ==
                        POS_ZCSIG) || localDW->icLoad_l);
  localZCE->MemoryX_Reset_ZCE = 0U;
  if (localDW->icLoad_l) {
    localDW->MemoryX_DSTATE[0] = localP->X0_Value[0];
    localDW->MemoryX_DSTATE[1] = localP->X0_Value[1];
    localDW->MemoryX_DSTATE[2] = localP->X0_Value[2];
    localDW->MemoryX_DSTATE[3] = localP->X0_Value[3];
  }

  rtb_MemoryX_n[0] = localDW->MemoryX_DSTATE[0];
  rtb_MemoryX_n[1] = localDW->MemoryX_DSTATE[1];
  rtb_MemoryX_n[2] = localDW->MemoryX_DSTATE[2];
  rtb_MemoryX_n[3] = localDW->MemoryX_DSTATE[3];

  // MATLAB Function: '<S285>/FixedHorizonOptimizer' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'
  //   Constant: '<S255>/Constant9'
  //   Constant: '<S256>/Constant'
  //   Delay: '<S352>/MemoryX'
  //   Product: '<S355>/Product1'
  //   RandomNumber: '<S4>/Measurement Noise'
  //   Sum: '<S256>/Sum2'
  //   UnitDelay: '<S257>/last_mv'
  //
  // MATLAB Function 'Adaptive MPC Controller/MPC/optimizer/FixedHorizonOptimizer': '<S286>:1' 
  // '<S286>:1:18' coder.extrinsic('mpcblock_optimizer_double_mex');
  // '<S286>:1:19' coder.extrinsic('mpcblock_optimizer_single_mex');
  // '<S286>:1:20' coder.extrinsic('mpcblock_refmd_double_mex');
  // '<S286>:1:21' coder.extrinsic('mpcblock_refmd_single_mex');
  //  Inputs (in BlockDataType except iA)
  //    xk:         current state (either x[k|k-1] from built-in KF or external x[k|k]) 
  // '<S286>:1:25' xk = convertDataType(xk0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    old_u:      last mv (calculated by MPC)
  // '<S286>:1:27' old_u = convertDataType(old_u0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    ym:         current measured output (used only with built-in KF)
  // '<S286>:1:29' ym = convertDataType(ym0,isDouble);
  //    ref:        output reference
  // '<S286>:1:31' ref = convertDataType(ref0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    md:         measured disturbance
  // '<S286>:1:33' md = convertDataType(md0,isDouble);
  //    umin:       run-time MV bound
  // '<S286>:1:35' umin = convertDataType(umin0,isDouble);
  //    umax:       run-time MV bound
  // '<S286>:1:37' umax = convertDataType(umax0,isDouble);
  //    ymin:       run-time OV bound
  // '<S286>:1:39' ymin = convertDataType(ymin0,isDouble);
  //    ymax:       run-time OV bound
  // '<S286>:1:41' ymax = convertDataType(ymax0,isDouble);
  //    E:          run-time mixed constraints
  // '<S286>:1:43' E = convertDataType(E0,isDouble);
  //    F:          run-time mixed constraints
  // '<S286>:1:45' F = convertDataType(F0,isDouble);
  //    G:          run-time mixed constraints
  // '<S286>:1:47' G = convertDataType(G0,isDouble);
  //    S:          run-time mixed constraints
  // '<S286>:1:49' S = convertDataType(S0,isDouble);
  //    switch_in:  if it matches "enable_value", MPC is active in control
  // '<S286>:1:51' switch_in = int32(switch_in0);
  //    ext_mv:     external last mv (actual)
  // '<S286>:1:53' ext_mv = convertDataType(ext_mv0,isDouble);
  //    MVtarget:   MV reference
  // '<S286>:1:55' MVtarget = convertDataType(MVtarget0,isDouble);
  //    ywt:        run-time OV weights
  // '<S286>:1:57' ywt = convertDataType(ywt0,isDouble);
  //    uwt:        run-time MV weights
  // '<S286>:1:59' uwt = convertDataType(uwt0,isDouble);
  //    duwt:       run-time DMV weights
  // '<S286>:1:61' duwt = convertDataType(duwt0,isDouble);
  //    rhoeps:     run-time Slack weights
  // '<S286>:1:63' ewt = convertDataType(ewt0,isDouble);
  //    a:          run-time A (must be in DT)
  // '<S286>:1:65' a = convertDataType(a0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    b:          run-time B (must be in DT)
  // '<S286>:1:67' b = convertDataType(b0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    c:          run-time C (must be in DT)
  // '<S286>:1:69' c = convertDataType(c0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    d:          run-time D (must be in DT)
  // '<S286>:1:71' d = convertDataType(d0,isDouble);
  //    U:          run-time nominal value
  // '<S286>:1:73' U = convertDataType(U0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    Y:          run-time nominal value
  // '<S286>:1:75' Y = convertDataType(Y0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    X:          run-time nominal value
  // '<S286>:1:77' X = convertDataType(X0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    DX:         run-time nominal value
  // '<S286>:1:79' DX = convertDataType(DX0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
  //    Pk:         covariance P[k|k-1] (used only with built-in KF)
  // '<S286>:1:81' Pk = convertDataType(Pk0,isDouble);
  // '<S286>:1:317' if isDouble
  //  convert an input signal to double precision when necessary
  // '<S286>:1:319' if isa(u,'double')
  // '<S286>:1:320' y = u;
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
  // '<S286>:1:115' isSimulation = coder.target('Sfun') && ~coder.target('RtwForRapid') && ~coder.target('RtwForSim'); 
  // '<S286>:1:116' isAdaptive = ~isLTV;
  //  isLTV=true forces isAdaptive=false
  // '<S286>:1:117' ZERO = zeros('like',ref);
  // '<S286>:1:118' ONE = ones('like',ref);
  // '<S286>:1:119' hasMD = nv>int32(1);
  //  Pre-allocate all the MEX block outputs for the simulation mode
  // '<S286>:1:123' if isSimulation
  //  Model update
  // '<S286>:1:137' nym = int32(length(myindex));
  // '<S286>:1:138' ai=zeros(nxp,nxp,'like',ref);
  // '<S286>:1:139' bi=zeros(nxp,nup,'like',ref);
  // '<S286>:1:140' ci=zeros(ny,nxp,'like',ref);
  // '<S286>:1:141' di=zeros(ny,nup,'like',ref);
  // '<S286>:1:143' ai(:,:)=a(:,:,1);
  // '<S286>:1:144' bi(:,:)=b(:,:,1);
  // '<S286>:1:145' ci(:,:)=c(:,:,1);
  // '<S286>:1:146' di(:,:)=d(:,:,1);
  //  Allocate matrices. Must allocate 3D matrix also in Adaptive case,
  //  otherwise EML code does not compile.
  // '<S286>:1:150' Bu=zeros(nx,nu,p+1,'like',ref);
  (void)std::memset(&Bu[0], 0, 252U * sizeof(real_T));

  // '<S286>:1:151' Bv=zeros(nx,nv,p+1,'like',ref);
  (void)std::memset(&Bv[0], 0, 84U * sizeof(real_T));

  // '<S286>:1:152' Dv=zeros(ny,nv,p+1,'like',ref);
  (void)std::memset(&Dv[0], 0, 42U * sizeof(real_T));

  // '<S286>:1:153' Dvm=zeros(nym,nv,p+1,'like',ref);
  // '<S286>:1:154' Cm=zeros(nym,nx,p+1,'like',ref);
  // '<S286>:1:155' [A(:,:,1),C(:,:,1),Bu(:,:,1),Bv(:,:,1),Cm(:,:,1),Dv(:,:,1),Dvm(:,:,1),Qk,Rk,Nk] = mpc_plantupdate(... 
  // '<S286>:1:156'     ai,bi,ci,di,A(:,:,1),B(:,:,1),C(:,:,1),D(:,:,1),mvindex,mdindex,unindex,nxp,nup,ny,nu,nv,nxid, ... 
  // '<S286>:1:157'     myindex,Uscale,Yscale,Cid,Did);
  (void)std::memcpy(&rtb_A_i[0], &b[0], sizeof(real_T) << 4UL);
  (void)std::memcpy(&b_B[0], &c[0], sizeof(real_T) << 5UL);
  for (i = 0; i < 8; i++) {
    rtb_C_k[i] = static_cast<real_T>(d[i]);
  }

  rtb_C_k[0] = rtb_C_fa[0];
  rtb_A_i[0] = rtb_A_p[0];
  rtb_C_k[1] = rtb_C_fa[1];
  rtb_A_i[1] = rtb_A_p[1];
  rtb_C_k[2] = rtb_C_fa[2];
  rtb_A_i[4] = rtb_A_p[2];
  rtb_C_k[3] = rtb_C_fa[3];
  rtb_A_i[5] = rtb_A_p[3];
  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    b_B[i] = rtb_B_em[i_0];
    b_B[i + 1] = rtb_B_em[i_0 + 1];
    Bu[i] = b_B[i];
    Bu[i + 1] = b_B[i + 1];
    Bu[i + 2] = b_B[i + 2];
    Bu[i + 3] = b_B[i + 3];
    i += 4;
    i_0 += 2;
  }

  Bv[2] = b_B[14];
  Bv[3] = b_B[15];
  Dv[0] = 0.0;
  Dv[1] = 0.0;

  // '<S286>:1:158' if isLTV
  //  Offset update together with Mlim, utarget, Bv and Dv values
  // '<S286>:1:174' [Mlim, utarget, uoff, voff, yoff, myoff, xoff, Bv, Dv] = ... 
  // '<S286>:1:175'     mpc_updateFromNominal(isAdaptive,isQP,Mlim,Mrows,...
  // '<S286>:1:176'        U,Uscale,uoff,mvindex,voff,mdindex,utarget,nu,nv-1,... 
  // '<S286>:1:177'        Y,Yscale,yoff,myoff,myindex,ny,...
  // '<S286>:1:178'        X,xoff,nxp,DX,A,Bu,Bv,C,Dv,nCC);
  for (i = 0; i < 86; i++) {
    b_Mlim[i] = static_cast<real_T>(e[i]);
  }

  (void)std::memset(&b_utarget[0], 0, 60U * sizeof(real_T));
  rtb_D_p[0] = rtu_u0[0];
  rtb_D_p[1] = rtu_u0[1];
  rtb_D_p[2] = rtu_u0[2];
  Y[0] = rtu_y0[0];
  Y[1] = rtu_y0[1];
  rtb_B_mk[0] = rtu_u0[0];
  rtb_B_mk[1] = rtu_u0[1];
  rtb_B_mk[2] = rtu_u0[2];
  for (int32_T b_ct{0}; b_ct < 86; b_ct++) {
    rtb_e_m = b_Mlim[b_ct];
    i = b_Mrows[b_ct];
    if (i <= 40) {
      rtb_e_m += 100.0 - Y[(i - (((i - 1) / static_cast<int32_T>(ny_i)) << 1UL))
        - 1];
    } else if (i <= 80) {
      rtb_e_m -= 100.0 - Y[(i - (((i - 41) >> 1UL) << 1UL)) - 41];
    } else if (i <= 140) {
      i = (i - div_nde_s32_floor(i - 81, static_cast<int32_T>(nu_h)) *
           static_cast<int32_T>(nu_h)) - 81;
      rtb_e_m += static_cast<real_T>(old_mvoff_0[i]) - rtb_D_p[i];
    } else {
      i = (i - div_nde_s32_floor(i - 141, static_cast<int32_T>(nu_h)) *
           static_cast<int32_T>(nu_h)) - 141;
      rtb_e_m -= static_cast<real_T>(old_mvoff_0[i]) - rtb_D_p[i];
    }

    b_Mlim[b_ct] = rtb_e_m;
  }

  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    old_mvoff = static_cast<int32_T>(old_mvoff_0[b_ct]);
    rtb_e_m = rtb_D_p[b_ct];
    i = 0;
    for (i_0 = 0; i_0 < 20; i_0++) {
      b_utarget_tmp = i + b_ct;
      b_utarget[b_utarget_tmp] = (b_utarget[b_utarget_tmp] + static_cast<real_T>
        (old_mvoff)) - rtb_e_m;
      i += 3;
    }
  }

  Bv[0] = localP->Constant9_Value[0];
  Bv[1] = localP->Constant9_Value[1];

  //  Remove last u offset
  // '<S286>:1:181' old_u = old_u - uoff;
  //  Get reference and MD signals -- accounting for previewing
  // '<S286>:1:184' if isSimulation
  // '<S286>:1:190' else
  //  When doing code generation, use M code directly
  // '<S286>:1:192' [rseq, vseq, v] = mpcblock_refmd(ref,md,nv,ny,p,yoff,voff,no_md,no_ref,openloopflag, RYscale, RMDscale); 
  for (int32_T b_ct{0}; b_ct < 21; b_ct++) {
    vseq[b_ct] = 1.0;
  }

  i = 0;
  for (int32_T b_ct{0}; b_ct < 20; b_ct++) {
    rseq[i] = rtu_r[0] - Y[0];
    rseq[i + 1] = rtu_r[1] - Y[1];
    i += 2;
  }

  //  External MV override.
  //  NOTE: old_u and ext_mv input signals are dimensionless and offset-free.
  // '<S286>:1:197' if no_mv
  //  no external mv: old_u is the optimal u[k-1] from last_mv
  // '<S286>:1:199' delmv = zeros(nu,1,'like',ref);
  //  Obtain x[k|k]
  // '<S286>:1:208' xk = xk - xoff;
  //  Remove offset
  // '<S286>:1:209' if CustomEstimation
  //  Input is x(k|k)
  // '<S286>:1:211' xest = xk;
  //  Real-time MV target override
  //  Note: utargetValue is a vector length p*nu.
  // '<S286>:1:231' if no_uref
  //  no external utarget
  // '<S286>:1:233' utargetValue = utarget;
  //  Real-time custom constraint override (scaled E/F/S)
  // '<S286>:1:241' if ~no_cc
  // '<S286>:1:250' return_sequence = return_mvseq || return_xseq || return_ovseq; 
  // '<S286>:1:251' if isSimulation
  // '<S286>:1:279' else
  //  When doing code generation, use M code directly
  // '<S286>:1:281' [u, cost, useq, status, iAout] = mpcblock_optimizer(...
  // '<S286>:1:282'             rseq, vseq, umin, umax, ymin, ymax, switch_in, xest, old_u, iA, ... 
  // '<S286>:1:283'             isQP, nu, ny, degrees, Hinv, Kx, Ku1, Kut, Kr, Kv, Mlim, ... 
  // '<S286>:1:284'             Mx, Mu1, Mv, utargetValue, p, uoff, voff, yoff, ... 
  // '<S286>:1:285'             false, CustomSolverCodeGen, UseSuboptimalSolution, ... 
  // '<S286>:1:286'             UseActiveSetSolver, ASOptions, IPOptions, MIQPOptions, nxQP, openloopflag, ... 
  // '<S286>:1:287'             no_umin, no_umax, no_ymin, no_ymax, no_cc, switch_inport, ... 
  // '<S286>:1:288'             no_switch, enable_value, return_cost, H, return_sequence, Linv, Ac, ... 
  // '<S286>:1:289'             ywt, uwt, duwt, ewt, no_ywt, no_uwt, no_duwt, no_rhoeps,... 
  // '<S286>:1:290'             Wy, Wdu, Jm, SuJm, Su1, Sx, Hv, Wu, I1, ...
  // '<S286>:1:291'             isAdaptive, isLTV, A, Bu, Bv, C, Dv, ...
  // '<S286>:1:292'             Mrows, nCC, Ecc, Fcc, Scc, Gcc, RYscale, RMVscale, m, isHyb, Mdis, Ndis, Vdis, numdis, maxdis); 
  tmp[0] = (localDW->MemoryX_DSTATE[0] + rtb_Add1_i[0]) - rtb_Add1_i[0];
  tmp[1] = (localDW->MemoryX_DSTATE[1] + rtb_Add1_i[1]) - rtb_Add1_i[1];
  tmp[2] = localP->Constant_Value[0] + localDW->MemoryX_DSTATE[2];
  tmp[3] = localP->Constant_Value[1] + localDW->MemoryX_DSTATE[3];
  rtb_D_a[0] = localDW->last_mv_DSTATE[0] - rtu_u0[0];
  rtb_D_a[1] = localDW->last_mv_DSTATE[1] - rtu_u0[1];
  rtb_D_a[2] = localDW->last_mv_DSTATE[2] - rtu_u0[2];
  (void)std::memset(&localDW->dv[0], 0, 1806U * sizeof(real_T));
  Y[0] = 0.018316915599999997;
  Y[1] = 0.018316915599999997;
  tmp_0[0] = 0.54594344475769718;
  tmp_1[0] = 0.0;
  tmp_0[1] = 0.54594344475769718;
  tmp_1[1] = 0.0;
  tmp_0[2] = 0.54594344475769718;
  tmp_1[2] = 0.0;

  // Memory: '<S257>/Memory'
  (void)std::memcpy(&tmp_6[0], &localDW->Memory_PreviousInput[0], 86U * sizeof
                    (boolean_T));

  // MATLAB Function: '<S285>/FixedHorizonOptimizer'
  (void)std::memcpy(&f_0[0], &f[0], 344U * sizeof(real_T));
  (void)std::memcpy(&g_0[0], &g[0], 258U * sizeof(real_T));
  (void)std::memcpy(&rtb_Q[0], &h[0], sizeof(real_T) << 4UL);
  (void)std::memcpy(&k_0[0], &k[0], 344U * sizeof(real_T));

  // Update for Memory: '<S257>/Memory' incorporates:
  //   MATLAB Function: '<S285>/FixedHorizonOptimizer'

  mpcblock_optimizer_l(rseq, vseq, tmp, rtb_D_a, tmp_6, b_Mlim, f_0, g_0,
                       localDW->dv, b_utarget, rtb_B_mk, rtb_Q, k_0, Y, tmp_0, l,
                       tmp_1, n, rtb_A_i, Bu, Bv, rtb_C_k, Dv, b_Mrows, rtb_D_p,
                       rtb_useq, &rtb_e_m, localDW->Memory_PreviousInput,
                       localDW);

  // Delay: '<S352>/MemoryP' incorporates:
  //   Constant: '<S352>/P0'
  //   DataTypeConversion: '<S352>/DataTypeConversionReset'

  // '<S286>:1:295' if return_xseq || return_ovseq
  // '<S286>:1:297' else
  // '<S286>:1:298' yseq = zeros(p+1,ny,'like',rseq);
  // '<S286>:1:299' xseq = zeros(p+1,nxQP,'like',rseq);
  // '<S286>:1:302' if CustomEstimation
  // '<S286>:1:303' xk1 = zeros(nx,1,'like',ref);
  // '<S286>:1:304' Pk1 = Pk;
  // '<S286>:1:311' xk1 = xk1 + xoff;
  //  Updated state must include offset
  //  return xest in original value
  // '<S286>:1:314' xest = xest + xoff;
  localDW->icLoad_h = ((static_cast<uint32_T>(localZCE->MemoryP_Reset_ZCE) ==
                        POS_ZCSIG) || localDW->icLoad_h);
  localZCE->MemoryP_Reset_ZCE = 0U;
  if (localDW->icLoad_h) {
    (void)std::memcpy(&localDW->MemoryP_DSTATE[0], &localP->P0_Value[0], sizeof
                      (real_T) << 4UL);
  }

  // DataTypeConversion: '<S352>/DataTypeConversionEnable' incorporates:
  //   Constant: '<S4>/Constant'

  rtb_DataTypeConversionEnable_a1 = (localP->Constant_Value_g != 0.0);

  // MATLAB Function: '<S256>/MATLAB Function' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'

  // MATLAB Function 'SupervisoryController/State2.controlLaw.AMPC2/State Estimator (KF)/MATLAB Function': '<S353>:1' 
  // '<S353>:1:2' [A, B, C, D, Q, R, N] = stateEstimator(Ap, Bp, Cp, Dp, Aod2, Bod2, Cod2, Dod2, Dmn2); 
  // 'stateEstimator:3' no = size(Cp,1);
  //  n_outputs
  // 'stateEstimator:4' ni = 3;
  //  n_inputs
  // 'stateEstimator:5' nsp = size(Ap,1);
  //  n_plant_states
  // 'stateEstimator:6' ns = nsp + no;
  //  n_states = n_plant_states + n_outputs
  // 'stateEstimator:8' A = zeros(ns);
  //  n_states x n_states
  // 'stateEstimator:9' B = zeros(ns,ni);
  //  n_states  x n_inputs
  // 'stateEstimator:10' C = zeros(no,ns);
  //  n_outputs x n_states
  // 'stateEstimator:11' D = zeros(no,ni);
  //  n_outputs x n_inputs
  // 'stateEstimator:12' Q = zeros(ns,ns);
  //  n_states  x n_states
  // 'stateEstimator:13' G = eye(ns);
  //  n_states  x n_states
  // 'stateEstimator:14' R = zeros(no,no);
  //  n_outputs x n_outputs
  // 'stateEstimator:15' N = zeros(ns,no);
  //  n_states  x n_outputs
  // 'stateEstimator:16' H = zeros(no,ns);
  //  n_outputs x n_states
  //  combine plant and output disturbance model
  // 'stateEstimator:19' A = [Ap, zeros(nsp,no);
  // 'stateEstimator:20'      zeros(no,nsp), Aod];
  rtb_A_i[0] = rtb_A_p[0];
  rtb_A_i[8] = 0.0;
  rtb_A_i[2] = 0.0;
  rtb_A_i[10] = rtP->Aod2[0];
  rtb_A_i[1] = rtb_A_p[1];
  rtb_A_i[9] = 0.0;
  rtb_A_i[3] = 0.0;
  rtb_A_i[11] = rtP->Aod2[1];
  rtb_A_i[4] = rtb_A_p[2];
  rtb_A_i[12] = 0.0;
  rtb_A_i[6] = 0.0;
  rtb_A_i[14] = rtP->Aod2[2];
  rtb_A_i[5] = rtb_A_p[3];
  rtb_A_i[13] = 0.0;
  rtb_A_i[7] = 0.0;
  rtb_A_i[15] = rtP->Aod2[3];

  // 'stateEstimator:21' B = [Bp;
  // 'stateEstimator:22'      zeros(no,ni)];
  // 'stateEstimator:23' C = [Cp Cod];
  rtb_C_k[0] = rtb_C_fa[0];
  rtb_C_k[4] = rtP->Cod2[0];
  rtb_C_k[1] = rtb_C_fa[1];
  rtb_C_k[5] = rtP->Cod2[1];
  rtb_C_k[2] = rtb_C_fa[2];
  rtb_C_k[6] = rtP->Cod2[2];
  rtb_C_k[3] = rtb_C_fa[3];
  rtb_C_k[7] = rtP->Cod2[3];

  // 'stateEstimator:24' D = Dp;
  // 'stateEstimator:26' B_est = [ [Bp; zeros(no,ni)] [zeros(size(Bp,1), size(Cp,1)); Bod] [zeros(size(Bp,1) + size(Bod,1), size(Cp,1))] ]; 
  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    B_est[i] = rtb_B_em[i_0];
    B_est[i + 2] = 0.0;
    B_est[i + 1] = rtb_B_em[i_0 + 1];
    B_est[i + 3] = 0.0;
    i += 4;
    i_0 += 2;
  }

  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 2; b_ct++) {
    B_est[i + 12] = 0.0;
    B_est[i + 14] = rtP->Bod2[i_0];
    B_est[i + 13] = 0.0;
    B_est[i + 15] = rtP->Bod2[i_0 + 1];
    B_est[i + 20] = 0.0;
    B_est[i + 21] = 0.0;
    B_est[i + 22] = 0.0;
    B_est[i + 23] = 0.0;
    i += 4;
    i_0 += 2;
  }

  // 'stateEstimator:27' D_est = [Dp Dod Dn];
  for (i = 0; i < 6; i++) {
    D_est[i] = rtb_D_m[i];
  }

  // 'stateEstimator:28' Q = B_est * B_est';
  for (i = 0; i < 4; i++) {
    D_est[i + 6] = rtP->Dod2[i];
    D_est[i + 10] = rtP->Dmn2[i];
    i_0 = 0;
    for (int32_T b_ct{0}; b_ct < 4; b_ct++) {
      int32_T tmp_2;
      b_utarget_tmp = i_0 + i;
      rtb_Q[b_utarget_tmp] = 0.0;
      tmp_2 = 0;
      for (old_mvoff = 0; old_mvoff < 7; old_mvoff++) {
        rtb_Q[b_utarget_tmp] += B_est[tmp_2 + i] * B_est[tmp_2 + b_ct];
        tmp_2 += 4;
      }

      i_0 += 4;
    }
  }

  // 'stateEstimator:29' R = D_est * D_est';
  i = 0;
  for (i_0 = 0; i_0 < 7; i_0++) {
    rtb_R_f_tmp[i_0] = D_est[i];
    rtb_R_f_tmp[i_0 + 7] = D_est[i + 1];
    i += 2;
  }

  // 'stateEstimator:30' N = B_est * D_est';
  for (i = 0; i < 2; i++) {
    for (i_0 = 0; i_0 < 2; i_0++) {
      old_mvoff = (i << 1UL) + i_0;
      rtb_C_fa[old_mvoff] = 0.0;
      for (int32_T b_ct{0}; b_ct < 7; b_ct++) {
        rtb_C_fa[old_mvoff] += D_est[(b_ct << 1UL) + i_0] * rtb_R_f_tmp[7 * i +
          b_ct];
      }
    }

    for (i_0 = 0; i_0 < 4; i_0++) {
      old_mvoff = (i << 2UL) + i_0;
      rtb_N[old_mvoff] = 0.0;
      for (int32_T b_ct{0}; b_ct < 7; b_ct++) {
        rtb_N[old_mvoff] += B_est[(b_ct << 2UL) + i_0] * rtb_R_f_tmp[7 * i +
          b_ct];
      }
    }
  }

  // Outputs for Atomic SubSystem: '<S352>/ScalarExpansionR'
  //  [k,L,~,Mx,~,My] = kalman(ss(A,[B G],C,[D H],dt), Q, R, N);
  //  [k,L,~,Mx,~,My] = kalman(ss(A,B_est,C,D_est,dt), Q, R, N);
  //  xhat = A*xhat_prev + B*u + L*(y - C*xhat_prev);
  //  yhat = C*xhat + D*u;
  ScalarExpansionR(rtb_C_fa, rtb_A_p);

  // End of Outputs for SubSystem: '<S352>/ScalarExpansionR'

  // Outputs for Atomic SubSystem: '<S352>/ScalarExpansionQ'
  ScalarExpansionQ(rtb_Q, rtb_Z_b);

  // End of Outputs for SubSystem: '<S352>/ScalarExpansionQ'

  // Outputs for Atomic SubSystem: '<S352>/ReducedQRN'
  // Constant: '<S352>/G' incorporates:
  //   Constant: '<S352>/H'

  ReducedQRN(localP->G_Value, localP->H_Value, rtb_Z_b, rtb_A_p, rtb_N,
             rtb_Product_d, rtb_C_fa, rtb_Product2_k);

  // End of Outputs for SubSystem: '<S352>/ReducedQRN'

  // Outputs for Atomic SubSystem: '<S352>/CalculatePL'
  // DataTypeConversion: '<S352>/DataTypeConversionEnable' incorporates:
  //   Constant: '<S4>/Constant'
  //   Delay: '<S352>/MemoryP'

  CalculatePL(rtb_A_i, rtb_C_k, rtb_Product_d, rtb_C_fa, rtb_Product2_k,
              localP->Constant_Value_g != 0.0, localDW->MemoryP_DSTATE, rtb_N,
              rtb_L_g, rtb_Z_b, rtb_Q);

  // End of Outputs for SubSystem: '<S352>/CalculatePL'

  // MATLAB Function: '<S393>/SqrtUsedFcn' incorporates:
  //   Constant: '<S393>/isSqrtUsed'

  SqrtUsedFcn(rtb_Z_b, localP->isSqrtUsed_Value, rtb_Product_d);

  // RandomNumber: '<S4>/Measurement Noise'
  rty_u[0] = localDW->NextOutput[0];

  // Sum: '<S4>/Sum' incorporates:
  //   Gain: '<S257>/u_scale'
  //   Product: '<S4>/Product'

  rty_u[0] = localP->u_scale_Gain[0] * rtb_D_p[0] + rtu_paramEst * rty_u[0];

  // Sum: '<S256>/Sum1' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'

  rtb_B_mk[0] = rty_u[0] - rtu_u0[0];

  // RandomNumber: '<S4>/Measurement Noise'
  rty_u[1] = localDW->NextOutput[1];

  // Sum: '<S4>/Sum' incorporates:
  //   Gain: '<S257>/u_scale'
  //   Product: '<S4>/Product'

  rty_u[1] = localP->u_scale_Gain[1] * rtb_D_p[1] + rtu_paramEst * rty_u[1];

  // Sum: '<S256>/Sum1' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'

  rtb_B_mk[1] = rty_u[1] - rtu_u0[1];

  // RandomNumber: '<S4>/Measurement Noise'
  rty_u[2] = localDW->NextOutput[2];

  // Sum: '<S4>/Sum' incorporates:
  //   Gain: '<S257>/u_scale'
  //   Product: '<S4>/Product'

  rty_u[2] = localP->u_scale_Gain[2] * rtb_D_p[2] + rtu_paramEst * rty_u[2];

  // Sum: '<S256>/Sum1' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'
  //   Sum: '<S255>/Add4'

  rtb_e_m = rty_u[2] - rtu_u0[2];
  rtb_B_mk[2] = rtb_e_m;

  // Sum: '<S256>/Sum6' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'
  //   Sum: '<S255>/Add2'

  rtb_C = rtu_y[0] - rtu_y0[0];
  rtb_Add1_i[0] = rtb_C;
  rtb_e_h = rtu_y[1] - rtu_y0[1];
  rtb_Add1_i[1] = rtb_e_h;

  // Outputs for Enabled SubSystem: '<S371>/MeasurementUpdate'
  // BusCreator: '<S255>/Bus Creator2'
  MeasurementUpdate(rtb_DataTypeConversionEnable_a1, rtb_L_g, rtb_Add1_i,
                    rtb_C_k, rtb_MemoryX_n, rtb_D_m, rtb_B_mk, localDW->Product3,
                    &localDW->MeasurementUpdate_a, &localP->MeasurementUpdate_a);

  // End of Outputs for SubSystem: '<S371>/MeasurementUpdate'

  // MATLAB Function: '<S256>/MATLAB Function' incorporates:
  //   BusCreator: '<S255>/Bus Creator2'

  i = 0;
  i_0 = 0;
  for (int32_T b_ct{0}; b_ct < 3; b_ct++) {
    rtb_B_e[i] = rtb_B_em[i_0];
    rtb_B_e[i + 2] = 0.0;
    rtb_B_e[i + 1] = rtb_B_em[i_0 + 1];
    rtb_B_e[i + 3] = 0.0;
    i += 4;
    i_0 += 2;
  }

  for (i = 0; i <= 2; i += 2) {
    // Product: '<S371>/B[k]*u[k]' incorporates:
    //   Sum: '<S256>/Sum1'

    tmp_5 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_A_p[i], tmp_5);
    tmp_3 = _mm_loadu_pd(&rtb_B_e[i]);
    tmp_4 = _mm_loadu_pd(&rtb_A_p[i]);
    (void)_mm_storeu_pd(&rtb_A_p[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_B_mk[0])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_B_e[i + 4]);
    tmp_4 = _mm_loadu_pd(&rtb_A_p[i]);
    (void)_mm_storeu_pd(&rtb_A_p[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_B_mk[1])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_B_e[i + 8]);
    tmp_4 = _mm_loadu_pd(&rtb_A_p[i]);
    (void)_mm_storeu_pd(&rtb_A_p[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_e_m)), tmp_4));

    // Product: '<S371>/A[k]*xhat[k|k-1]' incorporates:
    //   Product: '<S371>/B[k]*u[k]'

    (void)_mm_storeu_pd(&rtb_C_fa[i], tmp_5);
    tmp_5 = _mm_loadu_pd(&rtb_A_i[i]);
    tmp_3 = _mm_loadu_pd(&rtb_C_fa[i]);
    (void)_mm_storeu_pd(&rtb_C_fa[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_n[0])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_i[i + 4]);
    tmp_3 = _mm_loadu_pd(&rtb_C_fa[i]);
    (void)_mm_storeu_pd(&rtb_C_fa[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_n[1])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_i[i + 8]);
    tmp_3 = _mm_loadu_pd(&rtb_C_fa[i]);
    (void)_mm_storeu_pd(&rtb_C_fa[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_n[2])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_A_i[i + 12]);
    tmp_3 = _mm_loadu_pd(&rtb_C_fa[i]);
    (void)_mm_storeu_pd(&rtb_C_fa[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_MemoryX_n[3])), tmp_3));

    // Product: '<S371>/B[k]*u[k]'
    tmp_5 = _mm_loadu_pd(&rtb_A_p[i]);

    // Product: '<S371>/A[k]*xhat[k|k-1]' incorporates:
    //   Product: '<S371>/B[k]*u[k]'

    tmp_3 = _mm_loadu_pd(&rtb_C_fa[i]);

    // Sum: '<S371>/Add' incorporates:
    //   Product: '<S371>/B[k]*u[k]'

    tmp_4 = _mm_loadu_pd(&localDW->Product3[i]);

    // Update for Delay: '<S352>/MemoryX' incorporates:
    //   Product: '<S371>/B[k]*u[k]'

    (void)_mm_storeu_pd(&localDW->MemoryX_DSTATE[i], _mm_add_pd(_mm_add_pd(tmp_5,
      tmp_3), tmp_4));
  }

  for (i = 0; i <= 0; i += 2) {
    // Product: '<S355>/Product'
    tmp_5 = _mm_set1_pd(0.0);
    (void)_mm_storeu_pd(&rtb_Add1_i[i], tmp_5);
    tmp_3 = _mm_loadu_pd(&rtb_C_k[i]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_i[i]);
    (void)_mm_storeu_pd(&rtb_Add1_i[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_MemoryX_n[0])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_k[i + 2]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_i[i]);
    (void)_mm_storeu_pd(&rtb_Add1_i[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_MemoryX_n[1])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_k[i + 4]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_i[i]);
    (void)_mm_storeu_pd(&rtb_Add1_i[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_MemoryX_n[2])), tmp_4));
    tmp_3 = _mm_loadu_pd(&rtb_C_k[i + 6]);
    tmp_4 = _mm_loadu_pd(&rtb_Add1_i[i]);
    (void)_mm_storeu_pd(&rtb_Add1_i[i], _mm_add_pd(_mm_mul_pd(tmp_3, _mm_set1_pd
      (rtb_MemoryX_n[3])), tmp_4));

    // Product: '<S355>/Product1' incorporates:
    //   BusCreator: '<S255>/Bus Creator2'
    //   Product: '<S355>/Product'
    //   Sum: '<S256>/Sum1'

    (void)_mm_storeu_pd(&Y[i], tmp_5);
    tmp_5 = _mm_loadu_pd(&rtb_D_m[i]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_B_mk[0])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_D_m[i + 2]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd
      (rtb_B_mk[1])), tmp_3));
    tmp_5 = _mm_loadu_pd(&rtb_D_m[i + 4]);
    tmp_3 = _mm_loadu_pd(&Y[i]);
    (void)_mm_storeu_pd(&Y[i], _mm_add_pd(_mm_mul_pd(tmp_5, _mm_set1_pd(rtb_e_m)),
      tmp_3));

    // Product: '<S355>/Product'
    tmp_5 = _mm_loadu_pd(&rtb_Add1_i[i]);

    // Product: '<S355>/Product1' incorporates:
    //   Product: '<S355>/Product'

    tmp_3 = _mm_loadu_pd(&Y[i]);

    // BusCreator: '<S255>/Bus Creator2' incorporates:
    //   Product: '<S355>/Product'

    tmp_4 = _mm_loadu_pd(&rtu_y0[i]);

    // Sum: '<S256>/Sum3' incorporates:
    //   Product: '<S355>/Product'

    (void)_mm_storeu_pd(&rty_yhat[i], _mm_add_pd(_mm_add_pd(tmp_5, tmp_3), tmp_4));
  }

  // Update for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[0] = rtb_D_p[0];

  // Update for UnitDelay: '<S255>/Unit Delay5' incorporates:
  //   Sum: '<S255>/Add4'

  localDW->UnitDelay5_DSTATE[0] = rty_u[0] - rtu_u0[0];

  // Update for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[1] = rtb_D_p[1];

  // Update for UnitDelay: '<S255>/Unit Delay5' incorporates:
  //   Sum: '<S255>/Add4'

  localDW->UnitDelay5_DSTATE[1] = rty_u[1] - rtu_u0[1];

  // Update for UnitDelay: '<S257>/last_mv'
  localDW->last_mv_DSTATE[2] = rtb_D_p[2];

  // Update for UnitDelay: '<S255>/Unit Delay5'
  localDW->UnitDelay5_DSTATE[2] = rtb_e_m;

  // Update for UnitDelay: '<S255>/Unit Delay4'
  localDW->UnitDelay4_DSTATE[0] = rtb_C;
  localDW->UnitDelay4_DSTATE[1] = rtb_e_h;

  // Update for Delay: '<S290>/delayTheta'
  localDW->icLoad = false;

  // Update for Delay: '<S290>/delayPhi'
  localDW->icLoad_o = false;

  // Update for Delay: '<S290>/delayBufferU'
  localDW->icLoad_e = false;

  // Update for Delay: '<S290>/delayBufferY'
  localDW->icLoad_d = false;

  // Update for Delay: '<S290>/delayL'
  localDW->icLoad_e2 = false;

  // Update for Delay: '<S291>/delayTheta'
  localDW->icLoad_f = false;

  // Update for Delay: '<S291>/delayPhi'
  localDW->icLoad_ej = false;

  // Update for Delay: '<S291>/delayBufferU'
  localDW->icLoad_j = false;

  // Update for Delay: '<S291>/delayBufferY'
  localDW->icLoad_i = false;

  // Update for Delay: '<S291>/delayL'
  localDW->icLoad_jb = false;

  // Update for Delay: '<S352>/MemoryX'
  localDW->icLoad_l = false;

  // Update for Delay: '<S352>/MemoryP'
  localDW->icLoad_h = false;
  (void)std::memcpy(&localDW->MemoryP_DSTATE[0], &rtb_Q[0], sizeof(real_T) <<
                    4UL);

  // Update for RandomNumber: '<S4>/Measurement Noise'
  localDW->NextOutput[0] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[0]) *
    localP->MeasurementNoise_StdDev[0] + localP->MeasurementNoise_Mean[0];
  localDW->NextOutput[1] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[1]) *
    localP->MeasurementNoise_StdDev[1] + localP->MeasurementNoise_Mean[1];
  localDW->NextOutput[2] = rt_nrand_Upu32_Yd_f_pw_snf(&localDW->RandSeed[2]) *
    localP->MeasurementNoise_StdDev[2] + localP->MeasurementNoise_Mean[2];
}

// Function for Chart: '<Root>/SupervisoryController'
boolean_T SupervisoryController::isequal(const event_bus varargin_1, const
  event_bus varargin_2)
{
  int32_T k;
  boolean_T exitg1;
  boolean_T out;
  boolean_T p_1;
  p_1 = false;
  out = true;
  k = 0;
  exitg1 = false;
  while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
          false) && (k < 3)) {
    if (varargin_1.nextChs[k] != varargin_2.nextChs[k]) {
      out = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  if (out) {
    k = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (k < 3)) {
      if (varargin_1.chs[k] != varargin_2.chs[k]) {
        out = false;
        exitg1 = true;
      } else {
        k++;
      }
    }

    if (out) {
      if (varargin_1.holdTime == varargin_2.holdTime) {
        if (varargin_1.moveTime == varargin_2.moveTime) {
          k = 0;
          exitg1 = false;
          while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U))
                  == false) && (k < 3)) {
            if (!(varargin_1.destPos[k] == varargin_2.destPos[k])) {
              out = false;
              exitg1 = true;
            } else {
              k++;
            }
          }

          if (out) {
            out = ((varargin_1.destState == varargin_2.destState) &&
                   (varargin_1.srcState == varargin_2.srcState));
          }
        } else {
          out = false;
        }
      } else {
        out = false;
      }
    }
  }

  if (out) {
    p_1 = true;
  }

  return p_1;
}

// Function for Chart: '<Root>/SupervisoryController'
real_T SupervisoryController::minimum_f(const real_T x_data[], const int32_T
  *x_size)
{
  real_T ex;
  int32_T last;
  last = *x_size;
  if (static_cast<int32_T>(static_cast<uint8_T>(static_cast<int32_T>(*x_size - 1)))
      + 1 <= 2) {
    if (static_cast<int32_T>(static_cast<uint8_T>(static_cast<int32_T>(*x_size -
           1))) + 1 == 1) {
      ex = x_data[0];
    } else {
      ex = x_data[*x_size - 1];
      if (x_data[0] > ex) {
      } else if (std::isnan(x_data[0])) {
        if (!std::isnan(ex)) {
        } else {
          ex = x_data[0];
        }
      } else {
        ex = x_data[0];
      }
    }
  } else {
    int32_T idx;
    int32_T k;
    if (!std::isnan(x_data[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= *x_size)) {
        if (!std::isnan(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      for (k = idx + 1; k <= last; k++) {
        real_T tmp;
        tmp = x_data[k - 1];
        if (ex > tmp) {
          ex = tmp;
        }
      }
    }
  }

  return ex;
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::computeProfileParams_n(real_T i, const real_T
  wayPoints_data[], const int32_T wayPoints_size[2], const real_T Vel_data[],
  const int32_T *Vel_size, const real_T TAc_data[], const int32_T *TAc_size,
  real_T *vParam, real_T *aParam, real_T *tAParam, real_T *tFParam)
{
  real_T TAcSwitch_data[3];
  real_T VelSwitch_data[3];
  real_T tACandidates[2];
  real_T b_sF;
  real_T s0;
  real_T sF;
  int32_T deltaSign;
  int32_T loop_ub;
  int8_T b_data[2];
  boolean_T inputCombo_idx_0;
  boolean_T inputCombo_idx_3;
  *aParam = 0.0;
  *vParam = 0.0;
  *tAParam = 0.0;
  *tFParam = 0.0;
  s0 = wayPoints_data[static_cast<int32_T>(i) - 1];
  sF = wayPoints_data[(static_cast<int32_T>(i) + wayPoints_size[0]) - 1];
  deltaSign = 1;
  if (sF < s0) {
    b_sF = s0;
    s0 = sF;
    sF = b_sF;
    deltaSign = -1;
  }

  inputCombo_idx_0 = (*Vel_size != 0);
  inputCombo_idx_3 = (*TAc_size != 0);
  if (inputCombo_idx_0) {
    if (*Vel_size - 1 >= 0) {
      (void)std::memcpy(&VelSwitch_data[0], &Vel_data[0], static_cast<uint32_T>(*
        Vel_size) * sizeof(real_T));
    }
  } else {
    loop_ub = static_cast<int32_T>(i);
    for (int32_T i_0{0}; i_0 < loop_ub; i_0++) {
      VelSwitch_data[i_0] = 1.0;
    }
  }

  if (inputCombo_idx_3) {
    if (*TAc_size - 1 >= 0) {
      (void)std::memcpy(&TAcSwitch_data[0], &TAc_data[0], static_cast<uint32_T>(*
        TAc_size) * sizeof(real_T));
    }
  } else {
    loop_ub = static_cast<int32_T>(i);
    for (int32_T i_0{0}; i_0 < loop_ub; i_0++) {
      TAcSwitch_data[i_0] = 1.0;
    }
  }

  switch ((inputCombo_idx_0 << 3UL) + (inputCombo_idx_3 ? static_cast<int32_T>(1)
           : static_cast<int32_T>(0))) {
   case 9:
    *vParam = VelSwitch_data[static_cast<int32_T>(i) - 1];
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *aParam = *vParam / *tAParam;
    *tFParam = ((*vParam * *tAParam + sF) - s0) / *vParam;
    break;

   case 8:
    *vParam = VelSwitch_data[static_cast<int32_T>(i) - 1];
    *tFParam = (sF - s0) * 1.5 / *vParam;
    *tAParam = (*vParam * *tFParam + (s0 - sF)) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 6:
    *aParam = 1.0;
    *tFParam = 1.0;
    b_sF = std::sqrt((4.0 * s0 + 1.0) - 4.0 * sF) / 2.0;
    tACandidates[0] = 0.5 - b_sF;
    tACandidates[1] = b_sF + 0.5;
    loop_ub = 0;
    if (0.5 - b_sF > 0.0) {
      b_data[0] = 1;
      loop_ub = 1;
    }

    if (b_sF + 0.5 > 0.0) {
      b_data[loop_ub] = 2;
    }

    *tAParam = tACandidates[b_data[0] - 1];
    *vParam = (sF - s0) / (1.0 - *tAParam);
    break;

   case 5:
    *aParam = 1.0;
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *vParam = *tAParam;
    *tFParam = ((*tAParam * *tAParam + sF) - s0) / *tAParam;
    break;

   case 4:
    *aParam = 1.0;
    b_sF = sF - s0;
    *tAParam = std::sqrt(b_sF / 2.0);
    *vParam = *tAParam;
    *tFParam = (*tAParam * *tAParam + b_sF) / *tAParam;
    break;

   case 3:
    *tFParam = 1.0;
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *vParam = (sF - s0) / (1.0 - *tAParam);
    *aParam = *vParam / *tAParam;
    break;

   case 2:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 1:
    {
      real_T aParam_tmp;
      *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
      b_sF = *tAParam * *tAParam;
      aParam_tmp = sF - s0;
      *aParam = aParam_tmp / (b_sF * 2.0);
      *vParam = *tAParam * *aParam;
      *tFParam = (b_sF * *aParam + aParam_tmp) / *vParam;
    }
    break;

   case 0:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   default:
    // no actions
    break;
  }

  if (s0 == sF) {
    *aParam = 0.0;
    *vParam = 0.0;
    if (std::isnan(*tFParam)) {
      *tFParam = 1.0;
    } else if (*tFParam == 0.0) {
      *tFParam = 1.0;
    } else {
      // no actions
    }

    *tAParam = *tFParam / 3.0;
  }

  *vParam *= static_cast<real_T>(deltaSign);
  *aParam *= static_cast<real_T>(deltaSign);
}

// Function for Chart: '<Root>/SupervisoryController'
boolean_T SupervisoryController::checkPolyForMultipleBreaks(const real_T
  breakMat_data[], const int32_T breakMat_size[2])
{
  real_T y[4];
  int32_T b;
  boolean_T hasMultipleBreaks;
  hasMultipleBreaks = false;
  b = breakMat_size[0];
  for (int32_T i{0}; i <= b - 2; i++) {
    int32_T k;
    boolean_T b_y;
    boolean_T exitg1;
    y[0] = std::abs(breakMat_data[i] - breakMat_data[i + 1]);
    y[1] = std::abs(breakMat_data[i + breakMat_size[0]] - breakMat_data[(i +
      breakMat_size[0]) + 1]);
    y[2] = std::abs(breakMat_data[(breakMat_size[0] << 1UL) + i] -
                    breakMat_data[((breakMat_size[0] << 1UL) + i) + 1]);
    y[3] = std::abs(breakMat_data[breakMat_size[0] * 3 + i] - breakMat_data
                    [(breakMat_size[0] * 3 + i) + 1]);
    b_y = false;
    k = 0;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (k < 4)) {
      if (y[k] > 2.2204460492503131E-16) {
        b_y = true;
        exitg1 = true;
      } else {
        k++;
      }
    }

    hasMultipleBreaks = (b_y || hasMultipleBreaks);
  }

  return hasMultipleBreaks;
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::processPolynomialResults(const real_T breakMat_data[],
  const int32_T breakMat_size[2], const real_T coeffMat_data[], const int32_T
  coeffMat_size[2], boolean_T hasMultipleBreaks, cell_wrap_5 breaksCell_data[],
  int32_T *breaksCell_size, cell_wrap_6 coeffCell_data[], int32_T
  *coeffCell_size)
{
  int32_T b;
  int8_T c_data[9];
  boolean_T coeffIndex_data[9];
  boolean_T nthCoeffIndex_data[3];
  b = breakMat_size[0];
  *coeffCell_size = breakMat_size[0];
  *breaksCell_size = breakMat_size[0];
  for (int32_T ii{0}; ii < b; ii++) {
    if (hasMultipleBreaks) {
      int32_T ibcol;
      int32_T itilerow;
      int32_T nrows;
      int32_T nthCoeffIndex_size_idx_0;
      nthCoeffIndex_size_idx_0 = breakMat_size[0];
      nrows = breakMat_size[0];
      if (nrows - 1 >= 0) {
        (void)std::memset(&nthCoeffIndex_data[0], 0, static_cast<uint32_T>(nrows)
                          * sizeof(boolean_T));
      }

      nthCoeffIndex_data[ii] = true;
      for (itilerow = 0; itilerow < 3; itilerow++) {
        ibcol = itilerow * nthCoeffIndex_size_idx_0;
        for (nrows = 0; nrows < nthCoeffIndex_size_idx_0; nrows++) {
          coeffIndex_data[ibcol + nrows] = nthCoeffIndex_data[nrows];
        }
      }

      nrows = static_cast<int32_T>(static_cast<int8_T>(breakMat_size[0] * 3)) -
        1;
      itilerow = 0;
      for (ibcol = 0; ibcol <= nrows; ibcol++) {
        if (coeffIndex_data[ibcol]) {
          itilerow++;
        }
      }

      nthCoeffIndex_size_idx_0 = itilerow;
      itilerow = 0;
      for (ibcol = 0; ibcol <= nrows; ibcol++) {
        if (coeffIndex_data[ibcol]) {
          c_data[itilerow] = static_cast<int8_T>(ibcol + 1);
          itilerow++;
        }
      }

      coeffCell_data[ii].f1.size[0] = nthCoeffIndex_size_idx_0;
      coeffCell_data[ii].f1.size[1] = 3;
      for (itilerow = 0; itilerow < 3; itilerow++) {
        for (ibcol = 0; ibcol < nthCoeffIndex_size_idx_0; ibcol++) {
          coeffCell_data[ii].f1.data[ibcol + coeffCell_data[ii].f1.size[0] *
            itilerow] = coeffMat_data[(coeffMat_size[0] * itilerow +
            static_cast<int32_T>(c_data[ibcol])) - 1];
        }
      }

      breaksCell_data[ii].f1[0] = breakMat_data[ii];
      breaksCell_data[ii].f1[1] = breakMat_data[ii + breakMat_size[0]];
      breaksCell_data[ii].f1[2] = breakMat_data[(breakMat_size[0] << 1UL) + ii];
      breaksCell_data[ii].f1[3] = breakMat_data[breakMat_size[0] * 3 + ii];
    } else {
      int32_T nrows;
      coeffCell_data[ii].f1.size[0] = coeffMat_size[0];
      coeffCell_data[ii].f1.size[1] = 3;
      nrows = coeffMat_size[0] * 3;
      for (int32_T itilerow{0}; itilerow < nrows; itilerow++) {
        coeffCell_data[ii].f1.data[itilerow] = coeffMat_data[itilerow];
      }

      breaksCell_data[ii].f1[0] = breakMat_data[0];
      breaksCell_data[ii].f1[1] = breakMat_data[breakMat_size[0]];
      breaksCell_data[ii].f1[2] = breakMat_data[breakMat_size[0] << 1UL];
      breaksCell_data[ii].f1[3] = breakMat_data[breakMat_size[0] * 3];
    }
  }
}

// Function for Chart: '<Root>/SupervisoryController'
real_T SupervisoryController::maximum_e(const real_T x_data[], const int32_T
  *x_size)
{
  real_T ex;
  int32_T last;
  last = *x_size;
  if (static_cast<int32_T>(static_cast<uint8_T>(static_cast<int32_T>(*x_size - 1)))
      + 1 <= 2) {
    if (static_cast<int32_T>(static_cast<uint8_T>(static_cast<int32_T>(*x_size -
           1))) + 1 == 1) {
      ex = x_data[0];
    } else {
      ex = x_data[*x_size - 1];
      if (x_data[0] < ex) {
      } else if (std::isnan(x_data[0])) {
        if (!std::isnan(ex)) {
        } else {
          ex = x_data[0];
        }
      } else {
        ex = x_data[0];
      }
    }
  } else {
    int32_T idx;
    int32_T k;
    if (!std::isnan(x_data[0])) {
      idx = 1;
    } else {
      boolean_T exitg1;
      idx = 0;
      k = 2;
      exitg1 = false;
      while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
              false) && (k <= *x_size)) {
        if (!std::isnan(x_data[k - 1])) {
          idx = k;
          exitg1 = true;
        } else {
          k++;
        }
      }
    }

    if (idx == 0) {
      ex = x_data[0];
    } else {
      ex = x_data[idx - 1];
      for (k = idx + 1; k <= last; k++) {
        real_T tmp;
        tmp = x_data[k - 1];
        if (ex < tmp) {
          ex = tmp;
        }
      }
    }
  }

  return ex;
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::linspace(real_T d2, uint16_T n, real_T y_data[],
  int32_T y_size[2])
{
  y_size[0] = 1;
  y_size[1] = static_cast<int32_T>(n);
  if (n >= 1UL) {
    y_data[static_cast<int32_T>(n) - 1] = d2;
    if (n >= 2UL) {
      y_data[0] = 0.0;
      if (n >= 3UL) {
        if (-d2 == 0.0) {
          real_T delta1;
          int32_T d;
          delta1 = d2 / (static_cast<real_T>(n) - 1.0);
          d = static_cast<int32_T>(n) - 1;
          for (int32_T c_k{2}; c_k <= d; c_k++) {
            y_data[c_k - 1] = (static_cast<real_T>(static_cast<int32_T>((c_k <<
              1UL) - static_cast<int32_T>(n))) - 1.0) * delta1;
          }

          if ((static_cast<int32_T>(n) & 1) == 1) {
            y_data[n >> 1UL] = 0.0;
          }
        } else {
          real_T delta1;
          int32_T d;
          boolean_T guard1{ false };

          guard1 = false;
          if (d2 < 0.0) {
            if (std::abs(d2) > 8.9884656743115785E+307) {
              delta1 = d2 / (static_cast<real_T>(n) - 1.0);
              d = static_cast<int32_T>(n);
              for (int32_T c_k{0}; c_k <= d - 3; c_k++) {
                y_data[c_k + 1] = (static_cast<real_T>(c_k) + 1.0) * delta1;
              }
            } else {
              guard1 = true;
            }
          } else {
            guard1 = true;
          }

          if (guard1) {
            delta1 = d2 / (static_cast<real_T>(n) - 1.0);
            d = static_cast<int32_T>(n);
            for (int32_T c_k{0}; c_k <= d - 3; c_k++) {
              y_data[c_k + 1] = (static_cast<real_T>(c_k) + 1.0) * delta1;
            }
          }
        }
      }
    }
  }
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::ppval(const s_vjEZ2dxatR8VOmLd9oOqoD *pp, const
  real_T x_data[], const int32_T x_size[2], real_T v_data[], int32_T v_size[2])
{
  real_T b_xloc;
  int32_T b_high_i;
  int32_T b_low_i;
  int32_T b_low_ip1;
  int32_T b_mid_i;
  int32_T coefStride;
  int32_T elementsPerPage;
  elementsPerPage = pp->coefs.size[0] - 1;
  coefStride = pp->coefs.size[0] * 5;
  v_size[0] = static_cast<int32_T>(static_cast<int16_T>(pp->coefs.size[0]));
  v_size[1] = static_cast<int32_T>(static_cast<int16_T>(x_size[1]));
  if (pp->coefs.size[0] == 1) {
    elementsPerPage = x_size[1];

#pragma omp parallel for num_threads(omp_get_max_threads()) private(b_xloc,b_low_i,b_low_ip1,b_high_i,b_mid_i)

    for (int32_T b_ix = 0; b_ix < elementsPerPage; b_ix++) {
      b_xloc = x_data[b_ix];
      if (std::isnan(b_xloc)) {
        v_data[b_ix] = b_xloc;
      } else {
        b_low_i = 0;
        b_low_ip1 = 2;
        b_high_i = 6;
        while (b_high_i > b_low_ip1) {
          b_mid_i = ((b_low_i + b_high_i) + 1) >> 1UL;
          if (b_xloc >= pp->breaks[b_mid_i - 1]) {
            b_low_i = b_mid_i - 1;
            b_low_ip1 = b_mid_i + 1;
          } else {
            b_high_i = b_mid_i;
          }
        }

        b_xloc -= pp->breaks[b_low_i];
        v_data[b_ix] = (pp->coefs.data[b_low_i] * b_xloc + pp->
                        coefs.data[b_low_i + coefStride]) * b_xloc +
          pp->coefs.data[(coefStride << 1UL) + b_low_i];
      }
    }
  } else {
    int32_T c;
    c = x_size[1];
    for (int32_T c_ix{0}; c_ix < c; c_ix++) {
      real_T xloc;
      int32_T iv0;
      iv0 = (elementsPerPage + 1) * c_ix;
      xloc = x_data[c_ix];
      if (std::isnan(xloc)) {
        for (int32_T high_i{0}; high_i <= elementsPerPage; high_i++) {
          v_data[iv0 + high_i] = xloc;
        }
      } else {
        int32_T high_i;
        int32_T ic0;
        int32_T low_i;
        int32_T low_ip1;
        low_i = 1;
        low_ip1 = 2;
        high_i = 6;
        while (high_i > low_ip1) {
          ic0 = (low_i + high_i) >> 1UL;
          if (xloc >= pp->breaks[ic0 - 1]) {
            low_i = ic0;
            low_ip1 = ic0 + 1;
          } else {
            high_i = ic0;
          }
        }

        low_ip1 = (low_i - 1) * (elementsPerPage + 1);
        xloc -= pp->breaks[low_i - 1];
        for (ic0 = 0; ic0 <= elementsPerPage; ic0++) {
          v_data[iv0 + ic0] = pp->coefs.data[low_ip1 + ic0];
        }

        for (high_i = 0; high_i < 2; high_i++) {
          int32_T tmp_0;
          int32_T vectorUB;
          ic0 = (high_i + 1) * coefStride + low_ip1;
          low_i = ((elementsPerPage + 1) / 2) << 1UL;
          vectorUB = low_i - 2;
          for (int32_T c_j{0}; c_j <= vectorUB; c_j += 2) {
            __m128d tmp;
            tmp_0 = iv0 + c_j;
            tmp = _mm_loadu_pd(&v_data[tmp_0]);
            (void)_mm_storeu_pd(&v_data[tmp_0], _mm_add_pd(_mm_mul_pd(tmp,
              _mm_set1_pd(xloc)), _mm_loadu_pd(&pp->coefs.data[ic0 + c_j])));
          }

          for (int32_T c_j{low_i}; c_j <= elementsPerPage; c_j++) {
            tmp_0 = iv0 + c_j;
            v_data[tmp_0] = v_data[tmp_0] * xloc + pp->coefs.data[ic0 + c_j];
          }
        }
      }
    }
  }
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::generateTrajectoriesFromCoefs(const real_T breaks[4],
  const real_T coeffs_data[], const int32_T coeffs_size[2], real_T dim, const
  real_T t_data[], const int32_T t_size[2], real_T q_data[], int32_T q_size[2],
  real_T qd_data[], int32_T qd_size[2], real_T qdd_data[], int32_T qdd_size[2],
  real_T pp_breaks[6], real_T pp_coefs_data[], int32_T pp_coefs_size[3])
{
  __m128d tmp;
  s_vjEZ2dxatR8VOmLd9oOqoD expl_temp;
  s_vjEZ2dxatR8VOmLd9oOqoD expl_temp_0;
  s_vjEZ2dxatR8VOmLd9oOqoD expl_temp_1;
  real_T dCoeffs_data[45];
  real_T ddCoeffs_data[45];
  real_T newCoefs_data[45];
  real_T coefsWithFlatStart_data[36];
  real_T valueAtEnd_data[12];
  real_T newSegmentCoeffs_data[9];
  real_T b_breaks[5];
  real_T valueAtStart_data[3];
  real_T coefsWithFlatStart_size_idx_0_0;
  real_T evalPointVector_idx_0;
  real_T evalPointVector_idx_1;
  real_T evalPointVector_idx_2;
  real_T s;
  int32_T b_k;
  int32_T c;
  int32_T coefsWithFlatStart_size_idx_0_t;
  int32_T e_i;
  int32_T m;
  int32_T newSegmentCoeffs_size_idx_0_tmp;
  int32_T o;
  int32_T valueAtEnd_size_idx_0;
  if (dim < 1.0) {
    c = 0;
  } else {
    c = static_cast<int32_T>(dim);
  }

  m = static_cast<int32_T>(static_cast<uint8_T>(c));
  for (e_i = 0; e_i < m; e_i++) {
    o = c + e_i;
    s = coeffs_data[e_i / c * coeffs_size[0] + e_i % c] * 0.0 + coeffs_data[o /
      c * coeffs_size[0] + o % c] * 0.0;
    o = (c << 1UL) + e_i;
    valueAtStart_data[e_i] = coeffs_data[o / c * coeffs_size[0] + o % c] + s;
  }

  newSegmentCoeffs_size_idx_0_tmp = static_cast<int32_T>(dim);
  o = static_cast<int32_T>(dim) * 3;
  if (o - 1 >= 0) {
    (void)std::memset(&newSegmentCoeffs_data[0], 0, static_cast<uint32_T>(o) *
                      sizeof(real_T));
  }

  for (m = 0; m < c; m++) {
    newSegmentCoeffs_data[m + (static_cast<int32_T>(dim) << 1UL)] =
      valueAtStart_data[m];
  }

  coefsWithFlatStart_size_idx_0_0 = static_cast<real_T>(coeffs_size[0]) + dim;
  coefsWithFlatStart_size_idx_0_t = static_cast<int32_T>
    (coefsWithFlatStart_size_idx_0_0);
  e_i = static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0) * 3;
  if (e_i - 1 >= 0) {
    (void)std::memset(&coefsWithFlatStart_data[0], 0, static_cast<uint32_T>(e_i)
                      * sizeof(real_T));
  }

  for (m = 0; m < 3; m++) {
    for (b_k = 0; b_k < newSegmentCoeffs_size_idx_0_tmp; b_k++) {
      coefsWithFlatStart_data[b_k + static_cast<int32_T>
        (coefsWithFlatStart_size_idx_0_0) * m] = newSegmentCoeffs_data[
        static_cast<int32_T>(dim) * m + b_k];
    }
  }

  if (dim + 1.0 > coefsWithFlatStart_size_idx_0_0) {
    c = 1;
  } else {
    c = static_cast<int32_T>(static_cast<real_T>(dim + 1.0));
  }

  e_i = coeffs_size[0];
  for (m = 0; m < 3; m++) {
    for (b_k = 0; b_k < e_i; b_k++) {
      coefsWithFlatStart_data[((c + b_k) + static_cast<int32_T>
        (coefsWithFlatStart_size_idx_0_0) * m) - 1] = coeffs_data[coeffs_size[0]
        * m + b_k];
    }
  }

  b_breaks[0] = breaks[0] - 1.0;
  b_breaks[1] = breaks[0];
  b_breaks[2] = breaks[1];
  b_breaks[3] = breaks[2];
  b_breaks[4] = breaks[3];
  s = breaks[3] - breaks[2];
  evalPointVector_idx_0 = rt_powd_snf(s, 2.0);
  evalPointVector_idx_1 = rt_powd_snf(s, 1.0);
  evalPointVector_idx_2 = rt_powd_snf(s, 0.0);
  s = (static_cast<real_T>(static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0))
       - dim) + 1.0;
  if (s > static_cast<real_T>(static_cast<int32_T>
       (coefsWithFlatStart_size_idx_0_0))) {
    c = 0;
    m = 0;
  } else {
    c = static_cast<int32_T>(s) - 1;
    m = static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0);
  }

  m -= c;
  valueAtEnd_size_idx_0 = m;
  e_i = static_cast<int32_T>(static_cast<uint8_T>(m));
  for (b_k = 0; b_k < e_i; b_k++) {
    int32_T p_0;
    p_0 = m + b_k;
    s = coefsWithFlatStart_data[(b_k % m + c) + b_k / m * static_cast<int32_T>
      (coefsWithFlatStart_size_idx_0_0)] * evalPointVector_idx_0 +
      coefsWithFlatStart_data[(p_0 % m + c) + p_0 / m * static_cast<int32_T>
      (coefsWithFlatStart_size_idx_0_0)] * evalPointVector_idx_1;
    p_0 = (m << 1UL) + b_k;
    valueAtEnd_data[b_k] = coefsWithFlatStart_data[(p_0 % m + c) + p_0 / m *
      static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0)] *
      evalPointVector_idx_2 + s;
  }

  if (o - 1 >= 0) {
    (void)std::memset(&newSegmentCoeffs_data[0], 0, static_cast<uint32_T>(o) *
                      sizeof(real_T));
  }

  for (m = 0; m < valueAtEnd_size_idx_0; m++) {
    newSegmentCoeffs_data[m + (static_cast<int32_T>(dim) << 1UL)] =
      valueAtEnd_data[m];
  }

  valueAtEnd_size_idx_0 = static_cast<int32_T>(static_cast<real_T>(static_cast<
    real_T>(static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0)) + dim));
  e_i = valueAtEnd_size_idx_0 * 3;
  if (e_i - 1 >= 0) {
    (void)std::memset(&newCoefs_data[0], 0, static_cast<uint32_T>(e_i) * sizeof
                      (real_T));
  }

  for (m = 0; m < 3; m++) {
    for (b_k = 0; b_k < coefsWithFlatStart_size_idx_0_t; b_k++) {
      newCoefs_data[b_k + valueAtEnd_size_idx_0 * m] = coefsWithFlatStart_data[
        static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0) * m + b_k];
    }
  }

  for (m = 0; m < 3; m++) {
    for (b_k = 0; b_k < newSegmentCoeffs_size_idx_0_tmp; b_k++) {
      newCoefs_data[((static_cast<int32_T>(coefsWithFlatStart_size_idx_0_0) +
                      (b_k + 1)) + valueAtEnd_size_idx_0 * m) - 1] =
        newSegmentCoeffs_data[static_cast<int32_T>(dim) * m + b_k];
    }
  }

  for (m = 0; m < 5; m++) {
    pp_breaks[m] = b_breaks[m];
  }

  pp_breaks[5] = breaks[3] + 1.0;
  newSegmentCoeffs_size_idx_0_tmp = static_cast<int32_T>(static_cast<int8_T>
    (valueAtEnd_size_idx_0));
  e_i = static_cast<int32_T>(static_cast<int8_T>(valueAtEnd_size_idx_0)) * 3;
  if (e_i - 1 >= 0) {
    (void)std::memset(&dCoeffs_data[0], 0, static_cast<uint32_T>(e_i) * sizeof
                      (real_T));
  }

  for (c = 0; c < 2; c++) {
    o = (valueAtEnd_size_idx_0 / 2) << 1UL;
    b_k = o - 2;
    for (m = 0; m <= b_k; m += 2) {
      tmp = _mm_loadu_pd(&newCoefs_data[valueAtEnd_size_idx_0 * c + m]);
      (void)_mm_storeu_pd(&dCoeffs_data[m + static_cast<int32_T>
                          (static_cast<int8_T>(valueAtEnd_size_idx_0)) * (c + 1)],
                          _mm_mul_pd(_mm_set1_pd(static_cast<real_T>(
        static_cast<int32_T>(2 - c))), tmp));
    }

    for (m = o; m < valueAtEnd_size_idx_0; m++) {
      dCoeffs_data[m + static_cast<int32_T>(static_cast<int8_T>
        (valueAtEnd_size_idx_0)) * (c + 1)] =
        newCoefs_data[valueAtEnd_size_idx_0 * c + m] * static_cast<real_T>(
        static_cast<int32_T>(2 - c));
    }
  }

  e_i = static_cast<int32_T>(static_cast<int8_T>(valueAtEnd_size_idx_0)) * 3;
  if (e_i - 1 >= 0) {
    (void)std::memset(&ddCoeffs_data[0], 0, static_cast<uint32_T>(e_i) * sizeof
                      (real_T));
  }

  for (c = 0; c < 2; c++) {
    o = (static_cast<int32_T>(static_cast<int8_T>(valueAtEnd_size_idx_0)) / 2) <<
      1UL;
    b_k = o - 2;
    for (m = 0; m <= b_k; m += 2) {
      tmp = _mm_loadu_pd(&dCoeffs_data[static_cast<int32_T>(static_cast<int8_T>
        (valueAtEnd_size_idx_0)) * c + m]);
      (void)_mm_storeu_pd(&ddCoeffs_data[m + static_cast<int32_T>
                          (static_cast<int8_T>(valueAtEnd_size_idx_0)) * (c + 1)],
                          _mm_mul_pd(_mm_set1_pd(static_cast<real_T>(
        static_cast<int32_T>(2 - c))), tmp));
    }

    for (m = o; m < newSegmentCoeffs_size_idx_0_tmp; m++) {
      ddCoeffs_data[m + static_cast<int32_T>(static_cast<int8_T>
        (valueAtEnd_size_idx_0)) * (c + 1)] = dCoeffs_data[static_cast<int32_T>(
        static_cast<int8_T>(valueAtEnd_size_idx_0)) * c + m] *
        static_cast<real_T>(static_cast<int32_T>(2 - c));
    }
  }

  pp_coefs_size[0] = static_cast<int32_T>(dim);
  pp_coefs_size[1] = 5;
  pp_coefs_size[2] = 3;
  o = static_cast<int32_T>(dim) * 5 * 3;
  expl_temp.coefs.size[0] = static_cast<int32_T>(dim);
  expl_temp.coefs.size[1] = 5;
  expl_temp.coefs.size[2] = 3;
  if (o - 1 >= 0) {
    (void)std::memcpy(&pp_coefs_data[0], &newCoefs_data[0], static_cast<uint32_T>
                      (o) * sizeof(real_T));
  }

  if (o - 1 >= 0) {
    (void)std::memcpy(&expl_temp.coefs.data[0], &newCoefs_data[0],
                      static_cast<uint32_T>(o) * sizeof(real_T));
  }

  for (m = 0; m < 6; m++) {
    expl_temp.breaks[m] = pp_breaks[m];
  }

  ppval(&expl_temp, t_data, t_size, q_data, q_size);
  expl_temp_0.coefs.size[0] = static_cast<int32_T>(dim);
  expl_temp_0.coefs.size[1] = 5;
  expl_temp_0.coefs.size[2] = 3;
  if (o - 1 >= 0) {
    (void)std::memcpy(&expl_temp_0.coefs.data[0], &dCoeffs_data[0], static_cast<
                      uint32_T>(o) * sizeof(real_T));
  }

  for (m = 0; m < 6; m++) {
    expl_temp_0.breaks[m] = pp_breaks[m];
  }

  ppval(&expl_temp_0, t_data, t_size, qd_data, qd_size);
  expl_temp_1.coefs.size[0] = static_cast<int32_T>(dim);
  expl_temp_1.coefs.size[1] = 5;
  expl_temp_1.coefs.size[2] = 3;
  if (o - 1 >= 0) {
    (void)std::memcpy(&expl_temp_1.coefs.data[0], &ddCoeffs_data[0],
                      static_cast<uint32_T>(o) * sizeof(real_T));
  }

  for (m = 0; m < 6; m++) {
    expl_temp_1.breaks[m] = pp_breaks[m];
  }

  ppval(&expl_temp_1, t_data, t_size, qdd_data, qdd_size);
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::trapveltraj_p(const real_T wayPoints_data[], const
  int32_T wayPoints_size[2], uint16_T numSamples, real_T varargin_2, real_T
  q_data[], int32_T q_size[2])
{
  cell_wrap_5 breaksCell_data[3];
  emxArray_cell_wrap_6_3 coeffsCell;
  emxArray_s_vjEZ2dxatR8VOmLd9oOq ppCell;
  real_T coeffMat_data[27];
  real_T parameterMat_data[18];
  real_T breakMat_data[12];
  real_T coefs[9];
  real_T tAc_data[3];
  real_T vel_data[3];
  real_T segATime;
  real_T segAcc;
  real_T segFTime;
  real_T segVel;
  int32_T breakMat_size[2];
  int32_T c_size[2];
  int32_T coeffMat_size[2];
  int32_T d_size[2];
  int32_T cellSelection;
  int32_T i;
  int32_T indivPolyDim;
  int32_T j_size_idx_0;
  int32_T loop_ub;
  int32_T n;
  int32_T numComputedPolynomials;
  int32_T parameterMat_size_idx_0;
  int32_T tAc_size;
  int32_T vel_size;
  int8_T f_data[9];
  int8_T j_data[9];
  int8_T lspbSegIndices_data[9];
  int8_T rowSelection_data[3];
  int8_T j[2];
  boolean_T coefIndex_data[9];
  boolean_T hasMultipleBreaks;
  n = wayPoints_size[0];
  vel_size = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  cellSelection = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  for (int32_T c_i{0}; c_i < cellSelection; c_i++) {
    vel_data[c_i] = varargin_2;
  }

  tAc_size = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  cellSelection = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  for (int32_T c_i{0}; c_i < cellSelection; c_i++) {
    tAc_data[c_i] = 0.01;
  }

  q_size[0] = wayPoints_size[0];
  q_size[1] = static_cast<int32_T>(numSamples);
  cellSelection = wayPoints_size[0] * static_cast<int32_T>(numSamples);
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&q_data[0], 0, static_cast<uint32_T>(cellSelection) *
                      sizeof(real_T));
  }

  parameterMat_size_idx_0 = wayPoints_size[0];
  cellSelection = wayPoints_size[0] * 6;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&parameterMat_data[0], 0, static_cast<uint32_T>
                      (cellSelection) * sizeof(real_T));
  }

  coeffMat_size[0] = 3 * wayPoints_size[0];
  coeffMat_size[1] = 3;
  cellSelection = 3 * wayPoints_size[0] * 3;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&coeffMat_data[0], 0, static_cast<uint32_T>(cellSelection)
                      * sizeof(real_T));
  }

  breakMat_size[0] = wayPoints_size[0];
  breakMat_size[1] = 4;
  cellSelection = wayPoints_size[0] << 2UL;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&breakMat_data[0], 0, static_cast<uint32_T>(cellSelection)
                      * sizeof(real_T));
  }

  indivPolyDim = wayPoints_size[0];
  if (wayPoints_size[0] - 1 >= 0) {
    loop_ub = 3 * wayPoints_size[0];
  }

  for (i = 0; i < indivPolyDim; i++) {
    real_T parameterMat_data_tmp;
    computeProfileParams_n(static_cast<real_T>(i) + 1.0, wayPoints_data,
      wayPoints_size, vel_data, &vel_size, tAc_data, &tAc_size, &segVel, &segAcc,
      &segATime, &segFTime);
    parameterMat_data_tmp = wayPoints_data[i];
    parameterMat_data[i] = parameterMat_data_tmp;
    parameterMat_data[i + parameterMat_size_idx_0] = wayPoints_data[i +
      wayPoints_size[0]];
    parameterMat_data[i + (parameterMat_size_idx_0 << 1UL)] = segVel;
    parameterMat_data[i + parameterMat_size_idx_0 * 3] = segAcc;
    parameterMat_data[i + (parameterMat_size_idx_0 << 2UL)] = segATime;
    parameterMat_data[i + parameterMat_size_idx_0 * 5] = segFTime;
    (void)std::memset(&coefs[0], 0, 9U * sizeof(real_T));
    if (segVel == 0.0) {
      coefs[6] = parameterMat_data_tmp;
      coefs[7] = parameterMat_data_tmp;
      coefs[8] = parameterMat_data_tmp;
    } else {
      real_T coefs_tmp;
      coefs[0] = segAcc / 2.0;
      coefs[3] = 0.0;
      coefs[6] = parameterMat_data_tmp;
      coefs[1] = 0.0;
      coefs[4] = segVel;
      coefs_tmp = segAcc / 2.0 * (segATime * segATime);
      coefs[7] = coefs_tmp + parameterMat_data_tmp;
      coefs[2] = -segAcc / 2.0;
      coefs[5] = segVel;
      coefs[8] = (wayPoints_data[i + wayPoints_size[0]] + coefs_tmp) - segVel *
        segATime;
    }

    if (loop_ub - 1 >= 0) {
      (void)std::memset(&coefIndex_data[0], 0, static_cast<uint32_T>(loop_ub) *
                        sizeof(boolean_T));
    }

    segVel = static_cast<real_T>(static_cast<int32_T>((n << 1UL) + i)) + 1.0;
    if ((n == 0) || (static_cast<int32_T>(segVel) < i + 1)) {
      cellSelection = 0;
    } else {
      numComputedPolynomials = static_cast<int32_T>(static_cast<real_T>((segVel
        - static_cast<real_T>(static_cast<int32_T>(i + 1))) / static_cast<real_T>
        (n)));
      cellSelection = numComputedPolynomials + 1;
      for (int32_T c_i{0}; c_i <= numComputedPolynomials; c_i++) {
        lspbSegIndices_data[c_i] = static_cast<int8_T>((n * c_i + i) + 1);
      }
    }

    if (cellSelection - 1 >= 0) {
      (void)std::memcpy(&f_data[0], &lspbSegIndices_data[0],
                        static_cast<uint32_T>(cellSelection) * sizeof(int8_T));
    }

    for (int32_T c_i{0}; c_i < cellSelection; c_i++) {
      coefIndex_data[f_data[c_i] - 1] = true;
    }

    numComputedPolynomials = 3 * n - 1;
    cellSelection = 0;
    for (int32_T c_i{0}; c_i <= numComputedPolynomials; c_i++) {
      if (coefIndex_data[c_i]) {
        cellSelection++;
      }
    }

    j_size_idx_0 = cellSelection;
    cellSelection = 0;
    for (int32_T c_i{0}; c_i <= numComputedPolynomials; c_i++) {
      if (coefIndex_data[c_i]) {
        j_data[cellSelection] = static_cast<int8_T>(c_i + 1);
        cellSelection++;
      }
    }

    j[0] = static_cast<int8_T>(j_size_idx_0);
    cellSelection = static_cast<int32_T>(static_cast<int8_T>(j_size_idx_0));
    for (int32_T c_i{0}; c_i < 3; c_i++) {
      for (j_size_idx_0 = 0; j_size_idx_0 < cellSelection; j_size_idx_0++) {
        coeffMat_data[(static_cast<int32_T>(j_data[j_size_idx_0]) +
                       coeffMat_size[0] * c_i) - 1] = coefs[static_cast<int32_T>
          (j[0]) * c_i + j_size_idx_0];
      }
    }

    segVel = breakMat_data[i];
    breakMat_data[i + breakMat_size[0]] = segATime + segVel;
    breakMat_data[i + (breakMat_size[0] << 1UL)] = (segFTime - segATime) +
      segVel;
    breakMat_data[i + breakMat_size[0] * 3] = segFTime + segVel;
  }

  hasMultipleBreaks = checkPolyForMultipleBreaks(breakMat_data, breakMat_size);
  processPolynomialResults(breakMat_data, breakMat_size, coeffMat_data,
    coeffMat_size, hasMultipleBreaks, breaksCell_data, &tAc_size,
    coeffsCell.data, &coeffsCell.size);
  if (wayPoints_size[0] == 0) {
    vel_size = 0;
  } else {
    vel_size = wayPoints_size[0];
    for (i = 0; i < parameterMat_size_idx_0; i++) {
      vel_data[i] = parameterMat_data[parameterMat_size_idx_0 * 5 + i];
    }
  }

  linspace(maximum_e(vel_data, &vel_size), numSamples, rtDW.t_data,
           coeffMat_size);
  if (hasMultipleBreaks) {
    numComputedPolynomials = wayPoints_size[0];
    indivPolyDim = 1;
  } else {
    numComputedPolynomials = 1;
    indivPolyDim = wayPoints_size[0];
  }

  i = numComputedPolynomials;
  for (numComputedPolynomials = 0; numComputedPolynomials < i;
       numComputedPolynomials++) {
    if (hasMultipleBreaks) {
      rowSelection_data[0] = static_cast<int8_T>(numComputedPolynomials + 1);
      cellSelection = numComputedPolynomials;
    } else {
      if (n < 1) {
      } else {
        cellSelection = n - 1;
        for (int32_T c_i{0}; c_i <= cellSelection; c_i++) {
          rowSelection_data[c_i] = static_cast<int8_T>(c_i + 1);
        }
      }

      cellSelection = 0;
    }

    generateTrajectoriesFromCoefs(breaksCell_data[cellSelection].f1,
      coeffsCell.data[cellSelection].f1.data, coeffsCell.data[cellSelection].
      f1.size, static_cast<real_T>(indivPolyDim), rtDW.t_data, coeffMat_size,
      rtDW.b_data, breakMat_size, rtDW.c_data, c_size, rtDW.d_data, d_size,
      ppCell.data[numComputedPolynomials].breaks,
      ppCell.data[numComputedPolynomials].coefs.data,
      ppCell.data[numComputedPolynomials].coefs.size);
    cellSelection = breakMat_size[1];
    for (int32_T c_i{0}; c_i < cellSelection; c_i++) {
      loop_ub = breakMat_size[0];
      for (j_size_idx_0 = 0; j_size_idx_0 < loop_ub; j_size_idx_0++) {
        q_data[(static_cast<int32_T>(rowSelection_data[j_size_idx_0]) + q_size[0]
                * c_i) - 1] = rtDW.b_data[breakMat_size[0] * c_i + j_size_idx_0];
      }
    }
  }
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::computeProfileParams(real_T i, const real_T
  wayPoints_data[], const int32_T wayPoints_size[2], const real_T Vel_data[],
  const int32_T *Vel_size, real_T *vParam, real_T *aParam, real_T *tAParam,
  real_T *tFParam)
{
  real_T VelSwitch_data[3];
  real_T tACandidates[2];
  real_T b_sF;
  real_T s0;
  real_T sF;
  int32_T deltaSign;
  int32_T loop_ub;
  int8_T b_data[2];
  boolean_T inputCombo_idx_0;
  *aParam = 0.0;
  *vParam = 0.0;
  *tAParam = 0.0;
  *tFParam = 0.0;
  s0 = wayPoints_data[static_cast<int32_T>(i) - 1];
  sF = wayPoints_data[(static_cast<int32_T>(i) + wayPoints_size[0]) - 1];
  deltaSign = 1;
  if (sF < s0) {
    b_sF = s0;
    s0 = sF;
    sF = b_sF;
    deltaSign = -1;
  }

  inputCombo_idx_0 = (*Vel_size != 0);
  if (inputCombo_idx_0) {
    if (*Vel_size - 1 >= 0) {
      (void)std::memcpy(&VelSwitch_data[0], &Vel_data[0], static_cast<uint32_T>(*
        Vel_size) * sizeof(real_T));
    }
  } else {
    loop_ub = static_cast<int32_T>(i);
    for (int32_T i_0{0}; i_0 < loop_ub; i_0++) {
      VelSwitch_data[i_0] = 1.0;
    }
  }

  switch (inputCombo_idx_0 << 3UL) {
   case 8:
    *vParam = VelSwitch_data[static_cast<int32_T>(i) - 1];
    *tFParam = (sF - s0) * 1.5 / *vParam;
    *tAParam = (*vParam * *tFParam + (s0 - sF)) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 6:
    *aParam = 1.0;
    *tFParam = 1.0;
    b_sF = std::sqrt((4.0 * s0 + 1.0) - 4.0 * sF) / 2.0;
    tACandidates[0] = 0.5 - b_sF;
    tACandidates[1] = b_sF + 0.5;
    loop_ub = 0;
    if (0.5 - b_sF > 0.0) {
      b_data[0] = 1;
      loop_ub = 1;
    }

    if (b_sF + 0.5 > 0.0) {
      b_data[loop_ub] = 2;
    }

    *tAParam = tACandidates[b_data[0] - 1];
    *vParam = (sF - s0) / (1.0 - *tAParam);
    break;

   case 5:
    *aParam = 1.0;
    *tAParam = 1.0;
    *vParam = 1.0;
    *tFParam = (sF + 1.0) - s0;
    break;

   case 4:
    *aParam = 1.0;
    b_sF = sF - s0;
    *tAParam = std::sqrt(b_sF / 2.0);
    *vParam = *tAParam;
    *tFParam = (*tAParam * *tAParam + b_sF) / *tAParam;
    break;

   case 3:
    *tFParam = 1.0;
    *tAParam = 1.0;
    *vParam = sF - s0 > 0.0 ? (rtInf) : sF - s0 < 0.0 ? (rtMinusInf) : (rtNaN);
    *aParam = *vParam;
    break;

   case 2:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 1:
    *tAParam = 1.0;
    b_sF = sF - s0;
    *aParam = b_sF / 2.0;
    *vParam = *aParam;
    *tFParam = (b_sF + *aParam) / *aParam;
    break;

   case 0:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   default:
    // no actions
    break;
  }

  if (s0 == sF) {
    *aParam = 0.0;
    *vParam = 0.0;
    if (std::isnan(*tFParam)) {
      *tFParam = 1.0;
    } else if (*tFParam == 0.0) {
      *tFParam = 1.0;
    } else {
      // no actions
    }

    *tAParam = *tFParam / 3.0;
  }

  *vParam *= static_cast<real_T>(deltaSign);
  *aParam *= static_cast<real_T>(deltaSign);
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::trapveltraj(const real_T wayPoints_data[], const
  int32_T wayPoints_size[2], uint16_T numSamples, real_T varargin_2, real_T
  q_data[], int32_T q_size[2])
{
  cell_wrap_5 breaksCell_data[3];
  emxArray_cell_wrap_6_3 coeffsCell;
  emxArray_s_vjEZ2dxatR8VOmLd9oOq ppCell;
  real_T coeffMat_data[27];
  real_T parameterMat_data[18];
  real_T breakMat_data[12];
  real_T coefs[9];
  real_T vel_data[3];
  real_T segATime;
  real_T segAcc;
  real_T segFTime;
  real_T segVel;
  int32_T breakMat_size[2];
  int32_T coeffMat_size[2];
  int32_T g_size[2];
  int32_T h_size[2];
  int32_T c_i;
  int32_T cellSelection;
  int32_T i;
  int32_T indivPolyDim;
  int32_T j_size_idx_0;
  int32_T loop_ub;
  int32_T n;
  int32_T numComputedPolynomials;
  int32_T parameterMat_size_idx_0;
  int32_T vel_size;
  int8_T c_data[9];
  int8_T j_data[9];
  int8_T lspbSegIndices_data[9];
  int8_T rowSelection_data[3];
  int8_T j[2];
  boolean_T coefIndex_data[9];
  boolean_T hasMultipleBreaks;
  n = wayPoints_size[0];
  vel_size = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  cellSelection = static_cast<int32_T>(static_cast<int8_T>(wayPoints_size[0]));
  for (c_i = 0; c_i < cellSelection; c_i++) {
    vel_data[c_i] = varargin_2;
  }

  q_size[0] = wayPoints_size[0];
  q_size[1] = static_cast<int32_T>(numSamples);
  cellSelection = wayPoints_size[0] * static_cast<int32_T>(numSamples);
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&q_data[0], 0, static_cast<uint32_T>(cellSelection) *
                      sizeof(real_T));
  }

  parameterMat_size_idx_0 = wayPoints_size[0];
  cellSelection = wayPoints_size[0] * 6;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&parameterMat_data[0], 0, static_cast<uint32_T>
                      (cellSelection) * sizeof(real_T));
  }

  coeffMat_size[0] = 3 * wayPoints_size[0];
  coeffMat_size[1] = 3;
  cellSelection = 3 * wayPoints_size[0] * 3;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&coeffMat_data[0], 0, static_cast<uint32_T>(cellSelection)
                      * sizeof(real_T));
  }

  breakMat_size[0] = wayPoints_size[0];
  breakMat_size[1] = 4;
  cellSelection = wayPoints_size[0] << 2UL;
  if (cellSelection - 1 >= 0) {
    (void)std::memset(&breakMat_data[0], 0, static_cast<uint32_T>(cellSelection)
                      * sizeof(real_T));
  }

  indivPolyDim = wayPoints_size[0];
  if (wayPoints_size[0] - 1 >= 0) {
    loop_ub = 3 * wayPoints_size[0];
  }

  for (i = 0; i < indivPolyDim; i++) {
    real_T parameterMat_data_tmp;
    computeProfileParams(static_cast<real_T>(i) + 1.0, wayPoints_data,
                         wayPoints_size, vel_data, &vel_size, &segVel, &segAcc,
                         &segATime, &segFTime);
    parameterMat_data_tmp = wayPoints_data[i];
    parameterMat_data[i] = parameterMat_data_tmp;
    parameterMat_data[i + parameterMat_size_idx_0] = wayPoints_data[i +
      wayPoints_size[0]];
    parameterMat_data[i + (parameterMat_size_idx_0 << 1UL)] = segVel;
    parameterMat_data[i + parameterMat_size_idx_0 * 3] = segAcc;
    parameterMat_data[i + (parameterMat_size_idx_0 << 2UL)] = segATime;
    parameterMat_data[i + parameterMat_size_idx_0 * 5] = segFTime;
    (void)std::memset(&coefs[0], 0, 9U * sizeof(real_T));
    if (segVel == 0.0) {
      coefs[6] = parameterMat_data_tmp;
      coefs[7] = parameterMat_data_tmp;
      coefs[8] = parameterMat_data_tmp;
    } else {
      real_T coefs_tmp;
      coefs[0] = segAcc / 2.0;
      coefs[3] = 0.0;
      coefs[6] = parameterMat_data_tmp;
      coefs[1] = 0.0;
      coefs[4] = segVel;
      coefs_tmp = segAcc / 2.0 * (segATime * segATime);
      coefs[7] = coefs_tmp + parameterMat_data_tmp;
      coefs[2] = -segAcc / 2.0;
      coefs[5] = segVel;
      coefs[8] = (wayPoints_data[i + wayPoints_size[0]] + coefs_tmp) - segVel *
        segATime;
    }

    if (loop_ub - 1 >= 0) {
      (void)std::memset(&coefIndex_data[0], 0, static_cast<uint32_T>(loop_ub) *
                        sizeof(boolean_T));
    }

    segVel = static_cast<real_T>(static_cast<int32_T>((n << 1UL) + i)) + 1.0;
    if ((n == 0) || (static_cast<int32_T>(segVel) < i + 1)) {
      cellSelection = 0;
    } else {
      numComputedPolynomials = static_cast<int32_T>(static_cast<real_T>((segVel
        - static_cast<real_T>(static_cast<int32_T>(i + 1))) / static_cast<real_T>
        (n)));
      cellSelection = numComputedPolynomials + 1;
      for (c_i = 0; c_i <= numComputedPolynomials; c_i++) {
        lspbSegIndices_data[c_i] = static_cast<int8_T>((n * c_i + i) + 1);
      }
    }

    if (cellSelection - 1 >= 0) {
      (void)std::memcpy(&c_data[0], &lspbSegIndices_data[0],
                        static_cast<uint32_T>(cellSelection) * sizeof(int8_T));
    }

    for (c_i = 0; c_i < cellSelection; c_i++) {
      coefIndex_data[c_data[c_i] - 1] = true;
    }

    numComputedPolynomials = 3 * n - 1;
    cellSelection = 0;
    for (c_i = 0; c_i <= numComputedPolynomials; c_i++) {
      if (coefIndex_data[c_i]) {
        cellSelection++;
      }
    }

    j_size_idx_0 = cellSelection;
    cellSelection = 0;
    for (c_i = 0; c_i <= numComputedPolynomials; c_i++) {
      if (coefIndex_data[c_i]) {
        j_data[cellSelection] = static_cast<int8_T>(c_i + 1);
        cellSelection++;
      }
    }

    j[0] = static_cast<int8_T>(j_size_idx_0);
    cellSelection = static_cast<int32_T>(static_cast<int8_T>(j_size_idx_0));
    for (c_i = 0; c_i < 3; c_i++) {
      for (j_size_idx_0 = 0; j_size_idx_0 < cellSelection; j_size_idx_0++) {
        coeffMat_data[(static_cast<int32_T>(j_data[j_size_idx_0]) +
                       coeffMat_size[0] * c_i) - 1] = coefs[static_cast<int32_T>
          (j[0]) * c_i + j_size_idx_0];
      }
    }

    segVel = breakMat_data[i];
    breakMat_data[i + breakMat_size[0]] = segATime + segVel;
    breakMat_data[i + (breakMat_size[0] << 1UL)] = (segFTime - segATime) +
      segVel;
    breakMat_data[i + breakMat_size[0] * 3] = segFTime + segVel;
  }

  hasMultipleBreaks = checkPolyForMultipleBreaks(breakMat_data, breakMat_size);
  processPolynomialResults(breakMat_data, breakMat_size, coeffMat_data,
    coeffMat_size, hasMultipleBreaks, breaksCell_data, &c_i, coeffsCell.data,
    &coeffsCell.size);
  if (wayPoints_size[0] == 0) {
    vel_size = 0;
  } else {
    vel_size = wayPoints_size[0];
    for (i = 0; i < parameterMat_size_idx_0; i++) {
      vel_data[i] = parameterMat_data[parameterMat_size_idx_0 * 5 + i];
    }
  }

  linspace(maximum_e(vel_data, &vel_size), numSamples, rtDW.l_data,
           coeffMat_size);
  if (hasMultipleBreaks) {
    numComputedPolynomials = wayPoints_size[0];
    indivPolyDim = 1;
  } else {
    numComputedPolynomials = 1;
    indivPolyDim = wayPoints_size[0];
  }

  i = numComputedPolynomials;
  for (numComputedPolynomials = 0; numComputedPolynomials < i;
       numComputedPolynomials++) {
    if (hasMultipleBreaks) {
      rowSelection_data[0] = static_cast<int8_T>(numComputedPolynomials + 1);
      cellSelection = numComputedPolynomials;
    } else {
      if (n < 1) {
      } else {
        cellSelection = n - 1;
        for (c_i = 0; c_i <= cellSelection; c_i++) {
          rowSelection_data[c_i] = static_cast<int8_T>(c_i + 1);
        }
      }

      cellSelection = 0;
    }

    generateTrajectoriesFromCoefs(breaksCell_data[cellSelection].f1,
      coeffsCell.data[cellSelection].f1.data, coeffsCell.data[cellSelection].
      f1.size, static_cast<real_T>(indivPolyDim), rtDW.l_data, coeffMat_size,
      rtDW.f_data, breakMat_size, rtDW.g_data, g_size, rtDW.h_data, h_size,
      ppCell.data[numComputedPolynomials].breaks,
      ppCell.data[numComputedPolynomials].coefs.data,
      ppCell.data[numComputedPolynomials].coefs.size);
    cellSelection = breakMat_size[1];
    for (c_i = 0; c_i < cellSelection; c_i++) {
      loop_ub = breakMat_size[0];
      for (j_size_idx_0 = 0; j_size_idx_0 < loop_ub; j_size_idx_0++) {
        q_data[(static_cast<int32_T>(rowSelection_data[j_size_idx_0]) + q_size[0]
                * c_i) - 1] = rtDW.f_data[breakMat_size[0] * c_i + j_size_idx_0];
      }
    }
  }
}

//
// Function for Chart: '<Root>/SupervisoryController'
// function [trajectory, trajectorySize] = trajGen(event, y_i)
//
void SupervisoryController::trajGen(const event_bus event, const real_T y_i[3],
  real_T trajectory[7200], uint16_T *trajectorySize)
{
  real_T y_i_data[6];
  real_T tmp_data[3];
  real_T tmp;
  int32_T tmp_size[2];
  int32_T y_i_size[2];
  int32_T g_size_idx_0;
  int32_T h_size_idx_0;
  int8_T g_data[3];
  int8_T h_data[3];
  int8_T i_data[3];

  // MATLAB Function 'trajGen': '<S1>:50'
  // '<S1>:50:2' [trajectory, trajectorySize] = trajectoryGen(event, y_i, y_max, dt); 
  // 'trajectoryGen:3' traj = zeros(3, 2400);
  (void)std::memset(&trajectory[0], 0, 7200U * sizeof(real_T));

  // 'trajectoryGen:4' chs = event.chs;
  // 'trajectoryGen:5' numWaypts = cast(event.moveTime/dt, "uint16");
  tmp = std::round(event.moveTime / rtP.dt);
  if (tmp < 65536.0) {
    if (tmp >= 0.0) {
      *trajectorySize = static_cast<uint16_T>(tmp);
    } else {
      *trajectorySize = 0U;
    }
  } else {
    *trajectorySize = MAX_uint16_T;
  }

  // 'trajectoryGen:6' assert(numWaypts < 1200);
  // 'trajectoryGen:8' if event.srcState == event.destState
  if (event.srcState == event.destState) {
    int32_T trueCount;

    // 'trajectoryGen:9' traj(chs, 1:numWaypts) = trapveltraj(...
    // 'trajectoryGen:10'         [y(chs), y_max(chs).*event.destPos(chs)],...
    // 'trajectoryGen:11'         numWaypts,...
    // 'trajectoryGen:12'         PeakVelocity=min(y_max(chs))/event.moveTime);
    trueCount = 0;
    if (event.chs[0]) {
      trueCount = 1;
    }

    if (event.chs[1]) {
      trueCount++;
    }

    if (event.chs[2]) {
      trueCount++;
    }

    g_size_idx_0 = trueCount;
    trueCount = 0;
    if (event.chs[0]) {
      g_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      g_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      g_data[trueCount] = 3;
    }

    trueCount = 0;
    if (event.chs[0]) {
      h_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      h_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      h_data[trueCount] = 3;
    }

    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      // Inport: '<Root>/y_max'
      tmp_data[trueCount] = rtU.y_max[g_data[trueCount] - 1];
    }

    y_i_size[0] = g_size_idx_0;
    y_i_size[1] = 2;
    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      y_i_data[trueCount] = y_i[g_data[trueCount] - 1];
    }

    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      h_size_idx_0 = static_cast<int32_T>(g_data[trueCount]) - 1;

      // Inport: '<Root>/y_max'
      y_i_data[trueCount + g_size_idx_0] = rtU.y_max[h_size_idx_0] *
        event.destPos[h_size_idx_0];
    }

    trapveltraj(y_i_data, y_i_size, *trajectorySize, minimum_f(tmp_data,
      &g_size_idx_0) / event.moveTime, rtDW.tmp_data, tmp_size);
    g_size_idx_0 = tmp_size[1];
    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      int32_T y_i_0;
      y_i_0 = tmp_size[0];
      for (h_size_idx_0 = 0; h_size_idx_0 < y_i_0; h_size_idx_0++) {
        trajectory[(static_cast<int32_T>(h_data[h_size_idx_0]) + 3 * trueCount)
          - 1] = rtDW.tmp_data[tmp_size[0] * trueCount + h_size_idx_0];
      }
    }
  } else {
    int32_T trueCount;
    uint16_T tmp_0;

    // 'trajectoryGen:13' else
    //  for state transition, use rectangular velocity profile
    // 'trajectoryGen:14' traj(chs, 1:2*numWaypts) = trapveltraj(...
    // 'trajectoryGen:15'         [y(chs), 2*y_max(chs).*event.destPos(chs)],... 
    // 'trajectoryGen:16'         2*numWaypts,...
    // 'trajectoryGen:17'         PeakVelocity=min(y_max(chs))/event.moveTime,... 
    // 'trajectoryGen:18'         AccelTime=0.01);
    trueCount = 0;
    if (event.chs[0]) {
      trueCount = 1;
    }

    if (event.chs[1]) {
      trueCount++;
    }

    if (event.chs[2]) {
      trueCount++;
    }

    g_size_idx_0 = trueCount;
    trueCount = 0;
    if (event.chs[0]) {
      g_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      g_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      g_data[trueCount] = 3;
    }

    trueCount = 0;
    if (event.chs[0]) {
      trueCount = 1;
    }

    if (event.chs[1]) {
      trueCount++;
    }

    if (event.chs[2]) {
      trueCount++;
    }

    h_size_idx_0 = trueCount;
    trueCount = 0;
    if (event.chs[0]) {
      h_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      h_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      h_data[trueCount] = 3;
    }

    trueCount = 0;
    if (event.chs[0]) {
      i_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      i_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      i_data[trueCount] = 3;
    }

    for (trueCount = 0; trueCount < h_size_idx_0; trueCount++) {
      // Inport: '<Root>/y_max'
      tmp_data[trueCount] = rtU.y_max[h_data[trueCount] - 1];
    }

    y_i_size[0] = h_size_idx_0;
    y_i_size[1] = 2;
    for (trueCount = 0; trueCount < h_size_idx_0; trueCount++) {
      y_i_data[trueCount] = y_i[h_data[trueCount] - 1];
    }

    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      // Inport: '<Root>/y_max'
      y_i_data[trueCount + h_size_idx_0] = rtU.y_max[g_data[trueCount] - 1] *
        2.0 * event.destPos[h_data[trueCount] - 1];
    }

    if (*trajectorySize > 32767UL) {
      tmp_0 = MAX_uint16_T;
    } else {
      tmp_0 = static_cast<uint16_T>(*trajectorySize << 1UL);
    }

    trapveltraj_p(y_i_data, y_i_size, tmp_0, minimum_f(tmp_data, &h_size_idx_0) /
                  event.moveTime, rtDW.tmp_data_m, tmp_size);
    g_size_idx_0 = tmp_size[1];
    for (trueCount = 0; trueCount < g_size_idx_0; trueCount++) {
      int32_T y_i_0;
      y_i_0 = tmp_size[0];
      for (h_size_idx_0 = 0; h_size_idx_0 < y_i_0; h_size_idx_0++) {
        trajectory[(static_cast<int32_T>(i_data[h_size_idx_0]) + 3 * trueCount)
          - 1] = rtDW.tmp_data_m[tmp_size[0] * trueCount + h_size_idx_0];
      }
    }
  }
}

//
// Function for Chart: '<Root>/SupervisoryController'
// function [inTransitionRegion, eventDone, waypoint, holdTime] = handleEvent(event)
//
void SupervisoryController::handleEvent(const event_bus event, boolean_T
  *inTransitionRegion, boolean_T *eventDone, uint16_T *waypoint, real_T
  *holdTime) const
{
  int32_T c_size_idx_0;
  int32_T trueCount;
  int8_T b_data[3];
  int8_T c_data[3];
  int8_T d_data[3];
  int8_T e_data[3];
  int8_T f_data[3];
  boolean_T x_data[3];
  boolean_T exitg1;
  boolean_T y;

  // MATLAB Function 'handleEvent': '<S1>:43'
  // '<S1>:43:2' [inTransitionRegion, eventDone, waypoint, holdTime] = processEvent(event, waypt, holdT, y, trajSize, y_max, dt); 
  //  initialize flags
  // 'processEvent:4' evDone = false;
  *eventDone = false;

  // 'processEvent:5' inTransRegion = false;
  // 'processEvent:6' chs = event.chs;
  // 'processEvent:7' nextChs = event.nextChs;
  //  increment waypoint
  // 'processEvent:10' if currWaypt < numWaypts
  if (rtDW.waypt < rtDW.trajSize) {
    // 'processEvent:11' nextWaypt = currWaypt + 1;
    *waypoint = static_cast<uint16_T>(static_cast<uint32_T>(rtDW.waypt) + 1U);

    // 'processEvent:12' nextHoldTime = currHoldTime;
    *holdTime = rtDW.holdT;
  } else {
    // 'processEvent:13' else
    // 'processEvent:14' nextWaypt = currWaypt;
    *waypoint = rtDW.waypt;

    // 'processEvent:15' nextHoldTime = currHoldTime + dt;
    *holdTime = rtDW.holdT + rtP.dt;
  }

  //  check if we're in transition region
  // 'processEvent:19' if all(y(chs) < 0.9*y_max(chs)) && all(y(chs) > 0.1*y_max(chs)) 
  trueCount = 0;
  if (event.chs[0]) {
    b_data[0] = 1;
    trueCount = 1;
  }

  if (event.chs[1]) {
    b_data[trueCount] = 2;
    trueCount++;
  }

  if (event.chs[2]) {
    b_data[trueCount] = 3;
  }

  trueCount = 0;
  if (event.chs[0]) {
    trueCount = 1;
  }

  if (event.chs[1]) {
    trueCount++;
  }

  if (event.chs[2]) {
    trueCount++;
  }

  c_size_idx_0 = trueCount;
  trueCount = 0;
  if (event.chs[0]) {
    c_data[0] = 1;
    trueCount = 1;
  }

  if (event.chs[1]) {
    c_data[trueCount] = 2;
    trueCount++;
  }

  if (event.chs[2]) {
    c_data[trueCount] = 3;
  }

  for (trueCount = 0; trueCount < c_size_idx_0; trueCount++) {
    // Inport: '<Root>/y' incorporates:
    //   Inport: '<Root>/y_max'

    x_data[trueCount] = (rtU.y[c_data[trueCount] - 1] <
                         rtU.y_max[b_data[trueCount] - 1] * 0.9);
  }

  y = true;
  trueCount = 1;
  exitg1 = false;
  while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
          false) && (trueCount <= c_size_idx_0)) {
    if (!x_data[trueCount - 1]) {
      y = false;
      exitg1 = true;
    } else {
      trueCount++;
    }
  }

  if (y) {
    trueCount = 0;
    if (event.chs[0]) {
      d_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      d_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      d_data[trueCount] = 3;
    }

    trueCount = 0;
    if (event.chs[0]) {
      trueCount = 1;
    }

    if (event.chs[1]) {
      trueCount++;
    }

    if (event.chs[2]) {
      trueCount++;
    }

    c_size_idx_0 = trueCount;
    trueCount = 0;
    if (event.chs[0]) {
      e_data[0] = 1;
      trueCount = 1;
    }

    if (event.chs[1]) {
      e_data[trueCount] = 2;
      trueCount++;
    }

    if (event.chs[2]) {
      e_data[trueCount] = 3;
    }

    for (trueCount = 0; trueCount < c_size_idx_0; trueCount++) {
      // Inport: '<Root>/y' incorporates:
      //   Inport: '<Root>/y_max'

      x_data[trueCount] = (rtU.y[e_data[trueCount] - 1] >
                           rtU.y_max[d_data[trueCount] - 1] * 0.1);
    }

    trueCount = 1;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (trueCount <= c_size_idx_0)) {
      if (!x_data[trueCount - 1]) {
        y = false;
        exitg1 = true;
      } else {
        trueCount++;
      }
    }

    if (y) {
      // 'processEvent:20' inTransRegion = false;
      *inTransitionRegion = false;
    } else {
      // 'processEvent:21' else
      // 'processEvent:22' inTransRegion = true;
      *inTransitionRegion = true;
    }
  } else {
    // 'processEvent:21' else
    // 'processEvent:22' inTransRegion = true;
    *inTransitionRegion = true;
  }

  //  remain in current state
  // 'processEvent:26' if event.srcState == event.destState
  if (event.srcState == event.destState) {
    // 'processEvent:27' if nextWaypt == currWaypt & nextHoldTime > event.holdTime 
    if ((*waypoint == rtDW.waypt) && (*holdTime > event.holdTime)) {
      // 'processEvent:28' evDone = true;
      *eventDone = true;
    }

    //  transition to new state
  } else {
    // 'processEvent:32' else
    //  event.srcState ~= event.destState
    // 'processEvent:33' if all(y(nextChs))
    trueCount = 0;
    if (event.nextChs[0]) {
      trueCount = 1;
    }

    if (event.nextChs[1]) {
      trueCount++;
    }

    if (event.nextChs[2]) {
      trueCount++;
    }

    c_size_idx_0 = trueCount;
    trueCount = 0;
    if (event.nextChs[0]) {
      f_data[0] = 1;
      trueCount = 1;
    }

    if (event.nextChs[1]) {
      f_data[trueCount] = 2;
      trueCount++;
    }

    if (event.nextChs[2]) {
      f_data[trueCount] = 3;
    }

    y = true;
    trueCount = 1;
    exitg1 = false;
    while (((exitg1 ? static_cast<uint32_T>(1U) : static_cast<uint32_T>(0U)) ==
            false) && (trueCount <= c_size_idx_0)) {
      if (rtU.y[f_data[trueCount - 1] - 1] == 0.0) {
        y = false;
        exitg1 = true;
      } else {
        trueCount++;
      }
    }

    if (y) {
      //  measurements are available in all state1 channels
      // 'processEvent:34' evDone = true;
      *eventDone = true;
    }
  }
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::enter_internal_State1(void)
{
  // Entry Internal 'State1': '<S1>:176'
  // Entry Internal 'EventHandler': '<S1>:177'
  // Transition: '<S1>:164'
  rtDW.is_EventHandler_n = IN_RequestEvent;

  // Entry 'RequestEvent': '<S1>:163'
  // '<S1>:163:3' evDone = false;
  rtDW.evDone = false;

  // '<S1>:163:4' if waypt == 1
  if (rtDW.waypt == 1UL) {
    // Inport: '<Root>/y'
    //  hold curr pos
    // '<S1>:163:6' traj(:, waypt) = y;
    rtDW.traj[0] = rtU.y[0];
    rtDW.traj[1] = rtU.y[1];
    rtDW.traj[2] = rtU.y[2];
  } else {
    int32_T traj_tmp;

    // '<S1>:163:7' else
    //  hold last waypoint pos
    // '<S1>:163:9' traj(:,1) = traj(:, waypt);
    traj_tmp = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
    rtDW.traj[0] = rtDW.traj[traj_tmp];
    rtDW.traj[1] = rtDW.traj[traj_tmp + 1];
    rtDW.traj[2] = rtDW.traj[traj_tmp + 2];

    // '<S1>:163:10' waypt = 1;
    rtDW.waypt = 1U;
  }

  // Outport: '<Root>/requestEvent'
  // '<S1>:163:12' requestEvent = true;
  rtY.requestEvent = true;

  // Inport: '<Root>/y_o'
  //  request new event
  // Entry 'controlLaw': '<S1>:162'
  // '<S1>:162:3' yhat1 = y_o(chs1);
  rtDW.yhat1[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[0]) - 1];
  rtDW.yhat1[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[1]) - 1];
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::enter_internal_State0(void)
{
  // Entry Internal 'State0': '<S1>:1'
  // Entry Internal 'EventHandler': '<S1>:61'
  // Transition: '<S1>:64'
  rtDW.is_EventHandler = IN_RequestEvent;

  // Entry 'RequestEvent': '<S1>:65'
  // '<S1>:65:3' evDone = false;
  rtDW.evDone = false;

  // '<S1>:65:4' if waypt == 1
  if (rtDW.waypt == 1UL) {
    // Inport: '<Root>/y'
    //  hold curr pos
    // '<S1>:65:6' traj(:, waypt) = y;
    rtDW.traj[0] = rtU.y[0];
    rtDW.traj[1] = rtU.y[1];
    rtDW.traj[2] = rtU.y[2];
  } else {
    int32_T traj_tmp;

    // '<S1>:65:7' else
    //  hold last waypoint pos
    // '<S1>:65:9' traj(:,1) = traj(:, waypt);
    traj_tmp = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
    rtDW.traj[0] = rtDW.traj[traj_tmp];
    rtDW.traj[1] = rtDW.traj[traj_tmp + 1];
    rtDW.traj[2] = rtDW.traj[traj_tmp + 2];

    // '<S1>:65:10' waypt = 1;
    rtDW.waypt = 1U;
  }

  // Outport: '<Root>/requestEvent'
  // '<S1>:65:12' requestEvent = true;
  rtY.requestEvent = true;

  // Inport: '<Root>/y_o'
  //  request new event
  // Entry 'controlLaw': '<S1>:59'
  // '<S1>:59:3' yhat0 = y_o(chs0);
  rtDW.yhat0 = rtU.y_o[static_cast<int32_T>(rtDW.chs0) - 1];
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::State1(const int32_T *sfEvent)
{
  real_T r_k[2];
  real_T y0_g[2];
  real_T y_i[2];
  real_T holdT;
  uint16_T waypt;
  int8_T b_data[2];

  // Outport: '<Root>/currEv' incorporates:
  //   Inport: '<Root>/nextEv'
  //   Inport: '<Root>/y'
  //   Inport: '<Root>/y_o'
  //   Outport: '<Root>/currTraj'
  //   Outport: '<Root>/inTransRegion'
  //   Outport: '<Root>/requestEvent'
  //   Outport: '<Root>/yhat'

  // During 'State1': '<S1>:176'
  // '<S1>:47:1' sf_internal_predicateOutput = currEv.destState == 0 & evDone;
  if ((rtY.currEv.destState == 0.0) && rtDW.evDone) {
    // Disable for Function Call SubSystem: '<S1>/State1.controlLaw.AMPC1'
    // Transition: '<S1>:47'
    // Exit Internal 'State1': '<S1>:176'
    // Exit 'controlLaw': '<S1>:162'
    State1controlLawAMPC1_Disable(&rtDW.State1controlLawAMPC1_g,
      &rtP.State1controlLawAMPC1_g);

    // End of Disable for SubSystem: '<S1>/State1.controlLaw.AMPC1'
    // Exit Internal 'EventHandler': '<S1>:177'
    if (static_cast<uint32_T>(rtDW.is_EventHandler_n) == IN_RequestEvent) {
      // Exit 'RequestEvent': '<S1>:163'
      // '<S1>:163:14' requestEvent = false;
      rtDW.is_EventHandler_n = IN_NO_ACTIVE_CHILD;
    } else {
      rtDW.is_EventHandler_n = IN_NO_ACTIVE_CHILD;
    }

    rtDW.is_c6_SupervisoryController = IN_State0;

    // Entry 'State0': '<S1>:1'
    // '<S1>:1:3' waypt = 1;
    rtDW.waypt = 1U;

    // '<S1>:1:4' traj = zeros(3, 2400);
    (void)std::memset(&rtDW.traj[0], 0, 7200U * sizeof(real_T));

    // '<S1>:1:5' inTransRegion = true;
    rtY.inTransRegion = true;
    enter_internal_State0();

    // '<S1>:204:1' sf_internal_predicateOutput = currEv.destState == 2 & evDone; 
  } else if ((rtY.currEv.destState == 2.0) && rtDW.evDone) {
    // Disable for Function Call SubSystem: '<S1>/State1.controlLaw.AMPC1'
    // Transition: '<S1>:204'
    // Exit Internal 'State1': '<S1>:176'
    // Exit 'controlLaw': '<S1>:162'
    State1controlLawAMPC1_Disable(&rtDW.State1controlLawAMPC1_g,
      &rtP.State1controlLawAMPC1_g);

    // End of Disable for SubSystem: '<S1>/State1.controlLaw.AMPC1'
    // Exit Internal 'EventHandler': '<S1>:177'
    if (static_cast<uint32_T>(rtDW.is_EventHandler_n) == IN_RequestEvent) {
      // Exit 'RequestEvent': '<S1>:163'
      // '<S1>:163:14' requestEvent = false;
      rtDW.is_EventHandler_n = IN_NO_ACTIVE_CHILD;
    } else {
      rtDW.is_EventHandler_n = IN_NO_ACTIVE_CHILD;
    }

    rtDW.is_c6_SupervisoryController = IN_State2;

    // Entry 'State2': '<S1>:193'
    // '<S1>:193:3' waypt = 1;
    rtDW.waypt = 1U;

    // '<S1>:193:4' traj = zeros(3, 2400);
    (void)std::memset(&rtDW.traj[0], 0, 7200U * sizeof(real_T));

    // '<S1>:193:5' inTransRegion = true;
    rtY.inTransRegion = true;

    // Entry Internal 'State2': '<S1>:193'
    // Entry Internal 'EventHandler': '<S1>:194'
    // Transition: '<S1>:189'
    rtDW.is_EventHandler_k = IN_RequestEvent;

    // Entry 'RequestEvent': '<S1>:191'
    // '<S1>:191:3' evDone = false;
    rtDW.evDone = false;

    // '<S1>:191:4' if waypt == 1
    //  hold curr pos
    // '<S1>:191:6' traj(:, waypt) = y;
    rtDW.traj[0] = rtU.y[0];
    rtDW.traj[1] = rtU.y[1];
    rtDW.traj[2] = rtU.y[2];

    // '<S1>:191:12' requestEvent = true;
    rtY.requestEvent = true;

    //  request new event
    // Entry 'controlLaw': '<S1>:203'
    // '<S1>:203:3' yhat2 = y_o(chs2);
    rtDW.yhat2[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs2[0]) - 1];
    rtDW.yhat2[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs2[1]) - 1];
  } else {
    int32_T i;

    // During 'EventHandler': '<S1>:177'
    if (static_cast<uint32_T>(rtDW.is_EventHandler_n) == IN_HandleEvent) {
      // During 'HandleEvent': '<S1>:166'
      // '<S1>:165:1' sf_internal_predicateOutput = evDone;
      if (rtDW.evDone) {
        // Transition: '<S1>:165'
        rtDW.is_EventHandler_n = IN_RequestEvent;

        // Entry 'RequestEvent': '<S1>:163'
        // '<S1>:163:3' evDone = false;
        rtDW.evDone = false;

        // '<S1>:163:4' if waypt == 1
        if (rtDW.waypt == 1UL) {
          //  hold curr pos
          // '<S1>:163:6' traj(:, waypt) = y;
          rtDW.traj[0] = rtU.y[0];
          rtDW.traj[1] = rtU.y[1];
          rtDW.traj[2] = rtU.y[2];
        } else {
          // '<S1>:163:7' else
          //  hold last waypoint pos
          // '<S1>:163:9' traj(:,1) = traj(:, waypt);
          i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
          rtDW.traj[0] = rtDW.traj[i];
          rtDW.traj[1] = rtDW.traj[i + 1];
          rtDW.traj[2] = rtDW.traj[i + 2];

          // '<S1>:163:10' waypt = 1;
          rtDW.waypt = 1U;
        }

        // '<S1>:163:12' requestEvent = true;
        rtY.requestEvent = true;

        //  request new event
      } else {
        // Outport: '<Root>/inTransRegion'
        // '<S1>:166:7' [inTransRegion, evDone, waypt, holdT]...
        // '<S1>:166:8'     = handleEvent(currEv);
        handleEvent(rtY.currEv, &rtY.inTransRegion, &rtDW.evDone, &waypt, &holdT);
        rtDW.waypt = waypt;
        rtDW.holdT = holdT;
      }

      // During 'RequestEvent': '<S1>:163'
      // '<S1>:167:1' sf_internal_predicateOutput = ~isequal(nextEv, nullEv);
    } else if (!isequal(rtU.nextEv, rtP.nullEv)) {
      // Transition: '<S1>:167'
      // '<S1>:167:1' evDone = false;
      rtDW.evDone = false;

      // Exit 'RequestEvent': '<S1>:163'
      // '<S1>:163:14' requestEvent = false;
      rtY.requestEvent = false;
      rtDW.is_EventHandler_n = IN_HandleEvent;

      // Entry 'HandleEvent': '<S1>:166'
      // '<S1>:166:3' currEv = nextEv;
      rtY.currEv = rtU.nextEv;

      // '<S1>:166:4' [traj, trajSize] = trajGen(currEv, y);
      trajGen(rtY.currEv, rtU.y, rtDW.traj, &rtDW.trajSize);

      // '<S1>:166:5' holdT = 0;
      rtDW.holdT = 0.0;
    } else {
      // no actions
    }

    // During 'controlLaw': '<S1>:162'
    // '<S1>:162:5' if trueMeasAvail
    if (*sfEvent == static_cast<int32_T>(event_trueMeasAvail)) {
      // '<S1>:162:6' [u, yhat1] = AMPC1(traj(chs1, waypt), y(chs1), y_o(chs1), u_o, false); 
      i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
      r_k[0] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs1[0])) - 1];
      y_i[0] = rtU.y[static_cast<int32_T>(rtDW.chs1[0]) - 1];
      y0_g[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[0]) - 1];
      r_k[1] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs1[1])) - 1];
      y_i[1] = rtU.y[static_cast<int32_T>(rtDW.chs1[1]) - 1];
      y0_g[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[1]) - 1];

      // Outputs for Function Call SubSystem: '<S1>/State1.controlLaw.AMPC1'
      // Outport: '<Root>/u' incorporates:
      //   Inport: '<Root>/u_o'
      //   Inport: '<Root>/y'
      //   Inport: '<Root>/y_o'

      // Simulink Function 'AMPC1': '<S1>:175'
      State1controlLawAMPC1(r_k, y_i, y0_g, rtU.u_o, 0.0, rtY.u, rtDW.yhat1,
                            &rtDW.State1controlLawAMPC1_g,
                            &rtP.State1controlLawAMPC1_g, &rtP,
                            &rtPrevZCX.State1controlLawAMPC1_g);

      // End of Outputs for SubSystem: '<S1>/State1.controlLaw.AMPC1'
    } else if ((*sfEvent == static_cast<int32_T>(event_simMeasAvail)) &&
               rtY.inTransRegion) {
      real_T y_i_tmp;
      int32_T c_size_idx_0;

      // '<S1>:162:7' elseif simMeasAvail & inTransRegion
      //  use yhat instead of y for ch's w/o measurements
      // '<S1>:162:9' ysim1 = y(chs1);
      // '<S1>:162:10' ysim1(find(ysim1 == 0)) = yhat1(find(ysim1 == 0));
      holdT = rtU.y[static_cast<int32_T>(rtDW.chs1[0]) - 1];
      y_i[0] = holdT;
      y_i_tmp = rtU.y[static_cast<int32_T>(rtDW.chs1[1]) - 1];
      y_i[1] = y_i_tmp;
      i = 0;
      if (holdT == 0.0) {
        b_data[0] = 1;
        i = 1;
      }

      if (y_i_tmp == 0.0) {
        b_data[i] = 2;
      }

      i = 0;
      if (holdT == 0.0) {
        i = 1;
      }

      if (y_i_tmp == 0.0) {
        i++;
      }

      c_size_idx_0 = i;
      i = 0;
      if (holdT == 0.0) {
        r_k[0] = rtDW.yhat1[0];
        i = 1;
      }

      if (y_i_tmp == 0.0) {
        r_k[i] = rtDW.yhat1[1];
      }

      for (i = 0; i < c_size_idx_0; i++) {
        y_i[b_data[i] - 1] = r_k[i];
      }

      // '<S1>:162:11' [u, yhat1] = AMPC1(traj(chs1, waypt), ysim1, y_o(chs1), u_o, false); 
      i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
      r_k[0] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs1[0])) - 1];
      y0_g[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[0]) - 1];
      r_k[1] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs1[1])) - 1];
      y0_g[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[1]) - 1];

      // Outputs for Function Call SubSystem: '<S1>/State1.controlLaw.AMPC1'
      // Outport: '<Root>/u' incorporates:
      //   Inport: '<Root>/u_o'
      //   Inport: '<Root>/y'
      //   Inport: '<Root>/y_o'

      // Simulink Function 'AMPC1': '<S1>:175'
      State1controlLawAMPC1(r_k, y_i, y0_g, rtU.u_o, 0.0, rtY.u, rtDW.yhat1,
                            &rtDW.State1controlLawAMPC1_g,
                            &rtP.State1controlLawAMPC1_g, &rtP,
                            &rtPrevZCX.State1controlLawAMPC1_g);

      // End of Outputs for SubSystem: '<S1>/State1.controlLaw.AMPC1'
    } else {
      // no actions
    }

    // '<S1>:162:13' yhat = [0; yhat1];
    rtY.yhat[0] = 0.0;
    rtY.yhat[1] = rtDW.yhat1[0];
    rtY.yhat[2] = rtDW.yhat1[1];

    // '<S1>:162:14' currTraj = traj(:, waypt);
    rtY.currTraj[0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 3];
    i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
    rtY.currTraj[1] = rtDW.traj[i + 1];
    rtY.currTraj[2] = rtDW.traj[i + 2];
  }

  // End of Outport: '<Root>/currEv'
}

// Function for Chart: '<Root>/SupervisoryController'
void SupervisoryController::chartstep_c6_SupervisoryControl(const int32_T
  *sfEvent)
{
  real_T r_ki[2];
  real_T y0_i[2];
  real_T y_j[2];
  real_T ysim0;
  uint16_T waypt;
  int8_T b_data[2];

  // Chart: '<Root>/SupervisoryController' incorporates:
  //   Inport: '<Root>/nextEv'
  //   Inport: '<Root>/y'
  //   Inport: '<Root>/y_o'
  //   Outport: '<Root>/currEv'
  //   Outport: '<Root>/currTraj'
  //   Outport: '<Root>/inTransRegion'
  //   Outport: '<Root>/requestEvent'
  //   Outport: '<Root>/yhat'

  // During: SupervisoryController
  if (static_cast<uint32_T>(rtDW.is_active_c6_SupervisoryControl) == 0U) {
    // Entry: SupervisoryController
    rtDW.is_active_c6_SupervisoryControl = 1U;

    // Entry Internal: SupervisoryController
    // Transition: '<S1>:2'
    // '<S1>:2:1' currEv = struct('srcState', 0, 'destState', 0, 'destPos', [0; 0; 0], 'moveTime', 0, 'holdTime', 0, 'chs', logical([1; 0; 0]), 'nextChs', logical([1; 0; 0])); 
    rtY.currEv.srcState = 0.0;
    rtY.currEv.destState = 0.0;
    rtY.currEv.moveTime = 0.0;
    rtY.currEv.holdTime = 0.0;
    rtY.currEv.destPos[0] = 0.0;
    rtY.currEv.chs[0] = true;
    rtY.currEv.nextChs[0] = true;
    rtY.currEv.destPos[1] = 0.0;
    rtY.currEv.chs[1] = false;
    rtY.currEv.nextChs[1] = false;
    rtY.currEv.destPos[2] = 0.0;
    rtY.currEv.chs[2] = false;
    rtY.currEv.nextChs[2] = false;
    rtDW.is_c6_SupervisoryController = IN_State0;

    // Entry 'State0': '<S1>:1'
    // '<S1>:1:3' waypt = 1;
    rtDW.waypt = 1U;

    // '<S1>:1:4' traj = zeros(3, 2400);
    (void)std::memset(&rtDW.traj[0], 0, 7200U * sizeof(real_T));

    // Outport: '<Root>/inTransRegion' incorporates:
    //   Outport: '<Root>/currEv'

    // '<S1>:1:5' inTransRegion = true;
    rtY.inTransRegion = true;
    enter_internal_State0();
  } else {
    switch (rtDW.is_c6_SupervisoryController) {
     case IN_State0:
      {
        // During 'State0': '<S1>:1'
        // '<S1>:46:1' sf_internal_predicateOutput = currEv.destState == 1 & evDone; 
        if ((rtY.currEv.destState == 1.0) && rtDW.evDone) {
          // Disable for Function Call SubSystem: '<S1>/State0.controlLaw.AMPC0' 
          // Transition: '<S1>:46'
          // Exit Internal 'State0': '<S1>:1'
          // Exit 'controlLaw': '<S1>:59'
          State0controlLawAMPC0_Disable(&rtDW.State0controlLawAMPC0_n,
            &rtP.State0controlLawAMPC0_n);

          // End of Disable for SubSystem: '<S1>/State0.controlLaw.AMPC0'
          // Exit Internal 'EventHandler': '<S1>:61'
          if (static_cast<uint32_T>(rtDW.is_EventHandler) == IN_RequestEvent) {
            // Exit 'RequestEvent': '<S1>:65'
            // '<S1>:65:14' requestEvent = false;
            rtDW.is_EventHandler = IN_NO_ACTIVE_CHILD;
          } else {
            rtDW.is_EventHandler = IN_NO_ACTIVE_CHILD;
          }

          rtDW.is_c6_SupervisoryController = IN_State1;

          // Entry 'State1': '<S1>:176'
          // '<S1>:176:3' waypt = 1;
          rtDW.waypt = 1U;

          // '<S1>:176:4' traj = zeros(3, 2400);
          (void)std::memset(&rtDW.traj[0], 0, 7200U * sizeof(real_T));

          // '<S1>:176:5' inTransRegion = true;
          rtY.inTransRegion = true;
          enter_internal_State1();
        } else {
          int32_T c_size_idx_0;

          // During 'EventHandler': '<S1>:61'
          if (static_cast<uint32_T>(rtDW.is_EventHandler) == IN_HandleEvent) {
            // During 'HandleEvent': '<S1>:66'
            // '<S1>:69:1' sf_internal_predicateOutput = evDone;
            if (rtDW.evDone) {
              // Transition: '<S1>:69'
              rtDW.is_EventHandler = IN_RequestEvent;

              // Entry 'RequestEvent': '<S1>:65'
              // '<S1>:65:3' evDone = false;
              rtDW.evDone = false;

              // '<S1>:65:4' if waypt == 1
              if (rtDW.waypt == 1UL) {
                //  hold curr pos
                // '<S1>:65:6' traj(:, waypt) = y;
                rtDW.traj[0] = rtU.y[0];
                rtDW.traj[1] = rtU.y[1];
                rtDW.traj[2] = rtU.y[2];
              } else {
                // '<S1>:65:7' else
                //  hold last waypoint pos
                // '<S1>:65:9' traj(:,1) = traj(:, waypt);
                rtDW.traj[0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                  3];
                rtDW.traj[1] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                  3 + 1];
                rtDW.traj[2] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) *
                  3 + 2];

                // '<S1>:65:10' waypt = 1;
                rtDW.waypt = 1U;
              }

              // '<S1>:65:12' requestEvent = true;
              rtY.requestEvent = true;

              //  request new event
            } else {
              // Outport: '<Root>/inTransRegion'
              // '<S1>:66:7' [inTransRegion, evDone, waypt, holdT]...
              // '<S1>:66:8'     = handleEvent(currEv);
              handleEvent(rtY.currEv, &rtY.inTransRegion, &rtDW.evDone, &waypt,
                          &ysim0);
              rtDW.waypt = waypt;
              rtDW.holdT = ysim0;
            }

            // During 'RequestEvent': '<S1>:65'
            // '<S1>:68:1' sf_internal_predicateOutput = ~isequal(nextEv, nullEv); 
          } else if (!isequal(rtU.nextEv, rtP.nullEv)) {
            // Transition: '<S1>:68'
            // '<S1>:68:1' evDone = false;
            rtDW.evDone = false;

            // Exit 'RequestEvent': '<S1>:65'
            // '<S1>:65:14' requestEvent = false;
            rtY.requestEvent = false;
            rtDW.is_EventHandler = IN_HandleEvent;

            // Entry 'HandleEvent': '<S1>:66'
            // '<S1>:66:3' currEv = nextEv;
            rtY.currEv = rtU.nextEv;

            // '<S1>:66:4' [traj, trajSize] = trajGen(currEv, y);
            trajGen(rtY.currEv, rtU.y, rtDW.traj, &rtDW.trajSize);

            // '<S1>:66:5' holdT = 0;
            rtDW.holdT = 0.0;
          } else {
            // no actions
          }

          // During 'controlLaw': '<S1>:59'
          // '<S1>:59:5' if trueMeasAvail
          if (*sfEvent == static_cast<int32_T>(event_trueMeasAvail)) {
            // Outputs for Function Call SubSystem: '<S1>/State0.controlLaw.AMPC0' 
            // Outport: '<Root>/u' incorporates:
            //   Inport: '<Root>/u_o'
            //   Inport: '<Root>/y'
            //   Inport: '<Root>/y_o'

            // '<S1>:59:6' [u, yhat0] = AMPC0(traj(chs0, waypt), y(chs0), y_o(chs0), u_o, false); 
            // Simulink Function 'AMPC0': '<S1>:5'
            State0controlLawAMPC0(rtDW.traj[(static_cast<int32_T>(rtDW.chs0) + 3
              * (static_cast<int32_T>(rtDW.waypt) - 1)) - 1], rtU.y
                                  [static_cast<int32_T>(rtDW.chs0) - 1],
                                  rtU.y_o[static_cast<int32_T>(rtDW.chs0) - 1],
                                  rtU.u_o, 0.0, rtY.u, &rtDW.yhat0,
                                  &rtDW.State0controlLawAMPC0_n,
                                  &rtP.State0controlLawAMPC0_n, &rtP,
                                  &rtPrevZCX.State0controlLawAMPC0_n);

            // End of Outputs for SubSystem: '<S1>/State0.controlLaw.AMPC0'
          } else if ((*sfEvent == static_cast<int32_T>(event_simMeasAvail)) &&
                     rtY.inTransRegion) {
            int32_T i;

            // '<S1>:59:7' elseif simMeasAvail & inTransRegion
            //  use yhat instead of y for ch's w/o measurements
            // '<S1>:59:9' ysim0 = y(chs0);
            // '<S1>:59:10' ysim0(find(ysim0 == 0)) = yhat0(find(ysim0 == 0));
            i = -1;
            ysim0 = rtU.y[static_cast<int32_T>(rtDW.chs0) - 1];
            if (ysim0 == 0.0) {
              i = 0;
            }

            if (i >= 0) {
              ysim0 = rtDW.yhat0;
            }

            // Outputs for Function Call SubSystem: '<S1>/State0.controlLaw.AMPC0' 
            // Outport: '<Root>/u' incorporates:
            //   Inport: '<Root>/u_o'
            //   Inport: '<Root>/y'
            //   Inport: '<Root>/y_o'

            // '<S1>:59:11' [u, yhat0] = AMPC0(traj(chs0, waypt), ysim0, y_o(chs0), u_o, false); 
            // Simulink Function 'AMPC0': '<S1>:5'
            State0controlLawAMPC0(rtDW.traj[(static_cast<int32_T>(rtDW.chs0) + 3
              * (static_cast<int32_T>(rtDW.waypt) - 1)) - 1], ysim0, rtU.y_o[
                                  static_cast<int32_T>(rtDW.chs0) - 1], rtU.u_o,
                                  0.0, rtY.u, &rtDW.yhat0,
                                  &rtDW.State0controlLawAMPC0_n,
                                  &rtP.State0controlLawAMPC0_n, &rtP,
                                  &rtPrevZCX.State0controlLawAMPC0_n);

            // End of Outputs for SubSystem: '<S1>/State0.controlLaw.AMPC0'
          } else {
            // no actions
          }

          // '<S1>:59:13' yhat = [yhat0; 0; 0];
          rtY.yhat[0] = rtDW.yhat0;
          rtY.yhat[1] = 0.0;
          rtY.yhat[2] = 0.0;

          // '<S1>:59:14' currTraj = traj(:, waypt);
          rtY.currTraj[0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 3];
          c_size_idx_0 = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
          rtY.currTraj[1] = rtDW.traj[c_size_idx_0 + 1];
          rtY.currTraj[2] = rtDW.traj[c_size_idx_0 + 2];
        }
      }
      break;

     case IN_State1:
      State1(sfEvent);
      break;

     default:
      {
        // During 'State2': '<S1>:193'
        // '<S1>:207:1' sf_internal_predicateOutput = currEv.destState == 1 & evDone; 
        if ((rtY.currEv.destState == 1.0) && rtDW.evDone) {
          // Disable for Function Call SubSystem: '<S1>/State2.controlLaw.AMPC2' 
          // Transition: '<S1>:207'
          // Exit Internal 'State2': '<S1>:193'
          // Exit 'controlLaw': '<S1>:203'
          State2controlLawAMPC2_Disable(&rtDW.State2controlLawAMPC2_l,
            &rtP.State2controlLawAMPC2_l);

          // End of Disable for SubSystem: '<S1>/State2.controlLaw.AMPC2'
          // Exit Internal 'EventHandler': '<S1>:194'
          if (static_cast<uint32_T>(rtDW.is_EventHandler_k) == IN_RequestEvent)
          {
            // Exit 'RequestEvent': '<S1>:191'
            // '<S1>:191:14' requestEvent = false;
            rtDW.is_EventHandler_k = IN_NO_ACTIVE_CHILD;
          } else {
            rtDW.is_EventHandler_k = IN_NO_ACTIVE_CHILD;
          }

          rtDW.is_c6_SupervisoryController = IN_State1;

          // Entry 'State1': '<S1>:176'
          // '<S1>:176:3' waypt = 1;
          rtDW.waypt = 1U;

          // '<S1>:176:4' traj = zeros(3, 2400);
          (void)std::memset(&rtDW.traj[0], 0, 7200U * sizeof(real_T));

          // '<S1>:176:5' inTransRegion = true;
          rtY.inTransRegion = true;
          enter_internal_State1();
        } else {
          int32_T i;

          // During 'EventHandler': '<S1>:194'
          if (static_cast<uint32_T>(rtDW.is_EventHandler_k) == IN_HandleEvent) {
            // During 'HandleEvent': '<S1>:190'
            // '<S1>:188:1' sf_internal_predicateOutput = evDone;
            if (rtDW.evDone) {
              // Transition: '<S1>:188'
              rtDW.is_EventHandler_k = IN_RequestEvent;

              // Entry 'RequestEvent': '<S1>:191'
              // '<S1>:191:3' evDone = false;
              rtDW.evDone = false;

              // '<S1>:191:4' if waypt == 1
              if (rtDW.waypt == 1UL) {
                //  hold curr pos
                // '<S1>:191:6' traj(:, waypt) = y;
                rtDW.traj[0] = rtU.y[0];
                rtDW.traj[1] = rtU.y[1];
                rtDW.traj[2] = rtU.y[2];
              } else {
                // '<S1>:191:7' else
                //  hold last waypoint pos
                // '<S1>:191:9' traj(:,1) = traj(:, waypt);
                i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
                rtDW.traj[0] = rtDW.traj[i];
                rtDW.traj[1] = rtDW.traj[i + 1];
                rtDW.traj[2] = rtDW.traj[i + 2];

                // '<S1>:191:10' waypt = 1;
                rtDW.waypt = 1U;
              }

              // '<S1>:191:12' requestEvent = true;
              rtY.requestEvent = true;

              //  request new event
            } else {
              // Outport: '<Root>/inTransRegion'
              // '<S1>:190:7' [inTransRegion, evDone, waypt, holdT]...
              // '<S1>:190:8'     = handleEvent(currEv);
              handleEvent(rtY.currEv, &rtY.inTransRegion, &rtDW.evDone, &waypt,
                          &ysim0);
              rtDW.waypt = waypt;
              rtDW.holdT = ysim0;
            }

            // During 'RequestEvent': '<S1>:191'
            // '<S1>:192:1' sf_internal_predicateOutput = ~isequal(nextEv, nullEv); 
          } else if (!isequal(rtU.nextEv, rtP.nullEv)) {
            // Transition: '<S1>:192'
            // '<S1>:192:1' evDone = false;
            rtDW.evDone = false;

            // Exit 'RequestEvent': '<S1>:191'
            // '<S1>:191:14' requestEvent = false;
            rtY.requestEvent = false;
            rtDW.is_EventHandler_k = IN_HandleEvent;

            // Entry 'HandleEvent': '<S1>:190'
            // '<S1>:190:3' currEv = nextEv;
            rtY.currEv = rtU.nextEv;

            // '<S1>:190:4' [traj, trajSize] = trajGen(currEv, y);
            trajGen(rtY.currEv, rtU.y, rtDW.traj, &rtDW.trajSize);

            // '<S1>:190:5' holdT = 0;
            rtDW.holdT = 0.0;
          } else {
            // no actions
          }

          // During 'controlLaw': '<S1>:203'
          // '<S1>:203:5' if trueMeasAvail
          if (*sfEvent == static_cast<int32_T>(event_trueMeasAvail)) {
            // '<S1>:203:6' [u, yhat2] = AMPC2(traj(chs2, waypt), y(chs2), y_o(chs2), u_o, false); 
            i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
            r_ki[0] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs2[0])) - 1];
            y_j[0] = rtU.y[static_cast<int32_T>(rtDW.chs2[0]) - 1];
            y0_i[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs2[0]) - 1];
            r_ki[1] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs2[1])) - 1];
            y_j[1] = rtU.y[static_cast<int32_T>(rtDW.chs2[1]) - 1];
            y0_i[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs2[1]) - 1];

            // Outputs for Function Call SubSystem: '<S1>/State2.controlLaw.AMPC2' 
            // Outport: '<Root>/u' incorporates:
            //   Inport: '<Root>/u_o'
            //   Inport: '<Root>/y'
            //   Inport: '<Root>/y_o'

            // Simulink Function 'AMPC2': '<S1>:202'
            State2controlLawAMPC2(r_ki, y_j, y0_i, rtU.u_o, 0.0, rtY.u,
                                  rtDW.yhat2, &rtDW.State2controlLawAMPC2_l,
                                  &rtP.State2controlLawAMPC2_l, &rtP,
                                  &rtPrevZCX.State2controlLawAMPC2_l);

            // End of Outputs for SubSystem: '<S1>/State2.controlLaw.AMPC2'
          } else if ((*sfEvent == static_cast<int32_T>(event_simMeasAvail)) &&
                     rtY.inTransRegion) {
            real_T y_j_tmp;
            int32_T c_size_idx_0;

            // '<S1>:203:7' elseif simMeasAvail & inTransRegion
            //  use yhat instead of y for ch's w/o measurements
            // '<S1>:203:9' ysim2 = y(chs2);
            // '<S1>:203:10' ysim2(find(ysim2 == 0)) = yhat1(find(ysim2 == 0));
            ysim0 = rtU.y[static_cast<int32_T>(rtDW.chs2[0]) - 1];
            y_j[0] = ysim0;
            y_j_tmp = rtU.y[static_cast<int32_T>(rtDW.chs2[1]) - 1];
            y_j[1] = y_j_tmp;
            i = 0;
            if (ysim0 == 0.0) {
              b_data[0] = 1;
              i = 1;
            }

            if (y_j_tmp == 0.0) {
              b_data[i] = 2;
            }

            i = 0;
            if (ysim0 == 0.0) {
              i = 1;
            }

            if (y_j_tmp == 0.0) {
              i++;
            }

            c_size_idx_0 = i;
            i = 0;
            if (ysim0 == 0.0) {
              r_ki[0] = rtDW.yhat1[0];
              i = 1;
            }

            if (y_j_tmp == 0.0) {
              r_ki[i] = rtDW.yhat1[1];
            }

            for (i = 0; i < c_size_idx_0; i++) {
              y_j[b_data[i] - 1] = r_ki[i];
            }

            // '<S1>:203:11' [u, yhat2] = AMPC2(traj(chs2, waypt), ysim2, y_o(chs1), u_o, false); 
            i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
            r_ki[0] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs2[0])) - 1];
            y0_i[0] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[0]) - 1];
            r_ki[1] = rtDW.traj[(i + static_cast<int32_T>(rtDW.chs2[1])) - 1];
            y0_i[1] = rtU.y_o[static_cast<int32_T>(rtDW.chs1[1]) - 1];

            // Outputs for Function Call SubSystem: '<S1>/State2.controlLaw.AMPC2' 
            // Outport: '<Root>/u' incorporates:
            //   Inport: '<Root>/u_o'
            //   Inport: '<Root>/y'
            //   Inport: '<Root>/y_o'

            // Simulink Function 'AMPC2': '<S1>:202'
            State2controlLawAMPC2(r_ki, y_j, y0_i, rtU.u_o, 0.0, rtY.u,
                                  rtDW.yhat2, &rtDW.State2controlLawAMPC2_l,
                                  &rtP.State2controlLawAMPC2_l, &rtP,
                                  &rtPrevZCX.State2controlLawAMPC2_l);

            // End of Outputs for SubSystem: '<S1>/State2.controlLaw.AMPC2'
          } else {
            // no actions
          }

          // '<S1>:203:13' yhat = [yhat2(1); 0; yhat2(2)];
          rtY.yhat[0] = rtDW.yhat2[0];
          rtY.yhat[1] = 0.0;
          rtY.yhat[2] = rtDW.yhat2[1];

          // '<S1>:203:14' currTraj = traj(:, waypt);
          rtY.currTraj[0] = rtDW.traj[(static_cast<int32_T>(rtDW.waypt) - 1) * 3];
          i = (static_cast<int32_T>(rtDW.waypt) - 1) * 3;
          rtY.currTraj[1] = rtDW.traj[i + 1];
          rtY.currTraj[2] = rtDW.traj[i + 2];
        }
      }
      break;
    }
  }

  // End of Chart: '<Root>/SupervisoryController'
}

// Model step function
void SupervisoryController::step()
{
  int32_T sfEvent;
  ZCEventType zcEvent_idx_0;
  ZCEventType zcEvent_idx_1;

  // Chart: '<Root>/SupervisoryController' incorporates:
  //   TriggerPort: '<S1>/input events'

  // Inport: '<Root>/input events'
  zcEvent_idx_0 = rt_ZCFcn(ANY_ZERO_CROSSING,
    &rtPrevZCX.SupervisoryController_Trig_ZCE[0],
    (rtU.inputevents[0]));
  zcEvent_idx_1 = rt_ZCFcn(ANY_ZERO_CROSSING,
    &rtPrevZCX.SupervisoryController_Trig_ZCE[1],
    (rtU.inputevents[1]));
  if ((zcEvent_idx_0 != NO_ZCEVENT) || (zcEvent_idx_1 != NO_ZCEVENT)) {
    // Gateway: SupervisoryController
    if (static_cast<int32_T>(zcEvent_idx_0) != 0) {
      // Event: '<S1>:16'
      sfEvent = static_cast<int32_T>(event_trueMeasAvail);
      chartstep_c6_SupervisoryControl(&sfEvent);
    }

    if (static_cast<int32_T>(zcEvent_idx_1) != 0) {
      // Event: '<S1>:28'
      sfEvent = static_cast<int32_T>(event_simMeasAvail);
      chartstep_c6_SupervisoryControl(&sfEvent);
    }
  }
}

// Model initialize function
void SupervisoryController::initialize()
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  {
    real_T Sum_b[3];
    real_T Sum3_m[2];
    real_T Sum3_l;
    rtPrevZCX.SupervisoryController_Trig_ZCE[0] = UNINITIALIZED_ZCSIG;
    rtPrevZCX.SupervisoryController_Trig_ZCE[1] = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State2controlLawAMPC2_l.MemoryX_Reset_ZCE = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State2controlLawAMPC2_l.MemoryP_Reset_ZCE = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State1controlLawAMPC1_g.MemoryX_Reset_ZCE_b = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State1controlLawAMPC1_g.MemoryP_Reset_ZCE_k = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State0controlLawAMPC0_n.MemoryX_Reset_ZCE_j = UNINITIALIZED_ZCSIG;
    rtPrevZCX.State0controlLawAMPC0_n.MemoryP_Reset_ZCE_a = UNINITIALIZED_ZCSIG;
    rtDW.chs0 = 1U;
    rtDW.chs1[0] = 2U;
    rtDW.chs2[0] = 1U;
    rtDW.chs1[1] = 3U;
    rtDW.chs2[1] = 3U;

    // SystemInitialize for Outport: '<Root>/currEv'
    rtY.currEv.srcState = 0.0;
    rtY.currEv.destState = 0.0;
    rtY.currEv.moveTime = 0.0;
    rtY.currEv.holdTime = 0.0;
    rtY.currEv.destPos[0] = 0.0;
    rtY.currEv.chs[0] = false;
    rtY.currEv.nextChs[0] = false;
    rtY.currEv.destPos[1] = 0.0;
    rtY.currEv.chs[1] = false;
    rtY.currEv.nextChs[1] = false;
    rtY.currEv.destPos[2] = 0.0;
    rtY.currEv.chs[2] = false;
    rtY.currEv.nextChs[2] = false;

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/State0.controlLaw.AMPC0'

    State0controlLawAMPC0_Init(Sum_b, &Sum3_l, &rtDW.State0controlLawAMPC0_n,
      &rtP.State0controlLawAMPC0_n);

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/State1.controlLaw.AMPC1'

    State1controlLawAMPC1_Init(Sum_b, Sum3_m, &rtDW.State1controlLawAMPC1_g,
      &rtP.State1controlLawAMPC1_g);

    // SystemInitialize for Chart: '<Root>/SupervisoryController' incorporates:
    //   SubSystem: '<S1>/State2.controlLaw.AMPC2'

    State2controlLawAMPC2_Init(Sum_b, Sum3_m, &rtDW.State2controlLawAMPC2_l,
      &rtP.State2controlLawAMPC2_l);
  }
}

// Model terminate function
void SupervisoryController::terminate()
{
  // (no terminate code required)
}

// Constructor
SupervisoryController::SupervisoryController():
  rtU(),
  rtY(),
  rtDW(),
  rtPrevZCX()
{
  // Currently there is no constructor body generated.
}

// Destructor
SupervisoryController::~SupervisoryController()
{
  // Currently there is no destructor body generated.
}

//
// File trailer for generated code.
//
// [EOF]
//
