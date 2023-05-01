//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: sqrtMeasurementUpdate_MDG7rao2.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "sqrtMeasurementUpdate_MDG7rao2.h"
#include <cmath>
#include "qrFactor_6LFrwaYF.h"
#include "trisolve_fMVd8sA1.h"
#include "xnrm2_iCHzl0Tr.h"
#include "rt_hypotd_snf.h"
#include "xgemv_qNPKFBrB.h"
#include "xgerc_a2Ua5hDt.h"
#include <emmintrin.h>

extern "C"
{

#include "rt_nonfinite.h"

}

// Function for MATLAB Function: '<S38>/RLS'
void sqrtMeasurementUpdate_MDG7rao2(real_T L[16], const real_T H[4], real_T a0,
  real_T K[4])
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
  Sy = qrFactor_6LFrwaYF(H, L, prodVal);
  B_0[0] = Pxy[0];
  B_0[1] = Pxy[1];
  B_0[2] = Pxy[2];
  B_0[3] = Pxy[3];
  trisolve_fMVd8sA1(Sy, B_0);
  K[0] = B_0[0];
  K[1] = B_0[1];
  K[2] = B_0[2];
  K[3] = B_0[3];
  trisolve_fMVd8sA1(Sy, K);
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
    prodVal = xnrm2_iCHzl0Tr(4 - b_i, c_A, ii + 2);
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

        prodVal = rt_hypotd_snf(Sy, xnrm2_iCHzl0Tr(4 - b_i, c_A, ii + 2));
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
        xgemv_qNPKFBrB(lastv, coffset, c_A, ii + 6, c_A, ii + 1, Pxy);
        xgerc_a2Ua5hDt(lastv, coffset, -B_0[b_i], ii + 1, Pxy, c_A, ii + 6);
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
// File trailer for generated code.
//
// [EOF]
//
