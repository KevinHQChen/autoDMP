//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: qrFactor_6LFrwaYF.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "qrFactor_6LFrwaYF.h"
#include "xnrm2_NhwXTuci.h"
#include "rt_hypotd_snf.h"
#include <cmath>

// Function for MATLAB Function: '<S38>/RLS'
real_T qrFactor_6LFrwaYF(const real_T A[4], const real_T S[16], real_T Ns)
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
  s = xnrm2_NhwXTuci(4, b_A, 2);
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

      s = rt_hypotd_snf(b_S, xnrm2_NhwXTuci(4, b_A, 2));
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

//
// File trailer for generated code.
//
// [EOF]
//
