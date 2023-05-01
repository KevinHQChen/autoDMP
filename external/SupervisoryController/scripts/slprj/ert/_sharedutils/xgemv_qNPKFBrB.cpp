//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: xgemv_qNPKFBrB.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "xgemv_qNPKFBrB.h"
#include <cstring>
#include "div_nde_s32_floor.h"

// Function for MATLAB Function: '<S38>/RLS'
void xgemv_qNPKFBrB(int32_T m, int32_T n, const real_T A[20], int32_T ia0, const
                    real_T x[20], int32_T ix0, real_T y[4])
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

//
// File trailer for generated code.
//
// [EOF]
//
