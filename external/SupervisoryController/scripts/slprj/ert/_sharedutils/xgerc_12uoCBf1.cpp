//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: xgerc_12uoCBf1.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "xgerc_12uoCBf1.h"

// Function for MATLAB Function: '<S36>/FixedHorizonOptimizer'
void xgerc_12uoCBf1(int32_T b_m, int32_T n, real_T alpha1, int32_T ix0, const
                    real_T y[4], real_T b_A[16], int32_T ia0)
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

//
// File trailer for generated code.
//
// [EOF]
//
