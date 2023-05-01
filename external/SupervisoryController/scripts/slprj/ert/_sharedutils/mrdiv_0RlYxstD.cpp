//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: mrdiv_0RlYxstD.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "mrdiv_0RlYxstD.h"
#include <cmath>

// Function for MATLAB Function: '<S210>/Discrete-Time KF - Calculate PLMZ'
void mrdiv_0RlYxstD(const real_T A[8], const real_T B_4[4], real_T Y[8])
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
// File trailer for generated code.
//
// [EOF]
//
