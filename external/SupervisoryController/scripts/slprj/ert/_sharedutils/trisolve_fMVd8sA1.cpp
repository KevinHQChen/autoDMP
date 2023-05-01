//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: trisolve_fMVd8sA1.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "trisolve_fMVd8sA1.h"

// Function for MATLAB Function: '<S38>/RLS'
void trisolve_fMVd8sA1(real_T A, real_T B_1[4])
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

//
// File trailer for generated code.
//
// [EOF]
//
