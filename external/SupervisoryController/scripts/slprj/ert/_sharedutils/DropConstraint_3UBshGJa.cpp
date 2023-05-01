//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: DropConstraint_3UBshGJa.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "DropConstraint_3UBshGJa.h"

// Function for MATLAB Function: '<S141>/FixedHorizonOptimizer'
void DropConstraint_3UBshGJa(int32_T kDrop, boolean_T iA[86], int32_T *nA,
  int32_T iC[86])
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

//
// File trailer for generated code.
//
// [EOF]
//
