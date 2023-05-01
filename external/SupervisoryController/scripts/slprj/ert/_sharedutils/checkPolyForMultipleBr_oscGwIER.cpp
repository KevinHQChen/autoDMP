//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: checkPolyForMultipleBr_oscGwIER.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "checkPolyForMultipleBr_oscGwIER.h"
#include <cmath>

// Function for Chart: '<Root>/SupervisoryController'
boolean_T checkPolyForMultipleBr_oscGwIER(const real_T breakMat_data[], const
  int32_T breakMat_size[2])
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

//
// File trailer for generated code.
//
// [EOF]
//
