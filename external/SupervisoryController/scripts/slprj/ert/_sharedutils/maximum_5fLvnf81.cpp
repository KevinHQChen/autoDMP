//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: maximum_5fLvnf81.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "maximum_5fLvnf81.h"
#include <cmath>

// Function for Chart: '<Root>/SupervisoryController'
real_T maximum_5fLvnf81(const real_T x_data[], const int32_T *x_size)
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

//
// File trailer for generated code.
//
// [EOF]
//
