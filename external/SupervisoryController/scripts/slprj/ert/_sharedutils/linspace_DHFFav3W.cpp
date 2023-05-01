//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: linspace_DHFFav3W.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "linspace_DHFFav3W.h"
#include <cmath>

// Function for Chart: '<Root>/SupervisoryController'
void linspace_DHFFav3W(real_T d2, uint16_T n, real_T y_data[], int32_T y_size[2])
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

//
// File trailer for generated code.
//
// [EOF]
//
