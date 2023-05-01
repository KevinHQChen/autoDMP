//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: computeProfileParams_OfI5TC7i.cpp
//
// Code generated for Simulink model 'SupervisoryController'.
//
// Model version                  : 1.713
// Simulink Coder version         : 9.8 (R2022b) 13-May-2022
// C/C++ source code generated on : Mon May  1 12:50:40 2023
//
#include "rtwtypes.h"
#include "computeProfileParams_OfI5TC7i.h"
#include <cstring>
#include <cmath>

// Function for Chart: '<Root>/SupervisoryController'
void computeProfileParams_OfI5TC7i(real_T i, const real_T wayPoints_data[],
  const int32_T wayPoints_size[2], const real_T Vel_data[], const int32_T
  *Vel_size, const real_T TAc_data[], const int32_T *TAc_size, real_T *vParam,
  real_T *aParam, real_T *tAParam, real_T *tFParam)
{
  real_T TAcSwitch_data[3];
  real_T VelSwitch_data[3];
  real_T tACandidates[2];
  real_T b_sF;
  real_T s0;
  real_T sF;
  int32_T deltaSign;
  int32_T loop_ub;
  int8_T b_data[2];
  boolean_T inputCombo_idx_0;
  boolean_T inputCombo_idx_3;
  *aParam = 0.0;
  *vParam = 0.0;
  *tAParam = 0.0;
  *tFParam = 0.0;
  s0 = wayPoints_data[static_cast<int32_T>(i) - 1];
  sF = wayPoints_data[(static_cast<int32_T>(i) + wayPoints_size[0]) - 1];
  deltaSign = 1;
  if (sF < s0) {
    b_sF = s0;
    s0 = sF;
    sF = b_sF;
    deltaSign = -1;
  }

  inputCombo_idx_0 = (*Vel_size != 0);
  inputCombo_idx_3 = (*TAc_size != 0);
  if (inputCombo_idx_0) {
    if (*Vel_size - 1 >= 0) {
      (void)std::memcpy(&VelSwitch_data[0], &Vel_data[0], static_cast<uint32_T>(*
        Vel_size) * sizeof(real_T));
    }
  } else {
    loop_ub = static_cast<int32_T>(i);
    for (int32_T i_0{0}; i_0 < loop_ub; i_0++) {
      VelSwitch_data[i_0] = 1.0;
    }
  }

  if (inputCombo_idx_3) {
    if (*TAc_size - 1 >= 0) {
      (void)std::memcpy(&TAcSwitch_data[0], &TAc_data[0], static_cast<uint32_T>(*
        TAc_size) * sizeof(real_T));
    }
  } else {
    loop_ub = static_cast<int32_T>(i);
    for (int32_T i_0{0}; i_0 < loop_ub; i_0++) {
      TAcSwitch_data[i_0] = 1.0;
    }
  }

  switch ((inputCombo_idx_0 << 3UL) + (inputCombo_idx_3 ? static_cast<int32_T>(1)
           : static_cast<int32_T>(0))) {
   case 9:
    *vParam = VelSwitch_data[static_cast<int32_T>(i) - 1];
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *aParam = *vParam / *tAParam;
    *tFParam = ((*vParam * *tAParam + sF) - s0) / *vParam;
    break;

   case 8:
    *vParam = VelSwitch_data[static_cast<int32_T>(i) - 1];
    *tFParam = (sF - s0) * 1.5 / *vParam;
    *tAParam = (*vParam * *tFParam + (s0 - sF)) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 6:
    *aParam = 1.0;
    *tFParam = 1.0;
    b_sF = std::sqrt((4.0 * s0 + 1.0) - 4.0 * sF) / 2.0;
    tACandidates[0] = 0.5 - b_sF;
    tACandidates[1] = b_sF + 0.5;
    loop_ub = 0;
    if (0.5 - b_sF > 0.0) {
      b_data[0] = 1;
      loop_ub = 1;
    }

    if (b_sF + 0.5 > 0.0) {
      b_data[loop_ub] = 2;
    }

    *tAParam = tACandidates[b_data[0] - 1];
    *vParam = (sF - s0) / (1.0 - *tAParam);
    break;

   case 5:
    *aParam = 1.0;
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *vParam = *tAParam;
    *tFParam = ((*tAParam * *tAParam + sF) - s0) / *tAParam;
    break;

   case 4:
    *aParam = 1.0;
    b_sF = sF - s0;
    *tAParam = std::sqrt(b_sF / 2.0);
    *vParam = *tAParam;
    *tFParam = (*tAParam * *tAParam + b_sF) / *tAParam;
    break;

   case 3:
    *tFParam = 1.0;
    *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
    *vParam = (sF - s0) / (1.0 - *tAParam);
    *aParam = *vParam / *tAParam;
    break;

   case 2:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   case 1:
    {
      real_T aParam_tmp;
      *tAParam = TAcSwitch_data[static_cast<int32_T>(i) - 1];
      b_sF = *tAParam * *tAParam;
      aParam_tmp = sF - s0;
      *aParam = aParam_tmp / (b_sF * 2.0);
      *vParam = *tAParam * *aParam;
      *tFParam = (b_sF * *aParam + aParam_tmp) / *vParam;
    }
    break;

   case 0:
    *tFParam = 1.0;
    *vParam = (sF - s0) * 1.5;
    *tAParam = ((s0 - sF) + *vParam) / *vParam;
    *aParam = *vParam / *tAParam;
    break;

   default:
    // no actions
    break;
  }

  if (s0 == sF) {
    *aParam = 0.0;
    *vParam = 0.0;
    if (std::isnan(*tFParam)) {
      *tFParam = 1.0;
    } else if (*tFParam == 0.0) {
      *tFParam = 1.0;
    } else {
      // no actions
    }

    *tAParam = *tFParam / 3.0;
  }

  *vParam *= static_cast<real_T>(deltaSign);
  *aParam *= static_cast<real_T>(deltaSign);
}

//
// File trailer for generated code.
//
// [EOF]
//
