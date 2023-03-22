#pragma once

#include "ctrl/state/state.hpp"

extern "C" {
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#include "mpcmoveCodeGeneration.h"
#include "mpcmoveCodeGeneration_terminate.h"
#include "mpcmoveCodeGeneration_types.h"
#include "rt_nonfinite.h"
}

void argInit_3x1_real_T(double result[3]);
void argInit_3x3_real_T(double result[9]);
void argInit_62x1_boolean_T(boolean_T result[62]);
boolean_T argInit_boolean_T(void);
double argInit_real_T(void);
void argInit_struct4_T(struct4_T *result);
struct5_T argInit_struct5_T(void);
struct6_T argInit_struct6_T(void);
struct8_T argInit_struct8_T(void);

class State0 : public State {
  // num of states, outputs, inputs
  static constexpr int state = 0, numX = 3, numY = 1, numU = 1;

public:
  // system matrices
  Vector1ui ch{Vector1ui(0)};
  // Eigen::Matrix<double, numX, 1> dxhat{Eigen::Matrix<double, numX, 1>::Zero()};
  MDL<state, numX, numY, numU> *mdl;

  struct10_T info_;
  struct4_T stateData;
  struct5_T onlineData;

  State0(Supervisor *sv, Eigen::Vector3d uref_ = Eigen::Vector3d(84, 41, 40));
  ~State0();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
