#pragma once

#include "ctrl/state/state.hpp"

extern "C" {
#include "_coder_mpcmoveCodeGeneration_api.h" // Replace with the appropriate header file(s)
}

class State0 : public State {
  // num of states, outputs, inputs
  static constexpr int state = 0, numX = 3, numY = 1, numU = 1;

public:
  // system matrices
  Vector1ui ch{Vector1ui(0)};
  // Eigen::Matrix<double, numX, 1> dxhat{Eigen::Matrix<double, numX, 1>::Zero()};
  MDL<state, numX, numY, numU> *mdl;

  struct10_T info;
  struct4_T stateData;
  struct5_T onlineData;

  State0(Supervisor *sv, Eigen::Vector3d uref_ = Eigen::Vector3d(84, 41, 40));
  ~State0();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
