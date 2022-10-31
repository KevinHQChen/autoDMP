#pragma once

#include "ctrl/state/state.hpp"

class State0 : public State {
  // num of states, outputs, inputs
  static constexpr int state = 0, numX = 1, numY = 1, numU = 1;

public:
  // system matrices
  Vector1ui ch{Vector1ui(0)};
  // Eigen::Matrix<double, numX, 1> dxhat{Eigen::Matrix<double, numX, 1>::Zero()};
  MDL<state, numX, numY, numU> *mdl;

  State0(Supervisor *sv, Eigen::Vector3d uref_);
  ~State0();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
