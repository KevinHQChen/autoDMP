#pragma once

#include "ctrl/state/state.hpp"

class State2 : public State {
  // num of states, outputs, inputs
  static constexpr int state = 2, numX = 4, numY = 2, numU = 2;

public:
  // system matrices
  Vector2ui ch{Vector2ui(0, 2)};
  // Eigen::Matrix<double, numX, 1> dxhat{Eigen::Matrix<double, numX, 1>::Zero()};
  MDL<state, numX, numY, numU> *mdl;

  bool settled;
  int settlingTime;

  State2(Supervisor *sv, Eigen::Vector3d uref_);
  ~State2();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
