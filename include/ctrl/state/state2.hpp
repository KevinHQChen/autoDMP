#pragma once

#include "ctrl/state/state.hpp"

class State2 : public State {
public:
  // system matrices
  Vector2ui ch;
  Eigen::Matrix2d Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
  // dynamic (changes based on observer estimation error)
  Eigen::Matrix2d P0, P, Ko, temp, tempInv;

  State2(Supervisor *sv, Eigen::Vector3d uref_);
  ~State2();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
