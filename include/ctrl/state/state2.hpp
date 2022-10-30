#pragma once

#include "ctrl/state/state.hpp"

#define NUM_X 4 // num states
#define NUM_Y 2 // num outputs
#define NUM_U 2 // num inputs

class State2 : public State {
public:
  // system matrices
  Vector2ui ch;
  Eigen::Matrix<double, NUM_X, NUM_X> Ad, Ad_, Qw;
  Eigen::Matrix<double, NUM_X, NUM_Y> Bd;
  Eigen::Matrix<double, NUM_Y, NUM_X> Cd, Cd_, CdInv; // TODO remove CdInv
  Eigen::Matrix<double, NUM_U, NUM_X> K1;
  Eigen::Matrix<double, NUM_Y, NUM_Y> K2, Rv;
  // dynamic (changes based on observer estimation error)
  Eigen::Matrix<double, NUM_X, NUM_X> P0, P;
  Eigen::Matrix<double, NUM_Y, NUM_X> Ko;
  Eigen::Matrix<double, NUM_Y, NUM_Y> temp, tempInv;

  bool settled;
  int settlingTime;

  State2(Supervisor *sv, Eigen::Vector3d uref_);
  ~State2();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
