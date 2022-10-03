#pragma once

#include "ctrl/state/state.hpp"

class State1 : public State {
public:
  // system matrices
  Vector2ui ch;
  Eigen::Matrix2d Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
  // dynamic (changes based on observer estimation error)
  Eigen::Matrix2d P0, P, Ko, temp, tempInv;

  State1(Supervisor *supervisor);
  ~State1();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;

  /*
   * @brief: performs the required tasks for state transition corresponding to the received event
   * @param: event
   */
  virtual void handleEvent(Event *event) override;

  /*
   * @brief: generates control signals at each time step
   */
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
