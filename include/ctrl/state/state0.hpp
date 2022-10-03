#pragma once

#include "ctrl/state/state.hpp"

class State0 : public State {
public:
  // system matrices
  Vector1ui ch;
  Vector1d Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
  // dynamic (changes based on observer estimation error)
  Vector1d P0, P, Ko, temp, tempInv;

  State0(Supervisor *supervisor);
  ~State0();

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
