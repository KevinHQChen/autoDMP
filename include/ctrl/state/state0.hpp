#pragma once

#include "ctrl/state/state.hpp"

class State0 : public State {
  // system matrices
  Eigen::Matrix<unsigned int, 1, 1> ch;
  Vector1d Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv;
  // dynamic (changes based on observer estimation error)
  Vector1d P0, P, Ko, temp, tempInv;

  // instantaneous trajectory vectors
  Eigen::Vector3d du, u, usat, uref; // control signal vectors
  Eigen::Vector3d z0, z; // integral error: z = int (r - y) dt
  time_point<steady_clock> initTime, prevCtrlTime[3];
  duration<double> dt[3];
  Eigen::Vector3d yrefScale;
  Eigen::Vector3d y, yref0, yref;
  Eigen::Vector3d dy, dyref;
  Eigen::Vector3d yhat, dxhat, dyhat, ytilde;

  bool firstMeasAvail[1] = {false};
  bool openLoop = true;

public:
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
