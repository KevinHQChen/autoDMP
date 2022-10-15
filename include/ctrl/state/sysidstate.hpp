#pragma once

#include "ctrl/state/state.hpp"

class SysIDState : public State {
public:
  bool simMeasAvail_ = true, trueMeasAvail_ = true;
  SysIDState(Supervisor *sv, Eigen::Vector3d uref);
  ~SysIDState();

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
