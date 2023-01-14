#pragma once

#include "ctrl/state/state.hpp"

class SysIDState : public State {
public:
  bool simMeasAvail_ = true, trueMeasAvail_ = true;
  SysIDState(Supervisor *sv, Eigen::Vector3d uref);
  ~SysIDState();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
