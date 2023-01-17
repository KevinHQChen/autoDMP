#pragma once

#include "ctrl/state/state.hpp"

namespace py = pybind11;

using namespace py::literals;

class SysIDState : public State {
  bool *selChs_;
  float *minVals_;
  float *maxVals_;
  unsigned int numSamples_;
  int stp = 0;
  py::object incFcn;

public:
  bool simMeasAvail_ = true, trueMeasAvail_ = true;
  SysIDState(Supervisor *sv, Eigen::Vector3d uref);
  SysIDState(Supervisor *sv, Eigen::Vector3d uref, bool *selChs, float *minVals, float *maxVals,
             unsigned int numSamples);
  ~SysIDState();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
