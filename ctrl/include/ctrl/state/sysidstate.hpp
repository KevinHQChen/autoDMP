#pragma once

#include "ctrl/state/state.hpp"

namespace py = pybind11;
using namespace py::literals;

class SysIDState : public State {
  bool *selChs_;
  std::vector<float> minVals_;
  std::vector<float> maxVals_;
  int stp = 0;
  py::object simMeas;
  Eigen::MatrixXd excitationSignal_;

public:
  bool simMeasAvail_ = true, trueMeasAvail_ = true;
  SysIDState(Supervisor *sv, Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
             std::vector<float> maxVals, Eigen::MatrixXd &data);
  ~SysIDState();

  virtual bool measurementAvailable() override;
  virtual void updateMeasurement() override;
  virtual void handleEvent(Event *event) override;
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
