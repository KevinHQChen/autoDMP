#pragma once

#include "ctrl/state/state.hpp"

class State0 : public State {
  static Model model;
  static Eigen::MatrixXd yref_traj;
  double yref;

  int rotAngle;
  cv::Rect chanBBox;

public:
  State0(Supervisor *supervisor);
  ~State0();
  virtual Eigen::Matrix<int16_t, 3, 1> step() override;
};
