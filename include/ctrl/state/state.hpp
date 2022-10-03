#pragma once

#include "util/util.hpp"

// pump layout
//   top-down view
//     P3      P4
//     oil     water
//     | ch2   | ch1
//      \     /
//        \ /
//         |
//         | ch3
//         |
//        / \
//      /     \
//     |       |
//     P2      P1
//     outlet  outlet

//   camera view
//      y(2)
//      u(2)
//      outlet
//      P1 + P2
//         |
//         | ch3
//         |
//        / \
//      /     \
//     | ch2   | ch1
//     oil     water
//     P3      P4
//     u(1)    u(0)
//     y(1)    y(0)

struct Event;
class Supervisor; // forward declaration

class State {
protected:
  Supervisor *sv_;

public:
  // instantaneous trajectory vectors
  // control signal vectors
  Eigen::Vector3d du{Eigen::Vector3d::Zero()}, u, usat, uref;
  // integral error vectors: z = int (r - y) dt
  Eigen::Vector3d z0{Eigen::Vector3d::Zero()}, z{z0};
  time_point<steady_clock> initTime{steady_clock::now()};
  time_point<steady_clock> prevCtrlTime[3] = {initTime, initTime, initTime};
  duration<double> dt[3] = {0s, 0s, 0s};
  // state/output vectors
  Eigen::Vector3d yrefScale;
  Eigen::Vector3d y, yref0, yref;
  Eigen::Vector3d dy, dyref;
  Eigen::Vector3d yhat, dxhat{Eigen::Vector3d::Zero()}, dyhat{Eigen::Vector3d::Zero()}, ytilde;

  bool firstMeasAvail[3] = {false, false, false};
  bool measAvail[3] = {true, true, true};
  bool trueMeasAvail[3] = {true, true, true};
  bool obsv[3] = {false, false, false};
  bool stateTransitionCondition = false;

  State(Supervisor *sv, Eigen::Vector3d uref, Eigen::Vector3d yrefScale, Eigen::Vector3d yref0);
  virtual ~State();

  virtual bool measurementAvailable() = 0;
  virtual void updateMeasurement() = 0;

  virtual void handleEvent(Event *event) = 0;

  virtual Eigen::Matrix<int16_t, 3, 1> step() = 0; // pure virtual - makes State an abstract class
};
