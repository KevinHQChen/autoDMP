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

struct Event;
class Supervisor; // forward declaration

class State {
protected:
  Supervisor *sv_;

public:
  State(Supervisor *sv);
  virtual ~State();

  virtual bool measurementAvailable() = 0;
  virtual void updateMeasurement() = 0;

  virtual void handleEvent(Event *event) = 0;

  virtual Eigen::Matrix<int16_t, 3, 1> step() = 0; // pure virtual - makes State an abstract class
};
