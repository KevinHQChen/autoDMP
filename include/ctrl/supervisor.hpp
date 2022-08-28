#pragma once

#include "util/util.hpp"

class State; // forward declaration

class Supervisor {
  State *currentState_;
  Eigen::Matrix currPos;
  Eigen::Matrix u, usat;

public:
  Supervisor();
  ~Supervisor();

  bool procDataAvail();
  Eigen::Matrix<int16_t, 3, 1> step();

  template <typename T> T getPos() {}

  template <typename T> void updateState() {
    delete currentState_;
    currentState = new T(this);
  }
};
