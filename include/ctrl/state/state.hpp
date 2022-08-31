#pragma once

#include "util/util.hpp"

class Supervisor; // forward declaration

class State {
protected:
  Supervisor *sv_;

public:
  State(Supervisor *sv);
  virtual ~State();

  virtual Eigen::Matrix<int16_t, 3, 1> step() = 0; // pure virtual - makes State an abstract class
};
