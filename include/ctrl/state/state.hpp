#pragma once

#include "util/util.hpp"

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
