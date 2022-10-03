#include "ctrl/state/state.hpp"

State::State(Supervisor *sv, Eigen::Vector3d uref, Eigen::Vector3d yrefScale, Eigen::Vector3d yref0)
    : sv_(sv), uref(uref), yrefScale(yrefScale), yref0(yref0) {
  yref = (yref0.array() * yrefScale.array()).matrix();
  dyref = yref - yref0;
}

State::~State() {
  // virtual base class destructor
  // no impl needed
}
