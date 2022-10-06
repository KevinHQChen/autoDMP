#include "ctrl/state/state.hpp"

State::State(Supervisor *sv, Eigen::Vector3d uref, Eigen::Vector3d yrefScale)
    : sv_(sv), uref(uref), yrefScale(yrefScale), yref0(yrefScale) {
  yref = yref0;
  dyref = Eigen::Vector3d::Zero();
}

State::~State() {
  // virtual base class destructor
  // no impl needed
}
