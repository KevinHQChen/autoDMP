#include "ctrl/state/state.hpp"

State::State(Supervisor *sv, Eigen::Vector3d uref_, Eigen::Vector3d yrefScale)
    : sv_(sv), uref(uref_), u(uref), usat(uref), yrefScale(yrefScale), yref0(yrefScale),
      yref(yref0) {}

State::~State() {
  // virtual base class destructor
  // no impl needed
}
