#include "ctrl/state/state.hpp"

State::State(Supervisor* sv) : sv_(sv) { }

State::~State() {
    // virtual base class destructor
    // no impl needed
}
