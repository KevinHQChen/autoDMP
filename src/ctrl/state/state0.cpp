#include "ctrl/state/state0.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state1.hpp"
#include "ctrl/supervisor.hpp"

State0::State0(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_, Eigen::Vector3d(sv->imProc->impConf.getChanBBox()[0].height, 0, 0)),
      mdl(new MDL<state, numX, numY, numU>(sv)) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
  stateTransitionCondition = true;
}

State0::~State0() {
  // clean up any resources used by current state here
  delete mdl;
}

bool State0::measurementAvailable() { return State::measurementAvailable<numY>(ch); }

void State0::updateMeasurement() { State::updateMeasurement<numY>(ch, -1); }

void State0::handleEvent(Event *event) {
  if (event->srcState != 0) {
    info("Invalid event! source state should be 0, but is actually {}", event->srcState);
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
    return;
  }
  if (!startEvent) {
    yDest = (event->destPos.array() * yrefScale.array()).matrix();
    startEvent = true;
    yref = y;
  }
  bool destReached = true;
  bool junctionReached = true;

  // generate next waypoint if destination is not reached
  for (int i = 0; i != ch.rows(); ++i) {
    if (y(ch(i)) < yDest(ch(i)))
      dyref(ch(i)) = event->vel(ch(i)) * dt[ch(i)].count();
    else if (y(ch(i)) > yDest(ch(i)))
      dyref(ch(i)) = -event->vel(ch(i)) * dt[ch(i)].count();
    yref(ch(i)) += dyref(ch(i));

    if (std::abs(yref(ch(i)) - yDest(ch(i))) < event->vel(ch(i)) * 50e-3) {
      dyref(ch(i)) = 0;
      yref(ch(i)) = yDest(ch(i));
    }

    destReached &= std::abs(y(ch(i)) - yDest(ch(i))) < 1; // event->vel(ch(i)) * 25e-3;

    junctionReached &= y(ch(i)) < 0.9 * yrefScale(ch(i));
  }

  if (stateTransitionCondition && junctionReached)
    stateTransitionCondition = false;

  // remain in State 0
  //   |  |        |  |
  //   |  |   =>   |  |
  //  / /\_\      / /\_\
  // / /  \ \    / /  \ \.
  if (event->destState == 0) {
    if (destReached) {
      yref = yDest;
      startEvent = false;
      // z = Eigen::Vector3d::Zero();
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
    }
  }

  // transition to State 1
  //   |__|        |  |
  //   |  |   =>   |  |
  //  / /\ \      /_/\_\
  // / /  \ \    / /  \ \.
  if (event->destState == 1) {
    // start state transition when ch0 is 90% to junction
    if (yref(0) > 0.9 * yrefScale(0) && !stateTransitionCondition) {
      sv_->imProc->clearProcDataQueues();
      stateTransitionCondition = true;
    }
    // we're in sim mode and ch0 is 95% to junction
    // or ch1 or ch2 are observable and we're not observing ch0
    if ((yref(0) > 0.95 * yrefScale(0) && sv_->simModeActive) ||
        (stateTransitionCondition &&
         (!sv_->imProc->procDataQArr[1]->empty() || !sv_->imProc->procDataQArr[2]->empty()))) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State1>(usat);
    }
  }
}

Eigen::Matrix<int16_t, 3, 1> State0::step() {
  return State::step<state, numX, numY, numU>(ch, mdl);
}
