#include "ctrl/state/state1.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state2.hpp"
#include "ctrl/supervisor.hpp"

State1::State1(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_,
            Eigen::Vector3d(0, sv->imProc->impConf.getChROIs()[1].chHeight,
                            sv->imProc->impConf.getChROIs()[2].chHeight)),
      mdl(new MDL<state, numX, numY, numU>(sv)) {
  // clear all improc queues
  sv_->imProc->clearProcData();
  stateTransitionCondition = true;
}

State1::~State1() {
  // clean up any resources used by current state here
  delete mdl;
}

bool State1::measurementAvailable() { return State::measurementAvailable<numY>(ch); }

void State1::updateMeasurement() {
  if (stateTransitionCondition)
    State::updateMeasurement<numY>(ch, 1);
  else
    State::updateMeasurement<numY>(ch, -1);
}

void State1::handleEvent(Event *event) {
  if (event->srcState != 1) {
    info("Invalid event! source state should be 1, but is actually {}", event->srcState);
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

    junctionReached &= y(ch(i)) < 0.85 * yrefScale(ch(i));
  }

  if (stateTransitionCondition && junctionReached)
    stateTransitionCondition = false;

  // remain in State 1
  //   |__|        |__|
  //   |  |   =>   |  |
  //  /_/\ \      /_/\ \
  // / /  \ \    / /  \ \.
  if (event->destState == 1) {
    if (destReached) {
      yref = yDest;
      startEvent = false;
      // z = Eigen::Vector3d::Zero();
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
    }
  }

  // transition to State 0
  //   |  |        |__|
  //   |  |   =>   |  |
  //  /_/\_\      / /\ \
  // / /  \ \    / /  \ \.
  if (event->destState == 0) {
    // ch1 & ch2 are 85% to junction
    if (yref(1) > 0.85 * yrefScale(1) && yref(2) > 0.85 * yrefScale(2) &&
        !stateTransitionCondition) {
      sv_->imProc->clearProcData();
      stateTransitionCondition = true;
    }
    // we're in sim mode and ch1 & ch2 are 95% to junction
    // or ch0 is observable
    if ((yref(1) > 0.95 * yrefScale(1) && yref(2) > 0.95 * yrefScale(2) && sv_->simModeActive) ||
        (stateTransitionCondition && !sv_->poses[0].found)) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State0>(usat);
    }
  }

  // transition to State 2
  //   |  |        |__|
  //   |  |   =>   |  |
  //  /_/\_\      / /\_\
  // / /  \ \    / /  \_\.
  if (event->destState == 2) {
    // ch1 is 85% to junction
    if (yref(1) > 0.85 * yrefScale(1) && !stateTransitionCondition) {
      sv_->imProc->clearProcData();
      stateTransitionCondition = true;
    }
    // we're in sim mode and ch1 is 95% to junction
    // or ch0 is observable
    if ((yref(1) > 0.95 * yrefScale(1) && sv_->simModeActive) ||
        (stateTransitionCondition && !sv_->poses[0].found)) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State2>(usat);
    }
  }
}

Eigen::Vector3d State1::step() {
  return State::step<state, numX, numY, numU>(ch, mdl);
}
