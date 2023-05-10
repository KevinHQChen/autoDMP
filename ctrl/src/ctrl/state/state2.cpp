#include "ctrl/state/state2.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state1.hpp"
#include "ctrl/supervisor.hpp"

State2::State2(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_,
            Eigen::Vector3d(sv->imProc->impConf.getChROIs()[0].chHeight, 0,
                            sv->imProc->impConf.getChROIs()[2].chHeight)),
      mdl(new MDL<state, numX, numY, numU>(sv)) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
  stateTransitionCondition = false;
  settled = true;
  settlingTime = 40; // 40 * 25ms = 1s
}

State2::~State2() {
  // clean up any resources used by current state here
  delete mdl;
}

bool State2::measurementAvailable() { return State::measurementAvailable<numY>(ch); }

void State2::updateMeasurement() {
  if (!settled) {
    sv_->imProc->clearProcDataQueues();
    for (int i = 0; i != ch.rows(); ++i)
      prevCtrlTime[ch(i)] = steady_clock::now();
    return;
  }

  State::updateMeasurement<numY>(ch, 0);
}

void State2::handleEvent(Event *event) {
  if (event->srcState != 2) {
    info("Invalid event! source state should be 2, but is actually {}", event->srcState);
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
    return;
  }

  if (!settled)
    return;

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

  // remain in State 2
  //   |__|        |__|
  //   |  |   =>   |  |
  //  / /\_\      / /\_\
  // / /  \_\    / /  \_\.
  if (event->destState == 2) {
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
  //  / /\_\      /_/\_\
  // / /  \_\    / /  \ \.
  if (event->destState == 1) {
    // ch0 90%, ch2 85% to junction
    if (yref(0) > 0.9 * yrefScale(0) && yref(2) > 0.85 * yrefScale(2) &&
        !stateTransitionCondition) {
      sv_->imProc->clearProcDataQueues();
      stateTransitionCondition = true;
    }
    // we're in sim mode and ch0 & ch2 are 95% to junction
    // or ch1 is observable
    if ((yref(0) > 0.95 * yrefScale(1) && yref(2) > 0.95 * yrefScale(2) && sv_->simModeActive) ||
        (stateTransitionCondition && !sv_->imProc->procDataQArr[1]->empty())) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State1>(usat);
    }
  }
}

Eigen::Matrix<int16_t, 3, 1> State2::step() {
  if (!settled) {
    --settlingTime;
    if (settlingTime == 0) {
      settled = true;
      info("State 2 settled");
    }
    return uref.cast<int16_t>();
  } else
    return State::step<state, numX, numY, numU>(ch, mdl);
}
