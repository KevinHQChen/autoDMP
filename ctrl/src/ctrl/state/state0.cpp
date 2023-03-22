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
  argInit_struct4_T(&stateData);
  onlineData = argInit_struct5_T();

  // TODO update onlineData.limits.ymax with yrefScale
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
  // if (event->destPos.isZero() && event->vel.isZero()) {
  //   info("resetting state 0");
  //   sv_->updateState<State0>(uref);
  //   sv_->addEvent(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0));
  //   return;
  // }
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
  // Update measured output in online data.
  onlineData.signals.ym = dy(ch.array(), Eigen::all);

  // Update reference in online data.
  onlineData.signals.ref = yref(ch.array(), Eigen::all);

  // Update constraints in online data using current uref
  onlineData.limits.umin = (uref(0) > 40) ? -40 : std::fmax(-uref(0), 0);
  onlineData.limits.umax = (uref(0) < 250-40) ? 40 : std::fmin(250-uref,250);
  onlineData.limits.ymin = -y(ch.array(), Eigen::all);
  onlineData.limits.ymax = yrefScale(ch.array(), Eigen::all) - y(ch.array(), Eigen::all);

  // Compute and store control action.
  double du_;
  mpcmoveCodeGeneration(&stateData, &onlineData, &du_, &info);
  du(0) = du_;
  u = uref + du;
  return u.cast<int16_t>();
}
