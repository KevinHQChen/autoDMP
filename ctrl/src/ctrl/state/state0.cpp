#include "ctrl/state/state0.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state1.hpp"
#include "ctrl/supervisor.hpp"

void argInit_3x1_real_T(double result[3]) {
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

void argInit_3x3_real_T(double result[9]) {
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 9; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

void argInit_62x1_boolean_T(boolean_T result[62]) {
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 62; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_boolean_T();
  }
}

boolean_T argInit_boolean_T(void) { return false; }

double argInit_real_T(void) { return 0.0; }

void argInit_struct4_T(struct4_T *result) {
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  argInit_3x1_real_T(result->Plant);
  result->LastMove = argInit_real_T();
  argInit_3x3_real_T(result->Covariance);
  argInit_62x1_boolean_T(result->iA);
}

struct5_T argInit_struct5_T(void) {
  struct5_T result;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result.signals = argInit_struct6_T();
  result.limits = argInit_struct8_T();
  return result;
}

struct6_T argInit_struct6_T(void) {
  struct6_T result;
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.ref = result_tmp;
  result.ym = result_tmp;
  return result;
}

struct8_T argInit_struct8_T(void) {
  struct8_T result;
  double result_tmp;
  /* Set the value of each structure field.
Change this value to the value that the application requires. */
  result_tmp = argInit_real_T();
  result.ymax = result_tmp;
  result.umin = result_tmp;
  result.umax = result_tmp;
  result.ymin = result_tmp;
  return result;
}

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
  onlineData.signals.ym = dy(0);

  // Update reference in online data.
  onlineData.signals.ref = yref(0);

  // Update constraints in online data using current uref
  onlineData.limits.umin = (uref(0) > 40) ? -40 : std::fmax(-uref(0), 0);
  onlineData.limits.umax = (uref(0) < 250 - 40) ? 40 : std::fmin(250 - uref(0), 250);
  onlineData.limits.ymin = -y(0);
  onlineData.limits.ymax = yrefScale(0) - y(0);

  // Compute and store control action.
  double du_;
  mpcmoveCodeGeneration(&stateData, &onlineData, &du_, &info_);
  du(0) = du_;
  u = uref + du;
  return u.cast<int16_t>();
}
