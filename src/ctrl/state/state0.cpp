#include "ctrl/supervisor.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state1.hpp"

State0::State0(Supervisor *sv)
    : State(sv), ch(Eigen::Matrix<int, 1, 1>(0)),
      // system matrices
      Ad(openData(sv->getDataPath() + "state0/Ad.txt")),
      Ad_(openData(sv->getDataPath() + "state0/Ad_.txt")),
      Bd(openData(sv->getDataPath() + "state0/Bd.txt")),
      Cd(openData(sv->getDataPath() + "state0/Cd.txt")),
      Cd_(openData(sv->getDataPath() + "state0/Cd_.txt")), CdInv(Cd.inverse()),
      K1(openData(sv->getDataPath() + "state0/K1.txt")),
      K2(openData(sv->getDataPath() + "state0/K2.txt")),
      Qw(openData(sv->getDataPath() + "state0/Qw.txt")),
      Rv(openData(sv->getDataPath() + "state0/Rv.txt")), P0(Vector1d::Identity(1, 1)), P(P0),

      // initial conditions
      du(Eigen::Vector3d::Zero()), uref(Eigen::Vector3d(60, 40, 60)),

      z0(Eigen::Vector3d::Zero()), z(z0),
      initTime(steady_clock::now()), prevCtrlTime{initTime, initTime, initTime}, dt{0s, 0s, 0s},

      yrefScale(Eigen::Vector3d(sv->imProc->impConf.getChanBBox()[0].height, 0, 0)),
      yref0(Eigen::Vector3d(1, 0, 0)), yref((yref0.array() * yrefScale.array()).matrix()),
      dyref(yref - yref0),

      dxhat(Eigen::Vector3d::Zero()), dyhat(Cd * dxhat) {}

State0::~State0() {
    // clean up any resources used by current state here
}

// check for new measurements on selected channels
bool State0::measurementAvailable() {
  bool measAvail = true, trueMeasAvail = true;
  for (int i = 0; i != ch.rows(); ++i) {
    if (!openLoop)
      measAvail &= !sv_->imProc->procDataQArr[ch(i)]->empty();
    else {
      measAvail &=
          duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[ch(i)]).count() >= 25;
      trueMeasAvail &= !sv_->imProc->procDataQArr[ch(i)]->empty();
    }
  }

  if (trueMeasAvail)
    openLoop = false;

  return measAvail;
}

// update instantaneous trajectory vectors
void State0::updateMeasurement() {
  for (int i = 0; i != ch.rows(); ++i) {
    // update measurement vectors dy, y
    if (!openLoop)
      dy(ch(i)) = sv_->imProc->procDataQArr[ch(i)]->get().y - yref(ch(i));
    else
      dy(ch(i)) = dyhat(ch(i));
    y(ch(i)) = dy(ch(i)) + yref(ch(i));

    // update time step dt for numerical integration
    if (!firstMeasAvail[ch(i)])
      firstMeasAvail[ch(i)] = true;
    else // measure the time difference between consecutive measurements
      dt[ch(i)] = steady_clock::now() - prevCtrlTime[ch(i)];
    prevCtrlTime[ch(i)] = steady_clock::now();

    // update integral error based on time step for each channel's data
    z(ch(i)) += -dy(ch(i)) * dt[ch(i)].count();
  }
}

void State0::handleEvent(Event *event) {
  if (event->srcState != 0) {
    info("Invalid event! source state should be 0, but is actually {}", event->srcState);
    return;
  }

  Eigen::Vector3d yDest = (event->destPos.array() * yrefScale.array()).matrix();
  bool destReached = true;

  // generate next waypoint if destination is not reached
  // TODO figure out what units dt is measured in (is it milliseconds or seconds?)
  for (int i = 0; i != ch.rows(); ++i) {
    if (!firstMeasAvail[ch(i)]) {
      if (y(ch(i)) < yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) + event->vel(ch(i)) * 25e-3;
      else if (y(ch(i)) > yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) - event->vel(ch(i)) * 25e-3;
    } else {
      if (y(ch(i)) < yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) + event->vel(ch(i)) * dt[ch(i)].count() * 1e-3;
      else if (y(ch(i)) > yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) - event->vel(ch(i)) * dt[ch(i)].count() * 1e-3;
    }
    destReached &= std::abs(yref(ch(i)) - yDest(ch(i))) < event->vel(ch(i)) * 25e-3;
  }

  if (destReached) {
    yref = yDest;
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
    if (event->destState == 1) {
      // clear all improc queues
      sv_->imProc->clearProcDataQueues();
      sv_->updateState<State1>();
    }
  }
}

Eigen::Matrix<int16_t, 3, 1> State0::step() {
  // update state x and control signal u based on new measurements
  // Kalman observer
  // prediction
  // dxhat = dxhat_prev - dxref, Cd*dxref = dyref
  dxhat = Ad * dxhat(ch.array(), Eigen::all) + Bd * du(ch.array(), Eigen::all) - CdInv * dyref(ch.array(), Eigen::all);
  P = Ad * P * Ad_ + Qw;
  // correction
  temp = Cd * P * Cd_ + Rv;
  tempInv = temp.inverse();
  Ko = P * Cd_ * tempInv;
  dyhat = Cd * dxhat;
  dxhat = dxhat + Ko * (dy - dyhat);

  // apply updated control signals (with saturation limits) to pump
  // LQI control law
  du(ch.array(), Eigen::all) = -K1 * dxhat - K2 * z;
  u = uref + du;

  // apply saturation (+/- 20)
  for (int i = 0; i != du.rows(); ++i) {
    if (du(i) < -20) du(i) = -20;
    else if (du(i) > 20) du(i) = 20;
  }
  usat = uref + du;

  return usat.cast<int16_t>();

}
