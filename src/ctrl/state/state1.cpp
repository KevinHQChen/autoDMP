#include "ctrl/supervisor.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state1.hpp"

State1::State1(Supervisor *sv)
    : State(sv), ch(Vector2ui(1, 2)),
      // system matrices
      Ad(openData(sv->getDataPath() + "state1/Ad.txt")),
      Ad_(openData(sv->getDataPath() + "state1/Ad_.txt")),
      Bd(openData(sv->getDataPath() + "state1/Bd.txt")),
      Cd(openData(sv->getDataPath() + "state1/Cd.txt")),
      Cd_(openData(sv->getDataPath() + "state1/Cd_.txt")), CdInv(Cd.inverse()),
      K1(openData(sv->getDataPath() + "state1/K1.txt")),
      K2(openData(sv->getDataPath() + "state1/K2.txt")),
      Qw(openData(sv->getDataPath() + "state1/Qw.txt")),
      Rv(openData(sv->getDataPath() + "state1/Rv.txt")), P0(Eigen::Matrix2d::Identity(2, 2)), P(P0),

      // initial conditions
      du(Eigen::Vector3d::Zero()), uref(Eigen::Vector3d(65, 45, 65)),

      z0(Eigen::Vector3d::Zero()), z(z0),
      initTime(steady_clock::now()), prevCtrlTime{initTime, initTime, initTime}, dt{0s, 0s, 0s},

      yrefScale(Eigen::Vector3d(0, sv->imProc->impConf.getChanBBox()[1].height, sv->imProc->impConf.getChanBBox()[2].height)),
      yref0(Eigen::Vector3d(0, 1, 1)), yref((yref0.array() * yrefScale.array()).matrix()),
      dyref(yref - yref0),

      dxhat(Eigen::Vector3d::Zero()), dyhat(Eigen::Vector3d::Zero()) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
}

State1::~State1() {
    // clean up any resources used by current state here
}

// check for new measurements on selected channels
bool State1::measurementAvailable() {
  measAvail = true;
  trueMeasAvail = true;
  for (int i = 0; i != ch.rows(); ++i) {
    trueMeasAvail &= !sv_->imProc->procDataQArr[ch(i)]->empty();
    measAvail &=
        duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[ch(i)]).count() >= 25;
  }

  if (trueMeasAvail)
    openLoop = false;
  else if (stateTransitionCondition)
    openLoop = true;

  if (openLoop)
    return measAvail;
  else
    return trueMeasAvail;

}

// update instantaneous trajectory vectors
void State1::updateMeasurement() {
  for (int i = 0; i != ch.rows(); ++i) {
    // update measurement vectors dy, y
    if (trueMeasAvail)
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

void State1::handleEvent(Event *event) {
  if (event->srcState != 0) {
    info("Invalid event! source state should be 0, but is actually {}", event->srcState);
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
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

  if ((event->destState == 2) && (yref(0) > 0.9 * yrefScale(0)))


  if (destReached) {
    yref = yDest;
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
    if (event->destState == 1)
      sv_->updateState<State1>();
  }
}

Eigen::Matrix<int16_t, 3, 1> State1::step() {
  // update state x and control signal u based on new measurements
  // Kalman observer
  // prediction
  // dxhat = dxhat_prev - dxref, Cd*dxref = dyref
  dxhat(ch.array(), Eigen::all) = Ad * dxhat(ch.array(), Eigen::all) +
                                  Bd * du(ch.array(), Eigen::all) -
                                  CdInv * dyref(ch.array(), Eigen::all);
  P = Ad * P * Ad_ + Qw;
  // correction
  temp = Cd * P * Cd_ + Rv;
  tempInv = temp.inverse();
  Ko = P * Cd_ * tempInv;
  dyhat(ch.array(), Eigen::all) = Cd * dxhat(ch.array(), Eigen::all);
  dxhat(ch.array(), Eigen::all) = dxhat(ch.array(), Eigen::all) +
                                  Ko * (dy(ch.array(), Eigen::all) - dyhat(ch.array(), Eigen::all));

  // apply updated control signals (with saturation limits) to pump
  // LQI control law
  du(ch.array(), Eigen::all) = -K1 * dxhat(ch.array(), Eigen::all) - K2 * z(ch.array(), Eigen::all);
  u = uref + du;

  // apply saturation (+/- 20)
  for (int i = 0; i != du.rows(); ++i) {
    if (du(i) < -20) du(i) = -20;
    else if (du(i) > 20) du(i) = 20;
  }
  usat = uref + du;

  return usat.cast<int16_t>();

}
