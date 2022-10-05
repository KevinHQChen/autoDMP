#include "ctrl/state/state1.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state0.hpp"
#include "ctrl/supervisor.hpp"

State1::State1(Supervisor *sv)
    : State(sv, Eigen::Vector3d(65, 45, 65),
            Eigen::Vector3d(0, sv->imProc->impConf.getChanBBox()[1].height,
                            sv->imProc->impConf.getChanBBox()[2].height)),
      ch(Vector2ui(1, 2)),
      // system matrices
      Ad(openData(sv->getConfPath() + "state1/Ad.txt")),
      Ad_(openData(sv->getConfPath() + "state1/Ad_.txt")),
      Bd(openData(sv->getConfPath() + "state1/Bd.txt")),
      Cd(openData(sv->getConfPath() + "state1/Cd.txt")),
      Cd_(openData(sv->getConfPath() + "state1/Cd_.txt")), CdInv(Cd.inverse()),
      K1(openData(sv->getConfPath() + "state1/K1.txt")),
      K2(openData(sv->getConfPath() + "state1/K2.txt")),
      Qw(openData(sv->getConfPath() + "state1/Qw.txt")),
      Rv(openData(sv->getConfPath() + "state1/Rv.txt")), P0(Eigen::Matrix2d::Identity(2, 2)),
      P(P0) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
}

State1::~State1() {
  // clean up any resources used by current state here
}

// check for new measurements on selected channels
bool State1::measurementAvailable() {
  bool tmpMeasAvail = true;
  for (int i = 0; i != ch.rows(); ++i) {
    trueMeasAvail[ch(i)] = !sv_->imProc->procDataQArr[ch(i)]->empty();
    measAvail[ch(i)] =
        duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[ch(i)]).count() >= 25;

    if (obsv[ch(i)] || (!obsv[ch(i)] && trueMeasAvail[ch(i)])) {
      tmpMeasAvail &= trueMeasAvail[ch(i)];
      if (!stateTransitionCondition)
        obsv[ch(i)] = true;
    }
    else
      tmpMeasAvail &= measAvail[ch(i)];
  }

  return tmpMeasAvail;
}

// update instantaneous trajectory vectors
void State1::updateMeasurement() {
  for (int i = 0; i != ch.rows(); ++i) {
    // update measurement vectors dy, y
    if (trueMeasAvail[ch(i)])
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
  if (event->srcState != 1) {
    info("Invalid event! source state should be 1, but is actually {}", event->srcState);
    delete sv_->currEvent_;
    sv_->currEvent_ = nullptr;
    return;
  }

  Eigen::Vector3d yDest = (event->destPos.array() * yrefScale.array()).matrix();
  bool destReached = true;

  // generate next waypoint if destination is not reached
  for (int i = 0; i != ch.rows(); ++i) {
    if (!firstMeasAvail[ch(i)]) {
      if (y(ch(i)) < yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) + event->vel(ch(i)) * 25e-3;
      else if (y(ch(i)) > yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) - event->vel(ch(i)) * 25e-3;
    } else {
      if (y(ch(i)) < yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) + event->vel(ch(i)) * dt[ch(i)].count();
      else if (y(ch(i)) > yDest(ch(i)))
        yref(ch(i)) = y(ch(i)) - event->vel(ch(i)) * dt[ch(i)].count();
    }
    destReached &= std::abs(yref(ch(i)) - yDest(ch(i))) < event->vel(ch(i)) * 25e-3;
  }

  // remain in State 1
  //   |__|        |__|
  //   |  |   =>   |  |
  //  /_/\ \      /_/\ \
  // / /  \ \    / /  \ \.
  if (event->destState == 1) {
    if (destReached) {
      yref = yDest;
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
    }
  }

  // transition to State 0
  //   |__|        |  |
  //   |  |   =>   |  |
  //  /_/\ \      / /\_\
  // / /  \ \    / /  \ \.
  if (event->destState == 0) {
    if (yref(1) > 0.9 * yrefScale(1) && yref(2) > 0.9 * yrefScale(2) && obsv[1] && obsv[2]) {
      obsv[1] = false;
      obsv[2] = false;
      sv_->imProc->clearProcDataQueues();
      stateTransitionCondition = true;
    }
    if (!obsv[1] && !obsv[2] && !sv_->imProc->procDataQArr[0]->empty()) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State0>();
    }
  }

  // transition to State 2
  //   |__|        |__|
  //   |  |   =>   |  |
  //  /_/\ \      / /\_\
  // / /  \ \    / /  \ \.
  if (event->destState == 2) {
    if (yref(1) > 0.9 * yrefScale(1) && obsv[1]) {
      obsv[1] = false;
      sv_->imProc->clearProcDataQueues();
      stateTransitionCondition = true;
    }
    if (!obsv[1] && !sv_->imProc->procDataQArr[0]->empty()) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      // sv_->updateState<State2>();
    }
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
    if (du(i) < -20)
      du(i) = -20;
    else if (du(i) > 20)
      du(i) = 20;
  }
  usat = uref + du;

  return usat.cast<int16_t>();
}
