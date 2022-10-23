#include "ctrl/state/state2.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state1.hpp"
#include "ctrl/supervisor.hpp"

State2::State2(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_,
            Eigen::Vector3d(sv->imProc->impConf.getRotChanBBox()[0].height, 0,
                            sv->imProc->impConf.getChanBBox()[2].height)),
      ch(Vector2ui(0, 2)),
      // system matrices
      Ad(openData(sv->getConfPath() + "state2/Ad.txt")),
      Ad_(openData(sv->getConfPath() + "state2/Ad_.txt")),
      Bd(openData(sv->getConfPath() + "state2/Bd.txt")),
      Cd(openData(sv->getConfPath() + "state2/Cd.txt")),
      Cd_(openData(sv->getConfPath() + "state2/Cd_.txt")), CdInv(Cd.inverse()),
      K1(openData(sv->getConfPath() + "state2/K1.txt")),
      K2(openData(sv->getConfPath() + "state2/K2.txt")),
      Qw(openData(sv->getConfPath() + "state2/Qw.txt")),
      Rv(openData(sv->getConfPath() + "state2/Rv.txt")), P0(Eigen::Matrix2d::Identity(2, 2)),
      P(P0) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
  stateTransitionCondition = true;
}

State2::~State2() {
  // clean up any resources used by current state here
}

bool State2::measurementAvailable() { return State::measurementAvailable<2>(ch); }

void State2::updateMeasurement() { State::updateMeasurement<2>(ch); }

void State2::handleEvent(Event *event) {
  if (event->srcState != 2) {
    info("Invalid event! source state should be 2, but is actually {}", event->srcState);
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

  // remain in State 2
  //   |__|        |__|
  //   |  |   =>   |  |
  //  / /\_\      / /\_\
  // / /  \ \    / /  \ \.
  if (event->destState == 2) {
    if (destReached) {
      yref = yDest;
      startEvent = false;
      // z = Eigen::Vector3d::Zero();
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
    }
  }
}

Eigen::Matrix<int16_t, 3, 1> State2::step() {
  return State::step<2>(ch, Ad, Ad_, Bd, Cd, Cd_, CdInv, K1, K2, Qw, Rv, P, Ko, temp, tempInv);
}
