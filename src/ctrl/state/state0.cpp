#include "ctrl/state/state0.hpp"
#include "ctrl/state/state.hpp"
#include "ctrl/state/state1.hpp"
#include "ctrl/supervisor.hpp"

State0::State0(Supervisor *sv)
    : State(sv, Eigen::Vector3d(60, 40, 60),
            Eigen::Vector3d(sv->imProc->impConf.getRotChanBBox()[0].height, 0, 0))
{
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();

  ch<1>(Vector1ui(0));
  Ad<1>(openData(sv->getConfPath() + "state0/Ad.txt"));
  Ad_<1>(openData(sv->getConfPath() + "state0/Ad_.txt"));
  Bd<1>(openData(sv->getConfPath() + "state0/Bd.txt"));
  Cd<1>(openData(sv->getConfPath() + "state0/Cd.txt"));
  Cd_<1>(openData(sv->getConfPath() + "state0/Cd_.txt"));
  CdInv<1>(Cd<1>.inverse());
  K1<1>(openData(sv->getConfPath() + "state0/K1.txt"));
  K2<1>(openData(sv->getConfPath() + "state0/K2.txt"));
  Qw<1>(openData(sv->getConfPath() + "state0/Qw.txt"));
  Rv<1>(openData(sv->getConfPath() + "state0/Rv.txt"));
  P0<1>(Vector1d::Identity(1, 1));
  P<1>(P0<1>);
}

State0::~State0() {
  // clean up any resources used by current state here
}

bool State0::measurementAvailable() { return State::measurementAvailable<1>(ch); }



// (only called when new measurements are available)
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

  // generate next waypoint if destination is not reached
  for (int i = 0; i != ch<1>.rows(); ++i) {
    if (y(ch<1>(i)) < yDest(ch<1>(i)))
      dyref(ch<1>(i)) = event->vel(ch<1>(i)) * dt[ch<1>(i)].count();
    else if (y(ch<1>(i)) > yDest(ch<1>(i)))
      dyref(ch<1>(i)) = -event->vel(ch<1>(i)) * dt[ch<1>(i)].count();
    if (stateTransitionCondition)
      dyref(ch<1>(i)) = event->vel(ch<1>(i)) * dt[ch<1>(i)].count();
    yref(ch<1>(i)) += dyref(ch<1>(i));

    if (std::abs(yref(ch<1>(i)) - yDest(ch<1>(i))) < event->vel(ch<1>(i)) * 50e-3) {
      dyref(ch<1>(i)) = 0;
      yref(ch<1>(i)) = yDest(ch<1>(i));
    }

    destReached &= std::abs(y(ch<1>(i)) - yDest(ch<1>(i))) < 1; // event->vel(ch(i)) * 25e-3;
  }

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
  //   |  |        |__|
  //   |  |   =>   |  |
  //  / /\_\      /_/\ \
  // / /  \ \    / /  \ \.
  if (event->destState == 1) {
    // ch0 is 85% to junction and we're observing ch0
    if (yref(0) > 0.85 * yrefScale(0) && obsv[0]) {
      obsv[0] = false;
      sv_->imProc->clearProcDataQueues();
      stateTransitionCondition = true;
    }

    if (stateTransitionCondition) {
    }

    // we're in sim mode and ch0 is 95% to junction
    // or ch1 or ch2 are observable and we're not observing ch0
    if ((yref(0) > 0.95 * yrefScale(0) && sv_->simModeActive) ||
        (!obsv[0] &&
         (!sv_->imProc->procDataQArr[1]->empty() || !sv_->imProc->procDataQArr[2]->empty()))) {
      delete sv_->currEvent_;
      sv_->currEvent_ = nullptr;
      sv_->updateState<State1>(usat);
    }
  }
}
