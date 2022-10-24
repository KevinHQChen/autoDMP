#include "ctrl/state/sysidstate.hpp"
#include "ctrl/supervisor.hpp"

SysIDState::SysIDState(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_,
            Eigen::Vector3d(sv->imProc->impConf.getChanBBox()[0].height,
                            sv->imProc->impConf.getRotChanBBox()[1].height,
                            sv->imProc->impConf.getRotChanBBox()[2].height)) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
}

SysIDState::~SysIDState() {
  // clean up any resources used by current state here
}

// check for new measurements on selected channels
bool SysIDState::measurementAvailable() {
  simMeasAvail_ = true, trueMeasAvail_ = true;
  for (int i = 0; i < 3; ++i) {
    if (std::strcmp(sv_->sysidCh[i], "") != 0) {
      simMeasAvail_ &=
          duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[i]).count() >= 25;
      trueMeasAvail_ &= !sv_->imProc->procDataQArr[i]->empty();
    }
  }
  if (sv_->simModeActive)
    return simMeasAvail_;
  return trueMeasAvail_;
}

// update instantaneous trajectory vectors
void SysIDState::updateMeasurement() {
  for (int i = 0; i < 3; ++i) {
    if (std::strcmp(sv_->sysidCh[i], "") != 0) {
      // update y
      if (trueMeasAvail_)
        y(i) = sv_->imProc->procDataQArr[i]->get().y;
      else // generate random simulated measurement
        y(i) = 50 + (std::rand() % (600 - 50 + 1));
      prevCtrlTime[i] = steady_clock::now();

      // update du
      if (stp <= sv_->sysidSamples) {
        // state1 sysID settings
        if (stp % 10 == 0) {
          sv_->sysidMin = (40 + (std::rand() % (60 - 40 + 1))) / 100.0;
          sv_->sysidMax = sv_->sysidMin + 0.2;
        }

        // state0 sysID settings
        // if (stp % 10 == 0) {
        //   sv_->sysidMin = (20 + (std::rand() % (40 - 20 + 1))) / 100.0;
        //   sv_->sysidMax = 1 - sv_->sysidMin;
        // }

        if (y(i) > sv_->sysidMax * yrefScale(i))
          du(i) = -sv_->sysidDu[i];
        else if (y(i) < sv_->sysidMin * yrefScale(i))
          du(i) = sv_->sysidDu[i];
      } else
        du(i) = 0;
    }
  }
  stp++;
}

void SysIDState::handleEvent(Event *event) { return; }

Eigen::Matrix<int16_t, 3, 1> SysIDState::step() {
  u = uref + du;
  sv_->ctrlDataQueuePtr->out << u(0) << ", " << u(1) << ", " << u(2) << ", ";
  sv_->ctrlDataQueuePtr->out << y(0) << ", " << y(1) << ", " << y(2) << "\n";
  return u.cast<int16_t>();
}
