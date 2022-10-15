#include "ctrl/state/sysidstate.hpp"
#include "ctrl/supervisor.hpp"

SysIDState::SysIDState(Supervisor *sv, Eigen::Vector3d uref)
    : State(sv, uref, Eigen::Vector3d(0, 0, 0)) {
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
    if (std::strcmp(sv_->prbsRows[i], "") != 0) {
      simMeasAvail_ &=
          duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime[i]).count() >= 25;
      trueMeasAvail_ &= !sv_->imProc->procDataQArr[i]->empty();
    }
  }
  if (toml::get<bool>(sv_->conf["ctrl"]["simMode"]))
    return simMeasAvail_;
  return trueMeasAvail_;
}

// update instantaneous trajectory vectors
void SysIDState::updateMeasurement() {
  // update du each time measurement is available
  if (stp < sv_->scaledPrbs.rows()) {
    du = sv_->scaledPrbs.row(stp);
    stp++;
  } else
    du = Eigen::Vector3d::Zero();
  // update y
  for (int i = 0; i < 3; ++i) {
    if (std::strcmp(sv_->prbsRows[i], "") != 0) {
      if (trueMeasAvail_)
        y(i) = sv_->imProc->procDataQArr[i]->get().y;
      else
        y(i) = 100 + (std::rand() % (300 - 100 + 1));
      prevCtrlTime[i] = steady_clock::now();
    }
  }
}

void SysIDState::handleEvent(Event *event) { return; }

Eigen::Matrix<int16_t, 3, 1> SysIDState::step() {
  u = uref + du;
  sv_->ctrlDataQueuePtr->out << u(0) << ", " << u(1) << ", " << u(2) << ", ";
  sv_->ctrlDataQueuePtr->out << y(0) << ", " << y(1) << ", " << y(2) << "\n";
  return u.cast<int16_t>();
}
