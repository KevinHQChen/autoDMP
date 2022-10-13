#include "ctrl/state/sysidstate.hpp"
#include "ctrl/supervisor.hpp"

SysIDState::SysIDState(Supervisor *sv, Eigen::Vector3d uref)
    : State(sv, uref,
            Eigen::Vector3d(0, 0, 0)) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
}

SysIDState::~SysIDState() {
  // clean up any resources used by current state here
}

// check for new measurements on selected channels
bool SysIDState::measurementAvailable() { return false; }

// update instantaneous trajectory vectors
void SysIDState::updateMeasurement() { return; }

void SysIDState::handleEvent(Event *event) { return; }

Eigen::Matrix<int16_t, 3, 1> SysIDState::step() { return Eigen::Matrix<int16_t, 3, 1>(0, 0, 0); }
