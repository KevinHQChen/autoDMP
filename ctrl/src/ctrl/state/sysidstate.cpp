#include "ctrl/state/sysidstate.hpp"
#include "ctrl/supervisor.hpp"

SysIDState::SysIDState(Supervisor *sv, Eigen::Vector3d uref_)
    : State(sv, uref_,
            Eigen::Vector3d(sv->imProc->impConf.getChanBBox()[0].height,
                            sv->imProc->impConf.getRotChanBBox()[1].height,
                            sv->imProc->impConf.getRotChanBBox()[2].height)) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
  py::initialize_interpreter();
  py::eval_file("scripts/sysid.py"); // import sysid functions

  simMeas = py::module::import("sim_meas").attr("sim_meas");
}

SysIDState::SysIDState(Supervisor *sv, Eigen::Vector3d uref, bool *selChs,
                       std::vector<float> minVals, std::vector<float> maxVals,
                       unsigned int numSamples)
    : State(sv, uref,
            Eigen::Vector3d(sv->imProc->impConf.getChanBBox()[0].height,
                            sv->imProc->impConf.getRotChanBBox()[1].height,
                            sv->imProc->impConf.getRotChanBBox()[2].height)),
      selChs_(selChs), minVals_(minVals), maxVals_(maxVals), numSamples_(numSamples) {
  // clear all improc queues
  sv_->imProc->clearProcDataQueues();
  py::initialize_interpreter();
  py::eval_file("scripts/sysid.py"); // import sysid functions

  simMeas = py::module::import("sim_meas").attr("sim_meas");
  // incFcn = py::module::import("__main__").attr("incrementfcn");
}

SysIDState::~SysIDState() {
  // clean up any resources used by current state here
  py::finalize_interpreter();
}

// check for new measurements on selected channels
bool SysIDState::measurementAvailable() {
  simMeasAvail_ = true, trueMeasAvail_ = true;
  for (int i = 0; i < 3; ++i) {
    if (selChs_[i]) {
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
  if (trueMeasAvail_) { // update true measurement y
    for (int i = 0; i < 3; ++i)
      if (selChs_[i])
        y(i) = sv_->imProc->procDataQArr[i]->get().loc.y;
  } else // use simulated measurement based on a noisy second order model
    y = simMeas(u, y).cast<Eigen::Vector3d>();

  for (int i = 0; i < 3; ++i)
    if (selChs_[i])
      prevCtrlTime[i] = steady_clock::now();
}

void SysIDState::handleEvent(Event *event) { return; }

Eigen::Matrix<int16_t, 3, 1> SysIDState::step() {
  // update control signal u based on new measurements

  // apply updated control signal to pump

  std::cout << "the c++ value before calling the python script: " << stp << std::endl;

  if (!simMeas.is_none())
    stp = simMeas(stp).cast<int>();

  std::cout << "the c++ value after calling the python script: " << stp << std::endl;

  u = uref + du;
  sv_->ctrlDataQueuePtr->out << u(0) << ", " << u(1) << ", " << u(2) << ", ";
  sv_->ctrlDataQueuePtr->out << y(0) << ", " << y(1) << ", " << y(2) << "\n";
  return u.cast<int16_t>();
}
