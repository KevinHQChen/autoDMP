#include "ctrl/supervisor.hpp" // ensures supervisor.hpp compiles in isolation
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state2.hpp"
#include "ctrl/state/sysidstate.hpp"

Supervisor::Supervisor(ImProc* imProc, Pump* pump)
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")),
      simModeActive(toml::get<bool>(conf["ctrl"]["simMode"])),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), pump(pump), imProc(imProc),
      currState_(new State0(this, Eigen::Vector3d(100, 70, 70))),
      currEvent_(new Event(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0))),
      eventQueue_(new QueueFPS<Event *>(dataPath + "eventQueue.txt")),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")) {}

Supervisor::~Supervisor() {
  delete eventQueue_;
  delete currEvent_;
  delete currState_;
}

void Supervisor::startThread() {
  if (!startedCtrl) {
    info("Starting Supervisor...");
    startedCtrl = true;
    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();
    updateState<State0>(Eigen::Vector3d(94, 97, 97));
    // updateState<State1>(Eigen::Vector3d(90, 60, 50));
    // updateState<State2>(Eigen::Vector3d(135, 101, 101));
#ifdef USEPIEZOPUMP
    if (!simModeActive)
      pump->setFreq(200);
#endif
    ctrlThread = std::thread(&Supervisor::start, this);
    ctrlThread.detach();
  }
}

void Supervisor::stopThread() {
  if (startedCtrl) {
    info("Stopping Supervisor...");
    startedCtrl = false;
    if (ctrlThread.joinable())
      ctrlThread.join();
    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();
    // this->clearCtrlDataQueue();
  }
}

void Supervisor::start() {
  while (startedCtrl) {
    if (currState_->measurementAvailable()) {
      currState_->updateMeasurement();

      // events are pushed to a FIFO event queue by GUI
      // get the first event in the queue and pop it
      if (currEvent_ == nullptr && !eventQueue_->empty())
        currEvent_ = eventQueue_->get();

      // call the current state's trajectory generation function corresponding to the event
      if (currEvent_ != nullptr)
        currState_->handleEvent(currEvent_);

      // generate optimal control signals at current time step
      // and save to ctrlDataQueue
      if (!simModeActive)
        pump->sendSigs(currState_->step());
      else
        info(currState_->step());
    }
  }
}

void Supervisor::startSysIDThread(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                                  std::vector<float> maxVals, Eigen::MatrixXd &data) {
  if (!startedSysIDFlag) {
    info("Starting SysID...");
    startedSysIDFlag = true;

#ifdef USEPIEZOPUMP
    if (!simModeActive)
      pump->setFreq(200);
#endif

    updateState<SysIDState>(uref, selChs, minVals, maxVals, data);
    sysIDThread = std::thread(&Supervisor::startSysID, this);
    sysIDThread.detach();
  }
}

void Supervisor::stopSysIDThread() {
  if (startedSysIDFlag) {
    info("Stopping SysID...");
    startedSysIDFlag = false;

    if (!simModeActive)
      pump->sendSigs(currState_->uref.cast<int16_t>());
    else
      info(currState_->uref);

    if (sysIDThread.joinable())
      sysIDThread.join();
    updateState<State0>(Eigen::Vector3d(110, 79, 126));
  }
}

void Supervisor::startSysID() {
  while (startedSysIDFlag) {
    if (currState_->measurementAvailable()) {

      currState_->updateMeasurement();
      // send excitation signal to pump and save to ctrlDataQueue
      if (!simModeActive)
        pump->sendSigs(currState_->step());
      else
        info(currState_->step());
    }
  }
}

void Supervisor::addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel) {
  eventQueue_->push(new Event(srcState, destState, pos, vel));
}
