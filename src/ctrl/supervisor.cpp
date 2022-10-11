#include "ctrl/supervisor.hpp" // ensures supervisor.hpp compiles in isolation
#include "ctrl/state/state0.hpp"

Supervisor::Supervisor(ImProc *imProc)
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), imProc(imProc),
      currState_(new State0(this)),
      currEvent_(new Event(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0))),
      eventQueue_(new QueueFPS<Event *>(dataPath + "eventQueue.txt")) {}

Supervisor::~Supervisor() {
  delete eventQueue_;
  delete currEvent_;
  delete currState_;
}

void Supervisor::startThread() {
  if (!started()) {
    info("Starting Supervisor...");
    startedCtrl = true;
    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();
    ctrlThread = std::thread(&Supervisor::start, this);
    ctrlThread.detach();
  }
}

void Supervisor::stopThread() {
  if (started()) {
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
  while (started()) {
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
      info(currState_->step());
      // setPump(currState_->step());
    }
  }
}

bool Supervisor::started() { return startedCtrl; }

void Supervisor::addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel) {
  eventQueue_->push(new Event(srcState, destState, pos, vel));
}

// Eigen::Matrix<int16_t, 3, 1> Supervisor::getCtrlData() {
//     if (!ctrlDataQueuePtr->empty())
//         return ctrlDataQueuePtr->get();
//     return Eigen::Matrix<int16_t, 3, 1>::Zero();
// }