#include "ctrl/supervisor.hpp" // ensures supervisor.hpp compiles in isolation
#include "ctrl/state/state0.hpp"
#include "ctrl/state/sysidstate.hpp"

Supervisor::Supervisor(ImProc *imProc)
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), imProc(imProc),
      currState_(new State0(this)),
      currEvent_(new Event(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0))),
      eventQueue_(new QueueFPS<Event *>(dataPath + "eventQueue.txt")),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")) {}

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

      // TODO save ctrl data to file
    }
  }
}

void Supervisor::startSysIDThread() {
  if (!startedSysID()) {
    info("Starting SysID...");
    startedSysIDFlag = true;
    delete currState_;
    currState_ = new SysIDState(this, prbsUref);
    sysIDThread = std::thread(&Supervisor::startSysID, this);
    sysIDThread.detach();
  }
}

void Supervisor::stopSysIDThread() {
  if (startedSysID()) {
    info("Stopping SysID...");
    startedSysIDFlag = false;
    if (sysIDThread.joinable())
      sysIDThread.join();
  }
}

void Supervisor::startSysID() {
  while (startedSysID()) {
    // check if measurement is available
    bool measAvail = true, simMeasAvail = true, trueMeasAvail = true;
    for (int i = 0; i < 3; ++i) {
      if (std::strcmp(prbsRows[i], "") != 0) {
        simMeasAvail &=
            duration_cast<milliseconds>(steady_clock::now() - currState_->prevCtrlTime[i])
                .count() >= 25;
        trueMeasAvail &= !imProc->procDataQArr[i]->empty();
      }
    }
    if (toml::get<bool>(conf["ctrl"]["simMode"]))
      measAvail = simMeasAvail;
    else
      measAvail = trueMeasAvail;

    // update measurement and send excitation signal
    if (measAvail) {
      // update du each time measurement is available
      if (currState_->stp < scaledPrbs.rows()) {
        currState_->du = scaledPrbs.row(currState_->stp);
        currState_->stp++;
      } else
        currState_->du = Eigen::Vector3d::Zero();
      // update y
      for (int i = 0; i < 3; ++i) {
        if (std::strcmp(prbsRows[i], "") != 0) {
          if (trueMeasAvail)
            currState_->y(i) = imProc->procDataQArr[i]->get().y;
          else
            currState_->y(i) = 100 + (std::rand() % (300 - 100 + 1));
          currState_->prevCtrlTime[i] = steady_clock::now();
        }
      }
      // send excitation signal
      currState_->u = currState_->uref + currState_->du;
      // setPump(currState_->u.cast<int16_t>());

      // TODO save ctrl data to file
    }
  }
}

bool Supervisor::started() { return startedCtrl; }

bool Supervisor::startedSysID() { return startedSysIDFlag; }

void Supervisor::addEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel) {
  eventQueue_->push(new Event(srcState, destState, pos, vel));
}

// Eigen::Matrix<int16_t, 3, 1> Supervisor::getCtrlData() {
//     if (!ctrlDataQueuePtr->empty())
//         return ctrlDataQueuePtr->get();
//     return Eigen::Matrix<int16_t, 3, 1>::Zero();
// }
