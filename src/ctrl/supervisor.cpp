#include "ctrl/supervisor.hpp" // ensures supervisor.hpp compiles in isolation
#include "ctrl/state/state0.hpp"

Supervisor::Supervisor(ImProc *imProc)
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      currState_(new State0(this)), imProc(imProc) {}

Supervisor::~Supervisor() {
    delete currState_;
}

void Supervisor::startThread() {
  if (!started()) {
    info("Starting Supervisor...");
    startedCtrl = true;
    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
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
    this->clearCtrlDataQueue();
  }
}

void Supervisor::start() {
  while (started()) {
    if (currState_->measurementAvailable()) {
      // events are pushed to a FIFO event queue by GUI
      // get the first event in the queue and pop it
      if (currEvent_ == nullptr && !eventQueue_->empty())
          currEvent_ = eventQueue_->get();

      // call the current state's trajectory generation function corresponding to the event
      // (e.g. generate, merge, etc)
      if (currEvent_ != nullptr)
        currState_->handleEvent(currEvent_);

      currState_->updateMeasurement();

      // call the current state's step function to generate optimal control signals at current time
      // step
      controlSignals = step();
      setPump(controlSignals);
    }
  }
}

bool Supervisor::started() {
  return startedCtrl;
}

template <typename T> void Supervisor::updateState() {
    delete currState_;
    currState_ = new T(this);
}




Eigen::Matrix<int16_t, 3, 1> Supervisor::getCtrlData() {
    if (!ctrlDataQueuePtr->empty())
        return ctrlDataQueuePtr->get();
    return Eigen::Matrix<int16_t, 3, 1>::Zero();
}

bool Supervisor::procDataAvail() {
    // TODO
    bool procDataAvail = !ProcDataQueue.empty();
    if(procDataAvail)
        currPos = ProcDataQueue->get();
    else if (25ms has passed)
        currPos = ProcDataQueue->get();
    else
        procDataAvail = false;

    return procDataAvail;
}

void Supervisor::step() {
    usat = currState_->step();
}

void Supervisor::hold() {
    usat = currState_->hold();
}



// state0
class State0 : BaseState {

    public:
        step() {
            return BaseState->step(Model model, Trajectory traj);
        }
};


// state1



// state2

// main
int main() {
    Supervisor *supervisor = new Supervisor();

    while (true) {
        if(Supervisor->procDataAvail()) {
            controlSignals = Supervisor->step();
            setPump(controlSignals);
        }
    }

    delete supervisor;
}
