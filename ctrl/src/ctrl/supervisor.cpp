#include "ctrl/supervisor.hpp" // ensures supervisor.hpp compiles in isolation
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state2.hpp"
#include "ctrl/state/sysidstate.hpp"

#define NUM_PUMPS 4

event_bus Supervisor::nullEv = {};

Supervisor::Supervisor(ImProc *imProc, Pump *pump)
    : conf(Config::conf), simModeActive(toml::get<bool>(conf["ctrl"]["simMode"])),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), pump(pump), imProc(imProc),
      sup(new SupervisoryController()), supIn({}), supOut({}),
      evQueue_(new QueueFPS<event_bus>(dataPath + "evQueue.txt")), currState_(new State0(this)),
      currEvent_(new Event(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0))),
      eventQueue_(new QueueFPS<Event *>(dataPath + "eventQueue.txt")),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")) {}

Supervisor::~Supervisor() {
  delete eventQueue_;
  delete evQueue_;
  delete currEvent_;
  delete currState_;
  delete sup;
  delete ctrlDataQueuePtr;
}

void Supervisor::startThread() {
  if (!startedCtrl) {
    info("Starting Supervisor...");
    startedCtrl = true;
    firstMeasAvail = false;

    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();

    // initialize SupervisoryController (y_range, y_max, y_o, u_o, yhat)
    // - assume y_o is always y_max, and yhat is at y_max initially (TODO is this a good
    // assumption?)
    supIn.y_range[0] = 0.1;
    supIn.y_range[1] = 0.9;
    for (int ch = 0; ch < Config::numChans_; ++ch) {
      supIn.y_max[ch] =
          imProc->impConf.getChanBBox()[ch].height; // TODO add support for rotChanBBoxes
      supIn.y_o[ch] = supIn.y_max[ch];
      supOut.yhat[ch] = supIn.y_max[ch];
    }
    for (int pp = 0; pp < NUM_PUMPS; ++pp)
      supIn.u_o[pp] = pump->pumpVoltages[pp];
    sup->initialize();

    if (!simModeActive)
      pump->setFreq(200);
    ctrlThread = std::thread(&Supervisor::start, this);
    ctrlThread.detach();
  }
}

void Supervisor::stopThread() {
  if (startedCtrl) {
    info("Stopping Supervisor...");
    startedCtrl = false;

    supIn = {};
    supOut = {};
    sup->terminate();

    if (ctrlThread.joinable())
      ctrlThread.join();
    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();
    // this->clearCtrlDataQueue();
  }
}

bool Supervisor::measAvail() {
  supIn.inputevents[0] = false;
  supIn.inputevents[1] = false;

  // first measurement should only check ch0
  if (!firstMeasAvail)
    supOut.currEv.chs[0] = true;

  bool currChMeasAvail = true;
  allMeasAvail = true;
  anyMeasAvail = false;
  for (int ch = 0; ch < Config::numChans_; ++ch) {
    if (supOut.currEv.chs[ch]) {
      currChMeasAvail =
          !simModeActive
              ? !imProc->procDataQArr[ch]->empty()
              : duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime).count() >= 25;
      allMeasAvail &= currChMeasAvail;
      anyMeasAvail |= currChMeasAvail;

      if (currChMeasAvail)
        supIn.y[ch] = !simModeActive ? imProc->procDataQArr[ch]->get().loc.y : supOut.yhat[ch];
      else
        supIn.y[ch] = 0;
    }
  }
  simMeasAvail = duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime).count() >= 25;

  if (allMeasAvail) {
    firstMeasAvail = true;
    supIn.inputevents[0] = true;
    return true;
  } else if (anyMeasAvail || simMeasAvail) {
    supIn.inputevents[1] = true;
    return true;
  } else
    return false;
}

bool Supervisor::updateInputs() {
  if (supOut.requestEvent && !eventQueue_->empty())
    supIn.nextEv = evQueue_->get();
  else
    supIn.nextEv = nullEv;

  return true;
}

void Supervisor::start() {
  while (startedCtrl) {
    if (measAvail()) {
      updateInputs();
      sup->rtU = supIn;
      sup->step();
      supOut = sup->rtY;

      // // generate optimal control signals at current time step
      // // and save to ctrlDataQueue
      // if (!simModeActive)
      //   pump->sendSigs(currState_->step());
      // else
      //   info(currState_->step());
    }
  }
}

void Supervisor::startSysIDThread(Eigen::Vector3d uref, bool *selChs, std::vector<float> minVals,
                                  std::vector<float> maxVals, Eigen::MatrixXd &data) {
  if (!startedSysIDFlag) {
    info("Starting SysID...");
    startedSysIDFlag = true;

    if (!simModeActive)
      pump->setFreq(200);

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
    updateState<State0>();
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
