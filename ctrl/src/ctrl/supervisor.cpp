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
      evQueue_(new QueueFPS<event_bus>(dataPath + "eventQueue.txt")), currState_(new State0(this)),
      currEvent_(new Event(0, 0, Eigen::Vector3d(0.5, 0, 0), Eigen::Vector3d(10, 0, 0))),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")) {}

Supervisor::~Supervisor() {
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

    imProc->clearProcFrameQueues();
    imProc->clearTempFrameQueues();
    imProc->clearProcDataQueues();

    // initialize SupervisoryController (y_range, y_max, y_o, u_o, yhat, enAdapt)
    // - assume y_o is always y_max, and yhat is at y_max initially (TODO is this a good
    // assumption?)
    supIn.enAdapt = false;
    supIn.y_range[0] = 0.1;
    supIn.y_range[1] = 0.9;
    for (int ch = 0; ch < Config::numChans_; ++ch) {
      // TODO add support for non-90-degree channels
      if (ch == 0)
        supIn.y_max[ch] = imProc->impConf.getChanBBox()[ch].height;
      else
        supIn.y_max[ch] = imProc->impConf.getChanBBox()[ch].width;
      supIn.y_o[ch] = supIn.y_max[ch];
      supOut.yhat[ch] = supIn.y_max[ch];
    }
    for (int pp = 0; pp < NUM_PUMPS; ++pp) {
      if (pp < 2)
        supIn.u_o[0] = pump->pumpVoltages[pp];
      else
        supIn.u_o[pp - 1] = pump->pumpVoltages[pp];
    }
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
  allMeasAvail = true;
  anyMeasAvail = false;
  simMeasAvail = duration_cast<milliseconds>(steady_clock::now() - prevCtrlTime).count() >= 25;
  for (int ch = 0; ch < Config::numChans_; ++ch) {
    if (!simModeActive)
      trueMeasAvail[ch] = !imProc->procDataQArr[ch]->empty();
    else {
      if (supOut.currEv.chs[ch] == supOut.currEv.nextChs[ch])
        trueMeasAvail[ch] = simMeasAvail && supOut.currEv.chs[ch];
      else
        trueMeasAvail[ch] = !supOut.inTransRegion ? (simMeasAvail && supOut.currEv.chs[ch])
                                                  : (simMeasAvail && supOut.currEv.nextChs[ch]);
    }

    if (supOut.currEv.chs[ch]) {
      anyMeasAvail |= trueMeasAvail[ch];
      allMeasAvail &= trueMeasAvail[ch];
    }
  }

  if (allMeasAvail) {
    supIn.inputevents[0] = !supIn.inputevents[0];
    return true;
  } else if (anyMeasAvail || (simMeasAvail && supOut.inTransRegion)) {
    supIn.inputevents[1] = !supIn.inputevents[1];
    return true;
  } else
    return false;
}

bool Supervisor::updateInputs() {
  for (int ch = 0; ch < Config::numChans_; ++ch) {
    if (trueMeasAvail[ch]) {
      if (!simModeActive)
        supIn.y[ch] = imProc->procDataQArr[ch]->get().loc.y;
      else {
        if (supOut.currEv.chs[ch] == supOut.currEv.nextChs[ch])
          supIn.y[ch] = supOut.yhat[ch];
        else
          supIn.y[ch] = !supOut.inTransRegion ? supOut.yhat[ch] : supIn.y_o[ch];
      }

    } else // simMeasAvail (set y to 0 to tell the observer to only predict for this channel)
      supIn.y[ch] = 0;
  }

  ctrlDataQueuePtr->out << "allMA: " << allMeasAvail << ", anyMA: " << anyMeasAvail
                        << ", simMA: " << simMeasAvail << ", inTR: " << (bool)supOut.inTransRegion
                        << ", y: " << (double)supIn.y[0] << ", " << (double)supIn.y[1] << ", "
                        << (double)supIn.y[2] << ", yhat: " << (double)supOut.yhat[0] << ", "
                        << (double)supOut.yhat[1] << ", " << (double)supOut.yhat[2]
                        << ", u: " << (double)supOut.u[0] << ", " << (double)supOut.u[1] << ", "
                        << (double)supOut.u[2] << ", chs: " << (bool)supOut.currEv.chs[0]
                        << (bool)supOut.currEv.chs[1] << (bool)supOut.currEv.chs[2]
                        << ", nextChs: " << (bool)supOut.currEv.nextChs[0]
                        << (bool)supOut.currEv.nextChs[1] << (bool)supOut.currEv.nextChs[2] << "\n";

  if (supOut.requestEvent && !evQueue_->empty())
    supIn.nextEv = evQueue_->get();
  else
    supIn.nextEv = nullEv;

  return true;
}

void Supervisor::start() {
  while (startedCtrl) {
    if (measAvail()) {
      prevCtrlTime = steady_clock::now();
      // info("Time: {}", duration_cast<milliseconds>(prevCtrlTime - initTime).count());
      updateInputs();
      sup->rtU = supIn;
      // info("y0: {}", sup->rtU.y[0]);
      // info("y1: {}", sup->rtU.y[1]);
      // info("y2: {}", sup->rtU.y[2]);
      // info("y_max0: {}", sup->rtU.y_max[0]);
      // info("y_max1: {}", sup->rtU.y_max[1]);
      // info("y_max2: {}", sup->rtU.y_max[2]);
      // info("y_o0: {}", sup->rtU.y_o[0]);
      // info("y_o1: {}", sup->rtU.y_o[1]);
      // info("y_o2: {}", sup->rtU.y_o[2]);
      // info("u_o0: {}", sup->rtU.u_o[0]);
      // info("u_o1: {}", sup->rtU.u_o[1]);
      // info("u_o2: {}", sup->rtU.u_o[2]);
      // info("y_range0: {}", sup->rtU.y_range[0]);
      // info("y_range1: {}", sup->rtU.y_range[1]);
      // info("enAdapt: {}", sup->rtU.enAdapt);
      // info("u0: {}", sup->rtY.u[0]);
      // info("u1: {}", sup->rtY.u[1]);
      // info("u2: {}", sup->rtY.u[2]);
      // info("yhat0: {}", sup->rtY.yhat[0]);
      // info("yhat1: {}", sup->rtY.yhat[1]);
      // info("yhat2: {}", sup->rtY.yhat[2]);
      // info("inTransRegion: {}", sup->rtY.inTransRegion);
      // info("requestEvent: {}", sup->rtY.requestEvent);
      // info("currEv srcState: {}", sup->rtY.currEv.srcState);
      // info("currEv destState: {}", sup->rtY.currEv.destState);
      // info("currEv destPos0: {}", sup->rtY.currEv.destPos[0]);
      // info("currEv destPos1: {}", sup->rtY.currEv.destPos[1]);
      // info("currEv destPos2: {}", sup->rtY.currEv.destPos[2]);
      // info("currEv moveTime: {}", sup->rtY.currEv.moveTime);
      // info("currEv holdTime: {}", sup->rtY.currEv.holdTime);
      // info("currEv chs0: {}", sup->rtY.currEv.chs[0]);
      // info("currEv chs1: {}", sup->rtY.currEv.chs[1]);
      // info("currEv chs2: {}", sup->rtY.currEv.chs[2]);
      // info("currEv nextChs0: {}", sup->rtY.currEv.nextChs[0]);
      // info("currEv nextChs1: {}", sup->rtY.currEv.nextChs[1]);
      // info("currEv nextChs2: {}", sup->rtY.currEv.nextChs[2]);

      sup->step();
      supOut = sup->rtY;

      !simModeActive ? pump->sendSigs(Eigen::Matrix<int16_t, 3, 1>(supOut.u[0], supOut.u[1], 0))
                     : info("Pump inputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);
      // TODO save data to ctrlDataQueue
    }
  }
}

void Supervisor::addEvent(event_bus e) { evQueue_->push_back(e); }

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

    !simModeActive ? pump->sendSigs(currState_->uref.cast<int16_t>()) : info(currState_->uref);

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
      !simModeActive ? pump->sendSigs(currState_->step()) : info(currState_->step());
    }
  }
}
