#include "ctrl/supervisor.hpp"
#include "ctrl/state/state0.hpp"
#include "ctrl/state/state2.hpp"
#include "ctrl/state/sysidstate.hpp"

#define NUM_PUMPS 4

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

    imProc->clearProcData();

    // initialize SupervisoryController (y_range, y_max, y_o, u_o, yhat)
    // - assume y_o is always y_max, and yhat is at y_max initially (TODO is this a good
    // assumption?)
    supIn.excitation = 5;
    for (int ch = 0; ch < imProc->impConf.numChs_; ++ch) {
      supIn.y_max[ch] = imProc->yMax[ch];
      supIn.y0[ch] = supIn.y_max[ch];
      supOut.yhat[ch] = supIn.y_max[ch];
    }
    for (int pp = 0; pp < NUM_PUMPS; ++pp) {
      if (pp < 2)
        supIn.u0[0] = pump->pumpVoltages[pp];
      else
        supIn.u0[pp - 1] = pump->pumpVoltages[pp];
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
    imProc->clearProcData();
    // this->clearCtrlDataQueue();
  }
}

void Supervisor::start() {
  while (startedCtrl) {
    if (!imProc->procData->empty()) {
      if (supOut.requestEvent && !evQueue_->empty()) {
        supIn.nextEv = evQueue_->get();
        currEv_ = supIn.nextEv;
      } else
        supIn.nextEv = nullEv;

      prevCtrlTime = steady_clock::now();
      p = imProc->procData->get();
      updateMeas();

      // info("Time: {}", duration_cast<milliseconds>(prevCtrlTime - initTime).count());
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

      !simModeActive ? pump->sendSigs(Eigen::Vector3d(supOut.u[0], supOut.u[1], supOut.u[2]))
                     : info("Pump inputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);

      ctrlDataQueuePtr->out << "y: " << (double)supIn.ymeas[0] << ", " << (double)supIn.ymeas[1]
                            << ", " << (double)supIn.ymeas[2]
                            << ", yhat: " << (double)supOut.yhat[0] << ", "
                            << (double)supOut.yhat[1] << ", " << (double)supOut.yhat[2]
                            << ", u: " << (double)supOut.u[0] << ", " << (double)supOut.u[1] << ", "
                            << (double)supOut.u[2] << ", chs: " << (bool)supOut.currEv.chs[0]
                            << (bool)supOut.currEv.chs[1] << (bool)supOut.currEv.chs[2]
                            << ", nextChs: " << (bool)supOut.currEv.nextChs[0]
                            << (bool)supOut.currEv.nextChs[1] << (bool)supOut.currEv.nextChs[2]
                            << "\n";
      ctrlDataQueuePtr->out << "params: " << supOut.B_b[0] << ", " << supOut.B_b[1] << ", "
                            << supOut.B_b[2] << ", " << supOut.B_b[3] << ", " << supOut.B_b[4]
                            << ", " << supOut.B_b[5] << ", " << supOut.B_b[6] << ", "
                            << supOut.B_b[7] << ", " << supOut.B_b[8] << "\n";
    }
  }
}

void Supervisor::updateMeas() {
  if (currEv_.srcState == 0) {
    if (p[0].found[0]) {
      supIn.ymeas[0] = p[0].p[0];
      supIn.ymeas[1] = supIn.ymeas[2] = 0;
    }
    if (p[1].found[1] && p[2].found[1]) { // state0 -> state1
      supIn.ymeas[0] = 0;
      supIn.ymeas[1] = p[1].p[1];
      supIn.ymeas[2] = p[2].p[1];
    }
  }

  if (currEv_.srcState == 1) {
    if (p[1].found[0] && p[2].found[0]) {
      supIn.ymeas[0] = 0;
      supIn.ymeas[1] = p[1].p[0];
      supIn.ymeas[2] = p[2].p[0];
    }
    if (p[1].found[1] && p[2].found[1]) {
      supIn.ymeas[0] = 0;
      supIn.ymeas[1] = p[1].p[1];
      supIn.ymeas[2] = p[2].p[1];
    }
    if (p[1].found[2] && p[2].found[2]) {
      supIn.ymeas[0] = 0;
      supIn.ymeas[1] = (p[1].p[2] > imProc->yMax[1]) ? p[1].p[1] : p[1].p[2];
      supIn.ymeas[2] = (p[2].p[2] > imProc->yMax[2]) ? p[2].p[1] : p[2].p[2];
    }
    if ((p[0].found[0] && p[0].p[0] < imProc->yMax[0]) &&
        (p[2].found[2] && p[2].p[2] < imProc->yMax[2])) { // state1 -> state2
      supIn.ymeas[0] = p[0].p[0];
      supIn.ymeas[1] = 0;
      supIn.ymeas[2] = p[2].p[2];
    }
  }

  if (currEv_.srcState == 2) {
    if (p[0].found[0] && p[2].found[2]) {
      supIn.ymeas[0] = p[0].p[0];
      supIn.ymeas[1] = 0;
      supIn.ymeas[2] = p[2].p[2];
    }
    if (p[1].found[1] && p[2].found[1]) { // state2 -> state1
      supIn.ymeas[0] = 0;
      supIn.ymeas[1] = p[1].p[1];
      supIn.ymeas[2] = p[2].p[1];
    }
  }
  supIn.measAvail = !supIn.measAvail;
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

    !simModeActive ? pump->sendSigs(currState_->uref) : info(currState_->uref);

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
