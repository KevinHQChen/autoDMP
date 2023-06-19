#include "ctrl/supervisor.hpp"

Supervisor::Supervisor(ImProc *imProc, Pump *pump)
    : conf(Config::conf), simModeActive(toml::get<bool>(conf["ctrl"]["simMode"])),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), pump(pump), imProc(imProc),
      sup(new SupervisoryController()), supIn({}), supOut({}),
      evQueue_(new QueueFPS<event_bus>(dataPath + "eventQueue.txt")), currEv_(nullEv),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")) {}

Supervisor::~Supervisor() {
  delete evQueue_;
  delete sup;
  delete ctrlDataQueuePtr;
}

void Supervisor::startThread() {
  if (!startedCtrl) {
    info("Starting Supervisor...");
    startedCtrl = true;

    imProc->clearProcData();

    // initialize SupervisoryController (y_range, y_max, y_o, u_o, yhat)
    supIn.excitation = 5;
    for (int ch = 0; ch < imProc->impConf.numChs_; ++ch) {
      // primary
      supIn.ymax[ch] = imProc->yMax[ch];
      supIn.y0[ch] = 0;
      supOut.yhat[ch] = 0;
      // secondary
      supIn.ymax[imProc->impConf.numChs_ + ch] = imProc->yMax[ch];
      supIn.y0[imProc->impConf.numChs_ + ch] = 0;
      supOut.yhat[imProc->impConf.numChs_ + ch] = 0;

      supIn.u0[ch] = pump->pumpVoltages[ch];
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
        setCurrEv(supIn.nextEv);
        imProc->currEv = getCurrEv();
      } else
        supIn.nextEv = nullEv;

      prevCtrlTime = steady_clock::now();
      y = imProc->procData->get();
      for (int ch = 0; ch < 2 * imProc->impConf.numChs_; ++ch)
        supIn.y[ch] = y[ch];
      supIn.measAvail = !supIn.measAvail;

      // info("Time: {}", duration_cast<milliseconds>(prevCtrlTime - initTime).count());
      sup->rtU = supIn;
      sup->step();
      supOut = sup->rtY;

      !simModeActive ? pump->sendSigs(Eigen::Vector3d(supOut.u[0], supOut.u[1], supOut.u[2]))
                     : info("Pump inputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);
      // info("Pump inputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);

      ctrlDataQueuePtr->out << "y: " << (double)supIn.y[0] << ", " << (double)supIn.y[1] << ", "
                            << (double)supIn.y[2] << ", " << (double)supIn.y[3] << ", "
                            << (double)supIn.y[4] << ", " << (double)supIn.y[5]
                            << ", yhat: " << (double)supOut.yhat[0] << ", "
                            << (double)supOut.yhat[1] << ", " << (double)supOut.yhat[2] << ", "
                            << (double)supOut.yhat[3] << ", " << (double)supOut.yhat[4] << ", "
                            << (double)supOut.yhat[5] << ", u: " << (double)supOut.u[0] << ", "
                            << (double)supOut.u[1] << ", " << (double)supOut.u[2] << "\n";
      ctrlDataQueuePtr->out << "params: " << supOut.theta[0] << ", " << supOut.theta[1] << ", "
                            << supOut.theta[2] << ", " << supOut.theta[3] << ", " << supOut.theta[4]
                            << ", " << supOut.theta[5] << ", " << supOut.theta[6] << ", "
                            << supOut.theta[7] << ", " << supOut.theta[8] << ", " << supOut.theta[9]
                            << ", " << supOut.theta[10] << ", " << supOut.theta[11] << ", "
                            << supOut.theta[12] << ", " << supOut.theta[13] << ", "
                            << supOut.theta[14] << ", " << supOut.theta[15] << ", "
                            << supOut.theta[16] << ", " << supOut.theta[17] << ", "
                            << supOut.theta[18] << ", " << supOut.theta[19] << ", "
                            << supOut.theta[20] << ", " << supOut.theta[21] << ", "
                            << supOut.theta[22] << ", " << supOut.theta[23] << "\n";
    }
  }
}
