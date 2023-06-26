#include "ctrl/supervisor.hpp"

Supervisor::Supervisor(ImProc *imProc, Pump *pump, std::shared_ptr<logger> log)
    : conf(Config::conf), simModeActive(toml::get<bool>(conf["ctrl"]["simMode"])),
      dataPath(toml::get<std::string>(conf["ctrl"]["dataPath"])),
      confPath(toml::get<std::string>(conf["ctrl"]["confPath"])), pump(pump), imProc(imProc),
      sup(new SupervisoryController()), supIn({}), supOut({}),
      evQueue_(new QueueFPS<event_bus>(dataPath + "eventQueue.txt")),
      ctrlDataQueuePtr(new QueueFPS<int>(dataPath + "ctrlDataQueue.txt")), lg(log) {
  lg->info("Initializing Supervisor...");
  sup->initialize();
  rtM = sup->getRTM();
  nullEv = *(event_bus *)getSupervisorParam(&rtM->DataMapInfo.mmi, 0);
  currEv_ = nullEv;
}

Supervisor::~Supervisor() {
  lg->info("Terminating Supervisor...");
  stopThread();
  delete evQueue_;
  delete sup;
  delete ctrlDataQueuePtr;
}

void Supervisor::startThread() {
  if (!startedCtrl) {
    lg->info("Starting Supervisor...");
    startedCtrl = true;

    imProc->clearData();
    no = imProc->impConf.getNumChs();
    lg->info("Supervisor found {} channels.", no);

    // initialize SupervisoryController (y_range, y_max, y_o, u_o, yhat)
    supIn.excitation = 5;
    for (int ch = 0; ch < no; ++ch) {
      // primary
      supIn.ymax[ch] = imProc->yMax[ch];
      supIn.y0[ch] = 0;
      supOut.yhat[ch] = 0;
      // secondary
      supIn.ymax[no + ch] = imProc->yMax[ch];
      supIn.y0[no + ch] = 0;
      supOut.yhat[no + ch] = 0;

      supIn.u0[ch] = pump->outputs[ch];
    }
    sup->initialize();

    if (!simModeActive && pump->getPumpType() == "BARTELS")
      pump->setFreq(200);
    ctrlThread = std::thread(&Supervisor::start, this);
    ctrlThread.detach();
  }
}

void Supervisor::stopThread() {
  if (startedCtrl) {
    lg->info("Stopping Supervisor...");
    startedCtrl = false;
    std::lock_guard<std::mutex> guard(ctrlMtx); // wait for thread to finish
    supIn = {};
    supOut = {};
    sup->terminate();
    imProc->clearData();
    // this->clearCtrlDataQueue();
  }
}

void Supervisor::start() {
  while (startedCtrl) {
    std::lock_guard<std::mutex> guard(ctrlMtx);
    if (!imProc->procData->empty()) {
      // Timer t("ImProc");
      if (supOut.requestEvent && !evQueue_->empty()) {
        supIn.nextEv = evQueue_->get();
        setCurrEv(supIn.nextEv);
        imProc->setR(getCurrEv().r);
      } else
        supIn.nextEv = nullEv;

      prevCtrlTime = steady_clock::now();
      y = imProc->procData->get();
      for (int ch = 0; ch < 2 * no; ++ch)
        supIn.y[ch] = y[ch];
      supIn.measAvail = !supIn.measAvail;

      // lg->info("Time: {}", duration_cast<milliseconds>(prevCtrlTime - initTime).count());
      sup->rtU = supIn;
      sup->step();
      supOut = sup->rtY;

      !simModeActive ? pump->setOutputs(std::vector<double>(supOut.u, supOut.u + no))
                     : lg->info("Pump outputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);
      // lg->info("Pump inputs: {}, {}, {}", supOut.u[0], supOut.u[1], supOut.u[2]);

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

bool Supervisor::started() { return startedCtrl; }

void *Supervisor::getSupervisorParam(rtwCAPI_ModelMappingInfo *mmi, uint_T paramIdx) {
  const rtwCAPI_ModelParameters *modelParams;
  void **dataAddrMap;

  uint_T addrIdx;

  void *paramAddress;

  /* Assert the parameter index is less than total number of parameters */
  lg->info("SupervisoryController variable block parameters: {}",
           rtwCAPI_GetNumModelParameters(mmi));
  assert(paramIdx < rtwCAPI_GetNumModelParameters(mmi));

  /* Get modelParams, an array of rtwCAPI_ModelParameters structure  */
  modelParams = rtwCAPI_GetModelParameters(mmi);
  if (modelParams == NULL) {
    lg->error("Model parameters not available");
    return nullptr;
  }

  /* Get the address to this parameter */
  dataAddrMap = rtwCAPI_GetDataAddressMap(mmi);
  addrIdx = rtwCAPI_GetModelParameterAddrIdx(modelParams, paramIdx);
  paramAddress = (void *)rtwCAPI_GetDataAddress(dataAddrMap, addrIdx);
  if (paramAddress == NULL) {
    lg->error("Model parameter address not available");
    return nullptr;
  }

  return paramAddress;
}
