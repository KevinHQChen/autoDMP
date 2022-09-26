#pragma once

#include "util/util.hpp"
#include "improc/improc.hpp"

struct StateData;
class State; // forward declaration

class Supervisor {
  ordered_value conf;
  std::string dataPath;
  ImProc *imProc = nullptr;

  State *currState_;
  StateData *currStateData_;
  QueueFPS<Eigen::Matrix<int16_t, 3, 1>> *ctrlDataQueuePtr;

  std::atomic<bool> startedCtrl{false};
  std::thread ctrlThread;

  // Called within ctrlThread context
  void start();

public:
  Supervisor(ImProc *imProc);
  ~Supervisor();

  void startThread();
  void stopThread();
  bool started();

  StateData getCurrStateData();
  void clearCtrlDataQueue();

  bool procDataAvail();
  Eigen::Matrix<int16_t, 3, 1> step();

  template <typename T> T getPos() {}

  // events
  void hold();
  void generate();
  void merge();
  void split();

  template <typename T> void updateState() {
    delete currState_;
    currState_ = new T(this);
  }
};
