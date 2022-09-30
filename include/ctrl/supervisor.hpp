#pragma once

#include "util/util.hpp"
#include "improc/improc.hpp"

struct Event;
struct StateData;
class State; // forward declaration

class Supervisor {
  ordered_value conf;
  std::string dataPath;

  State *currState_ = nullptr;
  StateData *currStateData_;
  QueueFPS<Event *> *eventQueue_;
  QueueFPS<Eigen::Matrix<int16_t, 3, 1>> *ctrlDataQueuePtr;

  std::atomic<bool> startedCtrl{false};
  std::thread ctrlThread;

  // Called within ctrlThread context
  void start();

public:
  Event *currEvent_ = nullptr;
  ImProc *imProc = nullptr;

  Supervisor(ImProc *imProc);
  ~Supervisor();

  std::string getDataPath() const { return dataPath; }

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

  template <typename T> void updateState();
};
