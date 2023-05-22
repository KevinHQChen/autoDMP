#pragma once

#include "imcap/imcap.hpp"
#include "util/util.hpp"

struct Pose {
  cv::Point loc[3];
  double p[3];
  bool found[3];
};

void rotateMat(cv::Mat &src, cv::Mat &dst, double angle);

class RotRect : public cv::Rect {
public:
  int angle;
  double chHeight;

  // Constructors
  RotRect() : cv::Rect(), angle(0) {}
  RotRect(int _x, int _y, int _width, int _height, double _chHeight = 0, double _angle = 0)
      : cv::Rect(_x, _y, _width, _height), angle(_angle) {
    switch (angle) {
    case 0:
    case 180:
      chHeight = _height;
      break;
    case 90:
    case -90:
      chHeight = _width;
      break;
    default:
      chHeight = sqrt(_width * _width + _height * _height);
      break;
    }
  }
  RotRect(const cv::Rect &r, double _chHeight = 0, int _angle = 0) : cv::Rect(r), angle(_angle) {
    switch (angle) {
    case 0:
    case 180:
      chHeight = height;
      break;
    case 90:
    case -90:
      chHeight = width;
      break;
    default:
      chHeight = sqrt(width * width + height * height);
      break;
    }
  }
};

class ImProcConfig {
public:
  mutable std::mutex chanROIMtx, chWidthMtx, numChsMtx, bgSubHistoryMtx, bgSubThresMtx;
  std::vector<RotRect> chROIs_;
  int chWidth_, numChs_, bgSubHistory_;
  double bgSubThres_;

  ImProcConfig();
  ImProcConfig(const ImProcConfig &other);
  ImProcConfig &operator=(const ImProcConfig &other) {
    chROIs_ = other.getChROIs();
    chWidth_ = other.getChWidth();
    numChs_ = other.getNumChs();
    bgSubHistory_ = other.getBgSubHistory();
    bgSubThres_ = other.getBgSubThres();
    return *this;
  }

  std::vector<RotRect> getChROIs() const;
  int getChWidth() const;
  int getNumChs() const;
  int getBgSubHistory() const;
  double getBgSubThres() const;

  void setChROIs(const std::vector<RotRect> &chROIs);
  void setChWidth(int chWidth);
  void setNumChs(int numChs);
  void setBgSubHistory(int bgSubHistory);
  void setBgSubThres(double bgSubThres);

  void from_toml(const ordered_value &v);
  ordered_value into_toml() const;
};

class ImProc {
  ordered_value conf;
  std::string confPath, dataPath;
  ImCap *imCap;
  SharedBuffer<std::vector<cv::Mat>> procFrameBuf;

  std::vector<int> compParams;

  cv::Mat preFrame, fgMask, tempFgMask, tempChFgMask, tempRotChFgMask;
  std::vector<cv::Point> fgLocs;
  cv::Point currLoc;

  cv::Ptr<cv::BackgroundSubtractor> pBackSub;

  std::vector<cv::Mat> procFrameArr;
  std::vector<Pose> p;

  std::atomic<bool> startedImProc{false};
  std::thread procThread;
  // Called within imProcThread context
  void start();
  void chImProc(int ch);

public:
  std::vector<int> yMax;
  ImProcConfig impConf;
  QueueFPS<std::vector<Pose>> *procData;

  ImProc(ImCap *imCap);
  ~ImProc();

  // load/save template images, channel bounding boxes into imProcConfig
  void loadConfig();
  void saveConfig();

  void startProcThread();
  void stopProcThread();
  bool started();

  cv::Mat getProcFrame(int idx);

  void clearProcData();
};
