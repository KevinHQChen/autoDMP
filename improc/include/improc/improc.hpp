#pragma once

#include "imcap/imcap.hpp"
#include "util/util.hpp"

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

  cv::Ptr<cv::BackgroundSubtractor> pBackSub;
  cv::Mat preFrame, fgMask, tempFgMask, tempChFgMask;
  std::vector<cv::Point> fgLocs;
  std::vector<std::vector<double>> fgClstrs, fgClstrsFull;

  std::vector<cv::Mat> procFrameArr;
  std::vector<double> y, y1, y2, yPrev1, yPrev2;
  double r[2 * NUM_CHANS];

  std::atomic<bool> startedImProc{false};
  std::thread procThread;

  void start(); // Called within imProcThread context
  bool started();

  double minDist(std::vector<double> &vec, double value);

  void segAndOrientCh(cv::Mat &srcImg, cv::Mat &tmpImg, cv::Mat &destImg, RotRect &chROI,
                      int &chWidth);
  void findClusters(const std::vector<cv::Point> &fgLocs, std::vector<double> &clusters);
  bool anyNonZero(std::size_t start, std::size_t end);
  bool anyZeroCross(const std::vector<double> &vec1, const std::vector<double> &vec2);
  void rstOnZeroCross();

  std::mutex dataMtx;

public:
  std::vector<int> yMax;
  ImProcConfig impConf;
  QueueFPS<std::vector<double>> *procData;

  ImProc(ImCap *imCap);
  ~ImProc();

  // load/save template images, channel bounding boxes into imProcConfig
  void loadConfig();
  void saveConfig();

  void startThread();
  void stopThread();

  void setR(double r[2 * NUM_CHANS]);

  cv::Mat getProcFrame(int idx);

  void clearProcData();
};
