#pragma once

#include "imcap/imcap.hpp"
#include "util/util.hpp"

#define MAX_NO 20

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
  mutable std::mutex chanROIMtx, chWidthMtx, numChsMtx, bgSubHistoryMtx, bgSubThresMtx;
  std::vector<RotRect> chROIs_;
  int chWidth_, numChs_, bgSubHistory_;
  double bgSubThres_;

public:
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

  friend class ImProc;
};

class ImProc {

public:
  ImProc(ImCap *imCap, std::shared_ptr<logger> log);
  ~ImProc();

  std::vector<int> yMax;
  ImProcConfig impConf;
  QueueFPS<std::vector<double>> *procData;

  bool started();

  // load/save channel bounding boxes into/from imProcConfig
  void loadConfig();
  void saveConfig();

  void startThread();
  void stopThread();

  void setR(double r[2 * MAX_NO]);

  std::vector<double> getY();

  cv::Mat getProcFrame(int idx);

  void clearData();

  std::vector<bool> directMeasAvail, yState1, yState2;
  bool zeroCross1, zeroCross2;
  std::vector<unsigned char> zeroCross;

private:
  std::shared_ptr<logger> lg;
  ordered_value conf;
  std::string confPath, dataPath;
  std::vector<int> compParams;
  ImCap *imCap;

  std::mutex imProcMtx, yMtx;
  std::atomic<bool> startedImProc{false};
  std::thread procThread;
  SharedBuffer<std::vector<cv::Mat>> procFrameBuf;

  cv::Ptr<cv::BackgroundSubtractor> pBackSub;
  cv::Mat rectElement, preFrame, fgMask, tempFgMask, tempChFgMask;
  std::vector<cv::Point> fgLocs;
  std::vector<cv::Mat> procFrameArr;
  std::vector<std::vector<double>> directFgClstrs, inferredFgClstrs;
  // state is true for direct, false for inferred
  // std::vector<bool> directMeasAvail, yState1, yState2;
  std::vector<double> yDirect1, yDirect2, yInferred1, yInferred2;
  std::vector<double> y, y1, y2, yPrev1, yPrev2;
  int txCooldown1, txCooldown2, i2dOccurrences1, i2dOccurrences2;
  // bool txOccurred1, txOccurred2;
  double r[2 * MAX_NO];

  void start(); // Called within imProcThread context

  double minDist(std::vector<double> &vec, double value);

  void segAndOrientCh(cv::Mat &srcImg, cv::Mat &tmpImg, cv::Mat &destImg, RotRect &chROI,
                      int &chWidth);
  void findClusters(const std::vector<cv::Point> &fgLocs, std::vector<double> &clusters);
  bool anyNonZeroR(std::size_t start, std::size_t end);
  void rstOnZeroCross();
};
