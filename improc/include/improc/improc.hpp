#pragma once

#include "imcap/imcap.hpp"
#include "util/util.hpp"

struct Pose {
  int rot;
  cv::Point loc;
};

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

void rotateMat(cv::Mat &src, cv::Mat &dst, double angle);

class ImProcConfig {
public:
  mutable std::mutex chanROIMtx, chWidthMtx, numTmplsMtx, numChsMtx, tmplThresMtx, tmplMtx;
  std::vector<RotRect> chROIs_;
  int chWidth_, numTmpls_, numChs_;
  double tmplThres_; // template matching threshold: pxIntensity/255 (8-bit) [double]
  std::vector<cv::Mat> tmplImg_;

  ImProcConfig();
  ImProcConfig(const ImProcConfig &other);
  ImProcConfig &operator=(const ImProcConfig &other) {
    chROIs_ = other.getChROIs();
    chWidth_ = other.getChWidth();
    numTmpls_ = other.numTmpls_;
    numChs_ = other.numChs_;
    tmplThres_ = other.tmplThres_;
    return *this;
  }

  std::vector<RotRect> getChROIs() const;
  int getChWidth() const;
  int getNumTmpls() const;
  int getNumChs() const;
  double getTmplThres() const;
  std::vector<cv::Mat> getTmplImg() const;

  void setChROIs(const std::vector<RotRect> &chROIs);
  void setChWidth(int chWidth);
  void setNumTmpls(int numTmpls);
  void setNumChs(int numChs);
  void setTmplThres(double tmplThres);
  void setTmplImg(cv::Mat tmplImg);

  void setTmplImg(const std::vector<cv::Mat> &tmplImg);
  void clearTmplImgs();

  void from_toml(const ordered_value &v);
  ordered_value into_toml() const;
};

class ImProc {
  ordered_value conf;
  std::string confPath, dataPath;
  ImCap *imCap;
  QueueFPS<cv::Mat> *procFrameQueuePtr;
  std::vector<QueueFPS<cv::Mat> *> tempResultQueueArr, procFrameQueueArr;
  std::vector<int> compParams;

  cv::Mat preFrame{0, 0, CV_16UC1}, tempFrame{0, 0, CV_16UC1}, tempPreFrame{0, 0, CV_16UC1},
      tempProcFrame{0, 0, CV_16UC1};

  // for tmpl matching
  std::vector<cv::Mat> tempResultFrame;
  cv::Point minLoc, maxLoc;
  std::optional<cv::Point> currMaxLoc;
  double minVal, maxVal;

  std::atomic<bool> startedImProc{false};
  std::thread procThread;

  // Called within imProcThread context
  void start();

public:
  ImProcConfig impConf;
  std::vector<QueueFPS<Pose> *> procDataQArr;

  ImProc(ImCap *imCap);
  ~ImProc();

  // load/save template images, channel bounding boxes into imProcConfig
  void loadConfig();
  void saveConfig();

  void startProcThread();
  void stopProcThread();
  bool started();

  std::vector<cv::Mat> getTempFrames();
  cv::Mat getProcFrame(int idx);
  cv::Mat getProcFrame();
  cv::Point getProcData(int idx);
  void clearTempFrameQueues();
  void clearProcFrameQueues();
  void clearProcDataQueues();
};
