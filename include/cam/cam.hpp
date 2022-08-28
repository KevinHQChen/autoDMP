#pragma once

#include "util/util.hpp"

#include "atcore.h"
#include "atutility.h"

// create a generic type queue class template for storing frames (typename and
// class are equivalent) you specify the type when you use it
template <typename T> class QueueFPS : public std::queue<T> {
public:
  // constructor:
  // chrono uses the concepts of timepoints and durations
  // now() provides a single timepoint (num of ticks since epoch) from a
  // specific clock this has no meaning by itself (unless you want ticks since
  // epoch, in which case use now().time_since_epoch().count()) a difference
  // between 2 timepoints returns a duration for which the number of ticks is
  // given by duration.count()
  QueueFPS(std::string fileName)
      : counter(0), out(fileName), startTime(std::chrono::steady_clock::now()) {
  }

  void push(const T &entry) {
    std::lock_guard<std::mutex> lockGuard(mutex);
    std::queue<T>::push(entry);
    out << std::fixed << std::setprecision(3);
    out << std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::steady_clock::now() - startTime)
               .count()
        << ", ";
    counter += 1;
    if (counter == 1) {
      // Start counting from a second frame (warmup).
      tm.reset();
      tm.start();
    }
    tmSinceLastPush.reset();
    tmSinceLastPush.start();
  }

  T get() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    T entry = this->front();
    this->pop();
    return entry;
  }

  float getFPS() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    tm.stop();
    double fps = counter / tm.getTimeSec();
    tm.start();
    return static_cast<float>(fps);
  }

  int getTimeSinceLastPush() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    tmSinceLastPush.stop();
    double time = tmSinceLastPush.getTimeMilli();
    tmSinceLastPush.start();
    return (int)std::round(time);
  }

  bool empty_() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    return this->empty();
  }

  void clear() {
    std::lock_guard<std::mutex> lockGuard(mutex);
    while (!this->empty())
      this->pop();
  }

  unsigned int counter_() {
    // std::lock_guard<std::mutex> lockGuard(mutex);
    return counter;
  }

  ffstream out;
  std::chrono::time_point<std::chrono::steady_clock> startTime;

private:
  cv::TickMeter tm, tmSinceLastPush;
  std::mutex mutex;
  unsigned int counter;
};

struct camSettings {
  // Bools
  AT_BOOL Overlap;
  AT_BOOL StaticBlemishCorrection;
  AT_BOOL FastAOIFrameRateEnable;
  AT_BOOL FullAOIControl;
  AT_BOOL AlternatingReadoutDirection;
  AT_BOOL RollingShutterGlobalClear;
  AT_BOOL ScanSpeedControlEnable;
  AT_BOOL SensorCooling;
  AT_BOOL MetadataEnable;
  AT_BOOL MetadataTimestamp;
  AT_BOOL MetadataFrameinfo;

  // Ints
  AT_64 FrameCount;
  AT_64 Accumulatecount;
  AT_64 ImageSizeBytes;
  AT_64 TimeStampClock;
  AT_64 TimeStampClockFrequency;

  // Floats
  double ExposureTime;
  double FrameRate;
  double ReadoutTime;
  double RowReadTime;
  double LineScanSpeed;
  double TargetSensorTemperature;
  double SensorTemperature;

  // Enums
  AT_WC *TriggerMode;
  AT_WC *CycleMode;
  AT_WC *ElectronicShutteringMode;
  AT_WC *PixelReadoutRate;
  AT_WC *TemperatureStatus;
  AT_WC *SensorReadoutMode;
  AT_WC *SimplePreAmpGainControl;
  AT_WC *PixelEncoding;
  AT_WC *BitDepth;

  // Area/Region of Interest (AOI/ROI)
  AT_WC *AOIBinning;
  AT_WC *AOILayout;
  AT_64 AOIHBin;
  AT_64 AOIWidth;
  AT_64 AOILeft;
  AT_64 AOIVBin;
  AT_64 AOIHeight;
  AT_64 AOITop;
  AT_64 AOIStride;
};

class cam {
public:
  /* cam constructor implements Andor SDK3 camera initialization process:
      - Library initialization: AT_InitializeLibrary()
          - sets up internal data structures
          - detect attached cameras
      - Opens camera handle: AT_Open(int DeviceIndex, AT_H* Handle)
          - each camera handle refers to a particular camera
          - to open a camera handle, pass in its DeviceIndex,
          and the address of Handle will be updated (pass-by-reference)
  */
  cam(int cameraIdx, toml::value conf);

  /* alternate constructor that uses AT_OpenDevice
   */
  //  cam(int cameraIdx, AT_WC* cameraDescriptor);

  /* cam destructor implements Andor SDK3 internal cleanup of data structures:
      - Closes camera handle: AT_Close(AT_H Handle)
      - Library cleanup: AT_FinalizeLibrary()
  */
  ~cam();

  void start(const int &Ts); // void singleAcq(cv::Mat &image);
  void stop();
  int process(cv::Mat &image);
  AT_H handle;

private:
  int cameraIndex;
  toml::value config;
  int returnCode;
  AT_64 imageSizeBytes;
  AT_64 imageStride;
  AT_64 imageWidth;
  AT_64 imageHeight;
  AT_64 imageTop;
  AT_64 imageLeft;
  AT_WC *imageEncode;
  int bufferSize;
  double frameRate;
  int samplePeriod;
  int queueLength;

  AT_64 accumNumFrames; // should last 1.8e17 seconds before overflow
  unsigned char **buffers;
  unsigned char **alignedBuffers;
  std::mutex mutex;
};
