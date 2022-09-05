#pragma once

#include "util/util.hpp"

#include "atcore.h"
#include "atutility.h"

class Cam {
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
  Cam(int cameraIdx, ordered_value conf);

  /* alternate constructor that uses AT_OpenDevice
   */
  //  cam(int cameraIdx, AT_WC* cameraDescriptor);

  /* cam destructor implements Andor SDK3 internal cleanup of data structures:
      - Closes camera handle: AT_Close(AT_H Handle)
      - Library cleanup: AT_FinalizeLibrary()
  */
  ~Cam();
  int setFeatures();

  void start(const int &Ts); // void singleAcq(cv::Mat &image);
  void stop();
  bool process(cv::Mat &image);
  AT_H handle;

private:
  int cameraIndex;
  ordered_value config;
  ordered_value &camConf;
  int returnCode;
  cv::VideoCapture *offlineCam;
  cv::Mat rawImage;
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
