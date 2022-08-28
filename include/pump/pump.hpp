#pragma once

#define USEFGTPUMP FALSE
#define USEPIEZOPUMP TRUE

#include "util/util.hpp"

#include <windows.h>

#if USEFGTPUMP == TRUE
#include "fgt_SDK_Cpp.h" // include wrapper to fgt_SDK dll, functions can also be accessed by loading the dll
#endif

class pump {
public:
  pump();
  ~pump();

  void setVoltage(unsigned int chanIdx, int16_t voltage);
  void setFreq(unsigned int freq);
  void sendCmd(std::string cmd, char *data, int len);
  void sendSigs(Eigen::Matrix<int16_t, 3, 1> u);

private:
#if USEFGTPUMP == TRUE
  // structures holding controller/instrument identification and details
  fgt_CHANNEL_INFO channelInfo[256];       // each idx represents one channel
                                           // (numPressureChannels)
  fgt_CONTROLLER_INFO controllerInfo[256]; // each idx represents one
                                           // controller (numControllers)
  fgt_INSTRUMENT_TYPE
      instrumentType[256]; // None for each idx w/o an instrument
  unsigned short SN[256];  // Zero for each idx w/o an instrument

  unsigned char numControllers = 0;
  unsigned char numPressureChannels = 0;
  // unsigned int pressureIdx = 0;				// pressure
  // channel index, 0 if first detected channel of list
  float minPressure;
  float maxPressure;
#endif

#if USEPIEZOPUMP == TRUE
  HANDLE serialHandle;
  DCB serialParams;
  COMMTIMEOUTS timeouts;
  DWORD bytesRead;
  DWORD bytesWritten;
  char readData[4];
  std::mutex mutex;
#endif
};
