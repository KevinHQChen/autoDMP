#pragma once

// #define USEFGTPUMP FALSE
#define USEPIEZOPUMP TRUE

#include "util/util.hpp"

// Linux headers
#include <errno.h>   // Error integer and strerror() function
#include <fcntl.h>   // Contains file controls like O_RDWR
#include <termios.h> // Contains POSIX terminal control definitions
#include <unistd.h>  // write(), read(), close()

// #if USEFGTPUMP == TRUE
// #include "fgt_SDK_Cpp.h" // include wrapper to fgt_SDK dll, functions can also be accessed by
// loading the dll #endif

class Pump {
public:
  // pump state (voltages are ints to work with imgui slider)
  std::vector<int> pumpVoltages{0, 0, 0, 0}, prevPumpVoltages{0, 0, 0, 0};
  std::vector<bool> valveOnOff{false, false, false, false}, prevValveOnOff{false, false, false, false};
  std::vector<bool> valveState{false, false, false, false}, prevValveState{false, false, false, false};
  int freq{0}, prevFreq{0};

  Pump();
  ~Pump();

  // pumpIdx, valveIdx are 1-indexed
  void enablePump(unsigned int pumpIdx);
  void disablePump(unsigned int pumpIdx);
  void setVoltage(unsigned int pumpIdx, int16_t voltage);
  int getVoltage(unsigned int pumpIdx);

  void setFreq(int freq);
  void getFreq();

  void enableValve(unsigned int valveIdx);
  void disableValve(unsigned int valveIdx);
  void setValve(unsigned int pumpIdx, bool state);
  void getValve(unsigned int pumpIdx);

  void sendCmd(std::string cmd, int len);
  void sendSigs(Eigen::Matrix<int16_t, 3, 1> u);

private:
  // #if USEFGTPUMP == TRUE
  //   // structures holding controller/instrument identification and details
  //   fgt_CHANNEL_INFO channelInfo[256];       // each idx represents one channel
  //                                            // (numPressureChannels)
  //   fgt_CONTROLLER_INFO controllerInfo[256]; // each idx represents one
  //                                            // controller (numControllers)
  //   fgt_INSTRUMENT_TYPE instrumentType[256]; // None for each idx w/o an instrument
  //   unsigned short SN[256];                  // Zero for each idx w/o an instrument

  //   unsigned char numControllers = 0;
  //   unsigned char numPressureChannels = 0;
  //   // unsigned int pressureIdx = 0;				// pressure
  //   // channel index, 0 if first detected channel of list
  //   float minPressure;
  //   float maxPressure;
  // #endif

#if USEPIEZOPUMP == TRUE
  int serialPort;
  // Create new termios struct, we call it 'tty' for convention
  // No need for "= {0}" at the end as we'll immediately write the existing
  // config to this struct
  termios tty;
  char *readData{nullptr};
  std::mutex mutex;
#endif
};
