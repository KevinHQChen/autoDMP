#pragma once

#include "util/util.hpp"

/* BARTELS */
// Linux headers
#include <errno.h>   // Error integer and strerror() function
#include <fcntl.h>   // Contains file controls like O_RDWR
#include <termios.h> // Contains POSIX terminal control definitions
#include <unistd.h>  // write(), read(), close()
/* BARTELS END */

/* FLUIGENT */
#include "fgt_SDK_Cpp.h" // include wrapper to fgt_SDK dll, functions can also be accessed by loading the dll
/* FLUIGENT END */

class Pump {

public:
  Pump(std::shared_ptr<logger> log); //constructor declaration
  ~Pump(); //destructor declaration 

  // generic variables declaration
  std::vector<float> outputs, prevOutputs; //defines two vectors of type float one being current outputs and other being previous outputs

  /* BARTELS */
  bool valveState[4]{true, true, true, true}, prevValveState[4]{true, true, true, true};
  int freq{0}, prevFreq{0};
  /* BARTELS END */

  // generic function declarations
  bool setOutput(unsigned int pumpIdx, float voltage);
  int getOutput(unsigned int pumpIdx);
  void setOutputs(std::vector<double> u);
  std::string getPumpType();
  int getNumPumps();

  /* BARTELS */
  void setFreq(int freq);
  void getFreq();
  void enableValve(unsigned int valveIdx);
  void disableValve(unsigned int valveIdx);
  void setValve(unsigned int pumpIdx, bool state);
  void getValve(unsigned int pumpIdx);
  bool sendCmd(std::string cmd, int len);
  /* BARTELS END */

private:
  std::shared_ptr<logger> lg;
  std::mutex mutex;
  const std::string pumpType_ = toml::get<std::string>(Config::conf["pump"]["type"]); //pulled from config file setup.toml
  const bool simModeActive = toml::get<bool>(Config::conf["ctrl"]["simMode"]); //pulled from config file setup.toml

  /* FLUIGENT */
  // structures holding controller/instrument identification and details
  fgt_CHANNEL_INFO channelInfo[256];       // each idx represents one channel
                                           // (numPressureChannels)
  fgt_CONTROLLER_INFO controllerInfo[256]; // each idx represents one
                                           // controller (numControllers)
  fgt_INSTRUMENT_TYPE instrumentType[256]; // None for each idx w/o an instrument
  unsigned short SN[256];                  // Zero for each idx w/o an instrument

  unsigned char numControllers = 0;
  unsigned char numPressureChannels = 0;
  // unsigned int pressureIdx = 0;         // pressure
  // channel index, 0 if first detected channel of list
  float minPressure;
  float maxPressure;
  /* FLUIGENT END */

  /* BARTELS */
  int serialPort;
  // Create new termios struct, we call it 'tty' for convention
  // No need for "= {0}" at the end as we'll immediately write the existing
  // config to this struct
  termios tty;
  char *readData{nullptr};
  unsigned char numPumpChannels = 4;
  /* BARTELS END */
};
