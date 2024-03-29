#include "pump/pump.hpp"

Pump::Pump(std::shared_ptr<logger> log) : lg(log) {
  lg->info("Initializing pump...");

  if (simModeActive)
    return;

  if (pumpType_ == "FLUIGENT") {
    // detect number/type of instrument controllers and their serial numbers
    numControllers = Fgt_detect(SN, instrumentType);
    lg->info("Number of controllers detected: {}", int(numControllers));

    // only initialize MFCS-EZ (SN is populated sequentially for each detected controller)
    for (unsigned char controllerIdx = 0; controllerIdx < numControllers; controllerIdx++) {
      if (instrumentType[controllerIdx] == fgt_INSTRUMENT_TYPE::MFCS_EZ)
        lg->info("MFCS-EZ instrument detected at index: {}, serial number: {}", int(controllerIdx),
             SN[controllerIdx]);
      else
        SN[controllerIdx] = 0;
    }
    Fgt_initEx(SN);

    // Get total number of initialized pressure channel(s)
    Fgt_get_pressureChannelCount(&numPressureChannels);
    lg->info("Total number of pressure channels: {}", int(numPressureChannels));

    // Get detailed info about all pressure channels
    Fgt_get_pressureChannelsInfo(channelInfo);

    for (unsigned char ch = 0; ch < numPressureChannels; ch++) {
      // initialize data structures
      outputs.push_back(0.0);
      prevOutputs.push_back(0.0);

      // Get pressure limits
      unsigned int idx = channelInfo[ch].index;
      Fgt_get_pressureRange(idx, &minPressure, &maxPressure);
      lg->info("Channel {} max pressure: {} mbar, min pressure: {} mbar", idx, maxPressure,
           minPressure);

      // Calibrate pressure channels (set pressure commands will not be accepted during this time)
      // if (ch == 0) {
      //   std::cout << "Beginning pressure channel calibration, unplug all tubing from pump.\n";
      //   // std::cout << "Press enter to continue...\n";
      //   // getchar();
      // }
      lg->info("Calibrating pressure channel {}", idx);
      Fgt_calibratePressure(idx);
      lg->info("Pressure channel {} successfully calibrated", idx);
    }
  } else if (pumpType_ == "BARTELS") {
    // open serial port and check for errors (refer to:
    // https://blog.mbedded.ninja/programming/operating-systems/linux/linux-serial-ports-using-c-cpp/#overview)
    serialPort = open("/dev/ttyACM0", O_RDWR);
    if (serialPort < 0)
     lg->error("Error {} opening {}: {}", errno, ttyname(serialPort), strerror(errno));

    /* Configure serial port by modifying termios struct */
    // Read in existing settings, and handle any error
    // NOTE: This is important! POSIX states that the struct passed to tcsetattr()
    // must have been initialized with a call to tcgetattr() overwise behaviour
    // is undefined
    if (tcgetattr(serialPort, &tty) != 0)
     lg->error("Error {} from tcgetattr: {}", errno, strerror(errno));

    // control modes
    // Set 8N1 (8 bits/byte, no parity, one stop bit)
    tty.c_cflag &= ~PARENB; // Clear parity bit, disabling parity (most common)
    tty.c_cflag &=
        ~CSTOPB; // Clear stop field, only one stop bit used in communication (most common)
    tty.c_cflag &= ~CSIZE; // Clear all the size bits, then use one of the statements below
    tty.c_cflag |= CS8;    // 8 bits per byte (most common)

    tty.c_cflag &= ~CRTSCTS;       // Disable RTS/CTS hardware flow control (most common)
    tty.c_cflag |= CREAD | CLOCAL; // Turn on READ & ignore ctrl lines (CLOCAL = 1)

    // local modes
    tty.c_lflag &= ~ICANON; // Disable canonical mode
    tty.c_lflag &= ~ECHO;   // Disable echo
    tty.c_lflag &= ~ECHOE;  // Disable erasure
    tty.c_lflag &= ~ECHONL; // Disable new-line echo
    tty.c_lflag &= ~ISIG;   // Disable interpretation of INTR, QUIT and SUSP

    // input modes
    tty.c_iflag &= ~(IXON | IXOFF | IXANY); // Turn off s/w flow ctrl
    tty.c_iflag &= ~(IGNBRK | BRKINT | PARMRK | ISTRIP | INLCR | IGNCR |
                     ICRNL); // Disable any special handling of received bytes

    // output modes
    tty.c_oflag &= ~OPOST; // Prevent special interpretation of output bytes (e.g. newline chars)
    tty.c_oflag &= ~ONLCR; // Prevent conversion of newline to carriage return/line feed
    // tty.c_oflag &= ~OXTABS; // Prevent conversion of tabs to spaces (NOT PRESENT IN LINUX)
    // tty.c_oflag &= ~ONOEOT; // Prevent removal of C-d chars (0x004) in output (NOT PRESENT IN
    // LINUX)

    // read() will block until data is available (exact amount can be set by VMIN)
    // or timeout occurs (set by VTIME)
    tty.c_cc[VMIN] = 0;   // return as soon as any data is received
    tty.c_cc[VTIME] = 10; // wait for up to 1s (10 deciseconds) for data

    // set in/out baud rate to 9600
    cfsetspeed(&tty, B9600);

    // Save tty settings, also checking for error
    if (tcsetattr(serialPort, TCSANOW, &tty) != 0)
     lg->error("Error {} from tcsetattr: {}", errno, strerror(errno));

    lg->info("Pump serial port {} successfully configured", ttyname(serialPort));

    for (unsigned char ch = 0; ch < numPumpChannels; ch++) {
      // initialize data structures
      outputs.push_back(0.0);
      prevOutputs.push_back(0.0);
    }
  }
}

Pump::~Pump() {
  lg->info("Terminating pump...");

  if (simModeActive)
    return;

  if (pumpType_ == "FLUIGENT") {
    // set all pressures to 0mbar before closing
    for (unsigned char ch = 0; ch < numPressureChannels; ch++)
      Fgt_set_pressure(channelInfo[ch].index, 0);

    Fgt_close();
  } else if (pumpType_ == "BARTELS") {
    lg->info("Closing pump serial port {}", ttyname(serialPort));
    close(serialPort);
  }
}

bool Pump::setOutput(unsigned int pumpIdx, float output) {
  std::lock_guard lock(mutex);

  if (pumpType_ == "FLUIGENT") {
    // check pump index
    if (pumpIdx >= numPressureChannels) {
     lg->error("Pump index {} out of range (max {})", pumpIdx, numPressureChannels - 1);
      return false;
    }

    // do nothing if pump output has not changed
    if (output == prevOutputs[pumpIdx])
      return true;

    if (!simModeActive)
      Fgt_set_pressure(pumpIdx, output);

    // sim mode is active or command was successful
    outputs[pumpIdx] = output;
    prevOutputs[pumpIdx] = output;
    lg->info("Pump {} set to {} V.", pumpIdx + 1, output);
    return true;
  } else if (pumpType_ == "BARTELS") {
    int16_t outputInt = static_cast<int16_t>(output);
    // check pump index
    if (pumpIdx >= numPumpChannels) {
     lg->error("Pump index {} out of range (max {})", pumpIdx, numPumpChannels - 1);
      return false;
    }

    // do nothing if pump output has not changed
    if (outputInt == static_cast<int16_t>(prevOutputs[pumpIdx]))
      return true;

    auto pumpCommand = "P" + std::to_string(pumpIdx + 1) + "V" + std::to_string(outputInt) + "\r\n";

    if (!simModeActive && (!sendCmd(pumpCommand, 4) || std::strncmp("OK", readData, 2) != 0)) {
     lg->error("Error setting pump {} to {} V.", pumpIdx + 1, outputInt);
      return false;
    }

    // sim mode is active or command was successful
    outputs[pumpIdx] = output;
    prevOutputs[pumpIdx] = output;
    lg->info("Pump {} set to {} V.", pumpIdx + 1, outputInt);
    return true;
  } else {
   lg->error("Pump type {} not supported", pumpType_);
    return false;
  }
}

void Pump::setFreq(int freq_) {
  if (pumpType_ != "BARTELS")
    return;

  if (freq_ == prevFreq)
    return;

  if (simModeActive) {
    freq = freq_;
    prevFreq = freq_;
    lg->info("Set pump freq to {} Hz.", freq_);
    return;
  }

  std::string freqCommand = "F" + std::to_string(freq_) + "\r\n";

  sendCmd(freqCommand, 4);

  if (std::strncmp("OK", readData, 2) != 0)
   lg->error("Error setting pump freq.");
  else {
    freq = freq_;
    prevFreq = freq_;
    lg->info("Set pump freq to {} Hz.", freq_);
  }
}

void Pump::setValve(unsigned int valveIdx, bool state) {
  if (pumpType_ != "BARTELS")
    return;

  if (simModeActive) {
    valveState[valveIdx] = state;
    lg->info("Valve {} set to {}.", valveIdx, state ? "ON" : "OFF");
    return;
  }

  std::string valveCommand;
  if (state)
    valveCommand = "V" + std::to_string(valveIdx + 1) + "ON\r\n";
  else
    valveCommand = "V" + std::to_string(valveIdx + 1) + "OFF\r\n";

  sendCmd(valveCommand, 4);
  if (std::strncmp("OK", readData, 2) != 0)
   lg->error("Error setting valve {} to {}.", valveIdx + 1, state ? "ON" : "OFF");
  else {
    valveState[valveIdx] = state;
    lg->info("Valve {} set to {}.", valveIdx + 1, state ? "ON" : "OFF");
  }
}

bool Pump::sendCmd(std::string cmd, int len) {
  if (pumpType_ != "BARTELS")
    return false;

  bool ret = true;
  delete readData;
  readData = new char[len];
  if (write(serialPort, cmd.c_str(), cmd.length()) != cmd.length()) {
   lg->error("Error {} from write: {}", errno, strerror(errno));
    ret = false;
  }
  std::this_thread::sleep_for(1ms);
  // tcdrain(serialPort); // delay for output

  if (read(serialPort, readData, len) != len) {
   lg->error("Error {} from read: {}", errno, strerror(errno));
    ret = false;
  }

  lg->info("readData: {}", readData);

  return ret;
}

void Pump::setOutputs(std::vector<double> u) {
  if (pumpType_ == "FLUIGENT")
    for (int ch = 0; ch < numPressureChannels; ++ch)
      setOutput(ch, u[ch]);
  else if (pumpType_ == "BARTELS")
    for (int ch = 0; ch < numPumpChannels; ++ch)
      setOutput(ch, u[ch]);
}

std::string Pump::getPumpType() { return pumpType_; }

int Pump::getNumPumps() {
  if (pumpType_ == "FLUIGENT")
    return numPressureChannels;
  else // BARTELS
    return numPumpChannels;
}
