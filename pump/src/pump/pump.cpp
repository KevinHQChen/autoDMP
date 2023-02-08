#include "pump/pump.hpp"

Pump::Pump(bool simModeActive) : simModeActive(simModeActive) {
  if (simModeActive)
    return;

    // #if USEFGTPUMP == TRUE
    //   // detect number/type of instrument controllers and their serial numbers
    //   numControllers = Fgt_detect(SN, instrumentType);
    //   std::cout << "Number of controllers detected: " << int(numControllers) << "\n";

    //   // only initialize MFCS-EZ (SN is populated sequentially for each detected controller)
    //   for (unsigned char controllerIdx = 0; controllerIdx < numControllers; controllerIdx++) {
    //     if (instrumentType[controllerIdx] == fgt_INSTRUMENT_TYPE::MFCS_EZ) {
    //       std::cout << "MFCS-EZ instrument detected at index: " << int(controllerIdx)
    //                 << ", serial number: " << SN[controllerIdx] << "\n";
    //     } else {
    //       SN[controllerIdx] = 0;
    //     }
    //   }
    //   Fgt_initEx(SN);

    //   // Get total number of initialized pressure channel(s)
    //   Fgt_get_pressureChannelCount(&numPressureChannels);
    //   std::cout << "Total number of pressure channels: " << int(numPressureChannels) << "\n";

    //   // Get detailed info about all pressure channels
    //   Fgt_get_pressureChannelsInfo(channelInfo);

    //   for (unsigned char chanIdx = 0; chanIdx < numPressureChannels; chanIdx++) {
    //     // Get pressure limits
    //     unsigned int idx = channelInfo[chanIdx].index;
    //     Fgt_get_pressureRange(idx, &minPressure, &maxPressure);
    //     std::cout << "Channel " << idx << " max pressure: " << maxPressure << " mbar\n";
    //     std::cout << "Channel " << idx << " min pressure: " << minPressure << " mbar\n";

    //     // Calibrate pressure channels (set pressure commands will not be accepted during this
    //     time) if (chanIdx == 0) {
    //       std::cout << "Beginning pressure channel calibration, unplug all tubing from pump.\n";
    //       std::cout << "Press enter to continue...\n";
    //       getchar();
    //     }
    //     std::cout << "Calibrating pressure channel " << idx << "\n";
    //     Fgt_calibratePressure(idx);
    //     std::cout << "Done.\n";
    //   }
    // #endif

#if USEPIEZOPUMP == TRUE
  // open serial port and check for errors (refer to:
  // https://blog.mbedded.ninja/programming/operating-systems/linux/linux-serial-ports-using-c-cpp/#overview)
  serialPort = open("/dev/ttyACM0", O_RDWR);
  if (serialPort < 0)
    error("Error {} opening {}: {}", errno, ttyname(serialPort), strerror(errno));

  /* Configure serial port by modifying termios struct */
  // Read in existing settings, and handle any error
  // NOTE: This is important! POSIX states that the struct passed to tcsetattr()
  // must have been initialized with a call to tcgetattr() overwise behaviour
  // is undefined
  if (tcgetattr(serialPort, &tty) != 0)
    error("Error {} from tcgetattr: {}", errno, strerror(errno));

  // control modes
  // Set 8N1 (8 bits/byte, no parity, one stop bit)
  tty.c_cflag &= ~PARENB; // Clear parity bit, disabling parity (most common)
  tty.c_cflag &= ~CSTOPB; // Clear stop field, only one stop bit used in communication (most common)
  tty.c_cflag &= ~CSIZE;  // Clear all the size bits, then use one of the statements below
  tty.c_cflag |= CS8;     // 8 bits per byte (most common)

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
    error("Error {} from tcsetattr: {}", errno, strerror(errno));

  info("Pump serial port {} successfully configured", ttyname(serialPort));
#endif
}

Pump::~Pump() {
  if (simModeActive)
    return;

    // #if USEFGTPUMP == TRUE
    //   // set all pressures to 0mbar before closing
    //   for (unsigned char chanIdx = 0; chanIdx < numPressureChannels; chanIdx++)
    //     Fgt_set_pressure(channelInfo[chanIdx].index, 0);

    //   Fgt_close();
    // #endif

#if USEPIEZOPUMP == TRUE
  info("Closing pump serial port {}", ttyname(serialPort));
  close(serialPort);
#endif
}

bool Pump::setVoltage(unsigned int pumpIdx, int16_t voltage) {
  std::lock_guard lock(mutex);

  // Check if voltage is different from current voltage, return early if it's the same
  if (voltage == prevPumpVoltages[pumpIdx - 1])
    return true;

  auto presCommand = "P" + std::to_string(pumpIdx) + "V" + std::to_string(voltage) + "\r\n";

  if (!simModeActive && (!sendCmd(presCommand, 4) || std::strncmp("OK", readData, 2) != 0)) {
    error("Error setting pump {} to {} V.", pumpIdx, voltage);
    return false;
  }

  // sim mode is active
  pumpVoltages[pumpIdx - 1] = voltage;
  prevPumpVoltages[pumpIdx - 1] = voltage;
  info("Pump {} set to {} V.", pumpIdx, voltage);
  return true;
}

void Pump::setFreq(int freq_) {
  if (freq_ == prevFreq)
    return;

  if (simModeActive) {
    freq = freq_;
    prevFreq = freq_;
    info("Set pump freq to {} Hz.", freq_);
    return;
  }

  std::string freqCommand = "F" + std::to_string(freq_) + "\r\n";

  sendCmd(freqCommand, 4);

  if (std::strncmp("OK", readData, 2) != 0)
    error("Error setting pump freq.");
  else {
    freq = freq_;
    prevFreq = freq_;
    info("Set pump freq to {} Hz.", freq_);
  }
}

void Pump::setValve(unsigned int valveIdx, bool state) {
  if (valveState[valveIdx - 1] == state)
    return;

  if (simModeActive) {
    valveState[valveIdx - 1] = false;
    info("Valve {} disabled.", valveIdx);
    return;
  }

  std::string valveCommand;
  if (state)
    valveCommand = "V" + std::to_string(valveIdx) + "ON\r\n";
  else
    valveCommand = "V" + std::to_string(valveIdx) + "OFF\r\n";

  sendCmd(valveCommand, 4);
  if (std::strncmp("OK", readData, 2) != 0)
    error("Error setting valve {} to {}.", valveIdx, state ? "ON" : "OFF");
  else {
    valveState[valveIdx - 1] = state;
    info("Valve {} set to {}.", valveIdx, state ? "ON" : "OFF");
  }
}

bool Pump::sendCmd(std::string cmd, int len) {
  bool ret = true;
  delete readData;
  readData = new char[len];
  if (write(serialPort, cmd.c_str(), cmd.length()) != cmd.length()) {
    error("Error {} from write: {}", errno, strerror(errno));
    ret = false;
  }
  std::this_thread::sleep_for(1ms);
  // tcdrain(serialPort); // delay for output

  if (read(serialPort, readData, len) != len) {
    error("Error {} from read: {}", errno, strerror(errno));
    ret = false;
  }

  info("readData: {}", readData);

  return ret;
}

void Pump::sendSigs(Eigen::Matrix<int16_t, 3, 1> u) {
  // channel -> pump mapping
  setVoltage(1, u(0)); // ch1 = P1, P2
  setVoltage(2, u(0)); // ch1 = P1, P2
  setVoltage(3, u(1)); // ch2 = P3
  setVoltage(4, u(2)); // ch3 = P4
}
