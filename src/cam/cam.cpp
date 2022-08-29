#include "cam/cam.hpp"

cam::cam(int cameraIdx, ordered_value conf)
    : cameraIndex(cameraIdx), config(conf) {
  // init libs
  returnCode = AT_InitialiseLibrary();
  if (returnCode != AT_SUCCESS) {
    error("FAIL: zyla can not initialise library, return {}\n", returnCode);
    return;
  } else
    info("zyla initialized library\n");
  returnCode = AT_InitialiseUtilityLibrary();
  if (returnCode != AT_SUCCESS) {
    error("FAIL: zyla can not initialise utility library, return {}\n",
          returnCode);
    return;
  } else
    info("zyla initialized utility library\n");

  // open handle
  returnCode = AT_Open(cameraIndex, &handle);
  if (returnCode != AT_SUCCESS) {
    error("FAIL: zyla can not open, return {}\n", returnCode);
    return;
  } else
    info("zyla opened with handle {}\n", handle);

  // set camera features
  auto &andorConf = toml::find(conf, "andor");
  auto &boolConf = toml::find(andorConf, "bool");
  for (const auto &[key, val] : boolConf.as_table()) {
    wstring wKey = wstring(key.begin(), key.end());
    const AT_WC *wcKey = wKey.c_str();

    AT_BOOL impl;
    returnCode = AT_IsImplemented(handle, wcKey, &impl);
    AT_BOOL readOnly;
    returnCode = AT_IsReadOnly(handle, wcKey, &readOnly);
    AT_BOOL readable;
    returnCode = AT_IsReadable(handle, wcKey, &readable);
    AT_BOOL writable;
    returnCode = AT_IsWritable(handle, wcKey, &writable);
    if (returnCode != AT_SUCCESS)
      error("FAIL: can not set {}\n", key);
    else {
      info("i{}o{}r{}w{}\t{}", impl, readOnly, readable, writable, key);
      if (writable) {
        returnCode = AT_SetBool(handle, wcKey, toml::get<AT_BOOL>(val));
        if (returnCode != AT_SUCCESS) {
          error("FAIL: zyla can not set bool feature {} to {}, return {}\n",
                key, val, returnCode);
        } else
          info("zyla set bool feature {} to {}\n", key, val);
      }
    }
  }

  auto &intConf = toml::find(andorConf, "int");
  for (const auto &[key, val] : intConf.as_table()) {
    wstring wKey = wstring(key.begin(), key.end());
    const AT_WC *wcKey = wKey.c_str();

    AT_BOOL impl;
    returnCode = AT_IsImplemented(handle, wcKey, &impl);
    AT_BOOL readOnly;
    returnCode = AT_IsReadOnly(handle, wcKey, &readOnly);
    AT_BOOL readable;
    returnCode = AT_IsReadable(handle, wcKey, &readable);
    AT_BOOL writable;
    returnCode = AT_IsWritable(handle, wcKey, &writable);
    if (returnCode != AT_SUCCESS)
      error("FAIL: can not set {}\n", key);
    else {
      info("i{}o{}r{}w{}\t{}", impl, readOnly, readable, writable, key);
      if (writable) {
        returnCode = AT_SetInt(handle, wcKey, toml::get<AT_64>(val));
        if (returnCode != AT_SUCCESS) {
          error("FAIL: zyla can not set int feature {} to {}, return {}\n", key,
                val, returnCode);
        } else
          info("zyla set int feature {} to {}\n", key, val);
      }
    }
  }

  auto &floatConf = toml::find(andorConf, "float");
  for (const auto &[key, val] : floatConf.as_table()) {
    wstring wKey = wstring(key.begin(), key.end());
    const AT_WC *wcKey = wKey.c_str();

    AT_BOOL impl;
    returnCode = AT_IsImplemented(handle, wcKey, &impl);
    AT_BOOL readOnly;
    returnCode = AT_IsReadOnly(handle, wcKey, &readOnly);
    AT_BOOL readable;
    returnCode = AT_IsReadable(handle, wcKey, &readable);
    AT_BOOL writable;
    returnCode = AT_IsWritable(handle, wcKey, &writable);
    if (returnCode != AT_SUCCESS)
      error("FAIL: can not set {}\n", key);
    else {
      info("i{}o{}r{}w{}\t{}", impl, readOnly, readable, writable, key);
      if (writable) {
        returnCode = AT_SetFloat(handle, wcKey, toml::get<double>(val));
        if (returnCode != AT_SUCCESS) {
          error("FAIL: zyla can not set float feature {} to {}, return {}\n",
                key, val, returnCode);
        } else
          info("zyla set float feature {} to {}\n", key, val);
      }
    }
  }

  auto &enumConf = toml::find(andorConf, "enum");
  for (const auto &[key, val] : enumConf.as_table()) {
    wstring wKey = wstring(key.begin(), key.end());
    const AT_WC *wcKey = wKey.c_str();

    string sVal = toml::get<string>(val);
    wstring wVal = wstring(sVal.begin(), sVal.end());
    const AT_WC *wcVal = wVal.c_str();

    AT_BOOL impl;
    returnCode = AT_IsImplemented(handle, wcKey, &impl);
    AT_BOOL readOnly;
    returnCode = AT_IsReadOnly(handle, wcKey, &readOnly);
    AT_BOOL readable;
    returnCode = AT_IsReadable(handle, wcKey, &readable);
    AT_BOOL writable;
    returnCode = AT_IsWritable(handle, wcKey, &writable);
    if (returnCode != AT_SUCCESS)
      error("FAIL: can not set {}\n", key);
    else {
      info("i{}o{}r{}w{}\t{}", impl, readOnly, readable, writable, key);
      if (writable) {
        returnCode = AT_SetEnumString(handle, wcKey, wcVal);
        if (returnCode != AT_SUCCESS) {
          error("FAIL: zyla can not set enum feature {} to {}, return {}\n",
                key, val, returnCode);
        } else
          info("zyla set enum feature {} to {}\n", key, val);
      }
    }
  }
}

cam::~cam() {
  // close handle
  returnCode = AT_Close(handle);
  if (returnCode != AT_SUCCESS)
    error("FAIL: zyla can not close, return {}\n", returnCode);
  else
    info("zyla closed\n");

  // clean up libs
  returnCode = AT_FinaliseLibrary();
  if (returnCode != AT_SUCCESS)
    error("FAIL: zyla can not finalise library, return {}\n", returnCode);
  else
    info("zyla finalised library\n");

  returnCode = AT_FinaliseUtilityLibrary();
  if (returnCode != AT_SUCCESS)
    error("FAIL: zyla can not finalise utility library, return {}\n",
          returnCode);
  else
    info("zyla finalised utility library\n");
}

void cam::start(const int &Ts) {
  //// flush queue and wait buffers;
  returnCode = AT_Flush(handle);

  //// get info
  samplePeriod = Ts;
  // get size of memory to be allocated for storing each acquired image
  returnCode = AT_GetInt(handle, L"ImageSizeBytes", &imageSizeBytes);
  // cast from AT_64 to int so this value can be used in AT_QueueBuffer
  bufferSize = (int)(imageSizeBytes);

  returnCode = AT_GetFloat(handle, L"FrameRate", &frameRate);
  // set the circular buffer size (number of buffers)
  queueLength = (int)(frameRate * samplePeriod + 1); // min length = 1;
  returnCode = AT_GetInt(handle, L"AOIStride", &imageStride);
  returnCode = AT_GetInt(handle, L"AOIWidth", &imageWidth);
  returnCode = AT_GetInt(handle, L"AOILeft", &imageLeft);
  returnCode = AT_GetInt(handle, L"AOIHeight", &imageHeight);
  returnCode = AT_GetInt(handle, L"AOITop", &imageTop);

  info("Frame Size (bytes): {}\n"
       "Framerate: {}\n"
       "ROI width (pixels): {}\n"
       "ROI left (pixels): {}\n"
       "ROI height (pixels): {}\n"
       "ROI top (pixels): {}\n"
       "ROI Stride (bytes): {}\n"
       "Queue Length: {}\n",
       imageSizeBytes, frameRate, imageWidth, imageLeft, imageHeight, imageTop,
       imageStride, queueLength);

  int index = 0; // [0,1,2,3]:[Mono16,Mono12,Mono12Packed,Mono32]
  // AT_WC = andor-specific 16-bit wide character type
  // (used to represent all feature names, enumerated options, and string
  // feature values)
  AT_WC *str = new AT_WC[256];
  // get pixel encoding value index
  // [0,1,2,3]:[Mono16,Mono12,Mono12Packed,Mono32]
  returnCode = AT_GetEnumIndex(handle, L"Pixel Encoding", &index);
  // convert index to string representation
  returnCode =
      AT_GetEnumStringByIndex(handle, L"Pixel Encoding", index, str, 256);
  imageEncode = str;

  //// allocate buffer
  buffers = new unsigned char *[queueLength];
  alignedBuffers = new unsigned char *[queueLength];
  // memory allocated should be aligned on an 8-byte (64-bit) boundary
  // this helps system performance and prevents alignment faults
  for (int i = 0; i < queueLength; i++) {
    // each index of buffers*[] is a uc pointer
    // so for each index, we'll point it to a newly allocated space
    // that holds a single image divided into uc-sized chunks
    buffers[i] =
        new unsigned char[bufferSize + 8]; // add 7 to allow data alignment
    // this magic somehow aligns it to 64-bit
    alignedBuffers[i] = reinterpret_cast<unsigned char *>(
        (reinterpret_cast<unsigned long long>(buffers[i % queueLength]) + 7) &
        ~7);
  }
  //// pass buffers to fifo queue
  // AT_QueueBuffer will configure the corresponding area of memory to store the
  // acquired images
  for (int i = 0; i < queueLength; i++) {
    returnCode = AT_QueueBuffer(handle, alignedBuffers[i], bufferSize);
  }
  //// start acquisition
  returnCode = AT_Command(handle, L"AcquisitionStart");
  std::cerr << "acquisition start returns " << returnCode << "\n";
  accumNumFrames = 0;
}

void cam::stop() {
  // std::lock_guard<std::mutex> lockGuard(mutex);
  //// stop acquisition
  returnCode = AT_Command(handle, L"AcquisitionStop");
  std::cerr << "acquisition stop returns " << returnCode << "\n";
  returnCode = AT_Flush(handle);
  std::cerr << "queue/wait buffer flush returns " << returnCode << "\n";
  //// free buffer
  for (int i = 0; i < queueLength; i++) {
    delete[] buffers[i];
  }
  delete[] buffers;
  delete[] alignedBuffers;
}

int cam::process(cv::Mat &image) {
  // std::lock_guard<std::mutex> lockGuard(mutex);
  // grab buffer
  unsigned char *pointer;
  int size;
  // if a frame is available, store its address in pointer, size in size, and
  // timeout in ms AT_WaitBuffer will put the calling thread to sleep until
  // timeout elapses or an image becomes available (if timeout equals 0 it will
  // only get existing frames without waiting) (0 timeout causes AT_WaitBuffer
  // to hang quite often, setting to 15ms helps a lot) std::cerr << "zyla
  // process wait buffer returns ";
  returnCode = AT_WaitBuffer(handle, &pointer, &size, 15);
  // std::cerr << returnCode << "\n";

  if (returnCode == 0) {
    // re-queue circular buffer
    returnCode = AT_QueueBuffer(
        handle, alignedBuffers[accumNumFrames % queueLength], bufferSize);
    // std::cerr << "re-queue returns " << returnCode << "\n";
    accumNumFrames++;
    std::cerr << "accumNumFrames: " << accumNumFrames << "\n";
    // clean up buffer
    image = cv::Mat(imageHeight, imageWidth, CV_16UC1);
    returnCode = AT_ConvertBuffer(
        pointer, reinterpret_cast<unsigned char *>(image.data), imageWidth,
        imageHeight, imageStride, imageEncode, L"Mono16");
    // std::cerr << "convert returns " << returnCode << "\n";
  }
  return returnCode;
}
