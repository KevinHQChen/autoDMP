#include "improc/improc.hpp"

ImProc::ImProc(ImCap *imCap)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      procFrameQueuePtr(new QueueFPS<cv::Mat>(dataPath + "procFramesQueue.txt")),
      tempResultQueueArr({new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue1.txt"),
                          new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue2.txt"),
                          new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue3.txt")}),
      procFrameQueueArr({new QueueFPS<cv::Mat>(dataPath + "procFramesQueue1.txt"),
                         new QueueFPS<cv::Mat>(dataPath + "procFramesQueue2.txt"),
                         new QueueFPS<cv::Mat>(dataPath + "procFramesQueue3.txt")}),
      procData(new QueueFPS<std::vector<Pose>>(dataPath + "procDataQueue.txt")) {
  for (int ch = 0; ch < impConf.numChs_; ch++)
    procData->out << "time (ms), maxLoc.x (px), maxLoc.y (px)\n";

  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImProc::~ImProc() {
  stopProcThread();
  for (auto &q : procFrameQueueArr)
    delete q;
  for (auto &q : tempResultQueueArr)
    delete q;
}

// load channel/template images, bounding boxes from file
void ImProc::loadConfig() {
  // load bboxes from file
  ordered_value v = toml::parse<toml::discard_comments, tsl::ordered_map>(confPath + "config.toml");
  impConf.from_toml(v);

  // load template images from file
  info("Number of templates: {}", impConf.getNumTmpls());
  impConf.clearTmplImgs();
  for (int i = 0; i < impConf.getNumTmpls(); ++i)
    impConf.setTmplImg(
        cv::imread(confPath + "tmpl" + std::to_string(i) + ".png", cv::IMREAD_GRAYSCALE));
}

void ImProc::saveConfig() {
  // save bboxes to file
  ordered_value newImProcConfig(impConf);
  std::ofstream out(confPath + "config.toml");
  out << toml::format(newImProcConfig);
  out.close();

  // save template images to file
  int i = 0;
  for (auto &tmplImg : impConf.getTmplImg()) {
    cv::imwrite(confPath + "tmpl" + std::to_string(i) + ".png", tmplImg, compParams);
    ++i;
  }
}

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    for (int rot = 0; rot < impConf.getNumTmpls(); ++rot)
      tempResultFrame.push_back(cv::Mat());
    for (int ch = 0; ch < impConf.numChs_; ch++) {
      tempProcFrameArr.push_back(cv::Mat());
      poseData.push_back(Pose());
    }
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopProcThread() {
  if (started()) {
    info("Stopping image processing...");
    startedImProc = false;
    tempResultFrame.clear();
    tempProcFrameArr.clear();
    if (procThread.joinable())
      procThread.join();
    this->clearProcFrameQueues();
    this->clearTempFrameQueues();
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getFrame();
    try {
      auto startTime = high_resolution_clock::now();
      auto stopTime = high_resolution_clock::now();

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        tempPreFrame = preFrame(impConf.chROIs_[ch]);
        switch (impConf.chROIs_[ch].angle) {
        case 0:
          tempProcFrame = tempPreFrame;
          break;
        case 90:
          cv::rotate(tempPreFrame, tempProcFrame, cv::ROTATE_90_COUNTERCLOCKWISE);
          break;
        case -90:
          cv::rotate(tempPreFrame, tempProcFrame, cv::ROTATE_90_CLOCKWISE);
          break;
        default:
          rotateMat(tempPreFrame, tempProcFrame, impConf.chROIs_[ch].angle);
          int x = tempProcFrame.cols / 2 - impConf.chWidth_ / 2;
          int y = 0;
          int height = tempProcFrame.rows;
          tempProcFrame = tempProcFrame(cv::Rect(x, y, impConf.chWidth_, height));
          break;
        }

        // perform TM for each tmpl rotation, for each channel
        currMaxLoc.reset();
        int matchRot = 0;
        for (int rot = 0; rot < impConf.getNumTmpls(); ++rot) {
          // outputs a 32-bit float matrix to result (we're using normed cross-correlation)
          // info("rot angle: {}", rot);
          // info("tempProcFrame size: {}x{}", tempProcFrame.cols, tempProcFrame.rows);
          // info("tmpl size: {}x{}", impConf.getTmplImg()[rot].cols,
          // impConf.getTmplImg()[rot].rows);
          cv::matchTemplate(tempProcFrame, impConf.getTmplImg()[rot], tempResultFrame[rot],
                            cv::TM_CCOEFF_NORMED);
          cv::threshold(tempResultFrame[rot], tempResultFrame[rot], impConf.tmplThres_, 255,
                        cv::THRESH_TOZERO);
          cv::minMaxLoc(tempResultFrame[rot], &minVal, &maxVal, &minLoc, &maxLoc,
                        cv::Mat()); // we only need maxVal & maxLoc if we use correlation
          // keep only the maxLoc closest to junction (i.e. with the highest y value)
          if ((maxVal >= impConf.tmplThres_) && (maxLoc.y >= currMaxLoc.value_or(maxLoc).y)) {
            currMaxLoc = maxLoc;
            matchRot = rot;
          }
        }

        if (currMaxLoc.has_value()) {
          // save timestamp and maxLoc to file
          poseData[ch] = {matchRot, *currMaxLoc, true};
          procData->out << currMaxLoc->x << ", " << currMaxLoc->y << "\n";
          // draw tmpl match for each channel
          cv::rectangle(tempProcFrame, *currMaxLoc,
                        cv::Point(currMaxLoc->x + impConf.getTmplImg()[0].cols,
                                  currMaxLoc->y + impConf.getTmplImg()[0].rows),
                        cv::Scalar(0, 255, 0), 2, 8, 0);
        } else
          poseData[ch].found = false;

        // push processed frame to queue for display
        // tempResultQueueArr[ch]->push(tempResultFrame[i]);
        // procFrameQueueArr[ch]->push(tempProcFrame);
        tempProcFrameArr[ch] = tempProcFrame;

        // debug info
        // info("tempPreFrame size: {}", tempPreFrame.size());
        // info("tempProcFrame size: {}", tempProcFrame.size());
        // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
      } // iterate over all channels
      procData->push(poseData);
      procFrameBuf.set(tempProcFrameArr);
      // procFrameQueuePtr->push(tempFrame);
      stopTime = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stopTime - startTime);
      // info("imProc duration: {}", duration.count());
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
}

bool ImProc::started() { return startedImProc; }

std::vector<cv::Mat> ImProc::getTempFrames() {
  std::vector<cv::Mat> tempFrames;
  for (auto &q : tempResultQueueArr)
    tempFrames.push_back(q->get());
  return tempFrames;
}

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

cv::Mat ImProc::getProcFrame() {
  if (!procFrameQueuePtr->empty())
    return procFrameQueuePtr->get();
  return cv::Mat();
}

void ImProc::clearTempFrameQueues() {
  for (auto &q : tempResultQueueArr)
    q->clear();
}

void ImProc::clearProcFrameQueues() {
  for (auto &q : procFrameQueueArr)
    q->clear();
}

void ImProc::clearProcData() {
  procData->clear();
}
