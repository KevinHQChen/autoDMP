#include "improc/improc.hpp"

ImProc::ImProc(std::shared_ptr<ImCap> imCap)
    : conf(TOML11_PARSE_IN_ORDER("config/setup.toml")),
      confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])),
      numChans(toml::get<int>(conf["improc"]["numChans"])), imCap(imCap),
      procFrameQueuePtr(new QueueFPS<cv::Mat>(dataPath + "procFramesQueue.txt")),
      tempResultQueueArr({new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue1.txt"),
                          new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue2.txt"),
                          new QueueFPS<cv::Mat>(dataPath + "tempResultsQueue3.txt")}),
      procFrameQueueArr({new QueueFPS<cv::Mat>(dataPath + "procFramesQueue1.txt"),
                         new QueueFPS<cv::Mat>(dataPath + "procFramesQueue2.txt"),
                         new QueueFPS<cv::Mat>(dataPath + "procFramesQueue3.txt")}),
      tmplThres(toml::get<double>(conf["improc"]["tmplThres"])),
      procDataQArr({new QueueFPS<Pose>(dataPath + "procDataQueue1.txt"),
                    new QueueFPS<Pose>(dataPath + "procDataQueue2.txt"),
                    new QueueFPS<Pose>(dataPath + "procDataQueue3.txt")}) {
  for (int ch = 0; ch < numChans; ch++)
    procDataQArr[ch]->out << "time (ms), maxLoc.x (px), maxLoc.y (px)\n";

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
  for (int i = 0; i < NUM_TEMPLATES; ++i)
    impConf.tmplImg_[i] =
        cv::imread(confPath + "tmpl" + std::to_string(i) + ".png", cv::IMREAD_GRAYSCALE);
}

void ImProc::saveConfig() {
  // save bboxes to file
  ordered_value newImProcConfig(impConf);
  std::ofstream out(confPath + "config.toml");
  out << toml::format(newImProcConfig);
  out.close();

  // save template images to file
  for (int i = 0; i < NUM_TEMPLATES; ++i)
    cv::imwrite(confPath + "tmpl" + std::to_string(i) + ".png", impConf.tmplImg_[i], compParams);
}

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    imCap->clearPreFrameQueue();
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopProcThread() {
  if (started()) {
    info("Stopping image processing...");
    startedImProc = false;
    if (procThread.joinable())
      procThread.join();
    imCap->clearPreFrameQueue();
    this->clearProcFrameQueues();
    this->clearTempFrameQueues();
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getPreFrame();
    try {
      if (!preFrame.empty()) {
        auto startTime = high_resolution_clock::now();
        // grab most recent raw frame
        tempFrame = preFrame(impConf.getBBox());
        for (int ch = 0; ch < numChans; ++ch) {
          // use chanPose to crop preFrame
          tempPreFrame = tempFrame(impConf.getChanBBox()[ch]);
          tempProcFrame.release(); // a fresh cv::Mat is needed each time we call cv::rotate
          if (impConf.getRotAngle()[ch] == 90)
            cv::rotate(tempPreFrame, tempProcFrame, cv::ROTATE_90_COUNTERCLOCKWISE);
          else if (impConf.getRotAngle()[ch] == -90)
            cv::rotate(tempPreFrame, tempProcFrame, cv::ROTATE_90_CLOCKWISE);
          else
            tempProcFrame = tempPreFrame;
          cv::rectangle(tempFrame, impConf.getChanBBox()[ch], cv::Scalar::all(0));

          // if setup is currently active, use tmplBBox to update tmplImg
          if (startedSetup && ch == 0) {
            tmplFrames[0] = tempProcFrame(impConf.getTmplBBox());
            cv::flip(tmplFrames[0], tmplFrames[1], -1); // 180deg CCW (flip around x & y-axis)
            impConf.setTmplImg(0, tmplFrames[0].clone());
            impConf.setTmplImg(1, tmplFrames[1].clone());
          }
          if (startedSetup && ch == 1) {
            tmplFrames[2] = tempProcFrame(impConf.getTmplBBox());
            cv::flip(tmplFrames[2], tmplFrames[3], -1); // 180deg CCW (flip around x & y-axis)
            impConf.setTmplImg(2, tmplFrames[2].clone());
            impConf.setTmplImg(3, tmplFrames[3].clone());
          }

          // perform TM for each tmpl rotation, for each channel
          currMaxLoc.reset();
          int matchRot = 0;
          for (int rot = 0; rot < NUM_TEMPLATES; ++rot) {
            // outputs a 32-bit float matrix to result (we're using normed cross-correlation)
            cv::matchTemplate(tempProcFrame, impConf.getTmplImg()[rot], tempResultFrame[rot],
                              cv::TM_CCOEFF_NORMED);
            cv::threshold(tempResultFrame[rot], tempResultFrame[rot], tmplThres, 255,
                          cv::THRESH_TOZERO);
            cv::minMaxLoc(tempResultFrame[rot], &minVal, &maxVal, &minLoc, &maxLoc,
                          cv::Mat()); // we only need maxVal & maxLoc if we use correlation
            // keep only the maxLoc closest to junction (i.e. with the highest y value)
            if ((maxVal >= tmplThres) && (maxLoc.y >= currMaxLoc.value_or(maxLoc).y)) {
              currMaxLoc = maxLoc;
              matchRot = rot;
            }
          }

          if (currMaxLoc.has_value()) {
            // save timestamp and maxLoc to file
            Pose p = {matchRot, *currMaxLoc};
            procDataQArr[ch]->push(p);
            procDataQArr[ch]->out << currMaxLoc->x << ", " << currMaxLoc->y << "\n";
            // draw tmpl match for each channel
            cv::rectangle(tempProcFrame, *currMaxLoc,
                          cv::Point(currMaxLoc->x + impConf.getTmplImg()[0].cols,
                                    currMaxLoc->y + impConf.getTmplImg()[0].rows),
                          cv::Scalar(0, 255, 0), 2, 8, 0);
          }

          // push processed frame to queue for display
          // tempResultQueueArr[ch]->push(tempResultFrame[i]);
          procFrameQueueArr[ch]->push(tempProcFrame);

          // debug info
          // info("tempPreFrame size: {}", tempPreFrame.size());
          // info("tempProcFrame size: {}", tempProcFrame.size());
          // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
        } // iterate over all channels
        procFrameQueuePtr->push(tempFrame);
        auto stopTime = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stopTime - startTime);
        // info("imProc duration: {}", duration.count());
      }
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
}

bool ImProc::started() { return startedImProc; }

void ImProc::setSetupStatus(bool status) { startedSetup = status; }

std::vector<cv::Mat> ImProc::getTempFrames() {
  std::vector<cv::Mat> tempFrames;
  for (auto &q : tempResultQueueArr)
    tempFrames.push_back(q->get());
  return tempFrames;
}

cv::Mat ImProc::getProcFrame(int idx) {
  if (!procFrameQueueArr[idx]->empty())
    return procFrameQueueArr[idx]->get();
  return cv::Mat();
}

cv::Mat ImProc::getProcFrame() {
  if (!procFrameQueuePtr->empty())
    return procFrameQueuePtr->get();
  return cv::Mat();
}

cv::Point ImProc::getProcData(int idx) {
  if (!procDataQArr[idx]->empty())
    return procDataQArr[idx]->get().loc;
  return cv::Point();
}

void ImProc::clearTempFrameQueues() {
  for (auto &q : tempResultQueueArr)
    q->clear();
}

void ImProc::clearProcFrameQueues() {
  for (auto &q : procFrameQueueArr)
    q->clear();
}

void ImProc::clearProcDataQueues() {
  for (auto &q : procDataQArr)
    q->clear();
}
