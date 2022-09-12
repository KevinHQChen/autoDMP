#include "improc/improc.hpp"

ImProc::ImProc(ImCap *imCap)
    : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml")),
      confPath(toml::get<std::string>(conf["improc"]["path"])),
      numChans(toml::get<int>(conf["improc"]["numChans"])), imCap(imCap),
      procFrameQueuePtr(new QueueFPS<cv::Mat>("procFramesQueue.txt")),
      tempResultQueueArr({new QueueFPS<cv::Mat>("tempResultsQueue1.txt"),
                          new QueueFPS<cv::Mat>("tempResultsQueue2.txt"),
                          new QueueFPS<cv::Mat>("tempResultsQueue3.txt")}),
      procFrameQueueArr({new QueueFPS<cv::Mat>("procFramesQueue1.txt"),
                         new QueueFPS<cv::Mat>("procFramesQueue2.txt"),
                         new QueueFPS<cv::Mat>("procFramesQueue3.txt")}),
      procDataQArr({new QueueFPS<cv::Point>("procDataQueue1.txt"),
                    new QueueFPS<cv::Point>("procDataQueue2.txt"),
                    new QueueFPS<cv::Point>("procDataQueue3.txt")}) {
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
        // grab most recent raw frame
        tempFrame = preFrame(impConf.getBBox());
        cv::rectangle(tempFrame, cv::Point(impConf.getBBox().x, impConf.getBBox().y),
                      cv::Point(impConf.getBBox().x + impConf.getBBox().width,
                                impConf.getBBox().y + impConf.getBBox().height),
                      cv::Scalar(0, 255, 0), 2, 8, 0);
        for (int ch = 0; ch < numChans; ++ch) {
          // use chanPose to crop preFrame
          tempPreFrame = tempFrame(impConf.getChanBBox()[ch]);
          if (impConf.getRotAngle()[ch] != 0) {
            rotateMat(tempPreFrame, tempProcFrame, impConf.getRotAngle()[ch]);
            tempPreFrame = tempProcFrame(impConf.getRotChanBBox()[ch]);
          }
          tempProcFrame = tempPreFrame;

          // if setup is currently active, use tmplBBox to update tmplImg
          if (startedSetup && ch == 0) {
            tmplFrames[0] = tempProcFrame(impConf.getTmplBBox());
            cv::flip(tmplFrames[0], tmplFrames[1], -1); // 180deg CCW (flip around x & y-axis)
            impConf.setTmplImg(
                std::array<cv::Mat, NUM_TEMPLATES>{tmplFrames[0].clone(), tmplFrames[1].clone()});
          }

          // perform TM for each tmpl rotation, for each channel
          currMaxLoc.reset();
          for (int rot = 0; rot < NUM_TEMPLATES; ++rot) {
            // outputs a 32-bit float matrix to result (we're using normed cross-correlation)
            cv::matchTemplate(tempProcFrame, impConf.getTmplImg()[rot], tempResultFrame[rot],
                              cv::TM_CCOEFF_NORMED);
            cv::threshold(tempResultFrame[rot], tempResultFrame[rot], tmplThres, 255,
                          cv::THRESH_TOZERO);
            cv::minMaxLoc(tempResultFrame[rot], &minVal, &maxVal, &minLoc, &maxLoc,
                          cv::Mat()); // we only need maxVal & maxLoc if we use correlation
            // keep only the maxLoc closest to junction (i.e. with the highest y value)
            if ((maxVal >= tmplThres) && (maxLoc.y >= currMaxLoc.value_or(maxLoc).y))
              currMaxLoc = maxLoc;
          }

          if (currMaxLoc.has_value()) {
            // save timestamp and maxLoc to file
            procDataQArr[ch]->push(*currMaxLoc);
            procDataQArr[ch]->out << currMaxLoc->x << ", " << currMaxLoc->y << "\n";
            // draw tmpl match for each channel
            cv::rectangle(tempProcFrame, *currMaxLoc,
                          cv::Point(currMaxLoc->x + impConf.getTmplImg()[0].cols,
                                    currMaxLoc->y + impConf.getTmplImg()[0].rows),
                          cv::Scalar(0, 255, 0), 2, 8, 0);
          }

          // push processed frame to queue for display
          // tempResultQueueArr[ch]->push(tempResultFrame[i]);
          procFrameQueueArr[ch]->push(tempProcFrame.clone());

          // debug info
          // info("tempPreFrame size: {}", tempPreFrame.size());
          // info("tempProcFrame size: {}", tempProcFrame.size());
          // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
        } // iterate over all channels
        procFrameQueuePtr->push(preFrame.clone());
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


void ImProc::clearTempFrameQueues() {
  for (auto &q : tempResultQueueArr)
    q->clear();
}

void ImProc::clearProcFrameQueues() {
  for (auto &q : procFrameQueueArr)
    q->clear();
}
