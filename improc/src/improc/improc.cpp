#include "improc/improc.hpp"

ImProc::ImProc(ImCap *imCap)
    : conf(Config::conf), confPath(toml::get<std::string>(conf["improc"]["confPath"])),
      dataPath(toml::get<std::string>(conf["postproc"]["procDataPath"])), imCap(imCap),
      procData(new QueueFPS<std::vector<Pose>>(dataPath + "procDataQueue.txt")) {
  for (int ch = 0; ch < impConf.numChs_; ch++)
    procData->out << "time (ms), maxLoc.x (px), maxLoc.y (px)\n";

  // save images with proper format PNG, CV_16UC1
  compParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
  compParams.push_back(0);
}

ImProc::~ImProc() {
  stopProcThread();
  delete procData;
}

// load channel/template images, bounding boxes from file
void ImProc::loadConfig() {
  // load bboxes from file
  ordered_value v = toml::parse<toml::discard_comments, tsl::ordered_map>(confPath + "config.toml");
  impConf.from_toml(v);
}

void ImProc::saveConfig() {
  // save bboxes to file
  ordered_value newImProcConfig(impConf);
  std::ofstream out(confPath + "config.toml");
  out << toml::format(newImProcConfig);
  out.close();
}

void ImProc::startProcThread() {
  if (!started()) {
    info("Starting image processing...");
    startedImProc = true;
    pBackSub =
        cv::createBackgroundSubtractorMOG2(impConf.bgSubHistory_, impConf.bgSubThres_, false);

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
    tempProcFrameArr.clear();
    if (procThread.joinable())
      procThread.join();
  }
}

void ImProc::start() {
  while (started()) {
    preFrame = imCap->getFrame();
    try {
      auto startTime = high_resolution_clock::now();
      auto stopTime = high_resolution_clock::now();

      // update background model at preset learning rate,
      // apply dilation/erosion to remove noise
      pBackSub->apply(preFrame, fgMask, 0);
      cv::dilate(fgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 3);
      cv::erode(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 9);
      cv::dilate(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 3);

      // segment each channel
      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        tempChFgMask = tempFgMask(impConf.chROIs_[ch]);
        switch (impConf.chROIs_[ch].angle) {
        case 0:
          tempRotChFgMask = tempChFgMask;
          break;
        case 90:
          cv::rotate(tempChFgMask, tempRotChFgMask, cv::ROTATE_90_COUNTERCLOCKWISE);
          break;
        case -90:
          cv::rotate(tempChFgMask, tempRotChFgMask, cv::ROTATE_90_CLOCKWISE);
          break;
        case 180:
          cv::rotate(tempChFgMask, tempRotChFgMask, cv::ROTATE_180);
          break;
        default:
          rotateMat(tempChFgMask, tempRotChFgMask, impConf.chROIs_[ch].angle);
          int x = tempRotChFgMask.cols / 2 - impConf.chWidth_ / 2;
          int y = 0;
          int height = tempRotChFgMask.rows;
          tempRotChFgMask = tempRotChFgMask(cv::Rect(x, y, impConf.chWidth_, height));
          break;
        }

        // find maxLoc in each channel
        cv::findNonZero(tempRotChFgMask, fgLocs);
        if (!fgLocs.empty()) {
          currMaxLoc = *std::max_element(
              fgLocs.begin(), fgLocs.end(),
              [](const cv::Point &p1, const cv::Point &p2) { return p1.y < p2.y; });
          currMinLoc = *std::min_element(
              fgLocs.begin(), fgLocs.end(),
              [](const cv::Point &p1, const cv::Point &p2) { return p1.y > p2.y; });

          // save timestamp and maxLoc to file
          poseData[ch] = {currMaxLoc, currMinLoc, true};
          procData->out << currMaxLoc.x << ", " << currMaxLoc.y << "\n";
          // draw maxLoc/minLoc for each channel using black circle
          cv::circle(tempRotChFgMask, currMaxLoc, 1, cv::Scalar(100, 100, 100), 1);
          cv::circle(tempRotChFgMask, currMinLoc, 1, cv::Scalar(200, 200, 200), 1);
        } else
          poseData[ch].found = false;

        // push processed frame to queue for display
        tempProcFrameArr[ch] = tempRotChFgMask;

        // debug info
        // info("currPose rotChanBBox: {}", currPose.rotChanBBox[idx]);
      } // iterate over all channels
      procData->push(poseData);
      procFrameBuf.set(tempProcFrameArr);
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

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearProcData() { procData->clear(); }
