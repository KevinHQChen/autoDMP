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
      procFrameArr.push_back(cv::Mat());
      p.push_back(Pose());
      // TODO add support for non-90-degree channels
      if (ch == 0)
        yMax.push_back(impConf.getChROIs()[ch].height - impConf.getChWidth());
      else
        yMax.push_back(impConf.getChROIs()[ch].width - impConf.getChWidth());
    }
    procThread = std::thread(&ImProc::start, this);
    procThread.detach();
  }
}

void ImProc::stopProcThread() {
  if (started()) {
    info("Stopping image processing...");
    startedImProc = false;
    procFrameArr.clear();
    p.clear();
    yMax.clear();
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

      // update foreground mask based on background threshold,
      // (optionally update background model at preset learning rate),
      // then apply dilation/erosion to remove noise
      pBackSub->apply(preFrame, fgMask, 0);
      cv::dilate(fgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 3);
      cv::erode(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 9);
      cv::dilate(tempFgMask, tempFgMask, cv::Mat(), cv::Point(-1, -1), 3);

      for (int ch = 0; ch < impConf.numChs_; ++ch) {
        // segment each channel
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
        for (int pose = 0; pose < 3; ++pose)
          p[ch].found[pose] = false;
        procFrameArr[ch] = tempRotChFgMask;
        chImProc(ch);
      }
      procData->push(p);
      procFrameBuf.set(procFrameArr);
      stopTime = high_resolution_clock::now();
      auto duration = duration_cast<milliseconds>(stopTime - startTime);
      // info("imProc duration: {}", duration.count());
    } catch (cv::Exception &e) {
      error("Message: {}", e.what());
      error("Type: {}", type_name<decltype(e)>());
    }
  }
}

void ImProc::chImProc(int ch) {
  switch (ch) {
  case 0:
    // pose 0: fg pixel in center column of ch0 closest to junction
    cv::findNonZero(procFrameArr[ch].col(procFrameArr[ch].cols / 2), fgLocs);
    if (!fgLocs.empty()) {
      currLoc =
          *std::max_element(fgLocs.begin(), fgLocs.end(),
                            [](const cv::Point &p1, const cv::Point &p2) { return p1.y < p2.y; });
      p[ch].p[0] = currLoc.y;
      p[ch].found[0] = true;
      currLoc.x = procFrameArr[ch].cols / 2;
      p[ch].loc[0] = currLoc;
    }
    break;

  case 1:
  case 2:
    // pose 0: extension of ch1/2 to ch0
    if (p[0].found[0] && p[0].p[0] < yMax[0]) {
      p[ch].p[0] =
          yMax[ch] + std::sqrt(2 * std::pow(impConf.getChWidth() / 2.0, 2)) + yMax[0] - p[0].p[0];
      p[ch].found[0] = true;
      p[ch].loc[0] = cv::Point(0, 0);

      return;
    } else {
      // pose 1: fg pixel closest to yMaxPos
      cv::findNonZero(
          procFrameArr[ch](cv::Rect(0, yMax[ch], impConf.getChWidth(), impConf.getChWidth())),
          fgLocs);
      if (!fgLocs.empty()) {
        cv::Point yMaxPos = cv::Point(impConf.getChWidth() / 2.0, 0);
        currLoc = *std::min_element(fgLocs.begin(), fgLocs.end(),
                                    [&yMaxPos](const cv::Point &p1, const cv::Point &p2) {
                                      return cv::norm(p1 - yMaxPos) < cv::norm(p2 - yMaxPos);
                                    });
        p[ch].p[1] = yMax[ch] + cv::norm(currLoc - yMaxPos);
        p[ch].found[1] = true;
        currLoc.y = currLoc.y + yMax[ch];
        p[ch].loc[1] = currLoc;
      }

      // pose 2: fg pixel in center column of ch1/2 farthest/closest to junction
      cv::findNonZero(procFrameArr[ch].col(procFrameArr[ch].cols / 2), fgLocs);
      if (!fgLocs.empty()) {
        // fg pixel in center column of ch1/2 farthest from junction
        currLoc =
            *std::min_element(fgLocs.begin(), fgLocs.end(),
                              [](const cv::Point &p1, const cv::Point &p2) { return p1.y < p2.y; });
        p[ch].p[2] = currLoc.y;
        p[ch].found[2] = true;
        currLoc.x = procFrameArr[ch].cols / 2;
        p[ch].loc[2] = currLoc;
        // fg pixel in center column of ch1/2 closest to junction
        currLoc =
            *std::max_element(fgLocs.begin(), fgLocs.end(),
                              [](const cv::Point &p1, const cv::Point &p2) { return p1.y < p2.y; });
        if (currLoc.y < yMax[ch])
          p[ch].p[2] = currLoc.y;
        currLoc.x = procFrameArr[ch].cols / 2;
        p[ch].loc[2] = currLoc;
      }
    }
    break;
  }

  procData->out << currLoc.x << ", " << currLoc.y << "\n"; // save time/loc to file
  cv::circle(procFrameArr[ch], currLoc, 1, cv::Scalar(100, 100, 100), 1);
}

bool ImProc::started() { return startedImProc; }

cv::Mat ImProc::getProcFrame(int idx) { return procFrameBuf.get()[idx]; }

void ImProc::clearProcData() { procData->clear(); }
