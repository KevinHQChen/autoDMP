#include "improc/improc.hpp"

void imCap(ordered_value &conf,
          QueueFPS<cv::Mat>& rawFramesQueue,
          QueueFPS<cv::Mat>& preFramesQueue,
          cam* onlineCam,
          cv::VideoCapture* offlineCam,
          bool& run) {
}

void showImCap(ordered_value &conf, bool &startImCap, bool &needToQuit) {
  static cv::Mat currentImage = cv::Mat(0, 0, CV_16UC1);
  static Cam *cam;
  static GLuint textureID;
  static bool imCapSetup{false};

  // initialize camera libs, features, and start camera (allocate circular buffer to store frames)
  // timerInterval sets the size of buffers needed (min size of 1) multiplied by the
  // framerate(pretty sure). If timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is needed
  if (startImCap && !imCapSetup) {
    cam = new Cam(0, conf);
    cam->start((int)(100 / 1000)); // timerInterval of 100ms
    info("Starting image capture...");
    imCapSetup = true;
  }

  // set window to fullscreen
  static ImGuiWindowFlags imCapFlags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
  const ImGuiViewport* viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(viewport->WorkPos);
  ImGui::SetNextWindowSize(viewport->WorkSize);

  // start image capture
  if (startImCap && imCapSetup) {
    if (cam->process(currentImage))
      updateTexture(currentImage, textureID);
    else
      error("cannot read image\n");
    dear::Begin("Image Capture", &startImCap, imCapFlags) && []() {
      ImGui::Image((void *)(intptr_t)textureID, ImVec2(currentImage.cols, currentImage.rows));
    };
  }

  if (needToQuit)
    delete cam;
}

void imgAnnotation(QueueFPS<cv::Mat>& imgQueue, cv::Mat& img) {
            // adding basic frame annotation
            std::ostringstream label;
            label << std::fixed << std::setprecision(2) << imgQueue.getFPS();

            std::ostringstream label2;
            label2 << imgQueue.counter_();

            cv::putText(img, label.str(), cv::Point(0, 60), cv::FONT_HERSHEY_SIMPLEX, 0.25, cv::Scalar::all(0));
            cv::putText(img, label2.str(), cv::Point(0, 75), cv::FONT_HERSHEY_SIMPLEX, 0.25, cv::Scalar::all(0));
}
