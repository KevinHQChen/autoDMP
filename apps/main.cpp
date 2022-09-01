// config.hpp will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "cam/cam.hpp"
#include "gui/gui.hpp"
#include "util/util.hpp"

QueueFPS<cv::Mat> rawFramesQueue("rawFramesQueue.txt");
QueueFPS<cv::Mat> preFramesQueue("preFramesQueue.txt");
std::vector<QueueFPS<cv::Mat>*> tempResultsQueues{new QueueFPS<cv::Mat>("tempResultsQueue1.txt"), new QueueFPS<cv::Mat>("tempResultsQueue2.txt"), new QueueFPS<cv::Mat>("tempResultsQueue3.txt")};
std::vector<QueueFPS<cv::Mat>*> procFramesQueues{new QueueFPS<cv::Mat>("procFramesQueue1.txt"), new QueueFPS<cv::Mat>("procFramesQueue2.txt"), new QueueFPS<cv::Mat>("procFramesQueue3.txt")};
std::vector<QueueFPS<cv::Point>*> procDataQueues{new QueueFPS<cv::Point>("procDataQueue1.txt"), new QueueFPS<cv::Point>("procDataQueue2.txt"), new QueueFPS<cv::Point>("procDataQueue3.txt")};

void updateTexture(const cv::Mat &img, GLuint &textureID) {
  // mimic opencv grayscale image in opengl by making each color channel the same
  if (img.empty()) {
    error("Image is empty");
    return;
  }
  cv::Mat tmp;
  cv::merge(std::vector<cv::Mat>{img, img, img}, tmp);

  // update texture
  // glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  // create opengl texture identifier
  glGenTextures(1, &textureID);
  glBindTexture(GL_TEXTURE_2D, textureID);
  // setup filtering parameters for display
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  // upload pixels into texture
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glTexImage2D(
      GL_TEXTURE_2D,      // texture type
      0,                  // pyramid level (for mip-mapping), 0 is the top level
      GL_RGB,             // internal color format to convert to
                          // (https://www.khronos.org/opengl/wiki/Image_Format)
      tmp.cols, tmp.rows, // image width, height
      0,                  // border width in pixels (can be 1 or 0)
      GL_BGR,             // input image format
      GL_UNSIGNED_BYTE,   // image data type (https://www.khronos.org/opengl/wiki/OpenGL_Type)
      tmp.ptr());         // pointer to data
  glGenerateMipmap(GL_TEXTURE_2D);
}

void showImCap(ordered_value &conf, bool &startImCap) {
  static cv::Mat rawFrame = cv::Mat(0, 0, CV_16UC1);
  static GLuint textureID;

  // set window to fullscreen
  static ImGuiWindowFlags imCapFlags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoSavedSettings;
  const ImGuiViewport* viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(viewport->WorkPos);
  ImGui::SetNextWindowSize(viewport->WorkSize);

  // set up video saving
  static int codec = cv::VideoWriter::fourcc('M', 'J', 'P', 'G');    // allows us to save as avi
  // int codec = cv::VideoWriter::fourcc('R', 'G', 'B', 'A');    // allows us to save as avi
  // int codec = cv::VideoWriter::fourcc('a', 'v', 'c', '1');    // allows us to save as mp4
  // int codec = cv::VideoWriter::fourcc('H', '2', '6', '4');    // allows us to save as avi
  // int codec = cv::VideoWriter::fourcc('M', 'S', 'V', 'C');    // allows us to save as avi (David's codec)

  if (!rawFramesQueue.empty_()) {
    rawFrame = rawFramesQueue.get();    // get most recent frame
    if(conf.saveRaw) {
      // open video writer on first frame received to get its type and size
      // TODO: Change framerate at the end of video to be more accurate
      if(!firstRawFrameRcvd) {
        bool isColor = (rawFrame.type() == CV_8UC3);
        rawVideoWriter.open("rawVideo.avi", codec, conf.rawFPS, rawFrame.size(), isColor);
        if (!rawVideoWriter.isOpened()) {
          mOut << "Could not open the output video file for write\n";
          break;
        }
        firstRawFrameRcvd = true;
      }
    }
    imgAnnotation(rawFramesQueue, rawFrame);
    if(conf.saveRaw) {
      if (firstRawFrameRcvd) {
        rawFrame.convertTo(rawFrame8UC1, CV_8UC1, 255.0/65535.0);
        rawVideoWriter.write(rawFrame8UC1);    // video file will be closed and released automatically (via VideoWriter destructor) once we go out of scope
      }
    }

    updateTexture(rawFrame, textureID);
    dear::Begin("Image Capture", &startImCap, imCapFlags) && []() {
      ImGui::Image((void *)(intptr_t)textureID, ImVec2(rawFrame.cols, rawFrame.rows));
    };
  }
}

ImGuiWrapperReturnType appFn() {
  // TODO put all our GUI-related state in conf
  static bool startImCap{false}, needToQuit{false}, debug{false};
  static ordered_value conf =
      toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");

  dear::Begin("Menu") && []() {
    ImGui::Text("Instructions: TODO");
    ImGui::Text("Help (this might be a button?)");
    dear::MainMenuBar() && []() {
      dear::Menu("File") && []() { needToQuit = ImGui::MenuItem("Quit"); };
      dear::Menu("Setup") && []() { ImGui::MenuItem("Image Capture", nullptr, &startImCap); };
      dear::Menu("Debug") && []() { ImGui::MenuItem("Show Demo Window", nullptr, &debug); };
    };
  };

  showImCap(conf, startImCap);

  if (debug)
    ImGui::ShowDemoWindow();
    // ImGui::ShowMetricsWindow();

  if (needToQuit)
    return 0;

  return {};
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int, const char **) {
  ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");
  info("Config type: {}", type_name<decltype(conf)>());
  info("Parsed config: {}", conf);

  QueueFPS<cv::Mat> rawFramesQueue(L"rawFramesQueue.txt");
  QueueFPS<cv::Mat> preFramesQueue(L"preFramesQueue.txt");
  std::thread captureThread(imCap,
                            conf,
                            std::ref(rawFramesQueue),
                            std::ref(preFramesQueue),
                            onlineCam,
                            offlineCam,
                            std::ref(run));

  std::thread renderThread(imgui_main, std::ref(toml::find(conf, "gui")), std::ref(appFn));


  renderThread.join();
}

#if 0
// NOLINTNEXTLINE(bugprone-exception-escape)
int main() {
#if 0
  ordered_value conf = toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml");
  info("Config type: {}", type_name<decltype(conf)>());
  info("Parsed config: {}", conf);

  // VIDEO SOURCE SETUP
  Cam *cam = new Cam(0, conf);
  // start camera (allocate circular buffer to store frames)
  // timerInterval sets the size of buffers needed (min size of 1)
  // multiplied by the framerate(pretty sure)
  // if timerInterval is less than 1000ms, only 1 buffer (i.e .40 frames) is
  // needed
  cam->start((int)(100 / 1000)); // timerInterval of 100ms
  info("Starting image capture...");
  cv::Mat currentImage = cv::Mat(0, 0, CV_16UC1);
  cv::namedWindow("Display window");
  bool imCapSuccess;
  char key;

  // GUI SETUP (replaces the ugliness below)
  // std::thread captureThread(imgui_main, toml::find(config, "gui"), windowFn);

  while (true) {
    imCapSuccess = cam->process(currentImage);
    if (imCapSuccess) {
      cv::imshow("Display window", currentImage);
      key = (char)cv::waitKey(1);
      if (key != -1) {
        cv::destroyAllWindows();
        break;
      }
    } else {
      error("cannot read image\n");
    }
  }
  delete cam;
#endif
}
#endif
