#include "gui/gui.hpp"
#include <cstdio>

GUI::GUI() : conf(toml::parse<toml::discard_comments, tsl::ordered_map>("config/setup.toml"))
           , guiConf(toml::find<guiConfig>(conf, "gui"))
           , imCap(new ImCap())
           // , imProc(new ImProc())
           // , supervisor(new Supervisor())
{
  info("Config type: {}", type_name<decltype(guiConf)>());
  info("Parsed config: {}", toml::find(conf, "gui"));
}

GUI::~GUI() {
  stopGUIThread();
  delete imCap;
}

void GUI::showRawImCap() {
  if(!imCap->started())
    imCap->startCaptureThread();

  // set window to fullscreen
  const ImGuiViewport* viewport = ImGui::GetMainViewport();
  ImGui::SetNextWindowPos(viewport->WorkPos);
  ImGui::SetNextWindowSize(viewport->WorkSize);

  // show most recent frame or repeat previous frame if no new frame available
  rawFrame = imCap->getRawFrame();
  if (!rawFrame.empty()) {
    updateTexture(rawFrame);
    rawWidth = rawFrame.cols;
    rawHeight = rawFrame.rows;
  }
  dear::Begin("Image Capture", &guiConf.startImCap, imCapFlags) && [this]() {
      ImGui::Image((void *)(intptr_t)textureID, ImVec2(rawWidth, rawHeight));
  };
}

ImGuiWrapperReturnType GUI::render() {
  dear::Begin("Menu") && [this]() {
    ImGui::Text("Instructions: TODO");
    ImGui::Text("Help (this might be a button?)");
    dear::MainMenuBar() && [this]() {
      dear::Menu("File") && [this]() { needToQuit = ImGui::MenuItem("Quit"); };
      dear::Menu("Setup") && [this]() { ImGui::MenuItem("Image Capture", nullptr, &guiConf.startImCap); };
      dear::Menu("Debug") && [this]() { ImGui::MenuItem("Show Demo Window", nullptr, &guiConf.showDebug); };
    };
  };

  if (guiConf.startImCap)
    showRawImCap();
  else
    imCap->stopCaptureThread();

  if (guiConf.showDebug)
    ImGui::ShowDemoWindow(&guiConf.showDebug);

  if (needToQuit)
    return 0;

  return {};
}

void GUI::startGUIThread() {
  guiThread = std::thread(&GUI::imguiMain, this);
  guiThread.detach();
}

void GUI::stopGUIThread() {
  info("Stopping GUI thread...");
  needToQuit = true;
  guiThread.join();
}

// void GUI::showProcImCap(bool &startImCap) {}
