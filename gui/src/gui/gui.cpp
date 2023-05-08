#include "gui/gui.hpp"
#include "ctrl/state/state.hpp"

void GUI::imguiConfig() {
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  if (guiConf.keyboardNav)
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; // Enable Keyboard Controls
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;       // Enable Docking
  // io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;     // Enable Multi-Viewport / Platform
  // Windows
}

void GUI::imguiStyle() {
  ImGuiIO &io = ImGui::GetIO();
  (void)io;
  guiConf.startDark ? ImGui::StyleColorsDark() : ImGui::StyleColorsLight();
  // ImGui::StyleColorsClassic();

  ImGuiStyle &style = ImGui::GetStyle();
  style = ImGuiTheme::ShadesOfGray(3.f, 1.f, 1.f); // GrayVariations
  // style = ImGuiTheme::ShadesOfGray(3.f, 1.136f, 0.865f); // GrayVariations_Darker
  style.ScaleAllSizes(guiConf.scale);
}

GUI::GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv)
    : conf(Config::conf), guiConf(Config::guiConf), imCap_(imCap), imProc_(imProc), pump_(pump),
      sv_(sv) {
  sysIDWindow_ = std::make_shared<gui::SysIdWindow>(sv);
  pumpWindow_ = std::make_shared<gui::PumpWindow>(pump);
  imProcWindow_ = std::make_shared<gui::ImProcWindow>(imCap, imProc);
  ctrlWindow_ = std::make_shared<gui::CtrlWindow>(sv);

  // for documentation on runnerParam members, go to the associated header file
  runnerParams.appWindowParams.windowTitle = "autoDMP";
  runnerParams.appWindowParams.windowGeometry.fullScreenMode =
      HelloImGui::FullScreenMode::FullMonitorWorkArea;

  // custom ImGui config/style/fonts
  //   default config: ImGuiDefaultSettings::SetupDefaultImGuiConfig()
  //   default style: ImGuiDefaultSettings::SetupDefaultImGuiStyle()
  runnerParams.callbacks.SetupImGuiConfig = [this]() { imguiConfig(); };
  runnerParams.callbacks.SetupImGuiStyle = [this]() { imguiStyle(); };
  ImGuiTheme::ImGuiTweakedTheme theme;
  theme.Theme = ImGuiTheme::ImGuiTheme_::ImGuiTheme_MaterialFlat;
  runnerParams.imGuiWindowParams.tweakedTheme = theme;
  // runnerParams.callbacks.LoadAdditionalFonts = [this]() { LoadFonts(); };

  // STATUS BAR (disable for now): enable default + custom status bar
  runnerParams.imGuiWindowParams.showStatusBar = false;
  // runnerParams.imGuiWindowParams.showStatusFps = false;
  runnerParams.fpsIdling.enableIdling = false;
  runnerParams.callbacks.ShowStatus = nullptr; // GUI::renderStatusBar;

  // MENU BAR: enable default (showMenu_App, showMenu_View) + custom menu bar
  runnerParams.imGuiWindowParams.showMenuBar = true;
  runnerParams.callbacks.ShowMenus = [this]() { renderMenu(); };

  // define main GUI rendering function
  runnerParams.callbacks.ShowGui = [this]() { render(); };

  // enable docking, create MainDockSpace
  runnerParams.imGuiWindowParams.defaultImGuiWindowType =
      HelloImGui::DefaultImGuiWindowType::ProvideFullScreenDockSpace;
  // runnerParams.imGuiWindowParams.enableViewports = true;

  // define any external addons
  addOnsParams.withNodeEditor = false;
  addOnsParams.withImplot = true;
}

GUI::~GUI() {
  sysIDWindow_.reset();
  pumpWindow_.reset();
  imProcWindow_.reset();
  ctrlWindow_.reset();
}

void GUI::renderMenu() {
  // called by imgui_bundle within BeginMenuBar()/EndMenuBar() context
  if (ImGui::BeginMenu("File")) {
    runnerParams.appShallExit = ImGui::MenuItem("Quit");
    ImGui::EndMenu();
  }
  if (ImGui::BeginMenu("Setup")) {
    // shortcuts are just for show, they are not implemented here (yet)
    ImGui::MenuItem("Start Image Capture", "c", &imProcWindow_->visible_);
    ImGui::MenuItem("Setup Image Processing", "s", &imProcWindow_->improcSetupVisible_);
    ImGui::MenuItem("Start Image Processing", "i", &imProcWindow_->improcVisible_);
    ImGui::MenuItem("Start Pump Setup", "p", &pumpWindow_->visible_);
    ImGui::MenuItem("Setup System ID", "y", &sysIDWindow_->visible_);
    ImGui::MenuItem("Setup Supervisory Control", "u", &ctrlWindow_->ctrlSetupVisible_);
    ImGui::EndMenu();
  }
  if (ImGui::BeginMenu("Debug")) {
    ImGui::MenuItem("Show Demo Window", nullptr, &guiConf.showDebug);
    ImGui::EndMenu();
  }
}

void GUI::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_Q))
    runnerParams.appShallExit = true;
  imProcWindow_->render();
  pumpWindow_->render();
  ctrlWindow_->render();
  sysIDWindow_->render();
  if (guiConf.showDebug) {
    ImGui::ShowDemoWindow(&guiConf.showDebug);
    ImPlot::ShowDemoWindow(&guiConf.showDebug);
  }
}

void GUI::startGUIThread() {
  guiThread = std::thread([this]() { ImmApp::Run(runnerParams, addOnsParams); });
  guiThread.join();
}
