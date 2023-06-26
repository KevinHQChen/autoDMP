#include "gui/gui.hpp"

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

GUI::GUI(ImCap *imCap, ImProc *imProc, Pump *pump, Supervisor *sv, ImGuiConsole &console)
    : conf(Config::conf), guiConf(Config::guiConf), imCap_(imCap), imProc_(imProc), pump_(pump),
      sv_(sv), console_(console) {
  info("Config type: {}", type_name<decltype(Config::guiConf)>());
  info("Parsed config: {}", toml::find(Config::conf, "gui"));

  pumpWindow_ = std::make_shared<gui::PumpWindow>(pump);
  imProcWindow_ = std::make_shared<gui::ImProcWindow>(imCap, imProc);
  ctrlWindow_ = std::make_shared<gui::CtrlWindow>(sv);
  plotWindow_ = std::make_shared<gui::PlotWindow>(imProc, sv);

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

  // MENU BAR: disable default menus (showMenu_App, showMenu_View), define custom menu bar
  // source (menu_statusbar.cpp): HelloImGui::Menu_StatusBar::ShowMenu(RunnerParams &runnerParams))
  runnerParams.imGuiWindowParams.showMenuBar = true;
  runnerParams.imGuiWindowParams.showMenu_App = false;  // default app menu (quit)
  runnerParams.imGuiWindowParams.showMenu_View = false; // default view menu (theme, etc)
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
  pumpWindow_.reset();
  imProcWindow_.reset();
  ctrlWindow_.reset();
  plotWindow_.reset();
}

void GUI::renderMenu() {
  // called by imgui_bundle within BeginMenuBar()/EndMenuBar() context
  if (ImGui::BeginMenu("File")) {
    runnerParams.appShallExit = ImGui::MenuItem("Quit", "q");
    ImGui::MenuItem("Show Demo Window", nullptr, &guiConf.showDebug);
    ImGui::EndMenu();
  }
  if (ImGui::BeginMenu("View")) {
    // shortcuts are just for show, they are not implemented here (yet)
    ImGui::MenuItem("Image Capture", "d", &imProcWindow_->visible_);
    ImGui::MenuItem("Image Processing Setup", "s", &imProcWindow_->improcSetupVisible_);
    ImGui::MenuItem("Image Processing", "j", &imProcWindow_->improcVisible_);
    ImGui::MenuItem("Pump Control", "p", &pumpWindow_->visible_);
    ImGui::MenuItem("Supervisor Setup", "u", &ctrlWindow_->ctrlSetupVisible_);
    ImGui::MenuItem("Real-Time Plot", "v", &plotWindow_->visible_);
    ImGui::EndMenu();
  }
}

void GUI::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_Q))
    runnerParams.appShallExit = true;
  imProcWindow_->render();
  pumpWindow_->render();
  ctrlWindow_->render();
  plotWindow_->render();
  console_.render();
  if (guiConf.showDebug) {
    ImGui::ShowDemoWindow(&guiConf.showDebug);
    ImPlot::ShowDemoWindow(&guiConf.showDebug);
  }
}

void GUI::startThread() {
  guiThread = std::thread([this]() { ImmApp::Run(runnerParams, addOnsParams); });
  guiThread.join();
}
