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
}

void GUI::showCtrlSetup() {
  if (guiConf.startCtrlSetup) {
    if (ImGui::Begin("Ctrl Setup", &guiConf.startCtrlSetup)) {
      // disable tree node indentation
      ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, 0.0f);

      openAction = -1;
      if (ImGui::Button("Open all"))
        openAction = 1;
      ImGui::SameLine();
      if (ImGui::Button("Close all"))
        openAction = 0;
      ImGui::SameLine();
      ImGui::Separator();

      if (openAction != -1)
        ImGui::SetNextItemOpen(openAction != 0);
      if (ImGui::TreeNode("Add Event")) {
        if (ImGui::BeginTable("eventTable", 2, tableFlags)) {
          ImGui::TableSetupColumn("Property");
          ImGui::TableSetupColumn("Value");
          ImGui::TableHeadersRow();

          for (int row = 0; row < 4; ++row) {
            ImGui::TableNextRow();
            if (row == 0) {
              // Setup ItemWidth once (more efficient)
              ImGui::TableSetColumnIndex(1);
              ImGui::PushItemWidth(-FLT_MIN); // Right-aligned
            }
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%s", GUIEvent::props[row].c_str());
            ImGui::TableSetColumnIndex(1);
            if (row < 2)
              ImGui::SliderInt(GUIEvent::props[row].c_str(), currEvent.data[row],
                               GUIEvent::min[row], GUIEvent::max[row]);
            else
              ImGui::SliderInt3(GUIEvent::props[row].c_str(), currEvent.data[row],
                                GUIEvent::min[row], GUIEvent::max[row]);
          }
          ImGui::EndTable();
        }

        Eigen::Vector3d posVec((double)currEvent.pos[0] / 100.0, (double)currEvent.pos[1] / 100.0,
                               (double)currEvent.pos[2] / 100.0);
        Eigen::Vector3d velVec((double)currEvent.vel[0], (double)currEvent.vel[1],
                               (double)currEvent.vel[2]);
        if (ImGui::Button("Add Event")) {
          sv_->addEvent(currEvent.srcState, currEvent.destState, posVec, velVec);
          guiEventQueue.push_back(currEvent);
        }

        ImGui::SliderInt("Droplet length", &dropletLength, 0, 85);
        if (ImGui::Button("Add Droplet Generation Event")) {
          sv_->addEvent(0, 1, Eigen::Vector3d(200 / 100.0, 0, 0), Eigen::Vector3d(10, 0, 0));
          guiEventQueue.push_back(
              GUIEvent(0, 1, Eigen::Vector3d(200, 0, 0), Eigen::Vector3d(10, 0, 0)));

          sv_->addEvent(1, 1, Eigen::Vector3d(0, 84 / 100.0, (85 - dropletLength) / 100.0),
                        Eigen::Vector3d(0, 10, 10));
          guiEventQueue.push_back(GUIEvent(1, 1, Eigen::Vector3d(0, 84, (85 - dropletLength)),
                                           Eigen::Vector3d(0, 10, 10)));

          sv_->addEvent(1, 2, Eigen::Vector3d(0, 200 / 100.0, (85 - dropletLength) / 100.0),
                        Eigen::Vector3d(0, 10, 10));
          guiEventQueue.push_back(GUIEvent(1, 2, Eigen::Vector3d(0, 200, (85 - dropletLength)),
                                           Eigen::Vector3d(0, 10, 10)));

          sv_->addEvent(2, 2, Eigen::Vector3d(80 / 100.0, 0, 80 / 100.0),
                        Eigen::Vector3d(10, 0, 10));
          guiEventQueue.push_back(
              GUIEvent(2, 2, Eigen::Vector3d(80, 0, 80), Eigen::Vector3d(10, 0, 10)));
        }

        ImGui::TreePop();
      }
      ImGui::Separator();

      if (openAction != -1)
        ImGui::SetNextItemOpen(openAction != 0);
      if (ImGui::TreeNode("Event Queue")) {
        // sync gui event queue with supervisor
        for (int i = 0; i < guiEventQueue.size() - sv_->eventQueue_->size(); ++i)
          guiEventQueue.pop_front();

        if (ImGui::BeginTable("eventQueueTable", 5, tableFlags)) {
          ImGui::TableSetupColumn("#");
          for (int col = 0; col < 4; ++col)
            ImGui::TableSetupColumn(GUIEvent::props[col].c_str());
          ImGui::TableHeadersRow();

          int eventNum = 0;
          for (auto &event : guiEventQueue) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", eventNum);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", event.srcState);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%d", event.destState);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("(%d, %d, %d)", event.pos[0], event.pos[1], event.pos[2]);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("(%d, %d, %d)", event.vel[0], event.vel[1], event.vel[2]);
            ++eventNum;
          }
          ImGui::EndTable();
        }
        ImGui::TreePop();
      }
      ImGui::Separator();

      if (openAction != -1)
        ImGui::SetNextItemOpen(openAction != 0);
      if (ImGui::TreeNode("Supervisor Status")) {
        ImGui::Text("Current Event");
        if (ImGui::BeginTable("currEventTable", 4, tableFlags)) {
          for (int col = 0; col < 4; ++col)
            ImGui::TableSetupColumn(GUIEvent::props[col].c_str());
          ImGui::TableHeadersRow();

          if (sv_->currEvent_ != nullptr) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", sv_->currEvent_->srcState);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%d", sv_->currEvent_->destState);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("(%d, %d, %d)", (int)(sv_->currEvent_->destPos[0] * 100),
                        (int)(sv_->currEvent_->destPos[1] * 100),
                        (int)(sv_->currEvent_->destPos[2] * 100));
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("(%d, %d, %d)", (int)sv_->currEvent_->vel[0], (int)sv_->currEvent_->vel[1],
                        (int)sv_->currEvent_->vel[2]);
          }
          ImGui::EndTable();
        }
        ImGui::Separator();

        ImGui::Text("Current State");
        if (ImGui::BeginTable("currStateTable", 4, tableFlags)) {
          ImGui::TableSetupColumn("Vector");
          ImGui::TableSetupColumn("ch0");
          ImGui::TableSetupColumn("ch1");
          ImGui::TableSetupColumn("ch2");
          ImGui::TableHeadersRow();

          if (sv_->currState_ != nullptr) {
            displayVector3d("u", sv_->currState_->u);
            displayVector3d("u_sat", sv_->currState_->usat);
            displayVector3d("u_ref", sv_->currState_->uref);
            displayVector3d("y", sv_->currState_->y);
            displayVector3d("y_ref", sv_->currState_->yref);
            displayVector3d("y_refScale", sv_->currState_->yrefScale);
            displayArray3b("obsv", sv_->currState_->obsv, "Boolean vector of observed channels");
          }
          ImGui::EndTable();
        }
        ImGui::TreePop();
      }
      ImGui::Separator();

      if (!sysIDWindow_->visible_)
        if (ImGui::Button("Start System ID Setup"))
          sysIDWindow_->visible_ = true;
      if (sysIDWindow_->visible_)
        if (ImGui::Button("Stop System ID Setup"))
          sysIDWindow_->visible_ = false;

      if (!guiConf.startCtrl)
        if (ImGui::Button("Start Controller"))
          guiConf.startCtrl = true;
      if (guiConf.startCtrl)
        if (ImGui::Button("Stop Controller"))
          guiConf.startCtrl = false;

      if (!guiConf.pauseCtrlDataViz)
        if (ImGui::Button("Pause Data Display"))
          guiConf.pauseCtrlDataViz = true;
      if (guiConf.pauseCtrlDataViz)
        if (ImGui::Button("Start Data Display"))
          guiConf.pauseCtrlDataViz = false;
      if (ImGui::Button("Erase Data")) {
        guiTime = 0;
        for (auto &vec :
             std::vector<ScrollingBuffer>{u0, u1, u2, du0, du1, du2, y0, y1, y2, yref0, yref1,
                                          yref2, dxhat0, dxhat1, dxhat2, z0, z1, z2})
          vec.Erase();
      }

      ImGui::PopStyleVar();
      ImGui::End();
    }
  }
}

void GUI::showCtrl() {
  if (guiConf.startCtrl) {
    sv_->startThread();
    if (ImGui::Begin("Ctrl Data", &guiConf.startCtrl)) {
      if (!guiConf.pauseCtrlDataViz) {
        guiTime += ImGui::GetIO().DeltaTime;
        u0.AddPoint(guiTime, sv_->currState_->u(0));
        u1.AddPoint(guiTime, sv_->currState_->u(1));
        u2.AddPoint(guiTime, sv_->currState_->u(2));
        du0.AddPoint(guiTime, sv_->currState_->du(0));
        du1.AddPoint(guiTime, sv_->currState_->du(1));
        du2.AddPoint(guiTime, sv_->currState_->du(2));
        y0.AddPoint(guiTime, sv_->currState_->y(0));
        y1.AddPoint(guiTime, sv_->currState_->y(1));
        y2.AddPoint(guiTime, sv_->currState_->y(2));
        yref0.AddPoint(guiTime, sv_->currState_->yref(0));
        yref1.AddPoint(guiTime, sv_->currState_->yref(1));
        yref2.AddPoint(guiTime, sv_->currState_->yref(2));
        // TODO set mdl->dxhat to state->dxhat at each state::step() so we can plot it
        // dxhat0.AddPoint(guiTime, sv_->currState_->dxhat(0));
        // dxhat1.AddPoint(guiTime, sv_->currState_->dxhat(1));
        // dxhat2.AddPoint(guiTime, sv_->currState_->dxhat(2));
        z0.AddPoint(guiTime, sv_->currState_->z(0));
        z1.AddPoint(guiTime, sv_->currState_->z(1));
        z2.AddPoint(guiTime, sv_->currState_->z(2));
      }

      ImGui::SliderFloat("History", &history, 1, 60, "%.1f s");
      plotVector3d("##Control Input", "time (s)", "voltage (V)", 0, 250, ctrlVecs);
      plotVector3d("##Measured Output", "time (s)", "position (px)", 0, 600, measVecs);
      plotVector3d("##State Error, Integral Error", "time (s)", "error (px)", -500, 500, errorVecs);
      ImGui::End();
    }
  } else
    sv_->stopThread();
}

void GUI::renderMenu() {
  // called by imgui_bundle within BeginMenuBar()/EndMenuBar() context
  if (ImGui::BeginMenu("File")) {
    runnerParams.appShallExit = ImGui::MenuItem("Quit");
    ImGui::EndMenu();
  }
  if (ImGui::BeginMenu("Setup")) {
    ImGui::MenuItem("Start Image Capture", "c", &guiConf.startImCap);
    ImGui::MenuItem("Setup Image Processing", "s", &guiConf.startImProcSetup);
    ImGui::MenuItem("Start Image Processing", "i", &guiConf.startImProc);
    ImGui::MenuItem("Start Pump Setup", "p", &pumpWindow_->visible_);
    ImGui::MenuItem("Start Controller Setup", nullptr, &guiConf.startCtrlSetup);
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
  showCtrlSetup();
  showCtrl();
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

void GUI::plotVector3d(const char *plotName, const char *xAx, const char *yAx, double yMin,
                       double yMax, std::vector<std::pair<ScrollingBuffer *, std::string>> &vecs) {
  if (ImPlot::BeginPlot(plotName, ImVec2(-1, 300))) {
    ImPlot::SetupAxes(xAx, yAx); //, implotFlags, implotFlags);
    ImPlot::SetupAxisLimits(ImAxis_X1, guiTime - history, guiTime, ImGuiCond_Always);
    ImPlot::SetupAxisLimits(ImAxis_Y1, yMin, yMax);
    for (auto &vec : vecs)
      ImPlot::PlotLine(vec.second.c_str(), &vec.first->Data[0].x, &vec.first->Data[0].y,
                       vec.first->Data.size(), 0, vec.first->Offset, 2 * sizeof(float));
    ImPlot::EndPlot();
  }
}
