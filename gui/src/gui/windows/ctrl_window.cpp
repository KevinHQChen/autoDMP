#include "gui/windows/ctrl_window.hpp"

namespace gui {

CtrlWindow::CtrlWindow(Supervisor *sv) : sv_(sv) {
}

CtrlWindow::~CtrlWindow() {}

void CtrlWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_U))
      ctrlSetupVisible_ = !ctrlSetupVisible_;

  if (ctrlSetupVisible_) {
    if (ImGui::Begin("Ctrl Setup", &ctrlSetupVisible_)) {
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

          for (int row = 0; row < sizeof(GUIEvent::props) / sizeof(GUIEvent::props[0]); ++row) {
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

} // namespace gui
