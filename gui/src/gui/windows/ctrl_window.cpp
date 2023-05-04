#include "gui/windows/ctrl_window.hpp"

namespace gui {

CtrlWindow::CtrlWindow(Supervisor *sv) : sv_(sv) {}

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

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("%s", "Src State");
          ImGui::TableSetColumnIndex(1);
          ImGui::PushItemWidth(-FLT_MIN); // Right-aligned (Setup ItemWidth once (more efficient))
          ImGui::SliderInt("##srcState", &srcState, 0, 2);

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("%s", "Dest State");
          ImGui::TableSetColumnIndex(1);
          ImGui::SliderInt("##destState", &destState, 0, 2);

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("%s", "Target Pos [%]");
          ImGui::TableSetColumnIndex(1);
          ImGui::SliderInt3("##targetPos", targetPos, 0, 100);

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("%s", "Move Time [s]");
          ImGui::TableSetColumnIndex(1);
          ImGui::SliderInt("##moveTime", &moveTime, 0, 20);

          ImGui::TableNextRow();
          ImGui::TableSetColumnIndex(0);
          ImGui::Text("%s", "Hold Time [s]");
          ImGui::TableSetColumnIndex(1);
          ImGui::SliderInt("##holdTime", &holdTime, 0, 20);

          ImGui::EndTable();
        }

        if (ImGui::Button("Add Event"))
          sv_->addEvent(getEvent(srcState, destState,
                                 std::array<int, 3>{targetPos[0], targetPos[1], targetPos[2]},
                                 moveTime, holdTime));

        ImGui::SliderInt("Droplet length", &dropletLength, 0, 85);
        if (ImGui::Button("Add Droplet Generation Event")) {
          sv_->addEvent(getEvent(0, 1, std::array<int, 3>{100, 0, 0}, 10, 0));
          sv_->addEvent(getEvent(1, 1, std::array<int, 3>{0, 84, 85 - dropletLength}, 10, 5));
          sv_->addEvent(getEvent(1, 2, std::array<int, 3>{0, 100, 85 - dropletLength}, 10, 0));
          sv_->addEvent(getEvent(2, 2, std::array<int, 3>{80, 0, 80}, 10, 0));
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
            ImGui::Text("(%d, %d, %d)", event.destPos[0], event.destPos[1], event.destPos[2]);
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
