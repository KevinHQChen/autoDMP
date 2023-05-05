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
          sv_->addEvent(getEvent(2, 2, std::array<int, 3>{80, 0, 80}, 10, 5));
        }

        ImGui::TreePop();
      }
      ImGui::Separator();

      if (openAction != -1)
        ImGui::SetNextItemOpen(openAction != 0);
      if (ImGui::TreeNode("Event Queue")) {
        if (ImGui::BeginTable("eventQueueTable", 8, tableFlags)) {
          ImGui::TableSetupColumn("#");
          ImGui::TableSetupColumn("Src State");
          ImGui::TableSetupColumn("Dest State");
          ImGui::TableSetupColumn("Pos [%]");
          ImGui::TableSetupColumn("Move Time [s]");
          ImGui::TableSetupColumn("Hold Time [s]");
          ImGui::TableSetupColumn("chs");
          ImGui::TableSetupColumn("nextChs");
          ImGui::TableHeadersRow();

          int eventNum = 0;
          for (auto &event : *(sv_->evQueue_)) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%d", eventNum);
            ImGui::TableSetColumnIndex(1);
            ImGui::Text("%f", event.srcState);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%f", event.destState);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("(%f, %f, %f)", event.destPos[0], event.destPos[1], event.destPos[2]);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%f", event.moveTime);
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%f", event.holdTime);
            ImGui::TableSetColumnIndex(6);
            ImGui::Text("%d, %d, %d", event.chs[0], event.chs[1], event.chs[2]);
            ImGui::TableSetColumnIndex(7);
            ImGui::Text("%d, %d, %d", event.nextChs[0], event.nextChs[1], event.nextChs[2]);
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
        if (ImGui::BeginTable("currEventTable", 7, tableFlags)) {
          ImGui::TableSetupColumn("srcState");
          ImGui::TableSetupColumn("destState");
          ImGui::TableSetupColumn("Pos [%]");
          ImGui::TableSetupColumn("moveTime [s]");
          ImGui::TableSetupColumn("holdTime [s]");
          ImGui::TableSetupColumn("chs");
          ImGui::TableSetupColumn("nextChs");
          ImGui::TableHeadersRow();

          if (std::memcmp(&sv_->supOut.currEv, &Supervisor::nullEv, sizeof(event_bus)) != 0) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::Text("%f", sv_->supOut.currEv.srcState);
            ImGui::TableSetColumnIndex(2);
            ImGui::Text("%f", sv_->supOut.currEv.destState);
            ImGui::TableSetColumnIndex(3);
            ImGui::Text("(%f, %f, %f)", sv_->supOut.currEv.destPos[0],
                        sv_->supOut.currEv.destPos[1], sv_->supOut.currEv.destPos[2]);
            ImGui::TableSetColumnIndex(4);
            ImGui::Text("%f", sv_->supOut.currEv.moveTime);
            ImGui::TableSetColumnIndex(5);
            ImGui::Text("%f", sv_->supOut.currEv.holdTime);
            ImGui::TableSetColumnIndex(6);
            ImGui::Text("%d, %d, %d", sv_->supOut.currEv.chs[0], sv_->supOut.currEv.chs[1],
                        sv_->supOut.currEv.chs[2]);
            ImGui::TableSetColumnIndex(7);
            ImGui::Text("%d, %d, %d", sv_->supOut.currEv.nextChs[0], sv_->supOut.currEv.nextChs[1],
                        sv_->supOut.currEv.nextChs[2]);
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

          displayArray3d("y", sv_->supIn.y);
          displayArray3d("y_max", sv_->supIn.y_max);
          displayArray3d("y0", sv_->supIn.y_o);
          displayArray3d("u0", sv_->supIn.u_o);
          displayArray3d("u", sv_->supOut.u);
          displayArray3d("yhat", sv_->supOut.yhat);
          displayArray3d("traj", sv_->supOut.currTraj);
          ImGui::EndTable();
        }
        ImGui::TreePop();
      }
      ImGui::Separator();

      if (!ctrlVisible_)
        if (ImGui::Button("Start Controller"))
          ctrlVisible_ = true;
      if (ctrlVisible_)
        if (ImGui::Button("Stop Controller"))
          ctrlVisible_ = false;

      if (!pauseCtrlDataViz)
        if (ImGui::Button("Pause Data Display"))
          pauseCtrlDataViz = true;
      if (pauseCtrlDataViz)
        if (ImGui::Button("Start Data Display"))
          pauseCtrlDataViz = false;
      if (ImGui::Button("Erase Data")) {
        guiTime = 0;
        for (auto &vec : std::vector<ScrollingBuffer>{u0, u1, u2, y0, y1, y2, yhat0, yhat1, yhat2,
                                                      yref0, yref1, yref2})
          vec.Erase();
      }

      ImGui::PopStyleVar();
      ImGui::End();
    }
  }

  if (ctrlVisible_) {
    sv_->startThread();
    if (ImGui::Begin("Ctrl Data", &ctrlVisible_)) {
      if (!pauseCtrlDataViz) {
        guiTime += ImGui::GetIO().DeltaTime;
        u0.AddPoint(guiTime, sv_->supOut.u[0]);
        u1.AddPoint(guiTime, sv_->supOut.u[1]);
        u2.AddPoint(guiTime, sv_->supOut.u[2]);
        y0.AddPoint(guiTime, sv_->supIn.y[0]);
        y1.AddPoint(guiTime, sv_->supIn.y[1]);
        y2.AddPoint(guiTime, sv_->supIn.y[2]);
        yhat0.AddPoint(guiTime, sv_->supOut.yhat[0]);
        yhat1.AddPoint(guiTime, sv_->supOut.yhat[1]);
        yhat2.AddPoint(guiTime, sv_->supOut.yhat[2]);
        yref0.AddPoint(guiTime, sv_->supOut.currTraj[0]);
        yref1.AddPoint(guiTime, sv_->supOut.currTraj[1]);
        yref2.AddPoint(guiTime, sv_->supOut.currTraj[2]);
      }

      ImGui::SliderFloat("History", &history, 1, 60, "%.1f s");
      plotVector3d("##Control Input", "time (s)", "voltage (V)", 0, 250, ctrlVecs, guiTime, history);
      plotVector3d("##Measured Output", "time (s)", "position (px)", 0, 600, measVecs, guiTime, history);
      ImGui::End();
    }
  } else
    sv_->stopThread();
}

} // namespace gui
