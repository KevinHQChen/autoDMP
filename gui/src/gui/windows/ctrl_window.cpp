#include "gui/windows/ctrl_window.hpp"

namespace gui {

CtrlWindow::CtrlWindow(Supervisor *sv) : sv_(sv) {}

CtrlWindow::~CtrlWindow() {}

void CtrlWindow::renderAddEventDialog() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Add Event")) {
    if (ImGui::BeginTable("eventTable", 2, tableFlags)) {
      ImGui::TableSetupColumn("Property");
      ImGui::TableSetupColumn("Value");
      ImGui::TableHeadersRow();

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%s", "Target Pos [%]");
      ImGui::TableSetColumnIndex(1);
      ImGui::PushItemWidth(-FLT_MIN); // Right-aligned (Setup ItemWidth once (more efficient))
      ImGui::SliderInt3("##targetPos1", targetPos, 0, 100);
      ImGui::SliderInt3("##targetPos2", targetPos + 3, 0, 100);

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
      sv_->addEvent(
          getEvent(std::array<int, 2 * NUM_CHANS>{targetPos[0], targetPos[1], targetPos[2]},
                   moveTime, holdTime));

    ImGui::SliderInt("Droplet length", &dropletLength, 0, 85);
    ImGui::SliderInt("Number of Droplets", &numDroplets, 0, 85);
    if (ImGui::Button("Generate Droplets")) {
      for (int i = 0; i < numDroplets; ++i) {
        sv_->addEvent(
            getEvent(std::array<int, 2 * NUM_CHANS>{0, -10, -dropletLength, 0, 0, 0}, 10, 10));
        sv_->addEvent(
            getEvent(std::array<int, 2 * NUM_CHANS>{0, 25, -dropletLength, 0, 0, 0}, 5, 0));
        sv_->addEvent(getEvent(
            std::array<int, 2 * NUM_CHANS>{0, 0, 0, -10, 0, dropletLength / 2 - 100}, 10, 10));
        sv_->addEvent(getEvent(
            std::array<int, 2 * NUM_CHANS>{0, 0, 0, 10, 0, dropletLength / 2 - 100}, 5, 0));
      }
    }

    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderEventQueueContents() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Event Queue")) {
    if (ImGui::BeginTable("eventQueueTable", 4, tableFlags)) {
      ImGui::TableSetupColumn("#");
      ImGui::TableSetupColumn("Pos [%]");
      ImGui::TableSetupColumn("Move Time [s]");
      ImGui::TableSetupColumn("Hold Time [s]");
      ImGui::TableHeadersRow();

      int eventNum = 0;
      for (auto &event : *(sv_->evQueue_)) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        ImGui::Text("%d", eventNum);
        ImGui::TableSetColumnIndex(1);
        ImGui::Text("(%f, %f, %f, %f, %f, %f)", event.r[0], event.r[1], event.r[2], event.r[3],
                    event.r[4], event.r[5]);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text("%f", event.moveT);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text("%f", event.postT);
        ++eventNum;
      }
      ImGui::EndTable();
    }
    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderSupervisorStatus() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Supervisor Status")) {
    ImGui::Text("Current Event");
    if (ImGui::BeginTable("currEventTable", 3, tableFlags)) {
      ImGui::TableSetupColumn("Pos [%]");
      ImGui::TableSetupColumn("moveTime [s]");
      ImGui::TableSetupColumn("holdTime [s]");
      ImGui::TableHeadersRow();

      currEv_ = sv_->getCurrEv();
      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("(%f, %f, %f, %f, %f, %f)", currEv_.r[0], currEv_.r[1], currEv_.r[2],
                  currEv_.r[3], currEv_.r[4], currEv_.r[5]);
      ImGui::TableSetColumnIndex(1);
      ImGui::Text("%f", currEv_.moveT);
      ImGui::TableSetColumnIndex(2);
      ImGui::Text("%f", currEv_.postT);
      ImGui::EndTable();
    }
    ImGui::Separator();

    ImGui::Text("Current State");
    if (ImGui::BeginTable("currStateTable", 7, tableFlags)) {
      ImGui::TableSetupColumn("Vector");
      ImGui::TableSetupColumn("ch0_1");
      ImGui::TableSetupColumn("ch1_1");
      ImGui::TableSetupColumn("ch2_1");
      ImGui::TableSetupColumn("ch0_2");
      ImGui::TableSetupColumn("ch1_2");
      ImGui::TableSetupColumn("ch2_2");
      ImGui::TableHeadersRow();

      displayArray6d("y", sv_->supIn.y);
      displayArray6d("yhat", sv_->supOut.yhat);
      displayArray6d("y_max", sv_->supIn.ymax);
      displayArray6d("y0", sv_->supIn.y0);
      displayArray3d("u0", sv_->supIn.u0);
      displayArray3d("u", sv_->supOut.u);
      displayArray6d("traj", sv_->supOut.currTraj);
      ImGui::EndTable();
    }
    ImGui::TreePop();
  }
  ImGui::Separator();
}

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

      renderAddEventDialog();
      renderEventQueueContents();
      renderSupervisorStatus();

      if (!ctrlVisible_)
        if (ImGui::Button("Start Controller"))
          ctrlVisible_ = true;
      if (ctrlVisible_)
        if (ImGui::Button("Stop Controller"))
          ctrlVisible_ = false;

      for (int ch = 0; ch < 2 * sv_->imProc->impConf.getNumChs(); ++ch) {
        if (!sv_->supIn.enAdapt[ch])
          if (ImGui::Button(("Start Online Param Est Ch" + std::to_string(ch)).c_str()))
            sv_->supIn.enAdapt[ch] = true;
        if (sv_->supIn.enAdapt[ch])
          if (ImGui::Button(("Stop Online Param Est Ch" + std::to_string(ch)).c_str()))
            sv_->supIn.enAdapt[ch] = false;
      }

      excitationAmp = sv_->supIn.excitation;
      ImGui::SliderFloat("Excitation Amplitude", &excitationAmp, 0.0f, 10.0f);
      sv_->supIn.excitation = excitationAmp;

      paramLowerBound = sv_->supIn.p_;
      ImGui::SliderFloat("Param Lower Bound", &paramLowerBound, 0.0005f, 0.1f);
      sv_->supIn.p_ = paramLowerBound;

      covModification = sv_->supIn.dPmod_;
      ImGui::SliderFloat("Cov Modification", &covModification, 0.0005f, 0.1f);
      sv_->supIn.dPmod_ = covModification;

      ImGui::PopStyleVar();
      ImGui::End();
    }
  }

  if (ctrlVisible_)
    sv_->startThread();
  else
    sv_->stopThread();
}

} // namespace gui
