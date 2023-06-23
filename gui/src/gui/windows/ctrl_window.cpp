#include "gui/windows/ctrl_window.hpp"

namespace gui {

CtrlWindow::CtrlWindow(Supervisor *sv) : sv_(sv) {
  currEv_ = sv_->nullEv;
  newEv_ = sv_->nullEv;
}

CtrlWindow::~CtrlWindow() {}

void CtrlWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_U))
    ctrlSetupVisible_ = !ctrlSetupVisible_;
  if (ctrlSetupVisible_ && ImGui::Begin("Ctrl Setup", &ctrlSetupVisible_)) {
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
    renderDropletGenDialog();
    renderEventQueueContents();
    renderSupervisorStatus();
    renderControllerTuningDialog();

    ImGui::PopStyleVar();
    ImGui::End();
  }

  if (ctrlVisible_)
    sv_->startThread();
  else
    sv_->stopThread();
}

void CtrlWindow::renderAddEventDialog() {
  int numCols = 2;
  double rMax = 1.0, rMid = 0.0, rMin = -1.0, moveTMax = 20.0, moveTMin = 0.0;

  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Add Event")) {
    if (ImGui::BeginTable("eventTable", numCols, tableFlags)) {
      ImGui::TableSetupColumn("Property");
      ImGui::TableSetupColumn("Value");
      ImGui::TableHeadersRow();

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%s", "r [%]");
      ImGui::TableSetColumnIndex(1);
      ImGui::PushItemWidth(-FLT_MIN); // Right-aligned (Setup ItemWidth once (more efficient))
      for (int i = 0; i < sv_->no; ++i)
        ImGui::DragScalarN("##r1", ImGuiDataType_Double, &newEv_.r, sv_->no, 0.1f, &rMin, &rMax,
                           dragFmt);
      ImGui::Text(" ");
      for (int i = 0; i < sv_->no; ++i)
        ImGui::DragScalarN("##r2", ImGuiDataType_Double, &newEv_.r + sv_->no, sv_->no, 0.1f, &rMin,
                           &rMax, dragFmt);

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%s", "preT [s]");
      ImGui::TableSetColumnIndex(1);
      ImGui::DragScalar("##preT", ImGuiDataType_Double, &newEv_.preT, 0.1f, &moveTMin, &moveTMax,
                        dragFmt);

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%s", "moveT [s]");
      ImGui::TableSetColumnIndex(1);
      ImGui::DragScalar("##moveT", ImGuiDataType_Double, &newEv_.moveT, 0.1f, &moveTMin, &moveTMax,
                        dragFmt);

      ImGui::TableNextRow();
      ImGui::TableSetColumnIndex(0);
      ImGui::Text("%s", "postT [s]");
      ImGui::TableSetColumnIndex(1);
      ImGui::DragScalar("##postT", ImGuiDataType_Double, &newEv_.postT, 0.1f, &moveTMin, &moveTMax,
                        dragFmt);

      ImGui::EndTable();
    }

    if (ImGui::Button("Add Event"))
      sv_->addEvent(newEv_);
    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderDropletGenDialog() {
  double max = 1.0, min = 0.0;

  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Droplet Generation")) {
    ImGui::DragScalar("Wch", ImGuiDataType_Double, &wCh, 0.1f, &min, &max, dragFmt);
    ImGui::DragScalar("Dneck", ImGuiDataType_Double, &dNeck, 0.1f, &min, &max, dragFmt);
    ImGui::DragScalar("Dplug", ImGuiDataType_Double, &dPlug, 0.1f, &min, &max, dragFmt);
    ImGui::DragScalar("Droplet length", ImGuiDataType_Double, &dropletLength, 0.1f, &min, &max,
                      dragFmt);
    ImGui::DragScalar("Droplet quantity", ImGuiDataType_S32, &numDroplets);
    if (ImGui::Button("Generate Droplets")) {
      for (int i = 0; i < numDroplets; ++i) {
        // 1. pre-gen: get in position for droplet generation
        sv_->addEvent(
            getEvent(std::vector<double>{0, -wCh / 2, -dropletLength, 0, 0, 0}, 0, 10, 10));
        // 2. gen: perform droplet generation
        sv_->addEvent(
            getEvent(std::vector<double>{0, dNeck * 2, -dropletLength, 0, 0, 0}, 0, 10, 0));
        // 3. post-gen: move droplet out of the way for next droplet
        sv_->addEvent(getEvent(std::vector<double>{0, 0, 0, -wCh / 2, 0, dropletLength / 2 - 1.0},
                               0, 10, 10));
        // post-gen cont.: get in position for next droplet generation
        sv_->addEvent(
            getEvent(std::vector<double>{0, 0, 0, dPlug, 0, dropletLength / 2 - 1.0}, 0, 10, 0));
      }
    }
    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderEventQueueContents() {
  int numCols = 4;

  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Event Queue")) {
    if (ImGui::BeginTable("eventQueueTable", numCols, tableFlags)) {
      ImGui::TableSetupColumn("r [%]");
      ImGui::TableSetupColumn("preT [s]");
      ImGui::TableSetupColumn("moveT [s]");
      ImGui::TableSetupColumn("postT [s]");
      ImGui::TableHeadersRow();

      for (auto &event : *(sv_->evQueue_)) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        for (int i = 0; i < 2 * sv_->no; ++i) {
          if (i != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, event.r[i]);
        }
        ImGui::TableSetColumnIndex(1);
        ImGui::Text(txtFmt, event.preT);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text(txtFmt, event.moveT);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text(txtFmt, event.postT);
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
    {
      int numCols = 4;
      ImGui::Text("Current Event");
      if (ImGui::BeginTable("currEventTable", numCols, tableFlags)) {
        ImGui::TableSetupColumn("r [%]");
        ImGui::TableSetupColumn("preT [s]");
        ImGui::TableSetupColumn("moveT [s]");
        ImGui::TableSetupColumn("postT [s]");
        ImGui::TableHeadersRow();

        currEv_ = sv_->getCurrEv();
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        for (int i = 0; i < 2 * sv_->no; ++i) {
          if (i != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, currEv_.r[i]);
        }
        ImGui::TableSetColumnIndex(1);
        ImGui::Text(txtFmt, currEv_.preT);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text(txtFmt, currEv_.moveT);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text(txtFmt, currEv_.postT);
        ImGui::EndTable();
      }
      ImGui::Separator();
    }

    {
      int numCols = 2 * sv_->no + 1;
      ImGui::Text("Current State");
      if (ImGui::BeginTable("currStateTable", numCols, tableFlags)) {
        ImGui::TableSetupColumn("Vector");
        for (int i = 0; i < 2 * sv_->no; ++i)
          ImGui::TableSetupColumn("Chan %d", i);
        ImGui::TableHeadersRow();

        displayArrayNd("y1", std::vector<double>(sv_->supIn.y, sv_->supIn.y + sv_->no));
        displayArrayNd("y2",
                       std::vector<double>(sv_->supIn.y + sv_->no, sv_->supIn.y + 2 * sv_->no));
        displayArrayNd("yhat1", std::vector<double>(sv_->supOut.yhat, sv_->supOut.yhat + sv_->no));
        displayArrayNd("yhat2", std::vector<double>(sv_->supOut.yhat + sv_->no,
                                                    sv_->supOut.yhat + 2 * sv_->no));
        displayArrayNd("y0_1", std::vector<double>(sv_->supIn.y0, sv_->supIn.y0 + sv_->no));
        displayArrayNd("y0_2",
                       std::vector<double>(sv_->supIn.y0 + sv_->no, sv_->supIn.y0 + 2 * sv_->no));
        displayArrayNd("y_max1", std::vector<double>(sv_->supIn.ymax, sv_->supIn.ymax + sv_->no));
        displayArrayNd("y_max2", std::vector<double>(sv_->supIn.ymax + sv_->no,
                                                     sv_->supIn.ymax + 2 * sv_->no));
        displayArrayNd("r1",
                       std::vector<double>(sv_->supOut.currTraj, sv_->supOut.currTraj + sv_->no));
        displayArrayNd("r2", std::vector<double>(sv_->supOut.currTraj + sv_->no,
                                                 sv_->supOut.currTraj + 2 * sv_->no));
        displayArrayNd("u", std::vector<double>(sv_->supOut.u, sv_->supOut.u + sv_->no));
        displayArrayNd("u0", std::vector<double>(sv_->supIn.u0, sv_->supIn.u0 + sv_->no));
        // displayArray3d("y1", sv_->supIn.y);
        // displayArray3d("y2", sv_->supIn.y + NUM_CHANS);
        // displayArray3d("yhat1", sv_->supOut.yhat);
        // displayArray3d("yhat2", sv_->supOut.yhat + NUM_CHANS);
        // displayArray3d("y0_1", sv_->supIn.y0);
        // displayArray3d("y0_2", sv_->supIn.y0 + NUM_CHANS);
        // displayArray3d("y_max1", sv_->supIn.ymax);
        // displayArray3d("y_max2", sv_->supIn.ymax + NUM_CHANS);
        // displayArray3d("r1", sv_->supOut.currTraj);
        // displayArray3d("r2", sv_->supOut.currTraj + NUM_CHANS);
        // displayArray3d("u", sv_->supOut.u);
        // displayArray3d("u0", sv_->supIn.u0);
        ImGui::EndTable();
      }
    }
    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderControllerTuningDialog() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Controller Tuning")) {
    if (!ctrlVisible_)
      if (ImGui::Button("Start Controller"))
        ctrlVisible_ = true;
    if (ctrlVisible_)
      if (ImGui::Button("Stop Controller"))
        ctrlVisible_ = false;

    for (int ch = 0; ch < 2 * sv_->no; ++ch) {
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

    ImGui::TreePop();
  }
  ImGui::Separator();
}

} // namespace gui
