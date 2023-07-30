#include "gui/windows/ctrl_window.hpp"

namespace gui {

void CtrlWindow::loadEventsFromFile(const std::string &filename) {
  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string token;
    event_bus ev;
    int i = 0;

    // Get the size of the r array
    int r_size = sizeof(ev.r) / sizeof(ev.r[0]);

    while (std::getline(ss, token, ',')) {
      double value = std::stod(token);
      if (i < r_size)
        ev.r[i] = value;
      else if (i == r_size)
        ev.preT = value;
      else if (i == r_size + 1)
        ev.moveT = value;
      else if (i == r_size + 2)
        ev.postT = value;
      i++;
    }

    // scale r
    for (int ch = 0; ch < 2 * no; ++ch) {
      ev.r[ch] = ev.r[ch] / yMax[ch % no];
      std::clamp(ev.r[ch], -1.0, 1.0);
    }

    sv_->addEvent(ev);
  }
}

CtrlWindow::CtrlWindow(Supervisor *sv, ImProc *imProc) : sv_(sv), imProc_(imProc) {
  currEv_ = sv_->nullEv;
  newEv_ = sv_->nullEv;
}

CtrlWindow::~CtrlWindow() {}

void CtrlWindow::render() {
  if (ImGui::GetIO().KeyAlt && ImGui::IsKeyPressed(ImGuiKey_U))
    ctrlSetupVisible_ = !ctrlSetupVisible_;
  if (ctrlSetupVisible_ && ImGui::Begin("Ctrl Setup", &ctrlSetupVisible_)) {
    no = imProc_->impConf.getNumChs();
    if (yMax.empty())
      yMax = imProc_->yMax;

    if (rPixels.empty())
      rPixels.assign(2 * no, 0.0);

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
    renderSysIdDialog();

    ImGui::PopStyleVar();
    ImGui::End();
  } else
    yMax.clear();

  if (ctrlVisible_)
    sv_->startThread();
  else
    sv_->stopThread();
}

void CtrlWindow::renderAddEventDialog() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("Add Events From File")) {
    static std::array<char, 128> filenameBuffer = {};
    ImGui::InputText("Filename", filenameBuffer.data(), filenameBuffer.size());

    if (ImGui::Button("Add Events")) {
      std::string filename(filenameBuffer.data());
      loadEventsFromFile(filename);
    }
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
      ImGui::TableSetupColumn("r [px]");
      // ImGui::TableSetupColumn("preT [s]");
      ImGui::TableSetupColumn("moveT [s]");
      ImGui::TableSetupColumn("postT [s]");
      ImGui::TableHeadersRow();

      for (auto &event : *(sv_->evQueue_)) {
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        for (int ch = 0; ch < 2 * no; ++ch) {
          if (ch % no != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, event.r[ch]);
        }
        ImGui::TableSetColumnIndex(1);
        for (int ch = 0; ch < 2 * no; ++ch) {
          if (ch % no != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, event.r[ch] * yMax[ch % no]);
        }
        // ImGui::TableSetColumnIndex(2);
        // ImGui::Text(txtFmt, event.preT);
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
        ImGui::TableSetupColumn("r [px]");
        // ImGui::TableSetupColumn("preT [s]");
        ImGui::TableSetupColumn("moveT [s]");
        ImGui::TableSetupColumn("postT [s]");
        ImGui::TableHeadersRow();

        currEv_ = sv_->getCurrEv();
        ImGui::TableNextRow();
        ImGui::TableSetColumnIndex(0);
        for (int ch = 0; ch < 2 * no; ++ch) {
          if (ch % no != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, currEv_.r[ch]);
        }
        ImGui::TableSetColumnIndex(1);
        for (int ch = 0; ch < 2 * no; ++ch) {
          if (ch % no != 0)
            ImGui::SameLine();
          ImGui::Text(txtFmt, currEv_.r[ch] * yMax[ch % no]);
        }
        // ImGui::TableSetColumnIndex(2);
        // ImGui::Text(txtFmt, currEv_.preT);
        ImGui::TableSetColumnIndex(2);
        ImGui::Text(txtFmt, currEv_.moveT);
        ImGui::TableSetColumnIndex(3);
        ImGui::Text(txtFmt, currEv_.postT);
        ImGui::EndTable();
      }
      ImGui::Separator();
    }

    {
      int numCols = no + 1;
      ImGui::Text("Current State");
      if (ImGui::BeginTable("currStateTable", numCols, tableFlags)) {
        ImGui::TableSetupColumn("Vector");
        for (int i = 0; i < no; ++i) {
          std::string label = "Chan " + std::to_string(i);
          ImGui::TableSetupColumn(label.c_str());
        }

        ImGui::TableHeadersRow();

        displayArrayNd("y1", std::vector<double>(sv_->supIn.y, sv_->supIn.y + no));
        displayArrayNd("y2", std::vector<double>(sv_->supIn.y + no, sv_->supIn.y + 2 * no));
        // displayArrayNd("yPrev1", imProc_->yPrev1);
        // displayArrayNd("yPrev2", imProc_->yPrev2);
        displayArrayNd("yhat1", std::vector<double>(sv_->supOut.yhat, sv_->supOut.yhat + no));
        displayArrayNd("yhat2",
                       std::vector<double>(sv_->supOut.yhat + no, sv_->supOut.yhat + 2 * no));
        displayArrayNd("ywt1", std::vector<double>(sv_->supOut.ywt, sv_->supOut.ywt + no));
        displayArrayNd("ywt2", std::vector<double>(sv_->supOut.ywt + no, sv_->supOut.ywt + 2 * no));
        // displayArrayNd("y0_1", std::vector<double>(sv_->supIn.y0, sv_->supIn.y0 + no));
        // displayArrayNd("y0_2", std::vector<double>(sv_->supIn.y0 + no, sv_->supIn.y0 + 2 * no));
        displayArrayNd("y_max1", std::vector<double>(sv_->supIn.ymax, sv_->supIn.ymax + no));
        // displayArrayNd("y_max2",
        //                std::vector<double>(sv_->supIn.ymax + no, sv_->supIn.ymax + 2 * no));
        displayArrayNOptD("yDirect1", imProc_->yDirect1);
        displayArrayNOptD("yDirect2", imProc_->yDirect2);
        displayArrayNOptD("yInferred1", imProc_->yInferred1);
        displayArrayNOptD("yInferred2", imProc_->yInferred2);
        displayArrayNb("y1State", imProc_->yState1);
        displayArrayNb("y2State", imProc_->yState2);
        displayArrayNd("zeroCross1", std::vector<double>(sv_->supIn.yo, sv_->supIn.yo + no));
        displayArrayNd("zeroCross2",
                       std::vector<double>(sv_->supIn.yo + no, sv_->supIn.yo + 2 * no));
        displayArrayNd("currTraj1",
                       std::vector<double>(sv_->supOut.currTraj, sv_->supOut.currTraj + no));
        displayArrayNd("currTraj2", std::vector<double>(sv_->supOut.currTraj + no,
                                                        sv_->supOut.currTraj + 2 * no));
        displayArrayNd("u", std::vector<double>(sv_->supOut.u, sv_->supOut.u + no));
        // displayArrayNd("u0", std::vector<double>(sv_->supIn.u0, sv_->supIn.u0 + no));
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

    // for (int ch = 0; ch < 2 * no; ++ch) {
    //   if (!sv_->supIn.enAdapt[ch])
    //     if (ImGui::Button(("Start Online Param Est Ch" + std::to_string(ch)).c_str()))
    //       sv_->supIn.enAdapt[ch] = true;
    //   if (sv_->supIn.enAdapt[ch])
    //     if (ImGui::Button(("Stop Online Param Est Ch" + std::to_string(ch)).c_str()))
    //       sv_->supIn.enAdapt[ch] = false;
    // }

    umax = sv_->supIn.umax[0];
    ImGui::SliderScalar("uMax", ImGuiDataType_Double, &umax, &umaxMin, &umaxMax, "%.1f");
    for (int ch = 0; ch < no; ++ch)
      sv_->supIn.umax[ch] = umax;

    excitationAmp = sv_->supIn.excitation;
    ImGui::SliderScalar("Excitation Amplitude", ImGuiDataType_Double, &excitationAmp,
                        &excitationAmpMin, &excitationAmpMax, "%.1f");
    sv_->supIn.excitation = excitationAmp;

    paramLowerBound = sv_->supIn.p_;
    ImGui::SliderScalar("Param Lower Bound", ImGuiDataType_Double, &paramLowerBound,
                        &paramLowerBoundMin, &paramLowerBoundMax, "%.1f");
    sv_->supIn.p_ = paramLowerBound;

    covModification = sv_->supIn.dPmod_;
    ImGui::SliderScalar("Cov Modification", ImGuiDataType_Double, &covModification,
                        &covModificationMin, &covModificationMax, "%.1f");
    sv_->supIn.dPmod_ = covModification;

    uwt = sv_->supIn.uwt[0] / 0.025;
    ImGui::SliderScalar("uwt", ImGuiDataType_Double, &uwt, &uwtMin, &uwtMax, "%.1f");
    for (int ch = 0; ch < 2 * no; ++ch)
      sv_->supIn.uwt[ch] = 0.025 * uwt;

    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::renderSysIdDialog() {
  if (openAction != -1)
    ImGui::SetNextItemOpen(openAction != 0);
  if (ImGui::TreeNode("SysID Setup")) {
    // TODO add a button to set sysID to true/false (false by default)
    ImGui::Checkbox("Enable SysID", &sv_->sysID);
    ImGui::SliderScalar("Min Value", ImGuiDataType_Double, &minVal_, &minValMin, &minValMax, "%.1f");
    ImGui::SliderScalar("Max Value", ImGuiDataType_Double, &maxVal_, &maxValMin, &maxValMax, "%.1f");
    ImGui::SliderInt("Order", &order_, 1, 10);
    if (ImGui::Button("Generate Excitation Signal"))
      generateExcitationSignal(minVal_, maxVal_, order_);

    if (sv_->excitationSignal_.size() != 0) {
      timeVec_.clear();
      uVec_.clear();
      for (int i = 0; i < sv_->excitationSignal_.size(); ++i) {
        timeVec_.push_back(i * 0.1);
        uVec_.push_back(sv_->excitationSignal_(i));
      }

      if (ImPlot::BeginPlot("Preview")) {
        ImPlot::SetupAxes("time (s)", "Input");
        ImPlot::PlotLine("u", timeVec_.data(), uVec_.data(), sv_->excitationSignal_.size());
        ImPlot::EndPlot();
      }
    }
    ImGui::TreePop();
  }
  ImGui::Separator();
}

void CtrlWindow::generateExcitationSignal(double minVal, double maxVal, int order) {
    info("Generating excitation signal...");
    py::gil_scoped_acquire acquire;
    py::object prbs = py::module::import("prbs").attr("prbs");

    // Call the prbs function with the provided minVal, maxVal, and order parameters
    sv_->excitationSignal_ = prbs(minVal, maxVal, order).cast<Eigen::VectorXd>();

    info("Excitation signal dimensions: {}x{}", sv_->excitationSignal_.rows(), sv_->excitationSignal_.cols());
}

} // namespace gui
