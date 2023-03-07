#include "gui/windows/improc_window.hpp"

namespace gui {

ImProcWindow::ImProcWindow(std::shared_ptr<ImCap> imCap, std::shared_ptr<ImProc> imProc)
    : imCap_(imCap), imProc_(imProc) {
  imProcSetupToggle_ = std::make_unique<Toggle>("ImProc Setup", &improcSetupVisible_);
  rawImage_ = std::make_unique<IMMImage>("Raw Image", 0.25);
  procImage_ = std::make_unique<IMMImage>("Proc Image", 1);
  for (int i = 0; i < NUM_TEMPLATES; i++)
    tmplImages_[i] = std::make_unique<IMMImage>("TmpImage: " + std::to_string(i), 1, true);
  for (int i = 0; i < NUM_CHANS; i++)
    chImages_[i] = std::make_unique<IMMImage>("Channel " + std::to_string(i), 1);
}

ImProcWindow::~ImProcWindow() {
  imProcSetupToggle_.reset();
  rawImage_.reset();
  procImage_.reset();
  for (int i = 0; i < NUM_TEMPLATES; i++)
    tmplImages_[i].reset();
  for (int i = 0; i < NUM_CHANS; i++)
    chImages_[i].reset();
}

void ImProcWindow::render() {
  // raw image
  if (ImGui::IsKeyPressed(ImGuiKey_C))
    visible_ = !visible_;
  if (visible_) {
    imCap_->startCaptureThread();
    if (ImGui::Begin("Image Processing", &visible_)) {
      if (!improcVisible_)
        rawImage_->render(imCap_->getRawFrame());
      else
        imCap_->getRawFrame();
      imProcSetupToggle_->render();
      ImGui::End();
    }
  } else
    imCap_->stopCaptureThread();

  // improc setup
  if (ImGui::IsKeyPressed(ImGuiKey_S))
    improcSetupVisible_ = !improcSetupVisible_;
  imProc_->setSetupStatus(improcSetupVisible_);
  if (improcSetupVisible_) {
    if (ImGui::Begin("Image Processing Setup", &improcSetupVisible_)) {
      // update bbox
      int bbox[4] = {imProc_->impConf.getBBox().x, imProc_->impConf.getBBox().y,
                     imProc_->impConf.getBBox().width,
                     imProc_->impConf.getBBox().height}; // x, y, width, height
      ImGui::SliderInt("BBox.x", &bbox[0], 0, 1000);
      ImGui::SliderInt("BBox.y", &bbox[1], 0, 1000);
      ImGui::SliderInt("BBox.width", &bbox[2], 0, 1000);
      ImGui::SliderInt("BBox.height", &bbox[3], 0, 1000);
      imProc_->impConf.setBBox(cv::Rect(bbox[0], bbox[1], bbox[2], bbox[3]));

      // update junction
      int junc[2] = {imProc_->impConf.getJunction().x, imProc_->impConf.getJunction().y};
      junc[0] = bbox[2] / 2;
      ImGui::SliderInt("junction.y", &junc[1], 0, 1000);
      imProc_->impConf.setJunction(cv::Point(junc[0], junc[1]));

      // update chanWidth
      int chanWidth = imProc_->impConf.getChanWidth();
      ImGui::SliderInt("Channel Width", &chanWidth, 0, 100);
      imProc_->impConf.setChanWidth(chanWidth);

      // update rotAngles
      std::vector<int> rotAngles = imProc_->impConf.getRotAngle();
      for (int idx = 0; idx < NUM_CHANS; ++idx) {
        std::string chanWinName = "Channel " + std::to_string(idx);
        if (ImGui::CollapsingHeader(chanWinName.c_str())) {
          std::string rotWinName = "Rot Angle " + std::to_string(idx);
          ImGui::SliderInt(rotWinName.c_str(), &rotAngles[idx], -180, 180);
        }
      }
      imProc_->impConf.setRotAngle(rotAngles);

      // update chanBBox, rotChanBBox (using bbox, junction, chanWidth, rotAngle)
      std::vector<cv::Rect> chanBBoxes = imProc_->impConf.getChanBBox();
      chanBBoxes[0] = cv::Rect(junc[0] - chanWidth / 2, 0, chanWidth, junc[1]);
      chanBBoxes[1] = cv::Rect(0, junc[1] - chanWidth / 2, junc[0], chanWidth);
      chanBBoxes[2] = cv::Rect(junc[0], junc[1] - chanWidth / 2, junc[0], chanWidth);
      imProc_->impConf.setChanBBox(chanBBoxes);

      // update tmplBBox
      int tmplBBox[4] = {imProc_->impConf.getTmplBBox().x, imProc_->impConf.getTmplBBox().y,
                         imProc_->impConf.getTmplBBox().width,
                         imProc_->impConf.getTmplBBox().height}; // x, y, width, height
      ImGui::SliderInt("TmplBBox.x", &tmplBBox[0], 0, 100);
      ImGui::SliderInt("TmplBBox.y", &tmplBBox[1], 0, 1000);
      ImGui::SliderInt("TmplBBox.width", &tmplBBox[2], 0, 100);
      ImGui::SliderInt("TmplBBox.height", &tmplBBox[3], 0, 100);
      imProc_->impConf.setTmplBBox(cv::Rect(tmplBBox[0], tmplBBox[1], tmplBBox[2], tmplBBox[3]));

      // update tmplThres
      float tmplThres = imProc_->tmplThres;
      ImGui::SliderFloat("Tmpl Thres", &tmplThres, 0.0f, 1.0f, "ratio = %.3f");
      imProc_->tmplThres = tmplThres;

      // show tmplImgs
      for (int i = 0; i < NUM_TEMPLATES; ++i)
        tmplImages_[i]->render(imProc_->impConf.getTmplImg()[i]);

      // load from/save to file
      if (ImGui::Button("Load from file"))
        imProc_->loadConfig();
      if (ImGui::Button("Save to file"))
        imProc_->saveConfig();

      ImGui::End();
    }
  }

  if (ImGui::IsKeyPressed(ImGuiKey_I))
    improcVisible_ = !improcVisible_;
  if (improcVisible_) {
    imProc_->startProcThread();
    if (ImGui::Begin("Channels", &improcVisible_)) {
      for (int i = 0; i < NUM_CHANS; ++i) {
        if (i != 0)
          ImGui::SameLine();
        chImages_[i]->render(imProc_->getProcFrame(i));
      }
      ImGui::End();
    }
    if (ImGui::Begin("Proc Image", &improcVisible_)) {
      procImage_->render(imProc_->getProcFrame());
      ImGui::End();
    }
  } else
    imProc_->stopProcThread();
}

} // namespace gui
