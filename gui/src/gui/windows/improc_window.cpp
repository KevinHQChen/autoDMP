#include "gui/windows/improc_window.hpp"

namespace gui {

ImProcWindow::ImProcWindow(std::shared_ptr<ImCap> imCap, std::shared_ptr<ImProc> imProc)
    : imCap_(imCap), imProc_(imProc) {
  immvisionParams = ImmVision::ImageParams();
  immvisionParams.ImageDisplaySize = cv::Size(350, 0);
  immvisionParams.ZoomKey = "z";
  immvisionParams.RefreshImage = true;
}

ImProcWindow::~ImProcWindow() {}

void ImProcWindow::render() {
  if (ImGui::IsKeyPressed(ImGuiKey_C))
    imcapVisible_ = !imcapVisible_;
  if (imcapVisible_) {
    imCap_->startCaptureThread();

    if (ImGui::Begin("Raw Image Capture", &imcapVisible_)) {
      tmpMat = imCap_->getRawFrame();
      if (!tmpMat.empty())
        rawMat = tmpMat.clone();
      ImmVision::Image("##raw", rawMat, &immvisionParams);
      ImGui::End();
    }
  } else
    imCap_->stopCaptureThread();
}




} // namespace gui
