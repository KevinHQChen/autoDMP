#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

inline void HelpMarker(const char *desc) {
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}

inline void displayVector3d(const char *vecName, Eigen::Vector3d vec) {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", vecName);
  for (int i = 0; i < 3; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%f", vec(i));
  }
}

inline void displayArray3b(const char *arrName, bool arr[3], const char *helpText = "") {
  ImGui::TableNextRow();
  ImGui::TableSetColumnIndex(0);
  ImGui::Text("%s", arrName);
  ImGui::SameLine();
  HelpMarker(helpText);
  for (int i = 0; i < 3; ++i) {
    ImGui::TableSetColumnIndex(i + 1);
    ImGui::Text("%d", arr[i]);
  }
}

struct GUIEvent {
  int srcState = 0;
  int destState = 0;
  int pos[3] = {0, 0, 0};
  int vel[3] = {0, 0, 0};
  int *data[4];
  inline static const std::string props[4] = {"Src State", "Dest State", "Target Pos (ch0-2) [%]",
                                              "Target Vel (ch0-2) [px/s]"};
  inline static const int min[4] = {0, 0, 0, 0};
  inline static const int max[4] = {3, 3, 100, 20};

  GUIEvent() {
    int i = 0;
    data[i++] = &srcState;
    data[i++] = &destState;
    data[i++] = pos;
    data[i++] = vel;
  }

  GUIEvent(int srcState, int destState, Eigen::Vector3d pos, Eigen::Vector3d vel) {
    this->srcState = srcState;
    this->destState = destState;
    for (int i = 0; i < 3; i++) {
      this->pos[i] = pos(i);
      this->vel[i] = vel(i);
    }

    int i = 0;
    data[i++] = &srcState;
    data[i++] = &destState;
    data[i++] = this->pos;
    data[i++] = this->vel;
  }
};

namespace gui {

class CtrlWindow : public Window {
  const int numChans_ = toml::get<int>(Config::conf["improc"]["numChans"]);

  bool ctrlSetupVisible_{false}, ctrlVisible_{false};

  Supervisor *sv_;

  // displaying/modifying events, states
  ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg |
                               ImGuiTableFlags_NoBordersInBody;
  GUIEvent currEvent;
  std::deque<GUIEvent> guiEventQueue;
  int openAction = -1;
  int dropletLength = 0;

public:
  CtrlWindow(Supervisor *sv);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
