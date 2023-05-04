#pragma once

#include "gui/components/button.hpp"
#include "gui/components/slider.hpp"
#include "window.hpp"

#include "ctrl/supervisor.hpp"

#define NUM_CHANS 3

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

namespace gui {

class CtrlWindow : public Window {
  const int numChans_ = toml::get<int>(Config::conf["improc"]["numChans"]);

  bool ctrlSetupVisible_{false}, ctrlVisible_{false};

  Supervisor *sv_;

  // displaying/modifying events, states
  ImGuiTableFlags tableFlags = ImGuiTableFlags_BordersV | ImGuiTableFlags_BordersOuterH |
                               ImGuiTableFlags_Resizable | ImGuiTableFlags_RowBg |
                               ImGuiTableFlags_NoBordersInBody;

  int srcState, destState, moveTime, holdTime;
  int targetPos[NUM_CHANS];
  bool chs0[3] = {true, false, false}, chs1[3] = {false, true, true}, chs2[3] = {true, false, true};

  int openAction = -1;
  int dropletLength = 0;

  event_bus getEvent(int srcState, int destState, std::array<int, NUM_CHANS> targetPos,
                     int moveTime, int holdTime) {
    event_bus e;
    e.srcState = srcState;
    e.destState = destState;
    for (int i = 0; i < NUM_CHANS; ++i)
      e.destPos[i] = targetPos[i] / 100.0;
    e.moveTime = moveTime;
    e.holdTime = holdTime;

    switch (srcState) {
    case 0: // [1 0 0]
      std::memcpy(e.chs, chs0, sizeof(chs0));
      break;
    case 1: // [0 1 1]
      std::memcpy(e.chs, chs1, sizeof(chs1));
      break;
    case 2: // [1 0 1]
      std::memcpy(e.chs, chs2, sizeof(chs2));
      break;
    }
    switch (destState) {
    case 0: // [1 0 0]
      std::memcpy(e.chs, chs0, sizeof(chs0));
      break;
    case 1: // [0 1 1]
      std::memcpy(e.chs, chs1, sizeof(chs1));
      break;
    case 2: // [1 0 1]
      std::memcpy(e.chs, chs2, sizeof(chs2));
      break;
    }
    return e;
  }

public:
  CtrlWindow(Supervisor *sv);
  ~CtrlWindow();
  void render() override;
};

} // namespace gui
