#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class Dropdown {
public:
  Dropdown(std::string label, std::vector<std::string> items, std::function<void()> callback)
      : label_(label), items_(items), item_(items[0]), callback_(callback) {}

  void render() {
    ImGui::Text("%s", label_.c_str());
    ImGui::BeginGroup();
    std::string label = "##" + label_;
    if (ImGui::BeginCombo(label.c_str(), item_.c_str())) {
      for (auto &item : items_) {
        if (ImGui::Selectable(item.c_str())) {
          item_ = item;
        }
      }
      ImGui::EndCombo();
    }
    ImGui::EndGroup();
  }

private:
  std::string label_;
  std::vector<std::string> items_;
  std::string item_;
  std::function<void()> callback_;
};

} // namespace gui
