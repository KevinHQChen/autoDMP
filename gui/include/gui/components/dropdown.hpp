#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class Dropdown {
  std::string label_;
  std::vector<std::string> items_;
  std::string item_, prev_item_;
  std::function<void()> callback_{nullptr};

public:
  Dropdown(std::string label, std::vector<std::string> items,
           std::function<void()> callback = nullptr)
      : label_(label), items_(items), item_(items[0]), callback_(callback) {}

  void render() {
    ImGui::Text("%s", label_.c_str());
    ImGui::BeginGroup();
    std::string label = "##" + label_;
    if (ImGui::BeginCombo(label.c_str(), item_.c_str())) {
      for (auto &item : items_) {
        if (ImGui::Selectable(item.c_str())) {
          prev_item_ = item_;
          item_ = item;
          if (item_ != prev_item_ && callback_)
            callback_();
        }
      }
      ImGui::EndCombo();
    }
    ImGui::EndGroup();
  }

  std::string getItem() { return item_; }
};

} // namespace gui
