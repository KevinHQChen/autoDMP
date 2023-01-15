#pragma once

#include "imgui.h"
#include <functional>
#include <string>

namespace gui {

class Button {
public:
  Button(std::string label, std::function<void()> callback) : label_(label), callback_(callback) {}

  void render() {
    bool pressed = ImGui::Button(label_.c_str());
    if (pressed) {
      callback_();
    }
  }

private:
  std::string label_;
  std::function<void()> callback_;
};

} // namespace gui
