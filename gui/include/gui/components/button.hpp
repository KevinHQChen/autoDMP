#pragma once

#include "imgui.h"
#include "imgui_toggle/imgui_toggle.h"
#include "imgui_toggle/imgui_toggle_palette.h"
#include "imgui_toggle/imgui_toggle_presets.h"
#include "imgui_toggle/imgui_toggle_renderer.h"
#include <cassert>
#include <functional>
#include <string>

namespace gui {

class Toggle {
  std::string label_;
  bool *value_{nullptr};
  bool changed_{false};
  ImVec2 toggleSize_;
  std::function<void()> callback_;

public:
  Toggle(std::string label, bool *value = nullptr, std::function<void()> callback = nullptr,
         ImVec2 toggleSize = ImVec2(0, 0))
      : label_(label), value_(value), callback_(callback), toggleSize_(toggleSize) {
    if (value == nullptr)
      value_ = new bool(false);
    assert(value_ != nullptr);
  }

  void render() {
    ImGui::BeginGroup();
    changed_ = ImGui::Toggle(label_.c_str(), value_, toggleSize_);
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  bool changed() const { return changed_; }

  void setLabel(std::string label) { label_ = label; }

  bool get() const { return *value_; }
};

class ToggleArray {
  std::string label_;
  bool *values_{nullptr};
  std::vector<bool> changed_;
  int numToggles_;
  ImVec2 toggleSize_;
  std::function<void()> callback_;

public:
  ToggleArray(std::string label, bool *values, int numToggles,
              std::function<void()> callback = nullptr, ImVec2 toggleSize = ImVec2(0, 0))
      : label_(label), values_(values), changed_(std::vector<bool>(numToggles, false)),
        numToggles_(numToggles), callback_(callback), toggleSize_(toggleSize) {
    if (values_ == nullptr)
      values_ = new bool[numToggles_]{false};
  }

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < numToggles_; i++) {
      std::string label = "##" + label_ + std::to_string(i);
      changed_[i] = ImGui::Toggle(label.c_str(), &values_[i], toggleSize_);
      ImGui::SameLine();
    }
    ImGui::Text("%s", label_.c_str());
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  bool changed(int idx) const { return changed_[idx]; }

  void setLabel(std::string label) { label_ = label; }

  std::string getLabel() const { return label_; }

  bool get(int idx) const { return values_[idx]; }

  bool *get() const { return values_; }
};

class Button {
  std::string label_;
  bool *value_{nullptr};
  ImVec2 buttonSize_;
  std::function<void()> callback_;

public:
  Button(std::string label, bool *value = nullptr, std::function<void()> callback = nullptr,
         ImVec2 buttonSize = ImVec2(0, 0))
      : label_(label), value_(value), callback_(callback), buttonSize_(buttonSize) {
    if (value == nullptr)
      value_ = new bool(false);
    assert(value_ != nullptr);
  }

  void render() {
    ImGui::BeginGroup();
    *value_ = ImGui::Button(label_.c_str(), buttonSize_);
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  void setLabel(std::string label) { label_ = label; }

  bool get() const { return *value_; }
};

class ButtonArray {
  std::vector<std::string> labels_;
  std::vector<bool> *values_{nullptr};
  int numButtons_;
  ImVec2 buttonSize_;
  std::function<void()> callback_;

public:
  ButtonArray(std::vector<std::string> labels, std::vector<bool> *values,
              std::function<void()> callback = nullptr, ImVec2 buttonSize = ImVec2(0, 0))
      : labels_(labels), values_(values), callback_(callback), buttonSize_(buttonSize) {
    numButtons_ = labels_.size();
    if (values_ == nullptr)
      values_ = new std::vector<bool>(numButtons_, false);
  }

  void render() {
    ImGui::BeginGroup();
    for (int i = 0; i < numButtons_; i++) {
      (*values_)[i] = ImGui::Button(labels_[i].c_str(), buttonSize_);
      ImGui::SameLine();
    }
    ImGui::EndGroup();
    if (callback_)
      callback_();
  }

  void setLabel(std::string label, int idx = 0) { labels_[idx] = label; }

  std::string getLabel(int idx = 0) const { return labels_[idx]; }

  bool get(int index) const { return (*values_)[index]; }

  std::vector<bool> get() const { return *values_; }
};

} // namespace gui
