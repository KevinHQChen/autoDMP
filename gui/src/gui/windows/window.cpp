#include "gui/windows/window.hpp"

namespace gui {

void Window::registerCallback(std::function<void()> callback) { callbacks_.push_back(callback); }

void Window::handleEvents() {
  for (auto &callback : callbacks_) {
    callback();
  }
}

} // namespace gui
