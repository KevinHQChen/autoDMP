// config.hpp will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "gui/gui.hpp"

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int, const char **) {
  GUI gui;
  gui.startGUIThread();
}
