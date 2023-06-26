#include "gui/gui.hpp"

namespace py = pybind11;
using namespace py::literals;

int main(int, const char **) {
  // initialize logging
  ImGuiConsole console;

  // initialize python interpreter
  py::scoped_interpreter python;
  py::eval_file("ctrl/scripts/sysid.py"); // import sysid functions
  py::gil_scoped_release release;         // add this to release the GIL

  // initialize all subsystems
  auto imCap = new ImCap();
  auto imProc = new ImProc(imCap);
  auto pump = new Pump();
  auto sv = new Supervisor(imProc, pump);

  // start gui
  GUI gui(imCap, imProc, pump, sv, console);
  gui.startThread();

  delete sv;
  delete pump;
  delete imProc;
  delete imCap;
}
