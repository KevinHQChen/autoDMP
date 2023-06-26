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
  auto cam = new Cam(0, Config::conf, console.getLogger("Cam"));
  auto imCap = new ImCap(cam, console.getLogger("ImCap"));
  auto imProc = new ImProc(imCap, console.getLogger("ImProc"));
  auto pump = new Pump(console.getLogger("Pump"));
  auto sv = new Supervisor(imProc, pump, console.getLogger("Supervisor"));

  // start gui
  GUI gui(imCap, imProc, pump, sv, console, console.getLogger("GUI"));
  gui.startThread();

  delete sv;
  delete pump;
  delete imProc;
  delete imCap;
}
