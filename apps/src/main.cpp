#include "gui/gui.hpp"

namespace py = pybind11;
using namespace py::literals;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int, const char **) {
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
  info("Config type: {}", type_name<decltype(Config::guiConf)>());
  info("Parsed config: {}", toml::find(Config::conf, "gui"));
  GUI gui(imCap, imProc, pump, sv);
  gui.startGUIThread();

  delete sv;
  delete pump;
  delete imProc;
  delete imCap;
}
