// config.hpp will be generated automatically when you run the CMake
// configuration step. It creates a namespace called `autoDMP`. You can modify
// the source template at `configured_files/config.hpp.in`. #include
// <internal_use_only/config.hpp>

#include "gui/gui.hpp"

namespace py = pybind11;
using namespace py::literals;

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int, const char **) {
  py::scoped_interpreter python;
  py::eval_file("ctrl/scripts/sysid.py"); // import sysid functions
  py::gil_scoped_release release;         // add this to release the GIL

  ordered_value conf = TOML11_PARSE_IN_ORDER("config/setup.toml");
  guiConfig guiConf = toml::find<guiConfig>(conf, "gui");

  info("Config type: {}", type_name<decltype(guiConf)>());
  info("Parsed config: {}", toml::find(conf, "gui"));

  auto imCap = new ImCap();
  auto imProc = new ImProc(imCap);
  auto pump = new Pump(toml::get<bool>(conf["ctrl"]["simMode"]));
  auto sv = new Supervisor(imProc, pump);

  GUI gui(imCap, imProc, pump, sv);
  gui.startGUIThread();
}
