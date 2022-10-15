# Supervisory control library for autoDMP
set(HEADER_LIST "../../include/ctrl/supervisor.hpp"
    "../../include/ctrl/state/state.hpp"
    "../../include/ctrl/state/state0.hpp"
    "../../include/ctrl/state/state1.hpp"
    "../../include/ctrl/state/sysidstate.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(ctrl supervisor.cpp state/sysidstate.cpp
  state/state.cpp state/state0.cpp state/state1.cpp ${HEADER_LIST})   # target is created

target_link_libraries(
  ctrl
  PUBLIC util
  # PUBLIC improc # why don't we need this?
)

# We need this directory, and users of our library will need it too
target_include_directories(ctrl
  PUBLIC ../../include
)
