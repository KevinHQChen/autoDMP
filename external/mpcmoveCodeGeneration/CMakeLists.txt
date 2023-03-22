file(GLOB_RECURSE SOURCES "*.c" "*.h")

add_library(mpcmove ${SOURCES})

target_include_directories(mpcmove PUBLIC /autoDMP/external/MATLAB/extern/include)

set(mpcmove_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}"
  CACHE STRING "")
target_include_directories(mpcmove PUBLIC
  "$<BUILD_INTERFACE:${mpcmove_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")
