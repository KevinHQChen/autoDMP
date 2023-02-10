add_library(cam "./src/cam/cam.cpp")

# This depends on andor, opencv libraries
# https://www.pragmaticlinux.com/2022/03/how-to-link-a-shared-library-with-gcc-and-cmake/
# target_link_libraries(cam PRIVATE project_options project_warnings)

# includes
set(cam_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(cam PUBLIC
  "$<BUILD_INTERFACE:${cam_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# link dependencies
target_link_system_libraries(cam PRIVATE
  util
  atcore
  atutility)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
