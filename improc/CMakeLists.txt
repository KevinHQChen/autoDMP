add_library(improc
    "./src/improc/improc.cpp"
    "./src/improc/improc_helpers.cpp")

# includes
set(improc_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(improc PUBLIC
  "$<BUILD_INTERFACE:${improc_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

target_link_libraries(improc PRIVATE
  util
  imcap
  cam
  ctrl
  supervisor
  pump)

# tests
if(BUILD_TESTS)
  add_subdirectory("./test")
endif()
