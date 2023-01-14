add_library(ctrl
  "./src/ctrl/supervisor.cpp"
  "./src/ctrl/state/sysidstate.cpp"
  "./src/ctrl/state/state.cpp"
  "./src/ctrl/state/state0.cpp"
  "./src/ctrl/state/state1.cpp"
  "./src/ctrl/state/state2.cpp")
target_link_libraries(ctrl PRIVATE project_options project_warnings)

# includes
set(ctrl_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(ctrl PUBLIC
  "$<BUILD_INTERFACE:${ctrl_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# link dependencies
target_link_system_libraries(ctrl PRIVATE
  util
  improc
  imcap
  cam
  pump)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
