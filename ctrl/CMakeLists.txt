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

# find dependencies
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
set(ctrl_DEPENDENCIES_CONFIGURED
  pybind11
  CACHE STRING "")
foreach (dep ${ctrl_DEPENDENCIES_CONFIGURED})
  find_package(${dep} REQUIRED) # certain packages like boost may need different find_package calls, check the output of run_vcpkg
  if (DEFINED ${dep}_INCLUDE_DIRS)
    target_include_directories(ctrl PUBLIC ${${dep}_INCLUDE_DIRS})
  endif ()
endforeach ()


# link dependencies
target_link_system_libraries(ctrl PRIVATE
  pybind11::embed
  util
  improc
  imcap
  cam
  pump)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
