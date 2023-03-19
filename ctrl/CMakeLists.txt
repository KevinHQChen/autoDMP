add_library(ctrl
  "./src/ctrl/supervisor.cpp"
  "./src/ctrl/state/sysidstate.cpp"
  "./src/ctrl/state/state.cpp"
  "./src/ctrl/state/state0.cpp"
  "./src/ctrl/state/state1.cpp"
  "./src/ctrl/state/state2.cpp")

# includes
set(ctrl_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(ctrl PUBLIC
  "$<BUILD_INTERFACE:${ctrl_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")

# find dependencies
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
# print python version
message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python include dir: ${Python_INCLUDE_DIRS}")
message(STATUS "Python libraries: ${Python_LIBRARIES}")
message(STATUS "Python interpreter: ${Python_EXECUTABLE}")
message(STATUS "Python site packages: ${Python_SITELIB}")
# if a python virtual environment does not exist, create it, update pip, and use it to install dependencies
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/venv")
  execute_process(COMMAND ${Python_EXECUTABLE} -m venv venv WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
  # update pip and install dependencies
  execute_process(COMMAND venv/bin/python -m pip install --upgrade pip WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
  execute_process(COMMAND venv/bin/python -m pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/scripts/requirements.txt WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
endif()
# use virtual environment python interpreter
set(Python_EXECUTABLE ${CMAKE_CURRENT_SOURCE_DIR}/scripts/venv/bin/python)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/scripts/venv/bin/python)

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
target_link_libraries(ctrl PRIVATE
  ${Python_LIBRARY_RELEASE}
  pybind11::embed
  util
  improc
  imcap
  cam
  pump)

# tests
if(BUILD_TESTS)
  add_subdirectory("./test")
endif()
