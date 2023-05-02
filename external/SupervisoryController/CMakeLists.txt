file(GLOB_RECURSE SOURCES "*.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/scripts/SupervisoryController_ert_rtw/ert_main.cpp")

add_library(supervisor ${SOURCES})

set(supervisor_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}"
  CACHE STRING "")

# foreach(src ${SOURCES})
#     get_filename_component(dir ${src} DIRECTORY)
#     list(APPEND supervisor_INCLUDE_DIR ${dir})
# endforeach()
# list(REMOVE_DUPLICATES supervisor_INCLUDE_DIR)

target_include_directories(supervisor PUBLIC
  "$<BUILD_INTERFACE:${supervisor_INCLUDE_DIR}>"
  "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>"
  "${supervisor_INCLUDE_DIR}/scripts/SupervisoryController_ert_rtw"
  "${supervisor_INCLUDE_DIR}/scripts/slprj/ert/_sharedutils"
  "${supervisor_INCLUDE_DIR}/R2022b/simulink/include")
