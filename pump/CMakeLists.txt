add_library(pump "./src/pump/pump.cpp")
target_link_libraries(pump PRIVATE project_options project_warnings)

# includes
set(pump_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(pump PUBLIC
    "$<BUILD_INTERFACE:${pump_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>")
    # PUBLIC ../../external/fgt_SDK_Cpp           # for fgt_SDK_Cpp.h
    # PUBLIC ../../external/fgt_SDK_Cpp/dlls)     # for fgt_SDK.h

target_link_system_libraries(pump PRIVATE
    util)
    # fgt_SDK_cpp)

# tests
if(FEATURE_TESTS)
  add_subdirectory("./test")
endif()
