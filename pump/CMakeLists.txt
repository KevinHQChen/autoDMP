add_library(pump "./src/pump/pump.cpp")

# includes
include(FetchContent)
FetchContent_Declare(
    fgt_SDK_Cpp
    GIT_REPOSITORY https://github.com/Fluigent/fgt-SDK.git
    GIT_PROGRESS TRUE
    # Enter the desired git tag below (version 22.2)
    GIT_TAG "e158ea873f0ed7ac7bbbcc607496cfe19f8448aa"
    SOURCE_SUBDIR "C++"
)
FetchContent_MakeAvailable(fgt_SDK_Cpp)

set(pump_INCLUDE_DIR
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  CACHE STRING "")
target_include_directories(pump PUBLIC
    "$<BUILD_INTERFACE:${pump_INCLUDE_DIR}>"
    "$<INSTALL_INTERFACE:./${CMAKE_INSTALL_INCLUDEDIR}>"
    "${CMAKE_BINARY_DIR}/_deps/fgt_sdk_cpp-src/C++/fgt_SDK_Cpp"
    "${CMAKE_BINARY_DIR}/_deps/fgt_sdk_cpp-src/C++/fgt_SDK_Cpp/dlls")

# link dependencies
target_link_libraries(pump PRIVATE
    util
    fgt_SDK_cpp)

# tests
if(BUILD_TESTS)
  add_subdirectory("./test")
endif()
