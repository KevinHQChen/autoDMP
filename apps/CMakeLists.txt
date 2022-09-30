# add_executable(autoDMP autoDMP.cpp)
# target_link_libraries(robodrop PRIVATE util cam pump ctrl improc)

# add_executable(intro main.cpp)
# target_link_libraries(intro PRIVATE util cam gui imcap improc)
add_executable(testapp test.cpp)
target_link_libraries(testapp PRIVATE util)

target_include_directories(util
  PRIVATE "${CMAKE_BINARY_DIR}/configured_files/include"
)

# add_executable(sml_test sml_test.cpp)
# target_link_libraries(sml_test PRIVATE util)
