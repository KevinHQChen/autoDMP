add_executable(sinvec "/autoDMP/pump_pressure_inputs/test/class_test.cpp" "/autoDMP/pump_pressure_inputs/src/pump_pressure_inputs.cpp")
target_include_directories(sinvec PRIVATE "/autoDMP/pump_pressure_inputs/include")
target_link_libraries(sinvec PRIVATE util)