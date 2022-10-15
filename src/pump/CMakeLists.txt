# Pump library for RoboDrop2
set(HEADER_LIST "../../include/pump/pump.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(pump pump.cpp ${HEADER_LIST})

# This depends on fluigent library
# link each library where the target is created
target_link_libraries(pump
    PUBLIC util
    # PRIVATE fgt_SDK_cpp)
)

# We need this directory, and users of our library will need it too
target_include_directories(pump
    PUBLIC ../../include
    # PUBLIC ../../external/fgt_SDK_Cpp           # for fgt_SDK_Cpp.h
    # PUBLIC ../../external/fgt_SDK_Cpp/dlls)     # for fgt_SDK.h
)
