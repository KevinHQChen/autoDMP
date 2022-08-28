# Camera library for RoboDrop2
set(HEADER_LIST "../../include/cam/cam.hpp")

# Make an automatic library - will be static or dynamic based on user setting
add_library(cam cam.cpp ${HEADER_LIST})     # target is created

# This depends on andor, opencv libraries
# https://www.pragmaticlinux.com/2022/03/how-to-link-a-shared-library-with-gcc-and-cmake/
target_link_libraries(cam
    PUBLIC util
    PRIVATE atcore
    PRIVATE atutility)    # link each library where the target is created

# We need this directory, and users of our library will need it too
target_include_directories(cam
    PUBLIC ../../include)
