# project(ReaK_topologies)


find_package(OpenCV)

if(OpenCV_FOUND)
  message(STATUS "OpenCV library components were found: '${OpenCV_LIBS}'")
  add_definitions( "-DREAK_HAS_OPENCV" )
else()
  message(WARNING "OpenCV was not detected on this system, some test programs will not be built!")
endif()

add_subdirectory(interpolation)
add_subdirectory(spaces)

reak_cc_library(
  NAME
    topologies
  DEPS
    ReaK::spaces
    ReaK::interpolation
  PUBLIC
)

message("Configured: ReaK Topologies Library Components")
