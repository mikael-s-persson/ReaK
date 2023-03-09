# project(ReaK_geom)

find_package(Coin3D)
if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

add_subdirectory(shapes)
add_subdirectory(proximity)

reak_cc_library(
  NAME
    geom
  DEPS
    ReaK::shapes
    ReaK::proximity
  PUBLIC
)

message("Configured: ReaK Geometry Library Components")

