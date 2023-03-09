# project(ReaK_planning)

find_package(Coin3D)
if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
  include_directories(BEFORE SYSTEM ${COIN3D_INCLUDE_DIRS})
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

find_package( Qt4 )
if(COIN3D_FOUND AND QT_FOUND)
  find_library(SOQT4_LIBRARIES SoQt4)
endif()


add_subdirectory(graph_alg)
add_subdirectory(path_planning)
add_subdirectory(qt)

reak_cc_library(
  NAME
    planning
  DEPS
    ReaK::path_planning
    ReaK::graph_alg
  PUBLIC
)

message("Configured: ReaK Motion-Planning Library Components")


