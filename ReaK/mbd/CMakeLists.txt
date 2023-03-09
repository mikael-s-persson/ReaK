# project(ReaK_mbd)

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

add_subdirectory(kte)
add_subdirectory(models)
add_subdirectory(coin3D)
add_subdirectory(qt)

reak_cc_library(
  NAME
    mbd
  DEPS
    ReaK::kte
    ReaK::models
  PUBLIC
)

message("Configured: ReaK Multi-body Dynamics Library Components")

