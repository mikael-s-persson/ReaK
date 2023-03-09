
find_package( Qt4 )
if(COIN3D_FOUND AND QT_FOUND)
  find_library(SOQT4_LIBRARIES SoQt4)
endif()


add_subdirectory(systems)
add_subdirectory(integrators)
add_subdirectory(estimators)
add_subdirectory(controllers)
add_subdirectory(qt)

reak_cc_library(
  NAME
    control
  DEPS
    ReaK::controllers
    ReaK::estimators
    ReaK::integrators
    ReaK::systems
  PUBLIC
)

message("Configured: ReaK Control and Estimation Library Components")

