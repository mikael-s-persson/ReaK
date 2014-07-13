
set(SSSYSTEMS_SOURCES 
  "${SRCROOT}${RKSSSYSTEMSDIR}/quadrotor_system.cpp"
  "${SRCROOT}${RKSSSYSTEMSDIR}/satellite_basic_models.cpp"
  "${SRCROOT}${RKSSSYSTEMSDIR}/satellite_invar_models.cpp"
  "${SRCROOT}${RKSSSYSTEMSDIR}/satellite_modeling_options.cpp"
  "${SRCROOT}${RKSSSYSTEMSDIR}/near_buoyant_airship_models.cpp"
  "${SRCROOT}${RKSSSYSTEMSDIR}/airship_assembled_models.cpp"
)


set(SSSYSTEMS_HEADERS 
  "${RKSSSYSTEMSDIR}/airship_assembled_models.hpp"
  "${RKSSSYSTEMSDIR}/airship_basic_mixins.hpp"
  "${RKSSSYSTEMSDIR}/airship_drag_mixins.hpp"
  "${RKSSSYSTEMSDIR}/airship_IMU_mixins.hpp"
  "${RKSSSYSTEMSDIR}/airship_sonar_mixins.hpp"
  "${RKSSSYSTEMSDIR}/airship_thruster_mixins.hpp"
  "${RKSSSYSTEMSDIR}/quadrotor_system.hpp"
  "${RKSSSYSTEMSDIR}/satellite_basic_models.hpp"
  "${RKSSSYSTEMSDIR}/satellite_invar_models.hpp"
  "${RKSSSYSTEMSDIR}/satellite_modeling_options.hpp"
  "${RKSSSYSTEMSDIR}/satellite_modeling_po.hpp"
  "${RKSSSYSTEMSDIR}/near_buoyant_airship_models.hpp"
)

add_library(reak_ss_systems STATIC ${SSSYSTEMS_SOURCES})
setup_custom_target(reak_ss_systems)
setup_headers("${SSSYSTEMS_HEADERS}" "${RKSSSYSTEMSDIR}")
target_link_libraries(reak_ss_systems reak_topologies reak_mbd reak_core)

if(COIN3D_FOUND AND QT_FOUND)
  
#   include( ${QT_USE_FILE} )
  
#   add_executable(X8_test_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/X8_test_scene.cpp")
#   setup_custom_target(X8_test_scene)
#   target_link_libraries(X8_test_scene reak_examples reak_topologies reak_mbd_coin reak_core)
#   target_link_libraries(X8_test_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()
