
set(MISCEXAMPLES_SOURCES 
#   "${SRCROOT}${RKMISCEXAMPLESDIR}/IHAQR_topology.cpp"
)


set(MISCEXAMPLES_HEADERS 
  "${RKMISCEXAMPLESDIR}/quadrotor_system.hpp"
  "${RKMISCEXAMPLESDIR}/IHAQR_topology.hpp"
  "${RKMISCEXAMPLESDIR}/MEAQR_topology.hpp"
)

# add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
# setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
# setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")
# target_link_libraries(reak_robot_airship reak_topologies reak_interp reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/robot_kin_scene.cpp")
  setup_custom_target(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(robot_kin_scene reak_kte_coin reak_geom_coin reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
  target_link_libraries(robot_kin_scene ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(quadrotor_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/quadrotor_scene.cpp")
  setup_custom_target(quadrotor_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(quadrotor_scene reak_kte_coin reak_geom_coin reak_topologies reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
  target_link_libraries(quadrotor_scene ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
endif()

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")




