
set(MISCEXAMPLES_SOURCES 
  "${SRCROOT}${RKMISCEXAMPLESDIR}/X8_quadrotor_model.cpp"
  "${SRCROOT}${RKMISCEXAMPLESDIR}/X8_quadrotor_geom.cpp"
)


set(MISCEXAMPLES_HEADERS 
  "${RKMISCEXAMPLESDIR}/quadrotor_system.hpp"
  "${RKMISCEXAMPLESDIR}/IHAQR_topology.hpp"
  "${RKMISCEXAMPLESDIR}/MEAQR_topology.hpp"
  "${RKMISCEXAMPLESDIR}/X8_quadrotor_model.hpp"
  "${RKMISCEXAMPLESDIR}/X8_quadrotor_geom.hpp"
)

add_library(reak_examples STATIC ${MISCEXAMPLES_SOURCES})
setup_custom_target(reak_examples "${SRCROOT}${RKMISCEXAMPLESDIR}")
setup_headers("${MISCEXAMPLES_HEADERS}" "${RKMISCEXAMPLESDIR}")
target_link_libraries(reak_examples reak_topologies reak_interp reak_mbd reak_core)
  
add_executable(build_X8_obstacle_courses "${SRCROOT}${RKMISCEXAMPLESDIR}/build_X8_obstacle_courses.cpp")
setup_custom_target(build_X8_obstacle_courses "${SRCROOT}${RKMISCEXAMPLESDIR}")
target_link_libraries(build_X8_obstacle_courses reak_mbd reak_core)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/robot_kin_scene.cpp")
  setup_custom_target(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(robot_kin_scene reak_mbd_coin reak_core)
  target_link_libraries(robot_kin_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(quadrotor_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/quadrotor_scene.cpp")
  setup_custom_target(quadrotor_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(quadrotor_scene reak_examples reak_topologies reak_mbd_coin reak_core)
  target_link_libraries(quadrotor_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(build_quadrotor_space "${SRCROOT}${RKMISCEXAMPLESDIR}/build_quadrotor_space.cpp")
  setup_custom_target(build_quadrotor_space "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(build_quadrotor_space reak_examples reak_topologies reak_mbd_coin reak_core)
  
  add_executable(X8_test_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/X8_test_scene.cpp")
  setup_custom_target(X8_test_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(X8_test_scene reak_examples reak_topologies reak_mbd_coin reak_core)
  target_link_libraries(X8_test_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()

add_executable(X8_run_planner "${SRCROOT}${RKMISCEXAMPLESDIR}/X8_run_planner.cpp")
setup_custom_target(X8_run_planner "${SRCROOT}${RKMISCEXAMPLESDIR}")
target_link_libraries(X8_run_planner reak_examples reak_topologies reak_mbd reak_core)


