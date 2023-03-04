

add_executable(test_nlp_proximity "test_nlp_proximity.cpp")
ReaK_setup_target(test_nlp_proximity)
target_link_libraries(test_nlp_proximity reak_geom)

add_executable(test_gjk_proximity "test_gjk_proximity.cpp")
ReaK_setup_target(test_gjk_proximity)
target_link_libraries(test_gjk_proximity reak_geom)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(test_coin3d "test_coin3d.cpp")
  ReaK_setup_target(test_coin3d)
  target_link_libraries(test_coin3d ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
endif()

add_executable(build_X8_quadrotor_geom "build_X8_quadrotor_geom.cpp")
ReaK_setup_target(build_X8_quadrotor_geom)
target_link_libraries(build_X8_quadrotor_geom reak_geom)

add_executable(build_X8_obstacle_courses "build_X8_obstacle_courses.cpp")
ReaK_setup_target(build_X8_obstacle_courses)
target_link_libraries(build_X8_obstacle_courses reak_geom)

add_executable(build_MD148_lab "build_MD148_lab.cpp")
ReaK_setup_target(build_MD148_lab)
target_link_libraries(build_MD148_lab reak_geom)
target_link_libraries(build_MD148_lab absl::flags absl::flags_parse)



