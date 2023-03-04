
add_executable(build_quadrotor_space "build_quadrotor_space.cpp")
ReaK_setup_target(build_quadrotor_space)
target_link_libraries(build_quadrotor_space reak_control reak_topologies reak_mbd reak_math reak_core)

add_executable(quadrotor_scene "quadrotor_scene.cpp")
ReaK_setup_target(quadrotor_scene)
target_link_libraries(quadrotor_scene reak_control reak_topologies reak_mbd reak_math reak_core)

add_executable(estimate_airship "estimate_airship.cpp")
ReaK_setup_target(estimate_airship)
target_link_libraries(estimate_airship reak_control ${Boost_LIBRARIES})

add_executable(estimate_satellite3D "estimate_satellite3D.cpp")
ReaK_setup_target(estimate_satellite3D)
target_link_libraries(estimate_satellite3D reak_control)
target_link_libraries(estimate_satellite3D absl::flags absl::flags_parse)

add_executable(predict_satellite3D "predict_satellite3D.cpp")
ReaK_setup_target(predict_satellite3D)
target_link_libraries(predict_satellite3D reak_control)
target_link_libraries(predict_satellite3D absl::flags absl::flags_parse)

# add_executable(test_am "test_am.cpp")
# ReaK_setup_target(test_am)
# target_link_libraries(test_am reak_mbd)
# 
# add_executable(test_bm "test_bm.cpp")
# ReaK_setup_target(test_bm)
# target_link_libraries(test_bm reak_mbd)

