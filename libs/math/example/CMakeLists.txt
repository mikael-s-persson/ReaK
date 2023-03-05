
add_executable(test_mat_alg "test_mat_alg.cpp")
ReaK_setup_target(test_mat_alg)
target_link_libraries(test_mat_alg reak_math)

add_executable(test_mat_num "test_mat_num.cpp")
ReaK_setup_target(test_mat_num)
target_link_libraries(test_mat_num reak_math)

add_executable(test_mat_num_perf "test_mat_num_perf.cpp")
ReaK_setup_target(test_mat_num_perf)
target_link_libraries(test_mat_num_perf reak_math absl::log)

add_executable(test_mat_are "test_mat_are.cpp")
ReaK_setup_target(test_mat_are)
target_link_libraries(test_mat_are reak_math)

add_executable(test_mat_views "test_views.cpp")
ReaK_setup_target(test_mat_views)
target_link_libraries(test_mat_views reak_math)

add_executable(test_optim "test_optim.cpp")
ReaK_setup_target(test_optim)
target_link_libraries(test_optim reak_math)

add_executable(test_optim_nllsq "test_nllsq.cpp")
ReaK_setup_target(test_optim_nllsq)
target_link_libraries(test_optim_nllsq reak_math)

add_executable(test_optim_lp "test_lp.cpp")
ReaK_setup_target(test_optim_lp)
target_link_libraries(test_optim_lp reak_math)

add_executable(test_optim_qp "test_qp.cpp")
ReaK_setup_target(test_optim_qp)
target_link_libraries(test_optim_qp reak_math)

add_executable(test_optim_nlp "test_nlp.cpp")
ReaK_setup_target(test_optim_nlp)
target_link_libraries(test_optim_nlp reak_math)

add_executable(build_airship3D_transforms "build_airship3D_transforms.cpp")
ReaK_setup_target(build_airship3D_transforms)
target_link_libraries(build_airship3D_transforms reak_math)

add_executable(test_lie_grp_num "test_lie_grp_num.cpp")
ReaK_setup_target(test_lie_grp_num)
target_link_libraries(test_lie_grp_num reak_math)


