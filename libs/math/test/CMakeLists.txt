

add_executable(unit_test_vect_alg "unit_test_vect_alg.cpp")
ReaK_setup_test_program(unit_test_vect_alg)
target_link_libraries(unit_test_vect_alg reak_math GTest::gtest_main)

add_executable(unit_test_mat_alg "unit_test_mat_alg.cpp")
ReaK_setup_test_program(unit_test_mat_alg)
target_link_libraries(unit_test_mat_alg reak_math GTest::gtest_main)

add_executable(unit_test_mat_num "unit_test_mat_num.cpp")
ReaK_setup_test_program(unit_test_mat_num)
target_link_libraries(unit_test_mat_num reak_math GTest::gtest_main)

add_executable(unit_test_mat_are "unit_test_mat_are.cpp")
ReaK_setup_test_program(unit_test_mat_are)
add_custom_command(TARGET unit_test_mat_are PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  make_directory "${CMAKE_CURRENT_BINARY_DIR}/are_data")
add_custom_command(TARGET unit_test_mat_are PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  copy "${CMAKE_CURRENT_SOURCE_DIR}/carex_m/carex_data.txt" "${CMAKE_CURRENT_BINARY_DIR}/are_data/")
add_custom_command(TARGET unit_test_mat_are PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  copy "${CMAKE_CURRENT_SOURCE_DIR}/darex_m/darex_data.txt" "${CMAKE_CURRENT_BINARY_DIR}/are_data/")
target_link_libraries(unit_test_mat_are reak_math GTest::gtest_main)

add_executable(unit_test_rotations "unit_test_rotations.cpp")
ReaK_setup_test_program(unit_test_rotations)
target_link_libraries(unit_test_rotations reak_math GTest::gtest_main)

add_executable(unit_test_quat_alg "unit_test_quat_alg.cpp")
ReaK_setup_test_program(unit_test_quat_alg)
target_link_libraries(unit_test_quat_alg reak_math GTest::gtest_main)

add_executable(unit_test_root_finders "unit_test_root_finders.cpp")
ReaK_setup_test_program(unit_test_root_finders)
target_link_libraries(unit_test_root_finders reak_math GTest::gtest_main)




