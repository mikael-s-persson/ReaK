
  
add_executable(test_interp "test_interp.cpp")
ReaK_setup_target(test_interp)
target_link_libraries(test_interp reak_topologies)
target_link_libraries(test_interp absl::flags absl::flags_parse)

add_executable(test_interp_mc "test_interp_mc.cpp")
ReaK_setup_target(test_interp_mc)
target_link_libraries(test_interp_mc reak_topologies)
target_link_libraries(test_interp_mc absl::flags absl::flags_parse)

add_executable(test_interp_topo "test_interp_topo.cpp")
ReaK_setup_target(test_interp_topo)
target_link_libraries(test_interp_topo reak_topologies)
target_link_libraries(test_interp_topo absl::flags absl::flags_parse)

add_executable(test_interp_topo_rand "test_interp_topo_rand.cpp")
ReaK_setup_target(test_interp_topo_rand)
target_link_libraries(test_interp_topo_rand reak_topologies)
target_link_libraries(test_interp_topo_rand absl::flags absl::flags_parse)

add_executable(test_lim_interp "test_lim_interp.cpp")
ReaK_setup_target(test_lim_interp)
target_link_libraries(test_lim_interp reak_topologies)

add_executable(test_Ndof_interp "test_Ndof_interp.cpp")
ReaK_setup_target(test_Ndof_interp)
target_link_libraries(test_Ndof_interp reak_topologies)

if( OpenCV_FOUND )
  add_executable(test_sampling "test_sampling.cpp")
  ReaK_setup_target(test_sampling)
  target_link_libraries(test_sampling reak_topologies)
endif()

