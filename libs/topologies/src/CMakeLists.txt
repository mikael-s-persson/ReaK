
add_library(reak_topologies STATIC 
  "hyperball_topology_ext1.cpp"
  "hyperbox_topology_ext1.cpp"
  "hyperbox_topology_ext2.cpp"
  "hyperbox_topology_ext3.cpp"
  "joint_space_topologies_ext1.cpp"
  "joint_space_topologies_ext2.cpp"
  "joint_space_topologies_ext3.cpp"
  "joint_space_topologies_ext4.cpp"
  "joint_space_topologies_ext5.cpp"
  "line_topology_ext1.cpp"
  "manip_free_dynamic_workspace_ext1.cpp"
  "manip_free_workspace_ext1.cpp"
  "Ndof_cubic_spaces_ext.cpp"
  "Ndof_limits_ext1.cpp"
  "Ndof_linear_spaces_ext.cpp"
  "Ndof_quintic_spaces_ext.cpp"
  "Ndof_sap_spaces_ext.cpp"
  "Ndof_spaces_ext1.cpp"
  "Ndof_spaces_ext2.cpp"
  "Ndof_spaces_ext3.cpp"
  "Ndof_svp_spaces_ext.cpp"
  "ptrobot2D_test_world.cpp"
#   "sap_Ndof_old_detail.cpp"
  "se2_topologies_ext1.cpp"
  "se3_topologies_ext1.cpp"
  "so3_topologies_ext1.cpp"
  "sustained_acceleration_pulse_detail.cpp"
  "sustained_acceleration_pulse_Ndof_detail.cpp"
  "sustained_velocity_pulse_detail.cpp"
  "sustained_velocity_pulse_Ndof_detail.cpp"
  "vector_topology_ext1.cpp"
)
ReaK_setup_static_library(reak_topologies)
target_link_libraries(reak_topologies reak_mbd reak_geom reak_math)
if( OpenCV_FOUND )
  target_link_libraries(reak_topologies ${OpenCV_LIBS})
endif()

