
set(TOPOLOGIES_SOURCES 
  "${SRCROOT}${RKTOPOLOGIESDIR}/ptrobot2D_test_world.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/hyperball_topology_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/hyperbox_topology_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/hyperbox_topology_ext2.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/hyperbox_topology_ext3.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_limits_ext1.cpp"
#   "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_topologies_ext1.cpp"
#   "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_topologies_ext2.cpp"
#   "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_topologies_ext3.cpp"
#   "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_topologies_ext4.cpp"
#   "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/joint_space_topologies_ext5.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/line_topology_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/manip_free_workspace_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/manip_free_dynamic_workspace_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_limits_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_spaces_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_spaces_ext2.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_spaces_ext3.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_linear_spaces_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_cubic_spaces_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_quintic_spaces_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_svp_spaces_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/Ndof_sap_spaces_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/se2_topologies_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/se3_topologies_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/so3_topologies_ext1.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/extern_templates/vector_topology_ext1.cpp"
)

set(TOPOLOGIES_HEADERS 
  "${RKTOPOLOGIESDIR}/default_random_sampler.hpp"
  "${RKTOPOLOGIESDIR}/differentiable_space.hpp"
  "${RKTOPOLOGIESDIR}/direct_kinematics_topomap.hpp"
  "${RKTOPOLOGIESDIR}/direct_kinematics_topomap_detail.hpp"
  "${RKTOPOLOGIESDIR}/fixed_topology_random_sampler.hpp"
  "${RKTOPOLOGIESDIR}/generic_sampler_factory.hpp"
  "${RKTOPOLOGIESDIR}/hyperball_topology.hpp"
  "${RKTOPOLOGIESDIR}/hyperbox_topology.hpp"
  "${RKTOPOLOGIESDIR}/inverse_kinematics_topomap.hpp"
  "${RKTOPOLOGIESDIR}/inverse_kinematics_topomap_detail.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits.tpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits_detail.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_topologies.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_topologies_ext.hpp"
  "${RKTOPOLOGIESDIR}/line_topology.hpp"
  "${RKTOPOLOGIESDIR}/manip_3R3R_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_3R_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_ERA_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_dynamic_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_dynamic_workspace_ext.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_workspace_ext.hpp"
  "${RKTOPOLOGIESDIR}/manip_P3R3R_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_planning_traits.hpp"
  "${RKTOPOLOGIESDIR}/manip_SCARA_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_SSRMS_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/manip_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/metric_space_tuple.hpp"
  "${RKTOPOLOGIESDIR}/metric_space_tuple_fwd.hpp"
  "${RKTOPOLOGIESDIR}/Ndof_limits.hpp"
  "${RKTOPOLOGIESDIR}/Ndof_spaces.hpp"
  "${RKTOPOLOGIESDIR}/no_obstacle_space.hpp"
  "${RKTOPOLOGIESDIR}/proxy_model_updater.hpp"
  "${RKTOPOLOGIESDIR}/proxy_traj_applicator.hpp"
  "${RKTOPOLOGIESDIR}/ptrobot2D_test_world.hpp"
  "${RKTOPOLOGIESDIR}/rate_limited_space_metamaps.hpp"
  "${RKTOPOLOGIESDIR}/rate_limited_spaces.hpp"
  "${RKTOPOLOGIESDIR}/reachability_space.hpp"
  "${RKTOPOLOGIESDIR}/se2_topologies.hpp"
  "${RKTOPOLOGIESDIR}/se3_random_samplers.hpp"
  "${RKTOPOLOGIESDIR}/se3_topologies.hpp"
  "${RKTOPOLOGIESDIR}/so3_topologies.hpp"
  "${RKTOPOLOGIESDIR}/temporal_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/temporal_space.hpp"
  "${RKTOPOLOGIESDIR}/time_poisson_topology.hpp"
  "${RKTOPOLOGIESDIR}/time_topology.hpp"
  "${RKTOPOLOGIESDIR}/tuple_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/vect_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/vector_topology.hpp"
)


add_library(reak_topologies STATIC ${TOPOLOGIES_SOURCES})
setup_custom_target(reak_topologies "${SRCROOT}${RKTOPOLOGIESDIR}")
setup_headers("${TOPOLOGIES_HEADERS}" "${RKTOPOLOGIESDIR}")
target_link_libraries(reak_topologies reak_interp reak_core)
if( OpenCV_FOUND )
  target_link_libraries(reak_topologies ${OpenCV_LIBS})
endif()

if(NOT WIN32)

  if( OpenCV_FOUND )
    add_executable(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}/test_sampling.cpp")
    setup_custom_target(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}")
    target_link_libraries(test_sampling ${OpenCV_LIBS})
    target_link_libraries(test_sampling reak_topologies reak_core)
  endif()
endif()








