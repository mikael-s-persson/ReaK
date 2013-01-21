
set(TOPOLOGIES_SOURCES 
  "${SRCROOT}${RKTOPOLOGIESDIR}/hyperball_topology_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/hyperbox_topology_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/joint_space_limits_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/joint_space_topologies_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/line_topology_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_svp_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_sap_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_dynamic_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_svp_dynamic_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/manip_free_sap_dynamic_workspace_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/Ndof_limits_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/ptrobot2D_test_world.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/se2_topologies_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/se3_topologies_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/so3_topologies_ext.cpp"
  "${SRCROOT}${RKTOPOLOGIESDIR}/vector_topology_ext.cpp"
)

set(TOPOLOGIES_HEADERS 
  "${RKTOPOLOGIESDIR}/basic_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/default_random_sampler.hpp"
  "${RKTOPOLOGIESDIR}/differentiable_space.hpp"
  "${RKTOPOLOGIESDIR}/hyperball_topology.hpp"
  "${RKTOPOLOGIESDIR}/hyperbox_topology.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits_detail.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits_ext.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_topologies.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_topologies_ext.hpp"
  "${RKTOPOLOGIESDIR}/line_topology.hpp"
  "${RKTOPOLOGIESDIR}/manip_DK_proxy_updater.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_svp_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_sap_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_workspace_tsppf.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_dynamic_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_svp_dynamic_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_free_sap_dynamic_workspace.hpp"
  "${RKTOPOLOGIESDIR}/manip_workspaces.hpp"
  "${RKTOPOLOGIESDIR}/metric_space_tuple.hpp"
  "${RKTOPOLOGIESDIR}/Ndof_spaces.hpp"
  "${RKTOPOLOGIESDIR}/Ndof_limits.hpp"
  "${RKTOPOLOGIESDIR}/proxy_model_updater.hpp"
  "${RKTOPOLOGIESDIR}/ptrobot2D_test_world.hpp"
  "${RKTOPOLOGIESDIR}/rate_limited_spaces.hpp"
  "${RKTOPOLOGIESDIR}/reachability_space.hpp"
  "${RKTOPOLOGIESDIR}/se2_topologies.hpp"
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
target_link_libraries(reak_topologies reak_interp reak_core ${EXTRA_SYSTEM_LIBS})
if( OpenCV_FOUND )
  target_link_libraries(reak_topologies ${OpenCV_LIBS})
endif()

if(NOT WIN32)

  if( OpenCV_FOUND )
    add_executable(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}/test_sampling.cpp")
    setup_custom_target(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}")
    #target_link_libraries(test_adstar ${Boost_LIBRARIES})
    target_link_libraries(test_sampling ${OpenCV_LIBS})
    target_link_libraries(test_sampling reak_topologies reak_core ${EXTRA_SYSTEM_LIBS})
  endif()
endif()








