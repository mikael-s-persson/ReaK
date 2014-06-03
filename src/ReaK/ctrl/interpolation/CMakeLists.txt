
set(INTERPOLATION_SOURCES 
  "${SRCROOT}${RKINTERPOLATIONDIR}/sustained_velocity_pulse_detail.cpp"
  "${SRCROOT}${RKINTERPOLATIONDIR}/sustained_velocity_pulse_Ndof_detail.cpp"
  "${SRCROOT}${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_detail.cpp"
  "${SRCROOT}${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_Ndof_detail.cpp"
)

set(INTERPOLATION_HEADERS 
  "${RKINTERPOLATIONDIR}/constant_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/cubic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/discrete_point_path.hpp"
  "${RKINTERPOLATIONDIR}/discrete_point_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/extrapolator_concept.hpp"
  "${RKINTERPOLATIONDIR}/generic_interpolator_factory.hpp"
  "${RKINTERPOLATIONDIR}/interpolated_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/interpolator_concept.hpp"
  "${RKINTERPOLATIONDIR}/linear_interp.hpp"
  "${RKINTERPOLATIONDIR}/path_base.hpp"
  "${RKINTERPOLATIONDIR}/path_wrapper.hpp"
  "${RKINTERPOLATIONDIR}/point_to_point_path.hpp"
  "${RKINTERPOLATIONDIR}/predicted_trajectory_concept.hpp"
  "${RKINTERPOLATIONDIR}/quintic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/sap_metrics.hpp"
  "${RKINTERPOLATIONDIR}/sap_Ndof_metrics.hpp"
  "${RKINTERPOLATIONDIR}/sap_Ndof_reach_topologies.hpp"
  "${RKINTERPOLATIONDIR}/sap_Ndof_samplers.hpp"
  "${RKINTERPOLATIONDIR}/sap_reach_topologies.hpp"
  "${RKINTERPOLATIONDIR}/sap_samplers.hpp"
  "${RKINTERPOLATIONDIR}/seq_path_base.hpp"
  "${RKINTERPOLATIONDIR}/seq_path_wrapper.hpp"
  "${RKINTERPOLATIONDIR}/seq_trajectory_base.hpp"
  "${RKINTERPOLATIONDIR}/seq_trajectory_wrapper.hpp"
  "${RKINTERPOLATIONDIR}/sequential_path_concept.hpp"
  "${RKINTERPOLATIONDIR}/sequential_trajectory_concept.hpp"
  "${RKINTERPOLATIONDIR}/spatial_path_concept.hpp"
  "${RKINTERPOLATIONDIR}/spatial_trajectory_concept.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_detail.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_Ndof.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_Ndof_detail.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse_detail.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse_Ndof.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse_Ndof_detail.hpp"
  "${RKINTERPOLATIONDIR}/svp_metrics.hpp"
  "${RKINTERPOLATIONDIR}/svp_Ndof_metrics.hpp"
  "${RKINTERPOLATIONDIR}/svp_Ndof_reach_topologies.hpp"
  "${RKINTERPOLATIONDIR}/svp_Ndof_samplers.hpp"
  "${RKINTERPOLATIONDIR}/svp_reach_topologies.hpp"
  "${RKINTERPOLATIONDIR}/svp_samplers.hpp"
  "${RKINTERPOLATIONDIR}/trajectory_base.hpp"
  "${RKINTERPOLATIONDIR}/transformed_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/waypoint_container.hpp"
)

add_library(reak_interp STATIC ${INTERPOLATION_SOURCES})
setup_custom_target(reak_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(reak_interp reak_core)
setup_headers("${INTERPOLATION_HEADERS}" "${RKINTERPOLATIONDIR}")

add_executable(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp.cpp")
setup_custom_target(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp reak_topologies reak_interp reak_core)

add_executable(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_lim_interp.cpp")
setup_custom_target(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_lim_interp reak_topologies reak_interp reak_core)

add_executable(test_Ndof_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_Ndof_interp.cpp")
setup_custom_target(test_Ndof_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_Ndof_interp reak_topologies reak_interp reak_core)

add_executable(test_interp_mc "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp_mc.cpp")
setup_custom_target(test_interp_mc "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp_mc reak_topologies reak_interp reak_core)

add_executable(test_interp_topo "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp_topo.cpp")
setup_custom_target(test_interp_topo "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp_topo reak_topologies reak_interp reak_core)

add_executable(test_interp_topo_rand "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp_topo_rand.cpp")
setup_custom_target(test_interp_topo_rand "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp_topo_rand reak_topologies reak_interp reak_core)



