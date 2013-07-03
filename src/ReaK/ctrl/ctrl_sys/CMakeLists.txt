
set(CTRL_SYS_SOURCES 
)

set(CTRL_SYS_HEADERS 
  "${RKCTRLSYSDIR}/aggregate_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/belief_state_concept.hpp"
  "${RKCTRLSYSDIR}/belief_state_predictor.hpp"
  "${RKCTRLSYSDIR}/covar_topology.hpp"
  "${RKCTRLSYSDIR}/covariance_concept.hpp"
  "${RKCTRLSYSDIR}/covariance_info_matrix.hpp"
  "${RKCTRLSYSDIR}/covariance_matrix.hpp"
  "${RKCTRLSYSDIR}/decomp_covariance_matrix.hpp"
  "${RKCTRLSYSDIR}/discrete_linear_sss_concept.hpp"
  "${RKCTRLSYSDIR}/discrete_sss_concept.hpp"
  "${RKCTRLSYSDIR}/discretized_lti_sys.hpp"
  "${RKCTRLSYSDIR}/gaussian_belief_space.hpp"
  "${RKCTRLSYSDIR}/gaussian_belief_state.hpp"
  "${RKCTRLSYSDIR}/hybrid_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/invariant_aggregate_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/invariant_kalman_bucy_filter.hpp"
  "${RKCTRLSYSDIR}/invariant_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/invariant_symplectic_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/invariant_system_concept.hpp"
  "${RKCTRLSYSDIR}/kalman_bucy_filter.hpp"
  "${RKCTRLSYSDIR}/kalman_filter.hpp"
  "${RKCTRLSYSDIR}/kte_nl_system.hpp"
  "${RKCTRLSYSDIR}/linear_ss_system_concept.hpp"
  "${RKCTRLSYSDIR}/lti_discrete_sys.hpp"
  "${RKCTRLSYSDIR}/lti_ss_system.hpp"
  "${RKCTRLSYSDIR}/num_int_dtnl_system.hpp"
  "${RKCTRLSYSDIR}/sss_exceptions.hpp"
  "${RKCTRLSYSDIR}/state_estimator_concept.hpp"
  "${RKCTRLSYSDIR}/state_space_sys_concept.hpp"
  "${RKCTRLSYSDIR}/state_vector_concept.hpp"
  "${RKCTRLSYSDIR}/symplectic_kalman_filter.hpp"
  "${RKCTRLSYSDIR}/unscented_kalman_filter.hpp"
)

#add_library(reak_path_planning STATIC ${PATH_PLANNING_SOURCES})
#setup_custom_target(reak_path_planning "${SRCROOT}${RKPATHPLANNINGDIR}")
setup_headers("${CTRL_SYS_HEADERS}" "${RKCTRLSYSDIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")





