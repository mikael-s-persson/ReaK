
set(OPTIM_SOURCES 
)

set(OPTIM_HEADERS 
  "${RKOPTIMDIR}/augmented_lagrangian_methods.hpp"
  "${RKOPTIMDIR}/conjugate_gradient_methods.hpp"
  "${RKOPTIMDIR}/finite_diff_jacobians.hpp"
  "${RKOPTIMDIR}/finite_differences.hpp"
  "${RKOPTIMDIR}/function_types.hpp"
  "${RKOPTIMDIR}/gauss_newton_method.hpp"
  "${RKOPTIMDIR}/hessian_approx_update.hpp"
  "${RKOPTIMDIR}/jacobian_helper_functions.hpp"
  "${RKOPTIMDIR}/jacobian_transpose_method.hpp"
  "${RKOPTIMDIR}/levenberg_marquardt_method.hpp"
  "${RKOPTIMDIR}/limit_functions.hpp"
  "${RKOPTIMDIR}/line_search.hpp"
  "${RKOPTIMDIR}/mehrotra_method.hpp"
  "${RKOPTIMDIR}/nelder_mead_method.hpp"
  "${RKOPTIMDIR}/newton_methods.hpp"
  "${RKOPTIMDIR}/newton_search_directions.hpp"
  "${RKOPTIMDIR}/nl_interior_points_methods.hpp"
  "${RKOPTIMDIR}/optim_exceptions.hpp"
  "${RKOPTIMDIR}/quadratic_programs.hpp"
  "${RKOPTIMDIR}/quasi_newton_methods.hpp"
  "${RKOPTIMDIR}/sequential_qp_methods.hpp"
  "${RKOPTIMDIR}/simplex_method.hpp"
  "${RKOPTIMDIR}/trust_region_search.hpp"
)

setup_headers("${OPTIM_HEADERS}" "${RKOPTIMDIR}")

add_executable(test_optim "${SRCROOT}${RKOPTIMDIR}/test.cpp")
setup_custom_target(test_optim "${SRCROOT}${RKOPTIMDIR}")
target_link_libraries(test_optim reak_core)


add_executable(test_optim_nllsq "${SRCROOT}${RKOPTIMDIR}/test_nllsq.cpp")
setup_custom_target(test_optim_nllsq "${SRCROOT}${RKOPTIMDIR}")
target_link_libraries(test_optim_nllsq reak_core)


add_executable(test_optim_lp "${SRCROOT}${RKOPTIMDIR}/test_lp.cpp")
setup_custom_target(test_optim_lp "${SRCROOT}${RKOPTIMDIR}")
target_link_libraries(test_optim_lp reak_core)


add_executable(test_optim_qp "${SRCROOT}${RKOPTIMDIR}/test_qp.cpp")
setup_custom_target(test_optim_qp "${SRCROOT}${RKOPTIMDIR}")
target_link_libraries(test_optim_qp reak_core)


add_executable(test_optim_nlp "${SRCROOT}${RKOPTIMDIR}/test_nlp.cpp")
setup_custom_target(test_optim_nlp "${SRCROOT}${RKOPTIMDIR}")
target_link_libraries(test_optim_nlp reak_core)

