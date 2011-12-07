
set(OPTIM_SOURCES 
)

set(OPTIM_HEADERS 
  "${RKOPTIMDIR}/optim_exceptions.hpp"
  "${RKOPTIMDIR}/simplex_method.hpp"
  "${RKOPTIMDIR}/line_search.hpp"
  "${RKOPTIMDIR}/trust_region_search.hpp"
  "${RKOPTIMDIR}/finite_diff_jacobians.hpp"
  "${RKOPTIMDIR}/conjugate_gradient_methods.hpp"
  "${RKOPTIMDIR}/gauss_newton_method.hpp"
  "${RKOPTIMDIR}/jacobian_helper_functions.hpp"
  "${RKOPTIMDIR}/hessian_approx_update.hpp"
  "${RKOPTIMDIR}/limit_functions.hpp"
  "${RKOPTIMDIR}/nelder_mead_method.hpp"
  "${RKOPTIMDIR}/quasi_newton_methods.hpp"
  "${RKOPTIMDIR}/newton_search_directions.hpp"
  "${RKOPTIMDIR}/newton_methods.hpp"
  "${RKOPTIMDIR}/jacobian_transpose_method.hpp"
  "${RKOPTIMDIR}/levenberg_marquardt_method.hpp"
)

setup_headers("${OPTIM_HEADERS}" "${RKOPTIMDIR}")

add_executable(test_optim "${SRCROOT}${RKOPTIMDIR}/test.cpp")
setup_custom_target(test_optim "${SRCROOT}${RKOPTIMDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

target_link_libraries(test_optim reak_rtti)
target_link_libraries(test_optim ${BOOST_THREAD_LIB})


add_executable(test_optim_nllsq "${SRCROOT}${RKOPTIMDIR}/test_nllsq.cpp")
setup_custom_target(test_optim_nllsq "${SRCROOT}${RKOPTIMDIR}")

target_link_libraries(test_optim_nllsq reak_rtti)
target_link_libraries(test_optim_nllsq ${BOOST_THREAD_LIB})


add_executable(test_optim_lp "${SRCROOT}${RKOPTIMDIR}/test_lp.cpp")
setup_custom_target(test_optim_lp "${SRCROOT}${RKOPTIMDIR}")

target_link_libraries(test_optim_lp reak_core)
target_link_libraries(test_optim_lp ${BOOST_THREAD_LIB})

