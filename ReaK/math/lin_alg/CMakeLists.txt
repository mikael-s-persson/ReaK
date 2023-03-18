

reak_cc_library(
  NAME
    lin_alg
  HDRS
    "arithmetic_tuple.h"
    "complex_math.h"
    "mat_alg_diagonal.h"
    "mat_alg_general.h"
    "mat_alg.h"
    "mat_alg_identity.h"
    "mat_alg_lower_triangular.h"
    "mat_alg_nil.h"
    "mat_alg_orthogonal.h"
    "mat_alg_permutation.h"
    "mat_alg_rectangular.h"
    "mat_alg_scalar.h"
    "mat_alg_skew_symmetric.h"
    "mat_alg_square.h"
    "mat_alg_symmetric.h"
    "mat_alg_upper_triangular.h"
    "mat_are_solver.h"
    "mat_balance.h"
    "mat_cholesky.h"
    "mat_comparisons.h"
    "mat_composite_adaptor.h"
    "mat_concepts.h"
    "mat_ctrl_decomp.h"
    "mat_damped_matrix.h"
    "mat_exp_methods.h"
    "mat_gaussian_elim.h"
    "mat_givens_rot.h"
    "mat_hess_decomp.h"
    "mat_householder.h"
    "mat_jacobi_method.h"
    "mat_norms.h"
    "mat_num_exceptions.h"
    "mat_num.h"
    "mat_operators.h"
    "mat_op_results.h"
    "mat_qr_decomp.h"
    "mat_schur_decomp.h"
    "mat_slices.h"
    "mat_star_product.h"
    "mat_svd_method.h"
    "mat_traits.h"
    "mat_transpose_view.h"
    "mat_vector_adaptor.h"
    "mat_views.h"
    "named_compound_records.h"
    "stride_iterator.h"
    "vect_alg.h"
    "vect_concepts.h"
    "vect_index_iterator.h"
    "vect_traits.h"
    "vect_views.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
  PUBLIC
)

reak_cc_library(
  NAME
    lin_alg_matchers
  HDRS
    "mat_matchers.h"
    "vect_matchers.h"
  DEPS
    GTest::gmock
    ReaK::lin_alg
  PUBLIC
)

reak_cc_test(
  NAME
    vect_alg_test
  SRCS
    "vect_alg_test.cc"
  DEPS
    ReaK::lin_alg
    GTest::gtest_main
)

reak_cc_test(
  NAME
    mat_alg_test
  SRCS
    "mat_alg_test.cc"
  DEPS
    ReaK::lin_alg
    ReaK::lin_alg_matchers
    GTest::gtest_main
)

reak_cc_test(
  NAME
    mat_num_test
  SRCS
    "mat_num_test.cc"
  DEPS
    ReaK::lin_alg
    ReaK::lin_alg_matchers
    GTest::gtest_main
)

reak_cc_test(
  NAME
    mat_are_test
  SRCS
    "mat_are_test.cc"
  DEPS
    ReaK::lin_alg
    GTest::gtest_main
)
add_custom_command(TARGET reak_mat_are_test PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  make_directory "${CMAKE_CURRENT_BINARY_DIR}/are_data")
add_custom_command(TARGET reak_mat_are_test PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  copy "${CMAKE_CURRENT_SOURCE_DIR}/carex_m/carex_data.txt" "${CMAKE_CURRENT_BINARY_DIR}/are_data/")
add_custom_command(TARGET reak_mat_are_test PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E 
  copy "${CMAKE_CURRENT_SOURCE_DIR}/darex_m/darex_data.txt" "${CMAKE_CURRENT_BINARY_DIR}/are_data/")

reak_cc_binary(
  NAME
    mat_num_perf
  SRCS
    "mat_num_perf.cc"
  DEPS
    absl::log
    ReaK::lin_alg
)

if(Python3_FOUND)
  reak_cc_library(
    NAME
      py_lin_alg
    SRCS
      "py_vect_alg.cc"
    DEPS
      ReaK::lin_alg
      ${PYTHON_LIBRARIES}
      ${Boost_LIBRARIES}
      ${EXTRA_SYSTEM_LIBS}
    PUBLIC
  )
endif()
