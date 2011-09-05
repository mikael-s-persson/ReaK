
#set(LINALG_SOURCES "${SRCROOT}${RKLINALGDIR}/gen_mat_num.cpp")

set(LINALG_HEADERS 
                 "${RKLINALGDIR}/complex_math.hpp"
                 "${RKLINALGDIR}/mat_alg.hpp"
                 "${RKLINALGDIR}/mat_alg_diagonal.hpp"
                 "${RKLINALGDIR}/mat_alg_general.hpp"
                 "${RKLINALGDIR}/mat_alg_identity.hpp"
                 "${RKLINALGDIR}/mat_alg_lower_triangular.hpp"
                 "${RKLINALGDIR}/mat_alg_nil.hpp"
                 "${RKLINALGDIR}/mat_alg_orthogonal.hpp"
                 "${RKLINALGDIR}/mat_alg_rectangular.hpp"
                 "${RKLINALGDIR}/mat_alg_rectangular_fixed.hpp"
                 "${RKLINALGDIR}/mat_alg_skew_symmetric.hpp"
                 "${RKLINALGDIR}/mat_alg_square.hpp"
                 "${RKLINALGDIR}/mat_alg_symmetric.hpp"
                 "${RKLINALGDIR}/mat_alg_upper_triangular.hpp"
                 "${RKLINALGDIR}/mat_balance.hpp"
                 "${RKLINALGDIR}/mat_cholesky.hpp"
                 "${RKLINALGDIR}/mat_composite_adaptor.hpp"
                 "${RKLINALGDIR}/mat_concepts.hpp"
                 "${RKLINALGDIR}/mat_esn_expressions.hpp"
                 "${RKLINALGDIR}/mat_exp_methods.hpp"
                 "${RKLINALGDIR}/mat_gaussian_elim.hpp"
                 "${RKLINALGDIR}/mat_givens_rot.hpp"
                 "${RKLINALGDIR}/mat_hess_decomp.hpp"
                 "${RKLINALGDIR}/mat_householder.hpp"
                 "${RKLINALGDIR}/mat_jacobi_method.hpp"
                 "${RKLINALGDIR}/mat_norms.hpp"
                 "${RKLINALGDIR}/mat_num.hpp"
                 "${RKLINALGDIR}/mat_num_exceptions.hpp"
                 "${RKLINALGDIR}/mat_qr_decomp.hpp"
                 "${RKLINALGDIR}/mat_schur_decomp.hpp"
                 "${RKLINALGDIR}/mat_slices.hpp"
                 "${RKLINALGDIR}/mat_star_product.hpp"
                 "${RKLINALGDIR}/mat_svd_method.hpp"
                 "${RKLINALGDIR}/mat_traits.hpp"
                 "${RKLINALGDIR}/mat_vector_adaptor.hpp"
                 "${RKLINALGDIR}/mat_views.hpp"
                 "${RKLINALGDIR}/stride_iterator.hpp"
                 "${RKLINALGDIR}/vect_alg.hpp"
                 "${RKLINALGDIR}/vect_concepts.hpp"
                 "${RKLINALGDIR}/vect_traits.hpp"
)

#add_library(reak_lin_alg STATIC ${LINALG_SOURCES})
#setup_custom_target(reak_lin_alg "${SRCROOT}${RKLINALGDIR}")

setup_headers("${LINALG_HEADERS}" "${RKLINALGDIR}")

add_executable(test_lin_alg "${SRCROOT}${RKLINALGDIR}/test.cpp")
setup_custom_target(test_lin_alg "${SRCROOT}${RKLINALGDIR}")
target_link_libraries(test_lin_alg reak_rtti "${BOOST_DATE_TIME_LIB}")

add_executable(test_mat_views "${SRCROOT}${RKLINALGDIR}/test_views.cpp")
setup_custom_target(test_mat_views "${SRCROOT}${RKLINALGDIR}")
target_link_libraries(test_mat_views reak_rtti)






