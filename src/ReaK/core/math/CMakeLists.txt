
#set(MATH_SOURCES "${SRCROOT}${RKMATHDIR}/gen_mat_num.cpp")
#set(MATH_SOURCES ${MATH_SOURCES} "${SRCROOT}${RKMATHDIR}/gen_kinetostatics.cpp")
#set(MATH_SOURCES ${MATH_SOURCES} "${SRCROOT}${RKMATHDIR}/gen_time_tagged.cpp")

set(MATH_HEADERS "${RKMATHDIR}/vect_alg.hpp"
                 "${RKMATHDIR}/vect_traits.hpp"
                 "${RKMATHDIR}/vect_concepts.hpp"
                 "${RKMATHDIR}/mat_alg.hpp"
                 "${RKMATHDIR}/stride_iterator.hpp"
                 "${RKMATHDIR}/mat_alg_general.hpp"
                 "${RKMATHDIR}/mat_alg_diagonal.hpp"
                 "${RKMATHDIR}/mat_alg_identity.hpp"
                 "${RKMATHDIR}/mat_alg_lower_triangular.hpp"
                 "${RKMATHDIR}/mat_alg_nil.hpp"
                 "${RKMATHDIR}/mat_alg_orthogonal.hpp"
                 "${RKMATHDIR}/mat_alg_rectangular.hpp"
                 "${RKMATHDIR}/mat_alg_skew_symmetric.hpp"
                 "${RKMATHDIR}/mat_alg_square.hpp"
                 "${RKMATHDIR}/mat_alg_symmetric.hpp"
                 "${RKMATHDIR}/mat_alg_upper_triangular.hpp"
                 "${RKMATHDIR}/mat_vector_adaptor.hpp"
                 "${RKMATHDIR}/mat_views.hpp"
                 "${RKMATHDIR}/mat_slices.hpp"
                 "${RKMATHDIR}/mat_composite_adaptor.hpp"
                 "${RKMATHDIR}/mat_traits.hpp"
                 "${RKMATHDIR}/mat_concepts.hpp"
                 "${RKMATHDIR}/mat_num.hpp"
                 "${RKMATHDIR}/mat_gaussian_elim.hpp"
                 "${RKMATHDIR}/mat_cholesky.hpp"
                 "${RKMATHDIR}/mat_qr_decomp.hpp"
                 "${RKMATHDIR}/mat_jacobi_method.hpp"
                 "${RKMATHDIR}/mat_svd_method.hpp"
                 "${RKMATHDIR}/mat_norms.hpp"
                 "${RKMATHDIR}/mat_exp_methods.hpp"
                 "${RKMATHDIR}/mat_star_product.hpp"
                 "${RKMATHDIR}/mat_hess_decomp.hpp"
                 "${RKMATHDIR}/mat_givens_rot.hpp"
                 "${RKMATHDIR}/mat_householder.hpp"
                 "${RKMATHDIR}/mat_schur_decomp.hpp"
                 "${RKMATHDIR}/complex_math.hpp"
                 "${RKMATHDIR}/rotations.hpp"
                 "${RKMATHDIR}/rotations_2D.hpp"
                 "${RKMATHDIR}/rotations_3D.hpp"
                 "${RKMATHDIR}/kinetostatics.hpp"
                 "${RKMATHDIR}/gen_coord.hpp"
                 "${RKMATHDIR}/pose_2D.hpp"
                 "${RKMATHDIR}/frame_2D.hpp"
                 "${RKMATHDIR}/pose_3D.hpp"
                 "${RKMATHDIR}/frame_3D.hpp"
                 "${RKMATHDIR}/motion_jacobians.hpp"
)

#add_library(reak_math STATIC ${MATH_SOURCES})
#setup_custom_target(reak_math "${SRCROOT}${RKMATHDIR}")

setup_headers("${MATH_HEADERS}" "${RKMATHDIR}")

add_executable(test_math "${SRCROOT}${RKMATHDIR}/test.cpp")
setup_custom_target(test_math "${SRCROOT}${RKMATHDIR}")
target_link_libraries(test_math reak_rtti "${BOOST_DATE_TIME_LIB}")

add_executable(test_math_quat "${SRCROOT}${RKMATHDIR}/test_quat.cpp")
setup_custom_target(test_math_quat "${SRCROOT}${RKMATHDIR}")
target_link_libraries(test_math_quat reak_rtti)

add_executable(test_math_views "${SRCROOT}${RKMATHDIR}/test_views.cpp")
setup_custom_target(test_math_views "${SRCROOT}${RKMATHDIR}")
target_link_libraries(test_math_views reak_rtti)


include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




