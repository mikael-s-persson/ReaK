
#set(KINETOSTATICS_SOURCES "${SRCROOT}${RKMATHDIR}/gen_mat_num.cpp")

set(KINETOSTATICS_HEADERS 
                 "${RKKINETOSTATICSDIR}/frame_2D.hpp"
                 "${RKKINETOSTATICSDIR}/frame_3D.hpp"
                 "${RKKINETOSTATICSDIR}/gen_coord.hpp"
                 "${RKKINETOSTATICSDIR}/kinetostatics.hpp"
                 "${RKKINETOSTATICSDIR}/motion_jacobians.hpp"
                 "${RKKINETOSTATICSDIR}/pose_2D.hpp"
                 "${RKKINETOSTATICSDIR}/pose_3D.hpp"
                 "${RKKINETOSTATICSDIR}/quat_alg.hpp"
                 "${RKKINETOSTATICSDIR}/rotations.hpp"
                 "${RKKINETOSTATICSDIR}/rotations_2D.hpp"
                 "${RKKINETOSTATICSDIR}/rotations_3D.hpp"
)

#add_library(reak_kinetostatics STATIC ${KINETOSTATICS_SOURCES})
#setup_custom_target(reak_kinetostatics "${SRCROOT}${RKKINETOSTATICSDIR}")

setup_headers("${KINETOSTATICS_HEADERS}" "${RKKINETOSTATICSDIR}")

add_executable(test_kin "${SRCROOT}${RKKINETOSTATICSDIR}/test.cpp")
setup_custom_target(test_kin "${SRCROOT}${RKKINETOSTATICSDIR}")
target_link_libraries(test_kin reak_rtti)

add_executable(test_quat "${SRCROOT}${RKKINETOSTATICSDIR}/test_quat.cpp")
setup_custom_target(test_quat "${SRCROOT}${RKKINETOSTATICSDIR}")
target_link_libraries(test_quat reak_rtti)






