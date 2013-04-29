
set(KINETOSTATICS_SOURCES 
  "${SRCROOT}${RKKINETOSTATICSDIR}/kinetostatics_externs.cpp"
)

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

add_library(reak_kinetostatics STATIC ${KINETOSTATICS_SOURCES})
setup_custom_target(reak_kinetostatics "${SRCROOT}${RKKINETOSTATICSDIR}")
target_link_libraries(reak_kinetostatics reak_lin_alg reak_rtti ${EXTRA_SYSTEM_LIBS})

setup_headers("${KINETOSTATICS_HEADERS}" "${RKKINETOSTATICSDIR}")

add_executable(unit_test_rotations "${SRCROOT}${RKKINETOSTATICSDIR}/unit_test_rotations.cpp")
setup_custom_test_program(unit_test_rotations "${SRCROOT}${RKKINETOSTATICSDIR}")
target_link_libraries(unit_test_rotations reak_kinetostatics reak_lin_alg reak_rtti ${EXTRA_SYSTEM_LIBS})
target_link_libraries(unit_test_rotations ${Boost_LIBRARIES})

add_executable(unit_test_quat_alg "${SRCROOT}${RKKINETOSTATICSDIR}/unit_test_quat_alg.cpp")
setup_custom_test_program(unit_test_quat_alg "${SRCROOT}${RKKINETOSTATICSDIR}")
target_link_libraries(unit_test_quat_alg reak_kinetostatics reak_lin_alg reak_rtti ${EXTRA_SYSTEM_LIBS})
target_link_libraries(unit_test_quat_alg ${Boost_LIBRARIES})





