
set(MBDKTE_SOURCES 
  "${SRCROOT}${RKMBDKTEDIR}/inertia.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/damper.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/spring.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/revolute_joint.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/prismatic_joint.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/rigid_link.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/torsion_damper.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/torsion_spring.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/driving_actuator.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/force_actuator.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/dry_revolute_joint.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/joint_backlash.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/joint_friction.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/flexible_beam.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/inertial_beam.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/line_point_mindist.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/plane_point_mindist.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/mass_matrix_calculator.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/virtual_kte_interface.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/vmc_revolute_joint.cpp"
  "${SRCROOT}${RKMBDKTEDIR}/free_joints.cpp"
)


set(MBDKTE_HEADERS 
  "${RKMBDKTEDIR}/kte_map.hpp"
  "${RKMBDKTEDIR}/kte_map_chain.hpp"
  "${RKMBDKTEDIR}/kte_ext_mappings.hpp"
  "${RKMBDKTEDIR}/damper.hpp"
  "${RKMBDKTEDIR}/driving_actuator.hpp"
  "${RKMBDKTEDIR}/dry_revolute_joint.hpp"
  "${RKMBDKTEDIR}/flexible_beam.hpp"
  "${RKMBDKTEDIR}/force_actuator.hpp"
  "${RKMBDKTEDIR}/inertia.hpp"
  "${RKMBDKTEDIR}/inertial_beam.hpp"
  "${RKMBDKTEDIR}/jacobian_joint_map.hpp"
  "${RKMBDKTEDIR}/joint_backlash.hpp"
  "${RKMBDKTEDIR}/joint_friction.hpp"
  "${RKMBDKTEDIR}/line_point_mindist.hpp"
  "${RKMBDKTEDIR}/mass_matrix_calculator.hpp"
  "${RKMBDKTEDIR}/plane_point_mindist.hpp"
  "${RKMBDKTEDIR}/prismatic_joint.hpp"
  "${RKMBDKTEDIR}/reacting_kte.hpp"
  "${RKMBDKTEDIR}/revolute_joint.hpp"
  "${RKMBDKTEDIR}/rigid_link.hpp"
  "${RKMBDKTEDIR}/spring.hpp"
  "${RKMBDKTEDIR}/torsion_damper.hpp"
  "${RKMBDKTEDIR}/torsion_spring.hpp"
  "${RKMBDKTEDIR}/virtual_kte_interface.hpp"
  "${RKMBDKTEDIR}/vmc_revolute_joint.hpp"
  "${RKMBDKTEDIR}/kte_system_output.hpp"
  "${RKMBDKTEDIR}/kte_system_input.hpp"
  "${RKMBDKTEDIR}/state_measure.hpp"
  "${RKMBDKTEDIR}/state_controls.hpp"
  "${RKMBDKTEDIR}/free_joints.hpp"
)

add_library(reak_mbd_kte STATIC ${MBDKTE_SOURCES})
setup_custom_target(reak_mbd_kte "${SRCROOT}${RKMBDKTEDIR}")
setup_headers("${MBDKTE_HEADERS}" "${RKMBDKTEDIR}")

add_executable(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}/test_bm.cpp")
setup_custom_target(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}")

target_link_libraries(test_mbd_basic ${BOOST_THREAD_LIB})
target_link_libraries(test_mbd_basic reak_mbd_kte reak_core)

add_executable(test_mbd_adv "${SRCROOT}${RKMBDKTEDIR}/test_am.cpp")
setup_custom_target(test_mbd_adv "${SRCROOT}${RKMBDKTEDIR}")

target_link_libraries(test_mbd_adv ${BOOST_THREAD_LIB})
target_link_libraries(test_mbd_adv reak_mbd_kte reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




