

reak_cc_library(
  NAME
    kte
  SRCS
    "damper.cc"
    "driving_actuator.cc"
    "dry_revolute_joint.cc"
    "flexible_beam.cc"
    "force_actuator.cc"
    "free_joints.cc"
    "inertia.cc"
    "inertial_beam.cc"
    "joint_backlash.cc"
    "joint_friction.cc"
    "kte_chain_geometry.cc"
    "line_point_mindist.cc"
    "mass_matrix_calculator.cc"
    "plane_point_mindist.cc"
    "prismatic_joint.cc"
    "revolute_joint.cc"
    "rigid_link.cc"
    "spring.cc"
    "torsion_damper.cc"
    "torsion_spring.cc"
    "virtual_kte_interface.cc"
    "vmc_revolute_joint.cc"
  HDRS
    "damper.h"
    "driving_actuator.h"
    "dry_revolute_joint.h"
    "flexible_beam.h"
    "force_actuator.h"
    "free_joints.h"
    "inertia.h"
    "inertial_beam.h"
    "jacobian_joint_map.h"
    "joint_backlash.h"
    "joint_friction.h"
    "kte_chain_geometry.h"
    "kte_chain_visitation.h"
    "kte_ext_mappings.h"
    "kte_map_chain.h"
    "kte_map.h"
    "kte_system_input.h"
    "kte_system_output.h"
    "line_point_mindist.h"
    "mass_matrix_calculator.h"
    "plane_point_mindist.h"
    "prismatic_joint.h"
    "reacting_kte.h"
    "revolute_joint.h"
    "rigid_link.h"
    "spring.h"
    "state_controls.h"
    "state_measures.h"
    "torsion_damper.h"
    "torsion_spring.h"
    "virtual_kte_interface.h"
    "vmc_revolute_joint.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::kinetostatics
  PUBLIC
)

if(Python3_FOUND)
  reak_cc_library(
    NAME
      py_kte
    SRCS
      "py_kte.cc"
    DEPS
      ReaK::kte
      ${PYTHON_LIBRARIES}
      ${Boost_LIBRARIES}
      ${EXTRA_SYSTEM_LIBS}
    PUBLIC
  )
endif()

reak_cc_binary(
  NAME
    pendulum_advanced_test
  SRCS
    "pendulum_advanced_test.cc"
  DEPS
    ReaK::base
    ReaK::serialization
    ReaK::recorders
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::kte
)

reak_cc_binary(
  NAME
    pendulum_basic_test
  SRCS
    "pendulum_basic_test.cc"
  DEPS
    ReaK::base
    ReaK::serialization
    ReaK::recorders
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::kte
)

# reak_cc_binary(
#   NAME
#     pendulum_ctrl_test
#   SRCS
#     "pendulum_ctrl_test.cc"
#   DEPS
#     ReaK::base
#     ReaK::serialization
#     ReaK::recorders
#     ReaK::lin_alg
#     ReaK::kinetostatics
#     ReaK::kte
# )
