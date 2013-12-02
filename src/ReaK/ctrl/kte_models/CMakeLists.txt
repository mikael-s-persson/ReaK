
set(KTEMODELS_SOURCES 
  "${SRCROOT}${RKKTEMODELSDIR}/free_floating_platform.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_3R3R_arm.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_3R_arm.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_clik_calculator.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_dynamics_model.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_ERA_arm.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_kinematics_helper.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_kinematics_model.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_P3R3R_arm.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_SCARA_arm.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_SSRMS_arm.cpp"
)


set(KTEMODELS_HEADERS 
  "${RKKTEMODELSDIR}/direct_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/free_floating_platform.hpp"
  "${RKKTEMODELSDIR}/inverse_dynamics_model.hpp"
  "${RKKTEMODELSDIR}/inverse_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/manip_3R3R_arm.hpp"
  "${RKKTEMODELSDIR}/manip_3R_arm.hpp"
  "${RKKTEMODELSDIR}/manip_clik_calculator.hpp"
  "${RKKTEMODELSDIR}/manip_dynamics_model.hpp"
  "${RKKTEMODELSDIR}/manip_ERA_arm.hpp"
  "${RKKTEMODELSDIR}/manip_kinematics_helper.hpp"
  "${RKKTEMODELSDIR}/manip_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/manip_P3R3R_arm.hpp"
  "${RKKTEMODELSDIR}/manip_SCARA_arm.hpp"
  "${RKKTEMODELSDIR}/manip_SSRMS_arm.hpp"
)

add_library(reakobj_kte_models OBJECT ${KTEMODELS_SOURCES})

setup_headers("${KTEMODELS_HEADERS}" "${RKKTEMODELSDIR}")


set(KTEGEOMMODELS_SOURCES 
  "${SRCROOT}${RKKTEMODELSDIR}/chaser_target_model_data.cpp"
)


set(KTEGEOMMODELS_HEADERS 
  "${RKKTEMODELSDIR}/chaser_target_model_data.hpp"
)

add_library(reakobj_kte_geom_models OBJECT ${KTEGEOMMODELS_SOURCES})

setup_headers("${KTEGEOMMODELS_HEADERS}" "${RKKTEMODELSDIR}")


