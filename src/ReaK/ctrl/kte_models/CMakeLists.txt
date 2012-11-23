
set(KTEMODELS_SOURCES 
  "${SRCROOT}${RKKTEMODELSDIR}/manip_dynamics_model.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_kinematics_helper.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_kinematics_model.cpp"
  "${SRCROOT}${RKKTEMODELSDIR}/manip_clik_calculator.cpp"
)


set(KTEMODELS_HEADERS 
  "${RKKTEMODELSDIR}/direct_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/inverse_dynamics_model.hpp"
  "${RKKTEMODELSDIR}/inverse_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/manip_dynamics_model.hpp"
  "${RKKTEMODELSDIR}/manip_kinematics_helper.hpp"
  "${RKKTEMODELSDIR}/manip_kinematics_model.hpp"
  "${RKKTEMODELSDIR}/manip_clik_calculator.hpp"
)

add_library(reak_kte_models STATIC ${KTEMODELS_SOURCES})
setup_custom_target(reak_kte_models "${SRCROOT}${RKKTEMODELSDIR}")
setup_headers("${KTEMODELS_HEADERS}" "${RKKTEMODELSDIR}")
target_link_libraries(reak_kte_models reak_mbd_kte reak_core ${EXTRA_SYSTEM_LIBS})

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




