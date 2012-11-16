
set(INTERPOLATION_SOURCES 
)

set(INTERPOLATION_HEADERS 
  "${RKINTERPOLATIONDIR}/cubic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/generic_interpolator_factory.hpp"
  "${RKINTERPOLATIONDIR}/interpolated_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/linear_interp.hpp"
  "${RKINTERPOLATIONDIR}/quintic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse_detail.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse.hpp"
  "${RKINTERPOLATIONDIR}/sustained_acceleration_pulse_detail.hpp"
  "${RKINTERPOLATIONDIR}/waypoint_container.hpp"
)

#add_library(reak_interpolation STATIC ${INTERPOLATION_SOURCES})
#setup_custom_target(reak_interpolation "${SRCROOT}${RKINTERPOLATIONDIR}")
setup_headers("${INTERPOLATION_HEADERS}" "${RKINTERPOLATIONDIR}")

add_executable(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp.cpp")
setup_custom_target(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp ${Boost_LIBRARIES})
target_link_libraries(test_interp reak_topologies reak_core ${EXTRA_SYSTEM_LIBS})

add_executable(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_lim_interp.cpp")
setup_custom_target(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_lim_interp ${Boost_LIBRARIES})
target_link_libraries(test_lim_interp reak_topologies reak_core ${EXTRA_SYSTEM_LIBS})





