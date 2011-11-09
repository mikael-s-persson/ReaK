
set(INTERPOLATION_SOURCES 
)

set(INTERPOLATION_HEADERS 
  "${RKINTERPOLATIONDIR}/cubic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/linear_interp.hpp"
  "${RKINTERPOLATIONDIR}/quintic_hermite_interp.hpp"
  "${RKINTERPOLATIONDIR}/interpolated_trajectory.hpp"
  "${RKINTERPOLATIONDIR}/waypoint_container.hpp"
  "${RKINTERPOLATIONDIR}/sustained_velocity_pulse.hpp"
)

#add_library(reak_interpolation STATIC ${INTERPOLATION_SOURCES})
#setup_custom_target(reak_interpolation "${SRCROOT}${RKINTERPOLATIONDIR}")
setup_headers("${INTERPOLATION_HEADERS}" "${RKINTERPOLATIONDIR}")

add_executable(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp.cpp")
setup_custom_target(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp ${BOOST_THREAD_LIB})
target_link_libraries(test_interp reak_core)

add_executable(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_lim_interp.cpp")
setup_custom_target(test_lim_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_lim_interp ${BOOST_THREAD_LIB})
target_link_libraries(test_lim_interp reak_core)





