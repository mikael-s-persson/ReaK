
set(INTERPOLATION_SOURCES 
)

set(INTERPOLATION_HEADERS 
  "${RKINTERPOLATIONDIR}/.hpp"
)

#add_library(reak_interpolation STATIC ${INTERPOLATION_SOURCES})
#setup_custom_target(reak_interpolation "${SRCROOT}${RKINTERPOLATIONDIR}")
setup_headers("${INTERPOLATION_HEADERS}" "${RKINTERPOLATIONDIR}")

add_executable(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}/test_interp.cpp")
setup_custom_target(test_interp "${SRCROOT}${RKINTERPOLATIONDIR}")
target_link_libraries(test_interp reak_core)






