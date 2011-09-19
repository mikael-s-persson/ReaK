
set(INTERPOLATION_SOURCES 
)

set(INTERPOLATION_HEADERS 
  "${RKINTERPOLATIONDIR}/.hpp"
)

#add_library(reak_interpolation STATIC ${INTERPOLATION_SOURCES})
#setup_custom_target(reak_interpolation "${SRCROOT}${RKINTERPOLATIONDIR}")
setup_headers("${INTERPOLATION_HEADERS}" "${RKINTERPOLATIONDIR}")








