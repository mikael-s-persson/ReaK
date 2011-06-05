
set(INTEGRATORS_SOURCES )

set(INTEGRATORS_HEADERS 
  "${RKINTEGRATORSDIR}/integrator.hpp"
  "${RKINTEGRATORSDIR}/integration_exceptions.hpp"
  "${RKINTEGRATORSDIR}/fixed_step_integrators.hpp"
  "${RKINTEGRATORSDIR}/pred_corr_integrators.hpp"
  "${RKINTEGRATORSDIR}/variable_step_integrators.hpp"
)

#add_library(reak_integrators STATIC ${INTEGRATORS_SOURCES})
#setup_custom_target(reak_integrators "${SRCROOT}${RKINTEGRATORSDIR}")

setup_headers("${INTEGRATORS_HEADERS}" "${RKINTEGRATORSDIR}")

#add_executable(test_math "${SRCROOT}${RKMATHDIR}/test.cpp")
#setup_custom_target(test_math "${SRCROOT}${RKMATHDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




