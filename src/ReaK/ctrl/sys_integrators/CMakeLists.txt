
set(SYS_INTEGRATORS_SOURCES 
)


set(SYS_INTEGRATORS_HEADERS 
  "${RKSYSINTEGRATORSDIR}/adams_BM3_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/adams_BM5_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/dormand_prince45_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/euler_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/fehlberg45_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/hamming_iter_mod_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/hamming_mod_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/midpoint_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/runge_kutta4_integrator_sys.hpp"
  "${RKSYSINTEGRATORSDIR}/runge_kutta5_integrator_sys.hpp"
)


# add_library(reak_topologies STATIC ${TOPOLOGIES_SOURCES})
# setup_custom_target(reak_topologies "${SRCROOT}${RKTOPOLOGIESDIR}")
# target_link_libraries(reak_topologies reak_core)
setup_headers("${SYS_INTEGRATORS_HEADERS}" "${RKSYSINTEGRATORSDIR}")









