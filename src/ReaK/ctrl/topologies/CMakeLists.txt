
set(TOPOLOGIES_SOURCES 
)

set(TOPOLOGIES_HEADERS 
  "${RKTOPOLOGIESDIR}/line_topology.hpp"
  "${RKTOPOLOGIESDIR}/time_topology.hpp"
  "${RKTOPOLOGIESDIR}/reachability_space.hpp"
  "${RKTOPOLOGIESDIR}/temporal_space.hpp"
  "${RKTOPOLOGIESDIR}/differentiable_space.hpp"
)

#add_library(reak_topologies STATIC ${TOPOLOGIES_SOURCES})
#setup_custom_target(reak_topologies "${SRCROOT}${RKTOPOLOGIESDIR}")
setup_headers("${TOPOLOGIES_HEADERS}" "${RKTOPOLOGIESDIR}")








