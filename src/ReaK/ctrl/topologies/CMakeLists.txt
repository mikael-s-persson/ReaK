
set(TOPOLOGIES_SOURCES 
)

set(TOPOLOGIES_HEADERS 
  "${RKTOPOLOGIESDIR}/line_topology.hpp"
  "${RKTOPOLOGIESDIR}/time_topology.hpp"
  "${RKTOPOLOGIESDIR}/time_poisson_topology.hpp"
  "${RKTOPOLOGIESDIR}/reachability_space.hpp"
  "${RKTOPOLOGIESDIR}/temporal_space.hpp"
  "${RKTOPOLOGIESDIR}/basic_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/metric_space_tuple.hpp"
  "${RKTOPOLOGIESDIR}/tuple_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/differentiable_space.hpp"
  "${RKTOPOLOGIESDIR}/vector_topology.hpp"
  "${RKTOPOLOGIESDIR}/vect_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/hyperball_topology.hpp"
  "${RKTOPOLOGIESDIR}/hyperbox_topology.hpp"
  "${RKTOPOLOGIESDIR}/so3_topologies.hpp"
)

#add_library(reak_topologies STATIC ${TOPOLOGIES_SOURCES})
#setup_custom_target(reak_topologies "${SRCROOT}${RKTOPOLOGIESDIR}")
setup_headers("${TOPOLOGIES_HEADERS}" "${RKTOPOLOGIESDIR}")








