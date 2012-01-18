
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

if(NOT WIN32)

  if( OpenCV_FOUND )
    add_executable(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}/test_sampling.cpp")
    setup_custom_target(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}")
    #target_link_libraries(test_adstar ${Boost_LIBRARIES})
    target_link_libraries(test_sampling ${OpenCV_LIBS})
    target_link_libraries(test_sampling reak_core)
  endif()
endif()








