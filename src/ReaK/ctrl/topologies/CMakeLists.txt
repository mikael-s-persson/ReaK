
set(TOPOLOGIES_SOURCES 
)

if( OpenCV_FOUND )
  set(TOPOLOGIES_SOURCES 
    ${TOPOLOGIES_SOURCES}
    "${SRCROOT}${RKTOPOLOGIESDIR}/ptrobot2D_test_world.cpp"
  )
endif()


set(TOPOLOGIES_HEADERS 
  "${RKTOPOLOGIESDIR}/basic_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/default_random_sampler.hpp"
  "${RKTOPOLOGIESDIR}/differentiable_space.hpp"
  "${RKTOPOLOGIESDIR}/hyperball_topology.hpp"
  "${RKTOPOLOGIESDIR}/hyperbox_topology.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_limits.hpp"
  "${RKTOPOLOGIESDIR}/joint_space_topologies.hpp"
  "${RKTOPOLOGIESDIR}/line_topology.hpp"
  "${RKTOPOLOGIESDIR}/metric_space_tuple.hpp"
  "${RKTOPOLOGIESDIR}/ptrobot2D_test_world.hpp"
  "${RKTOPOLOGIESDIR}/rate_limited_spaces.hpp"
  "${RKTOPOLOGIESDIR}/reachability_space.hpp"
  "${RKTOPOLOGIESDIR}/se2_topologies.hpp"
  "${RKTOPOLOGIESDIR}/se3_topologies.hpp"
  "${RKTOPOLOGIESDIR}/so3_topologies.hpp"
  "${RKTOPOLOGIESDIR}/temporal_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/temporal_space.hpp"
  "${RKTOPOLOGIESDIR}/time_poisson_topology.hpp"
  "${RKTOPOLOGIESDIR}/time_topology.hpp"
  "${RKTOPOLOGIESDIR}/tuple_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/vect_distance_metrics.hpp"
  "${RKTOPOLOGIESDIR}/vector_topology.hpp"
)


add_library(reak_topologies STATIC ${TOPOLOGIES_SOURCES})
setup_custom_target(reak_topologies "${SRCROOT}${RKTOPOLOGIESDIR}")
setup_headers("${TOPOLOGIES_HEADERS}" "${RKTOPOLOGIESDIR}")
target_link_libraries(reak_topologies reak_core)
if( OpenCV_FOUND )
  target_link_libraries(reak_topologies ${OpenCV_LIBS})
endif()

if(NOT WIN32)

  if( OpenCV_FOUND )
    add_executable(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}/test_sampling.cpp")
    setup_custom_target(test_sampling "${SRCROOT}${RKTOPOLOGIESDIR}")
    #target_link_libraries(test_adstar ${Boost_LIBRARIES})
    target_link_libraries(test_sampling ${OpenCV_LIBS})
    target_link_libraries(test_sampling reak_core)
  endif()
endif()








