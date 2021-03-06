
set(GRAPH_ALG_SOURCES 
)

set(GRAPH_ALG_HEADERS 
  "${RKGRAPHALGDIR}/adaptive_star_neighborhood.hpp"
  "${RKGRAPHALGDIR}/adj_list_tree_overlay.hpp"
  "${RKGRAPHALGDIR}/adstar_search.hpp"
  "${RKGRAPHALGDIR}/any_graph.hpp"
  "${RKGRAPHALGDIR}/anytime_sbastar.hpp"
  "${RKGRAPHALGDIR}/avl_tree.hpp"
  "${RKGRAPHALGDIR}/avl_tree_detail.hpp"
  "${RKGRAPHALGDIR}/bgl_raw_property_graph.hpp"
  "${RKGRAPHALGDIR}/branch_and_bound_connector.hpp"
  "${RKGRAPHALGDIR}/fadprm.hpp"
  "${RKGRAPHALGDIR}/lazy_connector.hpp"
  "${RKGRAPHALGDIR}/lazy_sbastar.hpp"
  "${RKGRAPHALGDIR}/motion_graph_connector.hpp"
  "${RKGRAPHALGDIR}/neighborhood_functors.hpp"
  "${RKGRAPHALGDIR}/node_generators.hpp"
  "${RKGRAPHALGDIR}/prm_connector.hpp"
  "${RKGRAPHALGDIR}/probabilistic_roadmap.hpp"
  "${RKGRAPHALGDIR}/rr_graph.hpp"
  "${RKGRAPHALGDIR}/rr_tree.hpp"
  "${RKGRAPHALGDIR}/rrt_star.hpp"
  "${RKGRAPHALGDIR}/sbastar_rrtstar.hpp"
  "${RKGRAPHALGDIR}/sbastar_search.hpp"
  "${RKGRAPHALGDIR}/sbmp_visitor_concepts.hpp"
  "${RKGRAPHALGDIR}/tree_organizer_concept.hpp"
)

#add_library(reak_graph_alg STATIC ${GRAPH_ALG_SOURCES})
#setup_custom_target(reak_graph_alg "${SRCROOT}${RKGRAPHALGDIR}")
setup_headers("${GRAPH_ALG_HEADERS}" "${RKGRAPHALGDIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(test_graphs "${SRCROOT}${RKGRAPHALGDIR}/test_graphs.cpp")
setup_custom_target(test_graphs "${SRCROOT}${RKGRAPHALGDIR}")
target_link_libraries(test_graphs reak_topologies reak_core)

add_executable(unit_test_assoc_containers "${SRCROOT}${RKGRAPHALGDIR}/unit_test_assoc_containers.cpp")
setup_custom_test_program(unit_test_assoc_containers "${SRCROOT}${RKGRAPHALGDIR}")
target_link_libraries(unit_test_assoc_containers reak_core)







