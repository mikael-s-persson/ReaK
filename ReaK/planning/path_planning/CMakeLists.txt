
set(MANIP_PLANNERS_SOURCES 
  "fadprm_manip_planners_dynamic.cc"
  "fadprm_manip_planners_static.cc"
  "prm_manip_planners_dynamic.cc"
  "prm_manip_planners_static.cc"
#  "rrt_manip_planners_dynamic.cc"
  "rrt_manip_planners_static.cc"
  "rrtstar_manip_planners_dynamic.cc"
  "rrtstar_manip_planners_static.cc"
  "sbastar_manip_planners_dynamic.cc"
  "sbastar_manip_planners_static.cc"
)

# foreach(PLANNERS_DIM_COUNTER 0 1 2 3 4 5 6 7 8 9 10)
foreach(PLANNERS_DIM_COUNTER 7)
  add_library(reakobj_manip_planners_${PLANNERS_DIM_COUNTER} OBJECT ${MANIP_PLANNERS_SOURCES})
  set_target_properties(reakobj_manip_planners_${PLANNERS_DIM_COUNTER} PROPERTIES COMPILE_DEFINITIONS "RK_PLANNERS_DIM=${PLANNERS_DIM_COUNTER}")
  set(MANIP_PLANNERS_OBJECT_FILES ${MANIP_PLANNERS_OBJECT_FILES} $<TARGET_OBJECTS:reakobj_manip_planners_${PLANNERS_DIM_COUNTER}>)
endforeach(PLANNERS_DIM_COUNTER)

reak_cc_library(
  NAME
    path_planning
  SRCS
    ${MANIP_PLANNERS_OBJECT_FILES}
  HDRS
    "any_knn_synchro.h"
    "any_motion_graphs.h"
    "any_sbmp_reporter.h"
    "basic_sbmp_reporters.h"
    "density_calculators.h"
    "density_plan_visitors.h"
    "dvp_layout_adjacency_list.h"
    "dvp_tree_detail.h"
    "fadprm_manip_planners.h"
    "fadprm_path_planner.h"
    "frame_tracer_coin3d.h"
    "intercept_query.h"
    "MEAQR_rrtstar_planner.h"
    "MEAQR_sbastar_planner.h"
    "metric_space_search.h"
    "motion_graph_structures.h"
    "motion_planner_base.h"
    "multi_dvp_tree_search.h"
    "p2p_planning_query.h"
    "planner_exec_engines.h"
    "planner_exec_intercept.h"
    "planning_queries.h"
    "planning_visitors.h"
    "prm_manip_planners.h"
    "prm_path_planner.h"
    "reachability_sort.h"
    "rrt_manip_planners.h"
    "rrt_path_planner.h"
    "rrtstar_manip_planners.h"
    "rrtstar_path_planner.h"
    "sbastar_manip_planners.h"
    "sbastar_path_planner.h"
    "sbmp_point_recorder.h"
    "sbmp_reporter_concept.h"
    "solution_path_factories.h"
    "topological_search.h"
    "vlist_sbmp_report.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::kte
    ReaK::models
    ReaK::spaces
    ReaK::interpolation
    ReaK::graph_alg
  PUBLIC
)

reak_cc_library(
  NAME
    path_planning_po
  SRCS
    "path_planner_options_po.cc"
    "planning_space_options_po.cc"
  HDRS
    "path_planner_options.h"
    "path_planner_options_po.h"
    "planning_space_options.h"
    "planning_space_options_po.h"
  DEPS
    absl::flags
    ReaK::base
    ReaK::serialization
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::kte
    ReaK::models
    ReaK::spaces
    ReaK::interpolation
    ReaK::path_planning
  PUBLIC
)

reak_cc_binary(
  NAME
    dvp_tree_test
  SRCS
    "dvp_tree_test.cc"
  DEPS
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::spaces
    ReaK::interpolation
    ReaK::path_planning
)

reak_cc_binary(
  NAME
    dvp_adj_list_test
  SRCS
    "dvp_adj_list_test.cc"
  DEPS
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::spaces
    ReaK::interpolation
    ReaK::path_planning
)

reak_cc_binary(
  NAME
    planners_test
  SRCS
    "planners_test.cc"
  DEPS
    absl::flags
    absl::flags_parse
    ReaK::base
    ReaK::serialization
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::spaces
    ReaK::interpolation
    ReaK::path_planning
    ReaK::path_planning_po
)

add_custom_target(planners_hidim_test)
foreach(HIDIM_COUNTER 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
  reak_cc_binary(
    NAME
      planners_hidim_test_e${HIDIM_COUNTER}
    SRCS
      "planners_hidim_test.cc"
    DEPS
      absl::flags
      absl::flags_parse
      ReaK::base
      ReaK::serialization
      ReaK::lin_alg
      ReaK::kinetostatics
      ReaK::spaces
      ReaK::interpolation
      ReaK::path_planning
      ReaK::path_planning_po
    DEFINES
      "RK_HIDIM_PLANNER_N=${HIDIM_COUNTER}"
  )
  add_dependencies(planners_hidim_test planners_hidim_test_e${HIDIM_COUNTER})
endforeach(HIDIM_COUNTER)

add_subdirectory(tools)
