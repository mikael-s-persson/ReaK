
set(PATH_PLANNING_SOURCES 
)

set(PATH_PLANNING_HEADERS 
  "${RKPATHPLANNINGDIR}/any_knn_synchro.hpp"
  "${RKPATHPLANNINGDIR}/any_motion_graphs.hpp"
  "${RKPATHPLANNINGDIR}/any_sbmp_reporter.hpp"
  "${RKPATHPLANNINGDIR}/basic_sbmp_reporters.hpp"
  "${RKPATHPLANNINGDIR}/bounded_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/density_calculators.hpp"
  "${RKPATHPLANNINGDIR}/density_plan_visitors.hpp"
  "${RKPATHPLANNINGDIR}/dvp_layout_adjacency_list.hpp"
  "${RKPATHPLANNINGDIR}/dvp_tree_detail.hpp"
  "${RKPATHPLANNINGDIR}/extrapolator_concept.hpp"
  "${RKPATHPLANNINGDIR}/fadprm_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/fadprm_path_planner.tpp"
  "${RKPATHPLANNINGDIR}/fadprm_manip_planners.hpp"
  "${RKPATHPLANNINGDIR}/frame_tracer_coin3d.hpp"
  "${RKPATHPLANNINGDIR}/interpolator_concept.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_search.hpp"
  "${RKPATHPLANNINGDIR}/motion_graph_structures.hpp"
  "${RKPATHPLANNINGDIR}/motion_planner_base.hpp"
  "${RKPATHPLANNINGDIR}/path_base.hpp"
  "${RKPATHPLANNINGDIR}/path_planner_options.hpp"
  "${RKPATHPLANNINGDIR}/path_wrapper.hpp"
  "${RKPATHPLANNINGDIR}/planning_queries.hpp"
  "${RKPATHPLANNINGDIR}/planning_visitors.hpp"
  "${RKPATHPLANNINGDIR}/predicted_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/prm_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/prm_path_planner.tpp"
  "${RKPATHPLANNINGDIR}/prm_manip_planners.hpp"
  "${RKPATHPLANNINGDIR}/prob_distribution_concept.hpp"
  "${RKPATHPLANNINGDIR}/random_sampler_concept.hpp"
  "${RKPATHPLANNINGDIR}/reachability_sort.hpp"
  "${RKPATHPLANNINGDIR}/reachability_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/rrt_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/rrt_path_planner.tpp"
  "${RKPATHPLANNINGDIR}/rrt_manip_planners.hpp"
  "${RKPATHPLANNINGDIR}/rrtstar_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/rrtstar_path_planner.tpp"
  "${RKPATHPLANNINGDIR}/rrtstar_manip_planners.hpp"
  "${RKPATHPLANNINGDIR}/sbastar_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/sbastar_path_planner.tpp"
  "${RKPATHPLANNINGDIR}/sbastar_manip_planners.hpp"
  "${RKPATHPLANNINGDIR}/sbmp_point_recorder.hpp"
  "${RKPATHPLANNINGDIR}/sbmp_reporter_concept.hpp"
  "${RKPATHPLANNINGDIR}/seq_path_base.hpp"
  "${RKPATHPLANNINGDIR}/seq_path_wrapper.hpp"
  "${RKPATHPLANNINGDIR}/seq_trajectory_base.hpp"
  "${RKPATHPLANNINGDIR}/seq_trajectory_wrapper.hpp"
  "${RKPATHPLANNINGDIR}/sequential_path_concept.hpp"
  "${RKPATHPLANNINGDIR}/sequential_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_path_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/steerable_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/subspace_concept.hpp"
  "${RKPATHPLANNINGDIR}/tangent_bundle_concept.hpp"
  "${RKPATHPLANNINGDIR}/temporal_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/topological_map_concepts.hpp"
  "${RKPATHPLANNINGDIR}/topological_search.hpp"
  "${RKPATHPLANNINGDIR}/trajectory_base.hpp"
  "${RKPATHPLANNINGDIR}/transformed_trajectory.hpp"
  "${RKPATHPLANNINGDIR}/vlist_sbmp_report.hpp"
)



set(PLANNERS_SOURCES 
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/rrt_manip_planners_static.cpp"
#   "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/rrt_manip_planners_dynamic.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/rrtstar_manip_planners_static.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/rrtstar_manip_planners_dynamic.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/sbastar_manip_planners_static.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/sbastar_manip_planners_dynamic.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/prm_manip_planners_static.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/prm_manip_planners_dynamic.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/fadprm_manip_planners_static.cpp"
  "${SRCROOT}${RKPATHPLANNINGDIR}/extern_templates/fadprm_manip_planners_dynamic.cpp"
)


# foreach(PLANNERS_DIM_COUNTER 0 1 2 3 4 5 6 7 8 9 10)
foreach(PLANNERS_DIM_COUNTER 7)
  add_library(reakobj_planners_${PLANNERS_DIM_COUNTER} OBJECT ${PLANNERS_SOURCES})
  set_target_properties(reakobj_planners_${PLANNERS_DIM_COUNTER} PROPERTIES COMPILE_DEFINITIONS "RK_PLANNERS_DIM=${PLANNERS_DIM_COUNTER}")
  set(PLANNERS_OBJECT_FILES ${PLANNERS_OBJECT_FILES} $<TARGET_OBJECTS:reakobj_planners_${PLANNERS_DIM_COUNTER}>)
endforeach(PLANNERS_DIM_COUNTER)

add_library(reak_planners STATIC ${PLANNERS_OBJECT_FILES})
setup_custom_target(reak_planners "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(reak_planners reak_topologies reak_interp reak_core)



setup_headers("${PATH_PLANNING_HEADERS}" "${RKPATHPLANNINGDIR}")


add_executable(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}/test_vp_tree.cpp")
setup_custom_target(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_vp_tree reak_topologies reak_core)

add_executable(test_dvp_adj_list "${SRCROOT}${RKPATHPLANNINGDIR}/test_dvp_adj_list.cpp")
setup_custom_target(test_dvp_adj_list "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_dvp_adj_list reak_topologies reak_core)

add_executable(test_planners "${SRCROOT}${RKPATHPLANNINGDIR}/test_planners.cpp")
setup_custom_target(test_planners "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_planners reak_topologies reak_core)


add_custom_target(test_hidim_planners)
foreach(HIDIM_COUNTER 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
  add_executable(test_hidim_planners_e${HIDIM_COUNTER} EXCLUDE_FROM_ALL "${SRCROOT}${RKPATHPLANNINGDIR}/test_hidim_planners.cpp")
  setup_custom_target(test_hidim_planners_e${HIDIM_COUNTER} "${SRCROOT}${RKPATHPLANNINGDIR}")
  set_target_properties(test_hidim_planners_e${HIDIM_COUNTER} PROPERTIES COMPILE_DEFINITIONS "RK_HIDIM_PLANNER_N=${HIDIM_COUNTER}")
  target_link_libraries(test_hidim_planners_e${HIDIM_COUNTER} reak_topologies reak_core)
  add_dependencies(test_hidim_planners test_hidim_planners_e${HIDIM_COUNTER})
endforeach(HIDIM_COUNTER)

if()
#if(${QT_FOUND})
  include( ${QT_USE_FILE} )
  set( PATHPLANNING_PPTESTER_UIS
     "${SRCROOT}${RKPATHPLANNINGDIR}/qt_ui/pptester.ui" 
     "${SRCROOT}${RKPATHPLANNINGDIR}/qt_ui/ppresultview.ui"
     "${SRCROOT}${RKPATHPLANNINGDIR}/qt_ui/rrtproperties.ui"
  )
  QT4_WRAP_UI( PATHPLANNING_PPTESTER_UI_HDRS ${PATHPLANNING_PPTESTER_UIS} )
  set( PATHPLANNING_PPTESTER_MOC_HDRS
     "${SRCROOT}${RKPATHPLANNINGDIR}/qt_ui/pptester_impl.h"
  )
  QT4_WRAP_CPP( PATHPLANNING_PPTESTER_MOC_SRCS ${PATHPLANNING_PPTESTER_MOC_HDRS} )
  include_directories( AFTER ${CMAKE_BINARY_DIR}${RKPATHPLANNINGDIR} )
  SET (CMAKE_MODULE_PATH "/usr/share/cmake/Modules")

  INCLUDE (UsePlayerC++)

  include_directories(${PLAYERCPP_INCLUDE_DIRS})
  message(STATUS "PlayerC++ was found at '${PLAYERCPP_INCLUDE_DIRS}'")
  PLAYER_ADD_PLAYERCPP_CLIENT(pptester SOURCES "${SRCROOT}${RKPATHPLANNINGDIR}/qt_ui/pptester_impl.cpp" ${PATHPLANNING_PPTESTER_MOC_SRCS} ${PATHPLANNING_PPTESTER_UI_HDRS} LINKFLAGS ${replaceLib})
  target_link_libraries(pptester ${OpenCV_LIBS} ${QT_LIBRARIES})
  target_link_libraries(pptester reak_topologies reak_core)
  setup_custom_target(pptester "${SRCROOT}${RKPATHPLANNINGDIR}")
endif()








