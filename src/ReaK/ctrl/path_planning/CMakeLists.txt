
set(PATH_PLANNING_SOURCES 
)

set(PATH_PLANNING_HEADERS 
  "${RKPATHPLANNINGDIR}/bounded_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/dvp_layout_adjacency_list.hpp"
  "${RKPATHPLANNINGDIR}/dvp_tree_detail.hpp"
  "${RKPATHPLANNINGDIR}/global_rng.hpp"
  "${RKPATHPLANNINGDIR}/interpolator_concept.hpp"
  "${RKPATHPLANNINGDIR}/manipulator_topo_maps.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/metric_space_search.hpp"
  "${RKPATHPLANNINGDIR}/motion_planner_base.hpp"
  "${RKPATHPLANNINGDIR}/predicted_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/prob_distribution_concept.hpp"
  "${RKPATHPLANNINGDIR}/prm_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/random_sampler_concept.hpp"
  "${RKPATHPLANNINGDIR}/reachability_sort.hpp"
  "${RKPATHPLANNINGDIR}/reachability_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/rrt_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/rrtstar_path_planner.hpp"
  "${RKPATHPLANNINGDIR}/sbmp_reporter_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_path_concept.hpp"
  "${RKPATHPLANNINGDIR}/spatial_trajectory_concept.hpp"
  "${RKPATHPLANNINGDIR}/subspace_concept.hpp"
  "${RKPATHPLANNINGDIR}/tangent_bundle_concept.hpp"
  "${RKPATHPLANNINGDIR}/temporal_space_concept.hpp"
  "${RKPATHPLANNINGDIR}/topological_map_concepts.hpp"
  "${RKPATHPLANNINGDIR}/topological_search.hpp"
  "${RKPATHPLANNINGDIR}/trajectory_base.hpp"
  "${RKPATHPLANNINGDIR}/transformed_trajectory.hpp"
)

#add_library(reak_path_planning STATIC ${PATH_PLANNING_SOURCES})
#setup_custom_target(reak_path_planning "${SRCROOT}${RKPATHPLANNINGDIR}")
setup_headers("${PATH_PLANNING_HEADERS}" "${RKPATHPLANNINGDIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_executable(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}/test_vp_tree.cpp")
setup_custom_target(test_vp_tree "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_vp_tree reak_core)

add_executable(test_dvp_adj_list "${SRCROOT}${RKPATHPLANNINGDIR}/test_dvp_adj_list.cpp")
setup_custom_target(test_dvp_adj_list "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_dvp_adj_list reak_core)

add_executable(test_planners "${SRCROOT}${RKPATHPLANNINGDIR}/test_planners.cpp")
setup_custom_target(test_planners "${SRCROOT}${RKPATHPLANNINGDIR}")
target_link_libraries(test_planners reak_topologies reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")

find_package( Qt4 )
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
  target_link_libraries(pptester ${QT_LIBRARIES})
  target_link_libraries(pptester ${Boost_LIBRARIES})
  target_link_libraries(pptester ${OpenCV_LIBS})
  target_link_libraries(pptester reak_core)
  setup_custom_target(pptester "${SRCROOT}${RKPATHPLANNINGDIR}")
endif()







