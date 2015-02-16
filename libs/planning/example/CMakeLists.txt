
add_executable(test_graphs "test_graphs.cpp")
ReaK_setup_target(test_graphs)
target_link_libraries(test_graphs reak_topologies)

add_executable(test_vp_tree "test_vp_tree.cpp")
ReaK_setup_target(test_vp_tree)
target_link_libraries(test_vp_tree reak_topologies)

add_executable(test_dvp_adj_list "test_dvp_adj_list.cpp")
ReaK_setup_target(test_dvp_adj_list)
target_link_libraries(test_dvp_adj_list reak_topologies)

add_executable(test_planners "test_planners.cpp")
ReaK_setup_target(test_planners)
target_link_libraries(test_planners reak_topologies)

add_custom_target(test_hidim_planners)
foreach(HIDIM_COUNTER 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)
  add_executable(test_hidim_planners_e${HIDIM_COUNTER} EXCLUDE_FROM_ALL "test_hidim_planners.cpp")
  ReaK_setup_target(test_hidim_planners_e${HIDIM_COUNTER})
  set_target_properties(test_hidim_planners_e${HIDIM_COUNTER} PROPERTIES COMPILE_DEFINITIONS "RK_HIDIM_PLANNER_N=${HIDIM_COUNTER}")
  target_link_libraries(test_hidim_planners_e${HIDIM_COUNTER} reak_planning reak_topologies)
  add_dependencies(test_hidim_planners test_hidim_planners_e${HIDIM_COUNTER})
endforeach(HIDIM_COUNTER)

add_executable(run_CRS_planner "run_CRS_planner.cpp")
ReaK_setup_target(run_CRS_planner)
target_link_libraries(run_CRS_planner reak_planning reak_topologies)

add_executable(test_MEAQR_knn_methods "test_MEAQR_knn_methods.cpp")
ReaK_setup_target(test_MEAQR_knn_methods)
target_link_libraries(test_MEAQR_knn_methods reak_planning reak_control reak_topologies reak_mbd reak_core)

add_executable(X8_run_planner "X8_run_planner.cpp")
ReaK_setup_target(X8_run_planner)
target_link_libraries(X8_run_planner reak_planning reak_control reak_topologies)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  add_executable(X8_test_scene "X8_test_scene.cpp")
  ReaK_setup_target(X8_test_scene)
  target_link_libraries(X8_test_scene reak_planning reak_control reak_topologies reak_mbd_coin)
  target_link_libraries(X8_test_scene ${COIN3D_LIBRARIES} ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  
  set( PLANNER_3DVIEWER_UIS
    "planner_3Dview.ui" 
    "planner_space_config.ui"
    "../qt_ui/planner_alg_config.ui"
  )
  QT4_WRAP_UI( PLANNER_3DVIEWER_UI_HDRS ${PLANNER_3DVIEWER_UIS} )
  set( PLANNER_3DVIEWER_MOC_HDRS
    "planner_3Dviewer.hpp"
  )
  QT4_WRAP_CPP( PLANNER_3DVIEWER_MOC_SRCS ${PLANNER_3DVIEWER_MOC_HDRS} )

  set(PLANNER_3DVIEWER_SOURCES 
    "planner_3Dviewer.cpp"
  )
  
  add_executable(planner_3Dviewer ${PLANNER_3DVIEWER_SOURCES} ${PLANNER_3DVIEWER_MOC_SRCS} ${PLANNER_3DVIEWER_UI_HDRS})
  ReaK_setup_target(planner_3Dviewer)
  target_link_libraries(planner_3Dviewer reak_planning_qt reak_planning reak_mbd_coin)
  target_link_libraries(planner_3Dviewer ${COIN3D_LIBRARIES} ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  
  set( ROBOTAIRSHIP_CRSPLANNER_UIS
    "CRS_planner_window.ui" 
    "CRS_planner_config.ui"
    "CRS_run_dialog.ui"
  )
  QT4_WRAP_UI( ROBOTAIRSHIP_CRSPLANNER_UI_HDRS ${ROBOTAIRSHIP_CRSPLANNER_UIS} )
  set( ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS
    "CRS_run_dialog.hpp"
    "CRS_planner_impl.hpp"
  )
  QT4_WRAP_CPP( ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS ${ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS} )
  
  set( ROBOTAIRSHIP_CRSPLANNER_SRCS
    "CRS_run_dialog.cpp"
    "CRS_planner_impl.cpp"
    "CRS_planner_exec.cpp"
    "CRS_planner_dynexec.cpp"
  )
  
  add_executable(CRS_planner ${ROBOTAIRSHIP_CRSPLANNER_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_UI_HDRS})
  ReaK_setup_target(CRS_planner)
  target_link_libraries(CRS_planner reak_planning_qt reak_control_qt reak_planning reak_mbd_coin)
  target_link_libraries(CRS_planner ${COIN3D_LIBRARIES} ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()


if()
#if(${QT_FOUND})
  include( ${QT_USE_FILE} )
  set( PATHPLANNING_PPTESTER_UIS
     "pptester.ui" 
     "ppresultview.ui"
     "rrtproperties.ui"
  )
  QT4_WRAP_UI( PATHPLANNING_PPTESTER_UI_HDRS ${PATHPLANNING_PPTESTER_UIS} )
  set( PATHPLANNING_PPTESTER_MOC_HDRS
     "pptester_impl.h"
  )
  QT4_WRAP_CPP( PATHPLANNING_PPTESTER_MOC_SRCS ${PATHPLANNING_PPTESTER_MOC_HDRS} )
  include_directories( AFTER ${CMAKE_BINARY_DIR}${RKPATHPLANNINGDIR} )
  SET (CMAKE_MODULE_PATH "/usr/share/cmake/Modules")

  INCLUDE (UsePlayerC++)

  include_directories(${PLAYERCPP_INCLUDE_DIRS})
  message(STATUS "PlayerC++ was found at '${PLAYERCPP_INCLUDE_DIRS}'")
  PLAYER_ADD_PLAYERCPP_CLIENT(pptester SOURCES "pptester_impl.cpp" ${PATHPLANNING_PPTESTER_MOC_SRCS} ${PATHPLANNING_PPTESTER_UI_HDRS} LINKFLAGS ${replaceLib})
  target_link_libraries(pptester ${OpenCV_LIBS} ${QT_LIBRARIES})
  target_link_libraries(pptester reak_topologies reak_core)
  setup_custom_target(pptester)
endif()




