

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKPLANNINGQT_UIS
    "manip_space_config.ui"
    "planner_alg_config.ui"
  )
  QT4_WRAP_UI( RKPLANNINGQT_UI_HDRS ${RKPLANNINGQT_UIS} )
  set( RKPLANNINGQT_MOC_HDRS
    "manip_space_config_widget.h"
    "planner_alg_config_widget.h"
  )
  QT4_WRAP_CPP( RKPLANNINGQT_MOC_SRCS ${RKPLANNINGQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKPLANNINGQT_SOURCES 
    "manip_space_config_widget.cc"
    "planner_alg_config_widget.cc"
  )
  
  reak_cc_library(
    NAME
      planning_qt
    SRCS
      ${RKPLANNINGQT_SOURCES}
      ${RKPLANNINGQT_MOC_SRCS}
      ${RKPLANNINGQT_UI_HDRS}
    HDRS
      ReaK::path_planning
      ReaK::mbd_qt
      ${SOQT4_LIBRARIES}
      ${QT_LIBRARIES}
    PUBLIC
  )
  
  reak_cc_binary(
    NAME
      X8_test_scene
    SRCS
      "X8_test_scene.cc"
    HDRS
      ReaK::base
      ReaK::lin_alg
      ReaK::kinetostatics
      ReaK::proximity
      ReaK::path_planning
      ReaK::controllers
      ReaK::systems
      ReaK::spaces
      ReaK::interpolation
      ReaK::mbd_coin
      ${COIN3D_LIBRARIES}
      ${SOQT4_LIBRARIES}
      ${QT_LIBRARIES}
  )
  
  
  set( PLANNER_3DVIEWER_UIS
    "planner_3Dview.ui" 
    "planner_space_config.ui"
    "planner_alg_config.ui"
  )
  QT4_WRAP_UI( PLANNER_3DVIEWER_UI_HDRS ${PLANNER_3DVIEWER_UIS} )
  set( PLANNER_3DVIEWER_MOC_HDRS
    "planner_3Dviewer.h"
  )
  QT4_WRAP_CPP( PLANNER_3DVIEWER_MOC_SRCS ${PLANNER_3DVIEWER_MOC_HDRS} )

  set(PLANNER_3DVIEWER_SOURCES 
    "planner_3Dviewer.cc"
  )
  
  reak_cc_binary(
    NAME
      planner_3Dviewer
    SRCS
      ${PLANNER_3DVIEWER_SOURCES}
      ${PLANNER_3DVIEWER_MOC_SRCS}
      ${PLANNER_3DVIEWER_UI_HDRS}
    HDRS
      ReaK::base
      ReaK::lin_alg
      ReaK::kinetostatics
      ReaK::shapes
      ReaK::proximity
      ReaK::planning_qt
      ReaK::path_planning
      ReaK::mbd_coin
      ${COIN3D_LIBRARIES}
      ${SOQT4_LIBRARIES}
      ${QT_LIBRARIES}
  )
  
  
  set( ROBOTAIRSHIP_CRSPLANNER_UIS
    "CRS_planner_window.ui" 
    "CRS_planner_config.ui"
    "CRS_run_dialog.ui"
  )
  QT4_WRAP_UI( ROBOTAIRSHIP_CRSPLANNER_UI_HDRS ${ROBOTAIRSHIP_CRSPLANNER_UIS} )
  set( ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS
    "CRS_run_dialog.h"
    "CRS_planner_impl.h"
  )
  QT4_WRAP_CPP( ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS ${ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS} )
  
  set( ROBOTAIRSHIP_CRSPLANNER_SRCS
    "CRS_run_dialog.cc"
    "CRS_planner_impl.cc"
    "CRS_planner_exec.cc"
    "CRS_planner_dynexec.cc"
  )
  
  reak_cc_binary(
    NAME
      CRS_planner
    SRCS
      ${ROBOTAIRSHIP_CRSPLANNER_SRCS}
      ${ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS}
      ${ROBOTAIRSHIP_CRSPLANNER_UI_HDRS}
    HDRS
      absl::flags
      absl::flags_parse
      ReaK::base
      ReaK::lin_alg
      ReaK::kinetostatics
      ReaK::shapes
      ReaK::proximity
      ReaK::planning_qt
      ReaK::control_qt
      ReaK::systems
      ReaK::path_planning
      ReaK::mbd_coin
      ${COIN3D_LIBRARIES}
      ${SOQT4_LIBRARIES}
      ${QT_LIBRARIES}
  )
  
  if()
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

    include(UsePlayerC++)

    include_directories(${PLAYERCPP_INCLUDE_DIRS})
    message(STATUS "PlayerC++ was found at '${PLAYERCPP_INCLUDE_DIRS}'")
    PLAYER_ADD_PLAYERCPP_CLIENT(pptester SOURCES "pptester_impl.cc" ${PATHPLANNING_PPTESTER_MOC_SRCS} ${PATHPLANNING_PPTESTER_UI_HDRS} LINKFLAGS ${replaceLib})
    target_link_libraries(pptester ${OpenCV_LIBS} ${QT_LIBRARIES})
    target_link_libraries(pptester ReaK::spaces ReaK::interpolation ReaK::path_planning ReaK::base)
    setup_custom_target(pptester)
  endif()
  
endif()



