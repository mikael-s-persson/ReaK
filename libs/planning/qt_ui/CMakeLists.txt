

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKPLANNINGQT_UIS
    "manip_space_config.ui"
    "planner_alg_config.ui"
  )
  QT4_WRAP_UI( RKPLANNINGQT_UI_HDRS ${RKPLANNINGQT_UIS} )
  set( RKPLANNINGQT_MOC_HDRS
    "../include/ReaK/planning/qt/manip_space_config_widget.hpp"
    "../include/ReaK/planning/qt/planner_alg_config_widget.hpp"
  )
  QT4_WRAP_CPP( RKPLANNINGQT_MOC_SRCS ${RKPLANNINGQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKPLANNINGQT_SOURCES 
    "manip_space_config_widget.cpp"
    "planner_alg_config_widget.cpp"
  )
  
  add_library(reak_planning_qt STATIC ${RKPLANNINGQT_SOURCES} ${RKPLANNINGQT_MOC_SRCS} ${RKPLANNINGQT_UI_HDRS})
  ReaK_setup_static_library(reak_planning_qt)
  target_link_libraries(reak_planning_qt reak_planning reak_mbd_qt)
  target_link_libraries(reak_planning_qt ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()



