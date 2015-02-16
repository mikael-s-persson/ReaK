

if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKCTRLQT_UIS
    "target_predictor_config.ui"
  )
  QT4_WRAP_UI( RKCTRLQT_UI_HDRS ${RKCTRLQT_UIS} )
  set( RKCTRLQT_MOC_HDRS
    "../include/ReaK/control/qt/target_pred_config_widget.hpp"
  )
  QT4_WRAP_CPP( RKCTRLQT_MOC_SRCS ${RKCTRLQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKCTRLQT_SOURCES 
    "target_pred_config_widget.cpp"
  )
  
  add_library(reak_control_qt STATIC ${RKCTRLQT_SOURCES} ${RKCTRLQT_MOC_SRCS} ${RKCTRLQT_UI_HDRS})
  ReaK_setup_static_library(reak_control_qt)
  target_link_libraries(reak_control_qt reak_control reak_core_qt)
  target_link_libraries(reak_control_qt ${QT_LIBRARIES})
  
endif()



