

if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKCTRLQT_UIS
    "target_predictor_config.ui"
  )
  QT4_WRAP_UI( RKCTRLQT_UI_HDRS ${RKCTRLQT_UIS} )
  set( RKCTRLQT_MOC_HDRS
    "target_pred_config_widget.h"
  )
  QT4_WRAP_CPP( RKCTRLQT_MOC_SRCS ${RKCTRLQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKCTRLQT_SOURCES 
    "target_pred_config_widget.cc"
  )
  
  reak_cc_library(
    NAME
      control_qt
    SRCS
      ${RKCTRLQT_SOURCES}
      ${RKCTRLQT_MOC_SRCS}
      ${RKCTRLQT_UI_HDRS}
    DEPS
      ReaK::estimators
      ReaK::core_qt
      ${QT_LIBRARIES}
    PUBLIC
  )
endif()



