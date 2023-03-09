

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKMBDQT_UIS
    "kte_model_viewer.ui" 
    "chaser_target_mdl_config.ui"
    "chaser_target_interact.ui"
  )
  QT4_WRAP_UI( RKMBDQT_UI_HDRS ${RKMBDQT_UIS} )
  set( RKMBDQT_MOC_HDRS
    "kte_model_viewer.h"
    "rk_view3d_menu.h"
    "chaser_target_config_widget.h"
    "chaser_target_interact_widget.h"
  )
  QT4_WRAP_CPP( RKMBDQT_MOC_SRCS ${RKMBDQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKMBDQT_SOURCES 
    "kte_model_viewer.cc"
    "kte_model_viewer_loader.cc"
    "rk_view3d_menu.cc"
    "chaser_target_config_widget.cc"
    "chaser_target_interact_widget.cc"
  )
  
  reak_cc_library(
    NAME
      mbd_qt
    SRCS
      ${RKMBDQT_SOURCES}
      ${RKMBDQT_MOC_SRCS}
      ${RKMBDQT_UI_HDRS}
    DEPS
      ReaK::core_qt
      ReaK::mbd_coin
      ${QT_LIBRARIES}
    PUBLIC
  )
  
endif()



