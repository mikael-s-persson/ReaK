

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKMBDQT_UIS
#     "kte_model_viewer.ui" 
    "chaser_target_mdl_config.ui"
    "chaser_target_interact.ui"
  )
  QT4_WRAP_UI( RKMBDQT_UI_HDRS ${RKMBDQT_UIS} )
  set( RKMBDQT_MOC_HDRS
    "../include/ReaK/mbd/qt/rk_view3d_menu.hpp"
    "../include/ReaK/mbd/qt/chaser_target_config_widget.hpp"
    "../include/ReaK/mbd/qt/chaser_target_interact_widget.hpp"
  )
  QT4_WRAP_CPP( RKMBDQT_MOC_SRCS ${RKMBDQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKMBDQT_SOURCES 
    "rk_view3d_menu.cpp"
    "chaser_target_config_widget.cpp"
    "chaser_target_interact_widget.cpp"
  )
  
  add_library(reak_mbd_qt STATIC ${RKMBDQT_SOURCES} ${RKMBDQT_MOC_SRCS} ${RKMBDQT_UI_HDRS})
  ReaK_setup_static_library(reak_mbd_qt)
  target_link_libraries(reak_mbd_qt reak_core_qt reak_mbd_coin)
  target_link_libraries(reak_mbd_qt ${QT_LIBRARIES})
  
endif()



