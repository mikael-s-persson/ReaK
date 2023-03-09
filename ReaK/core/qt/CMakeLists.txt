
if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  set( RKCOREQT_UIS
    "rk_class_select.ui" 
    "rk_object_tree.ui" 
    "rk_prop_editor.ui" 
    "data_stream_opt.ui" 
  )
  QT4_WRAP_UI( RKCOREQT_UI_HDRS ${RKCOREQT_UIS} )
  set( RKCOREQT_MOC_HDRS
    "objtree_qtmodel.h"
    "obj_properties_qtmodel.h"
    "rk_object_tree_widget.h"
    "rk_prop_editor_widget.h"
    "data_stream_opt_widget.h"
  )
  QT4_WRAP_CPP( RKCOREQT_MOC_SRCS ${RKCOREQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKCOREQT_SOURCES 
    "objtree_qtmodel.cc"
    "obj_properties_qtmodel.cc"
    "rk_object_tree_widget.cc"
    "rk_prop_editor_widget.cc"
    "data_stream_opt_widget.cc"
  )

  reak_cc_library(
    NAME
      core_qt
    SRCS
      ${RKCOREQT_SOURCES}
      ${RKCOREQT_MOC_SRCS}
      ${RKCOREQT_UI_HDRS}
    DEPS
      ReaK::base
      ReaK::serialization
      ReaK::recorders
      ${QT_LIBRARIES}
    PUBLIC
  )
  
  reak_cc_binary(
    NAME
      rk_object_tree_widget_test
    SRCS
      "rk_object_tree_widget_test.cc"
    DEPS
      ReaK::core_qt
      ReaK::base
      ReaK::serialization
      ReaK::recorders
      ${QT_LIBRARIES}
  )
  
endif()



