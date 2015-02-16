
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
    "../include/ReaK/core/qt/objtree_qtmodel.hpp"
    "../include/ReaK/core/qt/obj_properties_qtmodel.hpp"
    "../include/ReaK/core/qt/rk_object_tree_widget.hpp"
    "../include/ReaK/core/qt/rk_prop_editor_widget.hpp"
    "../include/ReaK/core/qt/data_stream_opt_widget.hpp"
  )
  QT4_WRAP_CPP( RKCOREQT_MOC_SRCS ${RKCOREQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKCOREQT_SOURCES 
    "objtree_qtmodel.cpp"
    "obj_properties_qtmodel.cpp"
    "rk_object_tree_widget.cpp"
    "rk_prop_editor_widget.cpp"
    "data_stream_opt_widget.cpp"
  )
  
  add_library(reak_core_qt STATIC ${RKCOREQT_SOURCES} ${RKCOREQT_MOC_SRCS} ${RKCOREQT_UI_HDRS})
  ReaK_setup_static_library(reak_core_qt)
  target_link_libraries(reak_core_qt reak_core)
  target_link_libraries(reak_core_qt ${QT_LIBRARIES})
  
#   add_executable(test_objtreewidget "test_objtreewidget.cpp")
#   setup_custom_target(test_objtreewidget)
#   target_link_libraries(test_objtreewidget reak_core_qt reak_core)
#   target_link_libraries(test_objtreewidget ${QT_LIBRARIES})
  
endif()



