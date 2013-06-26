
if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_BINARY_DIR}${RKGUIQTDIR}" )
  
#   add_library(reak_gui_qt STATIC ${RKGUIQT_SOURCES})
#   setup_custom_target(reak_gui_qt "${SRCROOT}${RKGUIQTDIR}")
#   setup_headers("${RKGUIQT_HEADERS}" "${RKGUIQTDIR}")
#   target_link_libraries(reak_gui_qt reak_topologies reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
#   target_link_libraries(reak_gui_qt ${QT_LIBRARIES})
  
#   add_executable(test_CRS_planning "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_CRS_planning.cpp")
#   setup_custom_target(test_CRS_planning "${SRCROOT}${RKROBOTAIRSHIPDIR}")
#   target_link_libraries(test_CRS_planning reak_topologies reak_kte_coin reak_geom_coin reak_robot_airship reak_mbd_kte reak_geom_prox reak_geom reak_core)
#   target_link_libraries(test_CRS_planning ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
  set( RKGUIQT_UIS
     "${SRCROOT}${RKGUIQTDIR}/ui/rk_class_select.ui" 
     "${SRCROOT}${RKGUIQTDIR}/ui/rk_object_tree.ui" 
     "${SRCROOT}${RKGUIQTDIR}/ui/rk_prop_editor.ui" 
  )
  QT4_WRAP_UI( RKGUIQT_UI_HDRS ${RKGUIQT_UIS} )
  set( RKGUIQT_MOC_HDRS
     "${SRCROOT}${RKGUIQTDIR}/objtree_qtmodel.hpp"
     "${SRCROOT}${RKGUIQTDIR}/obj_properties_qtmodel.hpp"
     "${SRCROOT}${RKGUIQTDIR}/rk_object_tree_widget.hpp"
     "${SRCROOT}${RKGUIQTDIR}/rk_prop_editor_widget.hpp"
  )
  QT4_WRAP_CPP( RKGUIQT_MOC_SRCS ${RKGUIQT_MOC_HDRS} OPTIONS "-DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED")
  
  set(RKGUIQT_SOURCES 
    "${SRCROOT}${RKGUIQTDIR}/objtree_qtmodel.cpp"
    "${SRCROOT}${RKGUIQTDIR}/obj_properties_qtmodel.cpp"
    "${SRCROOT}${RKGUIQTDIR}/rk_object_tree_widget.cpp"
    "${SRCROOT}${RKGUIQTDIR}/rk_prop_editor_widget.cpp"
  )
  
  set(RKGUIQT_HEADERS 
#    "${RKGUIQTDIR}/objtree_qtmodel.hpp"
  )
    
  
  add_library(reak_gui_qt STATIC ${RKGUIQT_SOURCES} ${RKGUIQT_MOC_SRCS} ${RKGUIQT_UI_HDRS})
  setup_custom_target(reak_gui_qt "${SRCROOT}${RKGUIQTDIR}")
  setup_headers("${RKGUIQT_HEADERS}" "${RKGUIQTDIR}")
  target_link_libraries(reak_gui_qt reak_topologies reak_mbd_kte reak_core)
  target_link_libraries(reak_gui_qt ${QT_LIBRARIES})
  target_link_libraries(reak_gui_qt ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})
  
  add_executable(test_objtreewidget "${SRCROOT}${RKGUIQTDIR}/test_objtreewidget.cpp")
  setup_custom_target(test_objtreewidget "${SRCROOT}${RKGUIQTDIR}")
  target_link_libraries(test_objtreewidget reak_gui_qt reak_topologies reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
  target_link_libraries(test_objtreewidget ${QT_LIBRARIES})
  target_link_libraries(test_objtreewidget ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})


  if(COIN3D_FOUND)
    set( PLANNER_3DVIEWER_UIS
      "${SRCROOT}${RKGUIQTDIR}/ui/planner_3Dview.ui" 
      "${SRCROOT}${RKGUIQTDIR}/ui/planner_space_config.ui"
      "${SRCROOT}${RKGUIQTDIR}/ui/planner_alg_config.ui"
    )
    QT4_WRAP_UI( PLANNER_3DVIEWER_UI_HDRS ${PLANNER_3DVIEWER_UIS} )
    set( PLANNER_3DVIEWER_MOC_HDRS
      "${SRCROOT}${RKGUIQTDIR}/planner_3Dviewer.h"
    )
    QT4_WRAP_CPP( PLANNER_3DVIEWER_MOC_SRCS ${PLANNER_3DVIEWER_MOC_HDRS} )

    set(PLANNER_3DVIEWER_SOURCES 
      "${SRCROOT}${RKGUIQTDIR}/planner_3Dviewer.cpp"
    )
    
    add_executable(planner_3Dviewer ${PLANNER_3DVIEWER_SOURCES} ${PLANNER_3DVIEWER_MOC_SRCS} ${PLANNER_3DVIEWER_UI_HDRS})
    setup_custom_target(planner_3Dviewer "${SRCROOT}${RKGUIQTDIR}")
    target_link_libraries(planner_3Dviewer reak_gui_qt reak_kte_coin reak_geom_coin reak_topologies reak_interp reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
    target_link_libraries(planner_3Dviewer ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
    target_link_libraries(planner_3Dviewer ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})
    
    
    
    set( KTE_MODEL_VIEWER_UIS
      "${SRCROOT}${RKGUIQTDIR}/ui/kte_model_viewer.ui" 
    )
    QT4_WRAP_UI( KTE_MODEL_VIEWER_UI_HDRS ${KTE_MODEL_VIEWER_UIS} )
    set( KTE_MODEL_VIEWER_MOC_HDRS
      "${SRCROOT}${RKGUIQTDIR}/kte_model_viewer.h"
    )
    QT4_WRAP_CPP( KTE_MODEL_VIEWER_MOC_SRCS ${KTE_MODEL_VIEWER_MOC_HDRS} OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)

    set(KTE_MODEL_VIEWER_SOURCES 
      "${SRCROOT}${RKGUIQTDIR}/kte_model_viewer.cpp"
    )
    
    add_executable(kte_model_viewer ${KTE_MODEL_VIEWER_SOURCES} ${KTE_MODEL_VIEWER_MOC_SRCS} ${KTE_MODEL_VIEWER_UI_HDRS})
    setup_custom_target(kte_model_viewer "${SRCROOT}${RKGUIQTDIR}")
    target_link_libraries(kte_model_viewer reak_gui_qt reak_kte_coin reak_geom_coin reak_topologies reak_interp reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
    target_link_libraries(kte_model_viewer ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
    target_link_libraries(kte_model_viewer ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})
    

  endif()
  
endif()

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")
include_directories(AFTER "${SRCROOT}${RKGUIQTDIR}")




