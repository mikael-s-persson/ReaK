
if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(test_objtreewidget "test_objtreewidget.cpp")
  ReaK_setup_target(test_objtreewidget)
  target_link_libraries(test_objtreewidget reak_core_qt reak_core)
  target_link_libraries(test_objtreewidget ${QT_LIBRARIES})
  
endif()

