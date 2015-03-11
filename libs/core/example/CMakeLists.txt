
if(QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(test_objtreewidget "test_objtreewidget.cpp")
  ReaK_setup_target(test_objtreewidget)
  target_link_libraries(test_objtreewidget reak_core_qt reak_core)
  target_link_libraries(test_objtreewidget ${QT_LIBRARIES})
  
endif()

if(_COMPILER_SUPPORTS_CXX11)
  add_executable(test_rpc_calls "test_rpc_calls.cpp")
  ReaK_setup_target(test_rpc_calls)
  target_link_libraries(test_rpc_calls reak_core)

  add_executable(test_rpc_sv "test_rpc_sv.cpp")
  ReaK_setup_target(test_rpc_sv)
  target_link_libraries(test_rpc_sv reak_core)
endif()
