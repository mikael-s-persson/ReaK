

if(COIN3D_FOUND)
  
  add_executable(vrml2reak "vrml2reak.cpp")
  ReaK_setup_tool_program(vrml2reak)
  target_link_libraries(vrml2reak reak_mbd_coin)
  
endif()

