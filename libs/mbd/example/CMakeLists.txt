
  
add_executable(test_am "test_am.cpp")
ReaK_setup_target(test_am)
target_link_libraries(test_am reak_mbd)

add_executable(test_bm "test_bm.cpp")
ReaK_setup_target(test_bm)
target_link_libraries(test_bm reak_mbd)

add_executable(build_airship3D_target "build_airship3D_target.cpp")
ReaK_setup_target(build_airship3D_target)
target_link_libraries(build_airship3D_target reak_mbd)
target_link_libraries(build_airship3D_target absl::flags absl::flags_parse)

add_executable(build_P3R3R_model "build_P3R3R_model.cpp")
ReaK_setup_target(build_P3R3R_model)
target_link_libraries(build_P3R3R_model reak_mbd)
target_link_libraries(build_P3R3R_model absl::flags absl::flags_parse)

add_executable(CRS_jtctrl_to_EE "CRS_jtctrl_to_EE.cpp")
ReaK_setup_target(CRS_jtctrl_to_EE)
target_link_libraries(CRS_jtctrl_to_EE reak_mbd)
target_link_libraries(CRS_jtctrl_to_EE absl::flags absl::flags_parse)

add_executable(get_CRS_rel_world_transform "get_CRS_rel_world_transform.cpp")
ReaK_setup_target(get_CRS_rel_world_transform)
target_link_libraries(get_CRS_rel_world_transform reak_mbd)
target_link_libraries(get_CRS_rel_world_transform absl::flags_parse)

add_executable(pose_to_grapple "pose_to_grapple.cpp")
ReaK_setup_target(pose_to_grapple)
target_link_libraries(pose_to_grapple reak_mbd)
target_link_libraries(pose_to_grapple absl::flags absl::flags_parse)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  include_directories( AFTER "${CMAKE_CURRENT_BINARY_DIR}" )
  
  add_executable(test_oisg "test_oisg.cpp")
  ReaK_setup_target(test_oisg)
  target_link_libraries(test_oisg reak_mbd_coin)
  target_link_libraries(test_oisg ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(robot_kin_scene "robot_kin_scene.cpp")
  ReaK_setup_target(robot_kin_scene)
  target_link_libraries(robot_kin_scene reak_mbd_coin)
  target_link_libraries(robot_kin_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  
  set( KTE_MODEL_VIEWER_UIS
    "kte_model_viewer.ui" 
  )
  QT4_WRAP_UI( KTE_MODEL_VIEWER_UI_HDRS ${KTE_MODEL_VIEWER_UIS} )
  set( KTE_MODEL_VIEWER_MOC_HDRS
    "kte_model_viewer.hpp"
  )
  QT4_WRAP_CPP( KTE_MODEL_VIEWER_MOC_SRCS ${KTE_MODEL_VIEWER_MOC_HDRS} OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED)
  
  set(KTE_MODEL_VIEWER_SOURCES 
    "kte_model_viewer.cpp"
    "kte_mv_kte_loader.cpp"
  )
  
  add_executable(kte_model_viewer ${KTE_MODEL_VIEWER_SOURCES} ${KTE_MODEL_VIEWER_MOC_SRCS} ${KTE_MODEL_VIEWER_UI_HDRS})
  ReaK_setup_target(kte_model_viewer)
  target_link_libraries(kte_model_viewer reak_mbd_qt reak_mbd_coin)
  target_link_libraries(kte_model_viewer ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()

