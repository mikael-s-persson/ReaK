
set(ROBOTAIRSHIP_SOURCES 
#   "${SRCROOT}${RKROBOTAIRSHIPDIR}/satellite_invar_models.cpp"
)


set(ROBOTAIRSHIP_HEADERS 
#   "${RKROBOTAIRSHIPDIR}/satellite_invar_models.hpp"
)

add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")
target_link_libraries(reak_robot_airship reak_topologies reak_interp reak_mbd reak_core)


include_directories("${CMAKE_BINARY_DIR}/ReaK/gui_qt/")


add_executable(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D_transforms.cpp")
setup_custom_target(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_airship3D_transforms reak_topologies reak_interp reak_mbd reak_core)


add_executable(build_P3R3R_model "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_P3R3R_model.cpp")
setup_custom_target(build_P3R3R_model "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_P3R3R_model reak_topologies reak_interp reak_mbd reak_core)

add_executable(build_MD148_lab "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_MD148_lab.cpp")
setup_custom_target(build_MD148_lab "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_MD148_lab reak_geom_prox)

add_executable(build_airship3D_target "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D_target.cpp")
setup_custom_target(build_airship3D_target "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_airship3D_target reak_mbd)



if(COIN3D_FOUND AND QT_FOUND)

  include( ${QT_USE_FILE} )
  
  if(Boost_CHRONO_FOUND)

    include( ${QT_USE_FILE} )
    set( ROBOTAIRSHIP_CRSPLANNER_UIS
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/uis/CRS_planner_window.ui" 
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/uis/CRS_planner_config.ui"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/uis/CRS_run_dialog.ui"
    )
    QT4_WRAP_UI( ROBOTAIRSHIP_CRSPLANNER_UI_HDRS ${ROBOTAIRSHIP_CRSPLANNER_UIS} )
    set( ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_run_dialog.hpp"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_impl.hpp"
    )
    QT4_WRAP_CPP( ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS ${ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS} )
    include_directories( AFTER "${CMAKE_BINARY_DIR}${RKROBOTAIRSHIPDIR}" )
    
    set( ROBOTAIRSHIP_CRSPLANNER_SRCS
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_run_dialog.cpp"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planners_utility.cpp"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_impl.cpp"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_exec.cpp"
      "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_dynexec.cpp"
    )
    
    add_executable(CRS_planner ${ROBOTAIRSHIP_CRSPLANNER_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_UI_HDRS})
    setup_custom_target(CRS_planner "${SRCROOT}${RKROBOTAIRSHIPDIR}")
    target_link_libraries(CRS_planner reak_planners reak_ss_systems reak_topologies reak_interp reak_gui_qt reak_mbd_coin reak_core)
    target_link_libraries(CRS_planner ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
    
  endif()
  
endif()


add_executable(run_CRS_planner "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_CRS_planner.cpp")
setup_custom_target(run_CRS_planner "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(run_CRS_planner reak_planners reak_robot_airship reak_topologies reak_interp reak_mbd reak_core)


add_executable(estimate_satellite3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_satellite3D.cpp")
setup_custom_target(estimate_satellite3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_satellite3D reak_ss_systems reak_topologies reak_interp reak_robot_airship reak_core)

add_executable(estimate_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship.cpp")
setup_custom_target(estimate_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship reak_ss_systems reak_core)

add_executable(tryphon_state_estimator_stub "${SRCROOT}${RKROBOTAIRSHIPDIR}/tryphon_state_estimator_stub.cpp")
setup_custom_target(tryphon_state_estimator_stub "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(tryphon_state_estimator_stub reak_ss_systems reak_core)

add_executable(tryphon_ilqr_controller_stub "${SRCROOT}${RKROBOTAIRSHIPDIR}/tryphon_ilqr_controller_stub.cpp")
setup_custom_target(tryphon_ilqr_controller_stub "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(tryphon_ilqr_controller_stub reak_ss_systems reak_core)



