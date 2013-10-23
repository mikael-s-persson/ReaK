
set(ROBOTAIRSHIP_SOURCES 
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_2D_analog.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_geom_model.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_models.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrt_planners.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrt_planners_svp.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrt_planners_sap.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrtstar_planners.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrtstar_planners_svp.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_rrtstar_planners_sap.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_prm_planners.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_prm_planners_svp.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_prm_planners_sap.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_fadprm_planners.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_fadprm_planners_svp.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_fadprm_planners_sap.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_sbastar_planners.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_sbastar_planners_svp.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_sbastar_planners_sap.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/satellite_basic_models.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/satellite_invar_models.cpp"
)


set(ROBOTAIRSHIP_HEADERS 
  "${RKROBOTAIRSHIPDIR}/airship2D_lin_model.hpp"
  "${RKROBOTAIRSHIPDIR}/airship3D_lin_model.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_A465_2D_analog.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_A465_models.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_rrt_planners.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_rrtstar_planners.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_prm_planners.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_fadprm_planners.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_sbastar_planners.hpp"
  "${RKROBOTAIRSHIPDIR}/rigid_body_states.hpp"
  "${RKROBOTAIRSHIPDIR}/satellite_basic_models.hpp"
  "${RKROBOTAIRSHIPDIR}/satellite_invar_models.hpp"
)

add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")
target_link_libraries(reak_robot_airship reak_topologies reak_interp reak_mbd reak_core)



add_executable(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship2D.cpp")
setup_custom_target(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_airship2D reak_mbd reak_core)

add_executable(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship2D.cpp")
setup_custom_target(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(run_airship2D reak_topologies reak_interp reak_mbd reak_core)

add_executable(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship2D.cpp")
setup_custom_target(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship2D reak_topologies reak_interp reak_mbd reak_core)


  
add_executable(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D.cpp")
setup_custom_target(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_airship3D reak_topologies reak_interp reak_mbd reak_core)
  

add_executable(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D_transforms.cpp")
setup_custom_target(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_airship3D_transforms reak_topologies reak_interp reak_mbd reak_core)


add_executable(build_CRS_A465_model "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_CRS_A465_model.cpp")
setup_custom_target(build_CRS_A465_model "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_CRS_A465_model reak_topologies reak_interp reak_robot_airship reak_mbd reak_core)
#target_link_libraries(build_CRS_A465_model ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
add_executable(build_P3R3R_model "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_P3R3R_model.cpp")
setup_custom_target(build_P3R3R_model "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_P3R3R_model reak_topologies reak_interp reak_mbd reak_core)


if(COIN3D_FOUND AND QT_FOUND)

  include( ${QT_USE_FILE} )
  
  add_executable(test_lab_scene "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_lab_scene.cpp")
  setup_custom_target(test_lab_scene "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(test_lab_scene reak_topologies reak_interp reak_robot_airship reak_mbd_coin reak_core)
  target_link_libraries(test_lab_scene ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(test_CRS_planning "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_CRS_planning.cpp")
  setup_custom_target(test_CRS_planning "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(test_CRS_planning reak_topologies reak_interp reak_robot_airship reak_mbd_coin reak_core)
  target_link_libraries(test_CRS_planning ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  add_executable(test_CRS_plan_compiling "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_CRS_plan_compiling.cpp")
  setup_custom_target(test_CRS_plan_compiling "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(test_CRS_plan_compiling reak_topologies reak_interp reak_mbd_coin reak_core)
  target_link_libraries(test_CRS_plan_compiling ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
#if() #disabling CRS_planner target because it is too heavy to even compile.

if(Boost_CHRONO_FOUND)

  include( ${QT_USE_FILE} )
  set( ROBOTAIRSHIP_CRSPLANNER_UIS
     "${SRCROOT}${RKROBOTAIRSHIPDIR}/uis/CRS_planner_window.ui" 
     "${SRCROOT}${RKROBOTAIRSHIPDIR}/uis/CRS_planner_config.ui"
  )
  QT4_WRAP_UI( ROBOTAIRSHIP_CRSPLANNER_UI_HDRS ${ROBOTAIRSHIP_CRSPLANNER_UIS} )
  set( ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS
     "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_impl.h"
  )
  QT4_WRAP_CPP( ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS ${ROBOTAIRSHIP_CRSPLANNER_MOC_HDRS} )
  include_directories( AFTER "${CMAKE_BINARY_DIR}${RKROBOTAIRSHIPDIR}" )
  
  set( ROBOTAIRSHIP_CRSPLANNER_SRCS
    "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_impl.cpp"
    "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner_exec.cpp"
  )
  
  add_executable(CRS_planner ${ROBOTAIRSHIP_CRSPLANNER_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_UI_HDRS})
  setup_custom_target(CRS_planner "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(CRS_planner reak_robot_airship reak_topologies reak_interp reak_mbd_coin reak_core)
  target_link_libraries(CRS_planner ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
  
  set( ROBOTAIRSHIP_CRSPLANNER2_SRCS
    "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planners_utility.cpp"
    "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner2_impl.cpp"
    "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_planner2_exec.cpp"
  )
  
  add_executable(CRS_planner2 ${ROBOTAIRSHIP_CRSPLANNER2_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_MOC_SRCS} ${ROBOTAIRSHIP_CRSPLANNER_UI_HDRS})
  setup_custom_target(CRS_planner2 "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(CRS_planner2 reak_robot_airship reak_topologies reak_interp reak_mbd_coin reak_core)
  target_link_libraries(CRS_planner2 ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()

#endif()
  
endif()


add_executable(test_CRS_IK "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_CRS_IK.cpp")
setup_custom_target(test_CRS_IK "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(test_CRS_IK reak_robot_airship reak_topologies reak_interp reak_mbd reak_core)

add_executable(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship3D.cpp")
setup_custom_target(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(run_airship3D reak_topologies reak_interp reak_mbd reak_core)

add_executable(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D.cpp")
setup_custom_target(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship3D reak_topologies reak_interp reak_core)

add_executable(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_range.cpp")
setup_custom_target(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship3D_range reak_topologies reak_interp reak_core)

add_executable(estimate_airship3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_spins.cpp")
setup_custom_target(estimate_airship3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship3D_spins reak_topologies reak_interp reak_core)

add_executable(estimate_satellite3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_satellite3D_spins.cpp")
setup_custom_target(estimate_satellite3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_satellite3D_spins reak_topologies reak_interp reak_robot_airship reak_core)

add_executable(estimate_airship3D_mc "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_mc.cpp")
setup_custom_target(estimate_airship3D_mc "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(estimate_airship3D_mc reak_topologies reak_interp reak_mbd reak_core)


add_executable(run_satellite3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_satellite3D.cpp")
setup_custom_target(run_satellite3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(run_satellite3D reak_topologies reak_interp reak_robot_airship reak_core)



