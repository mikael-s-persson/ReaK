
set(ROBOTAIRSHIP_SOURCES 
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_2D_analog.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_geom_model.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/CRS_A465_models.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/satellite_basic_models.cpp"
  "${SRCROOT}${RKROBOTAIRSHIPDIR}/satellite_invar_models.cpp"
)


set(ROBOTAIRSHIP_HEADERS 
  "${RKROBOTAIRSHIPDIR}/airship2D_lin_model.hpp"
  "${RKROBOTAIRSHIPDIR}/airship3D_lin_model.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_A465_2D_analog.hpp"
  "${RKROBOTAIRSHIPDIR}/CRS_A465_models.hpp"
  "${RKROBOTAIRSHIPDIR}/rigid_body_states.hpp"
  "${RKROBOTAIRSHIPDIR}/satellite_basic_models.hpp"
  "${RKROBOTAIRSHIPDIR}/satellite_invar_models.hpp"
)

add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")
target_link_libraries(reak_robot_airship reak_mbd_kte reak_geom_prox reak_geom reak_core)

add_executable(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship2D.cpp")
setup_custom_target(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(build_airship2D ${Boost_LIBRARIES})
target_link_libraries(build_airship2D reak_mbd_kte reak_core)

add_executable(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship2D.cpp")
setup_custom_target(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(run_airship2D ${Boost_LIBRARIES})
target_link_libraries(run_airship2D reak_mbd_kte reak_core)

add_executable(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship2D.cpp")
setup_custom_target(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship2D ${Boost_LIBRARIES})
target_link_libraries(estimate_airship2D reak_core)


  
add_executable(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D.cpp")
setup_custom_target(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  
target_link_libraries(build_airship3D ${Boost_LIBRARIES})
target_link_libraries(build_airship3D reak_mbd_kte reak_geom_prox reak_geom reak_core)
  

add_executable(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D_transforms.cpp")
setup_custom_target(build_airship3D_transforms "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(build_airship3D_transforms ${Boost_LIBRARIES})
target_link_libraries(build_airship3D_transforms reak_mbd_kte reak_core)

add_executable(build_CRS_A465_model "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_CRS_A465_model.cpp")
setup_custom_target(build_CRS_A465_model "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(build_CRS_A465_model reak_robot_airship reak_mbd_kte reak_geom_prox reak_geom reak_core)
#target_link_libraries(build_CRS_A465_model ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  

if(COIN3D_FOUND AND QT_FOUND)

  include( ${QT_USE_FILE} )
  
  add_executable(test_lab_scene "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_lab_scene.cpp")
  setup_custom_target(test_lab_scene "${SRCROOT}${RKROBOTAIRSHIPDIR}")
  target_link_libraries(test_lab_scene reak_geom_coin reak_robot_airship reak_mbd_kte reak_geom_prox reak_geom reak_core)
  target_link_libraries(test_lab_scene ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
endif()


add_executable(test_CRS_IK "${SRCROOT}${RKROBOTAIRSHIPDIR}/test_CRS_IK.cpp")
setup_custom_target(test_CRS_IK "${SRCROOT}${RKROBOTAIRSHIPDIR}")
target_link_libraries(test_CRS_IK reak_robot_airship reak_mbd_kte reak_core)

add_executable(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship3D.cpp")
setup_custom_target(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(run_airship3D ${Boost_LIBRARIES})
target_link_libraries(run_airship3D reak_mbd_kte reak_core)

add_executable(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D.cpp")
setup_custom_target(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D ${Boost_LIBRARIES})
target_link_libraries(estimate_airship3D reak_core)

add_executable(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_range.cpp")
setup_custom_target(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D_range ${Boost_LIBRARIES})
target_link_libraries(estimate_airship3D_range reak_core)

add_executable(estimate_airship3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_spins.cpp")
setup_custom_target(estimate_airship3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D_spins ${Boost_LIBRARIES})
target_link_libraries(estimate_airship3D_spins reak_core)

add_executable(estimate_satellite3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_satellite3D_spins.cpp")
setup_custom_target(estimate_satellite3D_spins "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_satellite3D_spins ${Boost_LIBRARIES})
target_link_libraries(estimate_satellite3D_spins reak_robot_airship reak_core)

add_executable(estimate_airship3D_mc "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_mc.cpp")
setup_custom_target(estimate_airship3D_mc "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D_mc ${Boost_LIBRARIES})
target_link_libraries(estimate_airship3D_mc reak_mbd_kte reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")




