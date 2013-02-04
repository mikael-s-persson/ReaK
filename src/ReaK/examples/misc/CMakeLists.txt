
set(MISCEXAMPLES_SOURCES 
#   "${SRCROOT}${RKMISCEXAMPLESDIR}/CRS_A465_2D_analog.cpp"
)


set(MISCEXAMPLES_HEADERS 
#   "${RKMISCEXAMPLESDIR}/airship2D_lin_model.hpp"
)

# add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
# setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
# setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")
# target_link_libraries(reak_robot_airship reak_topologies reak_interp reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)

if(COIN3D_FOUND AND QT_FOUND)
  
  include( ${QT_USE_FILE} )
  
  add_executable(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}/robot_kin_scene.cpp")
  setup_custom_target(robot_kin_scene "${SRCROOT}${RKMISCEXAMPLESDIR}")
  target_link_libraries(robot_kin_scene reak_kte_coin reak_geom_coin reak_kte_models reak_mbd_kte reak_geom_prox reak_geom reak_core)
  target_link_libraries(robot_kin_scene ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  
endif()

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")




