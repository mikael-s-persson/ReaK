
set(ROBOTAIRSHIP_SOURCES 
)


set(ROBOTAIRSHIP_HEADERS 
  "${RKROBOTAIRSHIPDIR}/airship2D_lin_model.hpp"
)

#add_library(reak_robot_airship STATIC ${ROBOTAIRSHIP_SOURCES})
#setup_custom_target(reak_robot_airship "${SRCROOT}${RKROBOTAIRSHIPDIR}")
#setup_headers("${ROBOTAIRSHIP_HEADERS}" "${RKROBOTAIRSHIPDIR}")

add_executable(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship2D.cpp")
setup_custom_target(build_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(build_airship2D ${BOOST_THREAD_LIB})
target_link_libraries(build_airship2D reak_mbd_kte reak_core)

add_executable(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship2D.cpp")
setup_custom_target(run_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(run_airship2D ${BOOST_THREAD_LIB})
target_link_libraries(run_airship2D reak_mbd_kte reak_core)

add_executable(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship2D.cpp")
setup_custom_target(estimate_airship2D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship2D ${BOOST_THREAD_LIB})
target_link_libraries(estimate_airship2D reak_core)

add_executable(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/build_airship3D.cpp")
setup_custom_target(build_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(build_airship3D ${BOOST_THREAD_LIB})
target_link_libraries(build_airship3D reak_mbd_kte reak_core)

add_executable(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/run_airship3D.cpp")
setup_custom_target(run_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(run_airship3D ${BOOST_THREAD_LIB})
target_link_libraries(run_airship3D reak_mbd_kte reak_core)

add_executable(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D.cpp")
setup_custom_target(estimate_airship3D "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D ${BOOST_THREAD_LIB})
target_link_libraries(estimate_airship3D reak_core)

add_executable(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}/estimate_airship3D_range.cpp")
setup_custom_target(estimate_airship3D_range "${SRCROOT}${RKROBOTAIRSHIPDIR}")

target_link_libraries(estimate_airship3D_range ${BOOST_THREAD_LIB})
target_link_libraries(estimate_airship3D_range reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKCTRLDIR}")




