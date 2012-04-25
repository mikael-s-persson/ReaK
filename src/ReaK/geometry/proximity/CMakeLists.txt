
set(PROXIMITY_SOURCES 
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_record_2D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_record_3D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_finder_2D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_finder_3D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_circle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_ellipse.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ellipse_ellipse.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ellipse_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_rectangle_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_plane.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_sphere.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_ellipsoid.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_sphere.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_ellipsoid.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ellipsoid_ellipsoid.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ellipsoid_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ellipsoid_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_cylinder_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_cylinder_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_box_box.cpp"
)


set(PROXIMITY_HEADERS 
  "${RKPROXIMITYDIR}/proximity_record_2D.hpp"
  "${RKPROXIMITYDIR}/proximity_record_3D.hpp"
  "${RKPROXIMITYDIR}/proximity_finder_2D.hpp"
  "${RKPROXIMITYDIR}/proximity_finder_3D.hpp"
  "${RKPROXIMITYDIR}/prox_circle_circle.hpp"
  "${RKPROXIMITYDIR}/prox_circle_ellipse.hpp"
  "${RKPROXIMITYDIR}/prox_circle_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_ellipse_ellipse.hpp"
  "${RKPROXIMITYDIR}/prox_ellipse_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_rectangle_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_plane_plane.hpp"
  "${RKPROXIMITYDIR}/prox_plane_sphere.hpp"
  "${RKPROXIMITYDIR}/prox_plane_ellipsoid.hpp"
  "${RKPROXIMITYDIR}/prox_plane_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_plane_box.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_sphere.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_ellipsoid.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_box.hpp"
  "${RKPROXIMITYDIR}/prox_ellipsoid_ellipsoid.hpp"
  "${RKPROXIMITYDIR}/prox_ellipsoid_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_ellipsoid_box.hpp"
  "${RKPROXIMITYDIR}/prox_cylinder_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_cylinder_box.hpp"
  "${RKPROXIMITYDIR}/prox_box_box.hpp"
)

add_library(reak_geom_prox STATIC ${PROXIMITY_SOURCES})
setup_custom_target(reak_geom_prox "${SRCROOT}${RKPROXIMITYDIR}")
setup_headers("${PROXIMITY_HEADERS}" "${RKPROXIMITYDIR}")
target_link_libraries(reak_geom_prox reak_geom reak_core)

#add_executable(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}/test_bm.cpp")
#setup_custom_target(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}")

#target_link_libraries(test_mbd_basic ${Boost_LIBRARIES})
#target_link_libraries(test_mbd_basic reak_mbd_kte reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")




