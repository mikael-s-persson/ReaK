
set(PROXIMITY_SOURCES 
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_record_2D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_record_3D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_finder_2D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/proximity_finder_3D.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_circle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_crect.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_circle_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_crect_crect.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_crect_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_rectangle_rectangle.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_plane.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_sphere.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_ccylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_plane_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_sphere.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_ccylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_sphere_box.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ccylinder_ccylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ccylinder_cylinder.cpp"
  "${SRCROOT}${RKPROXIMITYDIR}/prox_ccylinder_box.cpp"
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
  "${RKPROXIMITYDIR}/prox_circle_crect.hpp"
  "${RKPROXIMITYDIR}/prox_circle_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_crect_crect.hpp"
  "${RKPROXIMITYDIR}/prox_crect_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_rectangle_rectangle.hpp"
  "${RKPROXIMITYDIR}/prox_plane_plane.hpp"
  "${RKPROXIMITYDIR}/prox_plane_sphere.hpp"
  "${RKPROXIMITYDIR}/prox_plane_ccylinder.hpp"
  "${RKPROXIMITYDIR}/prox_plane_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_plane_box.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_sphere.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_ccylinder.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_sphere_box.hpp"
  "${RKPROXIMITYDIR}/prox_ccylinder_ccylinder.hpp"
  "${RKPROXIMITYDIR}/prox_ccylinder_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_ccylinder_box.hpp"
  "${RKPROXIMITYDIR}/prox_cylinder_cylinder.hpp"
  "${RKPROXIMITYDIR}/prox_cylinder_box.hpp"
  "${RKPROXIMITYDIR}/prox_box_box.hpp"
)

add_library(reak_geom_prox STATIC ${PROXIMITY_SOURCES})
setup_custom_target(reak_geom_prox "${SRCROOT}${RKPROXIMITYDIR}")
setup_headers("${PROXIMITY_HEADERS}" "${RKPROXIMITYDIR}")
target_link_libraries(reak_geom_prox reak_geom reak_core)

add_executable(test_nlp_proximity "${SRCROOT}${RKPROXIMITYDIR}/test_nlp_proximity.cpp")
setup_custom_target(test_nlp_proximity "${SRCROOT}${RKPROXIMITYDIR}")
target_link_libraries(test_nlp_proximity reak_geom_prox reak_geom reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
include_directories(AFTER "${SRCROOT}${RKGEOMETRYDIR}")




