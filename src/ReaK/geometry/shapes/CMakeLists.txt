
set(SHAPES_SOURCES 
  "${SRCROOT}${RKSHAPESDIR}/box.cpp"
  "${SRCROOT}${RKSHAPESDIR}/capped_cylinder.cpp"
  "${SRCROOT}${RKSHAPESDIR}/capped_rectangle.cpp"
  "${SRCROOT}${RKSHAPESDIR}/circle.cpp"
  "${SRCROOT}${RKSHAPESDIR}/composite_shape_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/composite_shape_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/coord_arrows_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/coord_arrows_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/cylinder.cpp"
  "${SRCROOT}${RKSHAPESDIR}/geometry_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/geometry_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/grid_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/grid_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/line_seg_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/line_seg_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/plane.cpp"
  "${SRCROOT}${RKSHAPESDIR}/rectangle.cpp"
  "${SRCROOT}${RKSHAPESDIR}/shape_2D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/shape_3D.cpp"
  "${SRCROOT}${RKSHAPESDIR}/sphere.cpp"
)


set(SHAPES_HEADERS 
  "${RKSHAPESDIR}/base_render_pimpl.hpp"
  "${RKSHAPESDIR}/box.hpp"
  "${RKSHAPESDIR}/capped_cylinder.hpp"
  "${RKSHAPESDIR}/capped_rectangle.hpp"
  "${RKSHAPESDIR}/circle.hpp"
  "${RKSHAPESDIR}/composite_shape_2D.hpp"
  "${RKSHAPESDIR}/composite_shape_3D.hpp"
  "${RKSHAPESDIR}/coord_arrows_2D.hpp"
  "${RKSHAPESDIR}/coord_arrows_3D.hpp"
  "${RKSHAPESDIR}/cylinder.hpp"
  "${RKSHAPESDIR}/geometry_2D.hpp"
  "${RKSHAPESDIR}/geometry_3D.hpp"
  "${RKSHAPESDIR}/grid_2D.hpp"
  "${RKSHAPESDIR}/grid_3D.hpp"
  "${RKSHAPESDIR}/line_seg_2D.hpp"
  "${RKSHAPESDIR}/line_seg_3D.hpp"
  "${RKSHAPESDIR}/plane.hpp"
  "${RKSHAPESDIR}/rectangle.hpp"
  "${RKSHAPESDIR}/shape_2D.hpp"
  "${RKSHAPESDIR}/shape_3D.hpp"
  "${RKSHAPESDIR}/sphere.hpp"
)

add_library(reak_geom STATIC ${SHAPES_SOURCES})
setup_custom_target(reak_geom "${SRCROOT}${RKSHAPESDIR}")
setup_headers("${SHAPES_HEADERS}" "${RKSHAPESDIR}")
target_link_libraries(reak_geom reak_core)

add_library(reak_geom_nr STATIC "${SRCROOT}${RKSHAPESDIR}/no_render_pimples.cpp")
setup_custom_target(reak_geom_nr "${SRCROOT}${RKSHAPESDIR}")
target_link_libraries(reak_geom_nr reak_geom reak_core)


#add_executable(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}/test_bm.cpp")
#setup_custom_target(test_mbd_basic "${SRCROOT}${RKMBDKTEDIR}")

#target_link_libraries(test_mbd_basic ${Boost_LIBRARIES})
#target_link_libraries(test_mbd_basic reak_mbd_kte reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




