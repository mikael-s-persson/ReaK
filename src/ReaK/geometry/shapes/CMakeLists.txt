
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

#add_library(reak_geom_nr STATIC "${SRCROOT}${RKSHAPESDIR}/no_render_pimples.cpp")
#setup_custom_target(reak_geom_nr "${SRCROOT}${RKSHAPESDIR}")
#target_link_libraries(reak_geom_nr reak_geom reak_core)



find_package(Coin3D)

if(COIN3D_FOUND)
  message(STATUS "Coin3D library components were found: '${COIN3D_LIBRARIES}'")
else()
  message(WARNING "Coin3D was not detected on this system, some test programs will not be built!")
endif()

find_package( Qt4 )

if(COIN3D_FOUND AND QT_FOUND)
  FIND_LIBRARY(SOQT4_LIBRARIES SoQt4)

  add_library(reak_geom_coin STATIC "${SRCROOT}${RKSHAPESDIR}/oi_scene_graph.cpp")
  setup_custom_target(reak_geom_coin "${SRCROOT}${RKSHAPESDIR}")
  target_link_libraries(reak_geom_coin ${COIN3D_LIBRARIES} reak_geom reak_core)
  
  include( ${QT_USE_FILE} )
  
  add_executable(test_coin3d "${SRCROOT}${RKSHAPESDIR}/test_coin3d.cpp")
  setup_custom_target(test_coin3d "${SRCROOT}${RKSHAPESDIR}")
  
  target_link_libraries(test_coin3d ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
#  target_link_libraries(test_coin3d reak_mbd_kte reak_core)
  
  
  add_executable(test_oisg "${SRCROOT}${RKSHAPESDIR}/test_oisg.cpp")
  setup_custom_target(test_oisg "${SRCROOT}${RKSHAPESDIR}")
  
  target_link_libraries(test_oisg ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
  target_link_libraries(test_oisg reak_geom_coin reak_geom reak_core)
  
  include_directories(BEFORE ${COIN3D_INCLUDE_DIRS})
  include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
endif()



#target_link_libraries(test_mbd_basic ${Boost_LIBRARIES})
#target_link_libraries(test_mbd_basic reak_mbd_kte reak_core)

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




