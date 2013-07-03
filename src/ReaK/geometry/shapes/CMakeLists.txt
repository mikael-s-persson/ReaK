
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
  "${SRCROOT}${RKSHAPESDIR}/kte_chain_geometry.cpp"
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
  "${RKSHAPESDIR}/kte_chain_geometry.hpp"
  "${RKSHAPESDIR}/line_seg_2D.hpp"
  "${RKSHAPESDIR}/line_seg_3D.hpp"
  "${RKSHAPESDIR}/plane.hpp"
  "${RKSHAPESDIR}/rectangle.hpp"
  "${RKSHAPESDIR}/shape_2D.hpp"
  "${RKSHAPESDIR}/shape_3D.hpp"
  "${RKSHAPESDIR}/sphere.hpp"
)

add_library(reakobj_geom OBJECT ${SHAPES_SOURCES})


if(COIN3D_FOUND AND QT_FOUND)
  
  set(MBDCOIN_SOURCES
    "${SRCROOT}${RKSHAPESDIR}/oi_scene_graph.cpp"
    "${SRCROOT}${RKSHAPESDIR}/oi_reader.cpp"
    "${SRCROOT}${RKSHAPESDIR}/frame_tracer_coin3d_impl.cpp"
    "${SRCROOT}${RKSHAPESDIR}/oi_scene_graph_kte.cpp"
  )
  
  add_library(reak_mbd_coin STATIC ${MBDCOIN_SOURCES} $<TARGET_OBJECTS:reakobj_geom> $<TARGET_OBJECTS:reakobj_proximity> $<TARGET_OBJECTS:reakobj_mbd_kte> $<TARGET_OBJECTS:reakobj_kte_models>)
  setup_custom_target(reak_mbd_coin "${SRCROOT}${RKSHAPESDIR}")
  target_link_libraries(reak_mbd_coin ${COIN3D_LIBRARIES} reak_core ${Boost_LIBRARIES})
  
  
  include( ${QT_USE_FILE} )
  
  add_executable(test_coin3d "${SRCROOT}${RKSHAPESDIR}/test_coin3d.cpp")
  setup_custom_target(test_coin3d "${SRCROOT}${RKSHAPESDIR}")
  target_link_libraries(test_coin3d ${SOQT4_LIBRARIES} ${COIN3D_LIBRARIES} ${QT_LIBRARIES})
#  target_link_libraries(test_coin3d reak_mbd reak_core)
  
  
  add_executable(test_oisg "${SRCROOT}${RKSHAPESDIR}/test_oisg.cpp")
  setup_custom_target(test_oisg "${SRCROOT}${RKSHAPESDIR}")
  target_link_libraries(test_oisg reak_mbd_coin reak_geom_prox reak_core)
  target_link_libraries(test_oisg ${SOQT4_LIBRARIES} ${QT_LIBRARIES})
  
endif()


if(COIN3D_FOUND)
  
  add_executable(vrml2reak "${SRCROOT}${RKSHAPESDIR}/vrml2reak.cpp")
  setup_custom_target(vrml2reak "${SRCROOT}${RKSHAPESDIR}")
  target_link_libraries(vrml2reak reak_mbd_coin reak_geom_prox reak_core)
  
endif()




