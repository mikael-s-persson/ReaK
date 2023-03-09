

if(COIN3D_FOUND)
  
  reak_cc_library(
    NAME
      mbd_coin
    SRCS
      "oi_reader.cc"
      "oi_scene_graph.cc"
      "oi_scene_graph_kte.cc"
      "frame_tracer_coin3d_impl.cc"
    HDRS
      "oi_reader.h"
      "oi_scene_graph.h"
      "frame_tracer_coin3d_impl.h"
    DEPS
      ReaK::shapes
      ReaK::proximity
      ReaK::kte
      ${COIN3D_LIBRARIES}
    PUBLIC
  )

  if(QT_FOUND)
	  reak_cc_binary(
	    NAME
	      oi_scene_graph_test
	    SRCS
	      oi_scene_graph_test.cc
	    DEPS
	      ReaK::shapes
        ReaK::proximity
	      ReaK::mbd_coin
	      ${COIN3D_LIBRARIES}
	  )

	  reak_cc_binary(
	    NAME
	      robot_kin_scene
	    SRCS
	      robot_kin_scene.cc
	    DEPS
	      ReaK::shapes
        ReaK::proximity
	      ReaK::mbd_coin
	      ${COIN3D_LIBRARIES}
	  )
	endif()

  reak_cc_binary(
    NAME
      vrml2reak
    SRCS
      vrml2reak.cc
    DEPS
      absl::flags
      absl::flags_parse
      ReaK::serialization
      ReaK::shapes
      ReaK::proximity
      ReaK::mbd_coin
      ${COIN3D_LIBRARIES}
  )
  
endif()
