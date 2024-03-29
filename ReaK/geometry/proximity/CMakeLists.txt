

reak_cc_library(
  NAME
    proximity
  SRCS
    "prox_box_box.cc"
    "prox_ccylinder_box.cc"
    "prox_ccylinder_ccylinder.cc"
    "prox_ccylinder_cylinder.cc"
    "prox_circle_circle.cc"
    "prox_circle_crect.cc"
    "prox_circle_rectangle.cc"
    "prox_crect_crect.cc"
    "prox_crect_rectangle.cc"
    "prox_cylinder_box.cc"
    "prox_cylinder_cylinder.cc"
    "prox_fundamentals_2D.cc"
    "prox_fundamentals_3D.cc"
    "proximity_finder_2D.cc"
    "proximity_finder_3D.cc"
    "proximity_record_2D.cc"
    "proximity_record_3D.cc"
    "prox_plane_box.cc"
    "prox_plane_ccylinder.cc"
    "prox_plane_cylinder.cc"
    "prox_plane_plane.cc"
    "prox_plane_sphere.cc"
    "prox_rectangle_rectangle.cc"
    "prox_sphere_box.cc"
    "prox_sphere_ccylinder.cc"
    "prox_sphere_cylinder.cc"
    "prox_sphere_sphere.cc"
    "proxy_query_model.cc"
  HDRS
    "prox_box_box.h"
    "prox_ccylinder_box.h"
    "prox_ccylinder_ccylinder.h"
    "prox_ccylinder_cylinder.h"
    "prox_circle_circle.h"
    "prox_circle_crect.h"
    "prox_circle_rectangle.h"
    "prox_crect_crect.h"
    "prox_crect_rectangle.h"
    "prox_cylinder_box.h"
    "prox_cylinder_cylinder.h"
    "prox_fundamentals_2D.h"
    "prox_fundamentals_3D.h"
    "proximity_finder_2D.h"
    "proximity_finder_3D.h"
    "proximity_record_2D.h"
    "proximity_record_3D.h"
    "prox_plane_box.h"
    "prox_plane_ccylinder.h"
    "prox_plane_cylinder.h"
    "prox_plane_plane.h"
    "prox_plane_sphere.h"
    "prox_rectangle_rectangle.h"
    "prox_sphere_box.h"
    "prox_sphere_ccylinder.h"
    "prox_sphere_cylinder.h"
    "prox_sphere_sphere.h"
    "proxy_query_model.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
  PUBLIC
)

reak_cc_binary(
  NAME
    proximity_perf
  SRCS
    "proximity_perf.cc"
  DEPS
    ReaK::base
    ReaK::shapes
    ReaK::proximity
)

reak_cc_binary(
  NAME
    build_X8_quadrotor_geom
  SRCS
    "build_X8_quadrotor_geom.cc"
  DEPS
    ReaK::base
    ReaK::serialization
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
    ReaK::proximity
)

reak_cc_binary(
  NAME
    build_X8_obstacle_courses
  SRCS
    "build_X8_obstacle_courses.cc"
  DEPS
    ReaK::base
    ReaK::serialization
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
    ReaK::proximity
)

reak_cc_binary(
  NAME
    build_MD148_lab
  SRCS
    "build_MD148_lab.cc"
  DEPS
    absl::flags
    absl::flags_parse
    ReaK::base
    ReaK::serialization
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
    ReaK::proximity
)

reak_cc_binary(
  NAME
    proximity_nlp_test
  SRCS
    "proximity_nlp_test.cc"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
    ReaK::proximity
)

reak_cc_binary(
  NAME
    proximity_gjk_test
  SRCS
    "proximity_gjk_test.cc"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::shapes
    ReaK::proximity
)
