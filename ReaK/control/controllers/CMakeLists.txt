
reak_cc_library(
  NAME
    controllers
  HDRS
    "IHAQR_topology.h"
    "lqr_controllers.h"
    "MEAQR_topology.h"
    "ss_controller_concept.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
    ReaK::spaces
  PUBLIC
)
