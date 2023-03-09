

reak_cc_library(
  NAME
    math_integrators
  HDRS
    "fixed_step_integrators.h"
    "integration_exceptions.h"
    "integrator.h"
    "pred_corr_integrators.h"
    "variable_step_integrators.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
  PUBLIC
)

reak_cc_test(
  NAME
    integrators_test
  SRCS
    "integrators_test.cc"
  DEPS
    ReaK::serialization
    ReaK::lin_alg
    ReaK::math_integrators
    GTest::gtest_main
)
