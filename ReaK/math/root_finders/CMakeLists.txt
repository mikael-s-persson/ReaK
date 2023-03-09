

reak_cc_library(
  NAME
    root_finders
  HDRS
    "bisection_method.h"
    "broyden_method.h"
    "newton_raphson_method.h"
    "secant_method.h"
    "special_secant_method.h"
  DEPS
    ReaK::base
    ReaK::lin_alg
  PUBLIC
)

reak_cc_test(
  NAME
    root_finders_test
  SRCS
    "root_finders_test.cc"
  DEPS
    ReaK::lin_alg
    ReaK::root_finders
    GTest::gtest_main
)

reak_cc_binary(
  NAME
    root_finders_perf
  SRCS
    "root_finders_perf.cc"
  DEPS
    ReaK::lin_alg
    ReaK::root_finders
)
