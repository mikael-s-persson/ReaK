

reak_cc_library(
  NAME
    base_macros
  HDRS
    "defs.h"
  PUBLIC
)

reak_cc_library(
  NAME
    cnst_string
  HDRS
    "cnst_string.h"
  PUBLIC
)

reak_cc_test(
  NAME
    cnst_string_test
  SRCS
    "cnst_string_test.cc"
  DEPS
    ReaK::cnst_string
    GTest::gtest_main
)

reak_cc_library(
  NAME
    endian_conversions
  HDRS
    "endian_conversions.h"
  DEPS
    ReaK::base_macros
  PUBLIC
)

reak_cc_test(
  NAME
    endian_conversions_test
  SRCS
    "endian_conversions_test.cc"
  DEPS
    ReaK::endian_conversions
    GTest::gtest_main
)

reak_cc_library(
  NAME
    base
  HDRS
    "global_rng.h"
    "misc_math.h"
    "named_object.h"
    "scope_guard.h"
    "shared_object.h"
  DEPS
    ReaK::base_macros
    ReaK::rtti
    ReaK::serialization_base
  PUBLIC
)

if(Python3_FOUND)
  reak_cc_library(
    NAME
      py_base
    SRCS
      "py_base.cc"
    DEPS
      ReaK::base
      ${PYTHON_LIBRARIES}
      ${Boost_LIBRARIES}
      ${EXTRA_SYSTEM_LIBS}
    PUBLIC
  )
endif()
