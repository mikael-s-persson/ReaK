

reak_cc_library(
  NAME
    rtti
  SRCS
    "so_type.cc"
    "so_type_repo.cc"
  HDRS
    "rtti.h"
    "so_register_type.h"
    "so_type.h"
    "so_type_repo.h"
    "typed_containers.h"
    "typed_object.h"
    "typed_primitives.h"
  DEPS
    ReaK::cnst_string
  PUBLIC
)

reak_cc_test(
  NAME
    rtti_test
  SRCS
    "rtti_test.cc"
  DEPS
    ReaK::rtti
    GTest::gtest_main
)
