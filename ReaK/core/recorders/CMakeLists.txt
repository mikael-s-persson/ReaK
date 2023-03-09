

reak_cc_library(
  NAME
    recorders
  SRCS
    "ascii_recorder.cc"
    "bin_recorder.cc"
    "data_record.cc"
    "network_recorder.cc"
    "vector_recorder.cc"
  HDRS
    "ascii_recorder.h"
    "bin_recorder.h"
    "data_record.h"
    "network_recorder.h"
    "vector_recorder.h"
  DEPS
    ReaK::base_macros
    ReaK::rtti
    ReaK::endian_conversions
    ${Boost_LIBRARIES}   # For network recorder with boost::asio
    ${EXTRA_SYSTEM_LIBS}
  PUBLIC
)

reak_cc_test(
  NAME
    recorders_test
  SRCS
    "recorders_test.cc"
  DEPS
    ReaK::recorders
    ReaK::rtti
    GTest::gtest_main
)

reak_cc_library(
  NAME
    recorders_po
  SRCS
    "data_record_options.cc"
    "data_record_po.cc"
  HDRS
    "data_record_options.h"
    "data_record_po.h"
  DEPS
    absl::flags
    absl::strings
    ReaK::recorders
    ReaK::serialization
    ReaK::base_macros
  PUBLIC
)

reak_cc_binary(
  NAME
    convert_datastream
  SRCS
    "convert_datastream.cc"
  DEPS
    absl::flags
    absl::flags_parse
    ReaK::rtti
    ReaK::recorders
    ReaK::recorders_po
)
