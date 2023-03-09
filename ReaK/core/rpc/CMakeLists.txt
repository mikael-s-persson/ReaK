


reak_cc_library(
  NAME
    rpc
  SRCS
    "detail/disc_headers.cc"
    "detail/remote_function.cc"
    "detail/rpc_headers.cc"
    "rpc_exceptions.cc"
    "rpc_server.cc"
    "rpc_version.cc"
  HDRS
    "detail/disc_headers.h"
    "detail/remote_function.h"
    "detail/rpc_function_helpers.h"
    "detail/rpc_headers.h"
    "rpc_exceptions.h"
    "rpc_function.h"
    "rpc_server.h"
    "version.h"
  DEPS
    ReaK::base
    ReaK::serialization
    ${Boost_LIBRARIES}   # For network recorder with boost::asio
    ${EXTRA_SYSTEM_LIBS}
  PUBLIC
)

reak_cc_binary(
  NAME
    rpc_calls_test
  SRCS
    "rpc_calls_test.cc"
  DEPS
    ReaK::rpc
)

reak_cc_binary(
  NAME
    rpc_server_test
  SRCS
    "rpc_server_test.cc"
  DEPS
    ReaK::rpc
)
