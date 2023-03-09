# project(ReaK_core)

add_subdirectory(rtti)
add_subdirectory(base)
add_subdirectory(serialization)
add_subdirectory(recorders)
add_subdirectory(rpc)
add_subdirectory(sorting)
add_subdirectory(qt)

reak_cc_library(
  NAME
    core
  DEPS
    ReaK::base
    ReaK::rtti
    ReaK::serialization
    ReaK::recorders
    ReaK::rpc
    ReaK::sorting
  PUBLIC
)

message("Configured: ReaK Core Library Components")
