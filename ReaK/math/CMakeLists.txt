# project(ReaK_math)

add_subdirectory(lin_alg)
add_subdirectory(kinetostatics)
add_subdirectory(root_finders)
add_subdirectory(optimization)
add_subdirectory(integrators)

reak_cc_library(
  NAME
    math
  DEPS
    ReaK::lin_alg
    ReaK::kinetostatics
    ReaK::root_finders
    ReaK::optimization
    ReaK::math_integrators
  PUBLIC
)

message("Configured: ReaK Math Library Components")



