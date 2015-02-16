
add_library(reak_control STATIC 
  "airship_assembled_models.cpp"
  "near_buoyant_airship_models.cpp"
  "quadrotor_system.cpp"
  "satellite_invar_models.cpp"
  "satellite_basic_models.cpp"
  "satellite_modeling_options.cpp"
)
ReaK_setup_static_library(reak_control)
target_link_libraries(reak_control reak_topologies reak_math)

