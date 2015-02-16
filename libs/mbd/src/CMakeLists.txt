
add_library(reak_mbd STATIC 
  "chaser_target_model_data.cpp"
  "damper.cpp"
  "driving_actuator.cpp"
  "dry_revolute_joint.cpp"
  "flexible_beam.cpp"
  "force_actuator.cpp"
  "free_floating_platform.cpp"
  "free_joints.cpp"
  "inertia.cpp"
  "inertial_beam.cpp"
  "joint_backlash.cpp"
  "joint_friction.cpp"
  "kte_chain_geometry.cpp"
  "line_point_mindist.cpp"
  "manip_3R3R_arm.cpp"
  "manip_3R_arm.cpp"
  "manip_clik_calculator.cpp"
  "manip_dynamics_model.cpp"
  "manip_ERA_arm.cpp"
  "manip_kinematics_helper.cpp"
  "manip_kinematics_model.cpp"
  "manip_P3R3R_arm.cpp"
  "manip_SCARA_arm.cpp"
  "manip_SSRMS_arm.cpp"
  "manipulator_model.cpp"
  "mass_matrix_calculator.cpp"
  "navigation_model_data.cpp"
  "plane_point_mindist.cpp"
  "prismatic_joint.cpp"
  "revolute_joint.cpp"
  "rigid_link.cpp"
  "spring.cpp"
  "torsion_damper.cpp"
  "torsion_spring.cpp"
  "uav_kinematics.cpp"
  "virtual_kte_interface.cpp"
  "vmc_revolute_joint.cpp"
)
ReaK_setup_static_library(reak_mbd)
target_link_libraries(reak_mbd reak_geom reak_math)

if(COIN3D_FOUND)
  
  add_library(reak_mbd_coin STATIC 
    "oi_reader.cpp"
    "oi_scene_graph.cpp"
    "oi_scene_graph_kte.cpp"
    "frame_tracer_coin3d_impl.cpp"
  )
  ReaK_setup_static_library(reak_mbd_coin)
  target_link_libraries(reak_mbd_coin ${COIN3D_LIBRARIES} reak_mbd reak_geom)
  
endif()


