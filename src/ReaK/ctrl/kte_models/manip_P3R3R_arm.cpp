
/*
 *    Copyright 2013 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#include "manip_P3R3R_arm.hpp"

#include "mbd_kte/prismatic_joint.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/vect_alg.hpp"
#include "optimization/optim_exceptions.hpp"
#include <cmath>

namespace ReaK {

namespace kte {

//     shared_ptr< frame_3D<double> > m_base_frame;
//     shared_ptr< gen_coord<double> > m_track_coord;
//     shared_ptr< joint_dependent_frame_3D > m_output_frame;
//     vect<double,3> m_track_direction;
//     double m_track_lower_bound;
//     double m_track_upper_bound;
//     
//     manip_3R3R_kinematics m_arm_model;
//     shared_ptr< kte_map_chain > m_chain;

manip_P3R3R_kinematics::manip_P3R3R_kinematics(const std::string& aName,
                                               const shared_ptr< frame_3D<double> >& aBaseFrame,
                                               const vect<double,3>& aTrackDirection,
                                               double aBaseToShoulder, 
                                               double aShoulderToElbow,
                                               double aElbowToJoint4, 
                                               double aJoint4ToWrist,
                                               double aWristToFlange,
                                               const vect_n<double>& aJointLowerBounds,
                                               const vect_n<double>& aJointUpperBounds) :
                                               inverse_kinematics_model(aName),
                                               m_base_frame(aBaseFrame),
                                               m_track_coord(new gen_coord<double>(), scoped_deleter()),
                                               m_output_frame(new frame_3D<double>(), scoped_deleter()),
                                               m_track_direction(aTrackDirection),
                                               m_track_lower_bound(aJointLowerBounds[0]),
                                               m_track_upper_bound(aJointUpperBounds[0]),
                                               m_arm_model(
                                                 aName + "_arm",
                                                 m_output_frame,
                                                 aBaseToShoulder,
                                                 aShoulderToElbow,
                                                 aElbowToJoint4,
                                                 aJoint4ToWrist,
                                                 aWristToFlange,
                                                 vect_n<double>(aJointLowerBounds.begin()+1, aJointLowerBounds.end()),
                                                 vect_n<double>(aJointUpperBounds.begin()+1, aJointUpperBounds.end())),
                                               m_chain() {
  
  
  
  //declare all the joint jacobians.
  shared_ptr< jacobian_gen_3D<double> > track_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  
  //create prismatic joint
  shared_ptr< kte::prismatic_joint_3D > track_joint(new kte::prismatic_joint_3D(
      "manip_P3R3R_track",
      m_track_coord,
      m_track_direction,
      m_base_frame,
      m_output_frame,
      track_jacobian),
    scoped_deleter());
  
  m_chain = shared_ptr< kte_map_chain >(new kte_map_chain("manip_P3R3R_kin_model"), scoped_deleter());
  
  (*m_chain) << track_joint << m_arm_model.getKTEChain();
  
  m_arm_model.getDependentFrame3D(0)->add_joint(m_track_coord, track_jacobian);
  
};



void manip_P3R3R_kinematics::doDirectMotion() {
  m_chain->doMotion();
};


void manip_P3R3R_kinematics::doInverseMotion() {
  using std::sin; using std::cos; using std::fabs; 
  using std::atan2; using std::sqrt; using std::pow;
  
  m_arm_model.getDependentFrame3D(0)->mFrame
  frame_3D<double> EE_fr = m_arm_model.getDependentFrame3D(0)->mFrame->getFrameRelativeTo(m_base_frame);
  
  
  const double pos_epsilon = 1e-3;
  const double extend_epsilon = .02;
  
  
  m_arm_model.getBaseToShoulder();
  m_arm_model.getShoulderToElbow();
  m_arm_model.getElbowToJoint4();
  m_arm_model.getJoint4ToWrist();
  m_arm_model.getWristToFlange();
  
  
  quaternion<double>::zrot gl_to_track_rot(-0.5 * M_PI);
  vect<double,3> EE_z_axis = gl_to_track_rot * (EE_fr.Quat * vect_k);
  vect<double,3> EE_y_axis = gl_to_track_rot * (EE_fr.Quat * vect_j);
  vect<double,3> wrist_pos = gl_to_track_rot * EE_fr.Position
                           - m_arm_model.getWristToFlange() * EE_z_axis 
                           - m_arm_model.getBaseToShoulder() * vect_k;
  double elbow_to_wrist_dist = m_arm_model.getElbowToJoint4() + m_arm_model.getJoint4ToWrist();
  double shoulder_to_wrist = m_arm_model.getShoulderToElbow() + elbow_to_wrist_dist;
  double s2e_dist_sqr = m_arm_model.getShoulderToElbow() * m_arm_model.getShoulderToElbow();
  double e2w_dist_sqr = elbow_to_wrist_dist * elbow_to_wrist_dist;
  
  /*
   * find the maximum wrist to base distance, x_max and verifies if the required
   * position of the end-effector is within limits
   */
  /*Extended arm*/
  double c2_max = cos(joint_upper_bounds[1]);
  double s2_max = sin(joint_upper_bounds[1]);
  
  double x_max = 0.0;
  if (wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * c2_max) {
    if (wrist_pos[1] * wrist_pos[1] + wrist_pos[2] * wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * (shoulder_to_wrist - extend_epsilon))
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is outside the cylindrical workspace envelope (fully-extended arm).");
    x_max = sqrt( shoulder_to_wrist * shoulder_to_wrist - wrist_pos[1] * wrist_pos[1] - wrist_pos[2] * wrist_pos[2]);
  } else /*Bent arm*/ {
    // Verifies that the location can be reached, as far as height (z) is concerned
    double max_j23_angle = joint_upper_bounds[1] + joint_upper_bounds[2];
    if(max_j23_angle > M_PI)
      max_j23_angle = M_PI;
    double low_elbow_to_desired_height = wrist_pos[2] - A465_params.shoulder_to_elbow_dist * c2_max;
    double elbow_wrist_eps_dist = elbow_to_wrist_dist - extend_epsilon;
    
    if (low_elbow_to_desired_height < elbow_wrist_eps_dist * cos(max_j23_angle))
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too low for the manipulator to reach.");
    
    elbow_wrist_eps_dist *= elbow_wrist_eps_dist;
    low_elbow_to_desired_height *= low_elbow_to_desired_height;
    
    if ( low_elbow_to_desired_height + pow(fabs(wrist_pos[1]) - A465_params.shoulder_to_elbow_dist * s2_max, 2) > elbow_wrist_eps_dist)
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far below the manipulator.");
    
    double xy_dist = A465_params.shoulder_to_elbow_dist * s2_max + sqrt(elbow_wrist_eps_dist - low_elbow_to_desired_height);
    x_max = sqrt(xy_dist * xy_dist - wrist_pos[1] * wrist_pos[1]);
  };
  
  /* At this point, the range of joint values for the track is between +- x_max around wrist_pos[0] (to within track-limits). */
  
  /* Joint 0 (track) */
  
  // first, check that the limits of the track permit at least some solution:
  if( wrist_pos[0] - x_max > m_track_upper_bound )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far ahead of the track (beyond upper track limit).");
  if( wrist_pos[0] + x_max < m_track_lower_bound )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far behind the track (beyond lower track limit).");
  
  double x_desired = 0.0;
  /* solving for the best track position based on trying to get a right-angle at the elbow */ 
  double R0_rhs = s2e_dist_sqr + e2w_dist_sqr - wrist_pos[1] * wrist_pos[1] - wrist_pos[2] * wrist_pos[2];
  if(R0_rhs > 0) {
    // this means that it is possible to obtain an exact right angle at the elbow.
    if(EE_z_axis[0] >= 0.0)  // the EE needs to be pointing in the forward x direction:
      x_desired = wrist_pos[0] - sqrt(R0_rhs);
    else
      x_desired = wrist_pos[0] + sqrt(R0_rhs);
  } else {
    // this means that the best we can do is put the robot-base in x-alignment with the wrist position.
    x_desired = wrist_pos[0];
  };
  
  // clamp the x-solution to the track's range:
  if( x_desired < m_track_lower_bound )
    x_desired = m_track_lower_bound;
  else if( x_desired > m_track_upper_bound )
    x_desired = m_track_upper_bound;
  // update the wrist-position vector:
  wrist_pos[0] -= x_desired;
  
  /* reach forward solutions */
  x_desired;
  
  
  
  
  quaternion<double>::zrot gl_to_track_rot(-0.5 * M_PI);
  vect<double,3> EE_z_axis = EE_fr.Quat * vect_k;
  vect<double,3> EE_y_axis = EE_fr.Quat * vect_j;
  vect<double,3> wrist_pos = EE_fr.Position - wrist_to_flange * EE_z_axis - base_to_shoulder * vect_k;
  
  double elbow_to_wrist = elbow_to_joint_4 + joint_4_to_wrist;
  double shoulder_to_wrist = shoulder_to_elbow + elbow_to_wrist;
  double s2e_dist_sqr = shoulder_to_elbow * shoulder_to_elbow;
  double e2w_dist_sqr = elbow_to_wrist * elbow_to_wrist;
  
  const double pos_epsilon = 1e-4 * shoulder_to_wrist;
  const double extend_epsilon = .01 * shoulder_to_wrist;
  
  vect_n<double> solns[8];
  for(std::size_t i = 0; i < 8; ++i)
    solns[i].resize(6,0.0);
  
  /*
   * Verify if the required position of the end-effector is within the arm's reach.
   */
  /*Extended arm*/
  double c2_max, s2_max;
  if(fabs(joint_lower_bounds[1]) > fabs(joint_upper_bounds[1])) {
    c2_max = cos(joint_lower_bounds[1]);
    s2_max = sin(fabs(joint_lower_bounds[1]));
  } else {
    c2_max = cos(joint_upper_bounds[1]);
    s2_max = sin(fabs(joint_upper_bounds[1]));
  };
  
  if (wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * c2_max) {
    if (wrist_pos[0] * wrist_pos[0] + wrist_pos[1] * wrist_pos[1] + wrist_pos[2] * wrist_pos[2] > (shoulder_to_wrist - extend_epsilon) * (shoulder_to_wrist - extend_epsilon))
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is outside the spherical workspace envelope (fully-extended arm).");
  } else /*Bent arm*/ {
    // Verifies that the location can be reached, as far as height (z) is concerned
    double c23_max = -1.0;
    double max_j23_angle = joint_upper_bounds[1] + joint_upper_bounds[2];
    if(fabs(joint_lower_bounds[1] + joint_lower_bounds[2]) > max_j23_angle)
      max_j23_angle = fabs(joint_lower_bounds[1] + joint_lower_bounds[2]);
    if(max_j23_angle < M_PI)
      c23_max = cos(max_j23_angle);
    double low_elbow_to_desired_height = wrist_pos[2] - shoulder_to_elbow * c2_max;
    double elbow_wrist_eps = elbow_to_wrist - extend_epsilon;
    
    if (low_elbow_to_desired_height < elbow_wrist_eps * c23_max)
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too low for the manipulator to reach.");
    
    if ( low_elbow_to_desired_height * low_elbow_to_desired_height + pow(sqrt(wrist_pos[0] * wrist_pos[0] + wrist_pos[1] * wrist_pos[1]) - shoulder_to_elbow * s2_max, 2) > elbow_wrist_eps * elbow_wrist_eps)
      throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Desired wrist position is too far below the manipulator.");
  };
  
  /* Joint 1 */
  
  if( (fabs(wrist_pos[0]) < pos_epsilon) && (fabs(wrist_pos[1]) < pos_epsilon) ) {
    /* we're in the joint 1 singularity directly above the origin */
    solns[fun_posture][0] = m_joints[0]->q;
    solns[bun_posture][0] = clamp_to_pi_range(solns[fun_posture][0] + M_PI);
  } else {
    solns[fun_posture][0] = atan2(wrist_pos[1], wrist_pos[0]);
    solns[bun_posture][0] = clamp_to_pi_range(solns[fun_posture][0] + M_PI);
  };
  
  /* set up some variables for later */
  double s1_F = sin(solns[fun_posture][0]); /* forward */
  double c1_F = cos(solns[fun_posture][0]);
  
  /* Joint 3 - modified version (MP: just seems like a more straight forward way to do it) */ 
  
  double baseplane_dist_sqr = wrist_pos[0] * wrist_pos[0] + wrist_pos[1] * wrist_pos[1];
  double wrist_dist_sqr = baseplane_dist_sqr + wrist_pos[2] * wrist_pos[2];
  double j3tmp0 = s2e_dist_sqr + e2w_dist_sqr - wrist_dist_sqr;
  double j3tmp1 = 4.0 * s2e_dist_sqr * e2w_dist_sqr - j3tmp0 * j3tmp0;
    
  /* ensure we're within reach */
  if( j3tmp1 < 0.0 )
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! Cannot compute an elbow angle that would reach the desired wrist position.");
  
  /* now determine the Joint 3 angle (FUN) */
  solns[fun_posture][2] = -atan2(sqrt(j3tmp1), -j3tmp0);
  solns[bun_posture][2] = -solns[fun_posture][2];
  
  double inv_wrist_dist_sqr = 1.0 / wrist_dist_sqr;
  double baseplane_dist = sqrt(baseplane_dist_sqr);
  double s2e_projection = shoulder_to_elbow * sin(solns[fun_posture][2]);
  
  vect<double,2> to_wrist_along_l1(shoulder_to_elbow * cos(solns[fun_posture][2]) + elbow_to_wrist, s2e_projection);
  vect<double,2> to_wrist_global(baseplane_dist, wrist_pos[2]); // reach-forward vector.
  
  
  /* register all solutions for joints 1 and 3 */
  
  /* reach forward solutions */
  solns[fuf_posture][0] = solns[fun_posture][0];
  solns[fdn_posture][0] = solns[fun_posture][0];
  solns[fdf_posture][0] = solns[fun_posture][0];
  
  /* reach backward solutions */
  solns[buf_posture][0] = solns[bun_posture][0];
  solns[bdn_posture][0] = solns[bun_posture][0];
  solns[bdf_posture][0] = solns[bun_posture][0];
  
  /* backward/up and forward/down solutions */
  solns[buf_posture][2] = solns[bun_posture][2];
  solns[fdn_posture][2] = solns[bun_posture][2];
  solns[fdf_posture][2] = solns[bun_posture][2];
  
  /* foward/up and backward/down solutions */
  solns[fuf_posture][2] = solns[fun_posture][2];
  solns[bdn_posture][2] = solns[fun_posture][2];
  solns[bdf_posture][2] = solns[fun_posture][2];
  
  /* Prepare some variables that are constant throughout the iterations. */
  vect<double,2> EE_z_proj_F(c1_F * EE_z_axis[0] + s1_F * EE_z_axis[1], 
                            -s1_F * EE_z_axis[0] + c1_F * EE_z_axis[1]);
  vect<double,2> EE_y_proj_F(c1_F * EE_y_axis[0] + s1_F * EE_y_axis[1], 
                            -s1_F * EE_y_axis[0] + c1_F * EE_y_axis[1]);
  
  // solve for both elbow-cases (up or down).
  for(std::size_t i = 0; i < 2; ++i) {
    std::size_t e_o = 2 * i;
    double s2e_projection_i = (i ? -s2e_projection : s2e_projection);
    double c23 = (to_wrist_along_l1[0] * wrist_pos[2] + s2e_projection_i * baseplane_dist) * inv_wrist_dist_sqr;
    double s23 = /*reach_sign * */( s2e_projection_i * wrist_pos[2] - to_wrist_along_l1[0] * baseplane_dist) * inv_wrist_dist_sqr;
    
    solns[fuf_posture + e_o][1] = (solns[fun_posture + e_o][1] = atan2( s23, c23) - solns[fun_posture + e_o][2]);
    solns[buf_posture + e_o][1] = (solns[bun_posture + e_o][1] = atan2(-s23, c23) - solns[bun_posture + e_o][2]);
    
    double jt5_i = s23 * EE_z_axis[2] + c23 * EE_z_proj_F[0];
    
    /* first check if we're in the J5 singularity (i.e. J5 ~= 0)  */
    if( (fabs(jt5_i) < pos_epsilon) && (fabs(EE_z_proj_F[1]) < pos_epsilon) ) {
      double c4 = /*reach_sign * */cos(m_joints[3]->q);
      double s4 = /*reach_sign * */sin(m_joints[3]->q);
      double a4 = m_joints[3]->q;  /* F*N or B*F */
      
      solns[buf_posture + e_o][3] = (solns[fun_posture + e_o][3] = a4);
      solns[bun_posture + e_o][3] = (solns[fuf_posture + e_o][3] = clamp_to_pi_range(a4 + M_PI));
      
      solns[bun_posture + e_o][4] = (solns[fun_posture + e_o][4] = 0.0);
      solns[buf_posture + e_o][4] = (solns[fuf_posture + e_o][4] = 0.0);
      
      double c6 = -s4 * (c23 * EE_y_proj_F[0] + s23 * EE_y_axis[2]) + c4 * EE_y_proj_F[1];
      double s6 = -s4 * EE_y_proj_F[1] - c4 * (EE_y_axis[2] * s23 + EE_y_proj_F[0] * c23);
      
      double a6 = atan2(s6, c6);  /* wrist not flipped */
      solns[bun_posture + e_o][5] = (solns[fun_posture + e_o][5] = a6);
      solns[buf_posture + e_o][5] = (solns[fuf_posture + e_o][5] = clamp_to_pi_range(a6 + M_PI));
    } else {
      /* we're not singular in jt 5 */
      double a4 = atan2(EE_z_proj_F[1], jt5_i);  /* F*N or B*F */
      double s4 = /* reach_sign *  */ sin(a4);
      double c4 = /* reach_sign *  */ cos(a4);
      
      solns[buf_posture + e_o][3] = (solns[fun_posture + e_o][3] = a4);
      solns[bun_posture + e_o][3] = (solns[fuf_posture + e_o][3] = clamp_to_pi_range(a4 + M_PI));
      
      double c5 =  EE_z_axis[2] * c23 - EE_z_proj_F[0] * s23;
      double s5 = -c4 * jt5_i - s4 * EE_z_proj_F[1]; 
      
      double a5 = atan2(s5, c5);  /* wrist not flipped */
      solns[bun_posture + e_o][4] = (solns[fun_posture + e_o][4] =  a5);
      solns[buf_posture + e_o][4] = (solns[fuf_posture + e_o][4] = -a5);
      
      double c6 = -s4 * (c23 * EE_y_proj_F[0] + s23 * EE_y_axis[2])
                 + c4 * EE_y_proj_F[1];
      double s6 =  s4 * (EE_z_proj_F[1] * (EE_y_axis[2] * c23 - EE_y_proj_F[0] * s23) - c5 * EE_y_proj_F[1])
                 + c4 * (EE_y_axis[2] * EE_z_proj_F[0] - EE_y_proj_F[0] * EE_z_axis[2]);
      
      double a6 = atan2(s6, c6);  /* wrist not flipped */
      solns[bun_posture + e_o][5] = (solns[fun_posture + e_o][5] = a6);
      solns[buf_posture + e_o][5] = (solns[fuf_posture + e_o][5] = clamp_to_pi_range(a6 + M_PI));
    };
  };
  
  
  // Now, we just need to choose a best solution. 
  // We'll choose it based on the previous set of joint coordinates (whichever coordinates are currently in the joint pointers):
  double      best_cost = 100.0;
  std::size_t best_soln = 0;
  for(std::size_t i = 0; i < 8; ++i) {
    double cost = 0.0;
    for(std::size_t j = 0; j < 6; ++j) {
      if( ( solns[i][j] > joint_lower_bounds[j] ) && ( solns[i][j] < joint_upper_bounds[j] ) ) {
        cost += fabs( solns[i][j] - m_joints[j]->q ) / (joint_upper_bounds[j] - joint_lower_bounds[j]);
//         cost += fabs( 2.0 * solns[i][j] - joint_lower_bounds[j] - joint_upper_bounds[j] ) / (joint_upper_bounds[j] - joint_lower_bounds[j]);
      } else {
        cost = 100.0; // effectively makes this solution, one of the worst solutions (outside the bounds).
        break;
      };
    };
    if(cost < best_cost) {
      best_cost = cost;
      best_soln = i;
    };
  };
  if(best_cost > 99.0)
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! None of the inverse kinematics solutions respect the joint limits.");
  
  // set the joint coordinates to the computed values:
  m_joints[0]->q = solns[best_soln][0];
  m_joints[1]->q = solns[best_soln][1];
  m_joints[2]->q = solns[best_soln][2];
  m_joints[3]->q = solns[best_soln][3];
  m_joints[4]->q = solns[best_soln][4];
  m_joints[5]->q = solns[best_soln][5];
  
  // then, use the solution to compute the jacobian matrix:
  mat<double,mat_structure::rectangular> jac(6,6);
  getJacobianMatrix(jac);
  
  // finally, use the jacobian to find a solution for the joint velocities:
  mat<double,mat_structure::rectangular> x(6,1);
  mat<double,mat_structure::rectangular> b(6,1);
  b(0,0) = EE_fr.Velocity[0];    b(1,0) = EE_fr.Velocity[1];    b(2,0) = EE_fr.Velocity[2];
  b(3,0) = EE_fr.AngVelocity[0]; b(4,0) = EE_fr.AngVelocity[1]; b(5,0) = EE_fr.AngVelocity[2];
  linlsq_QR(jac, x, b, pos_epsilon);
  
  m_joints[0]->q_dot = x(0,0);
  m_joints[1]->q_dot = x(1,0);
  m_joints[2]->q_dot = x(2,0);
  m_joints[3]->q_dot = x(3,0);
  m_joints[4]->q_dot = x(4,0);
  m_joints[5]->q_dot = x(5,0);
  // acceleration is irrelevant (not part of start variables).
  m_joints[0]->q_ddot = 0.0;
  m_joints[1]->q_ddot = 0.0;
  m_joints[2]->q_ddot = 0.0;
  m_joints[3]->q_ddot = 0.0;
  m_joints[4]->q_ddot = 0.0;
  m_joints[5]->q_ddot = 0.0;
  
  m_chain->doMotion();
  
};

void manip_P3R3R_kinematics::getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const {
  /* calculate individual rotations */
  quaternion<double>::zrot q1( m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::yrot q3(-m_joints[2]->q);
  quaternion<double>::zrot q4( m_joints[3]->q);
  quaternion<double>::yrot q5(-m_joints[4]->q);
  quaternion<double>::zrot q6( m_joints[5]->q);
  
  quaternion<double> q_accum = q1.getQuaternion();
  vect<double,3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double,3> a2 = shoulder_to_elbow * (q_accum * vect_k);
  q_accum *= q3;
  vect<double,3> a3 = (elbow_to_joint_4 + joint_4_to_wrist) * (q_accum * vect_k);
  vect<double,3> e4 =  q_accum * vect_k;
  q_accum *= q4;
  vect<double,3> e5 = -(q_accum * vect_j);
  q_accum *= q5;
  vect<double,3> e6 =  q_accum * vect_k;
  q_accum *= q6;
  vect<double,3> a6 = wrist_to_flange * e6; 
  
  vect<double,3> s2f = a2 + a3 + a6; // shoulder to flange
  
  Jac.resize(std::make_pair(6,6));
  vect<double,3> v1 = vect_k % s2f;
  Jac(0,0) = v1[0];
  Jac(1,0) = v1[1];
  Jac(2,0) = v1[2];
  Jac(3,0) = 0.0;
  Jac(4,0) = 0.0;
  Jac(5,0) = 1.0;
  vect<double,3> v2 = e2 % s2f;
  Jac(0,1) = v2[0];
  Jac(1,1) = v2[1];
  Jac(2,1) = v2[2];
  Jac(3,1) = e2[0];
  Jac(4,1) = e2[1];
  Jac(5,1) = e2[2];
  vect<double,3> v3 = e2 % (s2f - a2);
  Jac(0,2) = v3[0];
  Jac(1,2) = v3[1];
  Jac(2,2) = v3[2];
  Jac(3,2) = e2[0];
  Jac(4,2) = e2[1];
  Jac(5,2) = e2[2];
  vect<double,3> v4 = e4 % a6;
  Jac(0,3) = v4[0];
  Jac(1,3) = v4[1];
  Jac(2,3) = v4[2];
  Jac(3,3) = e4[0];
  Jac(4,3) = e4[1];
  Jac(5,3) = e4[2];
  vect<double,3> v5 = e5 % a6;
  Jac(0,4) = v5[0];
  Jac(1,4) = v5[1];
  Jac(2,4) = v5[2];
  Jac(3,4) = e5[0];
  Jac(4,4) = e5[1];
  Jac(5,4) = e5[2];
  Jac(0,5) = 0.0;
  Jac(1,5) = 0.0;
  Jac(2,5) = 0.0;
  Jac(3,5) = e6[0];
  Jac(4,5) = e6[1];
  Jac(5,5) = e6[2];
};

void manip_P3R3R_kinematics::getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const {
  /* calculate individual rotations */
  quaternion<double>::zrot q1( m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::yrot q3(-m_joints[2]->q);
  quaternion<double>::zrot q4( m_joints[3]->q);
  quaternion<double>::yrot q5(-m_joints[4]->q);
  quaternion<double>::zrot q6( m_joints[5]->q);
  
  quaternion<double> q_accum = q1.getQuaternion();
  vect<double,3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double,3> a2 = shoulder_to_elbow * (q_accum * vect_k);
  q_accum *= q3;
  vect<double,3> a3 = (elbow_to_joint_4 + joint_4_to_wrist) * (q_accum * vect_k);
  vect<double,3> e4 =  q_accum * vect_k;
  q_accum *= q4;
  vect<double,3> e5 = -(q_accum * vect_j);
  q_accum *= q5;
  vect<double,3> e6 =  q_accum * vect_k;
  q_accum *= q6;
  vect<double,3> a6 = wrist_to_flange * e6; 
  
  vect<double,3> s2f = a2 + a3 + a6; // shoulder to flange
  
  Jac.resize(std::make_pair(6,6));
  vect<double,3> v1 = vect_k % s2f;
  Jac(0,0) = v1[0];
  Jac(1,0) = v1[1];
  Jac(2,0) = v1[2];
  Jac(3,0) = 0.0;
  Jac(4,0) = 0.0;
  Jac(5,0) = 1.0;
  vect<double,3> v2 = e2 % s2f;
  Jac(0,1) = v2[0];
  Jac(1,1) = v2[1];
  Jac(2,1) = v2[2];
  Jac(3,1) = e2[0];
  Jac(4,1) = e2[1];
  Jac(5,1) = e2[2];
  vect<double,3> v3 = e2 % (s2f - a2);
  Jac(0,2) = v3[0];
  Jac(1,2) = v3[1];
  Jac(2,2) = v3[2];
  Jac(3,2) = e2[0];
  Jac(4,2) = e2[1];
  Jac(5,2) = e2[2];
  vect<double,3> v4 = e4 % a6;
  Jac(0,3) = v4[0];
  Jac(1,3) = v4[1];
  Jac(2,3) = v4[2];
  Jac(3,3) = e4[0];
  Jac(4,3) = e4[1];
  Jac(5,3) = e4[2];
  vect<double,3> v5 = e5 % a6;
  Jac(0,4) = v5[0];
  Jac(1,4) = v5[1];
  Jac(2,4) = v5[2];
  Jac(3,4) = e5[0];
  Jac(4,4) = e5[1];
  Jac(5,4) = e5[2];
  Jac(0,5) = 0.0;
  Jac(1,5) = 0.0;
  Jac(2,5) = 0.0;
  Jac(3,5) = e6[0];
  Jac(4,5) = e6[1];
  Jac(5,5) = e6[2];
  
  
  //vect<double,3> q1_dot =  m_joints[0]->q_dot * vect_k;
  //vect<double,3> q2_dot = -m_joints[1]->q_dot * vect_j;
  //vect<double,3> q3_dot = -m_joints[2]->q_dot * vect_j;
  //vect<double,3> q4_dot =  m_joints[3]->q_dot * vect_k;
  //vect<double,3> q5_dot = -m_joints[4]->q_dot * vect_j;
  //vect<double,3> q6_dot =  m_joints[5]->q_dot * vect_k;
  vect<double,3> vi(1.0,0.0,0.0);
  vect<double,3> vj(0.0,1.0,0.0);
  vect<double,3> vk(0.0,0.0,1.0);
  
  vect<double,3> e2_dot =  q1 * ( m_joints[0]->q_dot * vi);
  
  vect<double,3> a2_dot = shoulder_to_elbow * (q1 * (m_joints[0]->q_dot * (vect_k % (q2 * vk))
                                             + q2 * (-m_joints[1]->q_dot * vi)));
  
  vect<double,3> e4_dot =  q1 * ( m_joints[0]->q_dot * (vect_k % (q2 * q3 * vk))
                         - q2 * ( m_joints[1]->q_dot * (vect_j % (q3 * vk))
                         + q3 * ( m_joints[2]->q_dot * vi)));
  vect<double,3> a3_dot = (elbow_to_joint_4 + joint_4_to_wrist) * e4_dot;

  vect<double,3> e5_dot =  q1 * (-m_joints[0]->q_dot * (vect_k % (q2 * q3 * q4 * vj))
                         + q2 * ( m_joints[1]->q_dot * (vect_j % (q3 * q4 * vj))
                         + q3 * ( m_joints[2]->q_dot * (vect_j % (q4 * vj))
                         + q4 * ( m_joints[3]->q_dot * vi))));
  
  vect<double,3> e6_dot =  q1 * ( m_joints[0]->q_dot * (vect_k % (q2 * q3 * q4 * q5 * vk))
                         - q2 * ( m_joints[1]->q_dot * (vect_j % (q3 * q4 * q5 * vk))
                         + q3 * ( m_joints[2]->q_dot * (vect_j % (q4 * q5 * vk))
                         - q4 * ( m_joints[3]->q_dot * (vect_k % (q5 * vk))
                         - q5 * ( m_joints[4]->q_dot * vi)))));
  
  vect<double,3> a6_dot = wrist_to_flange * e6_dot; 
  vect<double,3> s2f_dot = a2_dot + a3_dot + a6_dot; // shoulder to flange
  
  JacDot.resize(std::make_pair(6,6));
  vect<double,3> v1_dot = vect_k % s2f_dot;
  JacDot(0,0) = v1_dot[0];
  JacDot(1,0) = v1_dot[1];
  JacDot(2,0) = v1_dot[2];
  JacDot(3,0) = 0.0;
  JacDot(4,0) = 0.0;
  JacDot(5,0) = 0.0;
  vect<double,3> v2_dot = e2_dot % s2f + e2 % s2f_dot;
  JacDot(0,1) = v2_dot[0];
  JacDot(1,1) = v2_dot[1];
  JacDot(2,1) = v2_dot[2];
  JacDot(3,1) = e2_dot[0];
  JacDot(4,1) = e2_dot[1];
  JacDot(5,1) = e2_dot[2];
  vect<double,3> v3_dot = e2_dot % (s2f - a2) + e2 % (s2f_dot - a2_dot);
  JacDot(0,2) = v3_dot[0];
  JacDot(1,2) = v3_dot[1];
  JacDot(2,2) = v3_dot[2];
  JacDot(3,2) = e2_dot[0];
  JacDot(4,2) = e2_dot[1];
  JacDot(5,2) = e2_dot[2];
  vect<double,3> v4_dot = e4_dot % a6 + e4 % a6_dot;
  JacDot(0,3) = v4_dot[0];
  JacDot(1,3) = v4_dot[1];
  JacDot(2,3) = v4_dot[2];
  JacDot(3,3) = e4_dot[0];
  JacDot(4,3) = e4_dot[1];
  JacDot(5,3) = e4_dot[2];
  vect<double,3> v5_dot = e5_dot % a6 + e5 % a6_dot;
  JacDot(0,4) = v5_dot[0];
  JacDot(1,4) = v5_dot[1];
  JacDot(2,4) = v5_dot[2];
  JacDot(3,4) = e5_dot[0];
  JacDot(4,4) = e5_dot[1];
  JacDot(5,4) = e5_dot[2];
  JacDot(0,5) = 0.0;
  JacDot(1,5) = 0.0;
  JacDot(2,5) = 0.0;
  JacDot(3,5) = e6_dot[0];
  JacDot(4,5) = e6_dot[1];
  JacDot(5,5) = e6_dot[2];
};

vect_n<double> manip_P3R3R_kinematics::getJointPositions() const {
  vect_n<double> arm_jtpos = m_arm_model.getJointPositions();
  return vect_n<double>(m_track_coord->q, arm_jtpos[0], arm_jtpos[1], arm_jtpos[2], arm_jtpos[3], arm_jtpos[4], arm_jtpos[5]);
};

void manip_P3R3R_kinematics::setJointPositions(const vect_n<double>& aJointPositions) {
  m_track_coord->q = aJointPositions[0];
  m_arm_model.setJointPositions(vect_n<double>(aJointPositions.begin()+1, aJointPositions.end()));
};

vect_n<double> manip_P3R3R_kinematics::getJointVelocities() const {
  vect_n<double> arm_jtvel = m_arm_model.getJointVelocities();
  return vect_n<double>(m_track_coord->q_dot, arm_jtvel[0], arm_jtvel[1], arm_jtvel[2], arm_jtvel[3], arm_jtvel[4], arm_jtvel[5]);
};

void manip_P3R3R_kinematics::setJointVelocities(const vect_n<double>& aJointVelocities) {
  m_track_coord->q_dot = aJointVelocities[0];
  m_arm_model.setJointVelocities(vect_n<double>(aJointVelocities.begin()+1, aJointVelocities.end()));
};

vect_n<double> manip_P3R3R_kinematics::getJointAccelerations() const {
  vect_n<double> arm_jtacc = m_arm_model.getJointAccelerations();
  return vect_n<double>(m_track_coord->q_ddot, arm_jtacc[0], arm_jtacc[1], arm_jtacc[2], arm_jtacc[3], arm_jtacc[4], arm_jtacc[5]);
};

void manip_P3R3R_kinematics::setJointAccelerations(const vect_n<double>& aJointAccelerations) {
  m_track_coord->q_ddot = aJointAccelerations[0];
  m_arm_model.setJointAccelerations(vect_n<double>(aJointAccelerations.begin()+1, aJointAccelerations.end()));
};



void RK_CALL manip_P3R3R_kinematics::save(serialization::oarchive& A, unsigned int) const {
  inverse_kinematics_model::save(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(m_base_frame)
    & RK_SERIAL_SAVE_WITH_NAME(m_track_coord)
    & RK_SERIAL_SAVE_WITH_NAME(m_output_frame)
    & RK_SERIAL_SAVE_WITH_NAME(m_track_direction)
    & RK_SERIAL_SAVE_WITH_NAME(m_track_lower_bound)
    & RK_SERIAL_SAVE_WITH_NAME(m_track_upper_bound)
    & RK_SERIAL_SAVE_WITH_NAME(m_arm_model)
    & RK_SERIAL_SAVE_WITH_NAME(m_chain);
};

void RK_CALL manip_P3R3R_kinematics::load(serialization::iarchive& A, unsigned int) {
  inverse_kinematics_model::load(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(m_base_frame)
    & RK_SERIAL_LOAD_WITH_NAME(m_track_coord)
    & RK_SERIAL_LOAD_WITH_NAME(m_output_frame)
    & RK_SERIAL_LOAD_WITH_NAME(m_track_direction)
    & RK_SERIAL_LOAD_WITH_NAME(m_track_lower_bound)
    & RK_SERIAL_LOAD_WITH_NAME(m_track_upper_bound)
    & RK_SERIAL_LOAD_WITH_NAME(m_arm_model)
    & RK_SERIAL_LOAD_WITH_NAME(m_chain);
};


};

};








