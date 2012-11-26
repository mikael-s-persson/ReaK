
/*
 *    Copyright 2012 Sven Mikael Persson
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

#include "manip_SSRMS_arm.hpp"

#include "mbd_kte/revolute_joint.hpp"
#include "mbd_kte/rigid_link.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_qr_decomp.hpp"
#include "optimization/optim_exceptions.hpp"
#include <cmath>

namespace ReaK {

namespace kte {


//     shared_ptr< frame_3D<double> > m_base_frame;
//     std::vector< shared_ptr< gen_coord<double> > > m_joints;
//     shared_ptr< joint_dependent_frame_3D > m_EE;
//     vect_n<double> link_lengths;
//     vect_n<double> joint_offsets;
//     vect_n<double> joint_lower_bounds;
//     vect_n<double> joint_upper_bounds;

manip_SSRMS_kinematics::manip_SSRMS_kinematics(const std::string& aName,
                                               const shared_ptr< frame_3D<double> >& aBaseFrame,
                                               const vect_n<double>& aLinkLengths, 
                                               const vect_n<double>& aJointOffsets,
                                               const vect_n<double>& aJointLowerBounds,
                                               const vect_n<double>& aJointUpperBounds) :
                                               inverse_kinematics_model(aName),
                                               m_base_frame(aBaseFrame),
                                               link_lengths(aLinkLengths), 
                                               joint_offsets(aJointOffsets),
                                               joint_lower_bounds(aJointLowerBounds),
                                               joint_upper_bounds(aJointUpperBounds) {
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  
  
  //declare all the intermediate frames.
  shared_ptr< frame_3D<double> > joint_1_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_1_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_2_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_2_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_3_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_3_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_4_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_4_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_5_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_5_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_6_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_6_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_7_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_7_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > arm_EE(new frame_3D<double>(), scoped_deleter());
  
  //declare all the joint jacobians.
  shared_ptr< jacobian_gen_3D<double> > joint_1_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_2_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_3_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_4_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_5_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_6_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_7_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  
  
  //create link 
  shared_ptr< rigid_link_3D > link_0(new rigid_link_3D(
      "manip_SSRMS_link_0",
      m_base_frame,
      joint_1_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, joint_offsets[0]),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< revolute_joint_3D > joint_1(new revolute_joint_3D(
      "manip_SSRMS_joint_1",
      m_joints[0],
      vect<double,3>(0.0,0.0,1.0),
      joint_1_base,
      joint_1_end,
      joint_1_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_1(new rigid_link_3D(
      "manip_SSRMS_link_1",
      joint_1_end,
      joint_2_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[0], -joint_offsets[1], 0.0),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< revolute_joint_3D > joint_2(new revolute_joint_3D(
      "manip_SSRMS_joint_2",
      m_joints[1],
      vect<double,3>(0.0,-1.0,0.0),
      joint_2_base,
      joint_2_end,
      joint_2_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_2(new rigid_link_3D(
      "manip_SSRMS_link_2",
      joint_2_end,
      joint_3_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[1], 0.0, joint_offsets[2]),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< revolute_joint_3D > joint_3(new revolute_joint_3D(
      "manip_SSRMS_joint_3",
      m_joints[2],
      vect<double,3>(0.0,0.0,1.0),
      joint_3_base,
      joint_3_end,
      joint_3_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_3(new rigid_link_3D(
      "manip_SSRMS_link_3",
      joint_3_end,
      joint_4_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[2], 0.0, joint_offsets[3]),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< revolute_joint_3D > joint_4(new revolute_joint_3D(
      "manip_SSRMS_joint_4",
      m_joints[3],
      vect<double,3>(0.0,0.0,1.0),
      joint_4_base,
      joint_4_end,
      joint_4_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_4(new rigid_link_3D(
      "manip_SSRMS_link_4",
      joint_4_end,
      joint_5_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[3], 0.0, joint_offsets[4]),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< revolute_joint_3D > joint_5(new revolute_joint_3D(
      "manip_SSRMS_joint_5",
      m_joints[4],
      vect<double,3>(0.0,0.0,1.0),
      joint_5_base,
      joint_5_end,
      joint_5_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_5(new rigid_link_3D(
      "manip_SSRMS_link_5",
      joint_5_end,
      joint_6_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[4], 0.0, 0.0),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< kte::revolute_joint_3D > joint_6(new kte::revolute_joint_3D(
      "manip_SSRMS_joint_6",
      m_joints[5],
      vect<double,3>(0.0,-1.0,0.0),
      joint_6_base,
      joint_6_end,
      joint_6_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_6(new rigid_link_3D(
      "manip_SSRMS_link_6",
      joint_6_end,
      joint_7_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(link_lengths[5], -joint_offsets[5], 0.0),
        quaternion<double>::xrot(M_PI).getQuaternion())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< kte::revolute_joint_3D > joint_7(new kte::revolute_joint_3D(
      "manip_SSRMS_joint_7",
      m_joints[5],
      vect<double,3>(0.0,0.0,1.0),
      joint_7_base,
      joint_7_end,
      joint_7_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< rigid_link_3D > link_7(new rigid_link_3D(
      "manip_SSRMS_link_7",
      joint_7_end,
      arm_EE,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, joint_offsets[6]),
        quaternion<double>())),
    scoped_deleter());
  
  //create inertia
  m_EE = shared_ptr< joint_dependent_frame_3D >(new joint_dependent_frame_3D(arm_EE), scoped_deleter());
  
  m_EE->add_joint(m_joints[0], joint_1_jacobian);
  m_EE->add_joint(m_joints[1], joint_2_jacobian);
  m_EE->add_joint(m_joints[2], joint_3_jacobian);
  m_EE->add_joint(m_joints[3], joint_4_jacobian);
  m_EE->add_joint(m_joints[4], joint_5_jacobian);
  m_EE->add_joint(m_joints[5], joint_6_jacobian);
  m_EE->add_joint(m_joints[6], joint_7_jacobian);
  
  m_chain = shared_ptr< kte_map_chain >(new kte_map_chain("manip_SSRMS_kin_model"), scoped_deleter());
  
  *m_chain << link_0
           << joint_1
           << link_1
           << joint_2
           << link_2
           << joint_3
           << link_3
           << joint_4
           << link_4
           << joint_5
           << link_5
           << joint_6
           << link_6
           << joint_7
           << link_7;
  
};



void manip_SSRMS_kinematics::doDirectMotion() {
  m_chain->doMotion();
};



static double clamp_to_pi_range(double a) {
  return (a > M_PI ? a - 2.0 * M_PI : (a < -M_PI ? a + 2.0 * M_PI : a) );
};

static const int fun_posture = 0;
static const int fuf_posture = 1;
static const int fdn_posture = 2;
static const int fdf_posture = 3;
static const int bun_posture = 4;
static const int buf_posture = 5;
static const int bdn_posture = 6;
static const int bdf_posture = 7;


void manip_SSRMS_kinematics::doInverseMotion() {
#if 0
  using std::sin; using std::cos; using std::fabs; 
  using std::atan2; using std::sqrt; using std::pow;
  
  frame_3D<double> EE_fr = m_EE->mFrame->getFrameRelativeTo(m_base_frame);
  
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
#endif
};

void manip_SSRMS_kinematics::getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const {
  /* calculate individual rotations */
  quaternion<double>::zrot q1( m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::zrot q3( m_joints[2]->q);
  quaternion<double>::zrot q4( m_joints[3]->q);
  quaternion<double>::zrot q5( m_joints[4]->q);
  quaternion<double>::yrot q6(-m_joints[5]->q);
  quaternion<double>::xrot q67( M_PI);
  quaternion<double>::zrot q7( m_joints[6]->q);
  
  vect<double,3> a0(0.0, 0.0, joint_offsets[0]);
  quaternion<double> q_accum = q1.getQuaternion();
  vect<double,3> e1(0.0, 0.0, 1.0);
  vect<double,3> a1 = q_accum * vect<double,3>(link_lengths[0], -joint_offsets[1], 0.0);
  vect<double,3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double,3> a2 = q_accum * vect<double,3>(link_lengths[1], 0.0, joint_offsets[2]);
  vect<double,3> e3 = q_accum * vect_k;
  q_accum *= q3;
  vect<double,3> a3 = q_accum * vect<double,3>(link_lengths[2], 0.0, joint_offsets[3]);
  vect<double,3> e4 = q_accum * vect_k;
  q_accum *= q4;
  vect<double,3> a4 = q_accum * vect<double,3>(link_lengths[3], 0.0, joint_offsets[4]);
  vect<double,3> e5 = q_accum * vect_k;
  q_accum *= q5;
  vect<double,3> a5 = q_accum * vect<double,3>(link_lengths[4], 0.0, 0.0);
  vect<double,3> e6 = -(q_accum * vect_j);
  q_accum *= q6;
  vect<double,3> a6 = q_accum * vect<double,3>(link_lengths[5], -joint_offsets[5], 0.0);
  q_accum *= q67;
  vect<double,3> e7 = q_accum * vect_k;
  q_accum *= q7;
  vect<double,3> a7 = q_accum * vect<double,3>(0.0, 0.0, joint_offsets[6]);
  
  vect<double,3> a67 = a6 + a7;
  vect<double,3> a57 = a5 + a67;
  vect<double,3> a47 = a4 + a57;
  vect<double,3> a37 = a3 + a47;
  vect<double,3> a27 = a2 + a37;
  vect<double,3> a17 = a1 + a27;
  
  Jac.resize(std::make_pair(6,7));
  vect<double,3> v1 = vect_k % a17;
  Jac(0,0) = v1[0];
  Jac(1,0) = v1[1];
  Jac(2,0) = v1[2];
  Jac(3,0) = 0.0;
  Jac(4,0) = 0.0;
  Jac(5,0) = 1.0;
  vect<double,3> v2 = e2 % a27;
  Jac(0,1) = v2[0];
  Jac(1,1) = v2[1];
  Jac(2,1) = v2[2];
  Jac(3,1) = e2[0];
  Jac(4,1) = e2[1];
  Jac(5,1) = e2[2];
  vect<double,3> v3 = e3 % a37;
  Jac(0,2) = v3[0];
  Jac(1,2) = v3[1];
  Jac(2,2) = v3[2];
  Jac(3,2) = e3[0];
  Jac(4,2) = e3[1];
  Jac(5,2) = e3[2];
  vect<double,3> v4 = e4 % a47;
  Jac(0,3) = v4[0];
  Jac(1,3) = v4[1];
  Jac(2,3) = v4[2];
  Jac(3,3) = e4[0];
  Jac(4,3) = e4[1];
  Jac(5,3) = e4[2];
  vect<double,3> v5 = e5 % a57;
  Jac(0,4) = v5[0];
  Jac(1,4) = v5[1];
  Jac(2,4) = v5[2];
  Jac(3,4) = e5[0];
  Jac(4,4) = e5[1];
  Jac(5,4) = e5[2];
  vect<double,3> v6 = e6 % a67;
  Jac(0,5) = v6[0];
  Jac(1,5) = v6[1];
  Jac(2,5) = v6[2];
  Jac(3,5) = e6[0];
  Jac(4,5) = e6[1];
  Jac(5,5) = e6[2];
  Jac(0,6) = 0.0;
  Jac(1,6) = 0.0;
  Jac(2,6) = 0.0;
  Jac(3,6) = e7[0];
  Jac(4,6) = e7[1];
  Jac(5,6) = e7[2];
};

void manip_SSRMS_kinematics::getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const {
  /* calculate individual rotations */
  vect<double,3> ey(0.0,-1.0, 0.0);
  vect<double,3> ez(0.0, 0.0, 1.0);
  
  quaternion<double>::zrot q1( m_joints[0]->q);
  quaternion<double>::yrot q2(-m_joints[1]->q);
  quaternion<double>::zrot q3( m_joints[2]->q);
  quaternion<double>::zrot q4( m_joints[3]->q);
  quaternion<double>::zrot q5( m_joints[4]->q);
  quaternion<double>::yrot q6(-m_joints[5]->q);
  quaternion<double>::xrot q67( M_PI);
  quaternion<double>::zrot q7( m_joints[6]->q);
  
  vect<double,3> a1_p(link_lengths[0], -joint_offsets[1], 0.0);
  vect<double,3> a2_p(link_lengths[1], 0.0, joint_offsets[2]);
  vect<double,3> a3_p(link_lengths[2], 0.0, joint_offsets[3]);
  vect<double,3> a4_p(link_lengths[3], 0.0, joint_offsets[4]);
  vect<double,3> a5_p(link_lengths[4], 0.0, 0.0);
  vect<double,3> a6_p(link_lengths[5], -joint_offsets[5], 0.0);
  vect<double,3> a7_p(0.0, 0.0, joint_offsets[6]);
  
  quaternion<double> q_accum = q1.getQuaternion();
  vect<double,3> a1 = q_accum * a1_p;
  vect<double,3> e2 = -(q_accum * vect_j);
  q_accum *= q2;
  vect<double,3> a2 = q_accum * a2_p;
  vect<double,3> e3 = q_accum * vect_k;
  q_accum *= q3;
  vect<double,3> a3 = q_accum * a3_p;
  vect<double,3> e4 = q_accum * vect_k;
  q_accum *= q4;
  vect<double,3> a4 = q_accum * a4_p;
  vect<double,3> e5 = q_accum * vect_k;
  q_accum *= q5;
  vect<double,3> a5 = q_accum * a5_p;
  vect<double,3> e6 = -(q_accum * vect_j);
  q_accum *= q6;
  vect<double,3> a6 = q_accum * a6_p;
  q_accum *= q67;
  vect<double,3> e7 = q_accum * vect_k;
  q_accum *= q7;
  vect<double,3> a7 = q_accum * a7_p;
  
  vect<double,3> a67 = a6 + a7;
  vect<double,3> a57 = a5 + a67;
  vect<double,3> a47 = a4 + a57;
  vect<double,3> a37 = a3 + a47;
  vect<double,3> a27 = a2 + a37;
  vect<double,3> a17 = a1 + a27;
  
  Jac.resize(std::make_pair(6,7));
  vect<double,3> v1 = vect_k % a17;
  Jac(0,0) = v1[0];
  Jac(1,0) = v1[1];
  Jac(2,0) = v1[2];
  Jac(3,0) = 0.0;
  Jac(4,0) = 0.0;
  Jac(5,0) = 1.0;
  vect<double,3> v2 = e2 % a27;
  Jac(0,1) = v2[0];
  Jac(1,1) = v2[1];
  Jac(2,1) = v2[2];
  Jac(3,1) = e2[0];
  Jac(4,1) = e2[1];
  Jac(5,1) = e2[2];
  vect<double,3> v3 = e3 % a37;
  Jac(0,2) = v3[0];
  Jac(1,2) = v3[1];
  Jac(2,2) = v3[2];
  Jac(3,2) = e3[0];
  Jac(4,2) = e3[1];
  Jac(5,2) = e3[2];
  vect<double,3> v4 = e4 % a47;
  Jac(0,3) = v4[0];
  Jac(1,3) = v4[1];
  Jac(2,3) = v4[2];
  Jac(3,3) = e4[0];
  Jac(4,3) = e4[1];
  Jac(5,3) = e4[2];
  vect<double,3> v5 = e5 % a57;
  Jac(0,4) = v5[0];
  Jac(1,4) = v5[1];
  Jac(2,4) = v5[2];
  Jac(3,4) = e5[0];
  Jac(4,4) = e5[1];
  Jac(5,4) = e5[2];
  vect<double,3> v6 = e6 % a67;
  Jac(0,5) = v6[0];
  Jac(1,5) = v6[1];
  Jac(2,5) = v6[2];
  Jac(3,5) = e6[0];
  Jac(4,5) = e6[1];
  Jac(5,5) = e6[2];
  Jac(0,6) = 0.0;
  Jac(1,6) = 0.0;
  Jac(2,6) = 0.0;
  Jac(3,6) = e7[0];
  Jac(4,6) = e7[1];
  Jac(5,6) = e7[2];
  
  
  //vect<double,3> q1_dot = m_joints[0]->q_dot * ez;
  //vect<double,3> q2_dot = m_joints[1]->q_dot * ey;
  //vect<double,3> q3_dot = m_joints[2]->q_dot * ez;
  //vect<double,3> q4_dot = m_joints[3]->q_dot * ez;
  //vect<double,3> q5_dot = m_joints[4]->q_dot * ez;
  //vect<double,3> q6_dot = m_joints[5]->q_dot * ey;
  //vect<double,3> q7_dot = m_joints[6]->q_dot * ez;
  vect<double,3> a1_dot = q1 * (m_joints[0]->q_dot * ez % a1_p);
  
  vect<double,3> a2_dot = q2 * (m_joints[1]->q_dot * ey % a2_p);
  a2_p = q2 * a2_p;
  a2_dot = q1 * (m_joints[0]->q_dot * ez % a2_p + a2_dot);
  
  vect<double,3> a3_dot = q3 * (m_joints[2]->q_dot * ez % a3_p);
  a3_p = q3 * a3_p;
  a3_dot = q2 * (m_joints[1]->q_dot * ey % a3_p + a3_dot);
  a3_p = q2 * a3_p;
  a3_dot = q1 * (m_joints[0]->q_dot * ez % a3_p + a3_dot);
  
  vect<double,3> a4_dot = q4 * (m_joints[3]->q_dot * ez % a4_p);
  a4_p = q4 * a4_p;
  a4_dot = q3 * (m_joints[2]->q_dot * ez % a4_p + a4_dot);
  a4_p = q3 * a4_p;
  a4_dot = q2 * (m_joints[1]->q_dot * ey % a4_p + a4_dot);
  a4_p = q2 * a4_p;
  a4_dot = q1 * (m_joints[0]->q_dot * ez % a4_p + a4_dot);
  
  vect<double,3> a5_dot = q5 * (m_joints[4]->q_dot * ez % a5_p);
  a5_p = q5 * a5_p;
  a5_dot = q4 * (m_joints[3]->q_dot * ez % a5_p + a5_dot);
  a5_p = q4 * a5_p;
  a5_dot = q3 * (m_joints[2]->q_dot * ez % a5_p + a5_dot);
  a5_p = q3 * a5_p;
  a5_dot = q2 * (m_joints[1]->q_dot * ey % a5_p + a5_dot);
  a5_p = q2 * a5_p;
  a5_dot = q1 * (m_joints[0]->q_dot * ez % a5_p + a5_dot);
  
  vect<double,3> a6_dot = q6 * (m_joints[5]->q_dot * ey % a6_p);
  a6_p = q6 * a6_p;
  a6_dot = q5 * (m_joints[4]->q_dot * ez % a6_p + a6_dot);
  a6_p = q5 * a6_p;
  a6_dot = q4 * (m_joints[3]->q_dot * ez % a6_p + a6_dot);
  a6_p = q4 * a6_p;
  a6_dot = q3 * (m_joints[2]->q_dot * ez % a6_p + a6_dot);
  a6_p = q3 * a6_p;
  a6_dot = q2 * (m_joints[1]->q_dot * ey % a6_p + a6_dot);
  a6_p = q2 * a6_p;
  a6_dot = q1 * (m_joints[0]->q_dot * ez % a6_p + a6_dot);
  
  vect<double,3> a7_dot = q67 * q7 * (m_joints[6]->q_dot * ez % a7_p);
  a7_p = q67 * q7 * a7_p;
  a7_dot = q6 * (m_joints[5]->q_dot * ey % a7_p + a7_dot);
  a7_p = q6 * a7_p;
  a7_dot = q5 * (m_joints[4]->q_dot * ez % a7_p + a7_dot);
  a7_p = q5 * a7_p;
  a7_dot = q4 * (m_joints[3]->q_dot * ez % a7_p + a7_dot);
  a7_p = q4 * a7_p;
  a7_dot = q3 * (m_joints[2]->q_dot * ez % a7_p + a7_dot);
  a7_p = q3 * a7_p;
  a7_dot = q2 * (m_joints[1]->q_dot * ey % a7_p + a7_dot);
  a7_p = q2 * a7_p;
  a7_dot = q1 * (m_joints[0]->q_dot * ez % a7_p + a7_dot);
  
  
  vect<double,3> ey_p;
  vect<double,3> ez_p;
  
  vect<double,3> e2_dot = q1 * (m_joints[0]->q_dot * ez % ey);
  
  vect<double,3> e3_dot = q2 * (m_joints[1]->q_dot * ey % ez);
  ez_p = q2 * ez;
  e3_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e3_dot);
  
  ez_p = q3 * ez;
  vect<double,3> e4_dot = q2 * (m_joints[1]->q_dot * ey % ez_p);
  ez_p = q2 * ez_p;
  e4_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e4_dot);
  
  ez_p = q4 * ez;
  vect<double,3> e5_dot = q3 * (m_joints[2]->q_dot * ez % ez_p);
  ez_p = q3 * ez_p;
  e5_dot = q2 * (m_joints[1]->q_dot * ey % ez_p + e5_dot);
  ez_p = q2 * ez_p;
  e5_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e5_dot);
  
  vect<double,3> e6_dot = q5 * (m_joints[4]->q_dot * ez % ey);
  ey_p = q5 * ey;
  e6_dot = q4 * (m_joints[3]->q_dot * ez % ey_p + e6_dot);
  ey_p = q4 * ey_p;
  e6_dot = q3 * (m_joints[2]->q_dot * ez % ey_p + e6_dot);
  ey_p = q3 * ey_p;
  e6_dot = q2 * (m_joints[1]->q_dot * ey % ey_p + e6_dot);
  ey_p = q2 * ey_p;
  e6_dot = q1 * (m_joints[0]->q_dot * ez % ey_p + e6_dot);
  
  ez_p = q67 * ez;
  vect<double,3> e7_dot = q6 * (m_joints[5]->q_dot * ey % ez_p);
  ez_p = q6 * ez_p;
  e7_dot = q5 * (m_joints[4]->q_dot * ez % ez_p + e7_dot);
  ez_p = q5 * ez_p;
  e7_dot = q4 * (m_joints[3]->q_dot * ez % ez_p + e7_dot);
  ez_p = q4 * ez_p;
  e7_dot = q3 * (m_joints[2]->q_dot * ez % ez_p + e7_dot);
  ez_p = q3 * ez_p;
  e7_dot = q2 * (m_joints[1]->q_dot * ey % ez_p + e7_dot);
  ez_p = q2 * ez_p;
  e7_dot = q1 * (m_joints[0]->q_dot * ez % ez_p + e7_dot);
  
  
  vect<double,3> a67_dot = a6_dot + a7_dot;
  vect<double,3> a57_dot = a5_dot + a67_dot;
  vect<double,3> a47_dot = a4_dot + a57_dot;
  vect<double,3> a37_dot = a3_dot + a47_dot;
  vect<double,3> a27_dot = a2_dot + a37_dot;
  vect<double,3> a17_dot = a1_dot + a27_dot;
  
  
  JacDot.resize(std::make_pair(6,7));
  vect<double,3> v1_dot = vect_k % a17_dot;
  JacDot(0,0) = v1_dot[0];
  JacDot(1,0) = v1_dot[1];
  JacDot(2,0) = v1_dot[2];
  JacDot(3,0) = 0.0;
  JacDot(4,0) = 0.0;
  JacDot(5,0) = 0.0;
  vect<double,3> v2_dot = e2_dot % a27 + e2 % a27_dot;
  JacDot(0,1) = v2_dot[0];
  JacDot(1,1) = v2_dot[1];
  JacDot(2,1) = v2_dot[2];
  JacDot(3,1) = e2_dot[0];
  JacDot(4,1) = e2_dot[1];
  JacDot(5,1) = e2_dot[2];
  vect<double,3> v3_dot = e3_dot % a37 + e3 % a37_dot;
  JacDot(0,2) = v3_dot[0];
  JacDot(1,2) = v3_dot[1];
  JacDot(2,2) = v3_dot[2];
  JacDot(3,2) = e3_dot[0];
  JacDot(4,2) = e3_dot[1];
  JacDot(5,2) = e3_dot[2];
  vect<double,3> v4_dot = e4_dot % a47 + e4 % a47_dot;
  JacDot(0,3) = v4_dot[0];
  JacDot(1,3) = v4_dot[1];
  JacDot(2,3) = v4_dot[2];
  JacDot(3,3) = e4_dot[0];
  JacDot(4,3) = e4_dot[1];
  JacDot(5,3) = e4_dot[2];
  vect<double,3> v5_dot = e5_dot % a57 + e5 % a57_dot;
  JacDot(0,4) = v5_dot[0];
  JacDot(1,4) = v5_dot[1];
  JacDot(2,4) = v5_dot[2];
  JacDot(3,4) = e5_dot[0];
  JacDot(4,4) = e5_dot[1];
  JacDot(5,4) = e5_dot[2];
  vect<double,3> v6_dot = e6_dot % a67 + e6 % a67_dot;
  JacDot(0,5) = v6_dot[0];
  JacDot(1,5) = v6_dot[1];
  JacDot(2,5) = v6_dot[2];
  JacDot(3,5) = e6_dot[0];
  JacDot(4,5) = e6_dot[1];
  JacDot(5,5) = e6_dot[2];
  JacDot(0,6) = 0.0;
  JacDot(1,6) = 0.0;
  JacDot(2,6) = 0.0;
  JacDot(3,6) = e7_dot[0];
  JacDot(4,6) = e7_dot[1];
  JacDot(5,6) = e7_dot[2];
};

vect_n<double> manip_SSRMS_kinematics::getJointPositions() const {
  return vect_n<double>(m_joints[0]->q, 
                        m_joints[1]->q, 
                        m_joints[2]->q, 
                        m_joints[3]->q, 
                        m_joints[4]->q, 
                        m_joints[5]->q, 
                        m_joints[6]->q);
};

void manip_SSRMS_kinematics::setJointPositions(const vect_n<double>& aJointPositions) {
  m_joints[0]->q = aJointPositions[0];
  m_joints[1]->q = aJointPositions[1];
  m_joints[2]->q = aJointPositions[2];
  m_joints[3]->q = aJointPositions[3];
  m_joints[4]->q = aJointPositions[4];
  m_joints[5]->q = aJointPositions[5];
  m_joints[6]->q = aJointPositions[6];
};

vect_n<double> manip_SSRMS_kinematics::getJointVelocities() const {
  return vect_n<double>(m_joints[0]->q_dot, 
                        m_joints[1]->q_dot, 
                        m_joints[2]->q_dot, 
                        m_joints[3]->q_dot, 
                        m_joints[4]->q_dot, 
                        m_joints[5]->q_dot, 
                        m_joints[6]->q_dot);
};

void manip_SSRMS_kinematics::setJointVelocities(const vect_n<double>& aJointVelocities) {
  m_joints[0]->q_dot = aJointVelocities[0];
  m_joints[1]->q_dot = aJointVelocities[1];
  m_joints[2]->q_dot = aJointVelocities[2];
  m_joints[3]->q_dot = aJointVelocities[3];
  m_joints[4]->q_dot = aJointVelocities[4];
  m_joints[5]->q_dot = aJointVelocities[5];
  m_joints[6]->q_dot = aJointVelocities[6];
};

vect_n<double> manip_SSRMS_kinematics::getJointAccelerations() const {
  return vect_n<double>(m_joints[0]->q_ddot, 
                        m_joints[1]->q_ddot, 
                        m_joints[2]->q_ddot, 
                        m_joints[3]->q_ddot, 
                        m_joints[4]->q_ddot, 
                        m_joints[5]->q_ddot, 
                        m_joints[6]->q_ddot);
};

void manip_SSRMS_kinematics::setJointAccelerations(const vect_n<double>& aJointAccelerations) {
  m_joints[0]->q_ddot = aJointAccelerations[0];
  m_joints[1]->q_ddot = aJointAccelerations[1];
  m_joints[2]->q_ddot = aJointAccelerations[2];
  m_joints[3]->q_ddot = aJointAccelerations[3];
  m_joints[4]->q_ddot = aJointAccelerations[4];
  m_joints[5]->q_ddot = aJointAccelerations[5];
  m_joints[6]->q_ddot = aJointAccelerations[6];
};

vect_n<double> manip_SSRMS_kinematics::getDependentPositions() const {
  return vect_n<double>(
    m_EE->mFrame->Position[0], m_EE->mFrame->Position[1], m_EE->mFrame->Position[2],
    m_EE->mFrame->Quat[0], m_EE->mFrame->Quat[1], m_EE->mFrame->Quat[2], m_EE->mFrame->Quat[3]);
};

void manip_SSRMS_kinematics::setDependentPositions(const vect_n<double>& aDepPositions) {
  m_EE->mFrame->Position[0] = aDepPositions[0];
  m_EE->mFrame->Position[1] = aDepPositions[1];
  m_EE->mFrame->Position[2] = aDepPositions[2];
  m_EE->mFrame->Quat = quaternion<double>(vect<double,4>(aDepPositions[3], aDepPositions[4], aDepPositions[5], aDepPositions[6]));
};

vect_n<double> manip_SSRMS_kinematics::getDependentVelocities() const {
  return vect_n<double>(
    m_EE->mFrame->Velocity[0], m_EE->mFrame->Velocity[1], m_EE->mFrame->Velocity[2],
    m_EE->mFrame->AngVelocity[0], m_EE->mFrame->AngVelocity[1], m_EE->mFrame->AngVelocity[2]);
};

void manip_SSRMS_kinematics::setDependentVelocities(const vect_n<double>& aDepVelocities) {
  m_EE->mFrame->Velocity[0] = aDepVelocities[0];
  m_EE->mFrame->Velocity[1] = aDepVelocities[1];
  m_EE->mFrame->Velocity[2] = aDepVelocities[2];
  m_EE->mFrame->AngVelocity[0] = aDepVelocities[3];
  m_EE->mFrame->AngVelocity[1] = aDepVelocities[4];
  m_EE->mFrame->AngVelocity[2] = aDepVelocities[5];
};

vect_n<double> manip_SSRMS_kinematics::getDependentAccelerations() const {
  return vect_n<double>(
    m_EE->mFrame->Acceleration[0], m_EE->mFrame->Acceleration[1], m_EE->mFrame->Acceleration[2],
    m_EE->mFrame->AngAcceleration[0], m_EE->mFrame->AngAcceleration[1], m_EE->mFrame->AngAcceleration[2]);
};

void manip_SSRMS_kinematics::setDependentAccelerations(const vect_n<double>& aDepAccelerations) {
  m_EE->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_EE->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_EE->mFrame->Acceleration[2] = aDepAccelerations[2];
  m_EE->mFrame->AngAcceleration[0] = aDepAccelerations[3];
  m_EE->mFrame->AngAcceleration[1] = aDepAccelerations[4];
  m_EE->mFrame->AngAcceleration[2] = aDepAccelerations[5];
};


void RK_CALL manip_SSRMS_kinematics::save(serialization::oarchive& A, unsigned int) const {
  inverse_kinematics_model::save(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(m_base_frame)
    & RK_SERIAL_SAVE_WITH_NAME(m_joints)
    & RK_SERIAL_SAVE_WITH_NAME(m_EE)
    & RK_SERIAL_SAVE_WITH_NAME(link_lengths)
    & RK_SERIAL_SAVE_WITH_NAME(joint_offsets)
    & RK_SERIAL_SAVE_WITH_NAME(joint_lower_bounds)
    & RK_SERIAL_SAVE_WITH_NAME(joint_upper_bounds)
    & RK_SERIAL_SAVE_WITH_NAME(m_chain);
};

void RK_CALL manip_SSRMS_kinematics::load(serialization::iarchive& A, unsigned int) {
  inverse_kinematics_model::load(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(m_base_frame)
    & RK_SERIAL_LOAD_WITH_NAME(m_joints)
    & RK_SERIAL_LOAD_WITH_NAME(m_EE)
    & RK_SERIAL_LOAD_WITH_NAME(link_lengths)
    & RK_SERIAL_LOAD_WITH_NAME(joint_offsets)
    & RK_SERIAL_LOAD_WITH_NAME(joint_lower_bounds)
    & RK_SERIAL_LOAD_WITH_NAME(joint_upper_bounds)
    & RK_SERIAL_LOAD_WITH_NAME(m_chain);
};


};

};








