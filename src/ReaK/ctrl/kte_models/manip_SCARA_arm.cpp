
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

#include "manip_SCARA_arm.hpp"

#include "mbd_kte/revolute_joint.hpp"
#include "mbd_kte/rigid_link.hpp"

#include "mbd_kte/prismatic_joint.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_qr_decomp.hpp"
#include "optimization/optim_exceptions.hpp"
#include <cmath>

namespace ReaK {

namespace kte {


//     shared_ptr< frame_3D<double> > m_base_frame;
//     std::vector< shared_ptr< gen_coord<double> > > m_joints;
//     shared_ptr< joint_dependent_frame_3D > m_EE;
//     double m_link1_length;
//     double m_link2_length;
//     double m_link3_height;

manip_SCARA_kinematics::manip_SCARA_kinematics(const std::string& aName,
                                               const shared_ptr< frame_3D<double> >& aBaseFrame,
                                               double aLink1Length,
                                               double aLink2Length,
                                               double aLink3Height) :
                                               inverse_kinematics_model(aName),
                                               m_base_frame(aBaseFrame),
                                               m_link1_length(aLink1Length), 
                                               m_link2_length(aLink2Length),
                                               m_link3_height(aLink3Height) {
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  m_joints.push_back(shared_ptr< gen_coord<double> >(new gen_coord<double>(), scoped_deleter()));
  
  
  //declare all the intermediate frames.
  shared_ptr< frame_3D<double> > joint_1_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_2_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_2_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_3_base(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > joint_3_end(new frame_3D<double>(), scoped_deleter());
  shared_ptr< frame_3D<double> > arm_EE(new frame_3D<double>(), scoped_deleter());
  
  //declare all the joint jacobians.
  shared_ptr< jacobian_gen_3D<double> > joint_1_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_2_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  shared_ptr< jacobian_gen_3D<double> > joint_3_jacobian(new jacobian_gen_3D<double>(), scoped_deleter());
  
  
  //create revolute joint
  shared_ptr< kte::revolute_joint_3D > joint_1(new kte::revolute_joint_3D(
      "manip_SCARA_joint_1",
      m_joints[0],
      vect<double,3>(0.0,0.0,1.0),
      m_base_frame,
      joint_1_end,
      joint_1_jacobian),
    scoped_deleter());
  
  //create link from F to CM (note that this is very approximate!!!)
  shared_ptr< kte::rigid_link_3D > link_1(new kte::rigid_link_3D(
      "manip_SCARA_link_1",
      joint_1_end,
      joint_2_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(m_link1_length, 0.0, 0.0),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< kte::revolute_joint_3D > joint_2(new kte::revolute_joint_3D(
      "manip_SCARA_joint_2",
      m_joints[1],
      vect<double,3>(0.0,0.0,1.0),
      joint_2_base,
      joint_2_end,
      joint_2_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< kte::rigid_link_3D > link_2(new kte::rigid_link_3D(
      "manip_SCARA_link_2",
      joint_2_end,
      joint_3_base,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(m_link2_length, 0.0, 0.0),
        quaternion<double>())),
    scoped_deleter());
  
  //create revolute joint
  shared_ptr< kte::prismatic_joint_3D > joint_3(new kte::prismatic_joint_3D(
      "manip_SCARA_joint_3",
      m_joints[2],
      vect<double,3>(0.0,0.0,1.0),
      joint_3_base,
      joint_3_end,
      joint_3_jacobian),
    scoped_deleter());
  
  //create link 
  shared_ptr< kte::rigid_link_3D > link_3(new kte::rigid_link_3D(
      "manip_SCARA_link_3",
      joint_3_end,
      arm_EE,
      pose_3D<double>(
        weak_ptr<pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, m_link3_height),
        quaternion<double>())),
    scoped_deleter());
  
  //create inertia
  m_EE = shared_ptr< joint_dependent_frame_3D >(new joint_dependent_frame_3D(
      arm_EE),
    scoped_deleter());
  
  m_EE->add_joint(m_joints[0], joint_1_jacobian);
  m_EE->add_joint(m_joints[1], joint_2_jacobian);
  m_EE->add_joint(m_joints[2], joint_3_jacobian);
  
  m_chain = shared_ptr< kte_map_chain >(new kte_map_chain("manip_SCARA_kin_model"), scoped_deleter());
  
  *m_chain << joint_1
           << link_1
           << joint_2
           << link_2
           << joint_3
           << link_3;
  
};



void manip_SCARA_kinematics::doDirectMotion() {
  m_chain->doMotion();
};

void manip_SCARA_kinematics::doInverseMotion() {
  using std::acos;  using std::sin; using std::sqrt;
  using std::atan2; using std::fabs;
  
  frame_3D<double> EE_fr = m_EE->mFrame->getFrameRelativeTo(m_base_frame);
  vect<double,3> w_pos = EE_fr.Position;
  
  double w_dist = sqrt(w_pos[0] * w_pos[0] + w_pos[1] * w_pos[1]);
  if(w_dist > m_link1_length + m_link2_length)
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! Desired end-effector pose of the SCARA manipulator is outside the reachable workspace!");
  if(w_dist < 1.000001 * fabs(m_link1_length - m_link2_length))
    throw optim::infeasible_problem("Inverse kinematics problem is infeasible! Desired end-effector pose of the SCARA manipulator is too close to the center (inner work-space limit)!");
  
  double a2 = acos(-(m_link1_length * m_link1_length + m_link2_length * m_link2_length - w_dist * w_dist) / (2.0 * m_link1_length * m_link2_length));
  if(fabs(a2 - m_joints[0]->q) > fabs(a2 + m_joints[0]->q))
    a2 = -a2;
  
  double c1_p = (m_link1_length * m_link1_length + w_dist * w_dist - m_link2_length * m_link2_length) / (2.0 * m_link1_length * w_dist);  
  double s1_p = m_link2_length * sin(a2) / w_dist;
  double c1_c = w_pos[0] / w_dist;
  double s1_c = w_pos[1] / w_dist;
  double c1 = c1_c * c1_p - s1_c * s1_p;
  double s1 = s1_c * c1_p + c1_c * s1_p;
  double a1 = atan2(s1, c1);
  
  rot_mat_2D<double> R1 = rot_mat_2D<double>(a1);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(a1 + a2);
  
  vect<double,2> cl1 = double(1.0) % ( R1 * vect<double,2>(m_link1_length, 0.0) );
  vect<double,2> cl2 = double(1.0) % ( R12 * vect<double,2>(m_link2_length, 0.0) );
  
  mat<double,mat_structure::rectangular> A(2,2);
  A(0,0) = cl1[0]; A(1,0) = cl1[1];
  A(0,1) = cl2[0]; A(1,1) = cl2[1];
  mat<double,mat_structure::rectangular> b(2,1);
  b(0,0) = EE_fr.Velocity[0]; b(1,0) = EE_fr.Velocity[1];
  mat<double,mat_structure::rectangular> jt_vel(2,1);
  linlsq_QR(A, jt_vel, b); // solve for the joint velocities.
  
  m_joints[0]->q      = a1;
  m_joints[0]->q_dot  = jt_vel(0,0);
  m_joints[0]->q_ddot = 0.0;
  m_joints[1]->q      = a2;
  m_joints[1]->q_dot  = jt_vel(1,0);
  m_joints[1]->q_ddot = 0.0;
  m_joints[2]->q      = w_pos[2] - m_link3_height;
  m_joints[2]->q_dot  = EE_fr.Velocity[2];
  m_joints[2]->q_ddot = 0.0;
  
  m_chain->doMotion();
};

void manip_SCARA_kinematics::getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const {
  rot_mat_2D<double> R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  
  vect<double,2> cl1 = double(1.0) % ( R1 * vect<double,2>(m_link1_length, 0.0) );
  vect<double,2> cl2 = double(1.0) % ( R12 * vect<double,2>(m_link2_length, 0.0) );
  
  Jac.resize(std::make_pair(6,3));
  Jac(0,0) = cl1[0]; Jac(1,0) = cl1[1]; Jac(2,0) = 0.0; Jac(3,0) = 0.0; Jac(4,0) = 0.0; Jac(5,0) = 1.0;
  Jac(0,1) = cl2[0]; Jac(1,1) = cl2[1]; Jac(2,1) = 0.0; Jac(3,1) = 0.0; Jac(4,1) = 0.0; Jac(5,1) = 1.0;
  Jac(0,2) = 0.0;    Jac(1,2) = 0.0;    Jac(2,2) = 1.0; Jac(3,2) = 0.0; Jac(4,2) = 0.0; Jac(5,2) = 0.0;
};

void manip_SCARA_kinematics::getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const {
  rot_mat_2D<double> R1 = rot_mat_2D<double>(m_joints[0]->q);
  rot_mat_2D<double> R12 = rot_mat_2D<double>(m_joints[0]->q + m_joints[1]->q);
  
  vect<double,2> l1 = R1 * vect<double,2>(m_link1_length, 0.0);
  vect<double,2> l2 = R12 * vect<double,2>(m_link2_length, 0.0);
  
  Jac.resize(std::make_pair(6,3));
  Jac(0,0) = -l1[1]; Jac(1,0) = l1[0]; Jac(2,0) = 0.0; Jac(3,0) = 0.0; Jac(4,0) = 0.0; Jac(5,0) = 1.0;
  Jac(0,1) = -l2[1]; Jac(1,1) = l2[0]; Jac(2,1) = 0.0; Jac(3,1) = 0.0; Jac(4,1) = 0.0; Jac(5,1) = 1.0;
  Jac(0,2) = 0.0;    Jac(1,2) = 0.0;   Jac(2,2) = 1.0; Jac(3,2) = 0.0; Jac(4,2) = 0.0; Jac(5,2) = 0.0;
  
  JacDot.resize(std::make_pair(6,3));
  JacDot(0,0) = -m_joints[0]->q_dot * l1[0]; JacDot(1,0) = -m_joints[0]->q_dot * l1[1]; JacDot(2,0) = 0.0; JacDot(3,0) = 0.0; JacDot(4,0) = 0.0; JacDot(5,0) = 0.0;
  JacDot(0,1) = -m_joints[1]->q_dot * l2[0]; JacDot(1,1) = -m_joints[1]->q_dot * l2[1]; JacDot(2,1) = 0.0; JacDot(3,1) = 0.0; JacDot(4,1) = 0.0; JacDot(5,1) = 0.0;
  JacDot(0,2) = 0.0; JacDot(1,2) = 0.0; JacDot(2,2) = 0.0; JacDot(3,2) = 0.0; JacDot(4,2) = 0.0; JacDot(5,2) = 0.0;
};

vect_n<double> manip_SCARA_kinematics::getJointPositions() const {
  return vect_n<double>(m_joints[0]->q, m_joints[1]->q, m_joints[2]->q);
};

void manip_SCARA_kinematics::setJointPositions(const vect_n<double>& aJointPositions) {
  m_joints[0]->q = aJointPositions[0];
  m_joints[1]->q = aJointPositions[1];
  m_joints[2]->q = aJointPositions[2];
};

vect_n<double> manip_SCARA_kinematics::getJointVelocities() const {
  return vect_n<double>(m_joints[0]->q_dot, m_joints[1]->q_dot, m_joints[2]->q_dot);
};

void manip_SCARA_kinematics::setJointVelocities(const vect_n<double>& aJointVelocities) {
  m_joints[0]->q_dot = aJointVelocities[0];
  m_joints[1]->q_dot = aJointVelocities[1];
  m_joints[2]->q_dot = aJointVelocities[2];
};

vect_n<double> manip_SCARA_kinematics::getJointAccelerations() const {
  return vect_n<double>(m_joints[0]->q_ddot, m_joints[1]->q_ddot, m_joints[2]->q_ddot);
};

void manip_SCARA_kinematics::setJointAccelerations(const vect_n<double>& aJointAccelerations) {
  m_joints[0]->q_ddot = aJointAccelerations[0];
  m_joints[1]->q_ddot = aJointAccelerations[1];
  m_joints[2]->q_ddot = aJointAccelerations[2];
};

vect_n<double> manip_SCARA_kinematics::getDependentPositions() const {
  return vect_n<double>(
    m_EE->mFrame->Position[0], m_EE->mFrame->Position[1], m_EE->mFrame->Position[2],
    m_EE->mFrame->Quat[0], m_EE->mFrame->Quat[1], m_EE->mFrame->Quat[2], m_EE->mFrame->Quat[3]);
};

void manip_SCARA_kinematics::setDependentPositions(const vect_n<double>& aDepPositions) {
  m_EE->mFrame->Position[0] = aDepPositions[0];
  m_EE->mFrame->Position[1] = aDepPositions[1];
  m_EE->mFrame->Position[2] = aDepPositions[2];
  m_EE->mFrame->Quat = quaternion<double>(vect<double,4>(aDepPositions[3], aDepPositions[4], aDepPositions[5], aDepPositions[6]));
};

vect_n<double> manip_SCARA_kinematics::getDependentVelocities() const {
  return vect_n<double>(
    m_EE->mFrame->Velocity[0], m_EE->mFrame->Velocity[1], m_EE->mFrame->Velocity[2],
    m_EE->mFrame->AngVelocity[0], m_EE->mFrame->AngVelocity[1], m_EE->mFrame->AngVelocity[2]);
};

void manip_SCARA_kinematics::setDependentVelocities(const vect_n<double>& aDepVelocities) {
  m_EE->mFrame->Velocity[0] = aDepVelocities[0];
  m_EE->mFrame->Velocity[1] = aDepVelocities[1];
  m_EE->mFrame->Velocity[2] = aDepVelocities[2];
  m_EE->mFrame->AngVelocity[0] = aDepVelocities[3];
  m_EE->mFrame->AngVelocity[1] = aDepVelocities[4];
  m_EE->mFrame->AngVelocity[2] = aDepVelocities[5];
};

vect_n<double> manip_SCARA_kinematics::getDependentAccelerations() const {
  return vect_n<double>(
    m_EE->mFrame->Acceleration[0], m_EE->mFrame->Acceleration[1], m_EE->mFrame->Acceleration[2],
    m_EE->mFrame->AngAcceleration[0], m_EE->mFrame->AngAcceleration[1], m_EE->mFrame->AngAcceleration[2]);
};

void manip_SCARA_kinematics::setDependentAccelerations(const vect_n<double>& aDepAccelerations) {
  m_EE->mFrame->Acceleration[0] = aDepAccelerations[0];
  m_EE->mFrame->Acceleration[1] = aDepAccelerations[1];
  m_EE->mFrame->Acceleration[2] = aDepAccelerations[2];
  m_EE->mFrame->AngAcceleration[0] = aDepAccelerations[3];
  m_EE->mFrame->AngAcceleration[1] = aDepAccelerations[4];
  m_EE->mFrame->AngAcceleration[2] = aDepAccelerations[5];
};


void RK_CALL manip_SCARA_kinematics::save(serialization::oarchive& A, unsigned int) const {
  inverse_kinematics_model::save(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(m_base_frame)
    & RK_SERIAL_SAVE_WITH_NAME(m_joints)
    & RK_SERIAL_SAVE_WITH_NAME(m_EE)
    & RK_SERIAL_SAVE_WITH_NAME(m_link1_length)
    & RK_SERIAL_SAVE_WITH_NAME(m_link2_length)
    & RK_SERIAL_SAVE_WITH_NAME(m_link3_height)
    & RK_SERIAL_SAVE_WITH_NAME(m_chain);
};

void RK_CALL manip_SCARA_kinematics::load(serialization::iarchive& A, unsigned int) {
  inverse_kinematics_model::load(A,inverse_kinematics_model::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(m_base_frame)
    & RK_SERIAL_LOAD_WITH_NAME(m_joints)
    & RK_SERIAL_LOAD_WITH_NAME(m_EE)
    & RK_SERIAL_LOAD_WITH_NAME(m_link1_length)
    & RK_SERIAL_LOAD_WITH_NAME(m_link2_length)
    & RK_SERIAL_LOAD_WITH_NAME(m_link3_height)
    & RK_SERIAL_LOAD_WITH_NAME(m_chain);
};


};

};








