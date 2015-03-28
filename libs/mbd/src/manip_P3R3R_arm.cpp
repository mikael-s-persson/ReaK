
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

#include <ReaK/math/lin_alg/mat_alg.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/optimization/optim_exceptions.hpp>

#include <ReaK/mbd/models/manip_P3R3R_arm.hpp>
#include <ReaK/mbd/kte/prismatic_joint.hpp>

#include <cmath>
#include <algorithm>

namespace ReaK {

namespace kte {

//     shared_ptr< frame_3D<double> > m_base_frame;
//     shared_ptr< gen_coord<double> > m_track_coord;
//     shared_ptr< joint_dependent_frame_3D > m_output_frame;
//     double m_track_lower_bound;
//     double m_track_upper_bound;
//
//     manip_3R3R_kinematics m_arm_model;
//     shared_ptr< kte_map_chain > m_chain;

manip_P3R3R_kinematics::manip_P3R3R_kinematics( const std::string& aName,
                                                const shared_ptr< frame_3D< double > >& aBaseFrame,
                                                double aBaseToShoulder, double aShoulderToElbow, double aElbowToJoint4,
                                                double aJoint4ToWrist, double aWristToFlange,
                                                const vect_n< double >& aJointLowerBounds,
                                                const vect_n< double >& aJointUpperBounds )
    : inverse_kinematics_model( aName ), m_base_frame( aBaseFrame ),
      m_track_coord( new gen_coord< double >(), scoped_deleter() ),
      m_output_frame( new frame_3D< double >(), scoped_deleter() ), m_track_lower_bound( aJointLowerBounds[0] ),
      m_track_upper_bound( aJointUpperBounds[0] ),
      m_arm_model( aName + "_arm", m_output_frame, aBaseToShoulder, aShoulderToElbow, aElbowToJoint4, aJoint4ToWrist,
                   aWristToFlange, vect_n< double >( aJointLowerBounds.begin() + 1, aJointLowerBounds.end() ),
                   vect_n< double >( aJointUpperBounds.begin() + 1, aJointUpperBounds.end() ) ),
      m_chain() {

  if( !m_base_frame )
    m_base_frame = shared_ptr< frame_3D< double > >( new frame_3D< double >(), scoped_deleter() );


  // declare all the joint jacobians.
  shared_ptr< jacobian_gen_3D< double > > track_jacobian( new jacobian_gen_3D< double >(), scoped_deleter() );

  // create prismatic joint
  m_track_joint = shared_ptr< prismatic_joint_3D >(
    new prismatic_joint_3D( "manip_P3R3R_track", m_track_coord, vect< double, 3 >( 1.0, 0.0, 0.0 ), m_base_frame,
                            m_output_frame, track_jacobian ),
    scoped_deleter() );

  m_chain = shared_ptr< kte_map_chain >( new kte_map_chain( "manip_P3R3R_kin_model" ), scoped_deleter() );

  ( *m_chain ) << m_track_joint << m_arm_model.getKTEChain();

  m_arm_model.getDependentFrame3D( 0 )->add_joint( m_track_coord, track_jacobian );
};


void manip_P3R3R_kinematics::doDirectMotion() { m_chain->doMotion(); };


void manip_P3R3R_kinematics::doInverseMotion() {
  using std::sin;
  using std::cos;
  using std::fabs;
  using std::atan2;
  using std::sqrt;
  using std::pow;

  frame_3D< double > EE_fr = m_arm_model.getDependentFrame3D( 0 )->mFrame->getFrameRelativeTo( m_base_frame );

  const double extend_epsilon = .02;

  vect< double, 3 > EE_z_axis = EE_fr.Quat * vect_k;
  vect< double, 3 > wrist_pos = EE_fr.Position - m_arm_model.getWristToFlange() * EE_z_axis
                                - m_arm_model.getBaseToShoulder() * vect_k;
  double elbow_to_wrist_dist = m_arm_model.getElbowToJoint4() + m_arm_model.getJoint4ToWrist();

  double perp_wrist_dist_sqr = wrist_pos[2] * wrist_pos[2] + wrist_pos[1] * wrist_pos[1];

  /*
   * find the maximum wrist to base distance, x_max and verifies if the required
   * position of the end-effector is within limits
   */

  double c2_max = 1.0;
  if( m_arm_model.joint_upper_bounds[1] < 0.0 )
    c2_max = cos( m_arm_model.joint_upper_bounds[1] );
  if( m_arm_model.joint_lower_bounds[1] > 0.0 )
    c2_max = cos( m_arm_model.joint_lower_bounds[1] );
  double c2_min = -1.0;
  if( ( m_arm_model.joint_upper_bounds[1] < M_PI ) && ( m_arm_model.joint_lower_bounds[1] > -M_PI ) ) {
    if( M_PI - m_arm_model.joint_upper_bounds[1] > M_PI + m_arm_model.joint_lower_bounds[1] )
      c2_min = cos( m_arm_model.joint_lower_bounds[1] );
    else
      c2_min = cos( m_arm_model.joint_upper_bounds[1] );
  };

  double c3_max = 1.0;
  double c3_max_s3 = 0.0;
  if( m_arm_model.joint_upper_bounds[2] < 0.0 ) {
    c3_max = cos( m_arm_model.joint_upper_bounds[2] );
    c3_max_s3 = sin( m_arm_model.joint_upper_bounds[2] );
  };
  if( m_arm_model.joint_lower_bounds[2] > 0.0 ) {
    c3_max = cos( m_arm_model.joint_lower_bounds[2] );
    c3_max_s3 = sin( m_arm_model.joint_lower_bounds[2] );
  };
  double c3_min = -1.0;
  double c3_min_s3 = 0.0;
  if( ( m_arm_model.joint_upper_bounds[2] < M_PI ) && ( m_arm_model.joint_lower_bounds[2] > -M_PI ) ) {
    if( M_PI - m_arm_model.joint_upper_bounds[2] > M_PI + m_arm_model.joint_lower_bounds[2] ) {
      c3_min = cos( m_arm_model.joint_lower_bounds[2] );
      c3_min_s3 = sin( m_arm_model.joint_lower_bounds[2] );
    } else {
      c3_min = cos( m_arm_model.joint_upper_bounds[2] );
      c3_min_s3 = sin( m_arm_model.joint_upper_bounds[2] );
    };
  };

  double j23_lower_bound = m_arm_model.joint_lower_bounds[1] + m_arm_model.joint_lower_bounds[2];
  double j23_upper_bound = m_arm_model.joint_upper_bounds[1] + m_arm_model.joint_upper_bounds[2];
  double c23_max = 1.0;
  if( j23_upper_bound < 0.0 )
    c23_max = cos( j23_upper_bound );
  if( j23_lower_bound > 0.0 )
    c23_max = cos( j23_lower_bound );
  double c23_min = -1.0;
  if( ( j23_upper_bound < M_PI ) && ( j23_lower_bound > -M_PI ) ) {
    if( M_PI - j23_upper_bound > M_PI + j23_lower_bound )
      c3_min = cos( j23_lower_bound );
    else
      c3_min = cos( j23_upper_bound );
  };

  if( wrist_pos[2] >= 0.0 ) {
    double elbow_to_desired_height_min = wrist_pos[2] - m_arm_model.getShoulderToElbow() * c2_max;
    if( elbow_to_desired_height_min > ( elbow_to_wrist_dist - extend_epsilon ) * c23_max )
      throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! "
                                       "Desired wrist position is too high for the manipulator to reach." );

  } else if( wrist_pos[2] < 0.0 ) {
    double elbow_to_desired_height_min = wrist_pos[2] - m_arm_model.getShoulderToElbow() * c2_min;
    if( elbow_to_desired_height_min < ( elbow_to_wrist_dist - extend_epsilon ) * c23_min )
      throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! "
                                       "Desired wrist position is too low for the manipulator to reach." );
  };

  double shoulder_to_wrist_max_sqr = m_arm_model.getShoulderToElbow() + elbow_to_wrist_dist * c3_max;
  shoulder_to_wrist_max_sqr = shoulder_to_wrist_max_sqr * shoulder_to_wrist_max_sqr
                              + elbow_to_wrist_dist * elbow_to_wrist_dist * c3_max_s3 * c3_max_s3;
  if( perp_wrist_dist_sqr > shoulder_to_wrist_max_sqr ) {
    throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! "
                                     "Desired wrist position is outside the cylindrical workspace envelope "
                                     "(maximally-extended arm)." );
  };

  double x_max = sqrt( shoulder_to_wrist_max_sqr - perp_wrist_dist_sqr );

  double shoulder_to_wrist_min_sqr = m_arm_model.getShoulderToElbow() + elbow_to_wrist_dist * c3_min;
  shoulder_to_wrist_min_sqr = shoulder_to_wrist_min_sqr * shoulder_to_wrist_min_sqr
                              + elbow_to_wrist_dist * elbow_to_wrist_dist * c3_min_s3 * c3_min_s3;
  double x_min = 0.0;
  if( perp_wrist_dist_sqr < shoulder_to_wrist_min_sqr )
    x_min = sqrt( shoulder_to_wrist_min_sqr - perp_wrist_dist_sqr );

  /* At this point, the range of joint values for the track is between +- x_max around wrist_pos[0] and beyond +- x_min
   * around wrist_pos[0]. */

  /* Joint 0 (track) */

  // first, check that the limits of the track permit at least some solution:
  if( wrist_pos[0] - x_max > m_track_upper_bound )
    throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! "
                                     "Desired wrist position is too far ahead of the track (beyond upper track "
                                     "limit)." );
  if( wrist_pos[0] + x_max < m_track_lower_bound )
    throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! "
                                     "Desired wrist position is too far behind the track (beyond lower track limit)." );
  if( ( wrist_pos[0] + x_min > m_track_upper_bound ) && ( wrist_pos[0] - x_min < m_track_lower_bound ) )
    throw optim::infeasible_problem( "Inverse kinematics problem is infeasible! End-effector pose is out-of-reach! The "
                                     "track is too short to allow for a feasible solution." );

  double x_desired = 0.0;
  if( ( EE_z_axis[0] >= 0.0 ) && ( wrist_pos[0] - x_min > m_track_lower_bound ) ) {
    // if EE needs to point forward, then try to keep the manipulator in the smaller x-coord side.
    double x_try_min = wrist_pos[0] - x_min;
    double x_try_max = wrist_pos[0] - x_max;
    if( x_try_max < m_track_lower_bound )
      x_try_max = m_track_lower_bound;
    if( x_try_min > m_track_upper_bound )
      x_try_min = m_track_upper_bound;
    x_desired = ( x_try_min + x_try_max ) * 0.5;
  } else {
    // else, try to keep the manipulator in the bigger x-coord side.
    double x_try_min = wrist_pos[0] + x_min;
    double x_try_max = wrist_pos[0] + x_max;
    if( x_try_max > m_track_upper_bound )
      x_try_max = m_track_upper_bound;
    if( x_try_min < m_track_lower_bound )
      x_try_min = m_track_lower_bound;
    x_desired = ( x_try_min + x_try_max ) * 0.5;
  };

  m_track_coord->q = x_desired;
  m_track_coord->q_dot = EE_fr.Velocity[0];
  m_track_coord->q_ddot = 0.0;

  m_track_joint->doMotion();

  m_arm_model.doInverseMotion();
};

void manip_P3R3R_kinematics::getJacobianMatrix( mat< double, mat_structure::rectangular >& Jac ) const {
  mat< double, mat_structure::rectangular > Jac_manip( 6, 6 );
  m_arm_model.getJacobianMatrix( Jac_manip );

  Jac.resize( std::make_pair( 6, 7 ) );
  Jac( 0, 0 ) = 1.0;
  Jac( 1, 0 ) = 0.0;
  Jac( 2, 0 ) = 0.0;
  Jac( 3, 0 ) = 0.0;
  Jac( 4, 0 ) = 0.0;
  Jac( 5, 0 ) = 0.0;

  for( std::size_t i = 0; i < 5; ++i )
    for( std::size_t j = 0; j < 5; ++j )
      Jac( j, i + 1 ) = Jac_manip( j, i );
};

void manip_P3R3R_kinematics::getJacobianMatrixAndDerivative( mat< double, mat_structure::rectangular >& Jac,
                                                             mat< double, mat_structure::rectangular >& JacDot ) const {
  mat< double, mat_structure::rectangular > Jac_manip( 6, 6 );
  mat< double, mat_structure::rectangular > JacDot_manip( 6, 6 );
  m_arm_model.getJacobianMatrixAndDerivative( Jac_manip, JacDot_manip );

  Jac.resize( std::make_pair( 6, 7 ) );
  Jac( 0, 0 ) = 1.0;
  Jac( 1, 0 ) = 0.0;
  Jac( 2, 0 ) = 0.0;
  Jac( 3, 0 ) = 0.0;
  Jac( 4, 0 ) = 0.0;
  Jac( 5, 0 ) = 0.0;

  for( std::size_t i = 0; i < 6; ++i )
    for( std::size_t j = 0; j < 6; ++j )
      Jac( j, i + 1 ) = Jac_manip( j, i );

  JacDot.resize( std::make_pair( 6, 7 ) );
  JacDot( 0, 0 ) = 0.0;
  JacDot( 1, 0 ) = 0.0;
  JacDot( 2, 0 ) = 0.0;
  JacDot( 3, 0 ) = 0.0;
  JacDot( 4, 0 ) = 0.0;
  JacDot( 5, 0 ) = 0.0;

  for( std::size_t i = 0; i < 6; ++i )
    for( std::size_t j = 0; j < 6; ++j )
      JacDot( j, i + 1 ) = JacDot_manip( j, i );
};


vect_n< double > manip_P3R3R_kinematics::getJointPositionLowerBounds() const {
  vect_n< double > result = m_arm_model.getJointPositionLowerBounds();
  result.resize( 7 );
  std::rotate( result.begin(), result.begin() + 6, result.end() );
  result[0] = m_track_lower_bound;
  return result;
};

void manip_P3R3R_kinematics::setJointPositionLowerBounds( const vect_n< double >& aJointLowerBounds ) {
  m_track_lower_bound = aJointLowerBounds[0];
  m_arm_model.setJointPositionLowerBounds( vect_n< double >( aJointLowerBounds.begin() + 1, aJointLowerBounds.end() ) );
};

vect_n< double > manip_P3R3R_kinematics::getJointPositionUpperBounds() const {
  vect_n< double > result = m_arm_model.getJointPositionUpperBounds();
  result.resize( 7 );
  std::rotate( result.begin(), result.begin() + 6, result.end() );
  result[0] = m_track_upper_bound;
  return result;
};

void manip_P3R3R_kinematics::setJointPositionUpperBounds( const vect_n< double >& aJointUpperBounds ) {
  m_track_upper_bound = aJointUpperBounds[0];
  m_arm_model.setJointPositionUpperBounds( vect_n< double >( aJointUpperBounds.begin() + 1, aJointUpperBounds.end() ) );
};

vect_n< double > manip_P3R3R_kinematics::getJointPositions() const {
  vect_n< double > arm_jtpos = m_arm_model.getJointPositions();
  return vect_n< double >( m_track_coord->q, arm_jtpos[0], arm_jtpos[1], arm_jtpos[2], arm_jtpos[3], arm_jtpos[4],
                           arm_jtpos[5] );
};

void manip_P3R3R_kinematics::setJointPositions( const vect_n< double >& aJointPositions ) {
  m_track_coord->q = aJointPositions[0];
  m_arm_model.setJointPositions( vect_n< double >( aJointPositions.begin() + 1, aJointPositions.end() ) );
};

vect_n< double > manip_P3R3R_kinematics::getJointVelocities() const {
  vect_n< double > arm_jtvel = m_arm_model.getJointVelocities();
  return vect_n< double >( m_track_coord->q_dot, arm_jtvel[0], arm_jtvel[1], arm_jtvel[2], arm_jtvel[3], arm_jtvel[4],
                           arm_jtvel[5] );
};

void manip_P3R3R_kinematics::setJointVelocities( const vect_n< double >& aJointVelocities ) {
  m_track_coord->q_dot = aJointVelocities[0];
  m_arm_model.setJointVelocities( vect_n< double >( aJointVelocities.begin() + 1, aJointVelocities.end() ) );
};

vect_n< double > manip_P3R3R_kinematics::getJointAccelerations() const {
  vect_n< double > arm_jtacc = m_arm_model.getJointAccelerations();
  return vect_n< double >( m_track_coord->q_ddot, arm_jtacc[0], arm_jtacc[1], arm_jtacc[2], arm_jtacc[3], arm_jtacc[4],
                           arm_jtacc[5] );
};

void manip_P3R3R_kinematics::setJointAccelerations( const vect_n< double >& aJointAccelerations ) {
  m_track_coord->q_ddot = aJointAccelerations[0];
  m_arm_model.setJointAccelerations( vect_n< double >( aJointAccelerations.begin() + 1, aJointAccelerations.end() ) );
};


void RK_CALL manip_P3R3R_kinematics::save( serialization::oarchive& A, unsigned int ) const {
  inverse_kinematics_model::save( A, inverse_kinematics_model::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_SAVE_WITH_NAME( m_base_frame ) & RK_SERIAL_SAVE_WITH_NAME( m_track_coord )
    & RK_SERIAL_SAVE_WITH_NAME( m_output_frame ) & RK_SERIAL_SAVE_WITH_NAME( m_track_joint )
    & RK_SERIAL_SAVE_WITH_NAME( m_track_lower_bound ) & RK_SERIAL_SAVE_WITH_NAME( m_track_upper_bound )
    & RK_SERIAL_SAVE_WITH_NAME( m_arm_model ) & RK_SERIAL_SAVE_WITH_NAME( m_chain );
};

void RK_CALL manip_P3R3R_kinematics::load( serialization::iarchive& A, unsigned int ) {
  inverse_kinematics_model::load( A, inverse_kinematics_model::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_LOAD_WITH_NAME( m_base_frame ) & RK_SERIAL_LOAD_WITH_NAME( m_track_coord )
    & RK_SERIAL_LOAD_WITH_NAME( m_output_frame ) & RK_SERIAL_LOAD_WITH_NAME( m_track_joint )
    & RK_SERIAL_LOAD_WITH_NAME( m_track_lower_bound ) & RK_SERIAL_LOAD_WITH_NAME( m_track_upper_bound )
    & RK_SERIAL_LOAD_WITH_NAME( m_arm_model ) & RK_SERIAL_LOAD_WITH_NAME( m_chain );
};
};
};
