/**
 * \file manipulator_model_helper.hpp
 *
 * This library declares helper classes to use manipulator models, both kinematic only or
 * dynamic as well. Essentially, the model of the manipulator is only a KTE chain provided
 * by the user, but these manipulator-model classes take care of grouping the joints, their
 * limits, and their jacobian matrices.
 * The helper classes provided here are friends of the manipulator model classes and implement
 * some details of their functionalities. The helper classes can also be used elsewhere to access
 * lower-level functions.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2012
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_MANIPULATOR_MODEL_HELPER_HPP
#define REAK_MANIPULATOR_MODEL_HELPER_HPP

#include "manipulator_model.hpp"

#include <ReaK/math/optimization/nl_interior_points_methods.hpp>
#include <ReaK/math/optimization/function_types.hpp>

namespace ReaK {

namespace kte {


class manip_kin_mdl_joint_io {
public:
  /**
   * Default constructor.
   */
  manip_kin_mdl_joint_io( const manipulator_kinematics_model* aModel ) : model( aModel ){};

  /**
   * Default destructor.
   */
  ~manip_kin_mdl_joint_io(){};


private:
  const manipulator_kinematics_model* model;

public:
  template < typename Vector >
  void getJointPositions( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      result[j] = ( *it )->q;

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      result[j] = ( *it )->Position[0];
      ++j;
      result[j] = ( *it )->Position[1];
      ++j;
      result[j] = ( *it )->Rotation[0];
      ++j;
      result[j] = ( *it )->Rotation[1];
      ++j;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      result[j] = ( *it )->Position[0];
      ++j;
      result[j] = ( *it )->Position[1];
      ++j;
      result[j] = ( *it )->Position[2];
      ++j;
      result[j] = ( *it )->Quat[0];
      ++j;
      result[j] = ( *it )->Quat[1];
      ++j;
      result[j] = ( *it )->Quat[2];
      ++j;
      result[j] = ( *it )->Quat[3];
      ++j;
    };
  };

  template < typename Vector >
  void setJointPositions( const Vector& aJointPositions ) {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      ( *it )->q = aJointPositions[j];

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      ( *it )->Position[0] = aJointPositions[j];
      ++j;
      ( *it )->Position[1] = aJointPositions[j];
      ++j;
      ( *it )->Rotation = rot_mat_2D< double >( vect< double, 2 >( aJointPositions[j], aJointPositions[j + 1] ) );
      j += 2;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      ( *it )->Position[0] = aJointPositions[j];
      ++j;
      ( *it )->Position[1] = aJointPositions[j];
      ++j;
      ( *it )->Position[2] = aJointPositions[j];
      ++j;
      ( *it )->Quat = quaternion< double >( vect< double, 4 >( aJointPositions[j], aJointPositions[j + 1],
                                                               aJointPositions[j + 2], aJointPositions[j + 3] ) );
      j += 4;
    };
  };

  template < typename Vector >
  void getJointVelocities( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      result[j] = ( *it )->q_dot;

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      result[j] = ( *it )->Velocity[0];
      ++j;
      result[j] = ( *it )->Velocity[1];
      ++j;
      result[j] = ( *it )->AngVelocity;
      ++j;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      result[j] = ( *it )->Velocity[0];
      ++j;
      result[j] = ( *it )->Velocity[1];
      ++j;
      result[j] = ( *it )->Velocity[2];
      ++j;
      result[j] = ( *it )->AngVelocity[0];
      ++j;
      result[j] = ( *it )->AngVelocity[1];
      ++j;
      result[j] = ( *it )->AngVelocity[2];
      ++j;
    };
  };

  template < typename Vector >
  void setJointVelocities( const Vector& aJointVelocities ) {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      ( *it )->q_dot = aJointVelocities[j];

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      ( *it )->Velocity[0] = aJointVelocities[j];
      ++j;
      ( *it )->Velocity[1] = aJointVelocities[j];
      ++j;
      ( *it )->AngVelocity = aJointVelocities[j];
      ++j;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      ( *it )->Velocity[0] = aJointVelocities[j];
      ++j;
      ( *it )->Velocity[1] = aJointVelocities[j];
      ++j;
      ( *it )->Velocity[2] = aJointVelocities[j];
      ++j;
      ( *it )->AngVelocity[0] = aJointVelocities[j];
      ++j;
      ( *it )->AngVelocity[1] = aJointVelocities[j];
      ++j;
      ( *it )->AngVelocity[2] = aJointVelocities[j];
      ++j;
    };
  };

  template < typename Vector >
  void getJointAccelerations( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      result[j] = ( *it )->q_ddot;

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      result[j] = ( *it )->Acceleration[0];
      ++j;
      result[j] = ( *it )->Acceleration[1];
      ++j;
      result[j] = ( *it )->AngAcceleration;
      ++j;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      result[j] = ( *it )->Acceleration[0];
      ++j;
      result[j] = ( *it )->Acceleration[1];
      ++j;
      result[j] = ( *it )->Acceleration[2];
      ++j;
      result[j] = ( *it )->AngAcceleration[0];
      ++j;
      result[j] = ( *it )->AngAcceleration[1];
      ++j;
      result[j] = ( *it )->AngAcceleration[2];
      ++j;
    };
  };

  template < typename Vector >
  void setJointAccelerations( const Vector& aJointAccelerations ) {
    std::size_t j = 0;

    for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = model->mCoords.begin();
         it < model->mCoords.end(); ++it, ++j )
      ( *it )->q_ddot = aJointAccelerations[j];

    for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = model->mFrames2D.begin();
         it < model->mFrames2D.end(); ++it ) {
      ( *it )->Acceleration[0] = aJointAccelerations[j];
      ++j;
      ( *it )->Acceleration[1] = aJointAccelerations[j];
      ++j;
      ( *it )->AngAcceleration = aJointAccelerations[j];
      ++j;
    };

    for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = model->mFrames3D.begin();
         it < model->mFrames3D.end(); ++it ) {
      ( *it )->Acceleration[0] = aJointAccelerations[j];
      ++j;
      ( *it )->Acceleration[1] = aJointAccelerations[j];
      ++j;
      ( *it )->Acceleration[2] = aJointAccelerations[j];
      ++j;
      ( *it )->AngAcceleration[0] = aJointAccelerations[j];
      ++j;
      ( *it )->AngAcceleration[1] = aJointAccelerations[j];
      ++j;
      ( *it )->AngAcceleration[2] = aJointAccelerations[j];
      ++j;
    };
  };

  template < typename Vector >
  void getDependentPositions( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = model->mDependentGenCoords.begin();
         it < model->mDependentGenCoords.end(); ++it, ++j )
      result[j] = ( *it )->mFrame->q;

    for( std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = model->mDependent2DFrames.begin();
         it < model->mDependent2DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Position[0];
      ++j;
      result[j] = ( *it )->mFrame->Position[1];
      ++j;
      result[j] = ( *it )->mFrame->Rotation[0];
      ++j;
      result[j] = ( *it )->mFrame->Rotation[1];
      ++j;
    };

    for( std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = model->mDependent3DFrames.begin();
         it < model->mDependent3DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Position[0];
      ++j;
      result[j] = ( *it )->mFrame->Position[1];
      ++j;
      result[j] = ( *it )->mFrame->Position[2];
      ++j;
      result[j] = ( *it )->mFrame->Quat[0];
      ++j;
      result[j] = ( *it )->mFrame->Quat[1];
      ++j;
      result[j] = ( *it )->mFrame->Quat[2];
      ++j;
      result[j] = ( *it )->mFrame->Quat[3];
      ++j;
    };
  };

  template < typename Vector >
  void getDependentVelocities( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = model->mDependentGenCoords.begin();
         it < model->mDependentGenCoords.end(); ++it, ++j )
      result[j] = ( *it )->mFrame->q_dot;

    for( std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = model->mDependent2DFrames.begin();
         it < model->mDependent2DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Velocity[0];
      ++j;
      result[j] = ( *it )->mFrame->Velocity[1];
      ++j;
      result[j] = ( *it )->mFrame->AngVelocity;
      ++j;
    };

    for( std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = model->mDependent3DFrames.begin();
         it < model->mDependent3DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Velocity[0];
      ++j;
      result[j] = ( *it )->mFrame->Velocity[1];
      ++j;
      result[j] = ( *it )->mFrame->Velocity[2];
      ++j;
      result[j] = ( *it )->mFrame->AngVelocity[0];
      ++j;
      result[j] = ( *it )->mFrame->AngVelocity[1];
      ++j;
      result[j] = ( *it )->mFrame->AngVelocity[2];
      ++j;
    };
  };

  template < typename Vector >
  void getDependentAccelerations( Vector& result ) const {
    std::size_t j = 0;

    for( std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = model->mDependentGenCoords.begin();
         it < model->mDependentGenCoords.end(); ++it, ++j )
      result[j] = ( *it )->mFrame->q_ddot;

    for( std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = model->mDependent2DFrames.begin();
         it < model->mDependent2DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Acceleration[0];
      ++j;
      result[j] = ( *it )->mFrame->Acceleration[1];
      ++j;
      result[j] = ( *it )->mFrame->AngAcceleration;
      ++j;
    };

    for( std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = model->mDependent3DFrames.begin();
         it < model->mDependent3DFrames.end(); ++it ) {
      result[j] = ( *it )->mFrame->Acceleration[0];
      ++j;
      result[j] = ( *it )->mFrame->Acceleration[1];
      ++j;
      result[j] = ( *it )->mFrame->Acceleration[2];
      ++j;
      result[j] = ( *it )->mFrame->AngAcceleration[0];
      ++j;
      result[j] = ( *it )->mFrame->AngAcceleration[1];
      ++j;
      result[j] = ( *it )->mFrame->AngAcceleration[2];
      ++j;
    };
  };
};


/**
 * This class is a helper of the manipulator_kinematics_model class which is used to fill in
 * the Jacobian matrix (and its time-derivative). This class is useful because it is a friend
 * of the manipulator_kinematics_model class (thus, has access to its data members), but also
 * provides member function templates for filling in the matrices, meaning it can be used to
 * fill in any kind of matrix type (e.g. enabling the filling of matrix sub-blocks for example).
 */
class manip_kin_mdl_jac_calculator {
public:
  /**
   * Default constructor.
   */
  manip_kin_mdl_jac_calculator( const manipulator_kinematics_model* aModel ) : model( aModel ){};

  /**
   * Default destructor.
   */
  ~manip_kin_mdl_jac_calculator(){};

  /**
   * Get the Jacobian matrix for the system (or twist-shaping matrix). The Jacobian takes the velocity
   * information of the system coordinates and frames, and maps them to velocity information
   * of the system's dependent coordinates and frames.
   * \param Jac stores, as output, the calculated system's Jacobian matrix.
   */
  template < typename Matrix1 >
  void getJacobianMatrix( Matrix1& Jac ) const {
    getJacobianMatrixAndDerivativeImpl( Jac, static_cast< mat< double, mat_structure::rectangular >* >( nullptr ) );
  };

  /**
   * Get the Jacobian matrix for the system (or twist-shaping matrix), and its time-derivative.
   * The Jacobian takes the velocity information of the system coordinates and frames, and maps
   * them to velocity information of the system's dependent coordinates and frames. The time-derivative
   * of the Jacobian matrix will map the velocity information of the system coordinates and frames
   * to the acceleration information of the system's dependent coordinates and frames.
   * \param Jac stores, as output, the calculated system's Jacobian matrix.
   * \param JacDot stores, as output, the calculated time-derivative of the system's Jacobian matrix.
   */
  template < typename Matrix1, typename Matrix2 >
  void getJacobianMatrixAndDerivative( Matrix1& Jac, Matrix2& JacDot ) const {
    getJacobianMatrixAndDerivativeImpl( Jac, &JacDot );
  };

private:
  const manipulator_kinematics_model* model;

  template < typename Matrix1, typename Matrix2 >
  void getJacobianMatrixAndDerivativeImpl( Matrix1& Jac, Matrix2* JacDot ) const {
    std::size_t m = model->getDependentVelocitiesCount();
    std::size_t n = model->getJointVelocitiesCount();
    Jac = mat< double, mat_structure::nil >( m, n );
    if( JacDot )
      *JacDot = mat< double, mat_structure::nil >( m, n );

    std::size_t RowInd = 0;

    /****************************************************************************************
     *                             Gen Coords
     * *************************************************************************************/


    for (std::size_t i = 0; i < model->mCoords.size(); ++i) {
      RowInd = 0;

      for (std::size_t j = 0; j < model->mDependentGenCoords.size(); ++j) {
        if( model->mDependentGenCoords[j]->mUpStreamJoints.find( model->mCoords[i] )
            != model->mDependentGenCoords[j]->mUpStreamJoints.end() ) {
          mat_sub_block< Matrix1 > subJac = sub( Jac )( range( RowInd, RowInd + 1 ), range( i, i + 1 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot = sub ( *JacDot )( range( RowInd, RowInd + 1 ), range( i, i + 1 ) );
            model->mDependentGenCoords[j]->mUpStreamJoints[model->mCoords[i]]->write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependentGenCoords[j]->mUpStreamJoints[model->mCoords[i]]->write_to_matrices( subJac );
          };
        };
        RowInd++;
      };

      for (std::size_t j = 0; j < model->mDependent2DFrames.size(); ++j) {
        if( model->mDependent2DFrames[j]->mUpStreamJoints.find( model->mCoords[i] )
            != model->mDependent2DFrames[j]->mUpStreamJoints.end() ) {
          mat_sub_block< Matrix1 > subJac = sub( Jac )( range( RowInd, RowInd + 3 ), range( i, i + 1 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot = sub ( *JacDot )( range( RowInd, RowInd + 3 ), range( i, i + 1 ) );
            model->mDependent2DFrames[j]
              ->mUpStreamJoints[model->mCoords[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent2DFrames[j]
              ->mUpStreamJoints[model->mCoords[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 3;
      };

      for (std::size_t j = 0; j < model->mDependent3DFrames.size(); ++j) {
        if( model->mDependent3DFrames[j]->mUpStreamJoints.find( model->mCoords[i] )
            != model->mDependent3DFrames[j]->mUpStreamJoints.end() ) {
          mat_sub_block< Matrix1 > subJac = sub( Jac )( range( RowInd, RowInd + 6 ), range( i, i + 1 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot = sub ( *JacDot )( range( RowInd, RowInd + 6 ), range( i, i + 1 ) );
            model->mDependent3DFrames[j]
              ->mUpStreamJoints[model->mCoords[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent3DFrames[j]
              ->mUpStreamJoints[model->mCoords[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 6;
      };
    };


    /****************************************************************************************
     *                             2D Frames
     * *************************************************************************************/

    std::size_t base_i = model->mCoords.size();
    for (std::size_t i = 0; i < model->mFrames2D.size(); ++i) {
      RowInd = 0;

      for (std::size_t j = 0; j < model->mDependentGenCoords.size(); ++j) {
        if( model->mDependentGenCoords[j]->mUpStream2DJoints.find( model->mFrames2D[i] )
            != model->mDependentGenCoords[j]->mUpStream2DJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 1 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 1 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
            model->mDependentGenCoords[j]->mUpStream2DJoints[model->mFrames2D[i]]->write_to_matrices( subJac,
                                                                                                      subJacDot );
          } else {
            model->mDependentGenCoords[j]->mUpStream2DJoints[model->mFrames2D[i]]->write_to_matrices( subJac );
          };
        };
        RowInd++;
      };

      for (std::size_t j = 0; j < model->mDependent2DFrames.size(); ++j) {
        if( model->mDependent2DFrames[j]->mUpStream2DJoints.find( model->mFrames2D[i] )
            != model->mDependent2DFrames[j]->mUpStream2DJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 3 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 3 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
            model->mDependent2DFrames[j]
              ->mUpStream2DJoints[model->mFrames2D[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent2DFrames[j]
              ->mUpStream2DJoints[model->mFrames2D[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 3;
      };

      for (std::size_t j = 0; j < model->mDependent3DFrames.size(); ++j) {
        if( model->mDependent3DFrames[j]->mUpStream2DJoints.find( model->mFrames2D[i] )
            != model->mDependent3DFrames[j]->mUpStream2DJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 6 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 6 ), range( 3 * i + base_i, 3 * i + base_i + 3 ) );
            model->mDependent3DFrames[j]
              ->mUpStream2DJoints[model->mFrames2D[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent3DFrames[j]
              ->mUpStream2DJoints[model->mFrames2D[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 6;
      };
    };


    /****************************************************************************************
     *                             3D Frames
     * *************************************************************************************/

    base_i = model->mCoords.size() + 3 * model->mFrames2D.size();
    for (std::size_t i = 0; i < model->mFrames3D.size(); ++i) {
      RowInd = 0;

      for (std::size_t j = 0; j < model->mDependentGenCoords.size(); ++j) {
        if( model->mDependentGenCoords[j]->mUpStreamJoints.find( model->mCoords[i] )
            != model->mDependentGenCoords[j]->mUpStreamJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 1 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 1 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
            model->mDependentGenCoords[j]->mUpStream3DJoints[model->mFrames3D[i]]->write_to_matrices( subJac,
                                                                                                      subJacDot );
          } else {
            model->mDependentGenCoords[j]->mUpStream3DJoints[model->mFrames3D[i]]->write_to_matrices( subJac );
          };
        };
        RowInd++;
      };

      for (std::size_t j = 0; j < model->mDependent2DFrames.size(); ++j) {
        if( model->mDependent2DFrames[j]->mUpStream3DJoints.find( model->mFrames3D[i] )
            != model->mDependent2DFrames[j]->mUpStream3DJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 3 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 3 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
            model->mDependent2DFrames[j]
              ->mUpStream3DJoints[model->mFrames3D[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent2DFrames[j]
              ->mUpStream3DJoints[model->mFrames3D[i]]
              ->get_jac_relative_to( model->mDependent2DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 3;
      };

      for (std::size_t j = 0; j < model->mDependent3DFrames.size(); ++j) {
        if( model->mDependent3DFrames[j]->mUpStream3DJoints.find( model->mFrames3D[i] )
            != model->mDependent3DFrames[j]->mUpStream3DJoints.end() ) {
          mat_sub_block< Matrix1 > subJac
            = sub( Jac )( range( RowInd, RowInd + 6 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
          if( JacDot ) {
            mat_sub_block< Matrix2 > subJacDot
              = sub ( *JacDot )( range( RowInd, RowInd + 6 ), range( 6 * i + base_i, 6 * i + base_i + 6 ) );
            model->mDependent3DFrames[j]
              ->mUpStream3DJoints[model->mFrames3D[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac, subJacDot );
          } else {
            model->mDependent3DFrames[j]
              ->mUpStream3DJoints[model->mFrames3D[i]]
              ->get_jac_relative_to( model->mDependent3DFrames[j]->mFrame )
              .write_to_matrices( subJac );
          };
        };
        RowInd += 6;
      };
    };
  };
};


/**
 * This class is a helper of the manipulator_kinematics_model class which is used to perform
 * a closed-loop inverse kinematics (CLIK) calculation to compute the joint positions and
 * velocities necessary for a given set of dependent positions and velocities (end-effector).
 * The inverse kinematics calculation is done using a non-linear constrained optimization
 * method (the nlip_newton_tr_factory method, which is a non-linear interior-point Newton
 * method based on a trust-region search strategy, this method has shown rapid convergence for
 * inverse kinematics problems). This class automatically builds inequality constraints based on
 * the joint limits provided. Then, the class builds equality constraints such that the end-effector
 * pose and twist matches the desired values. Finally, one can specify a cost function to be minimized
 * if there is sufficient redundancy to permit such optimization of the resulting configuration.
 */
class manip_clik_calculator {
private:
  manipulator_kinematics_model* model;

public:
  shared_ptr< const optim::cost_evaluator > cost_eval;

  vect_n< double > lower_bounds;
  vect_n< double > upper_bounds;

  std::vector< gen_coord< double > > desired_gen_coords;
  std::vector< frame_2D< double > > desired_frame_2D;
  std::vector< frame_3D< double > > desired_frame_3D;


  struct ineq_function {
    const manip_clik_calculator* parent;

    ineq_function( const manip_clik_calculator* aParent ) : parent( aParent ){};

    vect_n< double > operator()( const vect_n< double >& x ) const {
      std::size_t l_size = 0;
      for( std::size_t i = 0; i < parent->lower_bounds.size(); ++i ) {
        if( parent->lower_bounds[i] != -std::numeric_limits< double >::infinity() )
          ++l_size;
      };
      std::size_t u_size = 0;
      for( std::size_t i = 0; i < parent->upper_bounds.size(); ++i ) {
        if( parent->upper_bounds[i] != std::numeric_limits< double >::infinity() )
          ++u_size;
      };
      vect_n< double > result( l_size + u_size );
      std::size_t j = 0;
      for( std::size_t i = 0; ( i < parent->lower_bounds.size() ) && ( i < x.size() ); ++i ) {
        if( parent->lower_bounds[i] != -std::numeric_limits< double >::infinity() ) {
          result[j] = x[i] - parent->lower_bounds[i];
          ++j;
        };
      };
      for( std::size_t i = 0; ( i < parent->upper_bounds.size() ) && ( i < x.size() ); ++i ) {
        if( parent->upper_bounds[i] != std::numeric_limits< double >::infinity() ) {
          result[j] = parent->upper_bounds[i] - x[i];
          ++j;
        };
      };
      return result;
    };
  };

  struct ineq_jac_filler {
    const manip_clik_calculator* parent;

    ineq_jac_filler( const manip_clik_calculator* aParent ) : parent( aParent ){};

    template < typename Matrix >
    void operator()( Matrix& J, const vect_n< double >& x, const vect_n< double >& h ) const {
      J = mat< double, mat_structure::nil >( h.size(), x.size() );
      std::size_t j = 0;
      for( std::size_t i = 0; ( i < parent->lower_bounds.size() ) && ( i < x.size() ); ++i ) {
        if( parent->lower_bounds[i] != -std::numeric_limits< double >::infinity() ) {
          J( j, i ) = 1.0;
          ++j;
        };
      };
      for( std::size_t i = 0; ( i < parent->upper_bounds.size() ) && ( i < x.size() ); ++i ) {
        if( parent->upper_bounds[i] != std::numeric_limits< double >::infinity() ) {
          J( j, i ) = -1.0;
          ++j;
        };
      };
    };
  };

  struct eq_function {
    const manip_clik_calculator* parent;

    eq_function( const manip_clik_calculator* aParent ) : parent( aParent ){};

    vect_n< double > operator()( const vect_n< double >& x ) const {

      const std::vector< shared_ptr< joint_dependent_gen_coord > >& dep_gen_coords = parent->model->DependentCoords();
      const std::vector< shared_ptr< joint_dependent_frame_2D > >& dep_frames_2D = parent->model->DependentFrames2D();
      const std::vector< shared_ptr< joint_dependent_frame_3D > >& dep_frames_3D = parent->model->DependentFrames3D();

      if( ( dep_gen_coords.size() != parent->desired_gen_coords.size() )
          || ( dep_frames_2D.size() != parent->desired_frame_2D.size() )
          || ( dep_frames_3D.size() != parent->desired_frame_3D.size() ) )
        throw std::range_error( "Improper inverse-kinematics problem, the number of desired frames does not match the "
                                "number of end-effector frames!" );

      manip_kin_mdl_joint_io( parent->model )
        .setJointPositions( x[range( 0, parent->model->getJointPositionsCount() )] );
      manip_kin_mdl_joint_io( parent->model ).setJointVelocities(
        x[range( parent->model->getJointPositionsCount(),
                 parent->model->getJointPositionsCount() + parent->model->getJointVelocitiesCount() )] );

      parent->model->doMotion();

      vect_n< double > result( parent->model->getDependentVelocitiesCount() * 2 + parent->model->Frames2D().size()
                               + parent->model->Frames3D().size() );


      // enforce the desired 'end-effector' frames.
      std::size_t j = 0;
      std::size_t k = parent->model->getDependentVelocitiesCount();

      for( std::size_t i = 0; i < dep_gen_coords.size(); ++i ) {
        result[j] = dep_gen_coords[i]->mFrame->q - parent->desired_gen_coords[i].q;
        ++j;
        result[k] = dep_gen_coords[i]->mFrame->q_dot - parent->desired_gen_coords[i].q_dot;
        ++k;
      };

      for( std::size_t i = 0; i < dep_frames_2D.size(); ++i ) {
        frame_2D< double > err = parent->desired_frame_2D[i].getFrameRelativeTo( dep_frames_2D[i]->mFrame );
        result[j] = -err.Position[0];
        ++j;
        result[j] = -err.Position[1];
        ++j;
        result[j] = -err.Rotation.getAngle();
        ++j;
        result[k] = -err.Velocity[0];
        ++k;
        result[k] = -err.Velocity[1];
        ++k;
        result[k] = -err.AngVelocity;
        ++k;
      };

      for( std::size_t i = 0; i < dep_frames_3D.size(); ++i ) {
        frame_3D< double > err = parent->desired_frame_3D[i].getFrameRelativeTo( dep_frames_3D[i]->mFrame );
        result[j] = -err.Position[0];
        ++j;
        result[j] = -err.Position[1];
        ++j;
        result[j] = -err.Position[2];
        ++j;
        axis_angle< double > aa = axis_angle< double >( err.Quat );
        vect< double, 3 > v = aa.angle() * aa.axis();
        result[j] = -v[0];
        ++j;
        result[j] = -v[1];
        ++j;
        result[j] = -v[2];
        ++j;
        result[k] = -err.Velocity[0];
        ++k;
        result[k] = -err.Velocity[1];
        ++k;
        result[k] = -err.Velocity[2];
        ++k;
        result[k] = -err.AngVelocity[0];
        ++k;
        result[k] = -err.AngVelocity[1];
        ++k;
        result[k] = -err.AngVelocity[2];
        ++k;
      };

      // enforce the normality of the rotation representation.
      j = parent->model->Coords().size();
      for( std::size_t i = 0; i < parent->model->Frames2D().size(); ++i ) {
        j += 2;
        result[k] = 1.0 - x[j] * x[j] - x[j + 1] * x[j + 1];
        ++k;
        j += 2;
      };

      for( std::size_t i = 0; i < parent->model->Frames3D().size(); ++i ) {
        j += 3;
        result[k] = 1.0 - x[j] * x[j] - x[j + 1] * x[j + 1] - x[j + 2] * x[j + 2] - x[j + 3] * x[j + 3];
        ++k;
        j += 4;
      };

      return result;
    };
  };

  struct eq_jac_filler {
    const manip_clik_calculator* parent;

    eq_jac_filler( const manip_clik_calculator* aParent ) : parent( aParent ){};

    template < typename Matrix >
    void operator()( Matrix& J, const vect_n< double >& x, const vect_n< double >& h ) const {

      const std::vector< shared_ptr< joint_dependent_gen_coord > >& dep_gen_coords = parent->model->DependentCoords();
      const std::vector< shared_ptr< joint_dependent_frame_2D > >& dep_frames_2D = parent->model->DependentFrames2D();
      const std::vector< shared_ptr< joint_dependent_frame_3D > >& dep_frames_3D = parent->model->DependentFrames3D();

      if( ( dep_gen_coords.size() != parent->desired_gen_coords.size() )
          || ( dep_frames_2D.size() != parent->desired_frame_2D.size() )
          || ( dep_frames_3D.size() != parent->desired_frame_3D.size() ) )
        throw std::range_error( "Improper inverse-kinematics problem, the number of desired frames does not match the "
                                "number of end-effector frames!" );

      manip_kin_mdl_joint_io( parent->model )
        .setJointPositions( x[range( 0, parent->model->getJointPositionsCount() )] );
      manip_kin_mdl_joint_io( parent->model ).setJointVelocities(
        x[range( parent->model->getJointPositionsCount(),
                 parent->model->getJointPositionsCount() + parent->model->getJointVelocitiesCount() )] );

      parent->model->doMotion();

      J = mat< double, mat_structure::nil >( h.size(), x.size() );

      mat_sub_block< Matrix > Jac_tmp = sub( J )(
        range( 0, parent->model->getDependentVelocitiesCount() ),
        range( 0, ( x.size() - parent->model->Frames2D().size() - parent->model->Frames3D().size() ) / 2 ) );
      mat_sub_block< Matrix > JacDot_tmp = sub( J )(
        range( parent->model->getDependentVelocitiesCount(), parent->model->getDependentVelocitiesCount() * 2 ),
        range( 0, ( x.size() - parent->model->Frames2D().size() - parent->model->Frames3D().size() ) / 2 ) );

      manip_kin_mdl_jac_calculator( parent->model ).getJacobianMatrixAndDerivative( Jac_tmp, JacDot_tmp );

      sub( J )( range( parent->model->getDependentVelocitiesCount(), parent->model->getDependentVelocitiesCount() * 2 ),
                range( ( x.size() + parent->model->Frames2D().size() + parent->model->Frames3D().size() ) / 2,
                       x.size() ) ) = Jac_tmp;

      // j is the index to the last position element of x.
      std::size_t j = ( x.size() + parent->model->Frames2D().size() + parent->model->Frames3D().size() ) / 2 - 1;
      // k is the index to the last valid column of J.
      std::size_t k = ( x.size() - parent->model->Frames2D().size() - parent->model->Frames3D().size() ) / 2 - 1;
      // l is the index to the last normality-constraint row of J.
      std::size_t l = h.size() - 1;

      for( std::size_t i = 0; i < parent->model->Frames3D().size(); ++i ) {
        mat< double, mat_structure::rectangular > H_inv(
          3, 4 ); // NOTE : Check this again, shouldn't there be a factor of 2 or 0.5 ???
        H_inv( 0, 0 ) = -x[j - 2];
        H_inv( 1, 0 ) = -x[j - 1];
        H_inv( 2, 0 ) = -x[j];
        H_inv( 0, 1 ) = x[j - 3];
        H_inv( 1, 1 ) = -x[j];
        H_inv( 2, 1 ) = x[j - 1];
        H_inv( 0, 2 ) = x[j];
        H_inv( 1, 2 ) = x[j - 3];
        H_inv( 2, 2 ) = -x[j - 2];
        H_inv( 0, 3 ) = -x[j - 1];
        H_inv( 1, 3 ) = x[j - 2];
        H_inv( 2, 3 ) = x[j - 3];

        // apply the transformation from quat_dot to omega:
        sub( J )( range( 0, parent->model->getDependentVelocitiesCount() * 2 ), range( j - 3, j + 1 ) )
          = sub( J )( range( 0, parent->model->getDependentVelocitiesCount() * 2 ), range( k - 2, k + 1 ) ) * H_inv;
        // fill in the normality-constraint jacobians:
        J( l, j - 3 ) = -2.0 * x[j - 3];
        J( l, j - 2 ) = -2.0 * x[j - 2];
        J( l, j - 1 ) = -2.0 * x[j - 1];
        J( l, j ) = -2.0 * x[j];
        k -= 3;
        j -= 4;
        --l;
        // copy the position row:
        for( std::size_t r = 0; r < 3; ++r )
          for( std::size_t s = 0; s < parent->model->getDependentVelocitiesCount() * 2; ++s )
            J( s, j - r ) = J( s, k - r );
        k -= 3;
        j -= 3;
      };

      for( std::size_t i = 0; i < parent->model->Frames2D().size(); ++i ) {
        // apply the transformation from quat_dot to omega:
        for( std::size_t s = 0; s < parent->model->getDependentVelocitiesCount() * 2; ++s ) {
          J( s, j - 1 ) = -J( s, k ) * x[j];
          J( s, j ) = J( s, k ) * x[j - 1];
        };
        // fill in the normality-constraint jacobians:
        J( l, j - 1 ) = -2.0 * x[j - 1];
        J( l, j ) = -2.0 * x[j];
        k -= 1;
        j -= 2;
        --l;
        // copy the position row:
        for( std::size_t r = 0; r < 2; ++r )
          for( std::size_t s = 0; s < parent->model->getDependentVelocitiesCount() * 2; ++s )
            J( s, j - r ) = J( s, k - r );
        k -= 2;
        i -= 2;
      };
    };
  };


  typedef optim::nlip_newton_tr_factory< optim::oop_cost_function, optim::oop_cost_grad, optim::oop_cost_hess, double,
                                         eq_function, eq_jac_filler, ineq_function,
                                         ineq_jac_filler > optim_factory_type;
  //     typedef optim::nlip_quasi_newton_tr_factory<optim::oop_cost_function,
  //                                                 optim::oop_cost_grad,
  //                                                 double,
  //                                                 eq_function,
  //                                                 eq_jac_filler,
  //                                                 ineq_function,
  //                                                 ineq_jac_filler> optim_factory_type;

  optim_factory_type optimizer;

  /**
   * Default constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  manip_clik_calculator( manipulator_kinematics_model* aModel,
                         const shared_ptr< const optim::cost_evaluator >& aCostEvaluator
                         = shared_ptr< const optim::cost_evaluator >(),
                         double aMaxRadius = 1.0, double aMu = 0.1, std::size_t aMaxIter = 300, double aTol = 1e-6,
                         double aEta = 1e-3, double aTau = 0.99 )
      : model( aModel ), cost_eval( aCostEvaluator ),
        optimizer( optim::oop_cost_function( aCostEvaluator ), optim::oop_cost_grad( aCostEvaluator ),
                   optim::oop_cost_hess( aCostEvaluator ), aMaxRadius, aMu, aMaxIter, eq_function( nullptr ),
                   eq_jac_filler( nullptr ), ineq_function( nullptr ), ineq_jac_filler( nullptr ), aTol, aEta, aTau ) {
    if( !model )
      throw optim::improper_problem( "CLIK error: The model pointer cannot be null!" );

    lower_bounds.resize( model->getJointPositionsCount() + model->getJointVelocitiesCount() );
    upper_bounds.resize( model->getJointPositionsCount() + model->getJointVelocitiesCount() );
    for( std::size_t i = 0; i < lower_bounds.size(); ++i ) {
      lower_bounds[i] = -std::numeric_limits< double >::infinity();
      upper_bounds[i] = std::numeric_limits< double >::infinity();
    };
  };


  /**
   * Default destructor.
   */
  ~manip_clik_calculator(){};

  void readDesiredFromModel() {

    desired_gen_coords.resize( model->DependentCoords().size() );
    for( std::size_t i = 0; i < desired_gen_coords.size(); ++i )
      desired_gen_coords[i] = *( model->DependentCoords()[i]->mFrame );

    desired_frame_2D.resize( model->DependentFrames2D().size() );
    for( std::size_t i = 0; i < desired_frame_2D.size(); ++i )
      desired_frame_2D[i] = *( model->DependentFrames2D()[i]->mFrame );

    desired_frame_3D.resize( model->DependentFrames3D().size() );
    for( std::size_t i = 0; i < desired_frame_3D.size(); ++i )
      desired_frame_3D[i] = *( model->DependentFrames3D()[i]->mFrame );
  };


  vect_n< double > readJointStatesFromModel() {
    vect_n< double > x( model->getJointPositionsCount() + model->getJointVelocitiesCount() );

    vect_ref_view< vect_n< double > > pos_x = x[range( 0, model->getJointPositionsCount() )];
    manip_kin_mdl_joint_io( model ).getJointPositions( pos_x );

    vect_ref_view< vect_n< double > > vel_x
      = x[range( model->getJointPositionsCount(), model->getJointPositionsCount() + model->getJointVelocitiesCount() )];
    manip_kin_mdl_joint_io( model ).getJointVelocities( vel_x );

    return x;
  };

  void writeJointStatesToModel( const vect_n< double >& x ) const {
    manip_kin_mdl_joint_io( model ).setJointPositions( x[range( 0, model->getJointPositionsCount() )] );
    manip_kin_mdl_joint_io( model ).setJointVelocities(
      x[range( model->getJointPositionsCount(), model->getJointPositionsCount() + model->getJointVelocitiesCount() )] );
  };

  void runOptimizer( vect_n< double >& x ) {
    shared_ptr< const optim::cost_evaluator > tmp_cost_eval = cost_eval;
    if( !cost_eval ) {
      tmp_cost_eval = shared_ptr< const optim::cost_evaluator >( new optim::quadratic_cost_evaluator(
        vect_n< double >( x.size(), double( 0.0 ) ),
        mat< double, mat_structure::symmetric >(
          ( mat< double, mat_structure::nil >( model->getJointPositionsCount(),
                                               model->getJointPositionsCount() + model->getJointVelocitiesCount() )
            | ( mat< double, mat_structure::nil >( model->getJointVelocitiesCount(), model->getJointPositionsCount() )
                & mat< double, mat_structure::identity >( model->getJointVelocitiesCount() ) ) ) ) ) );
    };

    optimizer.f = optim::oop_cost_function( tmp_cost_eval );
    optimizer.df = optim::oop_cost_grad( tmp_cost_eval );
    optimizer.fill_hessian = optim::oop_cost_hess( tmp_cost_eval );
    optimizer.g = eq_function( this );
    optimizer.fill_g_jac = eq_jac_filler( this );
    optimizer.h = ineq_function( this );
    optimizer.fill_h_jac = ineq_jac_filler( this );

    optimizer( x );
  };

  /**
   * This function takes its given desired 'end-effector' states and solves the inverse
   * kinematics problem to find the corresponding joint states.
   * \pre The joint states at which the manipulator model is before the function call
   *      will act as the initial guess for the inverse kinematics solution.
   *      Before calling this function, it is assumed that all the parameters of the
   *      optimization have been properly set.
   * \post The joint states which solve the inverse kinematics problem will be set in
   *       the manipulator model.
   */
  void solveInverseKinematics() {

    readDesiredFromModel();

    vect_n< double > x = readJointStatesFromModel();

    runOptimizer( x );

    writeJointStatesToModel( x );
  };
};
};
};

#endif
