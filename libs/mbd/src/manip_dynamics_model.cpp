
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


#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>

#include <ReaK/mbd/models/manip_dynamics_model.hpp>
#include <ReaK/mbd/models/manip_kinematics_helper.hpp>

namespace ReaK {

namespace kte {


manipulator_kinematics_model& manipulator_dynamics_model::
  operator<<( const shared_ptr< gen_coord< double > >& aCoord ) {
  if( aCoord ) {
    mMassCalc << aCoord;
    manipulator_kinematics_model::operator<<( aCoord );
  };
  return *this;
};

manipulator_kinematics_model& manipulator_dynamics_model::
  operator<<( const shared_ptr< frame_2D< double > >& aFrame2D ) {
  if( aFrame2D ) {
    mMassCalc << aFrame2D;
    manipulator_kinematics_model::operator<<( aFrame2D );
  };
  return *this;
};

manipulator_kinematics_model& manipulator_dynamics_model::
  operator<<( const shared_ptr< frame_3D< double > >& aFrame3D ) {
  if( aFrame3D ) {
    mMassCalc << aFrame3D;
    manipulator_kinematics_model::operator<<( aFrame3D );
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator<<( const shared_ptr< inertia_gen >& aInertiaGen ) {
  if( aInertiaGen ) {
    mMassCalc << aInertiaGen;
    *this << aInertiaGen->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator<<( const shared_ptr< inertia_2D >& aInertia2D ) {
  if( aInertia2D ) {
    mMassCalc << aInertia2D;
    *this << aInertia2D->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator<<( const shared_ptr< inertia_3D >& aInertia3D ) {
  if( aInertia3D ) {
    mMassCalc << aInertia3D;
    *this << aInertia3D->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator<<( const shared_ptr< system_input >& aInput ) {
  mInputs.push_back( aInput );
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator<<( const shared_ptr< system_output >& aOutput ) {
  mOutputs.push_back( aOutput );
  return *this;
};


vect_n< double > manipulator_dynamics_model::getJointStates() const {
  vect_n< double > result( getJointStatesCount() );

  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .getJointPositions( &result[0] );
  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .getJointVelocities( &result[getJointPositionsCount()] );

  return result;
};

void manipulator_dynamics_model::setJointStates( const vect_n< double >& aJointStates ) {
  if( aJointStates.size() != getJointStatesCount() )
    throw std::range_error( "Joint-state vector has incorrect dimensions!" );

  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .setJointPositions( &aJointStates[0] );
  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .setJointVelocities( &aJointStates[getJointPositionsCount()] );
};


void RK_CALL manipulator_dynamics_model::computeOutput( double aTime, const ReaK::vect_n< double >& aState,
                                                        ReaK::vect_n< double >& aOutput ) {
  setJointStates( aState );

  mModel->doMotion();
  mModel->clearForce();
  mModel->doForce();

  aOutput.resize( getOutputsCount() );

  unsigned int i = 0;
  for( std::vector< shared_ptr< system_output > >::iterator it = mOutputs.begin(); it != mOutputs.end(); ++it )
    for( unsigned int k = 0; k < ( *it )->getOutputCount(); ++k )
      aOutput[i++] = ( *it )->getOutput( k );
};

void RK_CALL manipulator_dynamics_model::setInput( const ReaK::vect_n< double >& aInput ) {
  if( aInput.size() != getInputsCount() )
    throw std::range_error( "The size of the input-vector to the manipulator model is not correct!" );

  unsigned int i = 0;
  for( std::vector< shared_ptr< system_input > >::iterator it = mInputs.begin(); it != mInputs.end(); ++it )
    for( unsigned int k = 0; k < ( *it )->getInputCount(); ++k )
      ( *it )->setInput( k, aInput[i++] );
};

vect_n< double > manipulator_dynamics_model::getInput() const {
  vect_n< double > result( getInputsCount() );

  unsigned int i = 0;
  for( std::vector< shared_ptr< system_input > >::const_iterator it = mInputs.begin(); it != mInputs.end(); ++it )
    for( unsigned int k = 0; k < ( *it )->getInputCount(); ++k )
      result[i++] = ( *it )->getInput( k );
  return result;
};

void RK_CALL manipulator_dynamics_model::computeStateRate( double aTime, const vect_n< double >& aState,
                                                           vect_n< double >& aStateRate ) {
  setJointStates( aState );

  mModel->doMotion();
  mModel->clearForce();
  mModel->doForce();

  aStateRate.resize( getJointStatesCount() );

  unsigned int j = 0;
  for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = mCoords.begin(); it < mCoords.end();
       ++it, ++j )
    aStateRate[j] = ( *it )->q_dot;

  for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = mFrames2D.begin(); it < mFrames2D.end();
       ++it ) {
    aStateRate[j] = ( *it )->Velocity[0];
    ++j;
    aStateRate[j] = ( *it )->Velocity[1];
    ++j;
    aStateRate[j] = -( *it )->Rotation[1] * ( *it )->AngVelocity;
    ++j;
    aStateRate[j] = ( *it )->Rotation[0] * ( *it )->AngVelocity;
    ++j;
  };

  for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = mFrames3D.begin(); it < mFrames3D.end();
       ++it ) {
    aStateRate[j] = ( *it )->Velocity[0];
    ++j;
    aStateRate[j] = ( *it )->Velocity[1];
    ++j;
    aStateRate[j] = ( *it )->Velocity[2];
    ++j;
    aStateRate[j] = ( *it )->QuatDot[0];
    ++j;
    aStateRate[j] = ( *it )->QuatDot[1];
    ++j;
    aStateRate[j] = ( *it )->QuatDot[2];
    ++j;
    aStateRate[j] = ( *it )->QuatDot[3];
    ++j;
  };

  for( std::vector< shared_ptr< gen_coord< double > > >::const_iterator it = mCoords.begin(); it < mCoords.end();
       ++it, ++j )
    aStateRate[j] = ( *it )->f;

  for( std::vector< shared_ptr< frame_2D< double > > >::const_iterator it = mFrames2D.begin(); it < mFrames2D.end();
       ++it ) {
    aStateRate[j] = ( *it )->Force[0];
    ++j;
    aStateRate[j] = ( *it )->Force[1];
    ++j;
    aStateRate[j] = ( *it )->Torque;
    ++j;
  };

  for( std::vector< shared_ptr< frame_3D< double > > >::const_iterator it = mFrames3D.begin(); it < mFrames3D.end();
       ++it ) {
    aStateRate[j] = ( *it )->Force[0];
    ++j;
    aStateRate[j] = ( *it )->Force[1];
    ++j;
    aStateRate[j] = ( *it )->Force[2];
    ++j;
    aStateRate[j] = ( *it )->Torque[0];
    ++j;
    aStateRate[j] = ( *it )->Torque[1];
    ++j;
    aStateRate[j] = ( *it )->Torque[2];
    ++j;
  };

  mat< double, mat_structure::symmetric > Msys( getJointAccelerationsCount() );
  getMassMatrix( Msys );
  try {
    mat_vect_adaptor< vect_n< double > > acc_as_mat( aStateRate, getJointAccelerationsCount(), 1,
                                                     getJointPositionsCount() );
    linsolve_Cholesky( Msys, acc_as_mat );
  } catch( singularity_error& e ) {
    RK_UNUSED( e );
    std::stringstream ss;
    ss << "Mass matrix is singular in the manipulator model '" << getName() << "' at time " << aTime << " seconds.";
    throw singularity_error( ss.str() );
  };
};

vect_n< double > manipulator_dynamics_model::getDependentStates() const {
  vect_n< double > result( getDependentStatesCount() );

  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .getDependentPositions( &result[0] );
  manip_kin_mdl_joint_io( shared_ptr< const direct_kinematics_model >( this, null_deleter() ) )
    .getDependentVelocities( &result[getDependentPositionsCount()] );

  return result;
};


void manipulator_dynamics_model::getMassMatrix( mat< double, mat_structure::symmetric >& M ) {
  mMassCalc.getMassMatrix( M );
};

void manipulator_dynamics_model::getMassMatrixAndDerivative( mat< double, mat_structure::symmetric >& M,
                                                             mat< double, mat_structure::square >& M_dot ) {
  mMassCalc.getMassMatrixAndDerivative( M, M_dot );
};

void manipulator_dynamics_model::get_TMT_TdMT( mat< double, mat_structure::rectangular >& Tcm,
                                               mat< double, mat_structure::symmetric >& Mcm,
                                               mat< double, mat_structure::rectangular >& Tcm_dot ) {
  mMassCalc.get_TMT_TdMT( Tcm, Mcm, Tcm_dot );
};
};
};
