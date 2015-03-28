
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

#include <ReaK/mbd/kte/revolute_joint.hpp>

namespace ReaK {

namespace kte {


void revolute_joint_2D::doMotion( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  mEnd->Parent = mBase->Parent;

  mEnd->Position = mBase->Position;
  mEnd->Velocity = mBase->Velocity;
  mEnd->Acceleration = mBase->Acceleration;

  if( !mAngle ) {
    mEnd->Rotation = mBase->Rotation;
    mEnd->AngVelocity = mBase->AngVelocity;
    mEnd->AngAcceleration = mBase->AngAcceleration;
  } else {
    mEnd->Rotation = mBase->Rotation * rot_mat_2D< double >( mAngle->q );
    mEnd->AngVelocity = mBase->AngVelocity + mAngle->q_dot;
    mEnd->AngAcceleration = mBase->AngAcceleration + mAngle->q_ddot;

    if( mJacobian ) {
      mJacobian->Parent = mEnd;

      mJacobian->qd_vel = vect< double, 2 >();
      mJacobian->qd_avel = 1.0;
      mJacobian->qd_acc = vect< double, 2 >();
      mJacobian->qd_aacc = 0.0;
    };
  };

  if( ( aFlag == store_kinematics ) && ( aStorage ) ) {
    if( !( aStorage->frame_2D_mapping[mBase] ) )
      aStorage->frame_2D_mapping[mBase]
        = shared_ptr< frame_2D< double > >( new frame_2D< double >( ( *mBase ) ), scoped_deleter() );
    else
      ( *( aStorage->frame_2D_mapping[mBase] ) ) = ( *mBase );
    if( !( aStorage->frame_2D_mapping[mEnd] ) )
      aStorage->frame_2D_mapping[mEnd]
        = shared_ptr< frame_2D< double > >( new frame_2D< double >( ( *mEnd ) ), scoped_deleter() );
    else
      ( *( aStorage->frame_2D_mapping[mEnd] ) ) = ( *mEnd );
    if( mAngle ) {
      if( !( aStorage->gen_coord_mapping[mAngle] ) )
        aStorage->gen_coord_mapping[mAngle]
          = shared_ptr< gen_coord< double > >( new gen_coord< double >( ( *mAngle ) ), scoped_deleter() );
      else
        ( *( aStorage->gen_coord_mapping[mAngle] ) ) = ( *mAngle );
    };
  };
};

void revolute_joint_2D::doForce( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  if( !mAngle ) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    mBase->Force += ReaK::rot_mat_2D< double >( mAngle->q ) * mEnd->Force;
    mAngle->f += mEnd->Torque;
  };

  if( ( aFlag == store_dynamics ) && ( aStorage ) ) {
    if( aStorage->frame_2D_mapping[mEnd] ) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};

void revolute_joint_2D::clearForce() {
  if( mEnd ) {
    mEnd->Force = vect< double, 2 >();
    mEnd->Torque = 0.0;
  };
  if( mBase ) {
    mBase->Force = vect< double, 2 >();
    mBase->Torque = 0.0;
  };
  if( mAngle ) {
    mAngle->f = 0.0;
  };
};

void revolute_joint_2D::applyReactionForce( double aForce ) {
  if( mAngle )
    mBase->Torque -= aForce;
};


void revolute_joint_3D::doMotion( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  mEnd->Parent = mBase->Parent;

  mEnd->Position = mBase->Position;
  mEnd->Velocity = mBase->Velocity;
  mEnd->Acceleration = mBase->Acceleration;

  if( !mAngle ) {
    mEnd->Quat = mBase->Quat;
    mEnd->AngVelocity = mBase->AngVelocity;
    mEnd->AngAcceleration = mBase->AngAcceleration;
  } else {
    quaternion< double > tmp_quat( axis_angle< double >( mAngle->q, mAxis ).getQuaternion() );
    rot_mat_3D< double > R2( tmp_quat.getRotMat() );
    mEnd->Quat = mBase->Quat * tmp_quat;
    mEnd->AngVelocity = ( mBase->AngVelocity * R2 ) + mAngle->q_dot * mAxis;
    mEnd->AngAcceleration = ( mBase->AngAcceleration * R2 )
                            + ( ( mBase->AngVelocity * R2 ) % ( mAngle->q_dot * mAxis ) ) + mAngle->q_ddot * mAxis;

    if( mJacobian ) {
      mJacobian->Parent = mEnd;

      mJacobian->qd_vel = vect< double, 3 >();
      mJacobian->qd_avel = mAxis;
      mJacobian->qd_acc = vect< double, 3 >();
      mJacobian->qd_aacc = vect< double, 3 >();
    };
  };

  mEnd->UpdateQuatDot();

  if( ( aFlag == store_kinematics ) && ( aStorage ) ) {
    if( !( aStorage->frame_3D_mapping[mBase] ) )
      aStorage->frame_3D_mapping[mBase]
        = shared_ptr< frame_3D< double > >( new frame_3D< double >( ( *mBase ) ), scoped_deleter() );
    else
      ( *( aStorage->frame_3D_mapping[mBase] ) ) = ( *mBase );
    if( !( aStorage->frame_3D_mapping[mEnd] ) )
      aStorage->frame_3D_mapping[mEnd]
        = shared_ptr< frame_3D< double > >( new frame_3D< double >( ( *mEnd ) ), scoped_deleter() );
    else
      ( *( aStorage->frame_3D_mapping[mEnd] ) ) = ( *mEnd );
    if( mAngle ) {
      if( !( aStorage->gen_coord_mapping[mAngle] ) )
        aStorage->gen_coord_mapping[mAngle]
          = shared_ptr< gen_coord< double > >( new gen_coord< double >( ( *mAngle ) ), scoped_deleter() );
      else
        ( *( aStorage->gen_coord_mapping[mAngle] ) ) = ( *mAngle );
    };
  };
};

void revolute_joint_3D::doForce( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  if( !mAngle ) {
    mBase->Force += mEnd->Force;
    mBase->Torque += mEnd->Torque;
  } else {
    rot_mat_3D< double > R( axis_angle< double >( mAngle->q, mAxis ).getRotMat() );
    mBase->Force += R * mEnd->Force;
    mAngle->f += mEnd->Torque * mAxis;
    mBase->Torque += R * ( mEnd->Torque - ( mEnd->Torque * mAxis ) * mAxis );
  };

  if( ( aFlag == store_dynamics ) && ( aStorage ) ) {
    if( aStorage->frame_3D_mapping[mEnd] ) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};


void revolute_joint_3D::clearForce() {
  if( mEnd ) {
    mEnd->Force = vect< double, 3 >();
    mEnd->Torque = vect< double, 3 >();
  };
  if( mBase ) {
    mBase->Force = vect< double, 3 >();
    mBase->Torque = vect< double, 3 >();
  };
  if( mAngle ) {
    mAngle->f = 0.0;
  };
};


void revolute_joint_3D::applyReactionForce( double aForce ) {
  if( mAngle )
    mBase->Torque -= aForce * mAxis;
};
};
};
