
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

#include <ReaK/mbd/kte/virtual_kte_interface.hpp>

namespace ReaK {

namespace kte {


void virtual_kte_interface_gen::doMotion( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  ( *mEnd ) = ( *mBase );

  if( ( aFlag == store_kinematics ) && ( aStorage ) ) {
    if( !( aStorage->gen_coord_mapping[mBase] ) )
      aStorage->gen_coord_mapping[mBase]
        = shared_ptr< gen_coord< double > >( new gen_coord< double >( ( *mBase ) ), scoped_deleter() );
    else
      ( *( aStorage->gen_coord_mapping[mBase] ) ) = ( *mBase );
    if( !( aStorage->gen_coord_mapping[mEnd] ) )
      aStorage->gen_coord_mapping[mEnd]
        = shared_ptr< gen_coord< double > >( new gen_coord< double >( ( *mEnd ) ), scoped_deleter() );
    else
      ( *( aStorage->gen_coord_mapping[mEnd] ) ) = ( *mEnd );
  };
};

void virtual_kte_interface_gen::doForce( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  mBase->f -= mEnd->f;

  if( ( aFlag == store_dynamics ) && ( aStorage ) ) {
    if( aStorage->gen_coord_mapping[mEnd] ) {
      aStorage->gen_coord_mapping[mEnd]->f = mEnd->f;
    };
  };
};

void virtual_kte_interface_gen::clearForce() {
  if( mEnd ) {
    mEnd->f = 0.0;
  };
  if( mBase ) {
    mBase->f = 0.0;
  };
};


void virtual_kte_interface_2D::doMotion( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  ( *mEnd ) = ( *mBase );

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
  };
};

void virtual_kte_interface_2D::doForce( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  mBase->Force -= mEnd->Force;
  mBase->Torque -= mEnd->Torque;

  if( ( aFlag == store_dynamics ) && ( aStorage ) ) {
    if( aStorage->frame_2D_mapping[mEnd] ) {
      aStorage->frame_2D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_2D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};

void virtual_kte_interface_2D::clearForce() {
  if( mEnd ) {
    mEnd->Force = vect< double, 2 >();
    mEnd->Torque = 0.0;
  };
  if( mBase ) {
    mBase->Force = vect< double, 2 >();
    mBase->Torque = 0.0;
  };
};


void virtual_kte_interface_3D::doMotion( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  ( *mEnd ) = ( *mBase );

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
  };
};

void virtual_kte_interface_3D::doForce( kte_pass_flag aFlag, const shared_ptr< frame_storage >& aStorage ) {
  if( ( !mEnd ) || ( !mBase ) )
    return;

  mBase->Force -= mEnd->Force;
  mBase->Torque -= mEnd->Torque;

  if( ( aFlag == store_dynamics ) && ( aStorage ) ) {
    if( aStorage->frame_3D_mapping[mEnd] ) {
      aStorage->frame_3D_mapping[mEnd]->Force = mEnd->Force;
      aStorage->frame_3D_mapping[mEnd]->Torque = mEnd->Torque;
    };
  };
};


void virtual_kte_interface_3D::clearForce() {
  if( mEnd ) {
    mEnd->Force = vect< double, 3 >();
    mEnd->Torque = vect< double, 3 >();
  };
  if( mBase ) {
    mBase->Force = vect< double, 3 >();
    mBase->Torque = vect< double, 3 >();
  };
};
};
};
