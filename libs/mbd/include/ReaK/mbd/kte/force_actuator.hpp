/**
 * \file force_actuator.hpp
 *
 * This library declares KTE models for a force actuator, that is, an actuator that applies some
 * given force. This model is needed to apply the reaction force that
 * the actuator will imbue on its base frame. This is a base class for other force actuator models,
 * but its default behaviour is to reverse the forces on the joints and apply the opposite reactions,
 * this can be used directly in a virtual model control scheme.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_FORCE_ACTUATOR_HPP
#define REAK_FORCE_ACTUATOR_HPP

#include <ReaK/core/base/defs.hpp>

#include "kte_map.hpp"

#include "reacting_kte.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements a force actuator acting on a joint with generalized coordinate. This base class'
 * default behaviour is to take the sum of forces on the joint and apply them as a reaction on the joint,
 * as if, the control would cancel all the forces on the joint, e.g., like in a virtual model control scheme.
 */
class force_actuator_gen : public kte_map {
protected:
  shared_ptr< gen_coord< double > > mFrame; ///< Holds the generalized coordinate on which the actuator acts.
  shared_ptr< reacting_kte_gen > mJoint;    ///< Holds the joint which will react to the actuator's force.

public:
  /**
   * Sets the frame on which the actuator acts.
   * \param aPtr The new frame on which the actuator acts.
   */
  void setFrame( const shared_ptr< gen_coord< double > >& aPtr ) { mFrame = aPtr; };
  /**
   * Returns the frame on which the actuator acts.
   * \return The frame on which the actuator acts.
   */
  shared_ptr< gen_coord< double > > Frame() const { return mFrame; };

  /**
   * Sets the joint which will react to the actuator's force.
   * \param aPtr The new joint which will react to the actuator's force.
   */
  void setJoint( const shared_ptr< reacting_kte_gen >& aPtr ) { mJoint = aPtr; };
  /**
   * Returns the joint which will react to the actuator's force.
   * \return The joint which will react to the actuator's force.
   */
  shared_ptr< reacting_kte_gen > Joint() const { return mJoint; };

  /**
  * Default constructor.
  */
  force_actuator_gen( const std::string& aName = "" ) : kte_map( aName ), mFrame(), mJoint(){};

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aFrame the generalized coordinate on which the actuator acts.
   * \param aJoint the joint which will react to the actuator's force.
   */
  force_actuator_gen( const std::string& aName, const shared_ptr< gen_coord< double > >& aFrame,
                      const shared_ptr< reacting_kte_gen >& aJoint )
      : kte_map( aName ), mFrame( aFrame ), mJoint( aJoint ){};

  /**
   * Default destructor.
   */
  virtual ~force_actuator_gen(){};

  virtual void doMotion( kte_pass_flag aFlag = nothing,
                         const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void doForce( kte_pass_flag aFlag = nothing,
                        const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void clearForce();

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    kte_map::save( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mJoint );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    kte_map::load( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mJoint );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( force_actuator_gen, 0xC2100019, 1, "force_actuator_gen", kte_map )
};


/**
 * This class implements a force actuator acting on a joint with a 2D frame. This base class'
 * default behaviour is to take the sum of forces on the joint and apply them as a reaction on the joint,
 * as if, the control would cancel all the forces on the joint, e.g., like in a virtual model control scheme.
 */
class force_actuator_2D : public kte_map {
protected:
  shared_ptr< frame_2D< double > > mFrame; ///< Holds the 2D frame on which the actuator acts.
  shared_ptr< reacting_kte_2D > mJoint;    ///< Holds the joint which will react to the actuator's force and torque.

public:
  /**
   * Sets the frame on which the actuator acts.
   * \param aPtr The new frame on which the actuator acts.
   */
  void setFrame( const shared_ptr< frame_2D< double > >& aPtr ) { mFrame = aPtr; };
  /**
   * Returns the frame on which the actuator acts.
   * \return The frame on which the actuator acts.
   */
  shared_ptr< frame_2D< double > > Frame() const { return mFrame; };

  /**
   * Sets the joint which will react to the actuator's force.
   * \param aPtr The new joint which will react to the actuator's force.
   */
  void setJoint( const shared_ptr< reacting_kte_2D >& aPtr ) { mJoint = aPtr; };
  /**
   * Returns the joint which will react to the actuator's force.
   * \return The joint which will react to the actuator's force.
   */
  shared_ptr< reacting_kte_2D > Joint() const { return mJoint; };

  /**
  * Default constructor.
  */
  force_actuator_2D( const std::string& aName = "" ) : kte_map( aName ), mFrame(), mJoint(){};

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aFrame the 2D frame on which the actuator acts.
   * \param aJoint the joint which will react to the actuator's force and torque.
   */
  force_actuator_2D( const std::string& aName, const shared_ptr< frame_2D< double > >& aFrame,
                     const shared_ptr< reacting_kte_2D >& aJoint )
      : kte_map( aName ), mFrame( aFrame ), mJoint( aJoint ){};

  /**
   * Default destructor.
   */
  virtual ~force_actuator_2D(){};

  virtual void doMotion( kte_pass_flag aFlag = nothing,
                         const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void doForce( kte_pass_flag aFlag = nothing,
                        const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void clearForce();

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    kte_map::save( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mJoint );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    kte_map::load( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mJoint );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( force_actuator_2D, 0xC210001A, 1, "force_actuator_2D", kte_map )
};

/**
 * This class implements a force actuator acting on a joint with a 3D frame. This base class'
 * default behaviour is to take the sum of forces on the joint and apply them as a reaction on the joint,
 * as if, the control would cancel all the forces on the joint, e.g., like in a virtual model control scheme.
 */
class force_actuator_3D : public kte_map {
protected:
  shared_ptr< frame_3D< double > > mFrame; ///< Holds the 3D frame on which the actuator acts.
  shared_ptr< reacting_kte_3D > mJoint;    ///< Holds the joint which will react to the actuator's force and torque.

public:
  /**
   * Sets the frame on which the actuator acts.
   * \param aPtr The new frame on which the actuator acts.
   */
  void setFrame( const shared_ptr< frame_3D< double > >& aPtr ) { mFrame = aPtr; };
  /**
   * Returns the frame on which the actuator acts.
   * \return The frame on which the actuator acts.
   */
  shared_ptr< frame_3D< double > > Frame() const { return mFrame; };

  /**
   * Sets the joint which will react to the actuator's force.
   * \param aPtr The new joint which will react to the actuator's force.
   */
  void setJoint( const shared_ptr< reacting_kte_3D >& aPtr ) { mJoint = aPtr; };
  /**
   * Returns the joint which will react to the actuator's force.
   * \return The joint which will react to the actuator's force.
   */
  shared_ptr< reacting_kte_3D > Joint() const { return mJoint; };

  /**
  * Default constructor.
  */
  force_actuator_3D( const std::string& aName = "" ) : kte_map( aName ), mFrame(), mJoint(){};

  /**
   * Parametrized constructor.
   * \param aName the name of the KTE model.
   * \param aFrame the 3D frame on which the actuator acts.
   * \param aJoint the joint which will react to the actuator's force and torque.
   */
  force_actuator_3D( const std::string& aName, const shared_ptr< frame_3D< double > >& aFrame,
                     const shared_ptr< reacting_kte_3D >& aJoint )
      : kte_map( aName ), mFrame( aFrame ), mJoint( aJoint ){};

  /**
   * Default destructor.
   */
  virtual ~force_actuator_3D(){};

  virtual void doMotion( kte_pass_flag aFlag = nothing,
                         const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void doForce( kte_pass_flag aFlag = nothing,
                        const shared_ptr< frame_storage >& aStorage = shared_ptr< frame_storage >() );

  virtual void clearForce();

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    kte_map::save( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mJoint );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    kte_map::load( A, kte_map::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mJoint );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( force_actuator_3D, 0xC210001B, 1, "force_actuator_3D", kte_map )
};
};
};

#endif
