/**
 * \file jacobian_joint_map.hpp
 *
 * This library declares types for jacobian mappings of generalized coordinates. These jacobian frames are
 * associated to the motion in a generalized coordinate by the joints that have these generalized coordinates
 * as input.
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

#ifndef REAK_JACOBIAN_JOINT_MAP_HPP
#define REAK_JACOBIAN_JOINT_MAP_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <ReaK/math/kinetostatics/motion_jacobians.hpp>
#include <map>

namespace ReaK {

namespace kte {


/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< shared_ptr< gen_coord< double > >, shared_ptr< jacobian_gen_gen< double > > > jacobian_joint_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< shared_ptr< gen_coord< double > >, shared_ptr< jacobian_gen_2D< double > > > jacobian_joint_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< shared_ptr< gen_coord< double > >, shared_ptr< jacobian_gen_3D< double > > > jacobian_joint_map_3D;


/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< shared_ptr< frame_2D< double > >, shared_ptr< jacobian_2D_gen< double > > > jacobian_joint2D_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< shared_ptr< frame_2D< double > >, shared_ptr< jacobian_2D_2D< double > > > jacobian_joint2D_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< shared_ptr< frame_2D< double > >, shared_ptr< jacobian_2D_3D< double > > > jacobian_joint2D_map_3D;


/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< shared_ptr< frame_3D< double > >, shared_ptr< jacobian_3D_gen< double > > > jacobian_joint3D_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< shared_ptr< frame_3D< double > >, shared_ptr< jacobian_3D_2D< double > > > jacobian_joint3D_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< shared_ptr< frame_3D< double > >, shared_ptr< jacobian_3D_3D< double > > > jacobian_joint3D_map_3D;


class joint_dependent_gen_coord : public shared_object {
public:
  shared_ptr< gen_coord< double > > mFrame; ///< Holds the generalized coordinate.

  jacobian_joint_map_gen mUpStreamJoints;     ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_gen mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_gen mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The generalized coordinate.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  joint_dependent_gen_coord( const shared_ptr< gen_coord< double > >& aFrame = shared_ptr< gen_coord< double > >(),
                             const jacobian_joint_map_gen& aUpStreamJoints = jacobian_joint_map_gen(),
                             const jacobian_joint2D_map_gen& aUpStream2DJoints = jacobian_joint2D_map_gen(),
                             const jacobian_joint3D_map_gen& aUpStream3DJoints = jacobian_joint3D_map_gen() )
      : mFrame( aFrame ), mUpStreamJoints( aUpStreamJoints ), mUpStream2DJoints( aUpStream2DJoints ),
        mUpStream3DJoints( aUpStream3DJoints ){};

  /**
   * Default destructor.
   */
  virtual ~joint_dependent_gen_coord(){};

  joint_dependent_gen_coord& add_joint( const shared_ptr< gen_coord< double > >& aJointFrame,
                                        const shared_ptr< jacobian_gen_gen< double > >& aJointJacobian ) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_gen_coord& add_joint( const shared_ptr< frame_2D< double > >& aJointFrame,
                                        const shared_ptr< jacobian_2D_gen< double > >& aJointJacobian ) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_gen_coord& add_joint( const shared_ptr< frame_3D< double > >& aJointFrame,
                                        const shared_ptr< jacobian_3D_gen< double > >& aJointJacobian ) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_gen_coord& remove_joint( const shared_ptr< gen_coord< double > >& aJointFrame ) {
    mUpStreamJoints.erase( mUpStreamJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_gen_coord& remove_joint( const shared_ptr< frame_2D< double > >& aJointFrame ) {
    mUpStream2DJoints.erase( mUpStream2DJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_gen_coord& remove_joint( const shared_ptr< frame_3D< double > >& aJointFrame ) {
    mUpStream3DJoints.erase( mUpStream3DJoints.find( aJointFrame ) );
    return *this;
  };


  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_SAVE_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_SAVE_WITH_NAME( mUpStream3DJoints );
  };

  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_LOAD_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_LOAD_WITH_NAME( mUpStream3DJoints );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( joint_dependent_gen_coord, 0xC2000002, 1, "joint_dependent_gen_coord", shared_object )
};


class joint_dependent_frame_2D : public shared_object {
public:
  shared_ptr< frame_2D< double > > mFrame; ///< Holds the generalized coordinate.

  jacobian_joint_map_2D mUpStreamJoints;     ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_2D mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_2D mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The generalized coordinate.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  joint_dependent_frame_2D( const shared_ptr< frame_2D< double > >& aFrame = shared_ptr< frame_2D< double > >(),
                            const jacobian_joint_map_2D& aUpStreamJoints = jacobian_joint_map_2D(),
                            const jacobian_joint2D_map_2D& aUpStream2DJoints = jacobian_joint2D_map_2D(),
                            const jacobian_joint3D_map_2D& aUpStream3DJoints = jacobian_joint3D_map_2D() )
      : mFrame( aFrame ), mUpStreamJoints( aUpStreamJoints ), mUpStream2DJoints( aUpStream2DJoints ),
        mUpStream3DJoints( aUpStream3DJoints ){};

  /**
   * Default destructor.
   */
  virtual ~joint_dependent_frame_2D(){};

  joint_dependent_frame_2D& add_joint( const shared_ptr< gen_coord< double > >& aJointFrame,
                                       const shared_ptr< jacobian_gen_2D< double > >& aJointJacobian ) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_2D& add_joint( const shared_ptr< frame_2D< double > >& aJointFrame,
                                       const shared_ptr< jacobian_2D_2D< double > >& aJointJacobian ) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_2D& add_joint( const shared_ptr< frame_3D< double > >& aJointFrame,
                                       const shared_ptr< jacobian_3D_2D< double > >& aJointJacobian ) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_2D& remove_joint( const shared_ptr< gen_coord< double > >& aJointFrame ) {
    mUpStreamJoints.erase( mUpStreamJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_frame_2D& remove_joint( const shared_ptr< frame_2D< double > >& aJointFrame ) {
    mUpStream2DJoints.erase( mUpStream2DJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_frame_2D& remove_joint( const shared_ptr< frame_3D< double > >& aJointFrame ) {
    mUpStream3DJoints.erase( mUpStream3DJoints.find( aJointFrame ) );
    return *this;
  };


  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_SAVE_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_SAVE_WITH_NAME( mUpStream3DJoints );
  };

  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_LOAD_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_LOAD_WITH_NAME( mUpStream3DJoints );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( joint_dependent_frame_2D, 0xC2000003, 1, "joint_dependent_frame_2D", shared_object )
};


class joint_dependent_frame_3D : public shared_object {
public:
  shared_ptr< frame_3D< double > > mFrame; ///< Holds the generalized coordinate.

  jacobian_joint_map_3D mUpStreamJoints;     ///< Holds the jacobian mappings of up-stream joints.
  jacobian_joint2D_map_3D mUpStream2DJoints; ///< Holds the jacobian mappings of up-stream 2D joints.
  jacobian_joint3D_map_3D mUpStream3DJoints; ///< Holds the jacobian mappings of up-stream 3D joints.

  /**
   * Parametrized and Default constructor.
   * \param aFrame The 3D frame.
   * \param aUpStreamJoints The jacobian mappings of up-stream joints.
   * \param aUpStream2DJoints The jacobian mappings of up-stream 2D joints.
   * \param aUpStream3DJoints The jacobian mappings of up-stream 3D joints.
   */
  joint_dependent_frame_3D( const shared_ptr< frame_3D< double > >& aFrame = shared_ptr< frame_3D< double > >(),
                            const jacobian_joint_map_3D& aUpStreamJoints = jacobian_joint_map_3D(),
                            const jacobian_joint2D_map_3D& aUpStream2DJoints = jacobian_joint2D_map_3D(),
                            const jacobian_joint3D_map_3D& aUpStream3DJoints = jacobian_joint3D_map_3D() )
      : mFrame( aFrame ), mUpStreamJoints( aUpStreamJoints ), mUpStream2DJoints( aUpStream2DJoints ),
        mUpStream3DJoints( aUpStream3DJoints ){};

  /**
   * Default destructor.
   */
  virtual ~joint_dependent_frame_3D(){};

  joint_dependent_frame_3D& add_joint( const shared_ptr< gen_coord< double > >& aJointFrame,
                                       const shared_ptr< jacobian_gen_3D< double > >& aJointJacobian ) {
    mUpStreamJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_3D& add_joint( const shared_ptr< frame_2D< double > >& aJointFrame,
                                       const shared_ptr< jacobian_2D_3D< double > >& aJointJacobian ) {
    mUpStream2DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_3D& add_joint( const shared_ptr< frame_3D< double > >& aJointFrame,
                                       const shared_ptr< jacobian_3D_3D< double > >& aJointJacobian ) {
    mUpStream3DJoints[aJointFrame] = aJointJacobian;
    return *this;
  };

  joint_dependent_frame_3D& remove_joint( const shared_ptr< gen_coord< double > >& aJointFrame ) {
    mUpStreamJoints.erase( mUpStreamJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_frame_3D& remove_joint( const shared_ptr< frame_2D< double > >& aJointFrame ) {
    mUpStream2DJoints.erase( mUpStream2DJoints.find( aJointFrame ) );
    return *this;
  };

  joint_dependent_frame_3D& remove_joint( const shared_ptr< frame_3D< double > >& aJointFrame ) {
    mUpStream3DJoints.erase( mUpStream3DJoints.find( aJointFrame ) );
    return *this;
  };


  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const {
    A& RK_SERIAL_SAVE_WITH_NAME( mFrame ) & RK_SERIAL_SAVE_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_SAVE_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_SAVE_WITH_NAME( mUpStream3DJoints );
  };

  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int ) {
    A& RK_SERIAL_LOAD_WITH_NAME( mFrame ) & RK_SERIAL_LOAD_WITH_NAME( mUpStreamJoints )
      & RK_SERIAL_LOAD_WITH_NAME( mUpStream2DJoints ) & RK_SERIAL_LOAD_WITH_NAME( mUpStream3DJoints );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( joint_dependent_frame_3D, 0xC2000004, 1, "joint_dependent_frame_3D", shared_object )
};
};
};

#endif // JACOBIAN_JOINT_MAP_HPP
