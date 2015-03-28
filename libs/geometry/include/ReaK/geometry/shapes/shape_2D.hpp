/**
 * \file shape_2D.hpp
 *
 * This library declares a base-class for 2D shapes (collidable primitives).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2012
 */

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

#ifndef REAK_SHAPE_2D_HPP
#define REAK_SHAPE_2D_HPP

#include "geometry_2D.hpp"
#include <ReaK/math/kinetostatics/pose_2D.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


class shape_2D; // forward-decl

class shape_2D_precompute_pack {
public:
  const shape_2D* parent;
  pose_2D< double > global_pose;

  shape_2D_precompute_pack() : parent(), global_pose(){};
};


/** This class is a base-class for all 2D shapes (collidable primitives). */
class shape_2D : public geometry_2D {
protected:
public:
  /**
   * This function returns the maximum radius of the shape (radius of the circle that bounds the shape).
   * \return The maximum radius of the shape.
   */
  virtual double getBoundingRadius() const = 0;

  virtual shape_2D_precompute_pack createPrecomputePack() const;

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   */
  shape_2D( const std::string& aName = "",
            const shared_ptr< pose_2D< double > >& aAnchor = shared_ptr< pose_2D< double > >(),
            const pose_2D< double >& aPose = pose_2D< double >() );

  /**
   * Default destructor.
   */
  virtual ~shape_2D(){};


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const;

  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int );

  RK_RTTI_MAKE_ABSTRACT_1BASE( shape_2D, 0xC3100008, 1, "shape_2D", geometry_2D )
};
};
};

#endif
