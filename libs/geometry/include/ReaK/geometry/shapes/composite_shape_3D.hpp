/**
 * \file composite_shape_3D.hpp
 *
 * This library declares a class for a composite of 3D shapes.
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

#ifndef REAK_COMPOSITE_SHAPE_3D_HPP
#define REAK_COMPOSITE_SHAPE_3D_HPP

#include "shape_3D.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents a composite of 3D shapes. */
class composite_shape_3D : public shape_3D {
protected:
  std::vector< shared_ptr< shape_3D > > mShapes;

public:
  /**
   * This function returns the maximum radius of the shape (radius of the sphere that bounds the shape).
   * \return The maximum radius of the shape.
   */
  virtual double getBoundingRadius() const;

  /**
   * This function returns a const-reference to the vector of shapes.
   * \return A const-reference to the vector of shapes.
   */
  const std::vector< shared_ptr< shape_3D > >& Shapes() const { return mShapes; };
  /**
   * This function returns a reference to the vector of shapes.
   * \return A reference to the vector of shapes.
   */
  std::vector< shared_ptr< shape_3D > >& Shapes() { return mShapes; };

  /**
   * Default constructor.
   */
  composite_shape_3D( const std::string& aName = "" );

  /**
   * Default destructor.
   */
  virtual ~composite_shape_3D(){};


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( ReaK::serialization::oarchive& A, unsigned int ) const;

  virtual void RK_CALL load( ReaK::serialization::iarchive& A, unsigned int );

  RK_RTTI_MAKE_CONCRETE_1BASE( composite_shape_3D, 0xC310000B, 1, "composite_shape_3D", shape_3D )
};
};
};

#endif
