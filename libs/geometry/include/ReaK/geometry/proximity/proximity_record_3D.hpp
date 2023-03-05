/**
 * \file proximity_record_3D.hpp
 *
 * This library declares a class to record the results of a proximity query between 3D shapes.
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

#ifndef REAK_PROXIMITY_RECORD_3D_HPP
#define REAK_PROXIMITY_RECORD_3D_HPP

#include "ReaK/core/base/shared_object.hpp"
#include "ReaK/math/lin_alg/vect_alg.hpp"

#include <limits>

namespace ReaK::geom {

/**
 * This class stores the data which results from a proximity query with 3D shapes.
 */
class proximity_record_3D : public shared_object {
 public:
  /** Holds the closest point (in global coordinates) on the first shape involved in the proximity query. */
  vect<double, 3> mPoint1;
  /** Holds the closest point (in global coordinates) on the second shape involved in the proximity query. */
  vect<double, 3> mPoint2;

  /** Holds the distance between the shapes, a negative value denotes penetration. */
  double mDistance;

  /** Default constructor. */
  proximity_record_3D() : mDistance(std::numeric_limits<double>::infinity()) {}

  /** Destructor. */
  ~proximity_record_3D() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override;

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(proximity_record_3D, 0xC3200002, 1,
                              "proximity_record_3D", shared_object)
};

}  // namespace ReaK::geom

#endif
