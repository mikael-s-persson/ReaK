/**
 * \file shape_3D.hpp
 *
 * This library declares a base-class for 3D shapes (collidable primitives).
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

#ifndef REAK_SHAPE_3D_HPP
#define REAK_SHAPE_3D_HPP

#include <ReaK/math/kinetostatics/pose_3D.hpp>
#include "geometry_3D.hpp"

namespace ReaK::geom {

class shape_3D;  // forward-decl

class shape_3D_precompute_pack {
 public:
  const shape_3D* parent{};
  pose_3D<double> global_pose;

  shape_3D_precompute_pack() = default;
};

/** This class is a base-class for all 3D shapes (collidable primitives). */
class shape_3D : public geometry_3D {
 protected:
 public:
  /**
   * This function returns the maximum radius of the shape (radius of the sphere that bounds the shape).
   * \return The maximum radius of the shape.
   */
  virtual double getBoundingRadius() const = 0;

  virtual shape_3D_precompute_pack createPrecomputePack() const;

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   */
  explicit shape_3D(const std::string& aName = "",
                    const std::shared_ptr<pose_3D<double>>& aAnchor =
                        std::shared_ptr<pose_3D<double>>(),
                    const pose_3D<double>& aPose = pose_3D<double>());

  /**
   * Default destructor.
   */
  ~shape_3D() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_ABSTRACT_1BASE(shape_3D, 0xC3100009, 1, "shape_3D", geometry_3D)
};

}  // namespace ReaK::geom

#endif
