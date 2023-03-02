/**
 * \file .hpp
 *
 * This library declares
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

#ifndef REAK_PROX_CYLINDER_BOX_HPP
#define REAK_PROX_CYLINDER_BOX_HPP

#include "proximity_finder_3D.hpp"

#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/cylinder.hpp>

namespace ReaK::geom {

proximity_record_3D compute_proximity(const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack1,
                                      const box& aBox,
                                      const shape_3D_precompute_pack& aPack2);

proximity_record_3D compute_proximity(const box& aBox,
                                      const shape_3D_precompute_pack& aPack1,
                                      const cylinder& aCylinder,
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a cylinder and a box.
 */
class prox_cylinder_box : public proximity_finder_3D {
 protected:
  const cylinder* mCylinder;
  const box* mBox;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aCCylinder The capped cylinder involved in the proximity query.
   * \param aBox The box involved in the proximity query.
   */
  explicit prox_cylinder_box(const cylinder* aCylinder = nullptr,
                             const box* aBox = nullptr);

  /** Destructor. */
  ~prox_cylinder_box() override = default;
};

}  // namespace ReaK::geom

#endif
