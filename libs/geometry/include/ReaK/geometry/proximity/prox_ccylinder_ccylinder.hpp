/**
 * \file prox_ccylinder_ccylinder.hpp
 *
 * This library declares a class for proximity queries between a capped cylinder and a capped cylinder.
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

#ifndef REAK_PROX_CCYLINDER_CCYLINDER_HPP
#define REAK_PROX_CCYLINDER_CCYLINDER_HPP

#include "ReaK/geometry/proximity/proximity_finder_3D.hpp"

#include "ReaK/geometry/shapes/capped_cylinder.hpp"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const capped_cylinder& aCCylinder1,
                                      const shape_3D_precompute_pack& aPack1,
                                      const capped_cylinder& aCCylinder2,
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between a capped cylinder and a capped cylinder.
 */
class prox_ccylinder_ccylinder : public proximity_finder_3D {
 protected:
  const capped_cylinder* mCCylinder1;
  const capped_cylinder* mCCylinder2;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aCCylinder1 The capped cylinder involved in the proximity query.
   * \param aCCylinder2 The capped cylinder involved in the proximity query.
   */
  explicit prox_ccylinder_ccylinder(
      const capped_cylinder* aCCylinder1 = nullptr,
      const capped_cylinder* aCCylinder2 = nullptr);

  /** Destructor. */
  ~prox_ccylinder_ccylinder() override = default;
};

}  // namespace ReaK::geom

#endif
