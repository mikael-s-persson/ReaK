/**
 * \file prox_sphere_sphere.hpp
 *
 * This library declares a class for proximity queries between spheres.
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

#ifndef REAK_PROX_SPHERE_SPHERE_HPP
#define REAK_PROX_SPHERE_SPHERE_HPP

#include "ReaK/geometry/proximity/proximity_finder_3D.hpp"

#include "ReaK/geometry/shapes/sphere.hpp"

namespace ReaK::geom {

proximity_record_3D compute_proximity(const sphere& aSphere1,
                                      const shape_3D_precompute_pack& aPack1,
                                      const sphere& aSphere2,
                                      const shape_3D_precompute_pack& aPack2);

/**
 * This class is for proximity queries between spheres.
 */
class prox_sphere_sphere : public proximity_finder_3D {
 protected:
  const sphere* mSphere1;
  const sphere* mSphere2;

 public:
  /** This function performs the proximity query on its associated shapes. */
  proximity_record_3D computeProximity(
      const shape_3D_precompute_pack& aPack1,
      const shape_3D_precompute_pack& aPack2) override;

  /**
   * Default constructor.
   * \param aSphere1 The first sphere involved in the proximity query.
   * \param aSphere2 The second sphere involved in the proximity query.
   */
  explicit prox_sphere_sphere(const sphere* aSphere1 = nullptr,
                              const sphere* aSphere2 = nullptr);

  /** Destructor. */
  ~prox_sphere_sphere() override = default;
};

}  // namespace ReaK::geom

#endif
