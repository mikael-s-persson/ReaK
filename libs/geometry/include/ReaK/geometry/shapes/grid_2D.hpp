/**
 * \file grid_2D.hpp
 *
 * This library declares a class for a 2D grid to be rendered.
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

#ifndef REAK_GRID_2D_HPP
#define REAK_GRID_2D_HPP

#include "ReaK/geometry/shapes/geometry_2D.hpp"

namespace ReaK::geom {

/** This class is a class for a 2D grid to be rendered. */
class grid_2D : public geometry_2D {
 protected:
  vect<double, 2> mDimensions;
  vect<std::size_t, 2> mSquareCounts;

 public:
  /**
   * This function returns the dimensions of the grid.
   * \return The dimensions of the grid.
   */
  const vect<double, 2>& getDimensions() const { return mDimensions; }
  /**
   * This function sets the new dimensions of the grid.
   * \param aDimensions The new dimensions of the grid.
   */
  void setDimensions(const vect<double, 2>& aDimensions) {
    mDimensions = aDimensions;
  }

  /**
   * This function returns the square-counts of the grid.
   * \return The square-counts.
   */
  const vect<std::size_t, 2>& getSquareCounts() const { return mSquareCounts; }
  /**
   * This function sets the new square-counts of the grid.
   * \param aSquareCounts The new square-counts.
   */
  void setSquareCounts(const vect<std::size_t, 2>& aSquareCounts) {
    mSquareCounts = aSquareCounts;
  }

  /**
   * Default constructor.
   * \param aName The name of the object.
   * \param aAnchor The anchor object for the geometry.
   * \param aPose The pose of the geometry (relative to the anchor).
   * \param aDimensions The dimensions.
   * \param aSquareCounts The square-counts.
   */
  explicit grid_2D(const std::string& aName = "",
                   const std::shared_ptr<pose_2D<double>>& aAnchor =
                       std::shared_ptr<pose_2D<double>>(),
                   const pose_2D<double>& aPose = pose_2D<double>(),
                   const vect<double, 2>& aDimensions = (vect<double, 2>()),
                   const vect<std::size_t, 2>& aSquareCounts =
                       (vect<std::size_t, 2>(10, 10)));

  /**
   * Default destructor.
   */
  ~grid_2D() override = default;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;

  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(grid_2D, 0xC3100006, 1, "grid_2D", geometry_2D)
};

}  // namespace ReaK::geom

#endif
