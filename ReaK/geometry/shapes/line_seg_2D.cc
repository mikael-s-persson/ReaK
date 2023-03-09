
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

#include "ReaK/geometry/shapes/line_seg_2D.h"

namespace ReaK::geom {

line_seg_2D::line_seg_2D(const std::string& aName,
                         const std::shared_ptr<pose_2D<double>>& aAnchor,
                         const pose_2D<double>& aPose,
                         const vect<double, 2>& aStart,
                         const vect<double, 2>& aEnd)
    : geometry_2D(aName, aAnchor, aPose), mStart(aStart), mEnd(aEnd) {}

void line_seg_2D::save(ReaK::serialization::oarchive& A,
                       unsigned int /*unused*/) const {
  geometry_2D::save(A, geometry_2D::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(mStart) & RK_SERIAL_SAVE_WITH_NAME(mEnd);
}

void line_seg_2D::load(ReaK::serialization::iarchive& A,
                       unsigned int /*unused*/) {
  geometry_2D::load(A, geometry_2D::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(mStart) & RK_SERIAL_LOAD_WITH_NAME(mEnd);
}

}  // namespace ReaK::geom
