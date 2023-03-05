
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

#include "ReaK/geometry/proximity/proximity_record_2D.hpp"

namespace ReaK::geom {

void proximity_record_2D::save(ReaK::serialization::oarchive& A,
                               unsigned int /*Version*/) const {
  shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(mPoint1) & RK_SERIAL_SAVE_WITH_NAME(mPoint2) &
      RK_SERIAL_SAVE_WITH_NAME(mDistance);
}

void proximity_record_2D::load(ReaK::serialization::iarchive& A,
                               unsigned int /*Version*/) {
  shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(mPoint1) & RK_SERIAL_LOAD_WITH_NAME(mPoint2) &
      RK_SERIAL_LOAD_WITH_NAME(mDistance);
}

}  // namespace ReaK::geom
