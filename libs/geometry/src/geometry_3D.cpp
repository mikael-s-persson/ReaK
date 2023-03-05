
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

#include "ReaK/geometry/shapes/geometry_3D.hpp"
#include <utility>

namespace ReaK::geom {

void geometry_3D::setAnchor(const std::shared_ptr<pose_3D<double>>& aAnchor) {
  mAnchor = aAnchor;
  mPose.Parent = mAnchor;
}

void geometry_3D::setPose(const pose_3D<double>& aPose) {
  mPose = aPose;
  mPose.Parent = mAnchor;
}

geometry_3D::geometry_3D(const std::string& aName,
                         std::shared_ptr<pose_3D<double>> aAnchor,
                         const pose_3D<double>& aPose)
    : mAnchor(std::move(aAnchor)), mPose(aPose) {
  mPose.Parent = mAnchor;
  this->setName(aName);
}

void geometry_3D::save(ReaK::serialization::oarchive& A,
                       unsigned int /*unused*/) const {
  named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(mAnchor) & RK_SERIAL_SAVE_WITH_NAME(mPose);
}

void geometry_3D::load(ReaK::serialization::iarchive& A,
                       unsigned int /*unused*/) {
  named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(mAnchor) & RK_SERIAL_LOAD_WITH_NAME(mPose);
}

}  // namespace ReaK::geom
