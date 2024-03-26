/**
 * \file calibrate_frames_3D.h
 *
 * This library declares functions to perform calibrations of relative poses and rotations,
 * between pose or point pairs expressed in different coordinate frames.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_MATH_KINETOSTATICS_CALIBRATE_FRAMES_3D_H_
#define REAK_MATH_KINETOSTATICS_CALIBRATE_FRAMES_3D_H_

#include "ReaK/math/kinetostatics/pose_3D.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_svd_method.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include <vector>

namespace ReaK {

template <typename T>
pose_3D<T> get_relative_pose_pointcloud(
    const std::vector<std::pair<vect<T, 3>, vect<T, 3>>>& aPointCloud) {

  vect<T, 3> centroid_1;
  vect<T, 3> centroid_2;
  for (const auto& [p1, p2] : aPointCloud) {
    centroid_1 += p1;
    centroid_2 += p2;
  }
  centroid_1 *= 1.0 / aPointCloud.size();
  centroid_2 *= 1.0 / aPointCloud.size();

  mat<T, mat_structure::rectangular> X(aPointCloud.size(), 3);
  mat<T, mat_structure::rectangular> B(aPointCloud.size(), 3);

  int i = 0;
  for (const auto& [p1, p2] : aPointCloud) {
    slice(X)(i, range(0, 3)) = p1 - centroid_1;
    slice(B)(i, range(0, 3)) = p2 - centroid_2;
    ++i;
  }

  // Solve the Orthogonal Procrustes problem: (Kabsch's algorithm)
  mat<T, mat_structure::square> C(transpose_view(X) * B);

  mat<T, mat_structure::square> U(3);
  mat<T, mat_structure::diagonal> E(3);
  mat<T, mat_structure::square> V(3);

  decompose_SVD(C, U, E, V);

  rot_mat_3D<T> Rt(U * transpose_view(V));

  return pose_3D<T>(std::weak_ptr<pose_3D<T>>(), centroid_1 - Rt * centroid_2,
                    quaternion<T>(Rt));
}

template <typename T>
rot_mat_3D<T> get_relative_rotation_vectcloud(
    const std::vector<std::pair<vect<T, 3>, vect<T, 3>>>& aVectCloud) {
  mat<T, mat_structure::rectangular> X(aVectCloud.size(), 3);
  mat<T, mat_structure::rectangular> B(aVectCloud.size(), 3);
  for (int i = 0; i < aVectCloud.size(); ++i) {
    for (const auto& [v1, v2] : aVectCloud) {
      slice(X)(i, range(0, 3)) = v1;
      slice(B)(i, range(0, 3)) = v2;
    }
  }

  // Solve the Orthogonal Procrustes problem: (Kabsch's algorithm)
  mat<T, mat_structure::square> C(transpose_view(X) * B);

  mat<T, mat_structure::square> U(3);
  mat<T, mat_structure::diagonal> E(3);
  mat<T, mat_structure::square> V(3);

  decompose_SVD(C, U, E, V);

  return rot_mat_3D<T>(U * transpose_view(V));
}

}  // namespace ReaK

#endif  // REAK_MATH_KINETOSTATICS_CALIBRATE_FRAMES_3D_H_
