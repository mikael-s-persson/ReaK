
/*
 *    Copyright 2011 Sven Mikael Persson
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

#include <ReaK/math/optimization/nl_interior_points_methods.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>

#include <iostream>

using namespace ReaK;

pose_3D<double> a1 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 0.0, 0.0),
    quaternion<double>(vect<double, 4>(0.8, 0.0, 0.6, 0.0)));
pose_3D<double> a2 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 3.0, 5.0),
    quaternion<double>(vect<double, 4>(0.8, -0.6, 0.0, 0.0)));
pose_3D<double> a3 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(10.0, -3.0, -2.0),
    quaternion<double>(vect<double, 4>(1.0, 0.0, 0.0, 0.0)));
pose_3D<double> a4 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(-3.0, -3.0, 6.0),
    quaternion<double>(vect<double, 4>(sqrt(3.0) / 3.0, 0.0, -sqrt(3.0) / 3.0,
                                       sqrt(3.0) / 3.0)));

pose_3D<double> a5 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 0.0, 0.0),
    quaternion<double>(vect<double, 4>(1.0, 0.0, 0.0, 0.0)));
pose_3D<double> a6 = pose_3D<double>(
    std::shared_ptr<pose_3D<double>>(), vect<double, 3>(0.5, 0.0, 4.0),
    quaternion<double>(vect<double, 4>(1.0, 0.0, 0.0, 0.0)));

auto cy1 = std::make_shared<geom::cylinder>(
    "cy1", std::shared_ptr<pose_3D<double>>(), a1, 5.0, 0.5);
auto cy2 = std::make_shared<geom::cylinder>(
    "cy2", std::shared_ptr<pose_3D<double>>(), a1, 10.0, 0.25);
auto cy3 = std::make_shared<geom::cylinder>(
    "cy3", std::shared_ptr<pose_3D<double>>(), a1, 1.0, 2.0);
auto cy4 = std::make_shared<geom::cylinder>(
    "cy4", std::shared_ptr<pose_3D<double>>(), a2, 5.0, 0.5);
auto cy5 = std::make_shared<geom::cylinder>(
    "cy5", std::shared_ptr<pose_3D<double>>(), a3, 5.0, 0.5);
auto cy6 = std::make_shared<geom::cylinder>(
    "cy6", std::shared_ptr<pose_3D<double>>(), a4, 5.0, 0.5);

auto cy7 = std::make_shared<geom::cylinder>(
    "cy7", std::shared_ptr<pose_3D<double>>(), a6, 5.0, 0.5);

auto bx1 =
    std::make_shared<geom::box>("bx1", std::shared_ptr<pose_3D<double>>(), a1,
                                vect<double, 3>(1.0, 2.0, 1.0));
auto bx2 =
    std::make_shared<geom::box>("bx2", std::shared_ptr<pose_3D<double>>(), a1,
                                vect<double, 3>(4.0, 1.0, 10.0));
auto bx3 =
    std::make_shared<geom::box>("bx3", std::shared_ptr<pose_3D<double>>(), a1,
                                vect<double, 3>(4.0, 4.0, 1.0));
auto bx4 =
    std::make_shared<geom::box>("bx4", std::shared_ptr<pose_3D<double>>(), a2,
                                vect<double, 3>(4.0, 2.0, 2.0));
auto bx5 =
    std::make_shared<geom::box>("bx5", std::shared_ptr<pose_3D<double>>(), a3,
                                vect<double, 3>(4.0, 2.0, 2.0));
auto bx6 =
    std::make_shared<geom::box>("bx6", std::shared_ptr<pose_3D<double>>(), a4,
                                vect<double, 3>(4.0, 2.0, 2.0));

auto bx7 =
    std::make_shared<geom::box>("bx7", std::shared_ptr<pose_3D<double>>(), a5,
                                vect<double, 3>(4.0, 2.0, 4.0));

template <typename BF1, typename BF2>
double get_brute_force_dist(BF1& bf1, BF2& bf2) {
  double min_dist_brute = std::numeric_limits<double>::infinity();
  for (double i = 0; i < 2.0 * M_PI; i += M_PI / 8.0) {
    for (double j = 0; j < M_PI; j += M_PI / 16.0) {
      vect<double, 3> u(cos(i) * sin(j), sin(i) * sin(j), cos(j));
      vect<double, 3> pt1 = bf1(u);
      for (double k = 0; k < 2.0 * M_PI; k += M_PI / 8.0) {
        for (double l = 0; l < M_PI; l += M_PI / 16.0) {
          vect<double, 3> v(cos(k) * sin(l), sin(k) * sin(l), cos(l));
          vect<double, 3> pt2 = bf2(v);
          double d = norm_2(pt2 - pt1);
          if (d < min_dist_brute)
            min_dist_brute = d;
        }
      }
    }
  }
  return min_dist_brute;
}

struct proximity_solver {
  std::shared_ptr<geom::shape_3D> mShape1;
  std::shared_ptr<geom::shape_3D> mShape2;

  proximity_solver(const std::shared_ptr<geom::shape_3D>& aShape1,
                   const std::shared_ptr<geom::shape_3D>& aShape2)
      : mShape1(aShape1), mShape2(aShape2) {}

  void operator()() {
    geom::proximity_record_3D result;
    vect<double, 3> v1, v2;
    using std::cos;
    using std::sin;

    if (mShape1->getObjectType() == geom::box::getStaticObjectType()) {
      std::shared_ptr<geom::box> bx1 =
          rtti::rk_dynamic_ptr_cast<geom::box>(mShape1);
      pose_3D<double> bx1_pose = bx1->getPose().getGlobalPose();
      if (mShape2->getObjectType() == geom::box::getStaticObjectType()) {
        // box-box case.
        auto bx2 = rtti::rk_dynamic_ptr_cast<geom::box>(mShape2);
        pose_3D<double> bx2_pose = bx2->getPose().getGlobalPose();

        std::cout << "Checking proximity between Box '" << bx1->getName()
                  << "' and Box '" << bx2->getName() << "'..." << std::endl;

        result =
            geom::findProximityByGJKEPA(geom::box_support_func(*bx1, bx1_pose),
                                        geom::box_support_func(*bx2, bx2_pose));

        v1 = mShape1->getPose().rotateToGlobal(
            mShape1->getPose().transformFromGlobal(result.mPoint1));
        v2 = mShape2->getPose().rotateToGlobal(
            mShape2->getPose().transformFromGlobal(result.mPoint2));

        geom::box_boundary_func bf1(*bx1, bx1_pose);
        geom::box_boundary_func bf2(*bx2, bx2_pose);
        v1 = bf1(v1);
        v2 = bf2(v2);

        std::cout << " which has brute-force approximate min-dist of: "
                  << get_brute_force_dist(bf1, bf2) << std::endl;

      } else {
        // box-cylinder case.
        auto cy2 = rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape2);
        pose_3D<double> cy2_pose = cy2->getPose().getGlobalPose();

        std::cout << "Checking proximity between Box '" << bx1->getName()
                  << "' and Cylinder '" << cy2->getName() << "'..."
                  << std::endl;

        result = geom::findProximityByGJKEPA(
            geom::box_support_func(*bx1, bx1_pose),
            geom::cylinder_support_func(*cy2, cy2_pose));

        v1 = mShape1->getPose().rotateToGlobal(
            mShape1->getPose().transformFromGlobal(result.mPoint1));
        v2 = mShape2->getPose().rotateToGlobal(
            mShape2->getPose().transformFromGlobal(result.mPoint2));

        geom::box_boundary_func bf1(*bx1, bx1_pose);
        geom::cylinder_boundary_func bf2(*cy2, cy2_pose);
        v1 = bf1(v1);
        v2 = bf2(v2);

        std::cout << " which has brute-force approximate min-dist of: "
                  << get_brute_force_dist(bf1, bf2) << std::endl;
      }
    } else {
      auto cy1 = rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape1);
      pose_3D<double> cy1_pose = cy1->getPose().getGlobalPose();
      if (mShape2->getObjectType() == geom::box::getStaticObjectType()) {
        // cylinder-box case.
        auto bx2 = rtti::rk_dynamic_ptr_cast<geom::box>(mShape2);
        pose_3D<double> bx2_pose = bx2->getPose().getGlobalPose();

        std::cout << "Checking proximity between Cylinder '" << cy1->getName()
                  << "' and Box '" << bx2->getName() << "'..." << std::endl;

        result = geom::findProximityByGJKEPA(
            geom::cylinder_support_func(*cy1, cy1_pose),
            geom::box_support_func(*bx2, bx2_pose));

        v1 = mShape1->getPose().rotateToGlobal(
            mShape1->getPose().transformFromGlobal(result.mPoint1));
        v2 = mShape2->getPose().rotateToGlobal(
            mShape2->getPose().transformFromGlobal(result.mPoint2));

        geom::cylinder_boundary_func bf1(*cy1, cy1_pose);
        geom::box_boundary_func bf2(*bx2, bx2_pose);
        v1 = bf1(v1);
        v2 = bf2(v2);

        std::cout << " which has brute-force approximate min-dist of: "
                  << get_brute_force_dist(bf1, bf2) << std::endl;

      } else {
        // cylinder-cylinder case.
        auto cy2 = rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape2);
        pose_3D<double> cy2_pose = cy2->getPose().getGlobalPose();

        std::cout << "Checking proximity between Cylinder '" << cy1->getName()
                  << "' and Cylinder '" << cy2->getName() << "'..."
                  << std::endl;

        result = geom::findProximityByGJKEPA(
            geom::cylinder_support_func(*cy1, cy1_pose),
            geom::cylinder_support_func(*cy2, cy2_pose));

        v1 = mShape1->getPose().rotateToGlobal(
            mShape1->getPose().transformFromGlobal(result.mPoint1));
        v2 = mShape2->getPose().rotateToGlobal(
            mShape2->getPose().transformFromGlobal(result.mPoint2));

        geom::cylinder_boundary_func bf1(*cy1, cy1_pose);
        geom::cylinder_boundary_func bf2(*cy2, cy2_pose);
        v1 = bf1(v1);
        v2 = bf2(v2);

        std::cout << " which has brute-force approximate min-dist of: "
                  << get_brute_force_dist(bf1, bf2) << std::endl;
      }
    }

    std::cout << "  -- Solution distance is: " << result.mDistance << std::endl;

    std::cout << "  -- Solution point-1 is: " << result.mPoint1 << std::endl;
    std::cout << "  -- Boundary at point-1 is: " << v1 << std::endl;

    std::cout << "  -- Solution point-2 is: " << result.mPoint2 << std::endl;
    std::cout << "  -- Boundary at point-2 is: " << v2 << std::endl;

    std::cout << "  -- Distance between pt-1 and pt-2 is: "
              << norm_2(result.mPoint1 - result.mPoint2) << std::endl;
  }
};

int main() {

  std::vector<proximity_solver> prox_tasks;
  prox_tasks.push_back(proximity_solver(cy1, cy4));
  prox_tasks.push_back(proximity_solver(cy1, cy5));
  prox_tasks.push_back(proximity_solver(cy1, cy6));
  prox_tasks.push_back(proximity_solver(cy2, cy4));
  prox_tasks.push_back(proximity_solver(cy3, cy4));

  prox_tasks.push_back(proximity_solver(bx1, bx4));
  prox_tasks.push_back(proximity_solver(bx1, bx5));
  prox_tasks.push_back(proximity_solver(bx1, bx6));
  prox_tasks.push_back(proximity_solver(bx2, bx4));
  prox_tasks.push_back(proximity_solver(bx3, bx4));

  prox_tasks.push_back(proximity_solver(bx1, cy4));
  prox_tasks.push_back(proximity_solver(bx2, cy4));
  prox_tasks.push_back(proximity_solver(bx3, cy4));

  prox_tasks.push_back(proximity_solver(cy1, bx4));
  prox_tasks.push_back(proximity_solver(cy2, bx4));
  prox_tasks.push_back(proximity_solver(cy3, bx4));

  prox_tasks.push_back(proximity_solver(cy7, bx7));

  for (std::size_t i = 0; i < prox_tasks.size(); ++i) {

    prox_tasks[i]();
  }

  return 0;
}
