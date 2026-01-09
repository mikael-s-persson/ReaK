
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

#include "ReaK/math/optimization/nl_interior_points_methods.h"

#include "ReaK/geometry/proximity/prox_fundamentals_3D.h"

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

struct proximity_solver {
  std::shared_ptr<geom::shape_3D> mShape1;
  std::shared_ptr<geom::shape_3D> mShape2;

  proximity_solver(const std::shared_ptr<geom::shape_3D>& aShape1,
                   const std::shared_ptr<geom::shape_3D>& aShape2)
      : mShape1(aShape1), mShape2(aShape2){};

  vect<double, 3> operator()() const {
    vect_n<double> x(0.0, 0.0, 0.0, 0.0);
    vect<double, 3> c1 =
        mShape1->getPose().transformToGlobal(vect<double, 3>(0.0, 0.0, 0.0));
    vect<double, 3> c2 =
        mShape2->getPose().transformToGlobal(vect<double, 3>(0.0, 0.0, 0.0));

    // pick the middle between centers as the starting guess at the solution.
    x[1] = (c1[0] + c2[0]) * 0.5;
    x[2] = (c1[1] + c2[1]) * 0.5;
    x[3] = (c1[2] + c2[2]) * 0.5;

    if (mShape1->get_object_type() == geom::box::get_static_object_type()) {
      std::shared_ptr<geom::box> bx_a =
          rtti::rk_dynamic_ptr_cast<geom::box>(mShape1);
      pose_3D<double> bx_a_pose = bx_a->getPose().getGlobalPose();
      double min_dim = 0.5 * bx_a->getDimensions()[0];
      if (min_dim > 0.5 * bx_a->getDimensions()[1]) {
        min_dim = 0.5 * bx_a->getDimensions()[1];
      }
      if (min_dim > 0.5 * bx_a->getDimensions()[2]) {
        min_dim = 0.5 * bx_a->getDimensions()[2];
      }
      if (mShape2->get_object_type() == geom::box::get_static_object_type()) {
        // box-box case.
        std::shared_ptr<geom::box> bx_b =
            rtti::rk_dynamic_ptr_cast<geom::box>(mShape2);
        pose_3D<double> bx_b_pose = bx_b->getPose().getGlobalPose();
        if (min_dim > 0.5 * bx_b->getDimensions()[0]) {
          min_dim = 0.5 * bx_b->getDimensions()[0];
        }
        if (min_dim > 0.5 * bx_b->getDimensions()[1]) {
          min_dim = 0.5 * bx_b->getDimensions()[1];
        }
        if (min_dim > 0.5 * bx_b->getDimensions()[2]) {
          min_dim = 0.5 * bx_b->getDimensions()[2];
        }

        x[0] = norm_2(c2 - c1) / min_dim;

        std::cout << "Checking proximity between Box '" << bx_a->get_name()
                  << "' and Box '" << bx_b->get_name() << "'..." << std::endl;

        try {
          optim::make_nlip_newton_tr(
              geom::slack_minimize_func(), geom::slack_minimize_grad(),
              geom::slack_minimize_hess(), min_dim, 0.1, 300, 1e-4, 1e-3, 0.9)
              .set_ineq_constraints(
                  geom::dual_slacking_func<geom::box_slacking_func,
                                           geom::box_slacking_func>(
                      geom::box_slacking_func(*bx_a, bx_a_pose),
                      geom::box_slacking_func(*bx_b, bx_b_pose)),
                  geom::dual_slacking_jac<geom::box_slacking_jac,
                                          geom::box_slacking_jac>(
                      geom::box_slacking_jac(*bx_a, bx_a_pose),
                      geom::box_slacking_jac(*bx_b, bx_b_pose)))(x);
        } catch (...) {};

        std::cout << "  -- The raw solution obtained was: " << x << std::endl;
        std::cout << "  -- The shape1 boundary functions give: "
                  << geom::box_slacking_func(*bx_a, bx_a_pose)(x) << std::endl;
        std::cout << "  -- The shape2 boundary functions give: "
                  << geom::box_slacking_func(*bx_b, bx_b_pose)(x) << std::endl;

      } else {
        // box-cylinder case.
        std::shared_ptr<geom::cylinder> cy_b =
            rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape2);
        pose_3D<double> cy_b_pose = cy_b->getPose().getGlobalPose();
        if (min_dim > 0.5 * cy_b->getLength()) {
          min_dim = 0.5 * cy_b->getLength();
        }
        if (min_dim > cy_b->getRadius()) {
          min_dim = cy_b->getRadius();
        }

        x[0] = norm_2(c2 - c1) / min_dim;

        std::cout << "Checking proximity between Box '" << bx_a->get_name()
                  << "' and Cylinder '" << cy_b->get_name() << "'..."
                  << std::endl;

        try {
          optim::make_nlip_newton_tr(
              geom::slack_minimize_func(), geom::slack_minimize_grad(),
              geom::slack_minimize_hess(), min_dim, 0.1, 300, 1e-4, 1e-3, 0.9)
              .set_ineq_constraints(
                  geom::dual_slacking_func<geom::box_slacking_func,
                                           geom::cylinder_slacking_func>(
                      geom::box_slacking_func(*bx_a, bx_a_pose),
                      geom::cylinder_slacking_func(*cy_b, cy_b_pose)),
                  geom::dual_slacking_jac<geom::box_slacking_jac,
                                          geom::cylinder_slacking_jac>(
                      geom::box_slacking_jac(*bx_a, bx_a_pose),
                      geom::cylinder_slacking_jac(*cy_b, cy_b_pose)))(x);
        } catch (...) {};

        std::cout << "  -- The raw solution obtained was: " << x << std::endl;
        std::cout << "  -- The shape1 boundary functions give: "
                  << geom::box_slacking_func(*bx_a, bx_a_pose)(x) << std::endl;
        std::cout << "  -- The shape2 boundary functions give: "
                  << geom::cylinder_slacking_func(*cy_b, cy_b_pose)(x)
                  << std::endl;
      };
    } else {
      std::shared_ptr<geom::cylinder> cy_a =
          rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape1);
      pose_3D<double> cy_a_pose = cy_a->getPose().getGlobalPose();
      double min_dim = 0.5 * cy_a->getLength();
      if (min_dim > cy_a->getRadius()) {
        min_dim = cy_a->getRadius();
      }
      if (mShape2->get_object_type() == geom::box::get_static_object_type()) {
        // cylinder-box case.
        std::shared_ptr<geom::box> bx_b =
            rtti::rk_dynamic_ptr_cast<geom::box>(mShape2);
        pose_3D<double> bx_b_pose = bx_b->getPose().getGlobalPose();
        if (min_dim > 0.5 * bx_b->getDimensions()[0]) {
          min_dim = 0.5 * bx_b->getDimensions()[0];
        }
        if (min_dim > 0.5 * bx_b->getDimensions()[1]) {
          min_dim = 0.5 * bx_b->getDimensions()[1];
        }
        if (min_dim > 0.5 * bx_b->getDimensions()[2]) {
          min_dim = 0.5 * bx_b->getDimensions()[2];
        }

        x[0] = norm_2(c2 - c1) / min_dim;

        std::cout << "Checking proximity between Cylinder '" << cy_a->get_name()
                  << "' and Box '" << bx_b->get_name() << "'..." << std::endl;

        try {
          optim::make_nlip_newton_tr(
              geom::slack_minimize_func(), geom::slack_minimize_grad(),
              geom::slack_minimize_hess(), min_dim, 0.1, 300, 1e-4, 1e-3, 0.9)
              .set_ineq_constraints(
                  geom::dual_slacking_func<geom::cylinder_slacking_func,
                                           geom::box_slacking_func>(
                      geom::cylinder_slacking_func(*cy_a, cy_a_pose),
                      geom::box_slacking_func(*bx_b, bx_b_pose)),
                  geom::dual_slacking_jac<geom::cylinder_slacking_jac,
                                          geom::box_slacking_jac>(
                      geom::cylinder_slacking_jac(*cy_a, cy_a_pose),
                      geom::box_slacking_jac(*bx_b, bx_b_pose)))(x);
        } catch (...) {};

        std::cout << "  -- The raw solution obtained was: " << x << std::endl;
        std::cout << "  -- The shape1 boundary functions give: "
                  << geom::cylinder_slacking_func(*cy_a, cy_a_pose)(x)
                  << std::endl;
        std::cout << "  -- The shape2 boundary functions give: "
                  << geom::box_slacking_func(*bx_b, bx_b_pose)(x) << std::endl;

      } else {
        // cylinder-cylinder case.
        std::shared_ptr<geom::cylinder> cy_b =
            rtti::rk_dynamic_ptr_cast<geom::cylinder>(mShape2);
        pose_3D<double> cy_b_pose = cy_b->getPose().getGlobalPose();
        if (min_dim > 0.5 * cy_b->getLength()) {
          min_dim = 0.5 * cy_b->getLength();
        }
        if (min_dim > cy_b->getRadius()) {
          min_dim = cy_b->getRadius();
        }

        x[0] = norm_2(c2 - c1) / min_dim;

        std::cout << "Checking proximity between Cylinder '" << cy_a->get_name()
                  << "' and Cylinder '" << cy_b->get_name() << "'..."
                  << std::endl;

        try {
          optim::make_nlip_newton_tr(
              geom::slack_minimize_func(), geom::slack_minimize_grad(),
              geom::slack_minimize_hess(), min_dim, 0.1, 300, 1e-4, 1e-3, 0.9)
              .set_ineq_constraints(
                  geom::dual_slacking_func<geom::cylinder_slacking_func,
                                           geom::cylinder_slacking_func>(
                      geom::cylinder_slacking_func(*cy_a, cy_a_pose),
                      geom::cylinder_slacking_func(*cy_b, cy_b_pose)),
                  geom::dual_slacking_jac<geom::cylinder_slacking_jac,
                                          geom::cylinder_slacking_jac>(
                      geom::cylinder_slacking_jac(*cy_a, cy_a_pose),
                      geom::cylinder_slacking_jac(*cy_b, cy_b_pose)))(x);
        } catch (...) {};

        std::cout << "  -- The raw solution obtained was: " << x << std::endl;
        std::cout << "  -- The shape1 boundary functions give: "
                  << geom::cylinder_slacking_func(*cy_a, cy_a_pose)(x)
                  << std::endl;
        std::cout << "  -- The shape2 boundary functions give: "
                  << geom::cylinder_slacking_func(*cy_b, cy_b_pose)(x)
                  << std::endl;
      };
    };

    vect<double, 3> result(x[1], x[2], x[3]);
    vect<double, 3> x1_rel = mShape1->getPose().transformFromGlobal(
        vect<double, 3>(x[1], x[2], x[3]));
    vect<double, 3> x2_rel = mShape2->getPose().transformFromGlobal(
        vect<double, 3>(x[1], x[2], x[3]));

    if (x[0] > 1e-6) {
      x1_rel *= (1.0 / x[0]);
      x2_rel *= (1.0 / x[0]);
      vect<double, 3> x1 = mShape1->getPose().transformToGlobal(x1_rel);
      vect<double, 3> x2 = mShape2->getPose().transformToGlobal(x2_rel);
      std::cout << "  -- The point on Shape1 is " << x1 << std::endl;
      std::cout << "  -- The point on Shape2 is " << x2 << std::endl;
      std::cout << "  -- The distance is "
                << (norm_2(x2 - x1) * (x[0] < 1.0 ? -1.0 : 1.0)) << std::endl;
    } else {
      std::cout << "  -- The center points of the shapes are coincident! At "
                << vect<double, 3>(x[1], x[2], x[3]) << std::endl;
    };

    return result;
  };
};

int main() {

  std::vector<proximity_solver> prox_tasks;
  prox_tasks.emplace_back(cy1, cy4);
  prox_tasks.emplace_back(cy1, cy5);
  prox_tasks.emplace_back(cy1, cy6);
  prox_tasks.emplace_back(cy2, cy4);
  prox_tasks.emplace_back(cy3, cy4);

  prox_tasks.emplace_back(bx1, bx4);
  prox_tasks.emplace_back(bx1, bx5);
  prox_tasks.emplace_back(bx1, bx6);
  prox_tasks.emplace_back(bx2, bx4);
  prox_tasks.emplace_back(bx3, bx4);

  prox_tasks.emplace_back(bx1, cy4);
  prox_tasks.emplace_back(bx2, cy4);
  prox_tasks.emplace_back(bx3, cy4);

  prox_tasks.emplace_back(cy1, bx4);
  prox_tasks.emplace_back(cy2, bx4);
  prox_tasks.emplace_back(cy3, bx4);

  for (auto& prox_task : prox_tasks) {

    prox_task();
  };

  return 0;
};
