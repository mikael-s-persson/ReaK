
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

#include "ReaK/geometry/proximity/proxy_query_model.h"

#include "ReaK/geometry/shapes/capped_rectangle.h"
#include "ReaK/geometry/shapes/circle.h"
#include "ReaK/geometry/shapes/composite_shape_2D.h"
#include "ReaK/geometry/shapes/coord_arrows_2D.h"
#include "ReaK/geometry/shapes/grid_2D.h"
#include "ReaK/geometry/shapes/line_seg_2D.h"
#include "ReaK/geometry/shapes/rectangle.h"

#include "ReaK/geometry/shapes/box.h"
#include "ReaK/geometry/shapes/capped_cylinder.h"
#include "ReaK/geometry/shapes/composite_shape_3D.h"
#include "ReaK/geometry/shapes/coord_arrows_3D.h"
#include "ReaK/geometry/shapes/cylinder.h"
#include "ReaK/geometry/shapes/grid_3D.h"
#include "ReaK/geometry/shapes/line_seg_3D.h"
#include "ReaK/geometry/shapes/plane.h"
#include "ReaK/geometry/shapes/sphere.h"

#include "ReaK/geometry/proximity/prox_circle_circle.h"
#include "ReaK/geometry/proximity/prox_circle_crect.h"
#include "ReaK/geometry/proximity/prox_circle_rectangle.h"
#include "ReaK/geometry/proximity/prox_crect_crect.h"
#include "ReaK/geometry/proximity/prox_crect_rectangle.h"
#include "ReaK/geometry/proximity/prox_rectangle_rectangle.h"

#include "ReaK/geometry/proximity/prox_box_box.h"
#include "ReaK/geometry/proximity/prox_ccylinder_box.h"
#include "ReaK/geometry/proximity/prox_ccylinder_ccylinder.h"
#include "ReaK/geometry/proximity/prox_ccylinder_cylinder.h"
#include "ReaK/geometry/proximity/prox_cylinder_box.h"
#include "ReaK/geometry/proximity/prox_cylinder_cylinder.h"
#include "ReaK/geometry/proximity/prox_plane_box.h"
#include "ReaK/geometry/proximity/prox_plane_ccylinder.h"
#include "ReaK/geometry/proximity/prox_plane_cylinder.h"
#include "ReaK/geometry/proximity/prox_plane_plane.h"
#include "ReaK/geometry/proximity/prox_plane_sphere.h"
#include "ReaK/geometry/proximity/prox_sphere_box.h"
#include "ReaK/geometry/proximity/prox_sphere_ccylinder.h"
#include "ReaK/geometry/proximity/prox_sphere_cylinder.h"
#include "ReaK/geometry/proximity/prox_sphere_sphere.h"

#include <algorithm>
#include <variant>

namespace ReaK::geom {

namespace {

struct convert_to_shape_2D_visitor {
  shape_2D& operator()(shape_2D& s) const { return s; }
  using result_type = shape_2D&;
};
}  // namespace

struct proxy_query_model_2D::variant_shape_cache {
  struct any_shape {
    std::variant<circle, capped_rectangle, rectangle> value;

    bool operator<(const any_shape& rhs) const {
      return this->value.index() < rhs.value.index();
    }
  };
  std::vector<any_shape> shapes;

  using iterator = std::vector<any_shape>::iterator;

  void addShape(const shape_2D& aShape) {
    any_shape tmp;
    // if the other is a circle..
    if (aShape.get_object_type() == circle::get_static_object_type()) {
      tmp.value = static_cast<const circle&>(aShape);
    }
    // if the other is a capped_rectangle..
    else if (aShape.get_object_type() ==
             capped_rectangle::get_static_object_type()) {
      tmp.value = static_cast<const capped_rectangle&>(aShape);
    }
    // if the other is a rectangle..
    else if (aShape.get_object_type() == rectangle::get_static_object_type()) {
      tmp.value = static_cast<const rectangle&>(aShape);
    };
    shapes.push_back(tmp);
    std::inplace_merge(shapes.begin(), shapes.end() - 1, shapes.end());
  };
};

proxy_query_model_2D::proxy_query_model_2D(const std::string& aName)
    : mShapeCache(new variant_shape_cache()) {
  set_name(aName);
}

proxy_query_model_2D::~proxy_query_model_2D() {
  delete mShapeCache;
}

proxy_query_model_2D& proxy_query_model_2D::addShape(
    const std::shared_ptr<shape_2D>& aShape) {
  mShapeCache->addShape(*aShape);
  return *this;
}

const shape_2D& proxy_query_model_2D::getShape(std::size_t i) const {
  return std::visit(convert_to_shape_2D_visitor(),
                    mShapeCache->shapes[i].value);
}

std::size_t proxy_query_model_2D::getShapeCount() const {
  return mShapeCache->shapes.size();
}

namespace {

struct precomputer_shape_2D_visitor {
  std::vector<shape_2D_precompute_pack>::iterator it_pre;
  explicit precomputer_shape_2D_visitor(
      std::vector<shape_2D_precompute_pack>::iterator aItPre)
      : it_pre(aItPre) {}
  void operator()(shape_2D& s) const { *it_pre = s.createPrecomputePack(); }
  using result_type = void;
};
}  // namespace

void proxy_query_model_2D::doPrecomputePass() {
  mPreComputePacks.resize(mShapeCache->shapes.size());

  auto it_pre = mPreComputePacks.begin();
  for (auto it = mShapeCache->shapes.begin(),
            it_end = mShapeCache->shapes.end();
       it != it_end; ++it_pre, ++it) {
    std::visit(precomputer_shape_2D_visitor(it_pre), it->value);
  }
}

void proxy_query_model_2D::save(ReaK::serialization::oarchive& A,
                                unsigned int /*unused*/) const {
  named_object::save(A, named_object::get_static_object_type()->version());

  // This is ugly and inefficient, but it keeps backward compatibility.
  std::vector<std::shared_ptr<shape_2D>> mShapeList;
  for (auto& shape : mShapeCache->shapes) {
    mShapeList.push_back(std::shared_ptr<shape_2D>(
        &std::visit(convert_to_shape_2D_visitor(), shape.value),
        null_deleter()));
  }

  A& RK_SERIAL_SAVE_WITH_NAME(mShapeList);
}

void proxy_query_model_2D::load(ReaK::serialization::iarchive& A,
                                unsigned int /*unused*/) {
  named_object::load(A, named_object::get_static_object_type()->version());
  std::vector<std::shared_ptr<shape_2D>> mShapeList;
  A& RK_SERIAL_LOAD_WITH_NAME(mShapeList);

  // This is ugly and inefficient, but it keeps backward compatibility.
  for (auto& it : mShapeList) {
    addShape(it);
  }
}

namespace {

struct compute_proximity_2D_visitor {
  const shape_2D_precompute_pack* p1;
  const shape_2D_precompute_pack* p2;
  compute_proximity_2D_visitor(const shape_2D_precompute_pack& aP1,
                               const shape_2D_precompute_pack& aP2)
      : p1(&aP1), p2(&aP2) {}
  template <typename T, typename U>
  proximity_record_2D operator()(const T& s1, const U& s2) const {
    return compute_proximity(s1, *p1, s2, *p2);
  }
  using result_type = proximity_record_2D;
};

struct get_bounding_radius_2D_visitor {
  double operator()(const shape_2D& s) const { return s.getBoundingRadius(); }
  using result_type = double;
};
}  // namespace

proximity_record_2D proxy_query_pair_2D::findMinimumDistance() const {
  proximity_record_2D result;

  if (!mModel1 || !mModel2) {
    return result;
  }

  mModel1->doPrecomputePass();
  mModel2->doPrecomputePass();

  for (std::size_t m1_i = 0; m1_i < mModel1->mPreComputePacks.size(); ++m1_i) {
    vect<double, 2> p1 = mModel1->mPreComputePacks[m1_i].global_pose.Position;
    double brad1 = std::visit(get_bounding_radius_2D_visitor(),
                              mModel1->mShapeCache->shapes[m1_i].value);
    for (std::size_t m2_i = 0; m2_i < mModel2->mPreComputePacks.size();
         ++m2_i) {
      vect<double, 2> p2 = mModel2->mPreComputePacks[m2_i].global_pose.Position;
      double brad2 = std::visit(get_bounding_radius_2D_visitor(),
                                mModel2->mShapeCache->shapes[m2_i].value);
      if (norm_2(p2 - p1) - brad1 - brad2 > result.mDistance) {
        continue;
      }

      proximity_record_2D tmp = std::visit(
          compute_proximity_2D_visitor(mModel1->mPreComputePacks[m1_i],
                                       mModel2->mPreComputePacks[m2_i]),
          mModel1->mShapeCache->shapes[m1_i].value,
          mModel2->mShapeCache->shapes[m2_i].value);

      if (result.mDistance > tmp.mDistance) {
        result = tmp;
      }
    }
  }

  return result;
}

bool proxy_query_pair_2D::gatherCollisionPoints(
    std::vector<proximity_record_2D>& aOutput) const {
  if (!mModel1 || !mModel2) {
    return false;
  }

  mModel1->doPrecomputePass();
  mModel2->doPrecomputePass();

  bool collision_found = false;

  for (std::size_t m1_i = 0; m1_i < mModel1->mPreComputePacks.size(); ++m1_i) {
    vect<double, 2> p1 = mModel1->mPreComputePacks[m1_i].global_pose.Position;
    double brad1 = std::visit(get_bounding_radius_2D_visitor(),
                              mModel1->mShapeCache->shapes[m1_i].value);
    for (std::size_t m2_i = 0; m2_i < mModel2->mPreComputePacks.size();
         ++m2_i) {
      vect<double, 2> p2 = mModel2->mPreComputePacks[m2_i].global_pose.Position;
      double brad2 = std::visit(get_bounding_radius_2D_visitor(),
                                mModel2->mShapeCache->shapes[m2_i].value);
      if (norm_2(p2 - p1) - brad1 - brad2 > 0.0) {
        continue;
      }

      proximity_record_2D tmp = std::visit(
          compute_proximity_2D_visitor(mModel1->mPreComputePacks[m1_i],
                                       mModel2->mPreComputePacks[m2_i]),
          mModel1->mShapeCache->shapes[m1_i].value,
          mModel2->mShapeCache->shapes[m2_i].value);

      if (tmp.mDistance < 0.0) {
        aOutput.push_back(tmp);
        collision_found = true;
      }
    }
  }

  return collision_found;
}

namespace {

struct convert_to_shape_3D_visitor {
  shape_3D& operator()(shape_3D& s) const { return s; }
  using result_type = shape_3D&;
};
};  // namespace

struct proxy_query_model_3D::variant_shape_cache {
  struct any_shape {
    std::variant<plane, sphere, box, cylinder, capped_cylinder> value;

    bool operator<(const any_shape& rhs) const {
      return this->value.index() < rhs.value.index();
    }
  };
  std::vector<any_shape> shapes;

  using iterator = std::vector<any_shape>::iterator;

  void addShape(const shape_3D& aShape) {
    any_shape tmp;
    // if the other is a plane..
    if (aShape.get_object_type() == plane::get_static_object_type()) {
      tmp.value = static_cast<const plane&>(aShape);
    }
    // if the other is a sphere..
    else if (aShape.get_object_type() == sphere::get_static_object_type()) {
      tmp.value = static_cast<const sphere&>(aShape);
    }
    // if the other is a ccylinder..
    else if (aShape.get_object_type() ==
             capped_cylinder::get_static_object_type()) {
      tmp.value = static_cast<const capped_cylinder&>(aShape);
    }
    // if the other is a cylinder..
    else if (aShape.get_object_type() == cylinder::get_static_object_type()) {
      tmp.value = static_cast<const cylinder&>(aShape);
    }
    // if the other is a box..
    else if (aShape.get_object_type() == box::get_static_object_type()) {
      tmp.value = static_cast<const box&>(aShape);
    }

    shapes.push_back(tmp);
    std::inplace_merge(shapes.begin(), shapes.end() - 1, shapes.end());
  }
};

proxy_query_model_3D::proxy_query_model_3D(const std::string& aName)
    : mShapeCache(new variant_shape_cache()) {
  set_name(aName);
}

proxy_query_model_3D::~proxy_query_model_3D() {
  delete mShapeCache;
}

proxy_query_model_3D& proxy_query_model_3D::addShape(
    const std::shared_ptr<shape_3D>& aShape) {
  mShapeCache->addShape(*aShape);
  return *this;
}

const shape_3D& proxy_query_model_3D::getShape(std::size_t i) const {
  return std::visit(convert_to_shape_3D_visitor(),
                    mShapeCache->shapes[i].value);
}

std::size_t proxy_query_model_3D::getShapeCount() const {
  return mShapeCache->shapes.size();
}

namespace {

struct precomputer_shape_3D_visitor {
  std::vector<shape_3D_precompute_pack>::iterator it_pre;
  explicit precomputer_shape_3D_visitor(
      std::vector<shape_3D_precompute_pack>::iterator aItPre)
      : it_pre(aItPre) {}
  void operator()(shape_3D& s) const { *it_pre = s.createPrecomputePack(); }
  using result_type = void;
};
};  // namespace

void proxy_query_model_3D::doPrecomputePass() {
  mPreComputePacks.resize(mShapeCache->shapes.size());

  auto it_pre = mPreComputePacks.begin();
  for (auto it = mShapeCache->shapes.begin(),
            it_end = mShapeCache->shapes.end();
       it != it_end; ++it_pre, ++it) {
    std::visit(precomputer_shape_3D_visitor(it_pre), it->value);
  }
}

void proxy_query_model_3D::save(ReaK::serialization::oarchive& A,
                                unsigned int /*unused*/) const {
  named_object::save(A, named_object::get_static_object_type()->version());

  // This is inefficient, but it keeps backward compatibility.
  std::vector<std::shared_ptr<shape_3D>> mShapeList;
  for (auto& shape : mShapeCache->shapes) {
    mShapeList.push_back(std::shared_ptr<shape_3D>(
        &std::visit(convert_to_shape_3D_visitor(), shape.value),
        null_deleter()));
  }

  A& RK_SERIAL_SAVE_WITH_NAME(mShapeList);
}

void proxy_query_model_3D::load(ReaK::serialization::iarchive& A,
                                unsigned int /*unused*/) {
  named_object::load(A, named_object::get_static_object_type()->version());
  std::vector<std::shared_ptr<shape_3D>> mShapeList;
  A& RK_SERIAL_LOAD_WITH_NAME(mShapeList);

  // This is inefficient, but it keeps backward compatibility.
  for (auto& it : mShapeList) {
    addShape(it);
  }
}

namespace {

struct compute_proximity_3D_visitor {
  const shape_3D_precompute_pack* p1;
  const shape_3D_precompute_pack* p2;
  compute_proximity_3D_visitor(const shape_3D_precompute_pack& aP1,
                               const shape_3D_precompute_pack& aP2)
      : p1(&aP1), p2(&aP2) {}
  template <typename T, typename U>
  proximity_record_3D operator()(const T& s1, const U& s2) const {
    return compute_proximity(s1, *p1, s2, *p2);
  }
  using result_type = proximity_record_3D;
};

struct get_bounding_radius_3D_visitor {
  double operator()(const shape_3D& s) const { return s.getBoundingRadius(); }
  using result_type = double;
};
}  // namespace

proximity_record_3D proxy_query_pair_3D::findMinimumDistance() const {
  proximity_record_3D result;

  if (!mModel1 || !mModel2) {
    return result;
  }

  mModel1->doPrecomputePass();
  mModel2->doPrecomputePass();

  for (std::size_t m1_i = 0; m1_i < mModel1->mPreComputePacks.size(); ++m1_i) {
    vect<double, 3> p1 = mModel1->mPreComputePacks[m1_i].global_pose.Position;
    double brad1 = std::visit(get_bounding_radius_3D_visitor(),
                              mModel1->mShapeCache->shapes[m1_i].value);
    for (std::size_t m2_i = 0; m2_i < mModel2->mPreComputePacks.size();
         ++m2_i) {
      vect<double, 3> p2 = mModel2->mPreComputePacks[m2_i].global_pose.Position;
      double brad2 = std::visit(get_bounding_radius_3D_visitor(),
                                mModel2->mShapeCache->shapes[m2_i].value);
      if (norm_2(p2 - p1) - brad1 - brad2 > result.mDistance) {
        continue;
      }

      proximity_record_3D tmp = std::visit(
          compute_proximity_3D_visitor(mModel1->mPreComputePacks[m1_i],
                                       mModel2->mPreComputePacks[m2_i]),
          mModel1->mShapeCache->shapes[m1_i].value,
          mModel2->mShapeCache->shapes[m2_i].value);

      if (result.mDistance > tmp.mDistance) {
        result = tmp;
      }
    }
  }

  return result;
}

bool proxy_query_pair_3D::gatherCollisionPoints(
    std::vector<proximity_record_3D>& aOutput) const {
  if (!mModel1 || !mModel2) {
    return false;
  }

  mModel1->doPrecomputePass();
  mModel2->doPrecomputePass();

  bool collision_found = false;

  for (std::size_t m1_i = 0; m1_i < mModel1->mPreComputePacks.size(); ++m1_i) {
    vect<double, 3> p1 = mModel1->mPreComputePacks[m1_i].global_pose.Position;
    double brad1 = std::visit(get_bounding_radius_3D_visitor(),
                              mModel1->mShapeCache->shapes[m1_i].value);
    for (std::size_t m2_i = 0; m2_i < mModel2->mPreComputePacks.size();
         ++m2_i) {
      vect<double, 3> p2 = mModel2->mPreComputePacks[m2_i].global_pose.Position;
      double brad2 = std::visit(get_bounding_radius_3D_visitor(),
                                mModel2->mShapeCache->shapes[m2_i].value);
      if (norm_2(p2 - p1) - brad1 - brad2 > 0.0) {
        continue;
      }

      proximity_record_3D tmp = std::visit(
          compute_proximity_3D_visitor(mModel1->mPreComputePacks[m1_i],
                                       mModel2->mPreComputePacks[m2_i]),
          mModel1->mShapeCache->shapes[m1_i].value,
          mModel2->mShapeCache->shapes[m2_i].value);

      if (tmp.mDistance < 0.0) {
        aOutput.push_back(tmp);
        collision_found = true;
      }
    }
  };

  return collision_found;
}

}  // namespace ReaK::geom
