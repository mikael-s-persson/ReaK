
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

#include "ReaK/mbd/coin3D/oi_scene_graph.h"

#include "ReaK/geometry/shapes/geometry_2D.h"
#include "ReaK/geometry/shapes/geometry_3D.h"
#include "ReaK/geometry/shapes/shape_2D.h"
#include "ReaK/geometry/shapes/shape_3D.h"

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

#include <Inventor/SbColor.h>           // for SbColor
#include <Inventor/SbVec3f.h>           // for SbVec3f
#include <Inventor/SbViewportRegion.h>  // for SbViewportRegion
#include <Inventor/SbXfBox3f.h>         // for SbXfBox3f
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/fields/SoMFColor.h>       // for SoMFColor
#include <Inventor/fields/SoMFInt32.h>       // for SoMFInt32
#include <Inventor/fields/SoMFVec3f.h>       // for SoMFVec3f
#include <Inventor/fields/SoSFFloat.h>       // for SoSFFloat
#include <Inventor/fields/SoSFRotation.h>    // for SoSFRotation
#include <Inventor/fields/SoSFVec3f.h>       // for SoSFVec3f
#include <Inventor/fields/SoSubField.h>      // for SoSFFloat::operator=, etc
#include <Inventor/nodes/SoBaseColor.h>      // for SoBaseColor
#include <Inventor/nodes/SoCoordinate3.h>    // for SoCoordinate3
#include <Inventor/nodes/SoCube.h>           // for SoCube
#include <Inventor/nodes/SoCylinder.h>       // for SoCylinder
#include <Inventor/nodes/SoLineSet.h>        // for SoLineSet
#include <Inventor/nodes/SoRotation.h>       // for SoRotation
#include <Inventor/nodes/SoSeparator.h>      // for SoSeparator
#include <Inventor/nodes/SoSphere.h>         // for SoSphere
#include <Inventor/nodes/SoSwitch.h>         // for SoSwitch
#include <Inventor/nodes/SoTransform.h>      // for SoTransform
#include <Inventor/nodes/SoTranslation.h>    // for SoTranslation
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor
#include <cmath>
#include <numbers>

namespace ReaK::geom {

void oi_scene_graph::update_anchors(void* aData, SoSensor* /*unused*/) {
  if (aData == nullptr) {
    return;
  }
  auto* SG = reinterpret_cast<oi_scene_graph*>(aData);

  using Iter2D = std::map<std::shared_ptr<pose_2D<double>>,
                          std::pair<SoSeparator*, SoTransform*>>::iterator;
  using Iter3D = std::map<std::shared_ptr<pose_3D<double>>,
                          std::pair<SoSeparator*, SoTransform*>>::iterator;

  std::unique_lock<std::recursive_mutex> lock_here(SG->mAnchorUpdatingMutex);

  for (auto& it : SG->mAnchor2DMap) {
    it.second.second->rotation.setValue(SbVec3f(0.0, 0.0, 1.0),
                                        it.first->Rotation.getAngle());
    it.second.second->translation.setValue(it.first->Position[0],
                                           it.first->Position[1], 0.0);
  }

  for (auto& it : SG->mAnchor3DMap) {
    it.second.second->rotation.setValue(it.first->Quat[1], it.first->Quat[2],
                                        it.first->Quat[3], it.first->Quat[0]);
    it.second.second->translation.setValue(
        it.first->Position[0], it.first->Position[1], it.first->Position[2]);
  }

  for (auto& mUpdateFunc : SG->mUpdateFuncs) {
    mUpdateFunc();
  }
}

void oi_scene_graph::setVisibility(bool aVisible) const {
  mRootSwitch->whichChild.setValue((aVisible ? SO_SWITCH_ALL : SO_SWITCH_NONE));
}

void oi_scene_graph::clearAll() {

  std::unique_lock<std::recursive_mutex> lock_here(mAnchorUpdatingMutex);

  mUpdateFuncs.clear();

  mAnchor2DMap.clear();
  mAnchor3DMap.clear();

  mRootSwitch->removeAllChildren();
}

void oi_scene_graph::enableAnchorUpdates() {
  mTimer->schedule();
}

void oi_scene_graph::disableAnchorUpdates() {
  mTimer->unschedule();
}

double oi_scene_graph::computeCharacteristicLength() {
  using std::sqrt;
  SbViewportRegion dummy_viewport;
  SoGetBoundingBoxAction bb_calc(dummy_viewport);
  bb_calc.apply(mRoot);
  float x = NAN;
  float y = NAN;
  float z = NAN;
  bb_calc.getXfBoundingBox().getSize(x, y, z);
  return mCharacteristicLength = sqrt((x * x + y * y + z * z) / 3.0);
}

oi_scene_graph::oi_scene_graph() {
  mRoot = new SoSeparator;
  mRoot->ref();
  mRootSwitch = new SoSwitch;
  mRoot->addChild(mRootSwitch);
  setVisibility(true);
  mTimer = new SoTimerSensor(oi_scene_graph::update_anchors, this);
}

oi_scene_graph::~oi_scene_graph() {
  mRoot->unref();
  delete mTimer;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const std::shared_ptr<pose_2D<double>>& aAnchor) {

  std::unique_lock<std::recursive_mutex> lock_here(aSG.mAnchorUpdatingMutex);

  if ((!aAnchor) ||
      (aSG.mAnchor2DMap.find(aAnchor) != aSG.mAnchor2DMap.end())) {
    return aSG;
  }
  if (!aAnchor->Parent.expired()) {
    aSG << aAnchor->Parent.lock();
  }

  auto* sep = new SoSeparator;
  auto* trans = new SoTransform;
  trans->translation.setValue(aAnchor->Position[0], aAnchor->Position[1], 0.0);
  trans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0),
                           aAnchor->Rotation.getAngle());
  sep->addChild(trans);
  aSG.mAnchor2DMap[aAnchor] = std::pair<SoSeparator*, SoTransform*>(sep, trans);

  if (aAnchor->Parent.expired()) {
    aSG.mRootSwitch->addChild(sep);
  } else {
    aSG.mAnchor2DMap[aAnchor->Parent.lock()].first->addChild(sep);
  }
  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const std::shared_ptr<pose_3D<double>>& aAnchor) {

  std::unique_lock<std::recursive_mutex> lock_here(aSG.mAnchorUpdatingMutex);

  if ((!aAnchor) ||
      (aSG.mAnchor3DMap.find(aAnchor) != aSG.mAnchor3DMap.end())) {
    return aSG;
  }
  if (!aAnchor->Parent.expired()) {
    aSG << aAnchor->Parent.lock();
  }

  auto* sep = new SoSeparator;
  auto* trans = new SoTransform;
  trans->translation.setValue(aAnchor->Position[0], aAnchor->Position[1],
                              aAnchor->Position[2]);
  trans->rotation.setValue(aAnchor->Quat[1], aAnchor->Quat[2], aAnchor->Quat[3],
                           aAnchor->Quat[0]);
  sep->addChild(trans);
  aSG.mAnchor3DMap[aAnchor] = std::pair<SoSeparator*, SoTransform*>(sep, trans);

  if (aAnchor->Parent.expired()) {
    aSG.mRootSwitch->addChild(sep);
  } else {
    aSG.mAnchor3DMap[aAnchor->Parent.lock()].first->addChild(sep);
  }
  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG, const color& aColor) {
  aSG.mCurrentColor = aColor;
  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_2D& aGeom2D) {

  auto* sep = new SoSeparator;

  auto* trans = new SoTransform;
  trans->translation.setValue(aGeom2D.getPose().Position[0],
                              aGeom2D.getPose().Position[1], 0.0);
  trans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0),
                           aGeom2D.getPose().Rotation.getAngle());
  sep->addChild(trans);

  if (aGeom2D.get_object_type() == line_seg_2D::get_static_object_type()) {
    const auto& ln_geom = static_cast<const line_seg_2D&>(aGeom2D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* coords = new SoCoordinate3;
    coords->point.set1Value(0, ln_geom.getStart()[0], ln_geom.getStart()[1],
                            0.0);
    coords->point.set1Value(1, ln_geom.getEnd()[0], ln_geom.getEnd()[1], 0.0);
    sep->addChild(coords);

    auto* ln_set = new SoLineSet;
    ln_set->numVertices.set1Value(0, 2);
    sep->addChild(ln_set);

  } else if (aGeom2D.get_object_type() == grid_2D::get_static_object_type()) {
    const auto& gd_geom = static_cast<const grid_2D&>(aGeom2D);
    double x_step =
        gd_geom.getDimensions()[0] / double(gd_geom.getSquareCounts()[0]);
    double y_step =
        gd_geom.getDimensions()[1] / double(gd_geom.getSquareCounts()[1]);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* coords = new SoCoordinate3;
    auto* ln_set = new SoLineSet;
    int j = 0;
    int k = 0;
    double x_value = -0.5 * gd_geom.getDimensions()[0];
    double y_value = -0.5 * gd_geom.getDimensions()[1];
    // First create the x-axis lines:
    for (std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
      coords->point.set1Value(j++, -0.5 * gd_geom.getDimensions()[0], y_value,
                              0.0);
      coords->point.set1Value(j++, 0.5 * gd_geom.getDimensions()[0], y_value,
                              0.0);
      y_value += y_step;
      ln_set->numVertices.set1Value(k++, 2);
    }
    // Second create the y-axis lines:
    for (std::size_t i = 0; i <= gd_geom.getSquareCounts()[0]; ++i) {
      coords->point.set1Value(j++, x_value, -0.5 * gd_geom.getDimensions()[1],
                              0.0);
      coords->point.set1Value(j++, x_value, 0.5 * gd_geom.getDimensions()[1],
                              0.0);
      x_value += x_step;
      ln_set->numVertices.set1Value(k++, 2);
    }

    sep->addChild(coords);
    sep->addChild(ln_set);

  } else if (aGeom2D.get_object_type() ==
             coord_arrows_2D::get_static_object_type()) {
    const auto& ca_geom = static_cast<const coord_arrows_2D&>(aGeom2D);

    auto* sep_x = new SoSeparator;

    auto* col_x = new SoBaseColor;
    col_x->rgb = SbColor(1, 0, 0);
    sep_x->addChild(col_x);

    auto* coords_x = new SoCoordinate3;
    coords_x->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_x->point.set1Value(1, ca_geom.getArrowLength(), 0.0, 0.0);
    sep_x->addChild(coords_x);

    auto* ln_set_x = new SoLineSet;
    ln_set_x->numVertices.set1Value(0, 2);
    sep_x->addChild(ln_set_x);

    sep->addChild(sep_x);

    auto* sep_y = new SoSeparator;

    auto* col_y = new SoBaseColor;
    col_y->rgb = SbColor(0, 1, 0);
    sep_y->addChild(col_y);

    auto* coords_y = new SoCoordinate3;
    coords_y->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_y->point.set1Value(1, 0.0, ca_geom.getArrowLength(), 0.0);
    sep_y->addChild(coords_y);

    auto* ln_set_y = new SoLineSet;
    ln_set_y->numVertices.set1Value(0, 2);
    sep_y->addChild(ln_set_y);

    sep->addChild(sep_y);

  } else if (aGeom2D.get_object_type() == circle::get_static_object_type()) {
    const auto& ci_geom = static_cast<const circle&>(aGeom2D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* ci_rot = new SoRotation;
    ci_rot->rotation.setValue(SbVec3f(1.0, 0.0, 0.0), std::numbers::pi / 2.0);
    sep->addChild(ci_rot);

    auto* ci_cyl = new SoCylinder;
    ci_cyl->radius = ci_geom.getRadius();
    ci_cyl->height = 1e-4;
    sep->addChild(ci_cyl);

  } else if (aGeom2D.get_object_type() == rectangle::get_static_object_type()) {
    const auto& re_geom = static_cast<const rectangle&>(aGeom2D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* re_cube = new SoCube;
    re_cube->width = re_geom.getDimensions()[0];
    re_cube->height = re_geom.getDimensions()[1];
    re_cube->depth = 1e-4;
    sep->addChild(re_cube);

  } else if (aGeom2D.get_object_type() ==
             capped_rectangle::get_static_object_type()) {
    const auto& re_geom = static_cast<const rectangle&>(aGeom2D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* re_cube = new SoCube;
    re_cube->width = re_geom.getDimensions()[0];
    re_cube->height = re_geom.getDimensions()[1];
    re_cube->depth = 1e-4;
    sep->addChild(re_cube);

    auto* re_rot = new SoRotation;
    re_rot->rotation.setValue(SbVec3f(1.0, 0.0, 0.0), std::numbers::pi / 2.0);
    sep->addChild(re_rot);

    auto* re_trans_right = new SoTranslation;
    re_trans_right->translation.setValue(0.5 * re_geom.getDimensions()[0], 0.0,
                                         0.0);
    sep->addChild(re_trans_right);

    auto* re_cyl_right = new SoCylinder;
    re_cyl_right->radius = 0.5 * re_geom.getDimensions()[1];
    re_cyl_right->height = 1e-4;
    sep->addChild(re_cyl_right);

    auto* re_trans_left = new SoTranslation;
    re_trans_left->translation.setValue(-re_geom.getDimensions()[0], 0.0, 0.0);
    sep->addChild(re_trans_left);

    auto* re_cyl_left = new SoCylinder;
    re_cyl_left->radius = 0.5 * re_geom.getDimensions()[1];
    re_cyl_left->height = 1e-4;
    sep->addChild(re_cyl_left);

  } else if (aGeom2D.get_object_type() ==
             composite_shape_2D::get_static_object_type()) {
    const auto& comp_geom = static_cast<const composite_shape_2D&>(aGeom2D);

    using Iter = std::vector<std::shared_ptr<shape_2D>>::const_iterator;
    for (const auto& it : comp_geom.Shapes()) {
      aSG << *it;
    }
  }

  if (!aGeom2D.getAnchor()) {
    aSG.mRootSwitch->addChild(sep);
  } else {

    std::unique_lock<std::recursive_mutex> lock_here(aSG.mAnchorUpdatingMutex);

    if (aSG.mAnchor2DMap.find(aGeom2D.getAnchor()) == aSG.mAnchor2DMap.end()) {
      aSG << aGeom2D.getAnchor();
    }
    aSG.mAnchor2DMap[aGeom2D.getAnchor()].first->addChild(sep);
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_3D& aGeom3D) {

  auto* sep = new SoSeparator;

  auto* trans = new SoTransform;
  trans->translation.setValue(aGeom3D.getPose().Position[0],
                              aGeom3D.getPose().Position[1],
                              aGeom3D.getPose().Position[2]);
  trans->rotation.setValue(aGeom3D.getPose().Quat[1], aGeom3D.getPose().Quat[2],
                           aGeom3D.getPose().Quat[3],
                           aGeom3D.getPose().Quat[0]);
  sep->addChild(trans);

  if (aGeom3D.get_object_type() == line_seg_3D::get_static_object_type()) {
    const auto& ln_geom = static_cast<const line_seg_3D&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* coords = new SoCoordinate3;
    coords->point.set1Value(0, ln_geom.getStart()[0], ln_geom.getStart()[1],
                            ln_geom.getStart()[2]);
    coords->point.set1Value(1, ln_geom.getEnd()[0], ln_geom.getEnd()[1],
                            ln_geom.getEnd()[2]);
    sep->addChild(coords);

    auto* ln_set = new SoLineSet;
    ln_set->numVertices.set1Value(0, 2);
    sep->addChild(ln_set);

  } else if (aGeom3D.get_object_type() == grid_3D::get_static_object_type()) {
    const auto& gd_geom = static_cast<const grid_3D&>(aGeom3D);
    double x_step =
        gd_geom.getDimensions()[0] / double(gd_geom.getSquareCounts()[0]);
    double y_step =
        gd_geom.getDimensions()[1] / double(gd_geom.getSquareCounts()[1]);
    double z_step =
        gd_geom.getDimensions()[2] / double(gd_geom.getSquareCounts()[2]);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* coords = new SoCoordinate3;
    auto* ln_set = new SoLineSet;
    int j = 0;
    int k = 0;
    double x_value = -0.5 * gd_geom.getDimensions()[0];
    double y_value = -0.5 * gd_geom.getDimensions()[1];
    double z_value = -0.5 * gd_geom.getDimensions()[2];
    for (std::size_t m = 0; m <= gd_geom.getSquareCounts()[2]; ++m) {
      // First create the x-axis lines:
      y_value = -0.5 * gd_geom.getDimensions()[1];
      for (std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
        coords->point.set1Value(j++, -0.5 * gd_geom.getDimensions()[0], y_value,
                                z_value);
        coords->point.set1Value(j++, 0.5 * gd_geom.getDimensions()[0], y_value,
                                z_value);
        y_value += y_step;
        ln_set->numVertices.set1Value(k++, 2);
      }
      // Second create the y-axis lines:
      x_value = -0.5 * gd_geom.getDimensions()[0];
      for (std::size_t i = 0; i <= gd_geom.getSquareCounts()[0]; ++i) {
        coords->point.set1Value(j++, x_value, -0.5 * gd_geom.getDimensions()[1],
                                z_value);
        coords->point.set1Value(j++, x_value, 0.5 * gd_geom.getDimensions()[1],
                                z_value);
        x_value += x_step;
        ln_set->numVertices.set1Value(k++, 2);
      }
      z_value += z_step;
    }

    x_value = -0.5 * gd_geom.getDimensions()[0];
    for (std::size_t m = 0; m <= gd_geom.getSquareCounts()[0]; ++m) {
      // First create the x-axis lines:
      y_value = -0.5 * gd_geom.getDimensions()[1];
      for (std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
        coords->point.set1Value(j++, x_value, y_value,
                                -0.5 * gd_geom.getDimensions()[2]);
        coords->point.set1Value(j++, x_value, y_value,
                                0.5 * gd_geom.getDimensions()[2]);
        y_value += y_step;
        ln_set->numVertices.set1Value(k++, 2);
      }
      x_value += x_step;
    }

    sep->addChild(coords);
    sep->addChild(ln_set);

  } else if (aGeom3D.get_object_type() ==
             coord_arrows_3D::get_static_object_type()) {
    const auto& ca_geom = static_cast<const coord_arrows_3D&>(aGeom3D);

    auto* sep_x = new SoSeparator;

    auto* col_x = new SoBaseColor;
    col_x->rgb = SbColor(1, 0, 0);
    sep_x->addChild(col_x);

    auto* coords_x = new SoCoordinate3;
    coords_x->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_x->point.set1Value(1, ca_geom.getArrowLength(), 0.0, 0.0);
    sep_x->addChild(coords_x);

    auto* ln_set_x = new SoLineSet;
    ln_set_x->numVertices.set1Value(0, 2);
    sep_x->addChild(ln_set_x);

    sep->addChild(sep_x);

    auto* sep_y = new SoSeparator;

    auto* col_y = new SoBaseColor;
    col_y->rgb = SbColor(0, 1, 0);
    sep_y->addChild(col_y);

    auto* coords_y = new SoCoordinate3;
    coords_y->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_y->point.set1Value(1, 0.0, ca_geom.getArrowLength(), 0.0);
    sep_y->addChild(coords_y);

    auto* ln_set_y = new SoLineSet;
    ln_set_y->numVertices.set1Value(0, 2);
    sep_y->addChild(ln_set_y);

    sep->addChild(sep_y);

    auto* sep_z = new SoSeparator;

    auto* col_z = new SoBaseColor;
    col_z->rgb = SbColor(0, 0, 1);
    sep_z->addChild(col_z);

    auto* coords_z = new SoCoordinate3;
    coords_z->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_z->point.set1Value(1, 0.0, 0.0, ca_geom.getArrowLength());
    sep_z->addChild(coords_z);

    auto* ln_set_z = new SoLineSet;
    ln_set_z->numVertices.set1Value(0, 2);
    sep_z->addChild(ln_set_z);

    sep->addChild(sep_z);

  } else if (aGeom3D.get_object_type() == plane::get_static_object_type()) {
    const auto& pl_geom = static_cast<const plane&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* pl_cube = new SoCube;
    pl_cube->width = pl_geom.getDimensions()[0];
    pl_cube->height = pl_geom.getDimensions()[1];
    pl_cube->depth = 1e-4;
    sep->addChild(pl_cube);

  } else if (aGeom3D.get_object_type() == sphere::get_static_object_type()) {
    const auto& sp_geom = static_cast<const sphere&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* sp_cyl = new SoSphere;
    sp_cyl->radius = sp_geom.getRadius();
    sep->addChild(sp_cyl);

  } else if (aGeom3D.get_object_type() == box::get_static_object_type()) {
    const box& bx_geom = static_cast<const box&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* bx_cube = new SoCube;
    bx_cube->width = bx_geom.getDimensions()[0];
    bx_cube->height = bx_geom.getDimensions()[1];
    bx_cube->depth = bx_geom.getDimensions()[2];
    sep->addChild(bx_cube);

  } else if (aGeom3D.get_object_type() == cylinder::get_static_object_type()) {
    const auto& cy_geom = static_cast<const cylinder&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* cy_rot = new SoRotation;
    cy_rot->rotation.setValue(SbVec3f(1.0, 0.0, 0.0), std::numbers::pi / 2.0);
    sep->addChild(cy_rot);

    auto* cy_cyl = new SoCylinder;
    cy_cyl->radius = cy_geom.getRadius();
    cy_cyl->height = cy_geom.getLength();
    sep->addChild(cy_cyl);

  } else if (aGeom3D.get_object_type() ==
             capped_cylinder::get_static_object_type()) {
    const auto& cy_geom = static_cast<const capped_cylinder&>(aGeom3D);

    auto* col = new SoBaseColor;
    col->rgb =
        SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);

    auto* cy_rot = new SoRotation;
    cy_rot->rotation.setValue(SbVec3f(1.0, 0.0, 0.0), std::numbers::pi / 2.0);
    sep->addChild(cy_rot);

    auto* cy_cyl = new SoCylinder;
    cy_cyl->radius = cy_geom.getRadius();
    cy_cyl->height = cy_geom.getLength();
    sep->addChild(cy_cyl);

    auto* cy_trans_top = new SoTranslation;
    cy_trans_top->translation.setValue(0.0, 0.5 * cy_geom.getLength(), 0.0);
    sep->addChild(cy_trans_top);

    auto* cy_sp_top = new SoSphere;
    cy_sp_top->radius = cy_geom.getRadius();
    sep->addChild(cy_sp_top);

    auto* cy_trans_bottom = new SoTranslation;
    cy_trans_bottom->translation.setValue(0.0, -cy_geom.getLength(), 0.0);
    sep->addChild(cy_trans_bottom);

    auto* cy_sp_bottom = new SoSphere;
    cy_sp_bottom->radius = cy_geom.getRadius();
    sep->addChild(cy_sp_bottom);

  } else if (aGeom3D.get_object_type() ==
             composite_shape_3D::get_static_object_type()) {
    const auto& comp_geom = static_cast<const composite_shape_3D&>(aGeom3D);

    using Iter = std::vector<std::shared_ptr<shape_3D>>::const_iterator;
    for (const auto& it : comp_geom.Shapes()) {
      aSG << *it;
    }
  }

  if (!aGeom3D.getAnchor()) {
    aSG.mRootSwitch->addChild(sep);
  } else {
    std::unique_lock<std::recursive_mutex> lock_here(aSG.mAnchorUpdatingMutex);

    if (aSG.mAnchor3DMap.find(aGeom3D.getAnchor()) == aSG.mAnchor3DMap.end()) {
      aSG << aGeom3D.getAnchor();
    }
    aSG.mAnchor3DMap[aGeom3D.getAnchor()].first->addChild(sep);
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const colored_model_2D& aModel) {

  for (const auto& i : aModel.mAnchorList) {
    aSG << i;
  }

  for (const auto& i : aModel.mGeomList) {
    aSG << i.mColor << *(i.mGeom);
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const colored_model_3D& aModel) {

  for (const auto& i : aModel.mAnchorList) {
    aSG << i;
  }

  for (const auto& i : aModel.mGeomList) {
    aSG << i.mColor << *(i.mGeom);
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const proxy_query_model_2D& aModel) {

  aSG << color(0.9, 0.9, 0.9);
  for (std::size_t i = 0; i < aModel.getShapeCount(); ++i) {
    aSG << aModel.getShape(i);
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const proxy_query_model_3D& aModel) {

  aSG << color(0.9, 0.9, 0.9);
  for (std::size_t i = 0; i < aModel.getShapeCount(); ++i) {
    aSG << aModel.getShape(i);
  }

  return aSG;
}

}  // namespace ReaK::geom
