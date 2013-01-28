
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

#include "oi_scene_graph.hpp"

#include "geometry_2D.hpp"
#include "geometry_3D.hpp"
#include "shape_2D.hpp"
#include "shape_3D.hpp"

#include "line_seg_2D.hpp"
#include "grid_2D.hpp"
#include "coord_arrows_2D.hpp"
#include "circle.hpp"
#include "rectangle.hpp"
#include "capped_rectangle.hpp"
#include "composite_shape_2D.hpp"

#include "line_seg_3D.hpp"
#include "grid_3D.hpp"
#include "coord_arrows_3D.hpp"
#include "plane.hpp"
#include "sphere.hpp"
#include "box.hpp"
#include "cylinder.hpp"
#include "capped_cylinder.hpp"
#include "composite_shape_3D.hpp"

#include <Inventor/SbColor.h>           // for SbColor
#include <Inventor/SbVec3f.h>           // for SbVec3f
#include <Inventor/SbViewportRegion.h>  // for SbViewportRegion
#include <Inventor/SbXfBox3f.h>         // for SbXfBox3f
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/fields/SoMFColor.h>  // for SoMFColor
#include <Inventor/fields/SoMFInt32.h>  // for SoMFInt32
#include <Inventor/fields/SoMFVec3f.h>  // for SoMFVec3f
#include <Inventor/fields/SoSFFloat.h>  // for SoSFFloat
#include <Inventor/fields/SoSFRotation.h>  // for SoSFRotation
#include <Inventor/fields/SoSFVec3f.h>  // for SoSFVec3f
#include <Inventor/fields/SoSubField.h>  // for SoSFFloat::operator=, etc
#include <Inventor/nodes/SoBaseColor.h>  // for SoBaseColor
#include <Inventor/nodes/SoCoordinate3.h>  // for SoCoordinate3
#include <Inventor/nodes/SoCube.h>      // for SoCube
#include <Inventor/nodes/SoCylinder.h>  // for SoCylinder
#include <Inventor/nodes/SoLineSet.h>   // for SoLineSet
#include <Inventor/nodes/SoRotation.h>  // for SoRotation
#include <Inventor/nodes/SoSeparator.h>  // for SoSeparator
#include <Inventor/nodes/SoSphere.h>    // for SoSphere
#include <Inventor/nodes/SoTransform.h>  // for SoTransform
#include <Inventor/nodes/SoTranslation.h>  // for SoTranslation
#include <Inventor/sensors/SoTimerSensor.h>  // for SoTimerSensor

#include <cmath>                        // for sqrt, M_PI

namespace ReaK {

namespace geom {



void oi_scene_graph::update_anchors(void* aData, SoSensor*) {
  if(!aData)
    return;
  oi_scene_graph* SG = reinterpret_cast<oi_scene_graph*>(aData);
  
  typedef std::map< shared_ptr< pose_2D<double> >, std::pair<SoSeparator*, SoTransform*> >::iterator Iter2D;
  typedef std::map< shared_ptr< pose_3D<double> >, std::pair<SoSeparator*, SoTransform*> >::iterator Iter3D;
  
  for(Iter2D it = SG->mAnchor2DMap.begin(); it != SG->mAnchor2DMap.end(); ++it) {
    it->second.second->rotation.setValue(SbVec3f(0.0,0.0,1.0), it->first->Rotation.getAngle());
    it->second.second->translation.setValue(it->first->Position[0],it->first->Position[1],0.0);
  };
  
  for(Iter3D it = SG->mAnchor3DMap.begin(); it != SG->mAnchor3DMap.end(); ++it) {
    it->second.second->rotation.setValue(it->first->Quat[1],it->first->Quat[2],it->first->Quat[3],it->first->Quat[0]);
    it->second.second->translation.setValue(it->first->Position[0],it->first->Position[1],it->first->Position[2]);
  };
  
  for(std::size_t i = 0; i < SG->mUpdateFuncs.size(); ++i) 
    SG->mUpdateFuncs[i]();
};


void oi_scene_graph::enableAnchorUpdates() {
  mTimer->schedule();
};


void oi_scene_graph::disableAnchorUpdates() {
  mTimer->unschedule();
};

double oi_scene_graph::computeCharacteristicLength() {
  using std::sqrt;
  SbViewportRegion dummy_viewport;
  SoGetBoundingBoxAction bb_calc(dummy_viewport);
  bb_calc.apply(mRoot);
  float x,y,z;
  bb_calc.getXfBoundingBox().getSize(x,y,z);
  return mCharacteristicLength = sqrt((x * x + y * y + z * z) / 3.0);
};



oi_scene_graph::oi_scene_graph() : mRoot(NULL), mTimer(NULL), mAnchor2DMap(), mAnchor3DMap() {
  mRoot = new SoSeparator;
  mRoot->ref();
  mTimer = new SoTimerSensor(oi_scene_graph::update_anchors, this);
};

oi_scene_graph::~oi_scene_graph() {
  mRoot->unref();
  delete mTimer;
};

oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_2D<double> >& aAnchor) {
  if((!aAnchor) || (aSG.mAnchor2DMap.find(aAnchor) != aSG.mAnchor2DMap.end()))
    return aSG;
  if(!aAnchor->Parent.expired())
    aSG << aAnchor->Parent.lock();
  
  SoSeparator* sep = new SoSeparator;
  SoTransform* trans = new SoTransform;
  trans->translation.setValue(aAnchor->Position[0], aAnchor->Position[1], 0.0);
  trans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0), aAnchor->Rotation.getAngle());
  sep->addChild(trans);
  aSG.mAnchor2DMap[aAnchor] = std::pair<SoSeparator*,SoTransform*>(sep,trans);
  
  if(aAnchor->Parent.expired()) {
    aSG.mRoot->addChild(sep);
  } else {
    aSG.mAnchor2DMap[aAnchor->Parent.lock()].first->addChild(sep);
  };
  return aSG;
};

oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_3D<double> >& aAnchor) {
  if((!aAnchor) || (aSG.mAnchor3DMap.find(aAnchor) != aSG.mAnchor3DMap.end()))
    return aSG;
  if(!aAnchor->Parent.expired())
    aSG << aAnchor->Parent.lock();
  
  SoSeparator* sep = new SoSeparator;
  SoTransform* trans = new SoTransform;
  trans->translation.setValue(aAnchor->Position[0], aAnchor->Position[1], aAnchor->Position[2]);
  trans->rotation.setValue(aAnchor->Quat[1], aAnchor->Quat[2], aAnchor->Quat[3], aAnchor->Quat[0]);
  sep->addChild(trans);
  aSG.mAnchor3DMap[aAnchor] = std::pair<SoSeparator*,SoTransform*>(sep,trans);
  
  if(aAnchor->Parent.expired()) {
    aSG.mRoot->addChild(sep);
  } else {
    aSG.mAnchor3DMap[aAnchor->Parent.lock()].first->addChild(sep);
  };
  return aSG;
};

oi_scene_graph& operator<<(oi_scene_graph& aSG, const color& aColor) {
  aSG.mCurrentColor = aColor;
  return aSG;
};

oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_2D& aGeom2D) {
  
  SoSeparator* sep = new SoSeparator;
  
  SoTransform* trans = new SoTransform;
  trans->translation.setValue(aGeom2D.getPose().Position[0], aGeom2D.getPose().Position[1], 0.0);
  trans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0), aGeom2D.getPose().Rotation.getAngle());
  sep->addChild(trans);
  
  if(aGeom2D.getObjectType() == line_seg_2D::getStaticObjectType()) {
    const line_seg_2D& ln_geom = static_cast<const line_seg_2D&>(aGeom2D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCoordinate3* coords = new SoCoordinate3;
    coords->point.set1Value(0, ln_geom.getStart()[0], ln_geom.getStart()[1], 0.0);
    coords->point.set1Value(1, ln_geom.getEnd()[0], ln_geom.getEnd()[1], 0.0);
    sep->addChild(coords);
    
    SoLineSet* ln_set = new SoLineSet;
    ln_set->numVertices.set1Value(0, 2);
    sep->addChild(ln_set);
    
  } else if(aGeom2D.getObjectType() == grid_2D::getStaticObjectType()) {
    const grid_2D& gd_geom = static_cast<const grid_2D&>(aGeom2D);
    double x_step = gd_geom.getDimensions()[0] / double(gd_geom.getSquareCounts()[0]);
    double y_step = gd_geom.getDimensions()[1] / double(gd_geom.getSquareCounts()[1]);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCoordinate3* coords = new SoCoordinate3;
    SoLineSet* ln_set = new SoLineSet;
    int j = 0;
    int k = 0;
    double x_value = -0.5 * gd_geom.getDimensions()[0];
    double y_value = -0.5 * gd_geom.getDimensions()[1];
    // First create the x-axis lines:
    for(std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
      coords->point.set1Value(j++, -0.5 * gd_geom.getDimensions()[0], y_value, 0.0);
      coords->point.set1Value(j++,  0.5 * gd_geom.getDimensions()[0], y_value, 0.0);
      y_value += y_step;
      ln_set->numVertices.set1Value(k++, 2);
    };
    // Second create the y-axis lines:
    for(std::size_t i = 0; i <= gd_geom.getSquareCounts()[0]; ++i) {
      coords->point.set1Value(j++, x_value, -0.5 * gd_geom.getDimensions()[1], 0.0);
      coords->point.set1Value(j++, x_value,  0.5 * gd_geom.getDimensions()[1], 0.0);
      x_value += x_step;
      ln_set->numVertices.set1Value(k++, 2);
    };
    
    sep->addChild(coords);
    sep->addChild(ln_set);
    
  } else if(aGeom2D.getObjectType() == coord_arrows_2D::getStaticObjectType()) {
    const coord_arrows_2D& ca_geom = static_cast<const coord_arrows_2D&>(aGeom2D);
    
    
    SoSeparator* sep_x = new SoSeparator;
    
    SoBaseColor * col_x = new SoBaseColor;
    col_x->rgb = SbColor(1, 0, 0);
    sep_x->addChild(col_x);
    
    SoCoordinate3* coords_x = new SoCoordinate3;
    coords_x->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_x->point.set1Value(1, ca_geom.getArrowLength(), 0.0, 0.0);
    sep_x->addChild(coords_x);
    
    SoLineSet* ln_set_x = new SoLineSet;
    ln_set_x->numVertices.set1Value(0, 2);
    sep_x->addChild(ln_set_x);
    
    sep->addChild(sep_x);
    
    
    SoSeparator* sep_y = new SoSeparator;
    
    SoBaseColor * col_y = new SoBaseColor;
    col_y->rgb = SbColor(0, 1, 0);
    sep_y->addChild(col_y);
    
    SoCoordinate3* coords_y = new SoCoordinate3;
    coords_y->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_y->point.set1Value(1, 0.0, ca_geom.getArrowLength(), 0.0);
    sep_y->addChild(coords_y);
    
    SoLineSet* ln_set_y = new SoLineSet;
    ln_set_y->numVertices.set1Value(0, 2);
    sep_y->addChild(ln_set_y);
    
    sep->addChild(sep_y);
    
  } else if(aGeom2D.getObjectType() == circle::getStaticObjectType()) {
    const circle& ci_geom = static_cast<const circle&>(aGeom2D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoRotation* ci_rot = new SoRotation;
    ci_rot->rotation.setValue(SbVec3f(1.0,0.0,0.0),M_PI / 2.0);
    sep->addChild(ci_rot);
    
    SoCylinder* ci_cyl = new SoCylinder;
    ci_cyl->radius = ci_geom.getRadius();
    ci_cyl->height = 1e-4;
    sep->addChild(ci_cyl);
    
  } else if(aGeom2D.getObjectType() == rectangle::getStaticObjectType()) {
    const rectangle& re_geom = static_cast<const rectangle&>(aGeom2D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCube* re_cube = new SoCube;
    re_cube->width = re_geom.getDimensions()[0];
    re_cube->height = re_geom.getDimensions()[1];
    re_cube->depth = 1e-4;
    sep->addChild(re_cube);
    
  } else if(aGeom2D.getObjectType() == capped_rectangle::getStaticObjectType()) {
    const rectangle& re_geom = static_cast<const rectangle&>(aGeom2D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCube* re_cube = new SoCube;
    re_cube->width = re_geom.getDimensions()[0];
    re_cube->height = re_geom.getDimensions()[1];
    re_cube->depth = 1e-4;
    sep->addChild(re_cube);
    
    SoRotation* re_rot = new SoRotation;
    re_rot->rotation.setValue(SbVec3f(1.0,0.0,0.0),M_PI / 2.0);
    sep->addChild(re_rot);
    
    SoTranslation* re_trans_right = new SoTranslation;
    re_trans_right->translation.setValue(0.5 * re_geom.getDimensions()[0],0.0,0.0);
    sep->addChild(re_trans_right);
    
    SoCylinder* re_cyl_right = new SoCylinder;
    re_cyl_right->radius = 0.5 * re_geom.getDimensions()[1];
    re_cyl_right->height = 1e-4;
    sep->addChild(re_cyl_right);
    
    SoTranslation* re_trans_left = new SoTranslation;
    re_trans_left->translation.setValue(-re_geom.getDimensions()[0],0.0,0.0);
    sep->addChild(re_trans_left);
    
    SoCylinder* re_cyl_left = new SoCylinder;
    re_cyl_left->radius = 0.5 * re_geom.getDimensions()[1];
    re_cyl_left->height = 1e-4;
    sep->addChild(re_cyl_left);
    
  } else if(aGeom2D.getObjectType() == composite_shape_2D::getStaticObjectType()) {
    const composite_shape_2D& comp_geom = static_cast<const composite_shape_2D&>(aGeom2D);
    
    typedef std::vector< shared_ptr< shape_2D > >::const_iterator Iter;
    for(Iter it = comp_geom.Shapes().begin(); it != comp_geom.Shapes().end(); ++it)
      aSG << *(*it);
    
  };
  
  if(!aGeom2D.getAnchor()) {
    aSG.mRoot->addChild(sep);
  } else {
    if(aSG.mAnchor2DMap.find(aGeom2D.getAnchor()) == aSG.mAnchor2DMap.end())
      aSG << aGeom2D.getAnchor();
    aSG.mAnchor2DMap[aGeom2D.getAnchor()].first->addChild(sep);
  };
  
  return aSG;
};

oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_3D& aGeom3D) {
  
  SoSeparator* sep = new SoSeparator;
  
  SoTransform* trans = new SoTransform;
  trans->translation.setValue(aGeom3D.getPose().Position[0], aGeom3D.getPose().Position[1], aGeom3D.getPose().Position[2]);
  trans->rotation.setValue(aGeom3D.getPose().Quat[1], aGeom3D.getPose().Quat[2], aGeom3D.getPose().Quat[3], aGeom3D.getPose().Quat[0]);
  sep->addChild(trans);
  
  if(aGeom3D.getObjectType() == line_seg_3D::getStaticObjectType()) {
    const line_seg_3D& ln_geom = static_cast<const line_seg_3D&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCoordinate3* coords = new SoCoordinate3;
    coords->point.set1Value(0, ln_geom.getStart()[0], ln_geom.getStart()[1], ln_geom.getStart()[2]);
    coords->point.set1Value(1, ln_geom.getEnd()[0], ln_geom.getEnd()[1], ln_geom.getEnd()[2]);
    sep->addChild(coords);
    
    SoLineSet* ln_set = new SoLineSet;
    ln_set->numVertices.set1Value(0, 2);
    sep->addChild(ln_set);
    
  } else if(aGeom3D.getObjectType() == grid_3D::getStaticObjectType()) {
    const grid_3D& gd_geom = static_cast<const grid_3D&>(aGeom3D);
    double x_step = gd_geom.getDimensions()[0] / double(gd_geom.getSquareCounts()[0]);
    double y_step = gd_geom.getDimensions()[1] / double(gd_geom.getSquareCounts()[1]);
    double z_step = gd_geom.getDimensions()[2] / double(gd_geom.getSquareCounts()[2]);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCoordinate3* coords = new SoCoordinate3;
    SoLineSet* ln_set = new SoLineSet;
    int j = 0;
    int k = 0;
    double x_value = -0.5 * gd_geom.getDimensions()[0];
    double y_value = -0.5 * gd_geom.getDimensions()[1];
    double z_value = -0.5 * gd_geom.getDimensions()[2];
    for(std::size_t m = 0; m <= gd_geom.getSquareCounts()[2]; ++m) {
      // First create the x-axis lines:
      y_value = -0.5 * gd_geom.getDimensions()[1];
      for(std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
        coords->point.set1Value(j++, -0.5 * gd_geom.getDimensions()[0], y_value, z_value);
        coords->point.set1Value(j++,  0.5 * gd_geom.getDimensions()[0], y_value, z_value);
        y_value += y_step;
        ln_set->numVertices.set1Value(k++, 2);
      };
      // Second create the y-axis lines:
      x_value = -0.5 * gd_geom.getDimensions()[0];
      for(std::size_t i = 0; i <= gd_geom.getSquareCounts()[0]; ++i) {
        coords->point.set1Value(j++, x_value, -0.5 * gd_geom.getDimensions()[1], z_value);
        coords->point.set1Value(j++, x_value,  0.5 * gd_geom.getDimensions()[1], z_value);
        x_value += x_step;
        ln_set->numVertices.set1Value(k++, 2);
      };
      z_value += z_step;
    };
    
    x_value = -0.5 * gd_geom.getDimensions()[0];
    for(std::size_t m = 0; m <= gd_geom.getSquareCounts()[0]; ++m) {
      // First create the x-axis lines:
      y_value = -0.5 * gd_geom.getDimensions()[1];
      for(std::size_t i = 0; i <= gd_geom.getSquareCounts()[1]; ++i) {
        coords->point.set1Value(j++, x_value, y_value, -0.5 * gd_geom.getDimensions()[2]);
        coords->point.set1Value(j++, x_value, y_value,  0.5 * gd_geom.getDimensions()[2]);
        y_value += y_step;
        ln_set->numVertices.set1Value(k++, 2);
      };
      x_value += x_step;
    };
    
    sep->addChild(coords);
    sep->addChild(ln_set);
    
  } else if(aGeom3D.getObjectType() == coord_arrows_3D::getStaticObjectType()) {
    const coord_arrows_3D& ca_geom = static_cast<const coord_arrows_3D&>(aGeom3D);
    
    SoSeparator* sep_x = new SoSeparator;
    
    SoBaseColor * col_x = new SoBaseColor;
    col_x->rgb = SbColor(1, 0, 0);
    sep_x->addChild(col_x);
    
    SoCoordinate3* coords_x = new SoCoordinate3;
    coords_x->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_x->point.set1Value(1, ca_geom.getArrowLength(), 0.0, 0.0);
    sep_x->addChild(coords_x);
    
    SoLineSet* ln_set_x = new SoLineSet;
    ln_set_x->numVertices.set1Value(0, 2);
    sep_x->addChild(ln_set_x);
    
    sep->addChild(sep_x);
    
    
    SoSeparator* sep_y = new SoSeparator;
    
    SoBaseColor * col_y = new SoBaseColor;
    col_y->rgb = SbColor(0, 1, 0);
    sep_y->addChild(col_y);
    
    SoCoordinate3* coords_y = new SoCoordinate3;
    coords_y->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_y->point.set1Value(1, 0.0, ca_geom.getArrowLength(), 0.0);
    sep_y->addChild(coords_y);
    
    SoLineSet* ln_set_y = new SoLineSet;
    ln_set_y->numVertices.set1Value(0, 2);
    sep_y->addChild(ln_set_y);
    
    sep->addChild(sep_y);
    
    
    SoSeparator* sep_z = new SoSeparator;
    
    SoBaseColor * col_z = new SoBaseColor;
    col_z->rgb = SbColor(0, 0, 1);
    sep_z->addChild(col_z);
    
    SoCoordinate3* coords_z = new SoCoordinate3;
    coords_z->point.set1Value(0, 0.0, 0.0, 0.0);
    coords_z->point.set1Value(1, 0.0, 0.0, ca_geom.getArrowLength());
    sep_z->addChild(coords_z);
    
    SoLineSet* ln_set_z = new SoLineSet;
    ln_set_z->numVertices.set1Value(0, 2);
    sep_z->addChild(ln_set_z);
    
    sep->addChild(sep_z);
    
  } else if(aGeom3D.getObjectType() == plane::getStaticObjectType()) {
    const plane& pl_geom = static_cast<const plane&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCube* pl_cube = new SoCube;
    pl_cube->width = pl_geom.getDimensions()[0];
    pl_cube->height = pl_geom.getDimensions()[1];
    pl_cube->depth = 1e-4;
    sep->addChild(pl_cube);
    
  } else if(aGeom3D.getObjectType() == sphere::getStaticObjectType()) {
    const sphere& sp_geom = static_cast<const sphere&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoSphere* sp_cyl = new SoSphere;
    sp_cyl->radius = sp_geom.getRadius();
    sep->addChild(sp_cyl);
    
  } else if(aGeom3D.getObjectType() == box::getStaticObjectType()) {
    const box& bx_geom = static_cast<const box&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoCube* bx_cube = new SoCube;
    bx_cube->width = bx_geom.getDimensions()[0];
    bx_cube->height = bx_geom.getDimensions()[1];
    bx_cube->depth = bx_geom.getDimensions()[2];
    sep->addChild(bx_cube);
    
  } else if(aGeom3D.getObjectType() == cylinder::getStaticObjectType()) {
    const cylinder& cy_geom = static_cast<const cylinder&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoRotation* cy_rot = new SoRotation;
    cy_rot->rotation.setValue(SbVec3f(1.0,0.0,0.0),M_PI / 2.0);
    sep->addChild(cy_rot);
    
    SoCylinder* cy_cyl = new SoCylinder;
    cy_cyl->radius = cy_geom.getRadius();
    cy_cyl->height = cy_geom.getLength();
    sep->addChild(cy_cyl);
    
  } else if(aGeom3D.getObjectType() == capped_cylinder::getStaticObjectType()) {
    const capped_cylinder& cy_geom = static_cast<const capped_cylinder&>(aGeom3D);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(aSG.mCurrentColor.R, aSG.mCurrentColor.G, aSG.mCurrentColor.B);
    sep->addChild(col);
    
    SoRotation* cy_rot = new SoRotation;
    cy_rot->rotation.setValue(SbVec3f(1.0,0.0,0.0),M_PI / 2.0);
    sep->addChild(cy_rot);
    
    SoCylinder* cy_cyl = new SoCylinder;
    cy_cyl->radius = cy_geom.getRadius();
    cy_cyl->height = cy_geom.getLength();
    sep->addChild(cy_cyl);
    
    SoTranslation* cy_trans_top = new SoTranslation;
    cy_trans_top->translation.setValue(0.0,0.5 * cy_geom.getLength(),0.0);
    sep->addChild(cy_trans_top);
    
    SoSphere* cy_sp_top = new SoSphere;
    cy_sp_top->radius = cy_geom.getRadius();
    sep->addChild(cy_sp_top);
    
    SoTranslation* cy_trans_bottom = new SoTranslation;
    cy_trans_bottom->translation.setValue(0.0,-cy_geom.getLength(),0.0);
    sep->addChild(cy_trans_bottom);
    
    SoSphere* cy_sp_bottom = new SoSphere;
    cy_sp_bottom->radius = cy_geom.getRadius();
    sep->addChild(cy_sp_bottom);
    
  } else if(aGeom3D.getObjectType() == composite_shape_3D::getStaticObjectType()) {
    const composite_shape_3D& comp_geom = static_cast<const composite_shape_3D&>(aGeom3D);
    
    typedef std::vector< shared_ptr< shape_3D > >::const_iterator Iter;
    for(Iter it = comp_geom.Shapes().begin(); it != comp_geom.Shapes().end(); ++it)
      aSG << *(*it);
    
  };
  
  if(!aGeom3D.getAnchor()) {
    aSG.mRoot->addChild(sep);
  } else {
    if(aSG.mAnchor3DMap.find(aGeom3D.getAnchor()) == aSG.mAnchor3DMap.end())
      aSG << aGeom3D.getAnchor();
    aSG.mAnchor3DMap[aGeom3D.getAnchor()].first->addChild(sep);
  };
  
  return aSG;
};



oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_2D& aModel) {
  
  for(std::size_t i = 0; i < aModel.mAnchorList.size(); ++i)
    aSG << aModel.mAnchorList[i];
  
  for(std::size_t i = 0; i < aModel.mGeomList.size(); ++i)
    aSG << aModel.mGeomList[i].mColor << *(aModel.mGeomList[i].mGeom);
  
  return aSG;
};
    
oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_3D& aModel) {
  
  for(std::size_t i = 0; i < aModel.mAnchorList.size(); ++i)
    aSG << aModel.mAnchorList[i];
  
  for(std::size_t i = 0; i < aModel.mGeomList.size(); ++i)
    aSG << aModel.mGeomList[i].mColor << *(aModel.mGeomList[i].mGeom);
  
  return aSG;
};



};


};





