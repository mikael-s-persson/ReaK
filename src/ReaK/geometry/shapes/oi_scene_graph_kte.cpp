
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

#include "line_seg_2D.hpp"
#include "grid_2D.hpp"
#include "coord_arrows_2D.hpp"
#include "circle.hpp"
#include "rectangle.hpp"
#include "capped_rectangle.hpp"
#include "composite_shape_2D.hpp"

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/free_joint.hpp"
#include "mbd_kte/revolute_joint.hpp"
#include "mbd_kte/prismatic_joint.hpp"
#include "mbd_kte/rigid_link.hpp"
#include "mbd_kte/damper.hpp"
#include "mbd_kte/spring.hpp"
#include "mbd_kte/torsion_damper.hpp"
#include "mbd_kte/torsion_spring.hpp"
#include "mbd_kte/line_point_mindist.hpp"
#include "mbd_kte/plane_point_mindist.hpp"


#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/sensors/SoTimerSensor.h>


namespace ReaK {

namespace geom {



oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map& aModel) {
  
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


oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map_chain& aModel) {
  
  const std::vector< shared_ptr< kte::kte_map > >& mdl_list = aModel.getKTEs();
  for(std::size_t i = 0; i < mdl_list.size(); ++i)
    aSG << *(mdl_list[i]);
  
  return aSG;
};



};


};





