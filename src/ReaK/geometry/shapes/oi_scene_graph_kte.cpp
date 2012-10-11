
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
#include "mbd_kte/free_joints.hpp"
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
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/sensors/SoTimerSensor.h>


namespace ReaK {

namespace geom {


static const float kte_gray_R = 0.8196;
static const float kte_gray_G = 0.8549;
static const float kte_gray_B = 0.8980;

static const float kte_red_R = 0.9922;
static const float kte_red_G = 0.0353;
static const float kte_red_B = 0.1922;


oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map& aModel) {
  
  // first check if the model is a kte-chain:
  const void* p_chain = aModel.castTo(kte::kte_map_chain::getStaticObjectType());
  if(p_chain)
    return aSG << static_cast<const kte::kte_map_chain&>(aModel);
  
  
  if(aModel.castTo(kte::revolute_joint_3D::getStaticObjectType())) {
    const kte::revolute_joint_3D& rev_joint = static_cast<const kte::revolute_joint_3D&>(aModel);
    
    SoSeparator* sep = new SoSeparator;
    
    vect<double,3> x_axis = rev_joint.Axis() % vect<double,3>(0.0,0.0,1.0);
    double x_axis_mag = norm_2(x_axis);
    if(x_axis_mag < 1e-5)
      x_axis = vect<double,3>(rev_joint.Axis()[2],0.0,0.0);
    else
      x_axis /= x_axis_mag;
    vect<double,3> z_axis = x_axis % rev_joint.Axis();
    SoTransform* trans = new SoTransform;
    SbMatrix rot_mat(
      x_axis[0], rev_joint.Axis()[0], z_axis[0], 0.0,
      x_axis[1], rev_joint.Axis()[1], z_axis[1], 0.0,
      x_axis[2], rev_joint.Axis()[2], z_axis[2], 0.0,
      0.0, 0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);
    
    SoCylinder* core_cyl = new SoCylinder;
    core_cyl->radius = aSG.mCharacteristicLength * 0.01;  // diameter about 2% of char-len
    core_cyl->height = aSG.mCharacteristicLength * 0.04;    // about 4% of char-len
    sep->addChild(core_cyl);
    
    SoBaseColor * white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);
    
    SoTexture2* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/red_arrow.tga");
    sep->addChild(red_arrow);
    
    SoCylinder* arrow_cyl = new SoCylinder;
    arrow_cyl->parts = SoCylinder::SIDES;
    arrow_cyl->radius = aSG.mCharacteristicLength * 0.015;
    arrow_cyl->height = aSG.mCharacteristicLength * 0.04;
    sep->addChild(arrow_cyl);
    
    if(!rev_joint.BaseFrame()) {
      aSG.mRoot->addChild(sep);
    } else {
      if(aSG.mAnchor3DMap.find(rev_joint.BaseFrame()) == aSG.mAnchor3DMap.end())
        aSG << rev_joint.BaseFrame();
      aSG.mAnchor3DMap[rev_joint.BaseFrame()].first->addChild(sep);
    };
    
  } else if(aModel.castTo(kte::prismatic_joint_3D::getStaticObjectType())) {
    const kte::prismatic_joint_3D& pri_joint = static_cast<const kte::prismatic_joint_3D&>(aModel);
    
    SoSeparator* sep = new SoSeparator;
    
    vect<double,3> x_axis = pri_joint.Axis() % vect<double,3>(0.0,0.0,1.0);
    double x_axis_mag = norm_2(x_axis);
    if(x_axis_mag < 1e-5)
      x_axis = vect<double,3>(pri_joint.Axis()[2],0.0,0.0);
    else
      x_axis /= x_axis_mag;
    vect<double,3> z_axis = x_axis % pri_joint.Axis();
    SoTransform* trans = new SoTransform;
    SbMatrix rot_mat(
      x_axis[0], pri_joint.Axis()[0], z_axis[0], 0.0,
      x_axis[1], pri_joint.Axis()[1], z_axis[1], 0.0,
      x_axis[2], pri_joint.Axis()[2], z_axis[2], 0.0,
      0.0, 0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);
    
    SoCube* core_cube = new SoCube;
    core_cube->width = aSG.mCharacteristicLength * 0.01;  // about 2% of char-len
    core_cube->height = aSG.mCharacteristicLength * 0.04;    // about 4% of char-len
    core_cube->depth = aSG.mCharacteristicLength * 0.01;  // about 2% of char-len
    sep->addChild(core_cube);
    
    SoBaseColor * white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);
    
    SoTexture2* red_arrows = new SoTexture2;
    red_arrows->filename.setValue("./models/textures/red_up_arrows.tga");
    sep->addChild(red_arrows);
    
    SoCylinder* arrow_cyl = new SoCylinder;
    arrow_cyl->parts = SoCylinder::SIDES;
    arrow_cyl->radius = aSG.mCharacteristicLength * 0.01;
    arrow_cyl->height = aSG.mCharacteristicLength * 0.04;
    sep->addChild(arrow_cyl);
    
    if(!pri_joint.BaseFrame()) {
      aSG.mRoot->addChild(sep);
    } else {
      if(aSG.mAnchor3DMap.find(pri_joint.BaseFrame()) == aSG.mAnchor3DMap.end())
        aSG << pri_joint.BaseFrame();
      aSG.mAnchor3DMap[pri_joint.BaseFrame()].first->addChild(sep);
    };
    
  } else if(aModel.castTo(kte::rigid_link_3D::getStaticObjectType())) {
    const kte::rigid_link_3D& lnk_obj = static_cast<const kte::rigid_link_3D&>(aModel);
    
    vect<double,3> y_axis = lnk_obj.PoseOffset().Position;
    double y_axis_mag = norm_2(y_axis);
    if(y_axis_mag < 1e-5)   // nothing to draw... 
      return aSG;
    else
      y_axis /= y_axis_mag;
    
    SoSeparator* sep = new SoSeparator;
    
    vect<double,3> x_axis = y_axis % vect<double,3>(0.0,0.0,1.0);
    double x_axis_mag = norm_2(x_axis);
    if(x_axis_mag < 1e-5)
      x_axis = vect<double,3>(y_axis[2],0.0,0.0);
    else
      x_axis /= x_axis_mag;
    vect<double,3> z_axis = x_axis % y_axis;
    SoTransform* trans = new SoTransform;
    trans->translation.setValue(
      lnk_obj.PoseOffset().Position[0] * 0.5,
      lnk_obj.PoseOffset().Position[1] * 0.5,
      lnk_obj.PoseOffset().Position[2] * 0.5);
    SbMatrix rot_mat(
      x_axis[0], y_axis[0], z_axis[0], 0.0,
      x_axis[1], y_axis[1], z_axis[1], 0.0,
      x_axis[2], y_axis[2], z_axis[2], 0.0,
      0.0, 0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);
    
    SoBaseColor * col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);
    
    SoCylinder* core_cyl = new SoCylinder;
    core_cyl->radius = aSG.mCharacteristicLength * 0.005;  // diameter about 1% of char-len (very thin rod).
    core_cyl->height = y_axis_mag;
    sep->addChild(core_cyl);
    
    if(!lnk_obj.BaseFrame()) {
      aSG.mRoot->addChild(sep);
    } else {
      if(aSG.mAnchor3DMap.find(lnk_obj.BaseFrame()) == aSG.mAnchor3DMap.end())
        aSG << lnk_obj.BaseFrame();
      aSG.mAnchor3DMap[lnk_obj.BaseFrame()].first->addChild(sep);
    };
    
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





