
/*
 *    Copyright 2013 Sven Mikael Persson
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

#include "oi_reader.hpp"

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
#include <Inventor/nodes/SoRotationXYZ.h>  // for SoRotationXYZ
#include <Inventor/nodes/SoSeparator.h>  // for SoSeparator
#include <Inventor/nodes/SoSphere.h>    // for SoSphere
#include <Inventor/nodes/SoMatrixTransform.h>  // for SoMatrixTransform
#include <Inventor/nodes/SoTransform.h>  // for SoTransform
#include <Inventor/nodes/SoTransformSeparator.h>  // for SoTransformSeparator
#include <Inventor/nodes/SoTransformation.h>  // for SoTransformation
#include <Inventor/nodes/SoTranslation.h>  // for SoTranslation
#include <Inventor/VRMLnodes/SoVRMLParent.h>  // for SoVRMLParent
#include <Inventor/VRMLnodes/SoVRMLGroup.h>  // for SoVRMLGroup
#include <Inventor/VRMLnodes/SoVRMLTransform.h>  // for SoVRMLTransform
#include <Inventor/VRMLnodes/SoVRMLColor.h>  // for SoVRMLColor
#include <Inventor/VRMLnodes/SoVRMLBox.h>  // for SoVRMLBox
#include <Inventor/VRMLnodes/SoVRMLCylinder.h>  // for SoVRMLCylinder
#include <Inventor/VRMLnodes/SoVRMLSphere.h>  // for SoVRMLSphere

#include <stack>
#include <cmath>                        // for sqrt, M_PI

namespace ReaK {

namespace geom {


double oi_reader::computeCharacteristicLength() {
  using std::sqrt;
  SbViewportRegion dummy_viewport;
  SoGetBoundingBoxAction bb_calc(dummy_viewport);
  bb_calc.apply(mRoot);
  float x,y,z;
  bb_calc.getXfBoundingBox().getSize(x,y,z);
  return sqrt((x * x + y * y + z * z) / 3.0);
};



oi_reader::oi_reader() : mRoot(NULL) {};

oi_reader::~oi_reader() {};

oi_reader& oi_reader::read_file(const std::string& aFileName) {
  SoInput file_in;
  if(!file_in.openFile(aFileName.c_str())) {
    return *this;
  };
  
  if(mRoot != NULL)
    mRoot->unref();
  
  mRoot = SoDB::readAll(&file_in);
  
  file_in.closeFile();
  
  return *this;
};

oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel) {
  if(aSG.mRoot == NULL)
    return aSG;
  
  std::stack< std::pair<SoGroup*,int> > so_group_stack;
  std::stack< shared_ptr< pose_3D<double> > > anchor_stack;
  std::stack< color > color_stack;
  
  so_group_stack.push(std::pair<SoGroup*,int>(aSG.mRoot,0));
  if(aModel.mAnchorList.empty()) {
    anchor_stack.push(shared_ptr< pose_3D<double> >(new pose_3D<double>()));
    aModel.addAnchor(anchor_stack.top());
  } else {
    anchor_stack.push(aModel.mAnchorList[0]);
  };
  color_stack.push(color(1.0,1.0,1.0)); // originally white.
  
  // so a depth-first traversal of the objects in the scene-tree:
  while(!so_group_stack.empty()) {
    SoGroup* current_group = so_group_stack.top().first;
    int current_child = so_group_stack.top().second;
    if(current_child >= current_group->getNumChildren()) {
      so_group_stack.pop();
      anchor_stack.pop();
      color_stack.pop();
      continue;
    };
    SoNode* current_node = current_group->getChild(current_child);
    so_group_stack.top().second = current_child + 1;
    
    if((current_node->getTypeId() == SoGroup::getClassTypeId()) ||
       (current_node->getTypeId() == SoSeparator::getClassTypeId()) ||
       (current_node->getTypeId() == SoTransformSeparator::getClassTypeId()) ||
       (current_node->getTypeId() == SoVRMLGroup::getClassTypeId()) ||
       (current_node->getTypeId() == SoVRMLTransform::getClassTypeId())) {
      
      so_group_stack.push(std::pair<SoGroup*,int>(static_cast<SoGroup*>(current_node),0));
      anchor_stack.push(shared_ptr< pose_3D<double> >(new pose_3D<double>(
        anchor_stack.top(),
        vect<double,3>(0.0,0.0,0.0),
        quaternion<double>()
      )));
      aModel.addAnchor(anchor_stack.top());
      
      color_stack.push(color_stack.top());
      
      if(current_node->getTypeId() == SoVRMLTransform::getClassTypeId()) {
        
        SoVRMLTransform* ptmp = static_cast<SoVRMLTransform*>(current_node);
        // set the transformation:
        SbVec3f rot_axis;
        float rot_angle;
        ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
        anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
          vect<double,3>(ptmp->translation.getValue()[0], ptmp->translation.getValue()[1], ptmp->translation.getValue()[2])
          + vect<double,3>(ptmp->center.getValue()[0], ptmp->center.getValue()[1], ptmp->center.getValue()[2]),
          axis_angle<double>(rot_angle, vect<double,3>(rot_axis[0],rot_axis[1],rot_axis[2])).getQuaternion()
        ));
        anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
          vect<double,3>(-ptmp->center.getValue()[0], -ptmp->center.getValue()[1], -ptmp->center.getValue()[2]),
          quaternion<double>()
        ));
      };
      
    } else if(current_node->getTypeId() == SoBaseColor::getClassTypeId()) {
      
      SoBaseColor* ptmp = static_cast<SoBaseColor*>(current_node);
      color_stack.top() = color(ptmp->rgb[0][0], ptmp->rgb[0][1], ptmp->rgb[0][2]);
      
    } else if(current_node->getTypeId() == SoVRMLColor::getClassTypeId()) {
      
      SoVRMLColor* ptmp = static_cast<SoVRMLColor*>(current_node);
      color_stack.top() = color(ptmp->color[0][0], ptmp->color[0][1], ptmp->color[0][2]);
      
    } else if(current_node->getTypeId() == SoCube::getClassTypeId()) {
      
      SoCube* ptmp = static_cast<SoCube*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new box(ptmp->getName().getString(),
                                          anchor_ptr, *anchor_stack.top(),
                                          vect<double,3>(ptmp->width.getValue(),
                                                         ptmp->height.getValue(),
                                                         ptmp->depth.getValue()))));
      
    } else if(current_node->getTypeId() == SoCylinder::getClassTypeId()) {
      
      SoCylinder* ptmp = static_cast<SoCylinder*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      pose_3D<double> final_pose = *anchor_stack.top();
      final_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, 0.0),
        quaternion<double>::xrot(M_PI * 0.5).getQuaternion()
      ));
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new cylinder(ptmp->getName().getString(),
                                               anchor_ptr, final_pose,
                                               ptmp->height.getValue(), 
                                               ptmp->radius.getValue())));
      
    } else if(current_node->getTypeId() == SoSphere::getClassTypeId()) {
      
      SoSphere* ptmp = static_cast<SoSphere*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new sphere(ptmp->getName().getString(),
                                             anchor_ptr, *anchor_stack.top(),
                                             ptmp->radius.getValue())));
      
    } else if(current_node->getTypeId() == SoVRMLBox::getClassTypeId()) {
      
      SoVRMLBox* ptmp = static_cast<SoVRMLBox*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new box(ptmp->getName().getString(),
                                          anchor_ptr, *anchor_stack.top(),
                                          vect<double,3>(ptmp->size.getValue()[0],
                                                         ptmp->size.getValue()[1],
                                                         ptmp->size.getValue()[2]))));
      
    } else if(current_node->getTypeId() == SoVRMLCylinder::getClassTypeId()) {
      
      SoVRMLCylinder* ptmp = static_cast<SoVRMLCylinder*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      pose_3D<double> final_pose = *anchor_stack.top();
      final_pose.addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, 0.0),
        quaternion<double>::xrot(M_PI * 0.5).getQuaternion()
      ));
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new cylinder(ptmp->getName().getString(),
                                               anchor_ptr, final_pose,
                                               ptmp->height.getValue(), 
                                               ptmp->radius.getValue())));
      
    } else if(current_node->getTypeId() == SoVRMLSphere::getClassTypeId()) {
      
      SoVRMLSphere* ptmp = static_cast<SoVRMLSphere*>(current_node);
      shared_ptr< pose_3D<double> > anchor_ptr;
      if(!anchor_stack.top()->Parent.expired())
        anchor_ptr = anchor_stack.top()->Parent.lock();
      aModel.addElement(color_stack.top(),
        shared_ptr< geometry_3D >(new sphere(ptmp->getName().getString(),
                                             anchor_ptr, *anchor_stack.top(),
                                             ptmp->radius.getValue())));
      
    } else if(current_node->getTypeId() == SoRotation::getClassTypeId()) {
      
      SoRotation* ptmp = static_cast<SoRotation*>(current_node);
      SbVec3f rot_axis;
      float rot_angle;
      ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(0.0, 0.0, 0.0),
        axis_angle<double>(rot_angle, vect<double,3>(rot_axis[0],rot_axis[1],rot_axis[2])).getQuaternion()
      ));
      
    } else if(current_node->getTypeId() == SoRotationXYZ::getClassTypeId()) {
      
      SoRotationXYZ* ptmp = static_cast<SoRotationXYZ*>(current_node);
      quaternion<double> rot;
      if(ptmp->axis.getValue() == SoRotationXYZ::X)
        rot = quaternion<double>::xrot(ptmp->angle.getValue()).getQuaternion();
      else if(ptmp->axis.getValue() == SoRotationXYZ::Y)
        rot = quaternion<double>::yrot(ptmp->angle.getValue()).getQuaternion();
      else
        rot = quaternion<double>::zrot(ptmp->angle.getValue()).getQuaternion();
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),vect<double,3>(0.0, 0.0, 0.0), rot));
      
    } else if(current_node->getTypeId() == SoTransform::getClassTypeId()) {
      
      SoTransform* ptmp = static_cast<SoTransform*>(current_node);
      // set the transformation:
      SbVec3f rot_axis;
      float rot_angle;
      ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(ptmp->translation.getValue()[0], ptmp->translation.getValue()[1], ptmp->translation.getValue()[2])
        + vect<double,3>(ptmp->center.getValue()[0], ptmp->center.getValue()[1], ptmp->center.getValue()[2]),
        axis_angle<double>(rot_angle, vect<double,3>(rot_axis[0],rot_axis[1],rot_axis[2])).getQuaternion()
      ));
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(-ptmp->center.getValue()[0], -ptmp->center.getValue()[1], -ptmp->center.getValue()[2]),
        quaternion<double>()
      ));
      
    } else if(current_node->getTypeId() == SoMatrixTransform::getClassTypeId()) {
      
      SoMatrixTransform* ptmp = static_cast<SoMatrixTransform*>(current_node);
      SbVec3f translation;
      SbRotation rotation;
      SbVec3f scaleFactor;
      SbRotation scaleOrientation;
      SbVec3f center(0.0,0.0,0.0);
      ptmp->matrix.getValue().getTransform(translation,rotation,scaleFactor,scaleOrientation,center);
      // set the transformation:
      SbVec3f rot_axis;
      float rot_angle;
      rotation.getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(translation[0], translation[1], translation[2]),
        axis_angle<double>(rot_angle, vect<double,3>(rot_axis[0],rot_axis[1],rot_axis[2])).getQuaternion()
      ));
      
    } else if(current_node->getTypeId() == SoTranslation::getClassTypeId()) {
      
      SoTranslation* ptmp = static_cast<SoTranslation*>(current_node);
      anchor_stack.top()->addBefore(pose_3D<double>(weak_ptr< pose_3D<double> >(),
        vect<double,3>(ptmp->translation.getValue()[0], ptmp->translation.getValue()[1], ptmp->translation.getValue()[2]),
        quaternion<double>()
      ));
      
    }; // else not supported
    
  };
  
  return aSG;
};



};


};





