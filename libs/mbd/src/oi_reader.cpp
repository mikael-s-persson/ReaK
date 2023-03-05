
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

#include "ReaK/mbd/coin3D/oi_reader.hpp"

#include "ReaK/geometry/shapes/geometry_3D.hpp"
#include "ReaK/geometry/shapes/shape_3D.hpp"

#include "ReaK/geometry/shapes/box.hpp"
#include "ReaK/geometry/shapes/composite_shape_3D.hpp"
#include "ReaK/geometry/shapes/cylinder.hpp"
#include "ReaK/geometry/shapes/plane.hpp"
#include "ReaK/geometry/shapes/sphere.hpp"

#include "ReaK/geometry/proximity/proxy_query_model.hpp"

#include <Inventor/SbColor.h>                         // for SbColor
#include <Inventor/SbVec3f.h>                         // for SbVec3f
#include <Inventor/SbViewportRegion.h>                // for SbViewportRegion
#include <Inventor/SbXfBox3f.h>                       // for SbXfBox3f
#include <Inventor/VRMLnodes/SoVRMLBox.h>             // for SoVRMLBox
#include <Inventor/VRMLnodes/SoVRMLColor.h>           // for SoVRMLColor
#include <Inventor/VRMLnodes/SoVRMLCylinder.h>        // for SoVRMLCylinder
#include <Inventor/VRMLnodes/SoVRMLGroup.h>           // for SoVRMLGroup
#include <Inventor/VRMLnodes/SoVRMLParent.h>          // for SoVRMLParent
#include <Inventor/VRMLnodes/SoVRMLSphere.h>          // for SoVRMLSphere
#include <Inventor/VRMLnodes/SoVRMLTransform.h>       // for SoVRMLTransform
#include <Inventor/actions/SoGetBoundingBoxAction.h>  // for SoGetBoundingBoxAction
#include <Inventor/fields/SoMFColor.h>                // for SoMFColor
#include <Inventor/fields/SoMFInt32.h>                // for SoMFInt32
#include <Inventor/fields/SoMFVec3f.h>                // for SoMFVec3f
#include <Inventor/fields/SoSFFloat.h>                // for SoSFFloat
#include <Inventor/fields/SoSFRotation.h>             // for SoSFRotation
#include <Inventor/fields/SoSFVec3f.h>                // for SoSFVec3f
#include <Inventor/fields/SoSubField.h>        // for SoSFFloat::operator=, etc
#include <Inventor/misc/SoChildList.h>         // for SoChildList
#include <Inventor/nodes/SoBaseColor.h>        // for SoBaseColor
#include <Inventor/nodes/SoCube.h>             // for SoCube
#include <Inventor/nodes/SoCylinder.h>         // for SoCylinder
#include <Inventor/nodes/SoMatrixTransform.h>  // for SoMatrixTransform
#include <Inventor/nodes/SoRotation.h>         // for SoRotation
#include <Inventor/nodes/SoRotationXYZ.h>      // for SoRotationXYZ
#include <Inventor/nodes/SoSeparator.h>        // for SoSeparator
#include <Inventor/nodes/SoSphere.h>           // for SoSphere
#include <Inventor/nodes/SoTransform.h>        // for SoTransform
#include <Inventor/nodes/SoTransformSeparator.h>  // for SoTransformSeparator
#include <Inventor/nodes/SoTransformation.h>      // for SoTransformation
#include <Inventor/nodes/SoTranslation.h>         // for SoTranslation
#include <cmath>

#include <memory>
#include <stack>

namespace ReaK::geom {

double oi_reader::computeCharacteristicLength() {
  if (mRoot == nullptr) {
    return 0.0;
  }
  using std::sqrt;
  SbViewportRegion dummy_viewport;
  SoGetBoundingBoxAction bb_calc(dummy_viewport);
  bb_calc.apply(mRoot);
  float x = NAN;
  float y = NAN;
  float z = NAN;
  bb_calc.getXfBoundingBox().getSize(x, y, z);
  return sqrt((x * x + y * y + z * z) / 3.0);
}

oi_reader::oi_reader() : mRoot(nullptr) {}

oi_reader::oi_reader(const std::string& aFileName) : mRoot(nullptr) {
  read_file(aFileName);
}

oi_reader::oi_reader(std::istream& aStream) : mRoot(nullptr) {
  read_stream(aStream);
}

oi_reader::oi_reader(const oi_reader& rhs) : mRoot(rhs.mRoot) {
  if (mRoot != nullptr) {
    mRoot->ref();
  }
}

oi_reader& oi_reader::operator=(const oi_reader& rhs) {
  if (mRoot != nullptr) {
    mRoot->unref();
  }
  mRoot = rhs.mRoot;
  if (mRoot != nullptr) {
    mRoot->ref();
  }
  return *this;
}

oi_reader::oi_reader(oi_reader&& rhs) noexcept : mRoot(rhs.mRoot) {
  rhs.mRoot = nullptr;
}

oi_reader& oi_reader::operator=(oi_reader&& rhs) noexcept {
  if (mRoot != nullptr) {
    mRoot->unref();
  }
  mRoot = rhs.mRoot;
  rhs.mRoot = nullptr;
  return *this;
}

oi_reader::~oi_reader() {
  if (mRoot != nullptr) {
    mRoot->unref();
  }
}

oi_reader& oi_reader::read_file(const std::string& aFileName) {
  SoInput file_in;
  if (file_in.openFile(aFileName.c_str()) == 0) {
    return *this;
  }

  if (mRoot != nullptr) {
    mRoot->unref();
  }

  mRoot = SoDB::readAll(&file_in);
  if (mRoot != nullptr) {
    mRoot->ref();
  }

  file_in.closeFile();

  return *this;
}

oi_reader& oi_reader::read_stream(std::istream& aStream) {
  if (!aStream) {
    return *this;
  }
  aStream.seekg(0, std::istream::end);
  int length = aStream.tellg();
  aStream.seekg(0, std::istream::beg);

  char* buffer = new char[length];
  aStream.read(buffer, length);

  SoInput file_in;
  file_in.setBuffer(buffer, length);

  if (mRoot != nullptr) {
    mRoot->unref();
  }

  mRoot = SoDB::readAll(&file_in);
  if (mRoot != nullptr) {
    mRoot->ref();
  }

  delete[] buffer;

  return *this;
}

static void read_sg_into_models(SoSeparator* aRoot, colored_model_3D* aModel,
                                proxy_query_model_3D* aProxy) {

  std::stack<std::pair<SoChildList*, int>> so_group_stack;
  std::stack<std::shared_ptr<pose_3D<double>>> anchor_stack;
  std::stack<color> color_stack;

  so_group_stack.push(std::pair<SoChildList*, int>(aRoot->getChildren(), 0));
  if (aModel != nullptr) {
    if (aModel->mAnchorList.empty()) {
      anchor_stack.push(std::make_shared<pose_3D<double>>());
      aModel->addAnchor(anchor_stack.top());
    } else {
      anchor_stack.push(aModel->mAnchorList[0]);
    }
  } else {
    anchor_stack.push(std::make_shared<pose_3D<double>>());
  }
  color_stack.push(color(1.0, 1.0, 1.0));  // originally white.

  // so a depth-first traversal of the objects in the scene-tree:
  while (!so_group_stack.empty()) {
    SoChildList* current_group = so_group_stack.top().first;
    int current_child = so_group_stack.top().second;
    if (current_child >= current_group->getLength()) {
      so_group_stack.pop();
      anchor_stack.pop();
      color_stack.pop();
      continue;
    }
    SoNode* current_node = (*current_group)[current_child];
    so_group_stack.top().second = current_child + 1;

    if ((current_node->getTypeId() == SoBaseColor::getClassTypeId()) != 0) {

      auto* ptmp = static_cast<SoBaseColor*>(current_node);
      color_stack.top() =
          color(ptmp->rgb[0][0], ptmp->rgb[0][1], ptmp->rgb[0][2]);

    } else if ((current_node->getTypeId() == SoVRMLColor::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoVRMLColor*>(current_node);
      color_stack.top() =
          color(ptmp->color[0][0], ptmp->color[0][1], ptmp->color[0][2]);

    } else if ((current_node->getTypeId() == SoCube::getClassTypeId()) != 0) {

      auto* ptmp = static_cast<SoCube*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      std::shared_ptr<box> new_box(new box(
          ptmp->getName().getString(), anchor_ptr, *anchor_stack.top(),
          vect<double, 3>(ptmp->width.getValue(), ptmp->height.getValue(),
                          ptmp->depth.getValue())));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_box);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_box);
      }

    } else if ((current_node->getTypeId() == SoCylinder::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoCylinder*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      pose_3D<double> final_pose = *anchor_stack.top();
      final_pose.addBefore(pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 0.0, 0.0),
          quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
      std::shared_ptr<cylinder> new_cyl(
          new cylinder(ptmp->getName().getString(), anchor_ptr, final_pose,
                       ptmp->height.getValue(), ptmp->radius.getValue()));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_cyl);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_cyl);
      }

    } else if ((current_node->getTypeId() == SoSphere::getClassTypeId()) != 0) {

      auto* ptmp = static_cast<SoSphere*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      std::shared_ptr<sphere> new_sph(
          new sphere(ptmp->getName().getString(), anchor_ptr,
                     *anchor_stack.top(), ptmp->radius.getValue()));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_sph);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_sph);
      }

    } else if ((current_node->getTypeId() == SoVRMLBox::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoVRMLBox*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      std::shared_ptr<box> new_box(new box(
          ptmp->getName().getString(), anchor_ptr, *anchor_stack.top(),
          vect<double, 3>(ptmp->size.getValue()[0], ptmp->size.getValue()[1],
                          ptmp->size.getValue()[2])));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_box);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_box);
      }

    } else if ((current_node->getTypeId() ==
                SoVRMLCylinder::getClassTypeId()) != 0) {

      auto* ptmp = static_cast<SoVRMLCylinder*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      pose_3D<double> final_pose = *anchor_stack.top();
      final_pose.addBefore(pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 0.0, 0.0),
          quaternion<double>::xrot(M_PI * 0.5).getQuaternion()));
      std::shared_ptr<cylinder> new_cyl(
          new cylinder(ptmp->getName().getString(), anchor_ptr, final_pose,
                       ptmp->height.getValue(), ptmp->radius.getValue()));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_cyl);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_cyl);
      }

    } else if ((current_node->getTypeId() == SoVRMLSphere::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoVRMLSphere*>(current_node);
      std::shared_ptr<pose_3D<double>> anchor_ptr;
      if (!anchor_stack.top()->Parent.expired()) {
        anchor_ptr = anchor_stack.top()->Parent.lock();
      }
      std::shared_ptr<sphere> new_sph(
          new sphere(ptmp->getName().getString(), anchor_ptr,
                     *anchor_stack.top(), ptmp->radius.getValue()));
      if (aModel != nullptr) {
        aModel->addElement(color_stack.top(), new_sph);
      }
      if (aProxy != nullptr) {
        aProxy->addShape(new_sph);
      }

    } else if ((current_node->getTypeId() == SoRotation::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoRotation*>(current_node);
      SbVec3f rot_axis;
      float rot_angle = NAN;
      ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 0.0, 0.0),
          axis_angle<double>(
              rot_angle, vect<double, 3>(rot_axis[0], rot_axis[1], rot_axis[2]))
              .getQuaternion()));

    } else if ((current_node->getTypeId() == SoRotationXYZ::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoRotationXYZ*>(current_node);
      quaternion<double> rot;
      if (ptmp->axis.getValue() == SoRotationXYZ::X) {
        rot = quaternion<double>::xrot(ptmp->angle.getValue()).getQuaternion();
      } else if (ptmp->axis.getValue() == SoRotationXYZ::Y) {
        rot = quaternion<double>::yrot(ptmp->angle.getValue()).getQuaternion();
      } else {
        rot = quaternion<double>::zrot(ptmp->angle.getValue()).getQuaternion();
      }
      anchor_stack.top()->addBefore(
          pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                          vect<double, 3>(0.0, 0.0, 0.0), rot));

    } else if ((current_node->getTypeId() == SoTransform::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoTransform*>(current_node);
      // set the transformation:
      SbVec3f rot_axis;
      float rot_angle = NAN;
      ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(),
          vect<double, 3>(ptmp->translation.getValue()[0],
                          ptmp->translation.getValue()[1],
                          ptmp->translation.getValue()[2]) +
              vect<double, 3>(ptmp->center.getValue()[0],
                              ptmp->center.getValue()[1],
                              ptmp->center.getValue()[2]),
          axis_angle<double>(
              rot_angle, vect<double, 3>(rot_axis[0], rot_axis[1], rot_axis[2]))
              .getQuaternion()));
      anchor_stack.top()->addBefore(
          pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                          vect<double, 3>(-ptmp->center.getValue()[0],
                                          -ptmp->center.getValue()[1],
                                          -ptmp->center.getValue()[2]),
                          quaternion<double>()));

    } else if ((current_node->getTypeId() ==
                SoMatrixTransform::getClassTypeId()) != 0) {

      auto* ptmp = static_cast<SoMatrixTransform*>(current_node);
      SbVec3f translation;
      SbRotation rotation;
      SbVec3f scaleFactor;
      SbRotation scaleOrientation;
      ptmp->matrix.getValue().getTransform(translation, rotation, scaleFactor,
                                           scaleOrientation);
      // set the transformation:
      SbVec3f rot_axis;
      float rot_angle = NAN;
      rotation.getValue(rot_axis, rot_angle);
      anchor_stack.top()->addBefore(pose_3D<double>(
          std::weak_ptr<pose_3D<double>>(),
          vect<double, 3>(translation[0], translation[1], translation[2]),
          axis_angle<double>(
              rot_angle, vect<double, 3>(rot_axis[0], rot_axis[1], rot_axis[2]))
              .getQuaternion()));

    } else if ((current_node->getTypeId() == SoTranslation::getClassTypeId()) !=
               0) {

      auto* ptmp = static_cast<SoTranslation*>(current_node);
      anchor_stack.top()->addBefore(
          pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                          vect<double, 3>(ptmp->translation.getValue()[0],
                                          ptmp->translation.getValue()[1],
                                          ptmp->translation.getValue()[2]),
                          quaternion<double>()));

    } else if (current_node->getTypeId().isDerivedFrom(
                   SoShape::getClassTypeId()) != 0) {

      // just register the bounding-box, cannot deal with the display.
      SbViewportRegion dummy_viewport;
      SoGetBoundingBoxAction bb_calc(dummy_viewport);
      bb_calc.apply(current_node);
      SbXfBox3f& bbox = bb_calc.getXfBoundingBox();
      if ((bbox.getVolume() > std::numeric_limits<float>::epsilon()) &&
          ((aProxy) != nullptr)) {
        float x = NAN;
        float y = NAN;
        float z = NAN;
        bbox.getSize(x, y, z);
        SbVec3f center = bbox.getCenter();

        SbVec3f translation;
        SbRotation rotation;
        SbVec3f scaleFactor;
        SbRotation scaleOrientation;
        bbox.getTransform().getTransform(translation, rotation, scaleFactor,
                                         scaleOrientation);
        // set the transformation:
        SbVec3f rot_axis;
        float rot_angle = NAN;
        rotation.getValue(rot_axis, rot_angle);

        std::shared_ptr<pose_3D<double>> anchor_ptr;
        if (!anchor_stack.top()->Parent.expired()) {
          anchor_ptr = anchor_stack.top()->Parent.lock();
        }
        pose_3D<double> final_pose = *anchor_stack.top();
        final_pose.addBefore(pose_3D<double>(
            std::weak_ptr<pose_3D<double>>(),
            vect<double, 3>(translation[0] + center[0],
                            translation[1] + center[1],
                            translation[2] + center[2]),
            axis_angle<double>(
                rot_angle,
                vect<double, 3>(rot_axis[0], rot_axis[1], rot_axis[2]))
                .getQuaternion()));

        std::shared_ptr<box> new_box(
            new box(current_node->getName().getString(), anchor_ptr, final_pose,
                    vect<double, 3>(x, y, z)));
        aProxy->addShape(new_box);
      }
    }

    SoChildList* next_childlist = current_node->getChildren();
    if ((next_childlist != nullptr) && (next_childlist->getLength() != 0)) {
      so_group_stack.push(std::pair<SoChildList*, int>(next_childlist, 0));
      anchor_stack.push(std::make_shared<pose_3D<double>>(
          anchor_stack.top(), vect<double, 3>(0.0, 0.0, 0.0),
          quaternion<double>()));
      if (aModel != nullptr) {
        aModel->addAnchor(anchor_stack.top());
      }

      color_stack.push(color_stack.top());

      if ((current_node->getTypeId() == SoVRMLTransform::getClassTypeId()) !=
          0) {

        auto* ptmp = static_cast<SoVRMLTransform*>(current_node);
        // set the transformation:
        SbVec3f rot_axis;
        float rot_angle = NAN;
        ptmp->rotation.getValue().getValue(rot_axis, rot_angle);
        anchor_stack.top()->addBefore(pose_3D<double>(
            std::weak_ptr<pose_3D<double>>(),
            vect<double, 3>(ptmp->translation.getValue()[0],
                            ptmp->translation.getValue()[1],
                            ptmp->translation.getValue()[2]) +
                vect<double, 3>(ptmp->center.getValue()[0],
                                ptmp->center.getValue()[1],
                                ptmp->center.getValue()[2]),
            axis_angle<double>(
                rot_angle,
                vect<double, 3>(rot_axis[0], rot_axis[1], rot_axis[2]))
                .getQuaternion()));
        anchor_stack.top()->addBefore(
            pose_3D<double>(std::weak_ptr<pose_3D<double>>(),
                            vect<double, 3>(-ptmp->center.getValue()[0],
                                            -ptmp->center.getValue()[1],
                                            -ptmp->center.getValue()[2]),
                            quaternion<double>()));
      }
    }
  }
}

oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel) {
  if (aSG.mRoot == nullptr) {
    return aSG;
  }

  read_sg_into_models(aSG.mRoot, &aModel, nullptr);

  return aSG;
}

oi_reader& operator>>(oi_reader& aSG, proxy_query_model_3D& aProxy) {
  if (aSG.mRoot == nullptr) {
    return aSG;
  }

  read_sg_into_models(aSG.mRoot, nullptr, &aProxy);

  return aSG;
}

oi_reader& oi_reader::translate_into(colored_model_3D& aModel,
                                     proxy_query_model_3D& aProxy) {
  if (mRoot == nullptr) {
    return *this;
  }

  read_sg_into_models(mRoot, &aModel, &aProxy);

  return *this;
}

}  // namespace ReaK::geom
