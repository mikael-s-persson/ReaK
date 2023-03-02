
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

#include <ReaK/mbd/coin3D/oi_scene_graph.hpp>

#include <ReaK/mbd/kte/damper.hpp>       // done.
#include <ReaK/mbd/kte/free_joints.hpp>  // done.
#include <ReaK/mbd/kte/inertia.hpp>      // done.
#include <ReaK/mbd/kte/kte_map.hpp>      // for kte_map
#include <ReaK/mbd/kte/kte_map_chain.hpp>
#include <ReaK/mbd/kte/prismatic_joint.hpp>  // done.
#include <ReaK/mbd/kte/revolute_joint.hpp>   // done.
#include <ReaK/mbd/kte/rigid_link.hpp>       // done.
#include <ReaK/mbd/kte/spring.hpp>           // done.
#include <ReaK/mbd/kte/torsion_damper.hpp>   // done.
#include <ReaK/mbd/kte/torsion_spring.hpp>   // done.
#include <utility>
//#include <ReaK/mbd/kte/line_point_mindist.hpp>
//#include <ReaK/mbd/kte/plane_point_mindist.hpp>

#include <Inventor/SbColor.h>              // for SbColor
#include <Inventor/SbMatrix.h>             // for SbMatrix
#include <Inventor/SbRotation.h>           // for SbRotation
#include <Inventor/SbVec3f.h>              // for SbVec3f
#include <Inventor/fields/SoMFColor.h>     // for SoMFColor
#include <Inventor/fields/SoMFInt32.h>     // for SoMFInt32
#include <Inventor/fields/SoMFVec2f.h>     // for SoMFVec2f
#include <Inventor/fields/SoMFVec3f.h>     // for SoMFVec3f
#include <Inventor/fields/SoSFBitMask.h>   // for SoSFBitMask
#include <Inventor/fields/SoSFFloat.h>     // for SoSFFloat
#include <Inventor/fields/SoSFRotation.h>  // for SoSFRotation
#include <Inventor/fields/SoSFString.h>    // for SoSFString
#include <Inventor/fields/SoSFVec3f.h>     // for SoSFVec3f
#include <Inventor/fields/SoSubField.h>    // for SoMFColor::operator=, etc
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoCube.h>
#include <Inventor/nodes/SoCylinder.h>
#include <Inventor/nodes/SoFaceSet.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoScale.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoSphere.h>
#include <Inventor/nodes/SoSwitch.h>
#include <Inventor/nodes/SoTexture2.h>
#include <Inventor/nodes/SoTextureCoordinate2.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoTranslation.h>

namespace ReaK::geom {

static const float kte_gray_R = 0.8196;
static const float kte_gray_G = 0.8549;
static const float kte_gray_B = 0.8980;

static const float kte_red_R = 0.9922;
static const float kte_red_G = 0.0353;
static const float kte_red_B = 0.1922;

namespace detail {

struct update_inter_frame_2D {
  std::shared_ptr<pose_2D<double>> Frame1;
  std::shared_ptr<pose_2D<double>> Frame2;
  SoTransform* pTrans;
  SoScale* pScale;

  update_inter_frame_2D(std::shared_ptr<pose_2D<double>> aFrame1,
                        std::shared_ptr<pose_2D<double>> aFrame2,
                        SoTransform* aTrans, SoScale* aScale)
      : Frame1(std::move(aFrame1)),
        Frame2(std::move(aFrame2)),
        pTrans(aTrans),
        pScale(aScale) {}

  void operator()() const {
    using std::atan2;

    vect<double, 2> p1 = Frame1->getGlobalPose().Position;
    vect<double, 2> p2 = Frame2->getGlobalPose().Position;

    vect<double, 2> midpoint = (p1 + p2) * 0.5;
    vect<double, 2> dp = (p2 - p1);

    pTrans->translation.setValue(midpoint[0], midpoint[1], 0.0);
    pTrans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0), atan2(-dp[0], dp[1]));

    pScale->scaleFactor.setValue(1.0, norm_2(dp) * 0.5, 1.0);
  }
};

struct update_inter_frame_3D {
  std::shared_ptr<pose_3D<double>> Frame1;
  std::shared_ptr<pose_3D<double>> Frame2;
  SoTransform* pTrans;
  SoScale* pScale;

  update_inter_frame_3D(std::shared_ptr<pose_3D<double>> aFrame1,
                        std::shared_ptr<pose_3D<double>> aFrame2,
                        SoTransform* aTrans, SoScale* aScale)
      : Frame1(std::move(aFrame1)),
        Frame2(std::move(aFrame2)),
        pTrans(aTrans),
        pScale(aScale) {}

  void operator()() const {
    using std::atan2;

    vect<double, 3> p1 = Frame1->getGlobalPose().Position;
    vect<double, 3> p2 = Frame2->getGlobalPose().Position;

    vect<double, 3> midpoint = (p1 + p2) * 0.5;
    vect<double, 3> dp = (p2 - p1);
    double dp_norm = norm_2(dp);
    if (dp_norm < 1e-5) {
      return;
    }
    dp /= dp_norm;

    pTrans->translation.setValue(midpoint[0], midpoint[1], midpoint[2]);

    vect<double, 3> x_axis = dp % vect<double, 3>(0.0, 0.0, 1.0);
    double x_axis_mag = norm_2(x_axis);
    if (x_axis_mag < 1e-5) {
      x_axis = vect<double, 3>(dp[2], 0.0, 0.0);
    } else {
      x_axis /= x_axis_mag;
    }
    vect<double, 3> z_axis = x_axis % dp;
    SbMatrix rot_mat(x_axis[0], x_axis[1], x_axis[2], 0.0, dp[0], dp[1], dp[2],
                     0.0, z_axis[0], z_axis[1], z_axis[2], 0.0, 0.0, 0.0, 0.0,
                     1.0);
    // SbMatrix rot_mat(
    //  x_axis[0], dp[0], z_axis[0], 0.0,
    //  x_axis[1], dp[1], z_axis[1], 0.0,
    //  x_axis[2], dp[2], z_axis[2], 0.0,
    //  0.0, 0.0, 0.0, 1.0);
    pTrans->rotation.setValue(SbRotation(rot_mat));

    pScale->scaleFactor.setValue(1.0, dp_norm * 0.5, 1.0);
  }
};
}  // namespace detail

oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map& aModel) {

  // first check if the model is a kte-chain:
  const void* p_chain =
      aModel.castTo(kte::kte_map_chain::getStaticObjectType());
  if (p_chain != nullptr) {
    return aSG << static_cast<const kte::kte_map_chain&>(aModel);
  }

  if (aModel.castTo(kte::revolute_joint_3D::getStaticObjectType()) != nullptr) {
    const auto& rev_joint = static_cast<const kte::revolute_joint_3D&>(aModel);

    auto* sep = new SoSeparator;

    vect<double, 3> x_axis = rev_joint.Axis() % vect<double, 3>(0.0, 0.0, 1.0);
    double x_axis_mag = norm_2(x_axis);
    if (x_axis_mag < 1e-5) {
      x_axis = vect<double, 3>(rev_joint.Axis()[2], 0.0, 0.0);
    } else {
      x_axis /= x_axis_mag;
    }
    vect<double, 3> z_axis = x_axis % rev_joint.Axis();
    auto* trans = new SoTransform;
    SbMatrix rot_mat(x_axis[0], x_axis[1], x_axis[2], 0.0, rev_joint.Axis()[0],
                     rev_joint.Axis()[1], rev_joint.Axis()[2], 0.0, z_axis[0],
                     z_axis[1], z_axis[2], 0.0, 0.0, 0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);

    auto* col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);

    auto* core_cyl = new SoCylinder;
    core_cyl->radius =
        aSG.mCharacteristicLength * 0.01;  // diameter about 2% of char-len
    core_cyl->height =
        aSG.mCharacteristicLength * 0.04;  // about 4% of char-len
    sep->addChild(core_cyl);

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/red_arrow.tga");
    sep->addChild(red_arrow);

    auto* arrow_cyl = new SoCylinder;
    arrow_cyl->parts = SoCylinder::SIDES;
    arrow_cyl->radius = aSG.mCharacteristicLength * 0.015;
    arrow_cyl->height = aSG.mCharacteristicLength * 0.04;
    sep->addChild(arrow_cyl);

    if (!rev_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(rev_joint.BaseFrame()) ==
          aSG.mAnchor3DMap.end()) {
        aSG << rev_joint.BaseFrame();
      }
      aSG.mAnchor3DMap[rev_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::revolute_joint_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& rev_joint = static_cast<const kte::revolute_joint_2D&>(aModel);

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/rev_joint_2D.tga");
    sep->addChild(red_arrow);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!rev_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(rev_joint.BaseFrame()) ==
          aSG.mAnchor2DMap.end()) {
        aSG << rev_joint.BaseFrame();
      }
      aSG.mAnchor2DMap[rev_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::prismatic_joint_3D::getStaticObjectType()) !=
             nullptr) {
    const auto& pri_joint = static_cast<const kte::prismatic_joint_3D&>(aModel);

    auto* sep = new SoSeparator;

    vect<double, 3> x_axis = pri_joint.Axis() % vect<double, 3>(0.0, 0.0, 1.0);
    double x_axis_mag = norm_2(x_axis);
    if (x_axis_mag < 1e-5) {
      x_axis = vect<double, 3>(pri_joint.Axis()[2], 0.0, 0.0);
    } else {
      x_axis /= x_axis_mag;
    }
    vect<double, 3> z_axis = x_axis % pri_joint.Axis();
    auto* trans = new SoTransform;
    SbMatrix rot_mat(x_axis[0], x_axis[1], x_axis[2], 0.0, pri_joint.Axis()[0],
                     pri_joint.Axis()[1], pri_joint.Axis()[2], 0.0, z_axis[0],
                     z_axis[1], z_axis[2], 0.0, 0.0, 0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);

    auto* col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);

    auto* core_cube = new SoCube;
    core_cube->width =
        aSG.mCharacteristicLength * 0.01;  // about 2% of char-len
    core_cube->height =
        aSG.mCharacteristicLength * 0.04;  // about 4% of char-len
    core_cube->depth =
        aSG.mCharacteristicLength * 0.01;  // about 2% of char-len
    sep->addChild(core_cube);

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* red_arrows = new SoTexture2;
    red_arrows->filename.setValue("./models/textures/red_up_arrows.tga");
    sep->addChild(red_arrows);

    auto* arrow_cyl = new SoCylinder;
    arrow_cyl->parts = SoCylinder::SIDES;
    arrow_cyl->radius = aSG.mCharacteristicLength * 0.01;
    arrow_cyl->height = aSG.mCharacteristicLength * 0.04;
    sep->addChild(arrow_cyl);

    if (!pri_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(pri_joint.BaseFrame()) ==
          aSG.mAnchor3DMap.end()) {
        aSG << pri_joint.BaseFrame();
      }
      aSG.mAnchor3DMap[pri_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::prismatic_joint_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& pri_joint = static_cast<const kte::prismatic_joint_2D&>(aModel);
    using std::atan2;

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    trans->rotation.setValue(
        SbVec3f(0.0, 0.0, 1.0),
        atan2(-(pri_joint.Axis()[0]), pri_joint.Axis()[1]));
    sep->addChild(trans);

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/prism_joint_2D.tga");
    sep->addChild(red_arrow);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!pri_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(pri_joint.BaseFrame()) ==
          aSG.mAnchor2DMap.end()) {
        aSG << pri_joint.BaseFrame();
      }
      aSG.mAnchor2DMap[pri_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::free_joint_3D::getStaticObjectType()) !=
             nullptr) {
    const auto& fr_joint = static_cast<const kte::free_joint_3D&>(aModel);

    auto* sep = new SoSeparator;

    auto* col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);

    auto* core_ball = new SoSphere;
    core_ball->radius =
        aSG.mCharacteristicLength * 0.01;  // diameter about 2% of char-len
    sep->addChild(core_ball);

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* rot_arrows = new SoTexture2;
    rot_arrows->filename.setValue("./models/textures/free_joint_3D_rot.tga");
    sep->addChild(rot_arrows);

    auto* arrow_ball = new SoSphere;
    arrow_ball->radius = aSG.mCharacteristicLength * 0.0125;
    sep->addChild(arrow_ball);

    auto* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/red_arrow.tga");
    sep->addChild(red_arrow);

    auto* rect_coords = new SoCoordinate3;
    auto* rect_texcoords = new SoTextureCoordinate2;

    // negative z face:
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 0.5, 1.0);
    rect_texcoords->point.set1Value(3, 0.5, 0.0);

    // positive x face:
    rect_coords->point.set1Value(4, aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(5, aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(6, aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(7, aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_texcoords->point.set1Value(4, 0.0, 0.0);
    rect_texcoords->point.set1Value(5, 0.0, 1.0);
    rect_texcoords->point.set1Value(6, 0.5, 1.0);
    rect_texcoords->point.set1Value(7, 0.5, 0.0);

    // negative y face:
    rect_coords->point.set1Value(8, -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(9, aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(10, aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_coords->point.set1Value(11, -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125,
                                 -aSG.mCharacteristicLength * 0.0125);
    rect_texcoords->point.set1Value(8, 0.0, 0.0);
    rect_texcoords->point.set1Value(9, 0.0, 1.0);
    rect_texcoords->point.set1Value(10, 0.5, 1.0);
    rect_texcoords->point.set1Value(11, 0.5, 0.0);

    sep->addChild(rect_coords);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    rect_face->numVertices.set1Value(1, 4);
    rect_face->numVertices.set1Value(2, 4);
    sep->addChild(rect_face);

    if (!fr_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(fr_joint.BaseFrame()) ==
          aSG.mAnchor3DMap.end()) {
        aSG << fr_joint.BaseFrame();
      }
      aSG.mAnchor3DMap[fr_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::free_joint_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& fr_joint = static_cast<const kte::free_joint_2D&>(aModel);

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* red_arrow = new SoTexture2;
    red_arrow->filename.setValue("./models/textures/free_joint_2D.tga");
    sep->addChild(red_arrow);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!fr_joint.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(fr_joint.BaseFrame()) ==
          aSG.mAnchor2DMap.end()) {
        aSG << fr_joint.BaseFrame();
      }
      aSG.mAnchor2DMap[fr_joint.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::rigid_link_3D::getStaticObjectType()) !=
             nullptr) {
    const auto& lnk_obj = static_cast<const kte::rigid_link_3D&>(aModel);

    vect<double, 3> y_axis = lnk_obj.PoseOffset().Position;
    double y_axis_mag = norm_2(y_axis);
    if (y_axis_mag < 1e-5) {  // nothing to draw...
      return aSG;
    }
    y_axis /= y_axis_mag;

    auto* sep = new SoSeparator;

    vect<double, 3> x_axis = y_axis % vect<double, 3>(0.0, 0.0, 1.0);
    double x_axis_mag = norm_2(x_axis);
    if (x_axis_mag < 1e-5) {
      x_axis = vect<double, 3>(y_axis[2], 0.0, 0.0);
    } else {
      x_axis /= x_axis_mag;
    }
    vect<double, 3> z_axis = x_axis % y_axis;
    auto* trans = new SoTransform;
    trans->translation.setValue(lnk_obj.PoseOffset().Position[0] * 0.5,
                                lnk_obj.PoseOffset().Position[1] * 0.5,
                                lnk_obj.PoseOffset().Position[2] * 0.5);
    SbMatrix rot_mat(x_axis[0], x_axis[1], x_axis[2], 0.0, y_axis[0], y_axis[1],
                     y_axis[2], 0.0, z_axis[0], z_axis[1], z_axis[2], 0.0, 0.0,
                     0.0, 0.0, 1.0);
    trans->rotation.setValue(SbRotation(rot_mat));
    sep->addChild(trans);

    auto* col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);

    auto* core_cyl = new SoCylinder;
    core_cyl->radius = aSG.mCharacteristicLength *
                       0.005;  // diameter about 1% of char-len (very thin rod).
    core_cyl->height = y_axis_mag;
    sep->addChild(core_cyl);

    if (!lnk_obj.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(lnk_obj.BaseFrame()) ==
          aSG.mAnchor3DMap.end()) {
        aSG << lnk_obj.BaseFrame();
      }
      aSG.mAnchor3DMap[lnk_obj.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::rigid_link_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& lnk_obj = static_cast<const kte::rigid_link_2D&>(aModel);
    using std::atan2;

    vect<double, 2> y_axis = lnk_obj.PoseOffset().Position;
    double y_axis_mag = norm_2(y_axis);

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    trans->rotation.setValue(SbVec3f(0.0, 0.0, 1.0),
                             atan2(-y_axis[0], y_axis[1]));
    sep->addChild(trans);

    auto* col = new SoBaseColor;
    col->rgb = SbColor(kte_gray_R, kte_gray_G, kte_gray_B);
    sep->addChild(col);

    auto* rod_coords = new SoCoordinate3;
    rod_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.005, 0.0,
                                0.0);
    rod_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.005,
                                y_axis_mag, 0.0);
    rod_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.005,
                                y_axis_mag, 0.0);
    rod_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.005, 0.0, 0.0);
    sep->addChild(rod_coords);

    auto* rod_face = new SoFaceSet;
    rod_face->numVertices.set1Value(0, 4);
    sep->addChild(rod_face);

    if (!lnk_obj.BaseFrame()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(lnk_obj.BaseFrame()) ==
          aSG.mAnchor2DMap.end()) {
        aSG << lnk_obj.BaseFrame();
      }
      aSG.mAnchor2DMap[lnk_obj.BaseFrame()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::spring_3D::getStaticObjectType()) != nullptr) {
    const auto& spr_obj = static_cast<const kte::spring_3D&>(aModel);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(spr_obj.Anchor1()) == aSG.mAnchor3DMap.end()) {
        aSG << spr_obj.Anchor1();
      }
      if (aSG.mAnchor3DMap.find(spr_obj.Anchor2()) == aSG.mAnchor3DMap.end()) {
        aSG << spr_obj.Anchor2();
      }
    }

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    sep->addChild(trans);

    auto* scal = new SoScale;
    sep->addChild(scal);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      aSG.mUpdateFuncs.emplace_back(detail::update_inter_frame_3D(
          spr_obj.Anchor1(), spr_obj.Anchor2(), trans, scal));
      aSG.mUpdateFuncs.back()();
    }

    auto* red_col = new SoBaseColor;
    red_col->rgb = SbColor(kte_red_R, kte_red_G, kte_red_B);
    sep->addChild(red_col);

    auto* sp_coords = new SoCoordinate3;
    {
      int i = 0;
      sp_coords->point.set1Value(i++, 0.0, -1.0, 0.0);
      sp_coords->point.set1Value(i++, 0.0, -0.9, 0.0);
      float h = -0.9;
      float a = 0.0;
      float r = aSG.mCharacteristicLength * 0.01;
      for (int j = 0; j < 181; ++j) {
        sp_coords->point.set1Value(i++, r * cos(a), h, r * sin(a));
        h += 0.01;
        a += 0.314159265359;  // PI / 10
      }
      sp_coords->point.set1Value(i++, 0.0, 0.9, 0.0);
      sp_coords->point.set1Value(i++, 0.0, 1.0, 0.0);
    }
    sep->addChild(sp_coords);

    auto* sp_set = new SoLineSet;
    sp_set->numVertices.set1Value(0, 185);
    sep->addChild(sp_set);

    /*
    SoBaseColor * white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    SoTexture2* twirls = new SoTexture2;
    twirls->filename.setValue("./models/textures/spring.tga");
    sep->addChild(twirls);

    SoCylinder* spr_cyl = new SoCylinder;
    spr_cyl->parts = SoCylinder::SIDES;
    spr_cyl->radius = aSG.mCharacteristicLength * 0.015;
    spr_cyl->height = 2.0;
    sep->addChild(spr_cyl);
    */
    aSG.mRootSwitch->addChild(
        sep);  // always at the root because "trans" is a global transformation.

  } else if (aModel.castTo(kte::spring_2D::getStaticObjectType()) != nullptr) {
    const auto& spr_obj = static_cast<const kte::spring_2D&>(aModel);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(spr_obj.Anchor1()) == aSG.mAnchor2DMap.end()) {
        aSG << spr_obj.Anchor1();
      }
      if (aSG.mAnchor2DMap.find(spr_obj.Anchor2()) == aSG.mAnchor2DMap.end()) {
        aSG << spr_obj.Anchor2();
      }
    }

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    sep->addChild(trans);

    auto* scal = new SoScale;
    sep->addChild(scal);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      aSG.mUpdateFuncs.emplace_back(detail::update_inter_frame_2D(
          spr_obj.Anchor1(), spr_obj.Anchor2(), trans, scal));
      aSG.mUpdateFuncs.back()();
    }

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* twirls = new SoTexture2;
    twirls->filename.setValue("./models/textures/spring.tga");
    sep->addChild(twirls);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01, -1.0,
                                 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01, 1.0,
                                 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01, 1.0, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01, -1.0,
                                 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    aSG.mRootSwitch->addChild(
        sep);  // always at the root because "trans" is a global transformation.

  } else if (aModel.castTo(kte::damper_3D::getStaticObjectType()) != nullptr) {
    const auto& dmp_obj = static_cast<const kte::damper_3D&>(aModel);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(dmp_obj.Anchor1()) == aSG.mAnchor3DMap.end()) {
        aSG << dmp_obj.Anchor1();
      }
      if (aSG.mAnchor3DMap.find(dmp_obj.Anchor2()) == aSG.mAnchor3DMap.end()) {
        aSG << dmp_obj.Anchor2();
      }
    }

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    sep->addChild(trans);

    auto* scal = new SoScale;
    sep->addChild(scal);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      aSG.mUpdateFuncs.emplace_back(detail::update_inter_frame_3D(
          dmp_obj.Anchor1(), dmp_obj.Anchor2(), trans, scal));
      aSG.mUpdateFuncs.back()();
    }

    // draw the center-line:

    auto* ln_sep = new SoSeparator;

    auto* ln_red_col = new SoBaseColor;
    ln_red_col->rgb = SbColor(kte_red_R, kte_red_G, kte_red_B);
    ln_sep->addChild(ln_red_col);

    auto* ln_coords = new SoCoordinate3;
    ln_coords->point.set1Value(0, 0.0, -1.0, 0.0);
    ln_coords->point.set1Value(1, 0.0, 1.0, 0.0);
    ln_sep->addChild(ln_coords);

    auto* ln_set = new SoLineSet;
    ln_set->numVertices.set1Value(0, 2);
    ln_sep->addChild(ln_set);

    sep->addChild(ln_sep);

    auto* red_col = new SoBaseColor;
    red_col->rgb = SbColor(kte_red_R, kte_red_G, kte_red_B);
    sep->addChild(red_col);

    // draw the outer hollow cylinder:
    auto* dmp_out_pos = new SoTranslation;
    dmp_out_pos->translation.setValue(0.0, 0.25, 0.0);
    sep->addChild(dmp_out_pos);

    auto* dmp_out_cyl = new SoCylinder;
    dmp_out_cyl->parts = SoCylinder::SIDES | SoCylinder::TOP;
    dmp_out_cyl->radius = aSG.mCharacteristicLength * 0.0075;
    dmp_out_cyl->height = 0.6;
    sep->addChild(dmp_out_cyl);

    // draw the inner solid cylinder:
    auto* dmp_in_pos = new SoTranslation;
    dmp_in_pos->translation.setValue(0.0, -0.5, 0.0);
    sep->addChild(dmp_in_pos);

    auto* dmp_in_cyl = new SoCylinder;
    dmp_in_cyl->radius = aSG.mCharacteristicLength * 0.005;
    dmp_in_cyl->height = 0.6;
    sep->addChild(dmp_in_cyl);

    aSG.mRootSwitch->addChild(
        sep);  // always at the root because "trans" is a global transformation.

  } else if (aModel.castTo(kte::damper_2D::getStaticObjectType()) != nullptr) {
    const auto& dmp_obj = static_cast<const kte::damper_2D&>(aModel);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(dmp_obj.Anchor1()) == aSG.mAnchor2DMap.end()) {
        aSG << dmp_obj.Anchor1();
      }
      if (aSG.mAnchor2DMap.find(dmp_obj.Anchor2()) == aSG.mAnchor2DMap.end()) {
        aSG << dmp_obj.Anchor2();
      }
    }

    auto* sep = new SoSeparator;

    auto* trans = new SoTransform;
    sep->addChild(trans);

    auto* scal = new SoScale;
    sep->addChild(scal);

    {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      aSG.mUpdateFuncs.emplace_back(detail::update_inter_frame_2D(
          dmp_obj.Anchor1(), dmp_obj.Anchor2(), trans, scal));
      aSG.mUpdateFuncs.back()();
    }

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* twirls = new SoTexture2;
    twirls->filename.setValue("./models/textures/dashpot.tga");
    sep->addChild(twirls);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01, -1.0,
                                 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01, 1.0,
                                 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01, 1.0, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01, -1.0,
                                 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    aSG.mRootSwitch->addChild(
        sep);  // always at the root because "trans" is a global transformation.

  } else if (aModel.castTo(kte::torsion_spring_3D::getStaticObjectType()) !=
             nullptr) {
    const auto& tor_spr = static_cast<const kte::torsion_spring_3D&>(aModel);
    using std::atan2;

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* spiral = new SoTexture2;
    spiral->filename.setValue("./models/textures/torsion_spring.tga");
    sep->addChild(spiral);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!tor_spr.Anchor1()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(tor_spr.Anchor1()) == aSG.mAnchor3DMap.end()) {
        aSG << tor_spr.Anchor1();
      }
      aSG.mAnchor3DMap[tor_spr.Anchor1()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::torsion_spring_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& tor_spr = static_cast<const kte::torsion_spring_2D&>(aModel);
    using std::atan2;

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* spiral = new SoTexture2;
    spiral->filename.setValue("./models/textures/torsion_spring.tga");
    sep->addChild(spiral);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!tor_spr.Anchor1()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(tor_spr.Anchor1()) == aSG.mAnchor2DMap.end()) {
        aSG << tor_spr.Anchor1();
      }
      aSG.mAnchor2DMap[tor_spr.Anchor1()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::torsion_damper_3D::getStaticObjectType()) !=
             nullptr) {
    const auto& tor_dmp = static_cast<const kte::torsion_damper_3D&>(aModel);
    using std::atan2;

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* spiral = new SoTexture2;
    spiral->filename.setValue("./models/textures/torsion_damper.tga");
    sep->addChild(spiral);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!tor_dmp.Anchor1()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(tor_dmp.Anchor1()) == aSG.mAnchor3DMap.end()) {
        aSG << tor_dmp.Anchor1();
      }
      aSG.mAnchor3DMap[tor_dmp.Anchor1()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::torsion_damper_2D::getStaticObjectType()) !=
             nullptr) {
    const auto& tor_dmp = static_cast<const kte::torsion_damper_2D&>(aModel);
    using std::atan2;

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* spiral = new SoTexture2;
    spiral->filename.setValue("./models/textures/torsion_damper.tga");
    sep->addChild(spiral);

    auto* rect_coords = new SoCoordinate3;
    rect_coords->point.set1Value(0, -aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(1, -aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(2, aSG.mCharacteristicLength * 0.01,
                                 aSG.mCharacteristicLength * 0.01, 0.0);
    rect_coords->point.set1Value(3, aSG.mCharacteristicLength * 0.01,
                                 -aSG.mCharacteristicLength * 0.01, 0.0);
    sep->addChild(rect_coords);

    auto* rect_texcoords = new SoTextureCoordinate2;
    rect_texcoords->point.set1Value(0, 0.0, 0.0);
    rect_texcoords->point.set1Value(1, 0.0, 1.0);
    rect_texcoords->point.set1Value(2, 1.0, 1.0);
    rect_texcoords->point.set1Value(3, 1.0, 0.0);
    sep->addChild(rect_texcoords);

    auto* rect_face = new SoFaceSet;
    rect_face->numVertices.set1Value(0, 4);
    sep->addChild(rect_face);

    if (!tor_dmp.Anchor1()) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(tor_dmp.Anchor1()) == aSG.mAnchor2DMap.end()) {
        aSG << tor_dmp.Anchor1();
      }
      aSG.mAnchor2DMap[tor_dmp.Anchor1()].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::inertia_3D::getStaticObjectType()) != nullptr) {
    const auto& cm_obj = static_cast<const kte::inertia_3D&>(aModel);

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* CoM_tex = new SoTexture2;
    CoM_tex->filename.setValue("./models/textures/CoM_gray_red.tga");
    sep->addChild(CoM_tex);

    auto* core_ball = new SoSphere;
    core_ball->radius = aSG.mCharacteristicLength *
                        0.01;  // diameter about 2% of char-len (very thin rod).
    sep->addChild(core_ball);

    if (!cm_obj.CenterOfMass()->mFrame) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor3DMap.find(cm_obj.CenterOfMass()->mFrame) ==
          aSG.mAnchor3DMap.end()) {
        aSG << cm_obj.CenterOfMass()->mFrame;
      }
      aSG.mAnchor3DMap[cm_obj.CenterOfMass()->mFrame].first->addChild(sep);
    }

  } else if (aModel.castTo(kte::inertia_2D::getStaticObjectType()) != nullptr) {
    const auto& cm_obj = static_cast<const kte::inertia_2D&>(aModel);

    auto* sep = new SoSeparator;

    auto* white_col = new SoBaseColor;
    white_col->rgb = SbColor(1.0, 1.0, 1.0);
    sep->addChild(white_col);

    auto* CoM_tex = new SoTexture2;
    CoM_tex->filename.setValue("./models/textures/CoM_gray_red.tga");
    sep->addChild(CoM_tex);

    auto* core_ball = new SoSphere;
    core_ball->radius = aSG.mCharacteristicLength *
                        0.01;  // diameter about 2% of char-len (very thin rod).
    sep->addChild(core_ball);

    if (!cm_obj.CenterOfMass()->mFrame) {
      aSG.mRootSwitch->addChild(sep);
    } else {
      std::unique_lock<std::recursive_mutex> lock_here(
          aSG.mAnchorUpdatingMutex);

      if (aSG.mAnchor2DMap.find(cm_obj.CenterOfMass()->mFrame) ==
          aSG.mAnchor2DMap.end()) {
        aSG << cm_obj.CenterOfMass()->mFrame;
      }
      aSG.mAnchor2DMap[cm_obj.CenterOfMass()->mFrame].first->addChild(sep);
    }
  }

  return aSG;
}

oi_scene_graph& operator<<(oi_scene_graph& aSG,
                           const kte::kte_map_chain& aModel) {

  const std::vector<std::shared_ptr<kte::kte_map>>& mdl_list = aModel.getKTEs();
  for (const auto& i : mdl_list) {
    aSG << *i;
  }

  return aSG;
}
}  // namespace ReaK::geom
