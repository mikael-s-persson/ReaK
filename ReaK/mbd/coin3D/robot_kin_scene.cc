/**
 * \file robot_kin_scene.cpp
 *
 * This application constructs the KTE-based kinematics models for a manipulator.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */

#include "ReaK/mbd/models/inverse_kinematics_model.h"
#include "ReaK/mbd/models/manip_ERA_arm.h"
#include "ReaK/mbd/models/manip_P3R3R_arm.h"
#include "ReaK/mbd/models/manip_SSRMS_arm.h"

#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/Qt/viewers/SoQtPlaneViewer.h>
#include <Inventor/events/SoKeyboardEvent.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoEventCallback.h>
#include <Inventor/nodes/SoLineSet.h>
#include <Inventor/nodes/SoSeparator.h>

#include "ReaK/mbd/coin3D/oi_scene_graph.h"

#include "ReaK/geometry/shapes/coord_arrows_3D.h"

#include "ReaK/mbd/kte/kte_map_chain.h"

#include "ReaK/core/serialization/xml_archiver.h"
#include "ReaK/math/optimization/optim_exceptions.h"

using namespace ReaK;

struct all_robot_info {
  //   kte::manip_ERA_kinematics builder;
  //   kte::manip_SSRMS_kinematics builder;
  kte::manip_P3R3R_kinematics builder;
  std::shared_ptr<kte::kte_map_chain> kin_chain;
  std::shared_ptr<frame_3D<double>> target_frame;
};

void keyboard_press_hdl(void* userData, SoEventCallback* eventCB) {
  all_robot_info* r_info = reinterpret_cast<all_robot_info*>(userData);
  const SoEvent* event = eventCB->getEvent();

  static bool IK_enabled = false;

  vect_n<double> j_pos = r_info->builder.getJointPositions();

  if (SO_KEY_PRESS_EVENT(event, Q)) {
    j_pos[0] += 0.01;
    //     if(j_pos[0] > r_info->builder.joint_upper_bounds[0])
    //       j_pos[0] = r_info->builder.joint_upper_bounds[0];
  } else if (SO_KEY_PRESS_EVENT(event, A)) {
    j_pos[0] -= 0.01;
    //     if(j_pos[0] < r_info->builder.joint_lower_bounds[0])
    //       j_pos[0] = r_info->builder.joint_lower_bounds[0];
  } else if (SO_KEY_PRESS_EVENT(event, W)) {
    j_pos[1] += 0.01;
    //     if(j_pos[1] > r_info->builder.joint_upper_bounds[1])
    //       j_pos[1] = r_info->builder.joint_upper_bounds[1];
  } else if (SO_KEY_PRESS_EVENT(event, S)) {
    j_pos[1] -= 0.01;
    //     if(j_pos[1] < r_info->builder.joint_lower_bounds[1])
    //       j_pos[1] = r_info->builder.joint_lower_bounds[1];
  } else if (SO_KEY_PRESS_EVENT(event, E)) {
    j_pos[2] += 0.01;
    //     if(j_pos[2] > r_info->builder.joint_upper_bounds[2])
    //       j_pos[2] = r_info->builder.joint_upper_bounds[2];
  } else if (SO_KEY_PRESS_EVENT(event, D)) {
    j_pos[2] -= 0.01;
    //     if(j_pos[2] < r_info->builder.joint_lower_bounds[2])
    //       j_pos[2] = r_info->builder.joint_lower_bounds[2];
  } else if (SO_KEY_PRESS_EVENT(event, R)) {
    j_pos[3] += 0.01;
    //     if(j_pos[3] > r_info->builder.joint_upper_bounds[3])
    //       j_pos[3] = r_info->builder.joint_upper_bounds[3];
  } else if (SO_KEY_PRESS_EVENT(event, F)) {
    j_pos[3] -= 0.01;
    //     if(j_pos[3] < r_info->builder.joint_lower_bounds[3])
    //       j_pos[3] = r_info->builder.joint_lower_bounds[3];
  } else if (SO_KEY_PRESS_EVENT(event, T)) {
    j_pos[4] += 0.01;
    //     if(j_pos[4] > r_info->builder.joint_upper_bounds[4])
    //       j_pos[4] = r_info->builder.joint_upper_bounds[4];
  } else if (SO_KEY_PRESS_EVENT(event, G)) {
    j_pos[4] -= 0.01;
    //     if(j_pos[4] < r_info->builder.joint_lower_bounds[4])
    //       j_pos[4] = r_info->builder.joint_lower_bounds[4];
  } else if (SO_KEY_PRESS_EVENT(event, Y)) {
    j_pos[5] += 0.01;
    //     if(j_pos[5] > r_info->builder.joint_upper_bounds[5])
    //       j_pos[5] = r_info->builder.joint_upper_bounds[5];
  } else if (SO_KEY_PRESS_EVENT(event, H)) {
    j_pos[5] -= 0.01;
    //     if(j_pos[5] < r_info->builder.joint_lower_bounds[5])
    //       j_pos[5] = r_info->builder.joint_lower_bounds[5];
  } else if (SO_KEY_PRESS_EVENT(event, U)) {
    j_pos[6] += 0.01;
    //     if(j_pos[6] > r_info->builder.joint_upper_bounds[6])
    //       j_pos[6] = r_info->builder.joint_upper_bounds[6];
  } else if (SO_KEY_PRESS_EVENT(event, J)) {
    j_pos[6] -= 0.01;
    //     if(j_pos[6] < r_info->builder.joint_lower_bounds[6])
    //       j_pos[6] = r_info->builder.joint_lower_bounds[6];
  } else if (SO_KEY_PRESS_EVENT(event, Z)) {
    r_info->target_frame->Position[2] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, X)) {
    r_info->target_frame->Position[2] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, UP_ARROW)) {
    r_info->target_frame->Position[1] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, DOWN_ARROW)) {
    r_info->target_frame->Position[1] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, LEFT_ARROW)) {
    r_info->target_frame->Position[0] -= 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, RIGHT_ARROW)) {
    r_info->target_frame->Position[0] += 0.01;
  } else if (SO_KEY_PRESS_EVENT(event, B)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(1.0, 0.0, 0.0));
  } else if (SO_KEY_PRESS_EVENT(event, N)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(0.0, 1.0, 0.0));
  } else if (SO_KEY_PRESS_EVENT(event, M)) {
    r_info->target_frame->Quat *= ReaK::axis_angle<double>(
        M_PI * 0.01, ReaK::vect<double, 3>(0.0, 0.0, 1.0));
  } else if (SO_KEY_PRESS_EVENT(event, P)) {
    IK_enabled = !IK_enabled;
  };

  if (IK_enabled) {
    try {
      *(r_info->builder.getDependentFrame3D(0)->mFrame) =
          *(r_info->target_frame);
      r_info->builder.doInverseMotion();
    } catch (ReaK::optim::infeasible_problem& e) {
      r_info->builder.setJointPositions(j_pos);
    };
  } else {
    r_info->builder.setJointPositions(j_pos);
  };

  r_info->kin_chain->doMotion();
};

int main(int argc, char** argv) {
  using namespace ReaK;

  all_robot_info r_info;

  r_info.kin_chain = r_info.builder.getKTEChain();

  r_info.target_frame = std::shared_ptr<frame_3D<double>>(new frame_3D<double>(
      std::weak_ptr<pose_3D<double>>(), vect<double, 3>(0.0, 5.0, 5.0),
      quaternion<double>(), vect<double, 3>(0.0, 0.0, 0.0),
      vect<double, 3>(0.0, 0.0, 0.0), vect<double, 3>(0.0, 0.0, 0.0),
      vect<double, 3>(0.0, 0.0, 0.0), vect<double, 3>(0.0, 0.0, 0.0),
      vect<double, 3>(0.0, 0.0, 0.0)));

  {
    QWidget* mainwin = SoQt::init(argc, argv, argv[0]);

    {

      r_info.kin_chain->doMotion();

      geom::oi_scene_graph sg;

      sg.setCharacteristicLength(10.0);

      sg << (*r_info.kin_chain);
      sg << geom::coord_arrows_3D("target_arrows", r_info.target_frame,
                                  pose_3D<double>(), 0.3);

      SoSeparator* root = new SoSeparator;
      root->ref();

      root->addChild(sg.getSceneGraph());

      SoEventCallback* keypressCB = new SoEventCallback;
      keypressCB->addEventCallback(SoKeyboardEvent::getClassTypeId(),
                                   keyboard_press_hdl, &r_info);
      root->addChild(keypressCB);

      sg.enableAnchorUpdates();

      // Use one of the convenient SoQt viewer classes.
      // SoQtPlaneViewer * eviewer = new SoQtPlaneViewer(mainwin);
      SoQtExaminerViewer* eviewer = new SoQtExaminerViewer(mainwin);
      eviewer->setSceneGraph(root);
      eviewer->show();

      // Pop up the main window.
      SoQt::show(mainwin);
      // Loop until exit.
      SoQt::mainLoop();

      sg.disableAnchorUpdates();

      // delete cone_rot_anim;
      // Clean up resources.
      delete eviewer;
      root->unref();
    };

    SoQt::done();
  };

  return 0;
};
