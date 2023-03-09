
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

#include "ReaK/mbd/kte/kte_map_chain.h"

#include "ReaK/mbd/kte/inertia.h"
#include "ReaK/mbd/kte/revolute_joint.h"
#include "ReaK/mbd/kte/rigid_link.h"

#include "ReaK/core/recorders/ascii_recorder.h"

#include "ReaK/core/serialization/xml_archiver.h"

#include "subarc/sa_root_node.hpp"
#include "subarc/sa_signal.hpp"

#include "ReaK/core/base/rk_typed_primitives.h"
#include "ReaK/mbd/kte/damper.h"
#include "ReaK/mbd/kte/spring.h"

#include <thread>

using namespace kte;
using namespace ReaK;
using namespace serialization;
using namespace sa;

class pendulum_simulation_node : public node {
 private:
  std::shared_ptr<signal<double>> input_force;
  std::shared_ptr<signal<double>> output_angle;
  std::shared_ptr<signal<double>> output_velocity;

  std::shared_ptr<gen_coord<double>> joint_coord;
  kte_map_chain pendulum;
  double sim_time;
  recorder::ascii_recorder output_rec;
  bool& is_done;

 protected:
  virtual void process() {

    double f_ext = input_force->getLast();

    joint_coord->q_ddot = 0.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl = joint_coord->f;

    joint_coord->q_ddot = 1.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl_in = joint_coord->f;

    f_nl_in = f_nl - f_nl_in;
    if (abs(f_nl_in) < 1E-10) {
      joint_coord->q_ddot = 0.0;
    } else {
      joint_coord->q_ddot = (f_nl + f_ext) / f_nl_in;
    }

    output_rec << sim_time << joint_coord->q << joint_coord->q_dot
               << joint_coord->q_ddot << f_nl << f_ext
               << recorder::data_recorder::end_value_row;

    joint_coord->q += joint_coord->q_dot * 0.001;
    joint_coord->q_dot += joint_coord->q_ddot * 0.001;
    sim_time += 0.001;
    std::this_thread::yield();

    if (sim_time > 20.0) {
      is_done = true;
    }

    std::cout << "\r" << sim_time;
    std::cout.flush();

    output_angle->setValue(joint_coord->q);
    output_velocity->setValue(joint_coord->q_dot);
  }

 public:
  pendulum_simulation_node(bool& aIsDone)
      : input_force(),
        output_angle(new signal<double>("joint_angle")),
        output_velocity(new signal<double>("joint_velocity")),
        joint_coord(),
        pendulum(),
        sim_time(0.0),
        output_rec("pendulum_ctrl_results.ssvdat"),
        is_done(aIsDone) {
    xml_iarchive pendulum_arc("pendulum.xml");
    iarchive& arc_ref = pendulum_arc;
    arc_ref >> joint_coord >> pendulum;

    output_rec << "time"
               << "q"
               << "qd"
               << "qdd"
               << "f_nl"
               << "f_ext" << recorder::data_recorder::end_name_row;
  }

  virtual ~pendulum_simulation_node() {
    output_rec << recorder::data_recorder::close;
  }

  virtual bool connectTo(std::shared_ptr<node> aNode) {
    std::weak_ptr<named_object> tmp = aNode->getOutputSignal("joint_torque");
    if (tmp.expired()) {
      return true;
    } else if (input_force =
                   rtti::rk_dynamic_ptr_cast<signal<double>>(tmp.lock())) {
      return true;
    }
    return false;
  }

  virtual bool isConnected() { return input_force; }

  virtual std::weak_ptr<named_object> getOutputSignal(
      const std::string& aSignalName) const {
    if (aSignalName == "joint_angle") {
      return output_angle;
    }
    if (aSignalName == "joint_velocity") {
      return output_velocity;
    }
    return std::weak_ptr<named_object>();
  }
};

class pendulum_control_node : public node {
 protected:
  std::shared_ptr<signal<double>> output_force;
  std::shared_ptr<signal<double>> input_angle;
  std::shared_ptr<signal<double>> input_velocity;

  std::shared_ptr<gen_coord<double>> joint_coord;
  kte_map_chain pendulum;

  double set_angle;
  const bool& is_done;

  virtual void process() {

    double prop_err = sin(set_angle) * cos(joint_coord->q) -
                      cos(set_angle) * cos(joint_coord->q);

    joint_coord->q_ddot = prop_err * 16.0 - joint_coord->q_dot * 8.0;

    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();

    output_force->setValue(-joint_coord->f);

    const unsigned int& w_a = input_angle->waitingListForNext();
    const unsigned int& w_v = input_velocity->waitingListForNext();
    while ((w_a) && (w_v) && (!is_done)) {
      std::this_thread::yield();
    }
    joint_coord->q = input_angle->getLast();
    joint_coord->q_dot = input_velocity->getLast();
  }

 public:
  pendulum_control_node(double aSetAngle, const bool& aIsDone)
      : output_force(new signal<double>("joint_torque")),
        input_angle(),
        input_velocity(),
        joint_coord(),
        pendulum(),
        set_angle(aSetAngle),
        is_done(aIsDone) {
    xml_iarchive pendulum_arc("pendulum.xml");
    iarchive& arc_ref = pendulum_arc;
    arc_ref >> joint_coord >> pendulum;
  }

  virtual ~pendulum_control_node() {}

  virtual bool connectTo(std::shared_ptr<node> aNode) {
    std::weak_ptr<named_object> tmp = aNode->getOutputSignal("joint_angle");
    if (!tmp.expired()) {
      input_angle = rtti::rk_dynamic_ptr_cast<signal<double>>(tmp.lock());
    }
    tmp = aNode->getOutputSignal("joint_velocity");
    if (!tmp.expired()) {
      input_velocity = rtti::rk_dynamic_ptr_cast<signal<double>>(tmp.lock());
    }
    return true;
  }

  virtual bool isConnected() { return ((input_angle) && (input_velocity)); }

  virtual std::weak_ptr<named_object> getOutputSignal(
      const std::string& aSignalName) const {
    if (aSignalName == "joint_torque") {
      return output_force;
    }
    return std::weak_ptr<named_object>();
  }
};

class pendulum_vmc_node : public pendulum_control_node {
 protected:
  std::shared_ptr<frame_2D<double>> end_frame;

  kte_map_chain spring_damper_system;

  virtual void process() {

    pendulum.doMotion();
    spring_damper_system.doMotion();
    pendulum.clearForce();
    spring_damper_system.clearForce();
    spring_damper_system.doForce();
    pendulum.doForce();

    output_force->setValue(-joint_coord->f);

    const unsigned int& w_a = input_angle->waitingListForNext();
    const unsigned int& w_v = input_velocity->waitingListForNext();
    while ((w_a) && (w_v) && (!is_done)) {
      std::this_thread::yield();
    }
    joint_coord->q = input_angle->getLast();
    joint_coord->q_dot = input_velocity->getLast();
  }

 public:
  pendulum_vmc_node(double aSetAngle, const bool& aIsDone)
      : pendulum_control_node(aSetAngle, aIsDone),
        end_frame(),
        spring_damper_system("spring_damper_system") {
    xml_iarchive pendulum_arc("pendulum.xml");
    iarchive& arc_ref = pendulum_arc;
    arc_ref >> joint_coord >> pendulum >> end_frame;

    auto fixed_anchor = std::make_shared<frame_2D<double>>(
        std::weak_ptr<pose_2D<double>>(), vect<double, 2>(0.0, 0.51),
        rot_mat_2D<double>(), vect<double, 2>(0.0, 0.0), 0.0,
        vect<double, 2>(0.0, 9.81), 0.0, vect<double, 2>(0.0, 0.0), 0.0);
    spring_damper_system << std::make_shared<spring_2D>("vmc_spring", end_frame,
                                                        fixed_anchor, 0.01,
                                                        -16.0)
                         << std::make_shared<damper_2D>("vmc_damper", end_frame,
                                                        fixed_anchor, -8.0);
    // spring_damper_system << std::make_shared< spring_2D>("vmc_spring",end_frame,fixed_anchor,0.01,-20.0)
    //                     << std::make_shared< damper_2D>("vmc_damper",end_frame,fixed_anchor,-15.0);
  };

  virtual ~pendulum_vmc_node(){};
};

int main() {

  root_node test_on_sim;

  bool done_flag;
  // auto ctrl_node = std::make_shared<pendulum_vmc_node>(M_PI * 0.5, done_flag);
  auto ctrl_node =
      std::make_shared<pendulum_control_node>(M_PI * 0.5, done_flag);
  auto sim_node = std::make_shared<pendulum_simulation_node>(done_flag);

  ctrl_node->connectTo(sim_node);
  sim_node->connectTo(ctrl_node);

  test_on_sim << sim_node << ctrl_node;

  test_on_sim.start();

  while (!done_flag)
    std::this_thread::yield();

  test_on_sim.stop();
};
