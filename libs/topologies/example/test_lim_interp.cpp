
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

#include <math.h>

#include <iostream>

using std::ptrdiff_t;
using std::size_t;

#include <ReaK/topologies/spaces/differentiable_space.hpp>
#include <ReaK/topologies/spaces/line_topology.hpp>
#include <ReaK/topologies/spaces/temporal_space.hpp>
#include <ReaK/topologies/spaces/time_poisson_topology.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>

#include <ReaK/core/recorders/ascii_recorder.hpp>
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse.hpp>
#include <ReaK/topologies/interpolation/sustained_velocity_pulse.hpp>
#include <memory>

int main(int argc, char** argv) {

  if (argc < 11) {

    std::cout << "Error: Arguments to the program were incorrect!" << std::endl
              << "Usage:" << std::endl
              << "\t\t./test_interp [time_step] [interp_time_step] [min_time] "
                 "[max_time] [max_vel] [max_accel] "
                 "[start_pt] [start_vel] [end_pt] [end_vel]"
              << std::endl;

    return 1;
  };

  double time_step = NAN;
  {
    std::stringstream ss(argv[1]);
    ss >> time_step;
  };

  double interp_time_step = NAN;
  {
    std::stringstream ss(argv[2]);
    ss >> interp_time_step;
  };

  double min_time = NAN;
  {
    std::stringstream ss(argv[3]);
    ss >> min_time;
  };

  double max_time = NAN;
  {
    std::stringstream ss(argv[4]);
    ss >> max_time;
  };

  double max_vel = NAN;
  {
    std::stringstream ss(argv[5]);
    ss >> max_vel;
  };

  double max_accel = NAN;
  {
    std::stringstream ss(argv[6]);
    ss >> max_accel;
  };

  double start_pt = NAN;
  {
    std::stringstream ss(argv[7]);
    ss >> start_pt;
  };

  double start_vel = NAN;
  {
    std::stringstream ss(argv[8]);
    ss >> start_vel;
  };

  double end_pt = NAN;
  {
    std::stringstream ss(argv[9]);
    ss >> end_pt;
  };

  double end_vel = NAN;
  {
    std::stringstream ss(argv[10]);
    ss >> end_vel;
  };

  using SpaceTupleType =
      ReaK::arithmetic_tuple<ReaK::pp::line_segment_topology<double>,
                             ReaK::pp::line_segment_topology<double>,
                             ReaK::pp::line_segment_topology<double>,
                             ReaK::pp::line_segment_topology<double>>;

  using TopoType =
      ReaK::pp::differentiable_space<ReaK::pp::time_poisson_topology,
                                     SpaceTupleType>;
  using PointType = ReaK::pp::topology_traits<TopoType>::point_type;
  using TempTopoType =
      ReaK::pp::temporal_space<TopoType, ReaK::pp::time_poisson_topology>;
  using TempPointType = ReaK::pp::topology_traits<TempTopoType>::point_type;

  std::shared_ptr<TempTopoType> topo = std::make_shared<TempTopoType>(
      "temporal_space",
      TopoType(SpaceTupleType(
          ReaK::pp::line_segment_topology<double>("pos_topo", -20.0, 20.0),
          ReaK::pp::line_segment_topology<double>("vel_topo", -max_vel,
                                                  max_vel),
          ReaK::pp::line_segment_topology<double>("acc_topo", -max_accel,
                                                  max_accel),
          ReaK::pp::line_segment_topology<double>("jerk_topo", -20.0, 20.0))));

  std::vector<TempPointType> pts;
  pts.emplace_back(0.0, PointType(start_pt, start_vel, 0.0, 0.0));
  pts.emplace_back(min_time, PointType(end_pt, end_vel, 0.0, 0.0));

  try {
    for (double current_end_time = min_time; current_end_time < max_time;
         current_end_time += interp_time_step) {
      std::stringstream ss;
      ss << "test_interp_results/svp_interp_" << current_end_time << ".ssv";
      RK_NOTICE(1, "Creating file: '" << ss.str() << "'");
      ReaK::recorder::ascii_recorder output_rec(ss.str());
      output_rec << "time"
                 << "pos"
                 << "vel"
                 << "acc"
                 << "jerk" << ReaK::recorder::data_recorder::end_name_row;

      pts[1].time = current_end_time;

      ReaK::pp::svp_interp_traj<TempTopoType> interp(pts.begin(), pts.end(),
                                                     topo);

      for (double t = 0.0; t <= current_end_time; t += time_step) {
        TempPointType p = interp.get_point_at_time(t);
        output_rec << p.time << ReaK::get<0>(p.pt) << ReaK::get<1>(p.pt)
                   << ReaK::get<2>(p.pt) << ReaK::get<3>(p.pt)
                   << ReaK::recorder::data_recorder::end_value_row;
      };
      output_rec << ReaK::recorder::data_recorder::flush;
    };

    for (double current_end_time = min_time; current_end_time < max_time;
         current_end_time += interp_time_step) {
      std::stringstream ss;
      ss << "test_interp_results/sap_interp_" << current_end_time << ".ssv";
      RK_NOTICE(1, "Creating file: '" << ss.str() << "'");
      ReaK::recorder::ascii_recorder output_rec(ss.str());
      output_rec << "time"
                 << "pos"
                 << "vel"
                 << "acc"
                 << "jerk" << ReaK::recorder::data_recorder::end_name_row;

      pts[1].time = current_end_time;

      ReaK::pp::sap_interp_traj<TempTopoType> interp(pts.begin(), pts.end(),
                                                     topo);

      for (double t = 0.0; t <= current_end_time; t += time_step) {
        TempPointType p = interp.get_point_at_time(t);
        output_rec << p.time << ReaK::get<0>(p.pt) << ReaK::get<1>(p.pt)
                   << ReaK::get<2>(p.pt) << ReaK::get<3>(p.pt)
                   << ReaK::recorder::data_recorder::end_value_row;
      };
      output_rec << ReaK::recorder::data_recorder::flush;
    };

  } catch (std::exception& e) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    std::cout << "Message: " << e.what() << std::endl;
    return 1;
  };

  using SpaceTupleType2 = ReaK::arithmetic_tuple<TopoType, TopoType>;
  using TopoType2 = ReaK::pp::metric_space_tuple<SpaceTupleType2>;
  using PointType2 = ReaK::pp::topology_traits<TopoType2>::point_type;
  using TempTopoType2 =
      ReaK::pp::temporal_space<TopoType2, ReaK::pp::time_poisson_topology>;
  using TempPointType2 = ReaK::pp::topology_traits<TempTopoType2>::point_type;

  std::shared_ptr<TempTopoType2> topo2 = std::make_shared<TempTopoType2>(
      "temporal_space_tuple",
      TopoType2(SpaceTupleType2(
          SpaceTupleType(
              ReaK::pp::line_segment_topology<double>("pos_topo", -20.0, 20.0),
              ReaK::pp::line_segment_topology<double>("vel_topo", -max_vel,
                                                      max_vel),
              ReaK::pp::line_segment_topology<double>("acc_topo", -max_accel,
                                                      max_accel),
              ReaK::pp::line_segment_topology<double>("jerk_topo", -20.0,
                                                      20.0)),
          SpaceTupleType(
              ReaK::pp::line_segment_topology<double>("pos_topo", -20.0, 20.0),
              ReaK::pp::line_segment_topology<double>("vel_topo", -max_vel,
                                                      max_vel),
              ReaK::pp::line_segment_topology<double>("acc_topo", -max_accel,
                                                      max_accel),
              ReaK::pp::line_segment_topology<double>("jerk_topo", -20.0,
                                                      20.0)))));

  std::vector<TempPointType2> pts2;
  pts2.emplace_back(0.0, PointType2(PointType(start_pt, start_vel, 0.0, 0.0),
                                    PointType(start_pt, start_vel, 0.0, 0.0)));
  pts2.emplace_back(min_time, PointType2(PointType(end_pt, end_vel, 0.0, 0.0),
                                         PointType(end_pt, end_vel, 0.0, 0.0)));

  ReaK::pp::svp_interp_traj<TempTopoType2> interp2(pts2.begin(), pts2.end(),
                                                   topo2);

  return 0;
};
