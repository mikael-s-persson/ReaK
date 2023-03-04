
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

#include <iostream>

using std::ptrdiff_t;
using std::size_t;

#include <ReaK/topologies/spaces/differentiable_space.hpp>
#include <ReaK/topologies/spaces/line_topology.hpp>
#include <ReaK/topologies/spaces/temporal_space.hpp>
#include <ReaK/topologies/spaces/time_poisson_topology.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>

#include <ReaK/topologies/spaces/Ndof_limits.hpp>
#include <ReaK/topologies/spaces/Ndof_spaces.hpp>

#include <ReaK/core/recorders/ascii_recorder.hpp>
#include <ReaK/topologies/interpolation/sustained_acceleration_pulse.hpp>
#include <ReaK/topologies/interpolation/sustained_velocity_pulse.hpp>

#include <ReaK/topologies/interpolation/sustained_acceleration_pulse_Ndof_detail.hpp>
#include <ReaK/topologies/interpolation/sustained_velocity_pulse_Ndof_detail.hpp>
#include <memory>

#define TEST_SIZE 2

int main(int argc, char** argv) {

  std::shared_ptr<std::istream> task_src_ptr;
  std::string in_filename = "piped";
  if (argc < 2) {
    task_src_ptr =
        std::shared_ptr<std::istream>(&std::cin, ReaK::null_deleter());
  } else {
    std::string tmp_fn = argv[1];
    in_filename =
        std::string(std::find(tmp_fn.rbegin(), tmp_fn.rend(), '/').base(),
                    std::find(tmp_fn.rbegin(), tmp_fn.rend(), '.').base() - 1);
    std::cout << "Filename is: " << in_filename << std::endl;
    task_src_ptr = std::shared_ptr<std::istream>(new std::ifstream(argv[1]));
  };
  std::istream& task_src = *task_src_ptr;

  double time_step;
  double min_time;
  double max_time;
  task_src >> time_step >> min_time >> max_time;

  ReaK::vect<double, 2> max_vel;
  ReaK::vect<double, 2> max_accel;
  ReaK::vect<double, 2> start_pt;
  ReaK::vect<double, 2> start_vel;
  ReaK::vect<double, 2> end_pt;
  ReaK::vect<double, 2> end_vel;
  task_src >> max_vel >> max_accel >> start_pt >> start_vel >> end_pt >>
      end_vel;

  using SpaceType = ReaK::pp::Ndof_rl_space<double, 2, 2>::type;
  using BoxTopoType =
      ReaK::pp::hyperbox_topology<ReaK::vect<double, 2>,
                                  ReaK::pp::inf_norm_distance_metric>;
  using SpaceTupleType =
      ReaK::arithmetic_tuple<BoxTopoType, BoxTopoType, BoxTopoType>;
  using PointType = ReaK::pp::topology_traits<SpaceType>::point_type;

  ReaK::vect<double, TEST_SIZE> max_pos;
  for (std::size_t i = 0; i < TEST_SIZE; ++i) {
    max_pos[i] = 100.0;
  }

  using DiffRule =
      ReaK::pp::Ndof_reach_time_differentiation<ReaK::vect<double, 2>>;
  using DiffTuple = ReaK::arithmetic_tuple<DiffRule, DiffRule>;

  std::shared_ptr<SpaceType> topo = std::make_shared<SpaceType>(
      SpaceTupleType(BoxTopoType("pos_topo", -max_pos, max_pos),
                     BoxTopoType("vel_topo", -max_vel, max_vel),
                     BoxTopoType("acc_topo", -max_accel, max_accel)),
      ReaK::pp::manhattan_tuple_distance(),
      DiffTuple(DiffRule(max_vel), DiffRule(max_accel)));
  ReaK::pp::time_topology t_topo;

  PointType start_point(start_pt, start_vel, ReaK::vect<double, TEST_SIZE>());
  PointType end_point(end_pt, end_vel, ReaK::vect<double, TEST_SIZE>());
  PointType result_point(start_pt, start_vel, ReaK::vect<double, TEST_SIZE>());

  try {
    {
      std::stringstream ss;
      ss << "test_interp_results/" << in_filename << "_svp.ssv";

      ReaK::recorder::ascii_recorder output_rec(ss.str());
      output_rec << "time";
      for (std::size_t i = 0; i < TEST_SIZE; ++i) {
        std::stringstream ss2;
        ss2 << i;
        std::string s = ss2.str();
        output_rec << ("pos_" + s) << ("vel_" + s) << ("acc_" + s);
      };
      output_rec << ReaK::recorder::data_recorder::end_name_row;

      ReaK::vect<double, TEST_SIZE> peak_velocity;

      double dt_total = max_time - min_time;

      double min_dt_final =
          ReaK::pp::detail::svp_compute_Ndof_interpolation_data_impl(
              start_point, end_point, peak_velocity, *topo, t_topo, dt_total);

      std::cout << "min dt final = \t" << min_dt_final << std::endl;
      std::cout << "peak velocity = \t" << peak_velocity << std::endl;

      for (double dt = 0.0; ((dt_total >= min_dt_final) && (dt <= dt_total));
           dt += time_step) {

        ReaK::pp::detail::svp_Ndof_interpolate_impl<2>(
            result_point, start_point, end_point, peak_velocity, *topo, t_topo,
            dt, dt_total);

        output_rec << (min_time + dt);
        for (std::size_t i = 0; i < TEST_SIZE; ++i) {
          output_rec << ReaK::get<0>(result_point)[i]
                     << ReaK::get<1>(result_point)[i]
                     << ReaK::get<2>(result_point)[i];
        };
        output_rec << ReaK::recorder::data_recorder::end_value_row;
      };
      output_rec << ReaK::recorder::data_recorder::flush;
    };

#if 1
    {
      std::stringstream ss;
      ss << "test_interp_results/" << in_filename << "_sap.ssv";

      ReaK::recorder::ascii_recorder output_rec(ss.str());
      output_rec << "time";
      for (std::size_t i = 0; i < TEST_SIZE; ++i) {
        std::stringstream ss2;
        ss2 << i;
        std::string s = ss2.str();
        output_rec << ("pos_" + s) << ("vel_" + s) << ("acc_" + s);
      };
      output_rec << ReaK::recorder::data_recorder::end_name_row;

      ReaK::vect<double, TEST_SIZE> peak_velocity;

      double dt_total = max_time - min_time;

      double min_dt_final =
          ReaK::pp::detail::sap_compute_Ndof_interpolation_data_impl(
              start_point, end_point, peak_velocity, *topo,
              ReaK::pp::time_topology(), dt_total);

      for (double dt = 0.0; ((dt_total >= min_dt_final) && (dt <= dt_total));
           dt += time_step) {

        ReaK::pp::detail::sap_Ndof_interpolate_impl<2>(
            result_point, start_point, end_point, peak_velocity, *topo, t_topo,
            dt, dt_total);

        output_rec << (min_time + dt);
        for (std::size_t i = 0; i < TEST_SIZE; ++i) {
          output_rec << ReaK::get<0>(result_point)[i]
                     << ReaK::get<1>(result_point)[i]
                     << ReaK::get<2>(result_point)[i];
        };
        output_rec << ReaK::recorder::data_recorder::end_value_row;
      };
      output_rec << ReaK::recorder::data_recorder::flush;
    };
#endif

  } catch (std::exception& e) {
    std::cout
        << "Error: An exception was thrown during the execution of this test!"
        << std::endl;
    std::cout << "Message: " << e.what() << std::endl;
    return 1;
  };

  return 0;
};
