/**
 * \file planner_exec_engines.hpp
 *
 * This library defines functions and classes useful to execute path-planners.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_PLANNER_EXEC_ENGINES_HPP
#define REAK_PLANNER_EXEC_ENGINES_HPP

#ifndef RK_DISABLE_RRT_PLANNER
#include "rrt_path_planner.hpp"
#endif
#ifndef RK_DISABLE_PRM_PLANNER
#include "prm_path_planner.hpp"
#endif
#ifndef RK_DISABLE_RRTSTAR_PLANNER
#include "rrtstar_path_planner.hpp"
#endif
#ifndef RK_DISABLE_FADPRM_PLANNER
#include "fadprm_path_planner.hpp"
#endif
#ifndef RK_DISABLE_SBASTAR_PLANNER
#include "sbastar_path_planner.hpp"
#endif

#ifndef RK_DISABLE_PLANNER_DEFINITIONS
#ifndef RK_DISABLE_RRT_PLANNER
#include "rrt_path_planner.tpp"
#endif
#ifndef RK_DISABLE_PRM_PLANNER
#include "prm_path_planner.tpp"
#endif
#ifndef RK_DISABLE_RRTSTAR_PLANNER
#include "rrtstar_path_planner.tpp"
#endif
#ifndef RK_DISABLE_FADPRM_PLANNER
#include "fadprm_path_planner.tpp"
#endif
#ifndef RK_DISABLE_SBASTAR_PLANNER
#include "sbastar_path_planner.tpp"
#endif
#endif

#include "path_planner_options.hpp"

#include "p2p_planning_query.hpp"

#include "basic_sbmp_reporters.hpp"
#include "vlist_sbmp_report.hpp"

#include <filesystem>

namespace ReaK::pp {

struct monte_carlo_mp_engine {

  std::size_t mc_run_count;

  std::ofstream timing_output;
  std::ofstream sol_events_output;

  std::stringstream time_ss;
  std::stringstream cost_ss;
  std::stringstream sol_ss;

  monte_carlo_mp_engine(std::size_t aMCRuns, const std::string& aPlannerName,
                        const std::string& aOutputPathStem)
      : mc_run_count(aMCRuns) {

    std::filesystem::create_directory(aOutputPathStem.c_str());

    timing_output.open(aOutputPathStem + "/" + aPlannerName + "_times.txt");
    sol_events_output.open(aOutputPathStem + "/" + aPlannerName +
                           "_solutions.txt");

    timing_output << aPlannerName << std::endl;
    sol_events_output << aPlannerName << ", Solutions" << std::endl;
    std::cout << "Running " << aPlannerName << std::endl;
  }

  template <typename Topology>
  std::shared_ptr<any_sbmp_reporter_chain<Topology>> create_reporter(
      std::shared_ptr<Topology> /*unused*/) {

    auto report_chain = std::make_shared<any_sbmp_reporter_chain<Topology>>();

    report_chain->add_reporter(timing_sbmp_report<>(time_ss));
    report_chain->add_reporter(least_cost_sbmp_report<>(cost_ss, &sol_ss));

    return report_chain;
  }

  template <typename Topology>
  void operator()(const planning_option_collection& plan_options,
                  std::shared_ptr<sample_based_planner<Topology>> planner,
                  planning_query<Topology>& mc_query) {
    std::size_t mc_num_records =
        plan_options.max_vertices / plan_options.prog_interval;

    std::vector<double> vertex_counts(mc_num_records, 0.0);
    std::vector<std::size_t> num_remaining_planners(mc_num_records, 0);
    std::vector<std::size_t> num_successful_planners(mc_num_records, 0);
    std::vector<double> time_values(mc_num_records, 0.0);
    std::vector<double> best_costs(mc_num_records, 1.0e10);
    std::vector<double> worst_costs(mc_num_records, 0.0);
    std::vector<double> avg_costs(mc_num_records, 0.0);

    cost_ss << std::fixed;
    sol_ss << std::fixed;

    for (std::size_t i = 0; i < mc_run_count; ++i) {
      time_ss.clear();
      cost_ss.clear();
      sol_ss.clear();

      std::cout << "\r" << std::setw(10) << i << std::flush;

      mc_query.reset_solution_records();
      planner->reset_internal_state();
      planner->solve_planning_query(mc_query);

      std::size_t v_count = 0;
      std::size_t t_val = 0;
      std::string tmp;
      std::size_t j = 0;
      while (std::getline(time_ss, tmp) && (!tmp.empty())) {
        std::stringstream ss_tmp(tmp);
        ss_tmp >> v_count >> t_val;
        vertex_counts[j] =
            (double(v_count) +
             double(num_remaining_planners[j]) * vertex_counts[j]) /
            double(num_remaining_planners[j] + 1);
        time_values[j] = (double(t_val) +
                          double(num_remaining_planners[j]) * time_values[j]) /
                         double(num_remaining_planners[j] + 1);
        num_remaining_planners[j] += 1;
        ++j;
      }

      double c_val = 1e10;
      j = 0;
      while (std::getline(cost_ss, tmp) && (!tmp.empty())) {
        std::stringstream ss_tmp(tmp);
        ss_tmp >> v_count >> c_val;
        if (c_val < best_costs[j]) {
          best_costs[j] = c_val;
        }
        if (c_val > worst_costs[j]) {
          worst_costs[j] = c_val;
        }
        if (c_val < 1.0e9) {
          avg_costs[j] = (double(c_val) +
                          double(num_successful_planners[j]) * avg_costs[j]) /
                         double(num_successful_planners[j] + 1);
          num_successful_planners[j] += 1;
        }
        ++j;
      }

      while (j < mc_num_records) {
        if (c_val < best_costs[j]) {
          best_costs[j] = c_val;
        }
        if (c_val > worst_costs[j]) {
          worst_costs[j] = c_val;
        }
        if (c_val < 1.0e9) {
          avg_costs[j] = (double(c_val) +
                          double(num_successful_planners[j]) * avg_costs[j]) /
                         double(num_successful_planners[j] + 1);
          num_successful_planners[j] += 1;
        }
        ++j;
      }

      std::string first_sol_event;
      std::getline(sol_ss, first_sol_event);
      if (!first_sol_event.empty()) {
        sol_events_output << first_sol_event << std::endl;
      }
    }
    for (std::size_t i = 0; i < mc_num_records; ++i) {
      timing_output << std::setw(9) << i << " " << std::setw(9)
                    << vertex_counts[i] << " " << std::setw(9)
                    << num_remaining_planners[i] << " " << std::setw(9)
                    << num_successful_planners[i] << " " << std::setw(9)
                    << time_values[i] << " " << std::setw(9) << best_costs[i]
                    << " " << std::setw(9) << worst_costs[i] << " "
                    << std::setw(9) << avg_costs[i] << std::endl;
    }

    std::cout << "Done!" << std::endl;
  }
};

struct vlist_print_mp_engine {

  std::string vlist_file_path;
  std::any vlist_report_te;

  std::ofstream cost_out;
  std::ofstream sol_out;

  vlist_print_mp_engine(const std::string& aPlannerName,
                        const std::string& aOutputPathStem)
      : vlist_file_path(aOutputPathStem + "/" + aPlannerName) {

    std::filesystem::create_directory(aOutputPathStem.c_str());

    std::string qualified_output_path = aOutputPathStem + "/" + aPlannerName;
    std::filesystem::create_directory(qualified_output_path.c_str());

    cost_out.open(vlist_file_path + "_times.txt");
    sol_out.open(vlist_file_path + "_solutions.txt");

    std::cout << "Outputting " << aPlannerName << std::endl;
  }

  template <typename Topology>
  std::shared_ptr<any_sbmp_reporter_chain<Topology>> create_reporter(
      std::shared_ptr<Topology> /*unused*/) {

    auto report_chain = std::make_shared<any_sbmp_reporter_chain<Topology>>();

    report_chain->add_reporter(least_cost_sbmp_report<>(cost_out, &sol_out));
    report_chain->add_reporter(print_sbmp_progress<>());

    vlist_report_te =
        std::any(vlist_sbmp_report<any_mg_vertex_printer<Topology>>(
            vlist_file_path + "_",
            any_mg_vertex_printer<Topology>(BASIC_MOTION_GRAPH_KIND)));

    report_chain->add_reporter(std::ref(
        std::any_cast<vlist_sbmp_report<any_mg_vertex_printer<Topology>>&>(
            vlist_report_te)));

    return report_chain;
  }

  template <typename Topology>
  void operator()(const planning_option_collection& plan_options,
                  std::shared_ptr<sample_based_planner<Topology>> planner,
                  planning_query<Topology>& pp_query) {

    std::any_cast<vlist_sbmp_report<any_mg_vertex_printer<Topology>>&>(
        vlist_report_te)
        .print_to_stream.graph_kind = planner->get_motion_graph_kind();

    pp_query.reset_solution_records();
    planner->reset_internal_state();
    planner->solve_planning_query(pp_query);

    std::cout << "Done!" << std::endl;
    std::cout << "The shortest distance is: "
              << pp_query.get_best_solution_distance() << std::endl;
  }
};

struct differ_report_mp_engine {

  differ_sbmp_report_to_space<> differed_report;

  differ_report_mp_engine(double segment_steps, const std::string& aPlannerName,
                          const std::string& aOutputPathStem)
      : differed_report("", segment_steps) {

    std::filesystem::create_directory(aOutputPathStem.c_str());

    std::string qualified_output_path = aOutputPathStem + "/" + aPlannerName;
    std::filesystem::create_directory(qualified_output_path.c_str());

    differed_report.file_path = qualified_output_path + "/";

    std::cout << "Outputting " << aPlannerName << std::endl;
  }

  template <typename Topology>
  std::shared_ptr<any_sbmp_reporter_chain<Topology>> create_reporter(
      std::shared_ptr<Topology> /*unused*/) {

    auto report_chain = std::make_shared<any_sbmp_reporter_chain<Topology>>();

    report_chain->add_reporter(differed_report);
    report_chain->add_reporter(print_sbmp_progress<>());

    return report_chain;
  }

  template <typename Topology>
  void operator()(const planning_option_collection& plan_options,
                  std::shared_ptr<sample_based_planner<Topology>> planner,
                  planning_query<Topology>& pp_query) {
    pp_query.reset_solution_records();
    planner->reset_internal_state();
    planner->solve_planning_query(pp_query);

    std::cout << "Done!" << std::endl;
    std::cout << "The shortest distance is: "
              << pp_query.get_best_solution_distance() << std::endl;
  }
};

template <typename Topology, typename PlanEngine>
void execute_p2p_planner(const std::shared_ptr<Topology>& world_topo,
                         const planning_option_collection& plan_options,
                         std::size_t world_dimensionality, PlanEngine& engine,
                         const topology_point_type_t<Topology>& p_start,
                         const topology_point_type_t<Topology>& p_goal) {

  // Create the reporter chain.
  std::shared_ptr<any_sbmp_reporter_chain<Topology>> p_report_chain =
      engine.create_reporter(world_topo);

  // Create the point-to-point query:
  path_planning_p2p_query<Topology> pp_query("pp_query", world_topo, p_start,
                                             p_goal, plan_options.max_results);

  // Create the planner:
  std::shared_ptr<sample_based_planner<Topology>> world_planner;

#ifndef RK_DISABLE_RRT_PLANNER
  if (plan_options.planning_algo == 0) {  // RRT

    world_planner = std::make_shared<rrt_planner<Topology>>(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options, 0.1, 0.05, *p_report_chain);

  } else
#endif
#ifndef RK_DISABLE_RRTSTAR_PLANNER
      if (plan_options.planning_algo == 1) {  // RRT*

    world_planner = std::make_shared<rrtstar_planner<Topology>>(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options, 0.1, 0.05, world_dimensionality,
        *p_report_chain);

  } else
#endif
#ifndef RK_DISABLE_PRM_PLANNER
      if (plan_options.planning_algo == 2) {  // PRM

    world_planner = std::make_shared<prm_planner<Topology>>(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method, 0.1, 0.05,
        plan_options.max_random_walk, world_dimensionality, *p_report_chain);

  } else
#endif
#ifndef RK_DISABLE_FADPRM_PLANNER
      if (plan_options.planning_algo == 4) {  // FADPRM

    auto tmp = std::make_shared<fadprm_planner<Topology>>(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method, 0.1, 0.05,
        plan_options.max_random_walk, world_dimensionality, *p_report_chain);

    tmp->set_initial_relaxation(plan_options.init_relax);

    world_planner = tmp;

  } else
#endif
#ifndef RK_DISABLE_SBASTAR_PLANNER
      if (plan_options.planning_algo == 3) {  // SBA*

    auto tmp = std::make_shared<sbastar_planner<Topology>>(
        world_topo, plan_options.max_vertices, plan_options.prog_interval,
        plan_options.store_policy | plan_options.knn_method,
        plan_options.planning_options, 0.1, 0.05, plan_options.max_random_walk,
        world_dimensionality, *p_report_chain);

    tmp->set_initial_density_threshold(0.0);
    tmp->set_initial_relaxation(plan_options.init_relax);
    tmp->set_initial_SA_temperature(plan_options.init_SA_temp);

    world_planner = tmp;

  } else
#endif
  {
  }

  if (!world_planner) {
    return;
  }

  // Solve the planning problem:
  engine(plan_options, world_planner, pp_query);
}

}  // namespace ReaK::pp

#endif
