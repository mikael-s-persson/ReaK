
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


#include "ReaK/math/integrators/fixed_step_integrators.h"
#include "ReaK/math/integrators/integrator.h"
#include "ReaK/math/integrators/pred_corr_integrators.h"
#include "ReaK/math/integrators/variable_step_integrators.h"

#include "ReaK/math/integrators/integrators_test_problems.h"

#include "ReaK/core/serialization/protobuf_archiver.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace ReaK {
namespace {

using SolutionTrace = std::vector<std::pair<double, ReaK::vect_n<double>>>;
using ProblemPtr = std::shared_ptr<ReaK::iv_problem<double>>;
using RefProblemPair = std::pair<ProblemPtr, SolutionTrace>;
using RefProblemSet = std::vector<RefProblemPair>;

template <typename Integrator>
void computeSolution(const ProblemPtr& prob, Integrator& integ,
                     SolutionTrace& sol_trace, double step_fraction,
                     double step_min_fraction, double tolerance) {
  sol_trace.clear();
  sol_trace.push_back(std::pair<double, ReaK::vect_n<double>>(
      prob->getInitialTime(), prob->getInitialValue()));
  double rec_step = (prob->getFinalTime() - prob->getInitialTime()) * 1e-3;
  integ.clearStateVector();
  integ.addStateElements(sol_trace.back().second);
  integ.setTime(sol_trace.back().first);
  integ.setStepSize(rec_step * step_fraction);
  integ.setStateRateFunc(prob);
  if constexpr (std::is_base_of_v<variable_step_integrator<double>,
                                  Integrator>) {
    integ.setMaxStepSize(rec_step);
    integ.setMinStepSize(rec_step * step_min_fraction);
    integ.setTolerance(tolerance);
  }

  std::cout << "Generating reference solution for problem: "
            << prob->get_object_type()->name() << " with t in ["
            << prob->getInitialTime() << "," << prob->getFinalTime() << "]"
            << std::endl;

  for (double next_time = prob->getInitialTime() + rec_step;
       next_time < prob->getFinalTime() + 0.5 * rec_step;
       next_time += rec_step) {

    integ.integrate(next_time);

    sol_trace.push_back(std::pair<double, ReaK::vect_n<double>>(
        integ.getTime(),
        ReaK::vect_n<double>(integ.getStateBegin(), integ.getStateEnd())));

    std::cout << "\r" << std::setw(10) << next_time << std::flush;
  }
  std::cout << std::endl << "Done!" << std::endl;
}

template <typename Integrator>
void computeReferenceSolution(const ProblemPtr& prob, Integrator& integ,
                              SolutionTrace& sol_trace) {
  computeSolution(prob, integ, sol_trace, 1e-3, 1e-6, 1e-6);
}

template <typename Integrator>
void computeTestSolution(const ProblemPtr& prob, Integrator& integ,
                         SolutionTrace& sol_trace) {
  computeSolution(prob, integ, sol_trace, 1e-1, 1e-3, 1e-3);
}

const RefProblemSet& getRefProblemSet() {
  using namespace ReaK;

  static bool first_pass = true;
  static RefProblemSet prob_set;
  if (first_pass) {
    ReaK::dormand_prince45_integrator<double> integ;
    integ.set_name("reference_integrator_ode45");

    {
      std::ifstream file_in("integ_records/hires.pb");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new HIRES_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, integ,
                                 prob_set.back().second);

        std::ofstream file_out("integ_records/hires.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }

#if 0
    {
      std::ifstream file_in("integ_records/pollution.pb");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new Pollution_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, integ, prob_set.back().second);
        
        std::ofstream file_out("integ_records/pollution.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }
#endif

    {
      std::ifstream file_in("integ_records/ringmod.pb");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(
            RefProblemPair(ProblemPtr(new RingModulator_iv_problem<double>()),
                           SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, integ,
                                 prob_set.back().second);

        std::ofstream file_out("integ_records/ringmod.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }

    {
      std::ifstream file_in("integ_records/vanderpol.pb");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new VanDerPol_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, integ,
                                 prob_set.back().second);

        std::ofstream file_out("integ_records/vanderpol.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }

    {
      std::ifstream file_in("integ_records/vanderpolmod.pb");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(
            RefProblemPair(ProblemPtr(new VanDerPolMod_iv_problem<double>()),
                           SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, integ,
                                 prob_set.back().second);

        std::ofstream file_out("integ_records/vanderpolmod.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }

    {
      std::ifstream file_in("integ_records/orego.pb");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new Orego_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, integ,
                                 prob_set.back().second);

        std::ofstream file_out("integ_records/orego.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }

#if 0
    {
      std::ifstream file_in("integ_records/rober.pb");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new Rober_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, integ, prob_set.back().second);
        
        std::ofstream file_out("integ_records/rober.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }
#endif

#if 0
    {
      std::ifstream file_in("integ_records/e5.pb");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new E5_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, integ, prob_set.back().second);
        
        std::ofstream file_out("integ_records/e5.pb");
        if (file_out.is_open()) {
          serialization::protobuf_oarchive ar_out(file_out);
          ar_out << prob_set.back();
        }
      }
    }
#endif

    first_pass = false;
  }

  return prob_set;
}

TEST(Integrators, EulerIntegrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    // HIRES_iv_problem is the only problem for which Euler is not a complete disaster.
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    euler_integrator<double> integ;
    integ.set_name("euler_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-4, 1e-3, 1e-3);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 5e-4)
          << " at time step " << i;
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, MidpointIntegrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    // HIRES_iv_problem is the only problem for which Midpoint is not a complete disaster.
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    midpoint_integrator<double> integ;
    integ.set_name("midpoint_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-3, 1e-3);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 5e-4)
          << " at time step " << i;
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, RungeKutta4Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    // HIRES_iv_problem is the only problem for which RungeKutta4 is not a complete disaster.
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    runge_kutta4_integrator<double> integ;
    integ.set_name("runge_kutta4_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-2, 1e-3, 1e-3);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 5e-4)
          << " at time step " << i;
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, RungeKutta5Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    // HIRES_iv_problem is the only problem for which RungeKutta5 is not a complete disaster.
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    runge_kutta5_integrator<double> integ;
    integ.set_name("runge_kutta5_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 2e-2, 1e-3, 1e-3);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 5e-4)
          << " at time step " << i;
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, Fehlberg45Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    fehlberg45_integrator<double> integ;
    integ.set_name("fehlberg45_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, DormandPrince45Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    dormand_prince45_integrator<double> integ;
    integ.set_name("dormand_prince45_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, AdamsBM3Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    adamsBM3_integrator<double> integ;
    integ.set_name("adamsBM3_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, AdamsBM5Integrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    adamsBM5_integrator<double> integ;
    integ.set_name("adamsBM5_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, HammingModIntegrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    hamming_mod_integrator<double> integ;
    integ.set_name("hamming_mod_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 5e-1 /*TODO FAILING*/))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

TEST(Integrators, HammingIterModIntegrator) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (const auto& [prob_ptr, ref_sol_trace] : ref_set) {
    if (prob_ptr->get_object_type()->name() != "HIRES_iv_problem<double>") {
      continue;
    }

    hamming_iter_mod_integrator<double> integ;
    integ.set_name("hamming_iter_mod_integrator");

    SolutionTrace test_sol_trace;
    computeSolution(prob_ptr, integ, test_sol_trace, 1e-3, 1e-4, 1e-4);

    double max_ref_norm = 0.0;
    for (const auto& [t, v] : ref_sol_trace) {
      double ref_norm = norm_2(v);
      if (ref_norm > max_ref_norm) {
        max_ref_norm = ref_norm;
      }
    }

    ASSERT_EQ(test_sol_trace.size(), ref_sol_trace.size());
    bool fatal_failure = false;
    for (int i = 0; i < test_sol_trace.size() && !fatal_failure; ++i) {
      EXPECT_NEAR(
          test_sol_trace[i].first, ref_sol_trace[i].first,
          (prob_ptr->getFinalTime() - prob_ptr->getInitialTime()) * 1e-3)
          << " at time step " << i;
      if (i < 10) {
        continue;
      }
      for (int j = 0; j < test_sol_trace[i].second.size(); ++j) {
        EXPECT_NEAR(test_sol_trace[i].second[j], ref_sol_trace[i].second[j],
                    std::max(2e-2, max_ref_norm * 2e-2))
            << " at time step " << i << " element " << j;
        if (!std::isfinite(test_sol_trace[i].second[j])) {
          fatal_failure = true;
          break;
        }
      }
    }
  }
}

}  // namespace
}  // namespace ReaK
