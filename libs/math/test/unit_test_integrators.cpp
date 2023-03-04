
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

#include <ReaK/core/base/defs.hpp>

#include <ReaK/math/integrators/fixed_step_integrators.hpp>
#include <ReaK/math/integrators/integrator.hpp>
#include <ReaK/math/integrators/pred_corr_integrators.hpp>
#include <ReaK/math/integrators/variable_step_integrators.hpp>

#include <ReaK/math/integrators/unit_test_integrators_problems.hpp>

#include <ReaK/core/serialization/protobuf_archiver.hpp>

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

void computeReferenceSolution(const ProblemPtr& prob,
                              SolutionTrace& sol_trace) {
  sol_trace.clear();
  sol_trace.push_back(std::pair<double, ReaK::vect_n<double>>(
      prob->getInitialTime(), prob->getInitialValue()));
  double rec_step = (prob->getFinalTime() - prob->getInitialTime()) * 1e-3;
  ReaK::dormand_prince45_integrator<double> integ(
      "reference_integrator_ode45", sol_trace.back().second,
      sol_trace.back().first, rec_step * 1e-3, prob, rec_step, rec_step * 1e-6,
      1e-6);

  std::cout << "Generating reference solution for problem: "
            << prob->getObjectType()->TypeName() << " with t in ["
            << prob->getInitialTime() << "," << prob->getFinalTime() << "]"
            << std::endl;

  for (double next_time = prob->getInitialTime() + rec_step;
       next_time < prob->getFinalTime() + 0.5 * rec_step;
       next_time += rec_step) {

    integ.integrate(next_time);

    sol_trace.push_back(std::pair<double, ReaK::vect_n<double>>(
        integ.getTime(),
        ReaK::vect_n<double>(integ.getStateBegin(), integ.getStateEnd())));

    std::cout << "\r" << std::setw(10) << next_time;
    for (std::size_t k = 0; k < sol_trace.back().second.size(); ++k) {
      std::cout << std::setw(10) << sol_trace.back().second[k];
    }
    std::cout << std::flush;
  }
  std::cout << std::endl << "Done!" << std::endl;
}

const RefProblemSet& getRefProblemSet() {
  using namespace ReaK;

  static bool first_pass = true;
  static RefProblemSet prob_set;
  if (first_pass) {

    {
      std::ifstream file_in("integ_records/hires.pbuf");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new HIRES_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);

        std::ofstream file_out("integ_records/hires.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }

#if 0
    {
      std::ifstream file_in("integ_records/pollution.pbuf");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new Pollution_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);
        
        std::ofstream file_out("integ_records/pollution.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }
#endif

    {
      std::ifstream file_in("integ_records/ringmod.pbuf");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(
            RefProblemPair(ProblemPtr(new RingModulator_iv_problem<double>()),
                           SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);

        std::ofstream file_out("integ_records/ringmod.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }

    {
      std::ifstream file_in("integ_records/vanderpol.pbuf");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new VanDerPol_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);

        std::ofstream file_out("integ_records/vanderpol.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }

    {
      std::ifstream file_in("integ_records/vanderpolmod.pbuf");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(
            RefProblemPair(ProblemPtr(new VanDerPolMod_iv_problem<double>()),
                           SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);

        std::ofstream file_out("integ_records/vanderpolmod.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }

    {
      std::ifstream file_in("integ_records/orego.pbuf");
      if (file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back(RefProblemPair(
            ProblemPtr(new Orego_iv_problem<double>()), SolutionTrace()));
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);

        std::ofstream file_out("integ_records/orego.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }

#if 0
    {
      std::ifstream file_in("integ_records/rober.pbuf");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new Rober_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);
        
        std::ofstream file_out("integ_records/rober.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }
#endif

#if 0
    {
      std::ifstream file_in("integ_records/e5.pbuf");
      if(file_in.is_open()) {
        serialization::protobuf_iarchive ar_in(file_in);
        RefProblemPair tmp;
        ar_in >> tmp;
        prob_set.push_back(tmp);
      } else {
        prob_set.push_back( RefProblemPair(ProblemPtr(new E5_iv_problem<double>()), SolutionTrace()) );
        computeReferenceSolution(prob_set.back().first, prob_set.back().second);
        
        std::ofstream file_out("integ_records/e5.pbuf");
        serialization::protobuf_oarchive ar_out(file_out);
        ar_out << prob_set.back();
      }
    }
#endif

    first_pass = false;
  }

  return prob_set;
}

TEST(Integrators, FixedStepIntegrators) {
  const RefProblemSet& ref_set = getRefProblemSet();

  for (std::size_t i = 0; i < ref_set.size(); ++i) {
    std::cout << "Reference solution is:" << std::endl;
    for (std::size_t j = 0; j < ref_set[i].second.size(); ++j) {
      std::cout << "\r" << std::setw(10) << ref_set[i].second[j].first;
      for (std::size_t k = 0; k < ref_set[i].second[j].second.size(); ++k)
        std::cout << std::setw(10) << ref_set[i].second[j].second[k];
      std::cout << std::flush;
    }
    std::cout << std::endl << "Done!" << std::endl;
  }
}

TEST(Integrators, VariableStepIntegrators) {}

TEST(Integrators, PredCorrIntegrators) {}

}  // namespace
}  // namespace ReaK
