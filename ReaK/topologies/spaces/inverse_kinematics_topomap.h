/**
 * \file inverse_kinematics_topomap.h
 *
 * This library provides classes that define topological mappings between a joint-space (generalized
 * coordinates) and the end-effector frame of a serial manipulator.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
 */

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

#ifndef REAK_TOPOLOGIES_SPACES_INVERSE_KINEMATICS_TOPOMAP_H_
#define REAK_TOPOLOGIES_SPACES_INVERSE_KINEMATICS_TOPOMAP_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/math/kinetostatics/gen_coord.h"
#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include "ReaK/mbd/models/inverse_kinematics_model.h"
#include "ReaK/mbd/models/manip_clik_calculator.h"

#include "ReaK/topologies/spaces/joint_space_limits.h"
#include "ReaK/topologies/spaces/joint_space_topologies.h"
#include "ReaK/topologies/spaces/ndof_spaces.h"
#include "ReaK/topologies/spaces/se2_topologies.h"
#include "ReaK/topologies/spaces/se3_topologies.h"

#include "ReaK/topologies/spaces/direct_kinematics_topomap.h"
#include "ReaK/topologies/spaces/inverse_kinematics_topomap_detail.h"

namespace ReaK::pp {

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_inverse_kin_map : public shared_object {
 public:
  using self = manip_inverse_kin_map;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::inverse_kinematics_model> model;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model which can do the inverse kinematics calculation.
   */
  explicit manip_inverse_kin_map(
      const std::shared_ptr<kte::inverse_kinematics_model>& aModel =
          std::shared_ptr<kte::inverse_kinematics_model>())
      : model(aModel) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam OutSpace The type of the output space (joint-space).
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename OutSpace>
  topology_point_type_t<OutSpace> extract_from_model(
      const OutSpace& space_out) const {

    model->doInverseMotion();

    topology_point_type_t<OutSpace> result;
    detail::read_joint_coordinates_impl(result, space_out, model);
    return result;
  }

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    return extract_from_model(space_out);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model);
  };
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400014, 1, "manip_inverse_kin_map",
                              shared_object)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates,
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_mapping<double>>
class manip_rl_inverse_kin_map : public shared_object {
 public:
  using self = manip_rl_inverse_kin_map<RateLimitMap>;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::inverse_kinematics_model> model;
  /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
  RateLimitMap joint_limits_map;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
   */
  explicit manip_rl_inverse_kin_map(
      const std::shared_ptr<kte::inverse_kinematics_model>& aModel =
          std::shared_ptr<kte::inverse_kinematics_model>(),
      const RateLimitMap& aJointLimitMap = RateLimitMap())
      : model(aModel), joint_limits_map(aJointLimitMap) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model, at the current end-effector configuration.
   * \tparam OutSpace The type of the output space (rate-limited joint-space).
   * \param space_out The output space, i.e. the rate-limited joint-space.
   * \return A point in the output space, i.e. the rate-limited joint coordinates.
   */
  template <typename OutSpace>
  topology_point_type_t<OutSpace> extract_from_model(
      const OutSpace& space_out) const {

    model->doInverseMotion();

    using NormalJointSpace = typename get_rate_illimited_space<OutSpace>::type;
    topology_point_type_t<NormalJointSpace> result_inter;
    detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(),
                                        model);

    return joint_limits_map.map_to_space(result_inter, NormalJointSpace(),
                                         space_out);
  }

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (rate-limited joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the rate-limited joint-space.
   * \return A point in the output space, i.e. the rate-limited joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    return extract_from_model(space_out);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model) &
        RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model) &
        RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400015, 1, "manip_rl_inverse_kin_map",
                              shared_object)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_clik_kin_map : public shared_object {
 public:
  using self = manip_clik_kin_map;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::direct_kinematics_model> model;
  /** This holds the inverse kinematics calculator factory. */
  mutable kte::manip_clik_calculator clik_calc;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  explicit manip_clik_kin_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel =
          std::shared_ptr<kte::direct_kinematics_model>(),
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator =
          std::shared_ptr<optim::cost_evaluator>(),
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99)
      : model(aModel),
        clik_calc(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta,
                  aTau) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    clik_calc.solveInverseKinematics();

    detail::read_joint_coordinates_impl(result, space_out, model);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model) & RK_SERIAL_SAVE_WITH_NAME(clik_calc);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model) & RK_SERIAL_LOAD_WITH_NAME(clik_calc);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400016, 1, "manip_clik_kin_map",
                              shared_object)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates,
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_mapping<double>>
class manip_rl_clik_kin_map : public shared_object {
 public:
  using self = manip_rl_clik_kin_map<RateLimitMap>;

  /** This data member points to a manipulator kinematics model to use for the mappings performed. */
  std::shared_ptr<kte::direct_kinematics_model> model;
  /** This holds the inverse kinematics calculator factory. */
  mutable kte::manip_clik_calculator clik_calc;
  /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
  RateLimitMap joint_limits_map;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  explicit manip_rl_clik_kin_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel =
          std::shared_ptr<kte::direct_kinematics_model>(),
      const RateLimitMap& aJointLimitMap = RateLimitMap(),
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator =
          std::shared_ptr<optim::cost_evaluator>(),
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99)
      : model(aModel),
        clik_calc(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta,
                  aTau),
        joint_limits_map(aJointLimitMap) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (rate-limited joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the rate-limited joint-space.
   * \return A point in the output space, i.e. the rate-limited joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    clik_calc.solveInverseKinematics();

    using NormalJointSpace = typename get_rate_illimited_space<OutSpace>::type;
    topology_point_type_t<NormalJointSpace> result_inter;
    detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(),
                                        model);
    result = joint_limits_map.map_to_space(result_inter, NormalJointSpace(),
                                           space_out);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(model) & RK_SERIAL_SAVE_WITH_NAME(clik_calc) &
        RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(model) & RK_SERIAL_LOAD_WITH_NAME(clik_calc) &
        RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400017, 1, "manip_rl_clik_kin_map",
                              shared_object)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename JointStateType>
class manip_clik_fig_kin_map : public manip_clik_kin_map {
 public:
  using self = manip_clik_fig_kin_map<JointStateType>;

  JointStateType initial_guess;

  /**
   * Parametrized Constructor.
   * \param aInitGuess The initial guess to act as the start for the inverse kinematics search.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  explicit manip_clik_fig_kin_map(
      const JointStateType& aInitGuess = JointStateType(),
      const std::shared_ptr<kte::direct_kinematics_model>& aModel =
          std::shared_ptr<kte::direct_kinematics_model>(),
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator =
          std::shared_ptr<optim::cost_evaluator>(),
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99)
      : manip_clik_kin_map(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter,
                           aTol, aEta, aTau),
        initial_guess(aInitGuess) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;

    detail::write_joint_coordinates_impl(initial_guess, space_out, model);

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    clik_calc.solveInverseKinematics();

    detail::read_joint_coordinates_impl(result, space_out, model);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    manip_clik_kin_map::save(
        A, manip_clik_kin_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(initial_guess);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    manip_clik_kin_map::load(
        A, manip_clik_kin_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(initial_guess);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240001A, 1, "manip_clik_fig_kin_map",
                              manip_clik_kin_map)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates,
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename JointStateType,
          typename RateLimitMap = joint_limits_mapping<double>>
class manip_rl_clik_fig_kin_map : public manip_rl_clik_kin_map<RateLimitMap> {
 public:
  using base_type = manip_rl_clik_kin_map<RateLimitMap>;
  using self = manip_rl_clik_fig_kin_map<JointStateType, RateLimitMap>;

  JointStateType initial_guess;

  /**
   * Parametrized Constructor.
   * \param aInitGuess The initial guess to act as the start for the inverse kinematics search.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   */
  explicit manip_rl_clik_fig_kin_map(
      const JointStateType& aInitGuess = {},
      const std::shared_ptr<kte::direct_kinematics_model>& aModel = {},
      const RateLimitMap& aJointLimitMap = {},
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator = {},
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99)
      : base_type(aModel, aJointLimitMap, aCostEvaluator, aMaxRadius, aMu,
                  aMaxIter, aTol, aEta, aTau),
        initial_guess(aInitGuess) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (rate-limited joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the rate-limited joint-space.
   * \return A point in the output space, i.e. the rate-limited joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;
    using NormalJointSpace = typename get_rate_illimited_space<OutSpace>::type;

    auto ip_inter = this->joint_limits_map.map_to_space(
        initial_guess, space_out, NormalJointSpace());
    detail::write_joint_coordinates_impl(ip_inter, NormalJointSpace(),
                                         this->model);

    detail::write_dependent_coordinates_impl(pt, space_in, this->model);

    this->clik_calc.solveInverseKinematics();

    topology_point_type_t<NormalJointSpace> result_inter;
    detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(),
                                        this->model);
    result = this->joint_limits_map.map_to_space(result_inter,
                                                 NormalJointSpace(), space_out);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(initial_guess);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(initial_guess);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240001B, 1, "manip_rl_clik_fig_kin_map",
                              base_type)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_clik_rnd_restart_map : public manip_clik_kin_map {
 public:
  using self = manip_clik_rnd_restart_map;

  std::size_t restart_count;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   * \param aRestartCount The number of restarts to try before giving up on the inverse kinematics search.
   */
  explicit manip_clik_rnd_restart_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel = {},
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator = {},
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99,
      std::size_t aRestartCount = 10)
      : manip_clik_kin_map(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter,
                           aTol, aEta, aTau),
        restart_count(aRestartCount) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;

    detail::write_dependent_coordinates_impl(pt, space_in, model);

    clik_calc.readDesiredFromModel();

    for (std::size_t i = 0; i < restart_count; ++i) {

      detail::write_joint_coordinates_impl(
          get(random_sampler, space_out)(space_out), space_out, model);

      vect_n<double> x = clik_calc.readJointStatesFromModel();

      try {
        clik_calc.runOptimizer(x);
      } catch (singularity_error& e) {
      } catch (maximum_iteration& e) {
      } catch (optim::infeasible_problem& e) {};

      if (norm_2(clik_calc.computeStatesError(x)) < clik_calc.tol * 10.0) {
        clik_calc.writeJointStatesToModel(x);
        break;
      }
      if (i == restart_count - 1) {
        throw optim::infeasible_problem(
            "The inverse kinematics problem cannot be solved!");
      }
    }

    detail::read_joint_coordinates_impl(result, space_out, model);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    manip_clik_kin_map::save(
        A, manip_clik_kin_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(restart_count);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    manip_clik_kin_map::load(
        A, manip_clik_kin_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(restart_count);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240001C, 1, "manip_clik_rnd_restart_map",
                              manip_clik_kin_map)
};

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates,
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_mapping<double>>
class manip_rl_clik_rnd_restart_map
    : public manip_rl_clik_kin_map<RateLimitMap> {
 public:
  using base_type = manip_rl_clik_kin_map<RateLimitMap>;
  using self = manip_rl_clik_rnd_restart_map<RateLimitMap>;

  std::size_t restart_count;

  /**
   * Parametrized Constructor.
   * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
   * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
   * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
   * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
   * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is
   * positive and should start with a rather large value (relative to the scale of the function) and will be
   * progressively decreased by the algorithm as it progresses).
   * \param aMaxIter The maximum number of iterations to perform.
   * \param aTol The tolerance on the norm of the gradient (and thus the step size).
   * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
   * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint
   * (barrier).
   * \param aRestartCount The number of restarts to try before giving up on the inverse kinematics search.
   */
  explicit manip_rl_clik_rnd_restart_map(
      const std::shared_ptr<kte::direct_kinematics_model>& aModel = {},
      const RateLimitMap& aJointLimitMap = {},
      const std::shared_ptr<optim::cost_evaluator>& aCostEvaluator = {},
      double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300,
      double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99,
      std::size_t aRestartCount = 10)
      : base_type(aModel, aJointLimitMap, aCostEvaluator, aMaxRadius, aMu,
                  aMaxIter, aTol, aEta, aTau),
        restart_count(aRestartCount){};

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (rate-limited joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the rate-limited joint-space.
   * \return A point in the output space, i.e. the rate-limited joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  topology_point_type_t<OutSpace> map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {
    topology_point_type_t<OutSpace> result;
    using NormalJointSpace = typename get_rate_illimited_space<OutSpace>::type;

    detail::write_dependent_coordinates_impl(pt, space_in, this->model);

    this->clik_calc.readDesiredFromModel();

    for (std::size_t i = 0; i < restart_count; ++i) {

      auto ip_inter = this->joint_limits_map.map_to_space(
          get(random_sampler, space_out)(space_out), space_out,
          NormalJointSpace());
      detail::write_joint_coordinates_impl(ip_inter, NormalJointSpace(),
                                           this->model);

      vect_n<double> x = this->clik_calc.readJointStatesFromModel();

      try {
        this->clik_calc.runOptimizer(x);
      } catch (singularity_error& e) {
      } catch (maximum_iteration& e) {
      } catch (optim::infeasible_problem& e) {};

      if (norm_2(this->clik_calc.computeStatesError(x)) <
          this->clik_calc.tol * 10.0) {
        this->clik_calc.writeJointStatesToModel(x);
        break;
      }
      if (i == restart_count - 1) {
        throw optim::infeasible_problem(
            "The inverse kinematics problem cannot be solved!");
      }
    }

    topology_point_type_t<NormalJointSpace> result_inter;
    detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(),
                                        this->model);
    result = this->joint_limits_map.map_to_space(result_inter,
                                                 NormalJointSpace(), space_out);

    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(restart_count);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(restart_count);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240001D, 1,
                              "manip_rl_clik_rnd_restart_map", base_type)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_INVERSE_KINEMATICS_TOPOMAP_H_
