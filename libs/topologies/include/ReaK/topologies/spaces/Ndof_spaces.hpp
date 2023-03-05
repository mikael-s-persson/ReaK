/**
 * \file Ndof_spaces.hpp
 *
 * This library provides classes to represent N-dof spaces of either static or dynamic (run-time) dimensions,
 * and of differentiation order 0, 1 or 2.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date December 2012
 */

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

#ifndef REAK_NDOF_SPACES_HPP
#define REAK_NDOF_SPACES_HPP

#include "ReaK/core/base/defs.hpp"
#include "ReaK/core/base/serializable.hpp"
#include "ReaK/math/kinetostatics/gen_coord.hpp"
#include "ReaK/math/lin_alg/arithmetic_tuple.hpp"

#include "ReaK/topologies/spaces/default_random_sampler.hpp"
#include "ReaK/topologies/spaces/hyperbox_topology.hpp"
#include "ReaK/topologies/spaces/rate_limited_spaces.hpp"
#include "ReaK/topologies/spaces/tuple_distance_metrics.hpp"

namespace ReaK {

namespace pp {

template <typename Topo>
struct is_Ndof_space : std::false_type {};

template <typename Topo>
static constexpr bool is_Ndof_space_v = is_Ndof_space<Topo>::value;

template <typename Topo>
struct is_Ndof_rl_space : std::false_type {};

template <typename Topo>
static constexpr bool is_Ndof_rl_space_v = is_Ndof_rl_space<Topo>::value;

/**
 * This class defines the differentiation rule to apply either to lift a
 * point-difference (e.g. finite-difference) to the tangent space, or to descend
 * a tangent vector to a point-difference, for topologies whose point-difference
 * vectors are expressed as reach-time values of a Ndof space.
 */
template <typename Vector>
struct Ndof_reach_time_differentiation : public serializable {
  using self = Ndof_reach_time_differentiation<Vector>;

  Vector max_rate_reach_time;

  explicit Ndof_reach_time_differentiation(const Vector& aMaxRateReachTime)
      : max_rate_reach_time(aMaxRateReachTime) {}

  Ndof_reach_time_differentiation()
      : Ndof_reach_time_differentiation(Vector{}) {}

  /**
   * This function will lift a point-difference vector into its corresponding tangent vector.
   * This function performs a simple division, dp * (max_rate_reach_time / dt).
   * \tparam Vector1 The destination type, a point in the tangent space.
   * \tparam Vector2 The source type, a point-difference in the base space.
   * \tparam U A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param v The resulting point in the tangent space.
   * \param dp The point-difference that is being lifted.
   * \param dt The time-difference value (i.e. the difference in the independent variable).
   */
  template <typename Vector1, typename Vector2, typename U, typename TSpace>
  void lift(Vector1& v, const Vector2& dp, const U& dt,
            const TSpace& /*unused*/) const {
    v = dp;
    for (std::size_t i = 0; i < v.size(); ++i) {
      v[i] *= max_rate_reach_time[i] / dt;
    }
  }
  /**
   * This function will descend a tangent vector into its corresponding point-difference vector.
   * This function performs a simple multiplication, v * (dt / max_rate_reach_time).
   * \tparam Vector1 The destination type, a point-difference in the base space.
   * \tparam Vector2 The source type, a point in the tangent space.
   * \tparam U A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param dp The resulting point-difference in the base space.
   * \param v The point in the tangent space that is being descended.
   * \param dt The time-difference value (i.e. the difference in the independent variable).
   */
  template <typename Vector1, typename Vector2, typename U, typename TSpace>
  void descend(Vector1& dp, const Vector2& v, const U& dt,
               const TSpace& /*unused*/) const {
    dp = v;
    for (std::size_t i = 0; i < dp.size(); ++i) {
      dp[i] *= dt / max_rate_reach_time[i];
    }
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    A& RK_SERIAL_SAVE_WITH_NAME(max_rate_reach_time);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    A& RK_SERIAL_LOAD_WITH_NAME(max_rate_reach_time);
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2420002, 1,
                              "Ndof_reach_time_differentiation", serializable)
};

/**
 * This meta-function generates a N-dof space type for a given value-type, dimension and order.
 * A N-dof space type has its coordinates represented in one vector for each differentiation level (position, velocity,
 * acceleration).
 * The underlying topologies are all hyper-boxes with a Manhattan distance metrics inside one differentiation level
 * and between differentiation levels.
 * \tparam T The value-type of the space.
 * \tparam Size The dimension of the N-dof space (number of degrees of freedom) (0 for dynamic dimension).
 * \tparam Order The differentiation order of the N-dof space.
 */
template <typename T, unsigned int Size = 0, unsigned int Order = 0>
struct Ndof_space {
  using type =
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<
                               vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T, unsigned int Size = 0, unsigned int Order = 0>
using Ndof_space_t = typename Ndof_space<T, Size, Order>::type;

template <typename T, unsigned int Size>
struct Ndof_space<T, Size, 1> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T, unsigned int Size>
struct Ndof_space<T, Size, 2> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T>
struct Ndof_space<T, 0, 0> {
  using type = differentiable_space<time_topology,
                                    arithmetic_tuple<hyperbox_topology<
                                        vect_n<T>, manhattan_distance_metric>>>;
};

template <typename T>
struct Ndof_space<T, 0, 1> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>;
};

template <typename T>
struct Ndof_space<T, 0, 2> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>;
};

template <typename Vector>
struct is_Ndof_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, manhattan_distance_metric>>>>
    : std::true_type {};

template <typename Vector>
struct is_Ndof_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, manhattan_distance_metric>,
                     hyperbox_topology<Vector, manhattan_distance_metric>>>>
    : std::true_type {};

template <typename Vector>
struct is_Ndof_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, manhattan_distance_metric>,
                     hyperbox_topology<Vector, manhattan_distance_metric>,
                     hyperbox_topology<Vector, manhattan_distance_metric>>>>
    : std::true_type {};

/**
 * This function creates a N-dof space of 0th-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \return A N-dof space of 0th-order of the given dimension.
 */
template <unsigned int Size, typename Vector1>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, Size, 0>{arithmetic_tuple<
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>>(
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect<ValueType, Size>(aLowerBound),
          vect<ValueType, Size>(aUpperBound)))};
}

/**
 * This function creates a N-dof space of 1st-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A readable vector type.
 * \tparam Vector2 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \return A N-dof space of 1st-order of the given dimension.
 */
template <unsigned int Size, typename Vector1, typename Vector2>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound,
                     const Vector2& aSpeedLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, Size, 1>{arithmetic_tuple<
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>,
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>>(
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect<ValueType, Size>(aLowerBound),
          vect<ValueType, Size>(aUpperBound)),
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_1st_order_space", vect<ValueType, Size>(-aSpeedLimit),
          vect<ValueType, Size>(aSpeedLimit)))};
}

/**
 * This function creates a N-dof space of 2nd-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A readable vector type.
 * \tparam Vector2 A readable vector type.
 * \tparam Vector3 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \return A N-dof space of 2nd-order of the given dimension.
 */
template <unsigned int Size, typename Vector1, typename Vector2,
          typename Vector3>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound,
                     const Vector2& aSpeedLimit, const Vector3& aAccelLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, Size, 2>{arithmetic_tuple<
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>,
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>,
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>>(
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect<ValueType, Size>(aLowerBound),
          vect<ValueType, Size>(aUpperBound)),
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_1st_order_space", vect<ValueType, Size>(-aSpeedLimit),
          vect<ValueType, Size>(aSpeedLimit)),
      hyperbox_topology<vect<ValueType, Size>, manhattan_distance_metric>(
          "Ndof_2nd_order_space", vect<ValueType, Size>(-aAccelLimit),
          vect<ValueType, Size>(aAccelLimit)))};
}

/**
 * This function creates a N-dof space of 0th-order of any dimension.
 * \tparam Vector1 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \return A N-dof space of 0th-order of any dimension.
 */
template <typename Vector1>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, 0, 0>{arithmetic_tuple<
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>>(
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect_n<ValueType>(aLowerBound),
          vect_n<ValueType>(aUpperBound)))};
}

/**
 * This function creates a N-dof space of 1st-order of any dimension.
 * \tparam Vector1 A readable vector type.
 * \tparam Vector2 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \return A N-dof space of 1st-order of any dimension.
 */
template <typename Vector1, typename Vector2>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound,
                     const Vector2& aSpeedLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, 0, 1>{arithmetic_tuple<
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>,
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>>(
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect_n<ValueType>(aLowerBound),
          vect_n<ValueType>(aUpperBound)),
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_1st_order_space", vect_n<ValueType>(-aSpeedLimit),
          vect_n<ValueType>(aSpeedLimit)))};
}

/**
 * This function creates a N-dof space of 2nd-order of any dimension.
 * \tparam Vector1 A readable vector type.
 * \tparam Vector2 A readable vector type.
 * \tparam Vector3 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \return A N-dof space of 2nd-order of any dimension.
 */
template <typename Vector1, typename Vector2, typename Vector3>
auto make_Ndof_space(const Vector1& aLowerBound, const Vector1& aUpperBound,
                     const Vector2& aSpeedLimit, const Vector3& aAccelLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  return Ndof_space_t<ValueType, 0, 2>{arithmetic_tuple<
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>,
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>,
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>>(
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_0th_order_space", vect_n<ValueType>(aLowerBound),
          vect_n<ValueType>(aUpperBound)),
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_1st_order_space", vect_n<ValueType>(-aSpeedLimit),
          vect_n<ValueType>(aSpeedLimit)),
      hyperbox_topology<vect_n<ValueType>, manhattan_distance_metric>(
          "Ndof_2nd_order_space", vect_n<ValueType>(-aAccelLimit),
          vect_n<ValueType>(aAccelLimit)))};
}

/**
 * This meta-function generates a rate-limited N-dof space type for a given value-type, dimension and order.
 * \tparam T The value-type of the space.
 * \tparam Size The dimension of the N-dof space (number of degrees of freedom) (0 for dynamic dimension).
 * \tparam Order The differentiation order of the N-dof space.
 */
template <typename T, unsigned int Size = 0, unsigned int Order = 0>
struct Ndof_rl_space {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T, unsigned int Size = 0, unsigned int Order = 0>
using Ndof_rl_space_t = typename Ndof_rl_space<T, Size, Order>::type;

template <typename T, unsigned int Size>
struct Ndof_rl_space<T, Size, 1> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T, unsigned int Size>
struct Ndof_rl_space<T, Size, 2> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T>
struct Ndof_rl_space<T, 0, 0> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

template <typename T>
struct Ndof_rl_space<T, 0, 1> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

template <typename T>
struct Ndof_rl_space<T, 0, 2> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

template <typename Vector>
struct is_Ndof_rl_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<Vector>>>
    : std::true_type {};

template <typename Vector>
struct is_Ndof_rl_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, inf_norm_distance_metric>,
                     hyperbox_topology<Vector, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<Vector>>>
    : std::true_type {};

template <typename Vector>
struct is_Ndof_rl_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<Vector, inf_norm_distance_metric>,
                     hyperbox_topology<Vector, inf_norm_distance_metric>,
                     hyperbox_topology<Vector, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<Vector>>>
    : std::true_type {};

/**
 * This function creates a rate-limited N-dof space of 0th-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \return A N-dof space of 0th-order of the given dimension.
 */
template <unsigned int Size, typename Vector1, typename Vector2>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        const Vector2& aSpeedLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
  }
  return Ndof_rl_space_t<ValueType, Size, 0>{arithmetic_tuple<
      hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>>(
      hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
          "Ndof_0th_order_rl_space", vect<ValueType, Size>(aLowerBound),
          vect<ValueType, Size>(aUpperBound)))};
}

/**
 * This function creates a rate-limited N-dof space of 1st-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A writable vector type.
 * \tparam Vector3 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \return A N-dof space of 1st-order of the given dimension.
 */
template <unsigned int Size, typename Vector1, typename Vector2,
          typename Vector3>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        Vector2 aSpeedLimit, const Vector3& aAccelLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
    aSpeedLimit[i] /= aAccelLimit[i];
  }
  return Ndof_rl_space_t<ValueType, Size, 1>{
      arithmetic_tuple<
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>>(
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
              "Ndof_0th_order_rl_space", vect<ValueType, Size>(aLowerBound),
              vect<ValueType, Size>(aUpperBound)),
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
              "Ndof_1st_order_rl_space", vect<ValueType, Size>(-aSpeedLimit),
              vect<ValueType, Size>(aSpeedLimit))),
      manhattan_tuple_distance(),
      arithmetic_tuple<Ndof_reach_time_differentiation<vect<ValueType, Size>>>(
          Ndof_reach_time_differentiation<vect<ValueType, Size>>(aSpeedLimit))};
}

/**
 * This function creates a rate-limited N-dof space of 2nd-order of a given dimension.
 * \tparam Size The dimension of the space (number of degrees of freedom).
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A writable vector type.
 * \tparam Vector3 A writable vector type.
 * \tparam Vector4 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \param aJerkLimit The limits of the jerk on the generalized coordinates.
 * \return A N-dof space of 2nd-order of the given dimension.
 */
template <unsigned int Size, typename Vector1, typename Vector2,
          typename Vector3, typename Vector4>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        Vector2 aSpeedLimit, Vector3 aAccelLimit,
                        const Vector4& aJerkLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
    aSpeedLimit[i] /= aAccelLimit[i];
    aAccelLimit[i] /= aJerkLimit[i];
  }
  return Ndof_rl_space_t<ValueType, Size, 2>{
      arithmetic_tuple<
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>>(
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
              "Ndof_0th_order_rl_space", vect<ValueType, Size>(aLowerBound),
              vect<ValueType, Size>(aUpperBound)),
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
              "Ndof_1st_order_rl_space", vect<ValueType, Size>(-aSpeedLimit),
              vect<ValueType, Size>(aSpeedLimit)),
          hyperbox_topology<vect<ValueType, Size>, inf_norm_distance_metric>(
              "Ndof_2nd_order_rl_space", vect<ValueType, Size>(-aAccelLimit),
              vect<ValueType, Size>(aAccelLimit))),
      manhattan_tuple_distance(),
      arithmetic_tuple<Ndof_reach_time_differentiation<vect<ValueType, Size>>,
                       Ndof_reach_time_differentiation<vect<ValueType, Size>>>(
          Ndof_reach_time_differentiation<vect<ValueType, Size>>(
              vect<ValueType, Size>(aSpeedLimit)),
          Ndof_reach_time_differentiation<vect<ValueType, Size>>(
              vect<ValueType, Size>(aAccelLimit)))};
}

/**
 * This function creates a rate-limited N-dof space of 0th-order of a given dimension.
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \return A N-dof space of 0th-order of the given dimension.
 */
template <typename Vector1, typename Vector2>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        const Vector2& aSpeedLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
  }
  return Ndof_rl_space_t<ValueType, 0, 0>{arithmetic_tuple<
      hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>>(
      hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
          "Ndof_0th_order_rl_space", vect_n<ValueType>(aLowerBound),
          vect_n<ValueType>(aUpperBound)))};
}

/**
 * This function creates a rate-limited N-dof space of 2nd-order of any dimension.
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A writable vector type.
 * \tparam Vector3 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \return A N-dof space of 2nd-order of any dimension.
 */
template <typename Vector1, typename Vector2, typename Vector3>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        Vector2 aSpeedLimit, const Vector3& aAccelLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
    aSpeedLimit[i] /= aAccelLimit[i];
  }
  return Ndof_rl_space_t<ValueType, 0, 1>{
      arithmetic_tuple<
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>,
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>>(
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
              "Ndof_0th_order_rl_space", vect_n<ValueType>(aLowerBound),
              vect_n<ValueType>(aUpperBound)),
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
              "Ndof_1st_order_rl_space", vect_n<ValueType>(-aSpeedLimit),
              vect_n<ValueType>(aSpeedLimit))),
      manhattan_tuple_distance(),
      arithmetic_tuple<Ndof_reach_time_differentiation<vect_n<ValueType>>>(
          Ndof_reach_time_differentiation<vect_n<ValueType>>(
              vect_n<ValueType>(aSpeedLimit)))};
}

/**
 * This function creates a rate-limited N-dof space of 2nd-order of any dimension.
 * \tparam Vector1 A writable vector type.
 * \tparam Vector2 A writable vector type.
 * \tparam Vector3 A writable vector type.
 * \tparam Vector4 A readable vector type.
 * \param aLowerBound The lower bounds of the positions of the generalized coordinates.
 * \param aUpperBound The upper bounds of the positions of the generalized coordinates.
 * \param aSpeedLimit The limits of the velocities of the generalized coordinates.
 * \param aAccelLimit The limits of the accelerations on the generalized coordinates.
 * \param aJerkLimit The limits of the jerk on the generalized coordinates.
 * \return A N-dof space of 2nd-order of any dimension.
 */
template <typename Vector1, typename Vector2, typename Vector3,
          typename Vector4>
auto make_Ndof_rl_space(Vector1 aLowerBound, Vector1 aUpperBound,
                        Vector2 aSpeedLimit, Vector3 aAccelLimit,
                        const Vector4& aJerkLimit) {
  using ValueType = vect_value_type_t<Vector1>;
  for (std::size_t i = 0; i < aSpeedLimit.size(); ++i) {
    aLowerBound[i] /= aSpeedLimit[i];
    aUpperBound[i] /= aSpeedLimit[i];
    aSpeedLimit[i] /= aAccelLimit[i];
    aAccelLimit[i] /= aJerkLimit[i];
  }
  return Ndof_rl_space_t<ValueType, 0, 2>{
      arithmetic_tuple<
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>,
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>,
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>>(
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
              "Ndof_0th_order_rl_space", vect_n<ValueType>(aLowerBound),
              vect_n<ValueType>(aUpperBound)),
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
              "Ndof_1st_order_rl_space", vect_n<ValueType>(-aSpeedLimit),
              vect_n<ValueType>(aSpeedLimit)),
          hyperbox_topology<vect_n<ValueType>, inf_norm_distance_metric>(
              "Ndof_2nd_order_rl_space", vect_n<ValueType>(-aAccelLimit),
              vect_n<ValueType>(aAccelLimit))),
      manhattan_tuple_distance(),
      arithmetic_tuple<Ndof_reach_time_differentiation<vect_n<ValueType>>,
                       Ndof_reach_time_differentiation<vect_n<ValueType>>>(
          Ndof_reach_time_differentiation<vect_n<ValueType>>(
              vect_n<ValueType>(aSpeedLimit)),
          Ndof_reach_time_differentiation<vect_n<ValueType>>(
              vect_n<ValueType>(aAccelLimit)))};
}

template <typename T, unsigned int Size>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>> {
  using type =
      differentiable_space<time_topology,
                           arithmetic_tuple<hyperbox_topology<
                               vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T, unsigned int Size>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T, unsigned int Size>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
        hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
          hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>;
};

template <typename T, unsigned int Size>
struct get_rate_limited_space<differentiable_space<
    time_topology, arithmetic_tuple<hyperbox_topology<
                       vect<T, Size>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T, unsigned int Size>
struct get_rate_limited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
        hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T, unsigned int Size>
struct get_rate_limited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
        hyperbox_topology<vect<T, Size>, manhattan_distance_metric>,
        hyperbox_topology<vect<T, Size>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>,
          hyperbox_topology<vect<T, Size>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect<T, Size>>>;
};

template <typename T>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>> {
  using type = differentiable_space<time_topology,
                                    arithmetic_tuple<hyperbox_topology<
                                        vect_n<T>, manhattan_distance_metric>>>;
};

template <typename T>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                     hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>;
};

template <typename T>
struct get_rate_illimited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                     hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                     hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
    manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
          hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>;
};

template <typename T>
struct get_rate_limited_space<differentiable_space<
    time_topology, arithmetic_tuple<hyperbox_topology<
                       vect_n<T>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

template <typename T>
struct get_rate_limited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

template <typename T>
struct get_rate_limited_space<differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<T>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<T>, manhattan_distance_metric>>>> {
  using type = differentiable_space<
      time_topology,
      arithmetic_tuple<hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>,
                       hyperbox_topology<vect_n<T>, inf_norm_distance_metric>>,
      manhattan_tuple_distance, Ndof_reach_time_differentiation<vect_n<T>>>;
};

}  // namespace pp

template <typename Vector>
gen_coord<vect_value_type_t<Vector>> get_gen_coord(
    const arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i) {
  return gen_coord<vect_value_type_t<Vector>>(get<0>(pt)[i], get<1>(pt)[i],
                                              get<2>(pt)[i], 0.0);
}

template <typename Vector>
gen_coord<vect_value_type_t<Vector>> get_gen_coord(
    const arithmetic_tuple<Vector, Vector>& pt, std::size_t i) {
  return gen_coord<vect_value_type_t<Vector>>(get<0>(pt)[i], get<1>(pt)[i], 0.0,
                                              0.0);
}

template <typename Vector>
gen_coord<vect_value_type_t<Vector>> get_gen_coord(
    const arithmetic_tuple<Vector>& pt, std::size_t i) {
  return gen_coord<vect_value_type_t<Vector>>(get<0>(pt)[i], 0.0, 0.0, 0.0);
}

template <typename Vector>
void set_gen_coord(arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i,
                   const gen_coord<vect_value_type_t<Vector>>& p) {
  get<0>(pt)[i] = p.q;
  get<1>(pt)[i] = p.q_dot;
  get<2>(pt)[i] = p.q_ddot;
}

template <typename Vector>
void set_gen_coord(arithmetic_tuple<Vector, Vector>& pt, std::size_t i,
                   const gen_coord<vect_value_type_t<Vector>>& p) {
  get<0>(pt)[i] = p.q;
  get<1>(pt)[i] = p.q_dot;
}

template <typename Vector>
void set_gen_coord(arithmetic_tuple<Vector>& pt, std::size_t i,
                   const gen_coord<vect_value_type_t<Vector>>& p) {
  get<0>(pt)[i] = p.q;
}

template <typename Vector>
const vect_value_type_t<Vector>& get_position(
    const arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i) {
  return get<0>(pt)[i];
}

template <typename Vector>
const vect_value_type_t<Vector>& get_position(
    const arithmetic_tuple<Vector, Vector>& pt, std::size_t i) {
  return get<0>(pt)[i];
}

template <typename Vector>
const vect_value_type_t<Vector>& get_position(
    const arithmetic_tuple<Vector>& pt, std::size_t i) {
  return get<0>(pt)[i];
}

template <typename Vector>
void set_position(arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i,
                  const vect_value_type_t<Vector>& p) {
  get<0>(pt)[i] = p;
}

template <typename Vector>
void set_position(arithmetic_tuple<Vector, Vector>& pt, std::size_t i,
                  const vect_value_type_t<Vector>& p) {
  get<0>(pt)[i] = p;
}

template <typename Vector>
void set_position(arithmetic_tuple<Vector>& pt, std::size_t i,
                  const vect_value_type_t<Vector>& p) {
  get<0>(pt)[i] = p;
}

template <typename Vector>
const vect_value_type_t<Vector>& get_velocity(
    const arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i) {
  return get<1>(pt)[i];
}

template <typename Vector>
const vect_value_type_t<Vector>& get_velocity(
    const arithmetic_tuple<Vector, Vector>& pt, std::size_t i) {
  return get<1>(pt)[i];
}

template <typename Vector>
void set_velocity(arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i,
                  const vect_value_type_t<Vector>& p) {
  get<1>(pt)[i] = p;
}

template <typename Vector>
void set_velocity(arithmetic_tuple<Vector, Vector>& pt, std::size_t i,
                  const vect_value_type_t<Vector>& p) {
  get<1>(pt)[i] = p;
}

template <typename Vector>
const vect_value_type_t<Vector>& get_acceleration(
    const arithmetic_tuple<Vector, Vector, Vector>& pt, std::size_t i) {
  return get<2>(pt)[i];
}

template <typename Vector>
void set_acceleration(arithmetic_tuple<Vector, Vector, Vector>& pt,
                      std::size_t i, const vect_value_type_t<Vector>& p) {
  get<2>(pt)[i] = p;
}

}  // namespace ReaK

#include "ReaK/topologies/spaces/joint_space_limits.hpp"

namespace ReaK::pp {

#define RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(NDOF)                          \
  extern template class differentiable_space<                                 \
      time_topology, arithmetic_tuple<hyperbox_topology<                      \
                         vect<double, NDOF>, manhattan_distance_metric>>>;    \
  extern template class differentiable_space<                                 \
      time_topology,                                                          \
      arithmetic_tuple<                                                       \
          hyperbox_topology<vect<double, NDOF>, manhattan_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, manhattan_distance_metric>>>; \
  extern template class differentiable_space<                                 \
      time_topology,                                                          \
      arithmetic_tuple<                                                       \
          hyperbox_topology<vect<double, NDOF>, manhattan_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, manhattan_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, manhattan_distance_metric>>>;

RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(1)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(2)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(3)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(4)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(5)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(6)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(7)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(8)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(9)
RK_NDOF_SPACES_MAKE_NORMAL_EXTERN_DECL(10)

extern template class differentiable_space<
    time_topology, arithmetic_tuple<hyperbox_topology<
                       vect_n<double>, manhattan_distance_metric>>>;
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<double>, manhattan_distance_metric>>>;
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<double>, manhattan_distance_metric>,
        hyperbox_topology<vect_n<double>, manhattan_distance_metric>>>;

#define RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(NDOF)                    \
  extern template class differentiable_space<                                \
      time_topology, arithmetic_tuple<hyperbox_topology<                     \
                         vect<double, NDOF>, inf_norm_distance_metric>>>;    \
  extern template class differentiable_space<                                \
      time_topology,                                                         \
      arithmetic_tuple<                                                      \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>>>; \
  extern template class differentiable_space<                                \
      time_topology,                                                         \
      arithmetic_tuple<                                                      \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>>>;

RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(1)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(2)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(3)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(4)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(5)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(6)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(7)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(8)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(9)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DECL(10)

extern template class differentiable_space<
    time_topology, arithmetic_tuple<hyperbox_topology<
                       vect_n<double>, inf_norm_distance_metric>>>;
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>>>;
extern template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>>>;

#define RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(NDOF)             \
                                                                          \
  extern template Ndof_rl_space_t<double, NDOF, 0>                        \
  joint_limits_mapping<double>::make_rl_joint_space(                      \
      const Ndof_space_t<double, NDOF, 0>&) const;                        \
  extern template Ndof_rl_space_t<double, NDOF, 1>                        \
  joint_limits_mapping<double>::make_rl_joint_space(                      \
      const Ndof_space_t<double, NDOF, 1>&) const;                        \
  extern template Ndof_rl_space_t<double, NDOF, 2>                        \
  joint_limits_mapping<double>::make_rl_joint_space(                      \
      const Ndof_space_t<double, NDOF, 2>&) const;                        \
                                                                          \
  extern template Ndof_space_t<double, NDOF, 0>                           \
  joint_limits_mapping<double>::make_normal_joint_space(                  \
      const Ndof_rl_space_t<double, NDOF, 0>&) const;                     \
  extern template Ndof_space_t<double, NDOF, 1>                           \
  joint_limits_mapping<double>::make_normal_joint_space(                  \
      const Ndof_rl_space_t<double, NDOF, 1>&) const;                     \
  extern template Ndof_space_t<double, NDOF, 2>                           \
  joint_limits_mapping<double>::make_normal_joint_space(                  \
      const Ndof_rl_space_t<double, NDOF, 2>&) const;                     \
                                                                          \
  extern template topology_point_type_t<Ndof_rl_space_t<double, NDOF, 0>> \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_space_t<double, NDOF, 0>>& pt,     \
      const Ndof_space_t<double, NDOF, 0>&,                               \
      const Ndof_rl_space_t<double, NDOF, 0>&) const;                     \
  extern template topology_point_type_t<Ndof_rl_space_t<double, NDOF, 1>> \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_space_t<double, NDOF, 1>>& pt,     \
      const Ndof_space_t<double, NDOF, 1>&,                               \
      const Ndof_rl_space_t<double, NDOF, 1>&) const;                     \
  extern template topology_point_type_t<Ndof_rl_space_t<double, NDOF, 2>> \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_space_t<double, NDOF, 2>>& pt,     \
      const Ndof_space_t<double, NDOF, 2>&,                               \
      const Ndof_rl_space_t<double, NDOF, 2>&) const;                     \
                                                                          \
  extern template topology_point_type_t<Ndof_space_t<double, NDOF, 0>>    \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_rl_space_t<double, NDOF, 0>>& pt,  \
      const Ndof_rl_space_t<double, NDOF, 0>&,                            \
      const Ndof_space_t<double, NDOF, 0>&) const;                        \
  extern template topology_point_type_t<Ndof_space_t<double, NDOF, 1>>    \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_rl_space_t<double, NDOF, 1>>& pt,  \
      const Ndof_rl_space_t<double, NDOF, 1>&,                            \
      const Ndof_space_t<double, NDOF, 1>&) const;                        \
  extern template topology_point_type_t<Ndof_space_t<double, NDOF, 2>>    \
  joint_limits_mapping<double>::map_to_space(                             \
      const topology_point_type_t<Ndof_rl_space_t<double, NDOF, 2>>& pt,  \
      const Ndof_rl_space_t<double, NDOF, 2>&,                            \
      const Ndof_space_t<double, NDOF, 2>&) const;

RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(0)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(1)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(2)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(3)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(4)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(5)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(6)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(7)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(8)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(9)
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC(10)

}  // namespace ReaK::pp

#endif
