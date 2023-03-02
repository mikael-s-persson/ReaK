/**
 * \file Ndof_limits.hpp
 *
 * This library provides classes to help create and manipulate N-dof topologies in over
 * a joint-space with limits (speed, acceleration, and jerk limits).
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

#ifndef REAK_NDOF_LIMITS_HPP
#define REAK_NDOF_LIMITS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>

#include "Ndof_spaces.hpp"

namespace ReaK::pp {

/**
 * This class template stores a set of vectors to represent the rate-limits on the Ndofs
 * of a manipulator. Basically, this class is just a POD class, but it also provides functions
 * to construct a rate-limited Ndof-space from a normal Ndof-space, or vice-versa. Also,
 * it can act as a mapping between rate-limited joint coordinates and normal joint coordinates.
 * \tparam T The value type of the underlying Ndof-space.
 */
template <typename T, unsigned int Size = 0>
struct Ndof_limits : public named_object {
  /** Holds the lower bounds for all generalized coordinates. */
  vect<T, Size> lower_bounds;
  /** Holds the upper bounds for all generalized coordinates. */
  vect<T, Size> upper_bounds;
  /** Holds the speed limit for all generalized coordinates. */
  vect<T, Size> speed_limits;
  /** Holds the acceleration limit for all generalized coordinates. */
  vect<T, Size> accel_limits;
  /** Holds the jerk limit for all generalized coordinates. */
  vect<T, Size> jerk_limits;

  using value_type = T;
  using self = Ndof_limits<T, Size>;

  /**
   * Default constructor.
   */
  explicit Ndof_limits(const std::string& aName) : named_object() {
    this->setName(aName);
  }

  Ndof_limits() : Ndof_limits("") {}

  /**
   * This function constructs a rate-limited N-dof space of 0th differentiation order.
   * \return A rate-limited N-dof space of 0th differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, Size, 0> make_rl_0th_space() const {
    return make_Ndof_rl_space<Size>(lower_bounds, upper_bounds, speed_limits);
  }

  /**
   * This function constructs a rate-limited N-dof space of 1st differentiation order.
   * \return A rate-limited N-dof space of 1st differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, Size, 1> make_rl_1st_space() const {
    return make_Ndof_rl_space<Size>(lower_bounds, upper_bounds, speed_limits,
                                    accel_limits);
  }

  /**
   * This function constructs a rate-limited N-dof space of 2nd differentiation order.
   * \return A rate-limited N-dof space of 2nd differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, Size, 2> make_rl_2nd_space() const {
    return make_Ndof_rl_space<Size>(lower_bounds, upper_bounds, speed_limits,
                                    accel_limits, jerk_limits);
  }

  /**
   * This function constructs a N-dof space of 0th differentiation order.
   * \return A N-dof space of 0th differentiation order from the stored limit values.
   */
  Ndof_space_t<T, Size, 0> make_0th_space() const {
    return make_Ndof_space<Size>(lower_bounds, upper_bounds);
  }

  /**
   * This function constructs a N-dof space of 1st differentiation order.
   * \return A N-dof space of 1st differentiation order from the stored limit values.
   */
  Ndof_space_t<T, Size, 1> make_1st_space() const {
    return make_Ndof_space<Size>(lower_bounds, upper_bounds, speed_limits);
  }

  /**
   * This function constructs a N-dof space of 2nd differentiation order.
   * \return A N-dof space of 2nd differentiation order from the stored limit values.
   */
  Ndof_space_t<T, Size, 2> make_2nd_space() const {
    return make_Ndof_space<Size>(lower_bounds, upper_bounds, speed_limits,
                                 accel_limits);
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>>& pt,
      const typename Ndof_space<T, Size, 0>::type& /*unused*/,
      const typename Ndof_rl_space<T, Size, 0>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>, vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>, vect<T, Size>>& pt,
      const typename Ndof_space<T, Size, 1>::type& /*unused*/,
      const typename Ndof_rl_space<T, Size, 1>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>, vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
    };
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>>& pt,
      const typename Ndof_space<T, Size, 2>::type& /*unused*/,
      const typename Ndof_rl_space<T, Size, 2>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
      get<2>(result)[i] /= jerk_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>>& pt,
      const typename Ndof_rl_space<T, Size, 0>::type& /*unused*/,
      const typename Ndof_space<T, Size, 0>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>, vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>, vect<T, Size>>& pt,
      const typename Ndof_rl_space<T, Size, 1>::type& /*unused*/,
      const typename Ndof_space<T, Size, 1>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>, vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of rate-limited joint coordinates into a set of normal joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of normal joint coordinates corresponding to given rate-limited joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>> map_to_space(
      const arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>>& pt,
      const typename Ndof_rl_space<T, Size, 2>::type& /*unused*/,
      const typename Ndof_space<T, Size, 2>::type& /*unused*/) const {
    arithmetic_tuple<vect<T, Size>, vect<T, Size>, vect<T, Size>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
      get<2>(result)[i] /= jerk_limits[i];
    }
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(lower_bounds) &
        RK_SERIAL_SAVE_WITH_NAME(upper_bounds) &
        RK_SERIAL_SAVE_WITH_NAME(speed_limits) &
        RK_SERIAL_SAVE_WITH_NAME(accel_limits) &
        RK_SERIAL_SAVE_WITH_NAME(jerk_limits);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(lower_bounds) &
        RK_SERIAL_LOAD_WITH_NAME(upper_bounds) &
        RK_SERIAL_LOAD_WITH_NAME(speed_limits) &
        RK_SERIAL_LOAD_WITH_NAME(accel_limits) &
        RK_SERIAL_LOAD_WITH_NAME(jerk_limits);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002E, 1, "Ndof_limits", named_object)
};

/**
 * This class template stores a set of vectors to represent the rate-limits on the Ndofs
 * of a manipulator. Basically, this class is just a POD class, but it also provides functions
 * to construct a rate-limited Ndof-space from a normal Ndof-space, or vice-versa. Also,
 * it can act as a mapping between rate-limited joint coordinates and normal joint coordinates.
 * \tparam T The value type of the underlying Ndof-space.
 */
template <typename T>
struct Ndof_limits<T, 0> : public named_object {
  /** Holds the lower bounds for all generalized coordinates. */
  vect_n<T> lower_bounds;
  /** Holds the upper bounds for all generalized coordinates. */
  vect_n<T> upper_bounds;
  /** Holds the speed limit for all generalized coordinates. */
  vect_n<T> speed_limits;
  /** Holds the acceleration limit for all generalized coordinates. */
  vect_n<T> accel_limits;
  /** Holds the jerk limit for all generalized coordinates. */
  vect_n<T> jerk_limits;

  using value_type = T;
  using self = Ndof_limits<T, 0>;

  /**
   * Default constructor.
   */
  explicit Ndof_limits(const std::string& aName) : named_object() {
    this->setName(aName);
  }

  Ndof_limits() : Ndof_limits("") {}

  /**
   * This function constructs a rate-limited N-dof space of 0th differentiation order.
   * \return A rate-limited N-dof space of 0th differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, 0, 0> make_rl_0th_space() const {
    return make_Ndof_rl_space(lower_bounds, upper_bounds, speed_limits);
  }

  /**
   * This function constructs a rate-limited N-dof space of 1st differentiation order.
   * \return A rate-limited N-dof space of 1st differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, 0, 1> make_rl_1st_space() const {
    return make_Ndof_rl_space(lower_bounds, upper_bounds, speed_limits,
                              accel_limits);
  }

  /**
   * This function constructs a rate-limited N-dof space of 2nd differentiation order.
   * \return A rate-limited N-dof space of 2nd differentiation order from the stored limit values.
   */
  Ndof_rl_space_t<T, 0, 2> make_rl_2nd_space() const {
    return make_Ndof_rl_space(lower_bounds, upper_bounds, speed_limits,
                              accel_limits, jerk_limits);
  }

  /**
   * This function constructs a N-dof space of 0th differentiation order.
   * \return A N-dof space of 0th differentiation order from the stored limit values.
   */
  Ndof_space_t<T, 0, 0> make_0th_space() const {
    return make_Ndof_space(lower_bounds, upper_bounds);
  }

  /**
   * This function constructs a N-dof space of 1st differentiation order.
   * \return A N-dof space of 1st differentiation order from the stored limit values.
   */
  Ndof_space_t<T, 0, 1> make_1st_space() const {
    return make_Ndof_space(lower_bounds, upper_bounds, speed_limits);
  }

  /**
   * This function constructs a N-dof space of 2nd differentiation order.
   * \return A N-dof space of 2nd differentiation order from the stored limit values.
   */
  Ndof_space_t<T, 0, 2> make_2nd_space() const {
    return make_Ndof_space(lower_bounds, upper_bounds, speed_limits,
                           accel_limits);
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>>& pt,
      const Ndof_space_t<T, 0, 0>& /*unused*/,
      const Ndof_rl_space_t<T, 0, 0>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>, vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>, vect_n<T>>& pt,
      const Ndof_space_t<T, 0, 1>& /*unused*/,
      const Ndof_rl_space_t<T, 0, 1>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>, vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>>& pt,
      const Ndof_space_t<T, 0, 2>& /*unused*/,
      const Ndof_rl_space_t<T, 0, 2>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
      get<2>(result)[i] /= jerk_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>>& pt,
      const Ndof_rl_space_t<T, 0, 0>& /*unused*/,
      const Ndof_space_t<T, 0, 0>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of normal joint coordinates into a set of rate-limited joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of rate-limited joint coordinates corresponding to given normal joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>, vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>, vect_n<T>>& pt,
      const Ndof_rl_space_t<T, 0, 1>& /*unused*/,
      const Ndof_space_t<T, 0, 1>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>, vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
    }
    return result;
  }

  /**
   * This function maps a set of rate-limited joint coordinates into a set of normal joint coordinates.
   * \param pt A point in the normal joint-space.
   * \return A set of normal joint coordinates corresponding to given rate-limited joint coordinates and the stored
   * limit values.
   */
  arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>> map_to_space(
      const arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>>& pt,
      const Ndof_rl_space_t<T, 0, 2>& /*unused*/,
      const Ndof_space_t<T, 0, 2>& /*unused*/) const {
    arithmetic_tuple<vect_n<T>, vect_n<T>, vect_n<T>> result = pt;
    for (std::size_t i = 0; i < get<0>(result).size(); ++i) {
      get<0>(result)[i] /= speed_limits[i];
      get<1>(result)[i] /= accel_limits[i];
      get<2>(result)[i] /= jerk_limits[i];
    }
    return result;
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(lower_bounds) &
        RK_SERIAL_SAVE_WITH_NAME(upper_bounds) &
        RK_SERIAL_SAVE_WITH_NAME(speed_limits) &
        RK_SERIAL_SAVE_WITH_NAME(accel_limits) &
        RK_SERIAL_SAVE_WITH_NAME(jerk_limits);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(lower_bounds) &
        RK_SERIAL_LOAD_WITH_NAME(upper_bounds) &
        RK_SERIAL_LOAD_WITH_NAME(speed_limits) &
        RK_SERIAL_LOAD_WITH_NAME(accel_limits) &
        RK_SERIAL_LOAD_WITH_NAME(jerk_limits);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC240002E, 1, "Ndof_limits", named_object)
};

extern template struct Ndof_limits<double>;

extern template struct Ndof_limits<double, 1>;
extern template struct Ndof_limits<double, 2>;
extern template struct Ndof_limits<double, 3>;
extern template struct Ndof_limits<double, 4>;
extern template struct Ndof_limits<double, 5>;
extern template struct Ndof_limits<double, 6>;
extern template struct Ndof_limits<double, 7>;
extern template struct Ndof_limits<double, 8>;
extern template struct Ndof_limits<double, 9>;
extern template struct Ndof_limits<double, 10>;

}  // namespace ReaK::pp

#endif
