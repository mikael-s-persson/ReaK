/**
 * \file seq_path_wrapper.hpp
 *
 * This library provides a sequential path-wrapper class template which makes an OOP-compatible path class
 * for a given sequential path (in the generic sense).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_SEQ_PATH_WRAPPER_HPP
#define REAK_SEQ_PATH_WRAPPER_HPP

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/interpolation/seq_path_base.hpp"
#include "ReaK/topologies/interpolation/sequential_path_concept.hpp"

#include "boost/concept_check.hpp"

namespace ReaK::pp {

/**
 * This class wraps a generic spatial path class into an OOP interface.
 * It, itself, also models the generic SpatialPathConcept, so this wrapper can
 * be used for both purposes.
 * \tparam SequentialPath The path type to be wrapped.
 */
template <typename SequentialPath>
class seq_path_wrapper
    : public seq_path_base<
          typename sequential_path_traits<SequentialPath>::topology> {
 public:
  using base_type =
      seq_path_base<typename sequential_path_traits<SequentialPath>::topology>;
  using self = seq_path_wrapper<SequentialPath>;

  using topology = typename base_type::topology;
  using point_type = typename base_type::point_type;

  BOOST_CONCEPT_ASSERT((SequentialPathConcept<SequentialPath, topology>));

  using wrapped_type = SequentialPath;

 protected:
  SequentialPath m_traj;

  using base_pt_dist_iterator_impl =
      typename base_type::point_distance_iterator_impl;
  using gen_pt_dist_iterator =
      typename sequential_path_traits<SequentialPath>::point_distance_iterator;

  struct point_distance_iterator_impl : public base_pt_dist_iterator_impl {

    gen_pt_dist_iterator base_it;

    explicit point_distance_iterator_impl(gen_pt_dist_iterator aBaseIt)
        : base_it(aBaseIt) {}

    ~point_distance_iterator_impl() override = default;

    void move_by_distance(double d) override { base_it += d; }

    bool is_equal_to(const base_pt_dist_iterator_impl* rhs) const override {
      return (base_it ==
              static_cast<const point_distance_iterator_impl*>(rhs)->base_it);
    }

    const point_type& get_point() const override { return *base_it; }

    base_pt_dist_iterator_impl* clone() const override {
      return new point_distance_iterator_impl(base_it);
    }
  };

  using base_pt_frac_iterator_impl =
      typename base_type::point_fraction_iterator_impl;
  using gen_pt_frac_iterator =
      typename sequential_path_traits<SequentialPath>::point_fraction_iterator;

  struct point_fraction_iterator_impl : public base_pt_frac_iterator_impl {

    gen_pt_frac_iterator base_it;

    explicit point_fraction_iterator_impl(gen_pt_frac_iterator aBaseIt)
        : base_it(aBaseIt) {}

    ~point_fraction_iterator_impl() override = default;

    void move_by_fraction(double f) override { base_it += f; }

    bool is_equal_to(const base_pt_frac_iterator_impl* rhs) const override {
      return (base_it ==
              static_cast<const point_fraction_iterator_impl*>(rhs)->base_it);
    }

    const point_type& get_point() const override { return *base_it; }

    base_pt_frac_iterator_impl* clone() const override {
      return new point_fraction_iterator_impl(base_it);
    }
  };

 public:
  using point_distance_iterator = typename base_type::point_distance_iterator;
  using point_fraction_iterator = typename base_type::point_fraction_iterator;

  wrapped_type& get_underlying_path() { return m_traj; }
  const wrapped_type& get_underlying_path() const { return m_traj; }

  wrapped_type& get_wrapped_object() { return m_traj; }
  const wrapped_type& get_wrapped_object() const { return m_traj; }

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aName The name for this object.
   * \param aTraj The wrapped path object to use.
   */
  explicit seq_path_wrapper(const std::string& aName = "",
                            const SequentialPath& aTraj = SequentialPath())
      : base_type(aName), m_traj(aTraj) {}

  ~seq_path_wrapper() override = default;

  /**
   * Returns the starting distance-iterator along the path.
   * \return The starting distance-iterator along the path.
   */
  point_distance_iterator begin_distance_travel() const override {
    return point_distance_iterator(
        new point_distance_iterator_impl(m_traj.begin_distance_travel()));
  }

  /**
   * Returns the end distance-iterator along the path.
   * \return The end distance-iterator along the path.
   */
  point_distance_iterator end_distance_travel() const override {
    return point_distance_iterator(
        new point_distance_iterator_impl(m_traj.end_distance_travel()));
  }

  /**
   * Returns the starting fraction-iterator along the path.
   * \return The starting fraction-iterator along the path.
   */
  point_fraction_iterator begin_fraction_travel() const override {
    return point_fraction_iterator(
        new point_fraction_iterator_impl(m_traj.begin_fraction_travel()));
  }

  /**
   * Returns the end fraction-iterator along the path.
   * \return The end fraction-iterator along the path.
   */
  point_fraction_iterator end_fraction_travel() const override {
    return point_fraction_iterator(
        new point_fraction_iterator_impl(m_traj.end_fraction_travel()));
  }

  double travel_distance(const point_type& a,
                         const point_type& b) const override {
    return m_traj.travel_distance(a, b);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    base_type::save(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(m_traj);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    base_type::load(A, base_type::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(m_traj);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2440012, 1, "seq_path_wrapper",
                              base_type)
};

}  // namespace ReaK::pp

#endif
