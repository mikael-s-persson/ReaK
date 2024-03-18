/**
 * \file seq_trajectory_base.h
 *
 * This library provides the base-class for sequential trajectories within a temporal topology.
 * This is a base-class that stems the object-oriented compatibility of other
 * sequential trajectory classes.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_TOPOLOGIES_INTERPOLATION_SEQ_TRAJECTORY_BASE_H_
#define REAK_TOPOLOGIES_INTERPOLATION_SEQ_TRAJECTORY_BASE_H_

#include "ReaK/core/base/defs.h"
#include "ReaK/core/base/named_object.h"

#include "ReaK/topologies/spaces/temporal_space_concept.h"

namespace ReaK::pp {

/**
 * This class defines the OOP interface for a sequential trajectory in a temporal topology.
 */
template <TemporalSpace Space>
class seq_trajectory_base : public named_object {
 public:
  using topology = Space;
  using point_type = typename topology_traits<topology>::point_type;
  using self = seq_trajectory_base<Space>;

 protected:
  struct point_time_iterator_impl {

    virtual ~point_time_iterator_impl() = default;
    ;

    virtual void move_by_time(double t) = 0;

    virtual bool is_equal_to(const point_time_iterator_impl* rhs) const = 0;

    virtual point_type get_point() const = 0;

    virtual point_time_iterator_impl* clone() const = 0;
  };

  struct point_fraction_iterator_impl {

    virtual ~point_fraction_iterator_impl() = default;
    ;

    virtual void move_by_fraction(double f) = 0;

    virtual bool is_equal_to(const point_fraction_iterator_impl* rhs) const = 0;

    virtual point_type get_point() const = 0;

    virtual point_fraction_iterator_impl* clone() const = 0;
  };

 public:
  class point_time_iterator {
   private:
    point_time_iterator_impl* p_impl;

   public:
    explicit point_time_iterator(point_time_iterator_impl* aPImpl = nullptr)
        : p_impl(aPImpl) {}

    point_time_iterator(const point_time_iterator& rhs)
        : p_impl(rhs.p_impl->clone()) {}
    point_time_iterator(point_time_iterator&& rhs) noexcept
        : p_impl(rhs.p_impl) {
      rhs.p_impl = nullptr;
    }
    friend void swap(point_time_iterator& rhs, point_time_iterator& lhs) {
      std::swap(rhs.p_impl, lhs.p_impl);
    }
    point_time_iterator& operator=(point_time_iterator rhs) {
      swap(*this, rhs);
      return *this;
    }
    ~point_time_iterator() { delete p_impl; }

    friend point_time_iterator operator+(point_time_iterator lhs, double rhs) {
      lhs.p_impl->move_by_time(rhs);
      return lhs;
    }

    friend point_time_iterator operator+(double lhs, point_time_iterator rhs) {
      rhs.p_impl->move_by_time(lhs);
      return rhs;
    }

    friend point_time_iterator& operator+=(point_time_iterator& lhs,
                                           double rhs) {
      lhs.p_impl->move_by_time(rhs);
      return lhs;
    }

    friend point_time_iterator operator-(point_time_iterator lhs, double rhs) {
      lhs.p_impl->move_by_time(-rhs);
      return lhs;
    }

    friend point_time_iterator& operator-=(point_time_iterator& lhs,
                                           double rhs) {
      lhs.p_impl->move_by_time(-rhs);
      return lhs;
    }

    friend bool operator==(const point_time_iterator& lhs,
                           const point_time_iterator& rhs) {
      return lhs.p_impl->is_equal_to(rhs.p_impl);
    }

    friend bool operator!=(const point_time_iterator& lhs,
                           const point_time_iterator& rhs) {
      return !(lhs.p_impl->is_equal_to(rhs.p_impl));
    }

    point_type operator*() const { return p_impl->get_point(); }
  };

  class point_fraction_iterator {
   private:
    point_fraction_iterator_impl* p_impl;

   public:
    explicit point_fraction_iterator(
        point_fraction_iterator_impl* aPImpl = nullptr)
        : p_impl(aPImpl) {}

    point_fraction_iterator(const point_fraction_iterator& rhs)
        : p_impl(rhs.p_impl->clone()) {}
    point_fraction_iterator(point_fraction_iterator&& rhs) noexcept
        : p_impl(rhs.p_impl) {
      rhs.p_impl = nullptr;
    }
    friend void swap(point_fraction_iterator& rhs,
                     point_fraction_iterator& lhs) {
      std::swap(rhs.p_impl, lhs.p_impl);
    }
    point_fraction_iterator& operator=(point_fraction_iterator rhs) {
      swap(*this, rhs);
      return *this;
    }
    ~point_fraction_iterator() { delete p_impl; }

    friend point_fraction_iterator operator+(point_fraction_iterator lhs,
                                             double rhs) {
      lhs.p_impl->move_by_fraction(rhs);
      return lhs;
    }

    friend point_fraction_iterator operator+(double lhs,
                                             point_fraction_iterator rhs) {
      rhs.p_impl->move_by_fraction(lhs);
      return rhs;
    }

    friend point_fraction_iterator& operator+=(point_fraction_iterator& lhs,
                                               double rhs) {
      lhs.p_impl->move_by_fraction(rhs);
      return lhs;
    }

    friend point_fraction_iterator operator-(point_fraction_iterator lhs,
                                             double rhs) {
      lhs.p_impl->move_by_fraction(-rhs);
      return lhs;
    }

    friend point_fraction_iterator& operator-=(point_fraction_iterator& lhs,
                                               double rhs) {
      lhs.p_impl->move_by_fraction(-rhs);
      return lhs;
    }

    friend bool operator==(const point_fraction_iterator& lhs,
                           const point_fraction_iterator& rhs) {
      return lhs.p_impl->is_equal_to(rhs.p_impl);
    }

    friend bool operator!=(const point_fraction_iterator& lhs,
                           const point_fraction_iterator& rhs) {
      return !(lhs.p_impl->is_equal_to(rhs.p_impl));
    }

    point_type operator*() const { return p_impl->get_point(); }
  };

  /**
   * Constructs the trajectory from a space, assumes the start and end are at the origin
   * of the space.
   * \param aName The name for this object.
   */
  explicit seq_trajectory_base(const std::string& aName = "") : named_object() {
    setName(aName);
  }

  ~seq_trajectory_base() override = default;

  /**
   * Returns the starting time-iterator along the path.
   * \return The starting time-iterator along the path.
   */
  virtual point_time_iterator begin_time_travel() const = 0;

  /**
   * Returns the end time-iterator along the path.
   * \return The end time-iterator along the path.
   */
  virtual point_time_iterator end_time_travel() const = 0;

  /**
   * Returns the starting fraction-iterator along the path.
   * \return The starting fraction-iterator along the path.
   */
  virtual point_fraction_iterator begin_fraction_travel() const = 0;

  /**
   * Returns the end fraction-iterator along the path.
   * \return The end fraction-iterator along the path.
   */
  virtual point_fraction_iterator end_fraction_travel() const = 0;

  /**
   * Computes the travel distance between two points, if traveling along the path.
   * \param a The first point.
   * \param b The second point.
   * \return The travel distance between two points if traveling along the path.
   */
  virtual double travel_distance(const point_type& a,
                                 const point_type& b) const = 0;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(self, 0xC2440014, 1, "seq_trajectory_base",
                              named_object)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_INTERPOLATION_SEQ_TRAJECTORY_BASE_H_
