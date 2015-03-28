/**
 * \file reversible_space_concept.hpp
 *
 * This library defines the traits and concepts that pertain to what can be considered
 * a reversible-space, as used in ReaK::pp. Reversible-spaces are based on the Topology concept
 * from the Boost.Graph library, but with additional requirements which are needed
 * in algorithms tailored for path-planning (mostly). Basically, the concept of a
 * reversible-space corresponds to spaces in which motions can
 * be done in reverse, whether steerable or not.
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

#ifndef REAK_REVERSIBLE_SPACE_CONCEPT_HPP
#define REAK_REVERSIBLE_SPACE_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include <cmath>
#include <boost/concept_check.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/tuple/tuple.hpp>

#include "metric_space_concept.hpp"
#include "steerable_space_concept.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {

namespace detail {

template < bool IsMetric >
struct move_backwards_concept_tester {
  template < typename Space, typename PointType >
  static void test( const Space&, const PointType&, PointType& ){};
};

template <>
struct move_backwards_concept_tester< true > {
  template < typename Space, typename PointType >
  static void test( const Space& s, const PointType& p1, PointType& p2 ) {
    p2 = s.move_position_back_to( p1, 1.0, p1 );
  };
};

template < bool IsSteerable >
struct steer_backwards_concept_tester {
  template < typename Space, typename PointType >
  static void test( const Space&, const PointType&, PointType& ){};
};

template <>
struct steer_backwards_concept_tester< true > {
  template < typename Space, typename PointType >
  static void test( const Space& s, const PointType& p1, PointType& p2 ) {
    typename steerable_space_traits< Space >::steer_record_type st_rec;
    boost::tie( p2, st_rec ) = s.steer_position_back_to( p1, 1.0, p1 );
  };
};
};


/**
 * This concept defines the requirements to fulfill in order to model a reversible-space
 * as used in ReaK::pp. A reversible-space is a special kind of topology in which
 * motions can be done in reverse, whether steerable or not.
 *
 * Valid expressions:
 *
 * if(is_metric_space<ReversibleSpace>)
 *   p2 = space.move_position_back_to(p1, d, p3);  A point (p2) can be obtained by moving a fraction (d) back from a
 *final point (p3) to previous point (p1).
 *
 * if(is_steerable_space<ReversibleSpace>)
 *   tie(p2, st_rec) = space.steer_position_back_to(p1, d, p3);  A point (p2) and a steer-record (st_rec) can be
 *obtained from attempting to steer a fraction (d) back from a final point (p3) to previous point (p1).
 *
 * \tparam ReversibleSpace The topology type to be checked for this concept.
 */
template < typename ReversibleSpace >
struct ReversibleSpaceConcept {
  typename topology_traits< ReversibleSpace >::point_type p1, p2;
  ReversibleSpace space;

  BOOST_CONCEPT_ASSERT( ( TopologyConcept< ReversibleSpace > ) );

  BOOST_CONCEPT_USAGE( ReversibleSpaceConcept ) {
    detail::move_backwards_concept_tester< is_metric_space< ReversibleSpace >::type::value >::test( space, p1, p2 );
    detail::steer_backwards_concept_tester< is_steerable_space< ReversibleSpace >::type::value >::test( space, p1, p2 );
  };
};


template < typename ReversibleSpace >
struct is_reversible_space : boost::mpl::false_ {};
};
};


#endif
