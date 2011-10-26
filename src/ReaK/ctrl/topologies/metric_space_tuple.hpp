/**
 * \file metric_space_tuple.hpp
 * 
 * This library provides classes that define a metric-space tuple class template. A metric-space tuple is 
 * a simple association of several topologies (metric-spaces) which, in turn, also models a metric-space
 * (conditional upon each underlying spaces being a metric-space as well).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_METRIC_SPACE_TUPLE_HPP
#define REAK_METRIC_SPACE_TUPLE_HPP

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "path_planning/metric_space_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "base/serializable.hpp"
#include "tuple_distance_metrics.hpp"

namespace ReaK {

namespace pp {
  
namespace detail {
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct point_type_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct point_type_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct point_type_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct point_type_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<9,SpaceTuple>::type >::point_type > type;
  };

#endif

  template <typename SpaceTuple>
  struct point_type_tuple : point_type_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct point_difference_type_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct point_difference_type_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_difference_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct point_difference_type_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_difference_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<9,SpaceTuple>::type >::point_difference_type > type;
  };

#endif

  template <typename SpaceTuple>
  struct point_difference_type_tuple : point_difference_type_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Order, typename SpaceTuple> 
  class metric_space_tuple_impl {
    public:
      typedef typename point_type_tuple< SpaceTuple >::type point_type;
      typedef typename point_difference_type_tuple< SpaceTuple >::type point_difference_type;
      
      static void random_point(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	metric_space_tuple_impl<Order-1,SpaceTuple>::random_point(s,p);
	get<Order>(p) = get<Order>(s).random_point();
      };
      
      static void difference(const SpaceTuple& s, point_difference_type& dp, const point_type& p1, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	metric_space_tuple_impl<Order-1,SpaceTuple>::difference(s,dp,p1,p2);
	get<Order>(dp) = get<Order>(s).difference(get<Order>(p1),get<Order>(p2));
      };
      
      static void move_position_toward(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	metric_space_tuple_impl<Order-1,SpaceTuple>::move_position_toward(s,pr,p1,d,p2);
	get<Order>(pr) = get<Order>(s).move_position_toward(get<Order>(p1),d,get<Order>(p2));
      };
      
      static void origin(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	metric_space_tuple_impl<Order-1,SpaceTuple>::origin(s,p);
	get<Order>(p) = get<Order>(s).origin();
      };
      
      static void adjust(const SpaceTuple& s, point_type& pr, const point_type& p, const point_difference_type& dp) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	metric_space_tuple_impl<Order-1,SpaceTuple>::adjust(s,pr,p,dp);
	get<Order>(pr) = get<Order>(s).adjust(get<Order>(p),get<Order>(dp));
      };
  };
  
  template <typename SpaceTuple> 
  class metric_space_tuple_impl<0,SpaceTuple> {
    public:
      typedef typename point_type_tuple< SpaceTuple >::type point_type;
      typedef typename point_difference_type_tuple< SpaceTuple >::type point_difference_type;
      
      static void random_point(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	get<0>(p) = get<0>(s).random_point();
      };
      
      static void difference(const SpaceTuple& s, point_difference_type& dp, const point_type& p1, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	get<0>(dp) = get<0>(s).difference(get<0>(p1),get<0>(p2));
      };
      
      static void move_position_toward(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	get<0>(pr) = get<0>(s).move_position_toward(get<0>(p1),d,get<0>(p2));
      };
      
      static void origin(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	get<0>(p) = get<0>(s).origin();
      };
      
      static void adjust(const SpaceTuple& s, point_type& pr, const point_type& p, const point_difference_type& dp) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	get<0>(pr) = get<0>(s).adjust(get<0>(p),get<0>(dp));
      };
  };
  
  
  
  
  
  
  
};



  
/**
 * This class template can be used to glue together a number of spaces into a tuple. This class template models 
 * the MetricSpaceConcept (if all underlying spaces do as well), and it is also a tuple class (meaning 
 * that the Boost.Tuple or std::tuple meta-functions, as well as the ReaK.Arithmetic-tuple meta-functions 
 * will work on this class, with the usual semantics).
 * 
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces to glue together.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a space-tuple (e.g. arithmetic_tuple).
 */
template <typename SpaceTuple, typename TupleDistanceMetric = manhattan_tuple_distance >
class metric_space_tuple : public serialization::serializable {
  protected:
    SpaceTuple m_spaces;
    TupleDistanceMetric m_dist;
    
  public:
    typedef metric_space_tuple< SpaceTuple, TupleDistanceMetric > self;
    
    typedef typename detail::point_type_tuple< SpaceTuple >::type point_type;
    typedef typename detail::point_difference_type_tuple< SpaceTuple >::type point_difference_type;
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     */
    metric_space_tuple(const SpaceTuple& aSpaces = SpaceTuple(), 
		       const TupleDistanceMetric& aDist = TupleDistanceMetric()) :
		       m_spaces(aSpaces), m_dist(aDist) { };
      
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      return m_dist(p1, p2, m_spaces);
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& dp) const {
      return m_dist(dp, m_spaces);
    };
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::random_point(m_spaces,result);
      return result;
    };
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      point_difference_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::difference(m_spaces,result,p1,p2);
      return result;
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& p1, double d, const point_type& p2) const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::move_position_toward(m_spaces,result,p1,d,p2);
      return result;
    };
    
    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::origin(m_spaces,result);
      return result;
    };
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& p1, const point_difference_type& dp) const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::adjust(m_spaces,result,p1,dp);
      return result;
    };
    
    
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space() const {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::tuples::get;
#endif
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space() {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::tuples::get;
#endif
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get() const {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::tuples::get;
#endif
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get() {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::tuples::get;
#endif
      return get<Idx>(m_spaces);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      A & RK_SERIAL_SAVE_WITH_NAME(m_spaces)
        & RK_SERIAL_SAVE_WITH_NAME(m_dist);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      A & RK_SERIAL_LOAD_WITH_NAME(m_spaces)
        & RK_SERIAL_LOAD_WITH_NAME(m_dist);
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC240000A,1,"metric_space_tuple",serialization::serializable)

};


/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(const metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_space<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_space<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get(const metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get(metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get<Idx>();
};
    




};

};


#ifdef RK_ENABLE_CXX0X_FEATURES

namespace std {
  
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_size< ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    public tuple_size< SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_size< const ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    public tuple_size< const SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_size< volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    public tuple_size< volatile SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_size< const volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    public tuple_size< const volatile SpaceTuple > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_element< Idx, ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename tuple_element< Idx, SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_element< Idx, const ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename tuple_element< Idx, const SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_element< Idx, volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename tuple_element< Idx, volatile SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class tuple_element< Idx, const volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename tuple_element< Idx, const volatile SpaceTuple >::type type;
  };
  
};

#else

namespace boost {
  
namespace tuples {
  
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  struct length< ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    ReaK::arithmetic_tuple_size< SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  struct length< const ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    ReaK::arithmetic_tuple_size< SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class length< volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    ReaK::arithmetic_tuple_size< SpaceTuple > { };
    
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  class length< const volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    ReaK::arithmetic_tuple_size< SpaceTuple > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class element< Idx, ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename element< Idx, SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class element< Idx, const ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename element< Idx, const SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class element< Idx, volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename element< Idx, volatile SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class element< Idx, const volatile ReaK::pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename element< Idx, const volatile SpaceTuple >::type type;
  };
  
};
  
};

#endif


namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_size< pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > : 
    arithmetic_tuple_size< SpaceTuple > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class arithmetic_tuple_element< Idx, pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename arithmetic_tuple_element< Idx, SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class arithmetic_tuple_element< Idx, const pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename arithmetic_tuple_element< Idx, const SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class arithmetic_tuple_element< Idx, volatile pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename arithmetic_tuple_element< Idx, volatile SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
  class arithmetic_tuple_element< Idx, const volatile pp::metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
    public:
      typedef typename arithmetic_tuple_element< Idx, const volatile SpaceTuple >::type type;
  };
  
};

#endif








