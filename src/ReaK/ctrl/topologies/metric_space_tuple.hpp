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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/serializable.hpp>
#include <ReaK/core/lin_alg/arithmetic_tuple.hpp>

#include "metric_space_concept.hpp"

#include "metric_space_tuple_fwd.hpp"
#include "tuple_distance_metrics.hpp"
#include "default_random_sampler.hpp"

namespace ReaK {

namespace pp {
  
namespace detail {
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct point_type_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
  
  template <std::size_t Size, typename... Spaces>
  struct point_type_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename topology_traits<Spaces>::point_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct point_type_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename topology_traits<Spaces>::point_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct point_type_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct point_type_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_type,
                              typename topology_traits< typename arithmetic_tuple_element<9,SpaceTuple>::type >::point_type > type;
  };

#endif

  template <typename SpaceTuple>
  struct point_type_tuple : point_type_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct point_difference_type_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
  
  template <std::size_t Size, typename... Spaces>
  struct point_difference_type_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename topology_traits<Spaces>::point_difference_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct point_difference_type_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename topology_traits<Spaces>::point_difference_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct point_difference_type_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<7,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<8,SpaceTuple>::type >::point_difference_type,
                              typename topology_traits< typename arithmetic_tuple_element<9,SpaceTuple>::type >::point_difference_type > type;
  };

#endif

  template <typename SpaceTuple>
  struct point_difference_type_tuple : point_difference_type_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct topo_dimensions_tuple_impl;
  
  template <typename SpaceTuple>
  struct topo_dimensions_tuple_impl< 0, SpaceTuple > {
    BOOST_STATIC_CONSTANT(std::size_t, value = 0);
  };
  
  template <std::size_t Size, typename SpaceTuple>
  struct topo_dimensions_tuple_impl {
    BOOST_STATIC_CONSTANT(std::size_t, value = (topology_traits< typename arithmetic_tuple_element<Size-1, SpaceTuple>::type >::dimensions + topo_dimensions_tuple_impl<Size-1, SpaceTuple>::value));
  };
  
  template <typename SpaceTuple>
  struct topo_dimensions_tuple : topo_dimensions_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Order, typename SpaceTuple> 
  class metric_space_tuple_impl {
    public:
      typedef typename point_type_tuple< SpaceTuple >::type point_type;
      typedef typename point_difference_type_tuple< SpaceTuple >::type point_difference_type;
      
      static void random_point(const SpaceTuple& s, point_type& p) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::random_point(s,p);
        get<Order>(p) = get(random_sampler,get<Order>(s))(get<Order>(s));
      };
      
      static void difference(const SpaceTuple& s, point_difference_type& dp, const point_type& p1, const point_type& p2) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::difference(s,dp,p1,p2);
        get<Order>(dp) = get<Order>(s).difference(get<Order>(p1),get<Order>(p2));
      };
      
      static void move_position_toward(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::move_position_toward(s,pr,p1,d,p2);
        get<Order>(pr) = get<Order>(s).move_position_toward(get<Order>(p1),d,get<Order>(p2));
      };
      
      static void move_position_back_to(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::move_position_back_to(s,pr,p1,d,p2);
        get<Order>(pr) = get<Order>(s).move_position_back_to(get<Order>(p1),d,get<Order>(p2));
      };
      
      static void origin(const SpaceTuple& s, point_type& p) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::origin(s,p);
        get<Order>(p) = get<Order>(s).origin();
      };
      
      static void adjust(const SpaceTuple& s, point_type& pr, const point_type& p, const point_difference_type& dp) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::adjust(s,pr,p,dp);
        get<Order>(pr) = get<Order>(s).adjust(get<Order>(p),get<Order>(dp));
      };
      
      static void bring_point_in_bounds(const SpaceTuple& s, point_type& p) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::bring_point_in_bounds(s,p);
        get<Order>(s).bring_point_in_bounds(get<Order>(p));
      };
      
      static void get_diff_to_boundary(const SpaceTuple& s, point_difference_type& dpr, const point_type& p) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::get_diff_to_boundary(s,dpr,p);
        get<Order>(dpr) = get<Order>(s).get_diff_to_boundary(get<Order>(p));
      };
      
      static void is_in_bounds(const SpaceTuple& s, bool& br, const point_type& p) {
        metric_space_tuple_impl<Order-1,SpaceTuple>::is_in_bounds(s,br,p);
        br = (br && get<Order>(s).is_in_bounds(get<Order>(p)));
      };
  };
  
  template <typename SpaceTuple> 
  class metric_space_tuple_impl<0,SpaceTuple> {
    public:
      typedef typename point_type_tuple< SpaceTuple >::type point_type;
      typedef typename point_difference_type_tuple< SpaceTuple >::type point_difference_type;
      
      static void random_point(const SpaceTuple& s, point_type& p) {
        get<0>(p) = get(random_sampler,get<0>(s))(get<0>(s));
      };
      
      static void difference(const SpaceTuple& s, point_difference_type& dp, const point_type& p1, const point_type& p2) {
        get<0>(dp) = get<0>(s).difference(get<0>(p1),get<0>(p2));
      };
      
      static void move_position_toward(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
        get<0>(pr) = get<0>(s).move_position_toward(get<0>(p1),d,get<0>(p2));
      };
      
      static void move_position_back_to(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
        get<0>(pr) = get<0>(s).move_position_back_to(get<0>(p1),d,get<0>(p2));
      };
      
      static void origin(const SpaceTuple& s, point_type& p) {
        get<0>(p) = get<0>(s).origin();
      };
      
      static void adjust(const SpaceTuple& s, point_type& pr, const point_type& p, const point_difference_type& dp) {
        get<0>(pr) = get<0>(s).adjust(get<0>(p),get<0>(dp));
      };
      
      static void bring_point_in_bounds(const SpaceTuple& s, point_type& p) {
        get<0>(s).bring_point_in_bounds(get<0>(p));
      };
      
      static void get_diff_to_boundary(const SpaceTuple& s, point_difference_type& dpr, const point_type& p) {
        get<0>(dpr) = get<0>(s).get_diff_to_boundary(get<0>(p));
      };
      
      static void is_in_bounds(const SpaceTuple& s, bool& br, const point_type& p) {
        br = get<0>(s).is_in_bounds(get<0>(p));
      };
  };
  
  
  
  
  
  
  
};



  
/**
 * This class template can be used to glue together a number of spaces into a tuple. Depending on the models 
 * supported by the underlying spaces included in the tuple, this class template models 
 * the TopologyConcept, the MetricSpaceConcept, the PointDistributionConcept, the LieGroupConcept, and 
 * the BoundedSpaceConcept (i.e. the metric-space tuple class will model all the concepts which are also 
 * modeled by all the spaces it includes). This class is also a tuple class (meaning 
 * that the Boost.Tuple or std::tuple meta-functions, as well as the ReaK.Arithmetic-tuple meta-functions 
 * will work on this class, with the usual semantics).
 * 
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces to glue together.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a space-tuple (e.g. arithmetic_tuple).
 */
template <typename SpaceTuple, typename TupleDistanceMetric >
class metric_space_tuple : public shared_object {
  protected:
    SpaceTuple m_spaces;
    TupleDistanceMetric m_dist;
    
  public:
    typedef metric_space_tuple< SpaceTuple, TupleDistanceMetric > self;
    
    typedef typename detail::point_type_tuple< SpaceTuple >::type point_type;
    typedef typename detail::point_difference_type_tuple< SpaceTuple >::type point_difference_type;
    
    typedef TupleDistanceMetric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = detail::topo_dimensions_tuple<SpaceTuple>::value);
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     */
    metric_space_tuple(const SpaceTuple& aSpaces = SpaceTuple(), 
                       const TupleDistanceMetric& aDist = TupleDistanceMetric()) :
                       m_spaces(aSpaces), m_dist(aDist) { };
    
    /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
      
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
    
    friend
    TupleDistanceMetric& get(distance_metric_t,self& space) {
      return space.m_dist;
    };
    
    friend
    const TupleDistanceMetric& get(distance_metric_t,const self& space) {
      return space.m_dist;
    };
    
   /*************************************************************************
    *                             PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::random_point(m_spaces,result);
      return result;
    };
    
   /*************************************************************************
    *                             TopologyConcept
    * **********************************************************************/
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      point_difference_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::difference(m_spaces,result,p1,p2);
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

    /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& p1, double d, const point_type& p2) const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::move_position_toward(m_spaces,result,p1,d,p2);
      return result;
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_back_to(const point_type& p1, double d, const point_type& p2) const {
      point_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::move_position_back_to(m_spaces,result,p1,d,p2);
      return result;
    };
    
    
    /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Brings a given point back with the bounds of the space.
     */
    void bring_point_in_bounds(point_type& p1) const {
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::bring_point_in_bounds(m_spaces, p1);
    };
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type get_diff_to_boundary(const point_type& p1) const {
      point_difference_type result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::get_diff_to_boundary(m_spaces,result,p1);
      return result;
    };
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    bool is_in_bounds(const point_type& p1) const {
      bool result;
      detail::metric_space_tuple_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, SpaceTuple>::is_in_bounds(m_spaces,result,p1);
      return result;
    };
    
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space_impl() const {
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space_impl() {
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_impl() const {
      return get<Idx>(m_spaces);
    };
    
    /**
     * This function returns the space at a given index.
     * \tparam Idx The index of the space.
     */
    template <int Idx>
    typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_impl() {
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
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240000A,1,"metric_space_tuple",shared_object)

};


/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(const metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_space_impl<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_space_impl<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get(const metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_impl<Idx>();
};
    
/**
 * This function returns the space at a given index.
 * \tparam Idx The index of the space.
 */
template <int Idx, typename SpaceTuple, typename TupleDistanceMetric>
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get(metric_space_tuple<SpaceTuple,TupleDistanceMetric>& s) {
  return s.template get_impl<Idx>();
};




/**
 * This meta-function can be used to glue together a number of spaces of the same type into a tuple. 
 * This class will generate a metric_space_tuple class, which has N spaces of type SpaceType.
 * 
 * \tparam SpaceType The type of the spaces to glue together as a metric-space tuple.
 * \tparam N The number of spaces to glue together as a metric-space tuple.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a space-tuple (e.g. arithmetic_tuple).
 */
template <typename SpaceType, std::size_t N, typename TupleDistanceMetric >
struct metric_space_array {
  char cannot_instantiation_the_general_template[0];
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,1,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,2,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,3,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,4,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,5,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,6,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,7,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,8,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,9,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};

template <typename SpaceType, typename TupleDistanceMetric>
struct metric_space_array<SpaceType,10,TupleDistanceMetric> {
  typedef metric_space_tuple< arithmetic_tuple<SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType,
                                               SpaceType>, TupleDistanceMetric > type;
};



};

};



#endif








