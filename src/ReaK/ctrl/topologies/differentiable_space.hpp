/**
 * \file differentiable_space.hpp
 * 
 * This library provides classes that define a differentiable-space. A differentiable-space is 
 * a simple association of two topologies to relate them by a derivative-integral relationship.
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

#ifndef REAK_DIFFERENTIABLE_SPACE_HPP
#define REAK_DIFFERENTIABLE_SPACE_HPP

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include <cmath>
#include "time_topology.hpp"
#include "path_planning/metric_space_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include <core/lin_alg/arithmetic_tuple.hpp>

namespace ReaK {

namespace pp {
  

struct default_differentiation_rule {
  template <typename T, typename U, typename V>
  static void lift(T& v, const U& dp, const V& dt) const {
    v = dp / dt;
  };
  template <typename T, typename U, typename V>
  static void descend(T& dp, const U& v, const V& dt) const {
    dp = v * dt;
  };
};


namespace detail {
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct differentiable_point_tuple_impl { 
    BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct differentiable_point_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct differentiable_point_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_tuple_impl< 8, SpaceTuple > { 
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
  struct differentiable_point_tuple_impl< 9, SpaceTuple > { 
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
  struct differentiable_point_tuple_impl< 10, SpaceTuple > { 
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
  struct differentiable_point_tuple : differentiable_point_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl { 
    BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct differentiable_point_difference_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_difference_type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct differentiable_point_difference_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename metric_topology_traits<Spaces>::point_difference_type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename metric_topology_traits< typename arithmetic_tuple_element<0,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<1,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<2,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<3,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<4,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<5,SpaceTuple>::type >::point_difference_type,
                              typename metric_topology_traits< typename arithmetic_tuple_element<6,SpaceTuple>::type >::point_difference_type > type;
  };

  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple_impl< 8, SpaceTuple > { 
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
  struct differentiable_point_difference_tuple_impl< 9, SpaceTuple > { 
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
  struct differentiable_point_difference_tuple_impl< 10, SpaceTuple > { 
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
  struct differentiable_point_difference_tuple : differentiable_point_difference_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  
  
  
  template <std::size_t Order, typename SpaceTuple> 
  class differentiable_space_impl {
    public:
      typedef typename differentiable_point_tuple< SpaceTuple >::type point_type;
      typedef typename differentiable_point_difference_tuple< SpaceTuple >::type point_difference_type;
      
      static double distance(const SpaceTuple& s, const point_type& p1, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	double result = differentiable_space_impl<Order-1,SpaceTuple>::distance(s,p1,p2);
	result += get<Order>(s).distance(get<Order>(p1),get<Order>(p2));
	return result;
      };
      
      static double norm(const SpaceTuple& s, const point_difference_type& dp) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	double result = differentiable_space_impl<Order-1,SpaceTuple>::norm(s,dp);
	result += get<Order>(s).norm(get<Order>(dp));
	return result;
      };
      
      static void random_point(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	differentiable_space_impl<Order-1,SpaceTuple>::random_point(s,p);
	get<Order>(p) = get<Order>(s).random_point();
      };
      
      static void difference(const SpaceTuple& s, point_difference_type& dp, const point_type& p1, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	differentiable_space_impl<Order-1,SpaceTuple>::difference(s,dp,p1,p2);
	get<Order>(dp) = get<Order>(s).difference(get<Order>(p1),get<Order>(p2));
      };
      
      static void move_position_toward(const SpaceTuple& s, point_type& pr, const point_type& p1, double d, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	differentiable_space_impl<Order-1,SpaceTuple>::move_position_toward(s,pr,p1,d,p2);
	get<Order>(pr) = get<Order>(s).move_position_toward(get<Order>(p1),d,get<Order>(p2));
      };
      
      static void origin(const SpaceTuple& s, point_type& p) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	differentiable_space_impl<Order-1,SpaceTuple>::origin(s,p);
	get<Order>(p) = get<Order>(s).origin();
      };
      
      static void adjust(const SpaceTuple& s, point_type& pr, const point_type& p, const point_difference_type& dp) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	differentiable_space_impl<Order-1,SpaceTuple>::adjust(s,pr,p,dp);
	get<Order>(pr) = get<Order>(s).adjust(get<Order>(p),get<Order>(dp));
      };
  };
  
  template <typename SpaceTuple> 
  class differentiable_space_impl<0,SpaceTuple> {
    public:
      typedef typename differentiable_point_tuple< SpaceTuple >::type point_type;
      typedef typename differentiable_point_difference_tuple< SpaceTuple >::type point_difference_type;
      
      static double distance(const SpaceTuple& s, const point_type& p1, const point_type& p2) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	return get<0>(s).distance(get<0>(p1),get<0>(p2));
      };
      
      static double norm(const SpaceTuple& s, const point_difference_type& dp) {
#ifdef RK_ENABLE_CXX0X_FEATURES
        using std::get;
#else
        using boost::tuples::get;
#endif
	return get<0>(s).norm(get<0>(dp));
      };
      
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



  
template <typename IndependentSpace, typename SpaceTuple, typename DiffRuleTuple>
class differentiable_space {
  protected:
    SpaceTuple m_space;
    
  public:
    template <int Idx>
    struct space {
      typedef typename arithmetic_tuple_element<Idx, SpaceTuple>::type type;
    };
      
    BOOST_STATIC_CONSTANT(std::size_t, differential_order = arithmetic_tuple_size<SpaceTuple>::type::value);
      
    typedef typename differentiable_point_tuple< SpaceTuple >::type point_type;
    typedef typename differentiable_point_difference_tuple< SpaceTuple >::type point_difference_type;
      
    double distance(const point_type& p1, const point_type& p2) const {
      return detail::differentiable_space_impl<differential_order, SpaceTuple>::distance(m_space, p1, p2);
    };
    
    double norm(const point_difference_type& dp) const {
      return detail::differentiable_space_impl<differential_order, SpaceTuple>::norm(m_space, dp);
    };
    
    point_type random_point() const {
      point_type result;
      detail::differentiable_space_impl<differential_order, SpaceTuple>::random_point(m_space,result);
      return result;
    };
    
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      point_difference_type result;
      detail::differentiable_space_impl<differential_order, SpaceTuple>::difference(m_space,result,p1,p2);
      return result;
    };
    
    point_type move_position_toward(const point_type& p1, double d, const point_type& p2) const {
      point_type result;
      detail::differentiable_space_impl<differential_order, SpaceTuple>::move_position_toward(m_space,result,p1,d,p2);
      return result;
    };
    
    point_type origin() const {
      point_type result;
      detail::differentiable_space_impl<differential_order, SpaceTuple>::origin(m_space,result);
      return result;
    };
    
    point_type adjust(const point_type& p1, const point_difference_type& dp) const {
      point_type result;
      detail::differentiable_space_impl<differential_order, SpaceTuple>::adjust(m_space,result,p1,dp);
      return result;
    };
      
    
    template <int Idx>
    const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(const IndependentSpace&) const {
#ifdef RK_ENABLE_CXX0X_FEATURES
      using std::get;
#else
      using boost::tuples::get;
#endif
      return get<Idx>(m_space);
    };
    
    template <int Idx>
    typename arithmetic_tuple_element<Idx, point_type>::type 
      lift_to_space(const typename arithmetic_tuple_element<Idx-1, point_difference_type>::type& dp,
		    const typename metric_topology_traits< IndependentSpace >::point_difference_type& dt) {
      typename arithmetic_tuple_element<Idx, point_type>::type result;
      arithmetic_tuple_element<Idx-1, DiffRuleTuple >::type::lift(result, dp, dt);
      return result;
    };
    
    template <int Idx>
    typename arithmetic_tuple_element<Idx, point_difference_type>::type 
      descend_to_space(const typename arithmetic_tuple_element<Idx+1, point_type>::type& v,
		       const typename metric_topology_traits< IndependentSpace >::point_difference_type& dt) {
      typename arithmetic_tuple_element<Idx, point_difference_type>::type result;
      arithmetic_tuple_element<Idx, DiffRuleTuple >::type::descend(result, v, dt);
      return result;
    };
    
};




};

};

#endif








