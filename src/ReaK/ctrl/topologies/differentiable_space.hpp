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

#ifdef RK_ENABLE_CXX0X_FEATURES
#include <tuple>
#else
#include <boost/tuple/tuple.hpp>
#endif

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
  
  
  template <typename SpaceTuple>
  struct differentiable_point_tuple { 
    BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  template <typename... Spaces>
  struct differentiable_point_tuple< std::tuple<Spaces...> > {
    typedef std::tuple< typename metric_topology_traits<Spaces>::point_type... > type;
  };
  
#else
  template <typename T1>
  struct differentiable_point_tuple< boost::tuples::tuple<T1> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type > type;
  };
  
  template <typename T1, typename T2>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type,
                                  typename metric_topology_traits<T6>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type,
                                  typename metric_topology_traits<T6>::point_type,
                                  typename metric_topology_traits<T7>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type,
                                  typename metric_topology_traits<T6>::point_type,
                                  typename metric_topology_traits<T7>::point_type,
                                  typename metric_topology_traits<T8>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type,
                                  typename metric_topology_traits<T6>::point_type,
                                  typename metric_topology_traits<T7>::point_type,
                                  typename metric_topology_traits<T8>::point_type,
                                  typename metric_topology_traits<T9>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  struct differentiable_point_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type,
                                  typename metric_topology_traits<T4>::point_type,
                                  typename metric_topology_traits<T5>::point_type,
                                  typename metric_topology_traits<T6>::point_type,
                                  typename metric_topology_traits<T7>::point_type,
                                  typename metric_topology_traits<T8>::point_type,
                                  typename metric_topology_traits<T9>::point_type,
                                  typename metric_topology_traits<T10>::point_type > type;
  };
#endif

  
  
  template <typename SpaceTuple>
  struct differentiable_point_difference_tuple { 
    BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  template <typename... Spaces>
  struct differentiable_point_difference_tuple< std::tuple<Spaces...> > {
    typedef std::tuple< typename metric_topology_traits<Spaces>::point_difference_type... > type;
  };
  
#else
  template <typename T1>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type > type;
  };
  
  template <typename T1, typename T2>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_type,
                                  typename metric_topology_traits<T2>::point_type,
                                  typename metric_topology_traits<T3>::point_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type,
                                  typename metric_topology_traits<T6>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type,
                                  typename metric_topology_traits<T6>::point_difference_type,
                                  typename metric_topology_traits<T7>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type,
                                  typename metric_topology_traits<T6>::point_difference_type,
                                  typename metric_topology_traits<T7>::point_difference_type,
                                  typename metric_topology_traits<T8>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type,
                                  typename metric_topology_traits<T6>::point_difference_type,
                                  typename metric_topology_traits<T7>::point_difference_type,
                                  typename metric_topology_traits<T8>::point_difference_type,
                                  typename metric_topology_traits<T9>::point_difference_type > type;
  };
  
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  struct differentiable_point_difference_tuple< boost::tuples::tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { 
    typedef boost::tuples::tuple< typename metric_topology_traits<T1>::point_difference_type,
                                  typename metric_topology_traits<T2>::point_difference_type,
                                  typename metric_topology_traits<T3>::point_difference_type,
                                  typename metric_topology_traits<T4>::point_difference_type,
                                  typename metric_topology_traits<T5>::point_difference_type,
                                  typename metric_topology_traits<T6>::point_difference_type,
                                  typename metric_topology_traits<T7>::point_difference_type,
                                  typename metric_topology_traits<T8>::point_difference_type,
                                  typename metric_topology_traits<T9>::point_difference_type,
                                  typename metric_topology_traits<T10>::point_difference_type > type;
  };
#endif

#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <typename IndependentSpace, typename SpaceTuple, typename DiffRuleTuple>
  class differentiable_space_impl {
    public:
      template <unsigned int Idx>
      struct space {
        typedef typename std::tuple_element<Idx, SpaceTuple>::type type;
      };
      
      BOOST_STATIC_CONSTANT(unsigned int, differential_order = std::tuple_size<SpaceTuple>::value);
      
      class point_type {
        public:
	  typedef typename differentiable_point_tuple< SpaceTuple >::type value_type;
        private:
	  value_type value;
        public:
	  point_type(const value_type& aValue) : value(aValue) { };
	  template <typename... Args>
	  point_type(Args&&... args) : value(value_type(std::forward<Args>(args)...)) { };
      };
      
      class point_difference_type {
        public:
          typedef typename differentiable_point_difference_tuple< SpaceTuple >::type value_type;
        private:
	  value_type value;
        public:
	  point_difference_type(const value_type& aValue) : value(aValue) { };
	  template <typename... Args>
	  point_difference_type(Args&&... args) : value(value_type(std::forward<Args>(args)...)) { };
	  
	  
      };
      
  };
  
  
#else
  
  template <typename IndependentSpace, typename SpaceTuple, typename DiffRuleTuple>
  class differentiable_space_impl {
    public:
      template <unsigned int Idx>
      struct space {
        typedef typename boost::tuples::element<Idx, SpaceTuple>::type type;
      };
      
      BOOST_STATIC_CONSTANT(unsigned int, differential_order = boost::tuples::length<SpaceTuple>::value);
      
      
      
  };
  
  
  
#endif
  
  
  template <typename IndependentSpace, typename SpaceTuple, typename DiffRuleTuple>
  class differentiable_space_impl {
    public:
      template <unsigned int Idx>
      struct space {
#ifdef RK_ENABLE_CXX0X_FEATURES
        typedef typename std::tuple_element<Idx, SpaceTuple>::type type;
#else
        typedef typename boost::tuples::element<Idx, SpaceTuple>::type type;
#endif
      };
      
#ifdef RK_ENABLE_CXX0X_FEATURES
      BOOST_STATIC_CONSTANT(unsigned int, differential_order = std::tuple_size<SpaceTuple>::value);
#else
      BOOST_STATIC_CONSTANT(unsigned int, differential_order = boost::tuples::length<SpaceTuple>::value);
#endif
      
      
      
  };
  
  
  
  
};



template <typename FirstOrderSpace, typename SecondOrderSpace, typename IndependentSpace = time_topology, typename DifferentiationRule = default_differentiation_rule>
class differentiable_space {
  public:
    
    
    
};




};

};

#endif








