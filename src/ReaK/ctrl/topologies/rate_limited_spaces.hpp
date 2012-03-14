/**
 * \file rate_limited_spaces.hpp
 * 
 * This library provides classes that define a number of reach-time-space class templates. A reach-time-space is 
 * a transformation on a tangent bundle such that points in the spaces and their differences
 * represent the time-to-reach (or reach-time) based on the N-order rate limitations on the tangent spaces.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_RATE_LIMITED_SPACES_HPP
#define REAK_RATE_LIMITED_SPACES_HPP

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/bounded_space_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "base/serializable.hpp"
#include "tuple_distance_metrics.hpp"
#include "time_topology.hpp"
#include "default_random_sampler.hpp"
#include "differentiable_space.hpp"
#include <lin_alg/vect_alg.hpp>

namespace ReaK {

namespace pp {


/**
 * This class defines the differentiation rule to apply either to lift a 
 * point-difference (e.g. finite-difference) to the tangent space, or to descend 
 * a tangent vector to a point-difference, for topologies whose point-difference 
 * vectors are expressed as reach-time values.
 */
struct reach_time_differentiation : public serialization::serializable {
  
  double max_rate_reach_time;
  
  reach_time_differentiation(double aMaxRateReachTime = 1.0) : max_rate_reach_time(aMaxRateReachTime) { };
  
  /**
   * This function will lift a point-difference vector into its corresponding tangent vector.
   * This function performs a simple division, dp * (max_rate_reach_time / dt).
   * \tparam T The destination type, a point in the tangent space.
   * \tparam U The source type, a point-difference in the base space.
   * \tparam V A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param v The resulting point in the tangent space.
   * \param dp The point-difference that is being lifted.
   * \param dt The time-difference value (i.e. the difference in the independent variable). 
   */
  template <typename T, typename U, typename V, typename TSpace>
  void lift(T& v, const U& dp, const V& dt, const TSpace&) const {
    v = dp * (max_rate_reach_time / dt);
  };
  /**
   * This function will descend a tangent vector into its corresponding point-difference vector.
   * This function performs a simple multiplication, v * (dt / max_rate_reach_time).
   * \tparam T The destination type, a point-difference in the base space.
   * \tparam U The source type, a point in the tangent space.
   * \tparam V A type representing the independent variable's difference (e.g. time-difference).
   * \tparam TSpace The type of the independent space (e.g. time-space).
   * \param dp The resulting point-difference in the base space.
   * \param v The point in the tangent space that is being descended.
   * \param dt The time-difference value (i.e. the difference in the independent variable). 
   */
  template <typename T, typename U, typename V, typename TSpace>
  void descend(T& dp, const U& v, const V& dt, const TSpace&) const {
    dp = v * (dt / max_rate_reach_time);
  };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { 
    A & RK_SERIAL_SAVE_WITH_NAME(max_rate_reach_time);
  };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { 
    A & RK_SERIAL_LOAD_WITH_NAME(max_rate_reach_time);
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE(reach_time_differentiation,0xC2420001,1,"reach_time_differentiation",serialization::serializable)
};




/**
 * This class defines a tuple of reach-time differentiation rules. This is useful for applying the 
 * reach-time differentiation rule to all the differentiation operations on a given differentiable_space.
 */
template <std::size_t Order>
struct reach_time_differentiation_tuple {
  //BOOST_STATIC_ASSERT(false);
};

template <>
struct reach_time_differentiation_tuple<0> {
  typedef arithmetic_tuple<reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<1> {
  typedef arithmetic_tuple<reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<2> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<3> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<4> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<5> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<6> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<7> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<8> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<9> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};

template <>
struct reach_time_differentiation_tuple<10> {
  typedef arithmetic_tuple<reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation,
                           reach_time_differentiation> type;
};


namespace detail {
  
  
  template <typename Idx, typename ForwardIter, typename ReachTimeDiffTuple>
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      arithmetic_tuple_size< ReachTimeDiffTuple >
    >,
  void >::type assign_reach_time_diff_rules(ForwardIter it, ReachTimeDiffTuple& diff_rule) {
    get< Idx::type::value >(diff_rule) = *it;
    assign_reach_time_diff_rules< boost::mpl::next<Idx> >(++it, diff_rule);
  };
  
  template <typename Idx, typename ForwardIter, typename ReachTimeDiffTuple>
  typename boost::disable_if<
    boost::mpl::less< 
      Idx,
      arithmetic_tuple_size< ReachTimeDiffTuple >
    >,
  void >::type assign_reach_time_diff_rules(ForwardIter, ReachTimeDiffTuple&) { 
    /* stop TMP-recursion */
  };
  
  
  template <typename ReachTimeDiffTuple, typename ForwardIter>
  ReachTimeDiffTuple construct_reach_time_diff_rules(ForwardIter first, ForwardIter last) {
    ReachTimeDiffTuple result;
    assign_reach_time_diff_rules< boost::mpl::size_t< 0 > >(first, result);
    return result;
  };
  
  
};




/**
 * This class template can be used to glue together a number of spaces by a reach-time differentiation 
 * relationship, where each differentiation / integration operation (or more formally speaking, each 
 * lift and descent through the tangent bundle) is governed by its own reach-time differentiation rule. 
 * This class template models the TopologyConcept (if all underlying spaces do as well), and models 
 * the TangentBundleConcept for as high an order as there are differentiation rules and spaces to support.
 * 
 * \tparam IndependentSpace The type of the independent-space against which the differentiation is 
 *                          taken (e.g. time_topology). There are no formal requirements on this type, 
 *                          it is merely used as a placeholder by this class (although the differentiation 
 *                          rules might require more of this type).
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces that are arranged 
 *                    in sequence of differentiation levels (e.g. space 0 -- diff --> space 1 -- diff --> space 2 ...).
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a 
 *                             space-tuple (e.g. arithmetic_tuple).
 */
template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric = manhattan_tuple_distance>
class reach_time_diff_space : public differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,typename reach_time_differentiation_tuple< arithmetic_tuple_size<SpaceTuple>::type::value - 1 >::type> {
  public:
    typedef reach_time_diff_space< IndependentSpace, SpaceTuple, TupleDistanceMetric > self;
    typedef typename reach_time_differentiation_tuple< arithmetic_tuple_size<SpaceTuple>::type::value - 1 >::type DiffRuleTuple;
    typedef differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> base_type;
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     * \param aDiffRules The differentiation rule tuple to initialize the diff-rule functors with.
     */
    reach_time_diff_space(const SpaceTuple& aSpaces = SpaceTuple(), 
			  const TupleDistanceMetric& aDist = TupleDistanceMetric(),
			  const DiffRuleTuple& aDiffRules = DiffRuleTuple()) :
			  base_type(aSpaces,aDist,aDiffRules) { };

    /**
     * Parametrized constructor which creates the reach-time differentiation rules from a list 
     * of rate-limits on spaces of the tangent bundle (element 0 should be the velocity limit, 
     * then the acceleration limit, and so on, so forth).
     * \param first An iterator to the first rate-limit used for differentiation rules.
     * \param last An iterator to the one-past-last rate-limit used for differentiation rules.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     */
    template <typename ForwardIter>
    reach_time_diff_space(ForwardIter first, ForwardIter last,
                          const SpaceTuple& aSpaces = SpaceTuple(), 
			  const TupleDistanceMetric& aDist = TupleDistanceMetric()) :
			  base_type(aSpaces, aDist, detail::construct_reach_time_diff_rules<DiffRuleTuple>(first, last)) { };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    /**
     * Parametrized constructor which creates the reach-time differentiation rules from a list 
     * of rate-limits on spaces of the tangent bundle (element 0 should be the velocity limit, 
     * then the acceleration limit, and so on, so forth).
     * \param aList An initializer-list to the rate-limits used for differentiation rules.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     */
    reach_time_diff_space(std::initializer_list<double> aList,
                          const SpaceTuple& aSpaces = SpaceTuple(), 
			  const TupleDistanceMetric& aDist = TupleDistanceMetric()) :
			  base_type(aSpaces, aDist, detail::construct_reach_time_diff_rules<DiffRuleTuple>(aList.begin(), aList.end())) { };
#endif
    
			  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400010,1,"reach_time_diff_space",base_type)

};






template <typename RateLimitedSpace>
struct get_rate_illimited_space { 
  typedef RateLimitedSpace type;
};

template <typename NormalSpace>
struct get_rate_limited_space {  
  typedef NormalSpace type;
};



namespace detail {
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct get_rate_illimited_space_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename get_rate_illimited_space<Spaces>::type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct get_rate_illimited_space_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename get_rate_illimited_space<Spaces>::type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<8,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_illimited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<8,SpaceTuple>::type >::type,
                              typename get_rate_illimited_space< typename arithmetic_tuple_element<9,SpaceTuple>::type >::type > type;
  };

#endif

  template <typename SpaceTuple>
  struct get_rate_illimited_space_tuple : get_rate_illimited_space_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
  template <std::size_t Size, typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, typename... Spaces>
  struct get_rate_limited_space_tuple_impl< Size, std::tuple<Spaces...> > {
    typedef arithmetic_tuple< typename get_rate_limited_space<Spaces>::type... > type;
  };
  
  template <std::size_t Size, typename... Spaces>
  struct get_rate_limited_space_tuple_impl< Size, arithmetic_tuple<Spaces...> > {
    typedef arithmetic_tuple< typename get_rate_limited_space<Spaces>::type... > type;
  };
  
#else
  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 1, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 2, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 3, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 4, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 5, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 6, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 7, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 8, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 9, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<8,SpaceTuple>::type >::type > type;
  };

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple_impl< 10, SpaceTuple > { 
    typedef arithmetic_tuple< typename get_rate_limited_space< typename arithmetic_tuple_element<0,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<1,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<2,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<3,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<4,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<5,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<6,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<7,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<8,SpaceTuple>::type >::type,
                              typename get_rate_limited_space< typename arithmetic_tuple_element<9,SpaceTuple>::type >::type > type;
  };

#endif

  template <typename SpaceTuple>
  struct get_rate_limited_space_tuple : get_rate_limited_space_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple > { };
  
  
};




template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_illimited_space< reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
  typedef differentiable_space<IndependentSpace,typename detail::get_rate_illimited_space_tuple<SpaceTuple>::type,TupleDistanceMetric> type;
};


template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_limited_space< differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
  typedef reach_time_diff_space<IndependentSpace,typename detail::get_rate_limited_space_tuple<SpaceTuple>::type,TupleDistanceMetric> type;
};



template <typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_illimited_space< metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
  typedef metric_space_tuple< typename detail::get_rate_illimited_space_tuple<SpaceTuple>::type, TupleDistanceMetric> type;
};


template <typename SpaceTuple, typename TupleDistanceMetric>
struct get_rate_limited_space< metric_space_tuple<SpaceTuple,TupleDistanceMetric> > {
  typedef metric_space_tuple< typename detail::get_rate_limited_space_tuple<SpaceTuple>::type, TupleDistanceMetric> type;
};








};

};


#endif








