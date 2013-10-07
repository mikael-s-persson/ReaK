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
#include "base/serializable.hpp"

#include "path_planning/metric_space_concept.hpp"
#include "path_planning/prob_distribution_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "tuple_distance_metrics.hpp"
#include "default_random_sampler.hpp"
#include "differentiable_space.hpp"

#include "rate_limited_space_metamaps.hpp"

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


namespace detail {
  
  template <typename Idx, typename ForwardIter, typename ReachTimeDiffTuple>
  typename boost::disable_if<
    boost::mpl::less< 
      Idx,
      arithmetic_tuple_size< ReachTimeDiffTuple >
    >,
  void >::type assign_reach_time_diff_rules(ForwardIter, ReachTimeDiffTuple&); // declaration only.
  
  template <typename Idx, typename ForwardIter, typename ReachTimeDiffTuple>
  typename boost::enable_if<
    boost::mpl::less<
      Idx,
      arithmetic_tuple_size< ReachTimeDiffTuple >
    >,
  void >::type assign_reach_time_diff_rules(ForwardIter it, ReachTimeDiffTuple& diff_rule) {
    get< Idx::type::value >(diff_rule) = *it;
    assign_reach_time_diff_rules< typename boost::mpl::next<Idx>::type >(++it, diff_rule);
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
class reach_time_diff_space : public differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, reach_time_differentiation > {
  public:
    typedef reach_time_diff_space< IndependentSpace, SpaceTuple, TupleDistanceMetric > self;
    typedef differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, reach_time_differentiation> base_type;
	typedef typename base_type::diff_rule_tuple diff_rule_tuple;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     * \param aDiffRules The differentiation rule tuple to initialize the diff-rule functors with.
     */
    reach_time_diff_space(const SpaceTuple& aSpaces = SpaceTuple(), 
			  const TupleDistanceMetric& aDist = TupleDistanceMetric(),
			  const diff_rule_tuple& aDiffRules = diff_rule_tuple()) :
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
			  base_type(aSpaces, aDist, detail::construct_reach_time_diff_rules<diff_rule_tuple>(first, last)) { };
    
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
			  base_type(aSpaces, aDist, detail::construct_reach_time_diff_rules<diff_rule_tuple>(aList.begin(), aList.end())) { };
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

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
struct is_metric_space< reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric> > : boost::mpl::true_ { };

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
struct is_point_distribution< reach_time_diff_space<IndependentSpace, SpaceTuple, TupleDistanceMetric> > : boost::mpl::true_ { };


template <typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace, typename IndependentSpace2, std::size_t Order>
struct derived_N_order_space< reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>, IndependentSpace2, Order > {
  typedef typename arithmetic_tuple_element<Order, SpaceTuple>::type type;
};

#if 1
/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(const reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& s, const IndependentSpace2& t) {
  return get_space<Idx>(static_cast<const differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,reach_time_differentiation>&>(s), t);
};
    
/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& get_space(reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& s, const IndependentSpace2& t) {
  return get_space<Idx>(static_cast<differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,reach_time_differentiation>&>(s), t);
};

/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
const reach_time_differentiation& get_diff_rule(const reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& s, const IndependentSpace2&) {
  return s.template get_diff_rule_impl<Idx>();
};
    
/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
reach_time_differentiation& get_diff_rule(reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& s, const IndependentSpace2&) {
  return s.template get_diff_rule_impl<Idx>();
};
    
/**
 * This function lifts a point-difference in space Idx-1 into a point in space Idx.
 * \tparam Idx The differential order of the destination space.
 * \param dp The point-difference in the space Idx-1.
 * \param dt The point-difference in the independent-space (e.g. time).
 * \param space The differentiable-space.
 * \param t_space The independent-space.
 * \return The point in space Idx which is the tangential lift of the point-difference in space Idx-1.
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
typename arithmetic_tuple_element<Idx, typename reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>::point_type>::type 
  lift_to_space(const typename arithmetic_tuple_element<Idx-1, typename reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>::point_difference_type>::type& dp,
                const typename topology_traits< IndependentSpace2 >::point_difference_type& dt,
                const reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& space,
                const IndependentSpace2& t_space) {
  return space.template lift_to_space<Idx>(dp,dt,t_space);
};
    
    /**
     * This function descends a point in space Idx+1 into a point-difference in space Idx.
     * \tparam Idx The differential order of the destination space.
     * \param v The point in the space Idx+1.
     * \param dt The point-difference in the independent-space (e.g. time).
     * \param space The differentiable-space.
     * \param t_space The independent-space.
     * \return The point-difference in space Idx which is the tangential descent of the point in space Idx+1.
     */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename IndependentSpace2>
typename arithmetic_tuple_element<Idx, typename reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>::point_difference_type>::type 
  descend_to_space(const typename arithmetic_tuple_element<Idx+1, typename reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>::point_type>::type& v,
                   const typename topology_traits< IndependentSpace2 >::point_difference_type& dt,
                   const reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric>& space,
                   const IndependentSpace2& t_space) {
  return space.template descend_to_space<Idx>(v,dt,t_space);
};

#endif


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



namespace ReaK {
  

/* Specialization, see general template docs. */
  template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_size< pp::reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > : 
    arithmetic_tuple_size< SpaceTuple > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_element< Idx, pp::reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
    typedef typename arithmetic_tuple_element< Idx, SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_element< Idx, const pp::reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
    typedef typename arithmetic_tuple_element< Idx, const SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_element< Idx, volatile pp::reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
    typedef typename arithmetic_tuple_element< Idx, volatile SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric>
  struct arithmetic_tuple_element< Idx, const volatile pp::reach_time_diff_space<IndependentSpace,SpaceTuple,TupleDistanceMetric> > {
    typedef typename arithmetic_tuple_element< Idx, const volatile SpaceTuple >::type type;
  };
  
};


#endif








