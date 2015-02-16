/**
 * \file differentiable_space.hpp
 * 
 * This library provides classes that define a differentiable-space. This class template can be used to 
 * glue together a number of spaces by a differentiation relationship,
 * where each differentiation / integration operation (or more formally speaking, each lift and descent 
 * through the tangent bundle) is governed by its own differentiation rule. This class template models 
 * the MetricSpaceConcept (if all underlying spaces do as well), and models the DifferentiableSpaceConcept
 * for as high an order as there are differentiation rules and spaces to support.
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

#ifndef REAK_DIFFERENTIABLE_SPACE_HPP
#define REAK_DIFFERENTIABLE_SPACE_HPP

#include <ReaK/core/base/defs.hpp>

#include "time_topology.hpp"
#include "metric_space_concept.hpp"
#include "tangent_bundle_concept.hpp"
#include "metric_space_tuple.hpp"

#include <cmath>

#include <boost/static_assert.hpp>

#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
#include <type_traits>
#else
#include <boost/type_traits.hpp>
#endif

namespace ReaK {

namespace pp {
  
/**
 * This class defines the default differentiation rule to apply either to lift a 
 * point-difference (e.g. finite-difference) to the tangent space, or to descend 
 * a tangent vector to a point-difference.
 */
struct default_differentiation_rule : public serializable {
  /**
   * This function will lift a point-difference vector into its corresponding tangent vector.
   * This function performs a simple division, dp / dt.
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
    v = dp * (1.0 / dt);
  };
  /**
   * This function will descend a tangent vector into its corresponding point-difference vector.
   * This function performs a simple multiplication, v * dt.
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
    dp = v * dt;
  };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
  virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const { };

  virtual void RK_CALL load(serialization::iarchive& A, unsigned int) { };

  RK_RTTI_MAKE_ABSTRACT_1BASE(default_differentiation_rule,0xC2420000,1,"default_differentiation_rule",serializable)
};

/**
 * This class defines a tuple of default differentiation rules. This is useful for applying the 
 * default differentiation rule to all the differentiation operations on a given differentiable_space.
 */
template <std::size_t Order, typename DiffRule>
struct differentiation_rule_array {
  //BOOST_STATIC_ASSERT(false);
};

template <typename DiffRule>
struct differentiation_rule_array<0,DiffRule> {
  typedef arithmetic_tuple<DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<1,DiffRule> {
  typedef arithmetic_tuple<DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<2,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<3,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<4,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<5,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<6,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<7,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<8,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<9,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

template <typename DiffRule>
struct differentiation_rule_array<10,DiffRule> {
  typedef arithmetic_tuple<DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule,
                           DiffRule> type;
};

  
/**
 * This class template can be used to glue together a number of spaces by a differentiation relationship,
 * where each differentiation / integration operation (or more formally speaking, each lift and descent 
 * through the tangent bundle) is governed by its own differentiation rule. This class template is based 
 * on the metric_space_tuple and will thus model all the concepts that metric_space_tuple models, and it 
 * models the TangentBundleConcept for as high an order as there are differentiation rules and spaces provided.
 * 
 * \tparam IndependentSpace The type of the independent-space against which the differentiation is 
 *                          taken (e.g. time_topology). There are no formal requirements on this type, 
 *                          it is merely used as a placeholder by this class (although the differentiation 
 *                          rules might require more of this type).
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces that are arranged 
 *                    in sequence of differentiation levels (e.g. space 0 -- diff --> space 1 -- diff --> space 2 ...).
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a 
 *                             space-tuple (e.g. arithmetic_tuple).
 * \tparam DiffRuleTuple A tuple type (e.g. arithmetic_tuple) which provides a set of differentiation rules 
 *                       to use in order to move vectors across the tangent bundle 
 */
template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric = manhattan_tuple_distance, 
          typename DiffRule = default_differentiation_rule >
class differentiable_space : public metric_space_tuple<SpaceTuple,TupleDistanceMetric> {
  public:
    typedef typename differentiation_rule_array< arithmetic_tuple_size<SpaceTuple>::type::value - 1, DiffRule >::type diff_rule_tuple;

  protected:
    diff_rule_tuple m_diff_rules;
    
  public:
    typedef differentiable_space< IndependentSpace, SpaceTuple, TupleDistanceMetric, DiffRule > self;
    typedef metric_space_tuple<SpaceTuple,TupleDistanceMetric> base_type;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef typename base_type::distance_metric_type distance_metric_type;
    typedef typename base_type::random_sampler_type random_sampler_type;
        
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
    /**
     * This nested class template is a meta-function to obtain the type of the space of a given 
     * differential order.
     * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
     * \tparam IndependentSpace2 The independent space against which the differentiation is done (e.g. time).
     */
    template <int Idx, typename IndependentSpace2 = IndependentSpace>
    struct space {
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
      BOOST_STATIC_ASSERT((std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value));
#else
      BOOST_STATIC_ASSERT((boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value));
#endif
      typedef typename arithmetic_tuple_element<Idx, SpaceTuple>::type type;
    };
      
    BOOST_STATIC_CONSTANT(std::size_t, differential_order = arithmetic_tuple_size<SpaceTuple>::type::value - 1);
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     * \param aDist The distance metric functor on the space-tuple.
     * \param aDiffRules The differentiation rule tuple to initialize the diff-rule functors with.
     */
    differentiable_space(const SpaceTuple& aSpaces = SpaceTuple(), 
                         const TupleDistanceMetric& aDist = TupleDistanceMetric(),
                         const diff_rule_tuple& aDiffRules = diff_rule_tuple()) :
                         base_type(aSpaces,aDist), 
                         m_diff_rules(aDiffRules) { };
    
    /**
     * This function returns the differentiation functor at a given order against a given independent-space.
     * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
     */
    template <int Idx>
    const DiffRule& get_diff_rule_impl(const IndependentSpace&) const {
      return get<Idx>(m_diff_rules);
    };
    
    /**
     * This function returns the differentiation functor at a given order against a given independent-space.
     * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
     */
    template <int Idx>
    DiffRule& get_diff_rule_impl(const IndependentSpace&) {
      return get<Idx>(m_diff_rules);
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
    template <int Idx>
    typename arithmetic_tuple_element<Idx, point_type>::type 
      lift_to_space(const typename arithmetic_tuple_element<Idx-1, point_difference_type>::type& dp,
                    const typename topology_traits< IndependentSpace >::point_difference_type& dt,
                    const IndependentSpace& t_space) const {
      typename arithmetic_tuple_element<Idx, point_type>::type result;
      get<Idx-1>(m_diff_rules).lift(result, dp, dt, t_space);
      return result;
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
    template <int Idx>
    typename arithmetic_tuple_element<Idx, point_difference_type>::type 
      descend_to_space(const typename arithmetic_tuple_element<Idx+1, point_type>::type& v,
                       const typename topology_traits< IndependentSpace >::point_difference_type& dt,
                       const IndependentSpace& t_space) const {
      typename arithmetic_tuple_element<Idx, point_difference_type>::type result;
      get<Idx>(m_diff_rules).descend(result, v, dt, t_space);
      return result;
    };
    

/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_diff_rules);
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_diff_rules);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400003,1,"differentiable_space",base_type)

};

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule>
struct is_metric_space< differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, DiffRule> > : boost::mpl::true_ { };

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule>
struct is_reversible_space< differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, DiffRule> > : 
  is_reversible_space< metric_space_tuple<SpaceTuple, TupleDistanceMetric> > { };

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule>
struct is_point_distribution< differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, DiffRule> > : boost::mpl::true_ { };

template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule>
struct is_metric_symmetric< differentiable_space<IndependentSpace, SpaceTuple, TupleDistanceMetric, DiffRule> > : 
  is_metric_symmetric< metric_space_tuple<SpaceTuple, TupleDistanceMetric> > { };


template <typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace, typename IndependentSpace2, std::size_t Order>
struct derived_N_order_space< differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>, IndependentSpace2, Order > {
  typedef typename arithmetic_tuple_element<Order, SpaceTuple>::type type;
};


/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& >::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
const typename arithmetic_tuple_element<Idx, SpaceTuple>::type& >::type
#endif
  get_space(const differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& s, const IndependentSpace2&) {
  return get_space<Idx>(s);
};
    
/**
 * This function returns the space at a given differential order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position, 1: velocity, 2: acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& >::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
typename arithmetic_tuple_element<Idx, SpaceTuple>::type& >::type
#endif
  get_space(differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& s, const IndependentSpace2&) {
  return get_space<Idx>(s);
};

/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
const DiffRule& >::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
const DiffRule& >::type
#endif
  get_diff_rule(const differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& s, const IndependentSpace2&) {
  return s.template get_diff_rule_impl<Idx>();
};
    
/**
 * This function returns the differentiation functor at a given order against a given independent-space.
 * \tparam Idx The differential order (e.g. 0: position/velocity, 1: velocity/acceleration).
 */
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
DiffRule& >::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
DiffRule& >::type
#endif
  get_diff_rule(differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& s, const IndependentSpace2&) {
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
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
arithmetic_tuple_element<Idx, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_type> >::type::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
arithmetic_tuple_element<Idx, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_type> >::type::type
#endif
  lift_to_space(const typename arithmetic_tuple_element<Idx-1, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_difference_type>::type& dp,
                const typename topology_traits< IndependentSpace2 >::point_difference_type& dt,
                const differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& space,
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
template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRule, typename IndependentSpace2>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
typename std::enable_if< std::is_convertible< const IndependentSpace2&, const IndependentSpace&>::value,
arithmetic_tuple_element<Idx, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_difference_type> >::type::type
#else
typename boost::enable_if< boost::is_convertible< const IndependentSpace2&, const IndependentSpace&>,
arithmetic_tuple_element<Idx, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_difference_type> >::type::type
#endif
  descend_to_space(const typename arithmetic_tuple_element<Idx+1, typename differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>::point_type>::type& v,
                   const typename topology_traits< IndependentSpace2 >::point_difference_type& dt,
                   const differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRule>& space,
                   const IndependentSpace2& t_space) {
  return space.template descend_to_space<Idx>(v,dt,t_space);
};






};

};



namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRuleTuple>
  struct arithmetic_tuple_size< pp::differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> > : 
    arithmetic_tuple_size< SpaceTuple > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRuleTuple>
  struct arithmetic_tuple_element< Idx, pp::differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> > {
    typedef typename arithmetic_tuple_element< Idx, SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRuleTuple>
  struct arithmetic_tuple_element< Idx, const pp::differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> > {
    typedef typename arithmetic_tuple_element< Idx, const SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRuleTuple>
  struct arithmetic_tuple_element< Idx, volatile pp::differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> > {
    typedef typename arithmetic_tuple_element< Idx, volatile SpaceTuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename IndependentSpace, typename SpaceTuple, typename TupleDistanceMetric, typename DiffRuleTuple>
  struct arithmetic_tuple_element< Idx, const volatile pp::differentiable_space<IndependentSpace,SpaceTuple,TupleDistanceMetric,DiffRuleTuple> > {
    typedef typename arithmetic_tuple_element< Idx, const volatile SpaceTuple >::type type;
  };
  
};

#endif








