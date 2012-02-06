/**
 * \file rate_limited_spaces.hpp
 * 
 * This library provides classes that define a number of reach-time-space class templates. A reach-time-space is 
 * a transformation on a differentiable space (tangent bundle) such that points in the spaces and their differences
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
 * This class template can be used to glue together a number of spaces into a tuple. This class template models 
 * the MetricSpaceConcept (if all underlying spaces do as well), and it is also a tuple class (meaning 
 * that the Boost.Tuple or std::tuple meta-functions, as well as the ReaK.Arithmetic-tuple meta-functions 
 * will work on this class, with the usual semantics).
 * 
 * \tparam SpaceTuple A tuple type (e.g. arithmetic_tuple) which provides a set of spaces to glue together.
 * \tparam TupleDistanceMetric A distance metric type which models the DistanceMetricConcept and operates on a space-tuple (e.g. arithmetic_tuple).
 */
template <typename DiffSpace, typename TimeSpace = time_topology >
class order1_rate_limited_space : public serialization::serializable {
  protected:
    DiffSpace m_space;
    TimeSpace t_space;
    
  public:
    typedef order1_rate_limited_space< DiffSpace, TimeSpace > self;
    
    typedef TimeSpace time_topology;
    typedef DiffSpace space_topology;
    
    typedef typename topology_traits<DiffSpace>::point_type point_type;
    typedef typename topology_traits<DiffSpace>::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_CONCEPT_ASSERT((TangentBundleConcept< space_topology, 1, time_topology >));
    BOOST_CONCEPT_ASSERT((BoundedSpaceConcept<typename derived_N_order_space< space_topology, time_topology, 0 >::type>));
    BOOST_CONCEPT_ASSERT((SphereBoundedSpaceConcept<typename derived_N_order_space< space_topology, time_topology, 1 >::type>));
    
    /**
     * Parametrized and default constructor.
     * \param aSpaces The space tuple to initialize the spaces with.
     */
    order1_rate_limited_space(const DiffSpace& aSpace = DiffSpace()) :
		              m_space(aSpace) { };
      
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      return get(distance_metric, get_space<0>(m_space,t_space))(get<0>(p1), get<0>(p2), get_space<0>(m_space,t_space))
           + get(distance_metric, get_space<1>(m_space,t_space))(get<1>(p1), get<1>(p2), get_space<1>(m_space,t_space));
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& dp) const {
      return get(distance_metric, get_space<0>(m_space,t_space))(get<0>(dp), get_space<0>(m_space,t_space))
           + get(distance_metric, get_space<1>(m_space,t_space))(get<1>(dp), get_space<1>(m_space,t_space));
    };
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      typedef typename derived_N_order_space<space_topology,time_topology,0>::type Space0;
      typedef typename derived_N_order_space<space_topology,time_topology,1>::type Space1;
    
      typedef typename topology_traits<Space0>::point_type PointType0;
      typedef typename topology_traits<Space1>::point_type PointType1;
    
      typedef typename topology_traits<Space0>::point_difference_type PointDiff0;
      typedef typename topology_traits<Space1>::point_difference_type PointDiff1;
      
      const Space0& s0 = get_space<0>(m_space,m_space.get_time_topology());
      const Space1& s1 = get_space<1>(m_space,m_space.get_time_topology());
      
      while(true) {
        //generate random point in underlying topology:
        point_type result = get(random_sampler, m_space)(m_space);
        
        //check if the position is possible based on the time it takes to stop the motion.
        double tv = get(distance_metric, s1)(get<1>(result.pt), s1.origin(), s1);
	PointDiff0 dp0 = descend_to_space<0>(get<1>(result.pt), 0.5 * tv, m_space.get_space_topology(), m_space.get_time_topology());
        if( s0.is_in_bounds(s0.adjust(get<0>(result.pt),  dp0)) && 
	    s0.is_in_bounds(s0.adjust(get<0>(result.pt), -dp0)) )
	  return result;
      };
    };
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1,p2);
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
    
    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC240000A,1,"order1_rate_limited_space",serialization::serializable)

};




};

};


#endif








