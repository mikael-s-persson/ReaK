/**
 * \file cubic_hermite_interp.hpp
 * 
 * This library provides an implementation of a trajectory within a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a cubic Hermite interpolation (cubic hermite spline, or cspline).
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

#ifndef REAK_CUBIC_HERMITE_INTERP_HPP
#define REAK_CUBIC_HERMITE_INTERP_HPP

#include "path_planning/spatial_trajectory_concept.hpp"

#include "path_planning/differentiable_space_concept.hpp"

#include "interpolated_trajectory.hpp"

#include "lin_alg/arithmetic_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

#include <list>
#include <map>
#include <limits>
#include "lin_alg/mat_num_exceptions.hpp"

namespace ReaK {

namespace pp {
  
  
namespace detail {
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type cubic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                              const DiffSpace& space, const TimeSpace& t_space,
					      double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<0>(result) = space.get_space<0>(t_space).adjust(get<0>(a), 
      ((3.0 - 2.0 * t_normal) * t_normal * t_normal) * space.get_space<0>(t_space).difference(get<0>(b),get<0>(a))
      + (t_normal * (1.0 - t_normal * (2.0 + t_normal))) * space.descend_to_space<0>(get<1>(a),t_factor, t_space)
      + (t_normal * t_normal * (t_normal - 1.0)) * space.descend_to_space<0>(get<1>(b),t_factor, t_space) );
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type cubic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                              const DiffSpace& space, const TimeSpace& t_space,
				 	      double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    cubic_hermite_interpolate_impl< boost::mpl::size_t<0>, PointType, DiffSpace, TimeSpace >(result,a,b,space,t_space,t_factor,t_normal);
    
    get<1>(result) = space.get_space<1>(t_space).adjust(get<1>(a),
      ((1.0 - t_normal) * 6.0 * t_normal) * space.get_space<1>(t_space).difference( space.lift_to_space<1>(space.get_space<0>(t_space).difference(get<0>(b),get<0>(a)), t_factor, t_space), get<1>(a))
      + (t_normal * (2.0 - 3.0 * t_normal)) * space.get_space<1>(t_space).difference(get<1>(b),get<1>(a)));
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<2> 
    >,
  void >::type cubic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                              const DiffSpace& space, const TimeSpace& t_space,
					      double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    cubic_hermite_interpolate_impl< boost::mpl::size_t<1>, PointType, DiffSpace, TimeSpace >(result,a,b,space,t_space,t_factor,t_normal);
    
    get<2>(result) = space.get_space<2>(t_space).adjust( space.lift_to_space<2>(space.get_space<1>(t_space).difference(get<1>(a),get<1>(b)), t_factor, t_space),
      (6.0 - 12.0 * t_normal) * space.get_space<2>(t_space).difference( space.lift_to_space<2>( space.get_space<1>(t_space).difference( space.lift_to_space<1>(space.get_space<0>(t_space).difference(get<0>(b),get<0>(a)), t_factor, t_space), get<1>(a)), t_factor, t_space), space.lift_to_space<2>(0.5 * space.get_space<1>(t_space).difference(get<1>(a),get<1>(b)), t_factor, t_space)));
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type cubic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                              const DiffSpace& space, const TimeSpace& t_space,
					      double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    cubic_hermite_interpolate_impl< boost::mpl::size_t<2>, PointType, DiffSpace, TimeSpace >(result,a,b,space,t_space,t_factor,t_normal);
    
    get<3>(result) = space.lift_to_space<3>(-12.0 * space.get_space<2>(t_space).difference( space.lift_to_space<2>( space.get_space<1>(t_space).difference( space.lift_to_space<1>(space.get_space<0>(t_space).difference(get<0>(b),get<0>(a)), t_factor, t_space), get<1>(a)), t_factor, t_space), space.lift_to_space<2>(0.5 * space.get_space<1>(t_space).difference(get<1>(a),get<1>(b)), t_factor, t_space)),t_factor, t_space);
  };
  
  template <typename Idx, typename PointType, typename DiffSpace, typename TimeSpace>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<3> 
    >,
  void >::type cubic_hermite_interpolate_impl(PointType& result, const PointType& a, const PointType& b,
                                              const DiffSpace& space, const TimeSpace& t_space,
					      double t_factor, double t_normal) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    cubic_hermite_interpolate_impl< boost::mpl::prior<Idx>, PointType, DiffSpace, TimeSpace >(result,a,b,space,t_space,t_factor,t_normal);
    
    get< Idx::type::value >(result) = space.get_space< Idx::type::value >(t_space).origin();
  };
  
};



template <typename PointType, typename Topology>
PointType cubic_hermite_interpolate(const PointType& a, const PointType& b, double t, const Topology& space) {
  BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
  BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 1, typename temporal_topology_traits<Topology>::time_topology >));
  double t_factor = b.time - a.time;
  if(std::fabs(t_factor) < std::numeric_limits<double>::epsilon())
    throw singularity_error("Normalizing factor in cubic Hermite spline is zero!");
  double t_normal = (t - a.time) / (b.time - a.time);
      
  PointType result;
  result.time = t;
  detail::cubic_hermite_interpolate_impl<boost::mpl::size_t<differentiable_space_traits< SpaceType >::order> >(result.pt, a.pt, b.pt, space.get_space_topology(), space.get_time_topology(), t_factor, t_normal);
      
  return result;      
};


struct cubic_hermite_interpolator {
  template <typename PointType, typename Topology>
  PointType operator()(const PointType& a, const PointType& b, double t, const Topology& space) const {
    return cubic_hermite_interpolate(a,b,t,space);
  };
};




  
/**
 * This class implements a trajectory in a temporal and once-differentiable topology.
 * The trajectory is represented by a set of waypoints and all intermediate points 
 * are computed with a linear interpolation. This class models the SpatialTrajectoryConcept.
 * \tparam Topology The topology type on which the points and the path can reside, should model the TemporalSpaceConcept and the DifferentiableSpaceConcept (order 1 with space against time).
 * \tparam DistanceMetric The distance metric used to assess the distance between points in the path, should model the DistanceMetricConcept.
 */
template <typename Topology, typename DistanceMetric = default_distance_metric>
class cubic_hermite_interp : public interpolated_trajectory<Topology,cubic_hermite_interpolator,DistanceMetric> {
  public:
    
    BOOST_CONCEPT_ASSERT((TemporalSpaceConcept<Topology>));
    BOOST_CONCEPT_ASSERT((DifferentiableSpaceConcept< typename temporal_topology_traits<Topology>::space_topology, 1, typename temporal_topology_traits<Topology>::time_topology >));
    
    typedef cubic_hermite_interp<Topology,DistanceMetric> self;
    typedef interpolated_trajectory<Topology,cubic_hermite_interpolator,DistanceMetric> base_class_type;
    
  public:
    /**
     * Constructs the path from a space, assumes the start and end are at the origin 
     * of the space.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    explicit cubic_hermite_interp(const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                                  base_class_type(aSpace, aDist) { };
    
    /**
     * Constructs the path from a space, the start and end points.
     * \param aSpace The space on which the path is.
     * \param aStart The start point of the path.
     * \param aEnd The end-point of the path.
     * \param aDist The distance metric functor that the path should use.
     */
    cubic_hermite_interp(const topology& aSpace, const point_type& aStart, const point_type& aEnd, const distance_metric& aDist = distance_metric()) :
                         base_class_type(aSpace, aStart, aEnd, aDist) { };
			
    /**
     * Constructs the path from a range of points and their space.
     * \tparam ForwardIter A forward-iterator type for getting points to initialize the path with.
     * \param aBegin An iterator to the first point of the path.
     * \param aEnd An iterator to the second point of the path.
     * \param aSpace The space on which the path is.
     * \param aDist The distance metric functor that the path should use.
     */
    template <typename ForwardIter>
    cubic_hermite_interp(ForwardIter aBegin, ForwardIter aEnd, const topology& aSpace, const distance_metric& aDist = distance_metric()) : 
                         base_class_type(aBegin, aEnd, aSpace, aDist) { };
    
};



};

};

#endif









