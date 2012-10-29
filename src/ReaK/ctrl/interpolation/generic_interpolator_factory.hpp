/**
 * \file generic_interpolator_factory.hpp
 * 
 * This library provides an implementation of a generic interpolator factory.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_GENERIC_INTERPOLATOR_FACTORY_HPP
#define REAK_GENERIC_INTERPOLATOR_FACTORY_HPP

#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/temporal_space_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "topologies/metric_space_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

namespace ReaK {

namespace pp {
  
namespace detail {

  

template < template <typename, typename> class InterpolatorImpl,
           typename SpaceType, typename TimeSpaceType>
class generic_interpolator_impl : public InterpolatorImpl< SpaceType, TimeSpaceType > { };


  
  template <std::size_t Size, typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace >
  struct generic_interpolator_impl_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  
  template <std::size_t Size, template <typename,typename> class InterpolatorImpl, typename TimeSpace, typename... Spaces>
  struct generic_interpolator_impl_tuple_impl< Size, std::tuple<Spaces...>, InterpolatorImpl, TimeSpace > {
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, Spaces, TimeSpace>... > type;
  };
  
  template <std::size_t Size, template <typename,typename> class InterpolatorImpl, typename TimeSpace, typename... Spaces>
  struct generic_interpolator_impl_tuple_impl< Size, arithmetic_tuple<Spaces...>, InterpolatorImpl, TimeSpace > {
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, Spaces, TimeSpace>... > type;
  };
  
#else
  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 1, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 2, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 3, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 4, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 5, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 6, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<5,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 7, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<5,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<6,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 8, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<5,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<6,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<7,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 9, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<5,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<6,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<7,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<8,SpaceTuple>::type, TimeSpace > > type;
  };

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple_impl< 10, SpaceTuple, InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<0,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<1,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<2,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<3,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<4,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<5,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<6,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<7,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<8,SpaceTuple>::type, TimeSpace >,
                              generic_interpolator_impl< InterpolatorImpl, typename arithmetic_tuple_element<9,SpaceTuple>::type, TimeSpace > > type;
  };

#endif

  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple : generic_interpolator_impl_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple, InterpolatorImpl, TimeSpace > { };
  


template <std::size_t Idx>
class gen_interpolator_recursion_impl {
  public:
    
    template <typename InterpolatorTuple, typename PointType, typename SpaceTuple, typename TimeSpace, typename Factory>
    static void initialize(InterpolatorTuple& interp,
                                  const PointType& start_point, const PointType& end_point, double dt,
		                  const SpaceTuple& space, const TimeSpace& t_space, const Factory& factory) {
      gen_interpolator_recursion_impl< Idx-1 >::initialize(interp,start_point,end_point,dt,space,t_space,factory);
      
      get<Idx>(interp).initialize(get<Idx>(start_point), get<Idx>(end_point), dt, 
				  get<Idx>(space), t_space, factory);
    };
    
    template <typename InterpolatorTuple, typename PointType, typename SpaceTuple, typename TimeSpace, typename Factory>
    static void compute_point(const InterpolatorTuple& interp,
                              PointType& result, const PointType& start_point, const PointType& end_point, 
		              const SpaceTuple& space, const TimeSpace& t_space, 
		              double dt, double dt_total, const Factory& factory) {
      gen_interpolator_recursion_impl< Idx-1 >::compute_point(interp,result,start_point,end_point,space,t_space,dt,dt_total,factory);
      
      get<Idx>(interp).compute_point(get<Idx>(result), get<Idx>(start_point), get<Idx>(end_point),
				     get<Idx>(space), t_space, dt, dt_total, factory);
    };
    
    template <typename InterpolatorTuple>
    static double get_minimum_travel_time(const InterpolatorTuple& interp) {
      double result = gen_interpolator_recursion_impl< Idx-1 >::get_minimum_travel_time(interp);
      
      double result_0 = get<Idx>(interp).get_minimum_travel_time();
      if( result > result_0 ) 
	return result;
      else
	return result_0;
    };
    
    
};


template <>
class gen_interpolator_recursion_impl<0> {
  public:
    
    template <typename InterpolatorTuple, typename PointType, typename SpaceTuple, typename TimeSpace, typename Factory>
    static void initialize(InterpolatorTuple& interp,
                                  const PointType& start_point, const PointType& end_point, double dt,
		                  const SpaceTuple& space, const TimeSpace& t_space, const Factory& factory) {
      get<0>(interp).initialize(get<0>(start_point), get<0>(end_point), dt, 
				  get<0>(space), t_space, factory);
    };
    
    template <typename InterpolatorTuple, typename PointType, typename SpaceTuple, typename TimeSpace, typename Factory>
    static void compute_point(const InterpolatorTuple& interp,
                              PointType& result, const PointType& start_point, const PointType& end_point, 
		              const SpaceTuple& space, const TimeSpace& t_space, 
		              double dt, double dt_total, const Factory& factory) {
      get<0>(interp).compute_point(get<0>(result), get<0>(start_point), get<0>(end_point),
				   get<0>(space), t_space, dt, dt_total, factory);
    };
    
    template <typename InterpolatorTuple>
    static double get_minimum_travel_time(const InterpolatorTuple& interp) {
      return get<0>(interp).get_minimum_travel_time();
    };
    
    
};



template < template <typename, typename> class InterpolatorImpl,
           typename TimeSpaceType, typename SpaceTuple, typename TupleDistMetric>
class generic_interpolator_impl< InterpolatorImpl, metric_space_tuple<SpaceTuple,TupleDistMetric>, TimeSpaceType > {
  public:
    typedef generic_interpolator_impl< InterpolatorImpl, metric_space_tuple<SpaceTuple,TupleDistMetric>, TimeSpaceType > self;
    
    typedef metric_space_tuple<SpaceTuple,TupleDistMetric> SpaceType;
    typedef typename topology_traits<SpaceType>::point_type point_type;
  private:
    typename generic_interpolator_impl_tuple< SpaceTuple, InterpolatorImpl, TimeSpaceType >::type interp;
    
  public:
    
    generic_interpolator_impl() { };
    
    template <typename Factory>
    void initialize(const point_type& start_point, const point_type& end_point, double dt,
		    const SpaceType& space, const TimeSpaceType& t_space, const Factory& factory) {
      gen_interpolator_recursion_impl< arithmetic_tuple_size< SpaceTuple >::value - 1 >::initialize(interp, start_point, end_point, dt, space, t_space, factory);
    };
    
    template <typename Factory>
    void compute_point(point_type& result, const point_type& start_point, const point_type& end_point,
		       const SpaceType& space, const TimeSpaceType& t_space, 
		       double dt, double dt_total, const Factory& factory) const {
      gen_interpolator_recursion_impl< arithmetic_tuple_size< SpaceTuple >::value - 1 >::compute_point(interp, result, start_point, end_point, space, t_space, dt, dt_total, factory);
    };
    
    double get_minimum_travel_time() const {
      return gen_interpolator_recursion_impl< arithmetic_tuple_size< SpaceTuple >::value - 1 >::get_minimum_travel_time(interp);
    };
    
    
};


};




/**
 * This functor class implements a generic interpolation in a temporal topology.
 * \tparam Factory The interpolator factory type which created this generic interpolator.
 * \tparam InterpolatorImpl The interpolator implementation template which is used to store and perform the specifics of the actual interpolation.
 */
template <typename Factory, 
          template<typename, typename> class InterpolatorImpl >
class generic_interpolator {
  public:
    typedef generic_interpolator<Factory, InterpolatorImpl> self;
    typedef typename Factory::point_type point_type;
    typedef typename Factory::topology topology;
    typedef typename temporal_space_traits<topology>::space_topology SpaceType;
    typedef typename temporal_space_traits<topology>::time_topology TimeSpaceType;
    
    typedef detail::generic_interpolator_impl< InterpolatorImpl, SpaceType, TimeSpaceType> interpolator_impl_type;
    
  private:
    const Factory* parent;
    const point_type* start_point;
    const point_type* end_point;
    interpolator_impl_type interp;
    
    void update_delta_value() {
      if(parent && start_point && end_point) {
	double delta_time = end_point->time - start_point->time;
	
	interp.initialize(start_point->pt, end_point->pt, delta_time,
			  parent->get_temporal_space()->get_space_topology(),
			  parent->get_temporal_space()->get_time_topology(),
			  *parent);
	
      };
    };
  
  public:
    
    
    /**
     * Default constructor.
     */
    generic_interpolator(const Factory* aParent = NULL, const point_type* aStart = NULL, const point_type* aEnd = NULL) :
                         parent(aParent), start_point(aStart), end_point(aEnd) {
      update_delta_value();
    };
    
    void set_segment(const point_type* aStart, const point_type* aEnd) {
      start_point = aStart;
      end_point = aEnd;
      update_delta_value();
    };
    
    const point_type* get_start_point() const { return start_point; };
    const point_type* get_end_point() const { return end_point; };
    
    template <typename DistanceMetric>
    double travel_distance_to(const point_type& pt, const DistanceMetric& dist) const {
      BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,topology>));
      if(parent && start_point)
	return dist(pt, *start_point, *(parent->get_temporal_space()));
      else
	return 0.0;
    };
    
    template <typename DistanceMetric>
    double travel_distance_from(const point_type& pt, const DistanceMetric& dist) const {
      BOOST_CONCEPT_ASSERT((DistanceMetricConcept<DistanceMetric,topology>));
      if(parent && end_point)
	return dist(*end_point, pt, *(parent->get_temporal_space()));
      else
	return 0.0;
    };
    
    point_type get_point_at_time(double t) const {
      if(!parent || !start_point || !end_point)
	return point_type();
      
      double dt_total = end_point->time - start_point->time;
      if(interp.get_minimum_travel_time() > dt_total)
	dt_total = interp.get_minimum_travel_time();
      double dt = t - start_point->time;
      
      point_type result;
      result.time = t;
      
      interp.compute_point(result.pt, start_point->pt, end_point->pt, parent->get_temporal_space()->get_space_topology(), parent->get_temporal_space()->get_time_topology(), dt, dt_total, *parent);
      
      return result;   
    };
    
    double get_minimum_travel_time() const {
      if(parent && start_point && end_point)
	return interp.get_minimum_travel_time();
      else 
	return std::numeric_limits<double>::infinity();
    };
    
    bool is_segment_feasible() const {
      if(parent && start_point && end_point)
        return (interp.get_minimum_travel_time() < end_point->time - start_point->time);
      else
	return false;
    };
    
};



/**
 * This meta-function is used to obtain the generic interpolator for a spatial interpolation (space and "time" 
 * topologies are separate), given an interpolator tag type.
 * \note The general template here is bogus. This class template needs to be specialized for each interpolator tag.
 */
template <typename InterpTag, typename SpaceType, typename TimeTopology = time_topology>
struct get_tagged_spatial_interpolator {
  typedef void type; // dummy, generalization is not possible.
  typedef void* pseudo_factory_type;
};

/**
 * This meta-function is used to obtain the generic interpolator for a temporal interpolation (space and "time" 
 * topologies are joined in one "temporal" topology), given an interpolator tag type.
 * \note The general template here is bogus. This class template needs to be specialized for each interpolator tag.
 */
template <typename InterpTag, typename TemporalSpaceType>
struct get_tagged_temporal_interpolator {
  typedef void type; // dummy, generalization is not possible.
};




};

};

#endif









