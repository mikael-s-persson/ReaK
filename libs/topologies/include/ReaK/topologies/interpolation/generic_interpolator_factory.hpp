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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

#include <ReaK/topologies/spaces/tangent_bundle_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>
#include <ReaK/topologies/spaces/metric_space_tuple.hpp>
#include <ReaK/topologies/spaces/time_topology.hpp>

#include <boost/concept_check.hpp>
#include <cmath>

namespace ReaK {

namespace pp {
  
namespace detail { namespace {

  

template < template <typename, typename> class InterpolatorImpl,
           typename SpaceType, typename TimeSpaceType>
class generic_interpolator_impl : public InterpolatorImpl< SpaceType, TimeSpaceType > { };


  
  template <typename SpaceTuple, template <typename,typename> class InterpolatorImpl, typename TimeSpace >
  struct generic_interpolator_impl_tuple { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
  
  template <template <typename,typename> class InterpolatorImpl, typename TimeSpace, typename... Spaces>
  struct generic_interpolator_impl_tuple< arithmetic_tuple<Spaces...>, InterpolatorImpl, TimeSpace > {
    typedef arithmetic_tuple< generic_interpolator_impl< InterpolatorImpl, Spaces, TimeSpace>... > type;
  };
  
#else
  
  template <typename S0, template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,void,void,void,void,void,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace > 
    > type;
  };
  
  template <typename S0, typename S1, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,void,void,void,void,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,void,void,void,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,void,void,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,void,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            typename S5, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,S5,void,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S5, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            typename S5, typename S6, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,S5,S6,void,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S5, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S6, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            typename S5, typename S6, typename S7, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,S5,S6,S7,void,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S5, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S6, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S7, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            typename S5, typename S6, typename S7, typename S8, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,S5,S6,S7,S8,void>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S5, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S6, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S7, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S8, TimeSpace >
    > type;
  };
  
  template <typename S0, typename S1, typename S2, typename S3, typename S4, 
            typename S5, typename S6, typename S7, typename S8, typename S9, 
            template <typename,typename> class InterpolatorImpl, typename TimeSpace>
  struct generic_interpolator_impl_tuple< 
      arithmetic_tuple<S0,S1,S2,S3,S4,S5,S6,S7,S8,S9>, 
    InterpolatorImpl, TimeSpace > { 
    typedef arithmetic_tuple< 
      generic_interpolator_impl< InterpolatorImpl, S0, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S1, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S2, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S3, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S4, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S5, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S6, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S7, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S8, TimeSpace >,
      generic_interpolator_impl< InterpolatorImpl, S9, TimeSpace >
    > type;
  };

#endif



template <typename TimeSpace, typename Factory>
struct gen_interp_impl_initializer {
  const TimeSpace* p_t_space;
  const Factory* p_factory;
  double dt;
  gen_interp_impl_initializer(const TimeSpace& t_space, const Factory& factory, double aDt) :
                              p_t_space(&t_space), p_factory(&factory), dt(aDt) {};
  template <typename Interpolator, typename Point, typename Space>
  void operator()(Interpolator& interp, const Point& start_point, 
                  const Point& end_point, const Space& space) const {
    interp.initialize(start_point, end_point, dt, space, *p_t_space, *p_factory);
  };
};

template <typename TimeSpace, typename Factory>
struct gen_interp_impl_computer {
  const TimeSpace* p_t_space;
  const Factory* p_factory;
  double dt, dt_total;
  gen_interp_impl_computer(const TimeSpace& t_space, const Factory& factory, 
                           double aDt, double aDtTotal) :
                           p_t_space(&t_space), p_factory(&factory), 
                           dt(aDt), dt_total(aDtTotal) {};
  template <typename Interpolator, typename Point, typename Space>
  void operator()(const Interpolator& interp, Point& result, 
                  const Point& start_point, const Point& end_point, const Space& space) const {
    interp.compute_point(result, start_point, end_point, space, *p_t_space, dt, dt_total, *p_factory);
  };
};

struct gen_interp_impl_mintime {
  double* p_result;
  gen_interp_impl_mintime(double& result) : p_result(&result) {};
  template <typename Interpolator>
  void operator()(const Interpolator& interp) const {
    double r0 = interp.get_minimum_travel_time();
    if( (*p_result) < r0 )
      (*p_result) = r0;
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
      tuple_for_each(interp, start_point, end_point, space, 
        gen_interp_impl_initializer<TimeSpaceType,Factory>(t_space, factory, dt)
      );
    };
    
    template <typename Factory>
    void compute_point(point_type& result, const point_type& start_point, const point_type& end_point,
                       const SpaceType& space, const TimeSpaceType& t_space, 
                       double dt, double dt_total, const Factory& factory) const {
      tuple_for_each(interp, result, start_point, end_point, space, 
        gen_interp_impl_computer<TimeSpaceType,Factory>(t_space, factory, dt, dt_total)
      );
    };
    
    double get_minimum_travel_time() const {
      double result = 0.0;
      tuple_for_each(interp, gen_interp_impl_mintime(result));
      return result;
    };
    
};


}; };



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









