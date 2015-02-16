/**
 * \file interpolated_topologies.hpp
 * 
 * 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_INTERPOLATED_TOPOLOGIES_HPP
#define REAK_INTERPOLATED_TOPOLOGIES_HPP

#include <ReaK/core/base/defs.hpp>

#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/proper_metric_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>
#include <ReaK/topologies/spaces/reversible_space_concept.hpp>
#include <ReaK/topologies/spaces/default_random_sampler.hpp>
#include <ReaK/topologies/spaces/rate_limited_space_metamaps.hpp>
#include <ReaK/topologies/interpolation/generic_interpolator_factory.hpp>

#include <boost/concept_check.hpp>

#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
#include <boost/function.hpp>
#else
#include <functional>
#endif

#include <boost/mpl/if.hpp>


namespace ReaK {

namespace pp {


namespace detail {
  
  template <typename BaseTopology, bool IsTemporal>
  struct interp_topo_base_get_temporal_topos {
    typedef typename temporal_space_traits< BaseTopology >::space_topology space_topology;
    typedef typename temporal_space_traits< BaseTopology >::time_topology time_topology;
  };
  
  template <typename BaseTopology>
  struct interp_topo_base_get_temporal_topos<BaseTopology, false> {
    typedef int space_topology;
    typedef int time_topology;
  };
  
  template <typename BaseTopology>
  struct interp_topo_base_temporal_topos : 
    interp_topo_base_get_temporal_topos<BaseTopology, is_temporal_space< BaseTopology >::type::value > { };
  
  
};


  

/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class interpolated_topology_base : public BaseTopology {
  public:
    BOOST_CONCEPT_ASSERT((TopologyConcept<BaseTopology>));
    
    typedef interpolated_topology_base<BaseTopology> self;
    
    typedef typename topology_traits< BaseTopology >::point_type point_type;
    typedef typename topology_traits< BaseTopology >::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef detail::interp_topo_base_temporal_topos<BaseTopology> TemporalToposSelect;
    typedef typename TemporalToposSelect::space_topology space_topology;
    typedef typename TemporalToposSelect::time_topology time_topology;
    
    typedef BaseTopology super_space_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< BaseTopology >::dimensions);
    
    
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
    typedef boost::function< bool(const point_type&) > validity_predicate_type;
#else
    typedef std::function< bool(const point_type&) > validity_predicate_type;
#endif
    
  protected:
    
    virtual point_type interp_topo_move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return static_cast<const BaseTopology&>(*this).move_position_toward(a, fraction, b);
    };
    
    virtual point_type interp_topo_move_position_back_to(const point_type& a, double fraction, const point_type& b) const {
      return static_cast<const BaseTopology&>(*this).move_position_back_to(a, fraction, b);
    };
    
    virtual double interp_topo_get_distance(const point_type& a, const point_type& b) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      return get(distance_metric, b_space)(a, b, b_space);
    };
    
    virtual double interp_topo_get_proper_distance(const point_type& a, const point_type& b) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      return get(proper_metric, b_space)(a, b, b_space);
    };
    
    virtual point_type interp_topo_move_position_toward_pred(const point_type& a, double fraction, const point_type& b,
                                                             double min_dist_interval, validity_predicate_type predicate) const {
      
      double dist_tot = this->interp_topo_get_distance(a, b);
      if(dist_tot == std::numeric_limits<double>::infinity())
        return a;
      if(dist_tot < min_dist_interval)
        return this->interp_topo_move_position_toward(a, fraction, b);
      
      double dist_inter = dist_tot * fraction;
      double dist_cur = min_dist_interval;
      point_type result = a;
      point_type last_result = a;
      while(dist_cur < dist_inter) {
        result = this->interp_topo_move_position_toward(a, dist_cur / dist_tot, b);
        if( !predicate(result) )
          return last_result;
        dist_cur += min_dist_interval;
        last_result = result;
      };
      if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
        return b;
      else if(fraction == 0.0)
        return a;
      
      return this->interp_topo_move_position_toward(a, fraction, b);
    };
    
    virtual point_type interp_topo_move_position_back_to_pred(const point_type& a, double fraction, const point_type& b,
                                                              double min_dist_interval, validity_predicate_type predicate) const {
      
      double dist_tot = this->interp_topo_get_distance(a, b);
      if(dist_tot == std::numeric_limits<double>::infinity())
        return b;
      if(dist_tot < min_dist_interval)
        return this->interp_topo_move_position_back_to(a, fraction, b);
      
      double dist_inter = dist_tot * fraction;
      double dist_cur = min_dist_interval;
      point_type result = b;
      point_type last_result = b;
      while(dist_cur < dist_inter) {
        result = this->interp_topo_move_position_back_to(a, dist_cur / dist_tot, b);
        if( !predicate(result) )
          return last_result;
        dist_cur += min_dist_interval;
        last_result = result;
      };
      if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
        return b;
      else if(fraction == 0.0)
        return a;
      
      return this->interp_topo_move_position_back_to(a, fraction, b);
    };
    
    virtual double interp_topo_get_distance_pred(const point_type& a, const point_type& b, double min_dist_interval, validity_predicate_type predicate) const {
      point_type b_tmp = this->interp_topo_move_position_toward_pred(a, 1.0, b, min_dist_interval, predicate);
      if(this->interp_topo_get_distance(b_tmp, b) < std::numeric_limits< double >::epsilon())
        return this->interp_topo_get_distance(a, b); //if b is reachable from a.
      else
        return std::numeric_limits<double>::infinity(); //b is not reachable from a.
    };
    
    virtual double interp_topo_get_norm(const point_difference_type& dp) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      return get(distance_metric, b_space)(dp, b_space);
    };
    
    virtual double interp_topo_get_proper_norm(const point_difference_type& dp) const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      return get(proper_metric, b_space)(dp, b_space);
    };
    
    virtual bool interp_topo_is_in_bounds(const point_type& a) const {
      return static_cast<const BaseTopology&>(*this).is_in_bounds(a);
    };
    
    virtual point_type interp_topo_random_point() const {
      const BaseTopology& b_space = static_cast<const BaseTopology&>(*this);
      return get(random_sampler, b_space)(b_space);
    };
    
    virtual point_type interp_topo_random_point_pred(validity_predicate_type predicate) const {
      point_type result;
      while(!predicate(result = this->interp_topo_random_point())) ; //output only free C-space points.
      return result;
    };
    
  private:
    
    // non-copyable
    interpolated_topology_base(const self& aTopo);
    self& operator=(const self& aTopo);
    
  public:
    
    explicit interpolated_topology_base(const BaseTopology& aTopo) : 
                                        BaseTopology(aTopo) { };
    
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
    template <typename... Args>
    interpolated_topology_base(Args&&... args) : BaseTopology(std::forward<Args>(args)...) { };
#else
    interpolated_topology_base() : BaseTopology() { };
    
    template <typename A1>
    explicit
    interpolated_topology_base(const A1& a1) : 
                               BaseTopology(a1) { };
    
    template <typename A1, typename A2>
    interpolated_topology_base(const A1& a1, const A2& a2) : 
                               BaseTopology(a1, a2) { };
    
    template <typename A1, typename A2, typename A3>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3) : 
                               BaseTopology(a1, a2, a3) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                               BaseTopology(a1, a2, a3, a4) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                               BaseTopology(a1, a2, a3, a4, a5) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                               BaseTopology(a1, a2, a3, a4, a5, a6) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                               BaseTopology(a1, a2, a3, a4, a5, a6, a7) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    interpolated_topology_base(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
                               BaseTopology(a1, a2, a3, a4, a5, a6, a7, a8) { };
#endif
    
    /**
     * Returns a const-reference to the super-space of this topology.
     * \note This function returns a const-reference to itself since the super-space is also 
     *       the base-class of this topology. The base class is not polymorphic, meaning that its
     *       distance metric and random-sampler are not overridden (non-virtual calls).
     */
    const super_space_type& get_super_space() const { return *this; };
    
    /*************************************************************************
     *                             MetricSpaceConcept
     * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const {
      return this->interp_topo_get_distance(a,b);
    };
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b, double min_dist_interval, validity_predicate_type predicate) const {
      return this->interp_topo_get_distance_pred(a, b, min_dist_interval, predicate);
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      return this->interp_topo_get_norm(delta);
    };
    
    
    /**
     * Returns the distance between two points.
     */
    double proper_distance(const point_type& a, const point_type& b) const {
      return this->interp_topo_get_proper_distance(a,b);
    };
    
    /**
     * Returns the norm of the difference between two points.
     */
    double proper_norm(const point_difference_type& delta) const {
      return this->interp_topo_get_proper_norm(delta);
    };
    
    /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      return this->interp_topo_is_in_bounds(a);
    };
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed within the reachable space.
     */
    point_type random_point() const {
      return this->interp_topo_random_point();
    };
    
    /**
     * Generates a random point in the space, uniformly distributed within the reachable space.
     */
    point_type random_point(validity_predicate_type predicate) const {
      return this->interp_topo_random_point_pred(predicate);
    };
    
   /*************************************************************************
    *                             LieGroupConcept
    * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return this->interp_topo_move_position_toward(a, fraction, b);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b,
                                    double min_dist_interval, validity_predicate_type predicate) const {
      return this->interp_topo_move_position_toward_pred(a, fraction, b, min_dist_interval, predicate);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_back_to(const point_type& a, double fraction, const point_type& b) const {
      return this->interp_topo_move_position_back_to(a, fraction, b);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     */
    point_type move_position_back_to(const point_type& a, double fraction, const point_type& b,
                                     double min_dist_interval, validity_predicate_type predicate) const {
      return this->interp_topo_move_position_back_to_pred(a, fraction, b, min_dist_interval, predicate);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      BaseTopology::save(A,BaseTopology::getStaticObjectType()->TypeVersion());
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      BaseTopology::load(A,BaseTopology::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400039,1,"interpolated_topology_base",BaseTopology)
    
};

template <typename BaseTopology>
struct is_metric_space< interpolated_topology_base<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_reversible_space< interpolated_topology_base<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_point_distribution< interpolated_topology_base<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct get_rate_illimited_space< interpolated_topology_base<BaseTopology> > : get_rate_illimited_space< BaseTopology > { };

template <typename BaseTopology>
struct is_temporal_space< interpolated_topology_base<BaseTopology> > : is_temporal_space< BaseTopology > { };

template <typename BaseTopology>
struct is_metric_symmetric< interpolated_topology_base<BaseTopology> > : is_metric_symmetric< BaseTopology > { };

template <typename BaseTopology>
struct get_proper_metric< interpolated_topology_base<BaseTopology> > {
  typedef default_proper_metric type;
};



/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology>
class wrapped_interp_topology : public named_object {
  public:
    typedef wrapped_interp_topology<BaseTopology> self;
    
    typedef interpolated_topology_base<BaseTopology> wrapped_base_type;
    
    typedef typename topology_traits< wrapped_base_type >::point_type point_type;
    typedef typename topology_traits< wrapped_base_type >::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef BaseTopology super_space_type;
    
    typedef typename wrapped_base_type::space_topology space_topology;
    typedef typename wrapped_base_type::time_topology time_topology;
    
    typedef typename wrapped_base_type::validity_predicate_type validity_predicate_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< wrapped_base_type >::dimensions);
    
  private:
    
    shared_ptr< wrapped_base_type > m_space;
    
  public:
    
    explicit wrapped_interp_topology(const shared_ptr< wrapped_base_type >& aSpace = shared_ptr< wrapped_base_type >()) : 
                                     m_space(aSpace) { 
      setName("wrapped_interp_topology");
    };
    
    operator wrapped_base_type&() { return *m_space; };
    operator const wrapped_base_type&() const { return *m_space; };
    
    /**
     * Returns a const-reference to the super-space of this topology.
     * \note This function returns a const-reference to itself since the super-space is also 
     *       the base-class of this topology. The base class is not polymorphic, meaning that its
     *       distance metric and random-sampler are not overridden (non-virtual calls).
     */
    const super_space_type& get_super_space() const { return m_space->get_super_space(); };
    
    /** Returns the underlying space topology. */
    const space_topology& get_space_topology() const { return m_space->get_space_topology(); };
    /** Returns the underlying time topology. */
    const time_topology& get_time_topology() const { return m_space->get_time_topology(); };
    
    /*************************************************************************
     *                             TopologyConcept
     * **********************************************************************/
    
    /**
     * Returns the difference between two points (a - b).
     * \param a The first point.
     * \param b The second point.
     * \return The difference between the two points.
     */
    point_difference_type difference(const point_type& a, const point_type& b) const { return m_space->difference(a, b); };
    
    /**
     * Returns the addition of a point-difference to a point.
     * \param a The starting point.
     * \param delta The point-difference.
     * \return The addition of a point-difference to a point.
     */
    point_type adjust(const point_type& a, const point_difference_type& delta) const { return m_space->adjust(a, delta); };
    
    /**
     * Returns the origin of the temporal-space.
     * \return The origin of the temporal-space.
     */
    point_type origin() const { return m_space->origin(); };
    
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const { return m_space->is_in_bounds(a); };
    
    /*************************************************************************
     *                             MetricSpaceConcept
     * **********************************************************************/
    
    /**
     * Computes the distance between two points in the temporal-space.
     * \param a The first point.
     * \param b The second point.
     * \return The distance between a and b.
     */
    double distance(const point_type& a, const point_type& b) const { return m_space->distance(a, b); };
    
    double distance(const point_type& a, const point_type& b, double min_dist_interval, validity_predicate_type predicate) const {
      return m_space->distance(a, b, min_dist_interval, predicate);
    };
    
    /**
     * Computes the norm of a difference between two points.
     * \param a The difference between two points.
     * \return The norm of a difference between two points.
     */
    double norm(const point_difference_type& a) const { return m_space->norm(a); };
    
    
    /**
     * Computes the distance between two points in the temporal-space.
     * \param a The first point.
     * \param b The second point.
     * \return The distance between a and b.
     */
    double proper_distance(const point_type& a, const point_type& b) const { return m_space->proper_distance(a, b); };
    
    /**
     * Computes the norm of a difference between two points.
     * \param a The difference between two points.
     * \return The norm of a difference between two points.
     */
    double proper_norm(const point_difference_type& a) const { return m_space->proper_norm(a); };
    
    /*************************************************************************
     *                             LieGroupConcept
     * **********************************************************************/
    
    /**
     * Returns a point which is at a fraction between two points a to b.
     * \param a The first point.
     * \param fraction The fraction between the two points (0 to 1).
     * \param b The second point.
     * \return The point which is at a fraction between two points.
     */
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return m_space->move_position_toward(a, fraction, b);
    };
    
    point_type move_position_toward(const point_type& a, double fraction, const point_type& b,
                                    double min_dist_interval, validity_predicate_type predicate) const {
      return m_space->move_position_toward(a, fraction, b, min_dist_interval, predicate);
    };
    
    /**
     * Returns a point which is at a backward fraction between two points a to b.
     * \param a The first point.
     * \param fraction The backward fraction between the two points (0 to 1).
     * \param b The second point.
     * \return The point which is at a backward fraction between two points.
     */
    point_type move_position_back_to(const point_type& a, double fraction, const point_type& b) const {
      return m_space->move_position_back_to(a, fraction, b);
    };
    
    point_type move_position_back_to(const point_type& a, double fraction, const point_type& b,
                                     double min_dist_interval, validity_predicate_type predicate) const {
      return m_space->move_position_back_to(a, fraction, b, min_dist_interval, predicate);
    };
    
    /*************************************************************************
     *                             PointDistributionConcept
     * **********************************************************************/
    
    /**
     * Returns a random point within the temporal-space.
     * \return A random point within the temporal-space.
     */
    point_type random_point() const { return m_space->random_point(); };
    
    point_type random_point(validity_predicate_type predicate) const {
      return m_space->random_point(predicate);
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(m_space);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(m_space);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240003D,1,"wrapped_interp_topology",named_object)
    
};

template <typename BaseTopology>
struct is_metric_space< wrapped_interp_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_reversible_space< wrapped_interp_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct is_point_distribution< wrapped_interp_topology<BaseTopology> > : boost::mpl::true_ { };

template <typename BaseTopology>
struct get_rate_illimited_space< wrapped_interp_topology<BaseTopology> > : get_rate_illimited_space< BaseTopology > { };

template <typename BaseTopology>
struct is_temporal_space< wrapped_interp_topology<BaseTopology> > : is_temporal_space< BaseTopology > { };

template <typename BaseTopology>
struct is_metric_symmetric< wrapped_interp_topology<BaseTopology> > : is_metric_symmetric< BaseTopology > { };

template <typename BaseTopology>
struct get_proper_metric< wrapped_interp_topology<BaseTopology> > {
  typedef default_proper_metric type;
};






namespace detail {
  
  
  template <typename InterpMethodTag, typename BaseTopology>
  typename boost::disable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot) {
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return a;
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return b;
    if(fraction == 0.0)
      return a;
    
    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(), InterpFactoryType());
    double dist_inter = dist_tot * fraction;
    PointType result = a;
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter, dist_tot, InterpFactoryType());
    return result;
  };
  
  template <typename InterpMethodTag, typename BaseTopology>
  typename boost::enable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot) {
    
    if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
      return a; //b is not reachable from a.
    
    typedef typename temporal_space_traits< BaseTopology >::space_topology SpaceTopoType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return a;
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return b;
    if(fraction == 0.0)
      return a;
    
    InterpType interp;
    double dt_total = (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * fraction;
    PointType result = a;
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, InterpFactoryType());
    result.time = a.time + dt;
    return result;
  };
  
  
  
  template <typename InterpMethodTag, typename BaseTopology, typename PredFunction>
  typename boost::disable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot, double min_interval, PredFunction predicate) {
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return a;
    if(dist_tot < min_interval)
      return interp_topo_move_pt_impl<InterpMethodTag>(b_space, a, fraction, b, dist_tot);
    
    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(), InterpFactoryType());
    double dist_inter = dist_tot * fraction;
    double dist_cur = min_interval;
    PointType result = a;
    PointType last_result = a;
    while(dist_cur < dist_inter) {
      interp.compute_point(result, a, b, b_space, time_topology(), dist_cur, dist_tot, InterpFactoryType());
      if(!predicate(result))
        return last_result;
      dist_cur += min_interval;
      last_result = result;
    };
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return b;
    if(fraction == 0.0)
      return a;
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter, dist_tot, InterpFactoryType());
    return result;
  };
  
  template <typename InterpMethodTag, typename BaseTopology, typename PredFunction>
  typename boost::enable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot, double min_interval, PredFunction predicate) {
    
    if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
      return a; //b is not reachable from a.
    
    typedef typename temporal_space_traits< BaseTopology >::space_topology SpaceTopoType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return a;
    if(dist_tot < min_interval)
      return interp_topo_move_pt_impl<InterpMethodTag>(b_space, a, fraction, b, dist_tot);
    
    InterpType interp;
    double dt_total = (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * fraction;
    double d = min_interval;
    PointType result = a;
    PointType last_result = a;
    while(d < dt) {
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), d, dt_total, InterpFactoryType());
      result.time = a.time + d;
      if(!predicate(result))
        return last_result;
      d += min_interval;
      last_result = result;
    };
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return b;
    if(fraction == 0.0)
      return a;
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, InterpFactoryType());
    result.time = a.time + dt;
    return result;
  };
  
  
  
  
  template <typename InterpMethodTag, typename BaseTopology>
  typename boost::disable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot) {
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return b;
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return a;
    if(fraction == 0.0)
      return b;
    
    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(), InterpFactoryType());
    double dist_inter = dist_tot * (1.0 - fraction);
    PointType result = a;
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter, dist_tot, InterpFactoryType());
    return result;
  };
  
  template <typename InterpMethodTag, typename BaseTopology>
  typename boost::enable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot) {
    
    if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
      return b; //a is not reachable from b.
    
    typedef typename temporal_space_traits< BaseTopology >::space_topology SpaceTopoType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return b;
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return a;
    if(fraction == 0.0)
      return b;
    
    InterpType interp;
    double dt_total = (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * (1.0 - fraction);
    PointType result = a;
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, InterpFactoryType());
    result.time = a.time + dt;
    return result;
  };
  
  
  
  template <typename InterpMethodTag, typename BaseTopology, typename PredFunction>
  typename boost::disable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot, double min_interval, PredFunction predicate) {
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, BaseTopology, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return b;
    if(dist_tot < min_interval)
      return interp_topo_move_pt_back_impl<InterpMethodTag>(b_space, a, fraction, b, dist_tot);
    
    InterpType interp;
    interp.initialize(a, b, dist_tot, b_space, time_topology(), InterpFactoryType());
    double dist_inter = dist_tot * (1.0 - fraction);
    double dist_cur = dist_tot - min_interval;
    PointType result = b;
    PointType last_result = b;
    while(dist_cur > dist_inter) {
      interp.compute_point(result, a, b, b_space, time_topology(), dist_cur, dist_tot, InterpFactoryType());
      if(!predicate(result))
        return last_result;
      dist_cur -= min_interval;
      last_result = result;
    };
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return a;
    if(fraction == 0.0)
      return b;
    interp.compute_point(result, a, b, b_space, time_topology(), dist_inter, dist_tot, InterpFactoryType());
    return result;
  };
  
  template <typename InterpMethodTag, typename BaseTopology, typename PredFunction>
  typename boost::enable_if< is_temporal_space< BaseTopology >,
  typename topology_traits< BaseTopology >::point_type >::type interp_topo_move_pt_back_impl(
    const BaseTopology& b_space, 
    const typename topology_traits< BaseTopology >::point_type& a, 
    double fraction, 
    const typename topology_traits< BaseTopology >::point_type& b,
    double dist_tot, double min_interval, PredFunction predicate) {
    
    if(a.time > b.time) // Am I trying to go backwards in time (impossible)?
      return b; //a is not reachable from b.
    
    typedef typename temporal_space_traits< BaseTopology >::space_topology SpaceTopoType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::type InterpType;
    typedef typename get_tagged_spatial_interpolator< InterpMethodTag, SpaceTopoType, time_topology>::pseudo_factory_type InterpFactoryType;
    typedef typename topology_traits< BaseTopology >::point_type PointType;
    
    if(dist_tot == std::numeric_limits<double>::infinity())
      return b;
    if(dist_tot < min_interval)
      return interp_topo_move_pt_back_impl<InterpMethodTag>(b_space, a, fraction, b, dist_tot);
    
    InterpType interp;
    double dt_total = (b.time - a.time);  // the free time that I have along the path.
    interp.initialize(a.pt, b.pt, dt_total, b_space.get_space_topology(), b_space.get_time_topology(), InterpFactoryType());
    double dt = dt_total * (1.0 - fraction);
    double d = dt_total - min_interval;
    PointType result = b;
    PointType last_result = b;
    while(d > dt) {
      interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), d, dt_total, InterpFactoryType());
      result.time = a.time + d;
      if(!predicate(result))
        return last_result;
      d -= min_interval;
      last_result = result;
    };
    if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
      return a;
    if(fraction == 0.0)
      return b;
    interp.compute_point(result.pt, a.pt, b.pt, b_space.get_space_topology(), b_space.get_time_topology(), dt, dt_total, InterpFactoryType());
    result.time = a.time + dt;
    return result;
  };
  
  
  
  
};





/**
 * This class wraps an interpolated topology which is a topology with a new travel function, distance metric and sampler.
 * \tparam BaseTopology The topology underlying this space, should model TopologyConcept.
 */
template <typename BaseTopology, typename InterpMethodTag>
class interpolated_topology : public interpolated_topology_base<BaseTopology> {
  public:
    
    typedef interpolated_topology_base<BaseTopology> base_type;
    typedef interpolated_topology<BaseTopology, InterpMethodTag> self;
    
    typedef typename base_type::point_type point_type;
    typedef typename base_type::point_difference_type point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    typedef typename base_type::space_topology space_topology;
    typedef typename base_type::time_topology time_topology;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = base_type::dimensions);
    
    typedef typename base_type::validity_predicate_type validity_predicate_type;
    
  protected:
    
    virtual point_type interp_topo_move_position_toward(const point_type& a, double fraction, const point_type& b) const {
      return detail::interp_topo_move_pt_impl<InterpMethodTag, BaseTopology>(
        static_cast<const BaseTopology&>(*this), a, fraction, b, this->interp_topo_get_distance(a, b));
    };
    
    virtual point_type interp_topo_move_position_toward_pred(const point_type& a, double fraction, const point_type& b,
                                                        double min_dist_interval, validity_predicate_type predicate) const {
      return detail::interp_topo_move_pt_impl<InterpMethodTag, BaseTopology, validity_predicate_type>(
        static_cast<const BaseTopology&>(*this), a, fraction, b, this->interp_topo_get_distance(a, b), min_dist_interval, predicate);
    };
    
    virtual point_type interp_topo_move_position_back_to(const point_type& a, double fraction, const point_type& b) const {
      return detail::interp_topo_move_pt_back_impl<InterpMethodTag, BaseTopology>(
        static_cast<const BaseTopology&>(*this), a, fraction, b, this->interp_topo_get_distance(a, b));
    };
    
    virtual point_type interp_topo_move_position_back_to_pred(const point_type& a, double fraction, const point_type& b,
                                                              double min_dist_interval, validity_predicate_type predicate) const {
      return detail::interp_topo_move_pt_back_impl<InterpMethodTag, BaseTopology, validity_predicate_type>(
        static_cast<const BaseTopology&>(*this), a, fraction, b, this->interp_topo_get_distance(a, b), min_dist_interval, predicate);
    };
    
    
  public:
    
    interpolated_topology(const BaseTopology& aTopo) : base_type(aTopo) { };
    
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
    template <typename... Args>
    interpolated_topology(Args&&... args) : base_type(std::forward<Args>(args)...) { };
#else
    interpolated_topology() : base_type() { };
    
    template <typename A1>
    interpolated_topology(const A1& a1) : 
                          base_type(a1) { };
    
    template <typename A1, typename A2>
    interpolated_topology(const A1& a1, const A2& a2) : 
                          base_type(a1, a2) { };
    
    template <typename A1, typename A2, typename A3>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3) : 
                          base_type(a1, a2, a3) { };
    
    template <typename A1, typename A2, typename A3, typename A4>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4) : 
                          base_type(a1, a2, a3, a4) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5) : 
                          base_type(a1, a2, a3, a4, a5) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6) : 
                          base_type(a1, a2, a3, a4, a5, a6) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7) { };
    
    template <typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
    interpolated_topology(const A1& a1, const A2& a2, const A3& a3, const A4& a4, const A5& a5, const A6& a6, const A7& a7, const A8& a8) : 
                          base_type(a1, a2, a3, a4, a5, a6, a7, a8) { };
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

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240003A,1,"interpolated_topology",base_type)
    
};

template <typename BaseTopology, typename InterpMethodTag>
struct is_metric_space< interpolated_topology<BaseTopology, InterpMethodTag> > : boost::mpl::true_ { };

template <typename BaseTopology, typename InterpMethodTag>
struct is_reversible_space< interpolated_topology<BaseTopology, InterpMethodTag> > : boost::mpl::true_ { };

template <typename BaseTopology, typename InterpMethodTag>
struct is_point_distribution< interpolated_topology<BaseTopology, InterpMethodTag> > : boost::mpl::true_ { };

template <typename BaseTopology, typename InterpMethodTag>
struct get_rate_illimited_space< interpolated_topology<BaseTopology, InterpMethodTag> > : get_rate_illimited_space< BaseTopology > { };

template <typename BaseTopology, typename InterpMethodTag>
struct is_temporal_space< interpolated_topology<BaseTopology, InterpMethodTag> > : is_temporal_space< BaseTopology > { };

template <typename BaseTopology, typename InterpMethodTag>
struct is_metric_symmetric< interpolated_topology<BaseTopology, InterpMethodTag> > : 
  boost::mpl::and_<
    is_metric_symmetric< InterpMethodTag >,
    is_metric_symmetric< BaseTopology >
  > { };

template <typename BaseTopology, typename InterpMethodTag>
struct get_proper_metric< interpolated_topology<BaseTopology, InterpMethodTag> > {
  typedef default_proper_metric type;
};








};

};


namespace ReaK {
  
  
/* Specialization, see general template docs. */
  template <typename BaseTopology>
  struct arithmetic_tuple_size< pp::interpolated_topology_base<BaseTopology> > : 
    arithmetic_tuple_size< BaseTopology > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, pp::interpolated_topology_base<BaseTopology> > :
    arithmetic_tuple_element< Idx, BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const pp::interpolated_topology_base<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, volatile pp::interpolated_topology_base<BaseTopology> > : 
    arithmetic_tuple_element< Idx, volatile BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology>
  struct arithmetic_tuple_element< Idx, const volatile pp::interpolated_topology_base<BaseTopology> > : 
    arithmetic_tuple_element< Idx, const volatile BaseTopology > { };
  
  
  
/* Specialization, see general template docs. */
  template <typename BaseTopology, typename InterpMethodTag>
  struct arithmetic_tuple_size< pp::interpolated_topology<BaseTopology, InterpMethodTag> > : 
    arithmetic_tuple_size< BaseTopology > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology, typename InterpMethodTag>
  struct arithmetic_tuple_element< Idx, pp::interpolated_topology<BaseTopology, InterpMethodTag> > :
    arithmetic_tuple_element< Idx, BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology, typename InterpMethodTag>
  struct arithmetic_tuple_element< Idx, const pp::interpolated_topology<BaseTopology, InterpMethodTag> > : 
    arithmetic_tuple_element< Idx, const BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology, typename InterpMethodTag>
  struct arithmetic_tuple_element< Idx, volatile pp::interpolated_topology<BaseTopology, InterpMethodTag> > : 
    arithmetic_tuple_element< Idx, volatile BaseTopology > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename BaseTopology, typename InterpMethodTag>
  struct arithmetic_tuple_element< Idx, const volatile pp::interpolated_topology<BaseTopology, InterpMethodTag> > : 
    arithmetic_tuple_element< Idx, const volatile BaseTopology > { };
  
};





#endif









