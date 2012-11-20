/**
 * \file manip_free_dynamic_workspace.hpp
 * 
 * This library defines a class for path-planning for a manipulator moving inside an environment 
 * with obstacles and physical limits (joint limits). This class is essentially just an assembly 
 * of many of the building blocks in the ReaK path-planning library. This class can also draw the 
 * elements of the motion graph (edges) as end-effector trajectories in a Coin3D scene-graph.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_MANIP_FREE_DYNAMIC_WORKSPACE_HPP
#define REAK_MANIP_FREE_DYNAMIC_WORKSPACE_HPP

#include "base/defs.hpp"
#include <boost/config.hpp>

#include "manip_free_workspace.hpp"

#include "path_planning/proxy_model_updater.hpp"       // needed by manip_dynamic_env
#include "temporal_space.hpp"         // needed by manip_dynamic_env
#include "time_poisson_topology.hpp"  // needed by manip_dynamic_env
#include "reachability_space.hpp"     // needed by manip_dynamic_env

namespace ReaK {

namespace pp {


/**
 * This class is used to represent the free-space of a manipulator in a dynamic environment.
 * Here, the term dynamic means that the space is a temporal space (time-space tuple) and proximity
 * queries are performed against an environment updated to the time-point of the current temporal 
 * point in question.
 */
template <typename RateLimitedJointSpace, typename InterpMethodTag>
class manip_dynamic_env : public named_object {
  public:
    typedef manip_dynamic_env<RateLimitedJointSpace,InterpMethodTag> self;
    typedef temporal_space<RateLimitedJointSpace, time_poisson_topology, reach_plus_time_metric> super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    typedef time_poisson_topology time_topology;
    typedef RateLimitedJointSpace space_topology;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    double min_interval;
    double max_edge_length;
    
    super_space_type m_space;
    typename metric_space_traits<RateLimitedJointSpace>::distance_metric_type m_distance;
    typename point_distribution_traits<RateLimitedJointSpace>::random_sampler_type m_rand_sampler;
    
    detail::manip_dk_proxy_env_impl m_prox_env;
    std::vector< shared_ptr< proxy_model_updater > > m_prox_updaters;
    
  public:
    
    /**
     * Returns a reference to the super-space in which this test-world is embedded.
     * \return A reference to the super-space in which this test-world is embedded.
     */
    super_space_type& get_super_space() { return m_space; };
    
    /**
     * Returns a const-reference to the super-space in which this test-world is embedded.
     * \return A const-reference to the super-space in which this test-world is embedded.
     */
    const super_space_type& get_super_space() const { return m_space; };
    
    /** Returns the underlying space topology. */
    const space_topology& get_space_topology() const { return m_space.get_space_topology(); };
    /** Returns the underlying time topology. */
    const time_topology& get_time_topology() const { return m_space.get_time_topology(); };
    
    /** Returns the underlying space topology. */
    space_topology& get_space_topology() { return m_space.get_space_topology(); };
    /** Returns the underlying time topology. */
    time_topology& get_time_topology() { return m_space.get_time_topology(); };
    
    /**
     * Checks if the given point is within the free-space.
     * \param p The point to be checked for being collision-free.
     * \return True if p is collision-free.
     */
    bool is_free(const point_type& p) const {
      for(std::size_t i = 0; i < m_prox_updaters.size(); ++i)
        m_prox_updaters[i]->synchronize_proxy_model(p.time);
      return m_prox_env.is_free(p.pt, m_space);
    };
    
    //Topology concepts:
    
    /**
     * Produces a random, collision-free point.
     * \return A random, collision-free point.
     */
    point_type random_point() const {
      point_type result;
      while(!is_free(result = point_type(m_space.get_time_topology().random_point(), m_rand_sampler(m_space.get_space_topology())))) ; //output only free C-space points.
      return result;
    };
    
    /**
     * Computes the distance between two points. If there is no collision-free line between
     * the two points, the distance is infinite.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return The collision-free distance between the two given points.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      using std::fabs;
      
      double actual_dist = get(distance_metric, m_space)(p1, p2, m_space);
      if(actual_dist == std::numeric_limits<double>::infinity())
        return actual_dist;
      
      if(fabs(p2.time - move_position_toward(p1, 1.0, p2).time) < std::numeric_limits< double >::epsilon())
        return actual_dist; //if p2 is reachable from p1, use Euclidean distance.
      else
        return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1, due to a collision.
    };
    
    /**
     * Computes the norm of the difference between two points. 
     * \param dp The point difference.
     * \return The norm of the difference between the two points.
     */
    double norm(const point_difference_type& dp) const {
      return get(distance_metric, m_space)(dp, m_space);
    };
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1, p2);
    };
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type origin() const {
      return m_space.origin();
    };
    
    /**
     * Returns the addition of a point-difference to a point.
     */
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return m_space.adjust(p, dp);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b, or as 
     * far as it can get before a collision.
     */
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      if(p1.time > p2.time) // Am I trying to go backwards in time (impossible)?
        return p1; //p2 is not reachable from p1.
      
      typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::type InterpType;
      typedef typename get_tagged_spatial_interpolator< InterpMethodTag, RateLimitedJointSpace, time_topology>::pseudo_factory_type InterpFactoryType;
      
      InterpType interp;
      double reach_time = m_distance(p1.pt, p2.pt, m_space.get_space_topology());
      double dt_total = (p2.time - p1.time);  // the free time that I have along the path.
      if(dt_total < reach_time) // There is not enough time to reach the end-point.
        return p1;
      interp.initialize(p1.pt, p2.pt, (p2.time - p1.time), m_space.get_space_topology(), m_space.get_time_topology(), InterpFactoryType());
      double dt = dt_total * fraction;
      dt = (dt < max_edge_length ? dt : max_edge_length);
      double d = min_interval;
      point_type result = p1;
      point_type last_result = p1;
      while(d < dt) {
        interp.compute_point(result.pt, p1.pt, p2.pt, m_space.get_space_topology(), m_space.get_time_topology(), d, dt_total, InterpFactoryType());
        result.time = p1.time + d;
        if(!is_free(result))
          return last_result;
        d += min_interval;
        last_result = result;
      };
      if((fraction == 1.0) && (dt_total < max_edge_length)) //these equal comparison are used for when exact end fractions are used.
        return p2;
      else if(fraction == 0.0)
        return p1;
      else {
        interp.compute_point(result.pt, p1.pt, p2.pt, m_space.get_space_topology(), m_space.get_time_topology(), dt, dt_total, InterpFactoryType());
        result.time = p1.time + dt;
        return result;
      };
    };
    
    /**
     * Returns a random point fairly near to the given point.
     */
    std::pair<point_type, bool> random_walk(const point_type& p_u) const {
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = point_type(0.0, m_rand_sampler(m_space.get_space_topology()));
        double reach_time = m_distance(p_u.pt, p_rnd.pt, m_space.get_space_topology());
        p_rnd.time = m_space.get_time_topology().random_point() + reach_time + p_u.time;
        p_v = move_position_toward(p_u, max_edge_length / reach_time, p_rnd);
        ++i;
      } while((p_v.time - p_u.time < min_interval) && (i <= 20));
      if(i > 20) {
        //could not expand vertex u, then just output a arbitrary C-free point.
        return std::make_pair(p_v, false);
      };
      return std::make_pair(p_v, true);
    };
    
    
    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aMaxEdgeLength The maximum length of an added edge, in time units (e.g., seconds).
     */
    manip_dynamic_env(const super_space_type& aSpace = super_space_type(),
                      const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >(),
                      const shared_ptr< joint_limits_collection<double> >& aJointLimitsMap = shared_ptr< joint_limits_collection<double> >(),
                      double aMinInterval = 0.1, 
                      double aMaxEdgeLength = 1.0) : 
                      min_interval(aMinInterval),
                      max_edge_length(aMaxEdgeLength),
                      m_space("manip_dynamic_env_underlying_space", 
                              aSpace, 
                              time_poisson_topology("time-poisson topology", aMinInterval, aMaxEdgeLength)),
                      m_distance(get(distance_metric, m_space.get_space_topology())),
                      m_rand_sampler(get(random_sampler, m_space.get_space_topology())), 
                      m_prox_env(aModel, aJointLimitsMap) { };
    
    virtual ~manip_dynamic_env() { };
    
    /**
     * Add a 2D proxy query pair to the collision environment.
     * \param aProxy The new 2D proxy query pair to add to the collision environment.
     * \return A reference back to 'this'.
     */
    self& operator<<(const shared_ptr< geom::proxy_query_pair_2D >& aProxy) {
      m_prox_env.m_proxy_env_2D.push_back(aProxy);
      return *this;
    };
    
    /**
     * Add a 3D proxy query pair to the collision environment.
     * \param aProxy The new 3D proxy query pair to add to the collision environment.
     * \return A reference back to 'this'.
     */
    self& operator<<(const shared_ptr< geom::proxy_query_pair_3D >& aProxy) {
      m_prox_env.m_proxy_env_3D.push_back(aProxy);
      return *this;
    };
    
    /**
     * Add a functor to update the proxy-query models for a given time value (parameter to functor call).
     * This builds a list of functors that will be called just before each proximity-queries in order to make sure 
     * the geometric models on which the query is performed are in their correct dynamic state.
     * Typically, such an updater would obtain the state from some trajectory (predicted or controlled)
     * and then place the components of the geometric model in the correct configuration.
     * \param aFunc A functor to add to the list of proxy model updaters.
     * \return A reference back to 'this'.
     */
    self& add_proxy_model_updater(const shared_ptr< proxy_model_updater >& aUpdater) {
      m_prox_updaters.push_back(aUpdater);
      return *this;
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(min_interval)
        & RK_SERIAL_SAVE_WITH_NAME(max_edge_length)
        & RK_SERIAL_SAVE_WITH_NAME(m_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_distance)
        & RK_SERIAL_SAVE_WITH_NAME(m_rand_sampler)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_model)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_joint_limits_map)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_3D)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_updaters);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(min_interval)
        & RK_SERIAL_LOAD_WITH_NAME(max_edge_length)
        & RK_SERIAL_LOAD_WITH_NAME(m_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_distance)
        & RK_SERIAL_LOAD_WITH_NAME(m_rand_sampler)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_model)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_joint_limits_map)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_updaters);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400028,1,"manip_dynamic_env",named_object)
    
    
};





#define RK_GENERATE_MDENV_REACHINTERP

// Create a manip_dynamic_env specialization for the SAP-interpolator
#define RK_REACHINTERP_TAG sap_interpolation_tag
#define RK_REACHINTERP_TOPOLOGY sap_reach_topology

#include "manip_free_workspace_tsppf.hpp"

#undef RK_REACHINTERP_TAG
#undef RK_REACHINTERP_TOPOLOGY

// Create a manip_dynamic_env specialization for the SVP-interpolator
#define RK_REACHINTERP_TAG svp_interpolation_tag
#define RK_REACHINTERP_TOPOLOGY svp_reach_topology

#include "manip_free_workspace_tsppf.hpp"

#undef RK_REACHINTERP_TAG
#undef RK_REACHINTERP_TOPOLOGY

#undef RK_GENERATE_MDENV_REACHINTERP


template <typename RateLimitedJointSpace, typename InterpMethodTag>
struct is_metric_space< manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag> > : boost::mpl::true_ { };
        
template <typename RateLimitedJointSpace, typename InterpMethodTag>
struct is_point_distribution< manip_dynamic_env<RateLimitedJointSpace, InterpMethodTag> > : boost::mpl::true_ { };





};

};

#endif

