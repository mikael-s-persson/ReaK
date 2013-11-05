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

#include "proxy_model_updater.hpp"       // needed by manip_dynamic_env
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
template <typename BaseJointSpace>
class manip_dynamic_env : public named_object {
  public:
    typedef manip_dynamic_env<BaseJointSpace> self;
    typedef temporal_space<BaseJointSpace, time_poisson_topology, reach_plus_time_metric> base_temporal_joint_space;
    typedef interpolated_topology_base<base_temporal_joint_space> super_space_base_type;
    typedef wrapped_interp_topology<base_temporal_joint_space> super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    typedef time_poisson_topology time_topology;
    typedef BaseJointSpace space_topology;
    
    typedef proxy_model_applicator<BaseJointSpace> applicator_type;
    typedef detail::manip_dk_proxy_env_impl<BaseJointSpace> dk_proxy_env_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    double min_interval;
    super_space_type m_space;
    
    dk_proxy_env_type m_prox_env;
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
    
    /**
     * Checks if the given point is within the free-space.
     * \param p The point to be checked for being collision-free.
     * \return True if p is collision-free.
     */
    bool is_free(const point_type& p) const {
      if(!m_space.is_in_bounds(p))
        return false;
      for(std::size_t i = 0; i < m_prox_updaters.size(); ++i)
        m_prox_updaters[i]->synchronize_proxy_model(p.time);
      return m_prox_env.is_free(p.pt, m_space.get_space_topology());
    };
    
    struct is_free_predicate {
      const self* parent;
      is_free_predicate(const self* aParent) : parent(aParent) { };
      bool operator()(const point_type& p) const { return parent->is_free(p); };
    };
    
    //Topology concepts:
    
    /**
     * Produces a random, collision-free point.
     * \return A random, collision-free point.
     */
    point_type random_point() const {
      return m_space.random_point( is_free_predicate(this) ); 
    };
    
    /**
     * Computes the distance between two points. If there is no collision-free line between
     * the two points, the distance is infinite.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return The collision-free distance between the two given points.
     */
    double distance(const point_type& p1, const point_type& p2) const {
      return m_space.distance(p1, p2, min_interval, is_free_predicate(this));
    };
    
    /**
     * Computes the norm of the difference between two points. 
     * \param dp The point difference.
     * \return The norm of the difference between the two points.
     */
    double norm(const point_difference_type& dp) const { return m_space.norm(dp); };
    
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
      return m_space.move_position_toward(p1, fraction, p2, min_interval, is_free_predicate(this));
    };
    
    
    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aMaxEdgeLength The maximum length of an added edge, in time units (e.g., seconds).
     */
    explicit
    manip_dynamic_env(const BaseJointSpace& aSpace = BaseJointSpace(),
                      const shared_ptr< applicator_type >& aApplicator = shared_ptr< applicator_type >(),
                      double aMinInterval = 0.1, double aMaxEdgeLength = 10.0) : 
                      min_interval(aMinInterval),
                      m_space(
                        shared_ptr<super_space_base_type>(new super_space_base_type(
                          base_temporal_joint_space(
                            "manip_dynamic_env_underlying_space", 
                            aSpace, time_poisson_topology("time-poisson topology", aMinInterval, aMaxEdgeLength)
                          )))),
                      m_prox_env(aApplicator) { };
    
    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aMinInterval The minimum length of the travel between collision detection calls.
     */
    template <typename InterpMethodTag>
    explicit
    manip_dynamic_env(InterpMethodTag aInterpTag,
                      const BaseJointSpace& aSpace = BaseJointSpace(),
                      const shared_ptr< applicator_type >& aApplicator = shared_ptr< applicator_type >(),
                      double aMinInterval = 0.1, double aMaxEdgeLength = 10.0) : 
                      min_interval(aMinInterval),
                      m_space(
                        shared_ptr<super_space_base_type>(new interpolated_topology<base_temporal_joint_space, InterpMethodTag>(
                          base_temporal_joint_space(
                            "manip_dynamic_env_underlying_space", 
                            aSpace, time_poisson_topology("time-poisson topology", aMinInterval, aMaxEdgeLength)
                          )))),
                      m_prox_env(aApplicator) { };
    
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
        & RK_SERIAL_SAVE_WITH_NAME(m_space)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_applicator)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_3D)
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_updaters);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(min_interval)
        & RK_SERIAL_LOAD_WITH_NAME(m_space)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_applicator)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_updaters);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400028,1,"manip_dynamic_env",named_object)
    
    
};


template <typename BaseJointSpace>
struct is_metric_space< manip_dynamic_env<BaseJointSpace> > : boost::mpl::true_ { };
        
template <typename BaseJointSpace>
struct is_point_distribution< manip_dynamic_env<BaseJointSpace> > : boost::mpl::true_ { };

template <typename BaseJointSpace>
struct is_temporal_space< manip_dynamic_env<BaseJointSpace> > : boost::mpl::true_ { };



};

};



#include "manip_free_dynamic_workspace_ext.hpp"


#endif

