/**
 * \file manip_free_workspace.hpp
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

#ifndef REAK_MANIP_FREE_WORKSPACE_HPP
#define REAK_MANIP_FREE_WORKSPACE_HPP


#include "path_planning/random_sampler_concept.hpp"
#include "path_planning/metric_space_concept.hpp"

#include "basic_distance_metrics.hpp"
#include "default_random_sampler.hpp"

#include "hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"

#include "interpolation/sustained_acceleration_pulse.hpp"
#include "interpolation/sustained_velocity_pulse.hpp"

#include "proximity/proxy_query_model.hpp"

#include "patj_planning/manipulator_topo_maps.hpp"
#include "joint_space_limits.hpp"

namespace ReaK {

namespace pp {


namespace detail {

template <typename RateLimitedJointSpace>
class manip_dk_proxy_env_impl {
  public:
    typedef manip_dk_proxy_env<RateLimitedJointSpace> self;
    typedef topology_traits< RateLimitedJointSpace >::point_type point_type;
    typedef topology_traits< RateLimitedJointSpace >::point_difference_type point_difference_type;
    
    shared_ptr< kte::manipulator_kinematics_model > m_model; 
    shared_ptr< joint_limits_collection<double> > m_joint_limits_map;
    
    std::vector< shared_ptr< geom::proxy_query_pair_2D > > m_proxy_env_2D;
    std::vector< shared_ptr< geom::proxy_query_pair_3D > > m_proxy_env_3D;
    
    manip_dk_proxy_env(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
                       const shared_ptr< joint_limits_collection<double> >& aJointLimitMap = shared_ptr< joint_limits_collection<double> >()) :
                       m_model(aModel), m_joint_limits_map(aJointLimitMap) { };
    
    bool is_free(const point_type& pt, const RateLimitedJointSpace& space) const {
      typedef typename get_rate_illimited_space< RateLimitedJointSpace >::type NormalJointSpace;
      NormalJointSpace normal_j_space;
      typename topology_traits<NormalJointSpace>::point_type pt_inter = m_joint_limits_map->map_to_space(pt, space, normal_j_space);
      detail::write_joint_coordinates_impl(pt_inter, normal_j_space, m_model);
      
      // update the kinematics model with the given joint states.
      m_model->doMotion();
      
      // NOTE: it is assumed that the proxy environments are properly linked with the manipulator kinematic model.
      
      for( std::vector< shared_ptr< geom::proxy_query_pair_2D > >::const_iterator it = m_proxy_env_2D.begin(); it != m_proxy_env_2D.end(); ++it) {
        double tmp = (*it)->findMinimumDistance();
        if(tmp < 0.0)
          return false;
      };
      for( std::vector< shared_ptr< geom::proxy_query_pair_3D > >::const_iterator it = m_proxy_env_3D.begin(); it != m_proxy_env_3D.end(); ++it) {
        double tmp = (*it)->findMinimumDistance();
        if(tmp < 0.0)
          return false;
      };
      
      return true;
    };
    
    
    
};

};


/**
 * This class is used to represent the free-space of a manipulator in a quasi-static environment.
 * Here, the term quasi-static refers to the fact that proximity queries are performed against an 
 * environment that is assumed to be unchanging (at least, not significantly while this topology is used).
 */
template <typename RateLimitedJointSpace, typename InterpMethodTag>
class manip_quasi_static_env : public named_object {
  public:
    typedef manip_quasi_static_env<RateLimitedJointSpace,InterpMethodTag> self;
    typedef RateLimitedJointSpace super_space_type;
    typedef topology_traits< super_space_type >::point_type point_type;
    typedef topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    double min_interval;
    double max_edge_length;
    
    super_space_type m_space;
    typename metric_space_traits<super_space_type>::distance_metric_type m_distance;
    typename point_distribution_traits<super_space_type>::random_sampler_type m_rand_sampler;
    
    detail::manip_dk_proxy_env_impl<RateLimitedJointSpace> m_prox_env;
    
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
    
    
    /**
     * Checks if the given point is within the free-space.
     * \param p The point to be checked for being collision-free.
     * \return True if p is collision-free.
     */
    bool is_free(const point_type& p) const {
      return m_prox_env.is_free(pt, m_space);
    };
    
    //Topology concepts:
    
    /**
     * Produces a random, collision-free point.
     * \return A random, collision-free point.
     */
    point_type random_point() const {
      point_type result;
      while(!m_prox_env.is_free(result = m_rand_sampler(m_space), m_space)) ; //output only free C-space points.
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
      if(m_distance(p2, move_position_toward(p1, 1.0, p2), m_space) < std::numeric_limits< double >::epsilon())
        return m_distance(p1, p2, m_space); //if p2 is reachable from p1, use Euclidean distance.
      else
        return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1.
    };
    
    /**
     * Computes the norm of the difference between two points. 
     * \param dp The point difference.
     * \return The norm of the difference between the two points.
     */
    double norm(const point_difference_type& dp) const {
      return m_distance(dp, m_space);
    };
    
    /**
     * Returns the difference between two points (a - b).
     */
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1,p2);
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
      return m_space.adjust(p,dp);
    };
    
    /**
     * Returns a point which is at a fraction between two points a to b, or as 
     * far as it can get before a collision.
     * \todo TODO: Write this using the interpolator corresponding to the given tag.
     */
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      double dist = m_distance(p1, p2, m_space);
      double d = min_interval;
      while(d < dist * fraction) {
        if (!m_prox_env.is_free(m_space.move_position_toward(p1, (d / dist), p2))) {
          return m_space.move_position_toward(p1,((d - min_interval) / dist), p2);
        };
        d += min_interval;
      };
      if(fraction == 1.0) //these equal comparison are used for when exact end fractions are used.
        return p2;
      else if(fraction == 0.0)
        return p1;
      else 
        return m_space.move_position_toward(p1, fraction, p2);
    };
    
    /**
     * Returns a random point fairly near to the given point.
     */
    std::pair<point_type, bool> random_walk(const point_type& p_u) const {
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = m_rand_sampler(m_space);
        double dist = m_distance(p_u, p_rnd, m_space);
        p_v = move_position_toward(p_u, max_edge_length / dist, p_rnd);
        ++i;
      } while((m_distance(p_u, p_v, m_space) < min_interval) && (i <= 10));
      if(i > 10) {
        //could not expand vertex u, then just output a random C-free point.
        return std::make_pair(p_v, false);
      };
      return std::make_pair(p_v, true);
    };
    
    
    /**
     * Parametrized constructor (this class is a RAII class).
     * \param aMaxEdgeLength The maximum length of an added edge, in time units (e.g., seconds).
     */
    manip_quasi_static_env(double aMinInterval = 0.1, 
                           double aMaxEdgeLength = 1.0) : 
                           min_interval(aMinInterval),
                           max_edge_length(aMaxEdgeLength),
                           m_space(),
                           m_distance(get(distance_metric, m_space)),
                           m_rand_sampler(get(random_sampler, m_space)), 
                           m_prox_env() { };
    
    virtual ~manip_quasi_static_env() { };
    
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
        & RK_SERIAL_SAVE_WITH_NAME(m_prox_env.m_proxy_env_3D);
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
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400027,1,"manip_quasi_static_env",named_object)
    
    
};


};

};

#endif

