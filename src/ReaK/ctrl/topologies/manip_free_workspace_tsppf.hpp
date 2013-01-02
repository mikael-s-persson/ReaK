
// Template-specialization pre-processor-based factories.


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


#ifdef RK_GENERATE_MQSENV_REACHINTERP

template <typename RateLimitedJointSpace>
class manip_quasi_static_env<RateLimitedJointSpace, RK_REACHINTERP_TAG> : public named_object {
  public:
    typedef manip_quasi_static_env<RateLimitedJointSpace, RK_REACHINTERP_TAG> self;
    typedef RK_REACHINTERP_TOPOLOGY<RateLimitedJointSpace> super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    double min_interval;
    double max_edge_length;
    
    super_space_type m_space;
    
    detail::manip_dk_proxy_env_impl m_prox_env;
    
  public:
    
    super_space_type& get_super_space() { return m_space; };
    const super_space_type& get_super_space() const { return m_space; };
    
    bool is_free(const point_type& p) const {
      return m_prox_env.is_free(p, m_space.get_super_space());
    };
    
    //Topology concepts:
    
    point_type random_point() const {
      point_type result;
      while(!m_prox_env.is_free(result = m_space.random_point(), m_space.get_super_space())) ; //output only free C-space points.
      return result;
    };
    
    double distance(const point_type& p1, const point_type& p2) const {
      if(m_space.distance(p2, move_position_toward(p1, 1.0, p2)) < std::numeric_limits< double >::epsilon())
        return m_space.distance(p1, p2); //if p2 is reachable from p1, use reach-time distance.
      else
        return std::numeric_limits<double>::infinity(); //p2 is not reachable from p1.
    };
    
    double norm(const point_difference_type& dp) const { return m_space.norm(dp); };
    
    point_difference_type difference(const point_type& p1, const point_type& p2) const { return m_space.difference(p1,p2); };
    
    point_type origin() const { return m_space.origin(); };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const { return m_space.adjust(p,dp); };
    
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      typedef typename get_tagged_spatial_interpolator< RK_REACHINTERP_TAG, RateLimitedJointSpace, time_topology>::type InterpType;
      InterpType interp;
      interp.initialize(p1, p2, 0.0, static_cast<const RateLimitedJointSpace&>(m_space), time_topology(), m_space.get_pseudo_factory());
      double dt_min = interp.get_minimum_travel_time();
      double dt = dt_min * fraction;
      dt = (dt < max_edge_length ? dt : max_edge_length);
      double d = min_interval;
      point_type result = p1;
      point_type last_result = p1;
      while(d < dt) {
        interp.compute_point(result, p1, p2, static_cast<const RateLimitedJointSpace&>(m_space), time_topology(), d, dt_min, m_space.get_pseudo_factory());
        if(!m_prox_env.is_free(result, m_space.get_super_space()))
          return last_result;
        d += min_interval;
        last_result = result;
      };
      if((fraction == 1.0) && (dt_min < max_edge_length)) //these equal comparison are used for when exact end fractions are used.
        return p2;
      else if(fraction == 0.0)
        return p1;
      else {
        interp.compute_point(result, p1, p2, static_cast<const RateLimitedJointSpace&>(m_space), time_topology(), dt, dt_min, m_space.get_pseudo_factory());
        return result;
      };
    };
    
    std::pair<point_type, bool> random_walk(const point_type& p_u) const {
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = m_space.random_point();
        double dist = m_space.distance(p_u, p_rnd);
        p_v = move_position_toward(p_u, max_edge_length / dist, p_rnd);
        ++i;
      } while((m_space.distance(p_u, p_v) < min_interval) && (i <= 10));
      if(i > 10) {
        //could not expand vertex u, then just output a random C-free point.
        return std::make_pair(p_v, false);
      };
      return std::make_pair(p_v, true);
    };
    
    manip_quasi_static_env(const RateLimitedJointSpace& aSpace = RateLimitedJointSpace(),
                           const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >(),
                           const shared_ptr< joint_limits_collection<double> >& aJointLimitsMap = shared_ptr< joint_limits_collection<double> >(),
                           double aMinInterval = 0.1, 
                           double aMaxEdgeLength = 1.0,
                           double aInterpTolerance = 1e-6, 
                           unsigned int aInterpMaxIter = 60) : 
                           min_interval(aMinInterval),
                           max_edge_length(aMaxEdgeLength),
                           m_space(aSpace, aInterpTolerance, aInterpMaxIter),
                           m_prox_env(aModel, aJointLimitsMap) { };
    
    virtual ~manip_quasi_static_env() { };
    
    self& operator<<(const shared_ptr< geom::proxy_query_pair_2D >& aProxy) {
      m_prox_env.m_proxy_env_2D.push_back(aProxy);
      return *this;
    };
    
    self& operator<<(const shared_ptr< geom::proxy_query_pair_3D >& aProxy) {
      m_prox_env.m_proxy_env_3D.push_back(aProxy);
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
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_model)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_joint_limits_map)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400027,1,"manip_quasi_static_env",named_object)
    
    
};

#endif







#ifdef RK_GENERATE_MDENV_REACHINTERP


template <typename RateLimitedJointSpace>
class manip_dynamic_env<RateLimitedJointSpace, RK_REACHINTERP_TAG> : public named_object {
  public:
    typedef manip_dynamic_env<RateLimitedJointSpace, RK_REACHINTERP_TAG> self;
    typedef temporal_space< RK_REACHINTERP_TOPOLOGY<RateLimitedJointSpace>, time_poisson_topology, reach_plus_time_metric> super_space_type;
    typedef typename topology_traits< super_space_type >::point_type point_type;
    typedef typename topology_traits< super_space_type >::point_difference_type point_difference_type;
    
    typedef time_poisson_topology time_topology;
    typedef RK_REACHINTERP_TOPOLOGY<RateLimitedJointSpace> space_topology;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = topology_traits< super_space_type >::dimensions);
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
  private:
    double min_interval;
    double max_edge_length;
    
    super_space_type m_space;
    
    detail::manip_dk_proxy_env_impl m_prox_env;
    std::vector< shared_ptr< proxy_model_updater > > m_prox_updaters;
    
  public:
    
    super_space_type& get_super_space() { return m_space; };
    
    const super_space_type& get_super_space() const { return m_space; };
    
    const space_topology& get_space_topology() const { return m_space.get_space_topology(); };
    const time_topology& get_time_topology() const { return m_space.get_time_topology(); };
    
    space_topology& get_space_topology() { return m_space.get_space_topology(); };
    time_topology& get_time_topology() { return m_space.get_time_topology(); };
    
    bool is_free(const point_type& p) const {
      for(std::size_t i = 0; i < m_prox_updaters.size(); ++i)
        m_prox_updaters[i]->synchronize_proxy_model(p.time);
      return m_prox_env.is_free(p.pt, m_space.get_space_topology());
    };
    
    //Topology concepts:
    
    point_type random_point() const {
      point_type result;
      while(!is_free(result = point_type(m_space.get_time_topology().random_point(), m_space.get_space_topology().random_point()))) ; //output only free C-space points.
      return result;
    };
    
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
    
    double norm(const point_difference_type& dp) const {
      return get(distance_metric, m_space)(dp, m_space);
    };
    
    point_difference_type difference(const point_type& p1, const point_type& p2) const {
      return m_space.difference(p1, p2);
    };
    
    point_type origin() const {
      return m_space.origin();
    };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return m_space.adjust(p, dp);
    };
    
    point_type move_position_toward(const point_type& p1, double fraction, const point_type& p2) const {
      if(p1.time > p2.time) // Am I trying to go backwards in time (impossible)?
        return p1; //p2 is not reachable from p1.
        
      typedef typename get_tagged_spatial_interpolator< RK_REACHINTERP_TAG, RateLimitedJointSpace, time_topology>::type InterpType;
      InterpType interp;
      double reach_time = m_space.get_space_topology().distance(p1.pt, p2.pt);
      double dt_total = (p2.time - p1.time);  // the free time that I have along the path.
      if(dt_total < reach_time) // There is not enough time to reach the end-point.
        return p1;
      interp.initialize(p1.pt, p2.pt, (p2.time - p1.time), static_cast<const RateLimitedJointSpace&>(m_space.get_space_topology()), m_space.get_time_topology(), m_space.get_space_topology().get_pseudo_factory());
      double dt = dt_total * fraction;
      dt = (dt < max_edge_length ? dt : max_edge_length);
      double d = min_interval;
      point_type result = p1;
      point_type last_result = p1;
      while(d < dt) {
        interp.compute_point(result.pt, p1.pt, p2.pt, static_cast<const RateLimitedJointSpace&>(m_space.get_space_topology()), m_space.get_time_topology(), d, dt_total, m_space.get_space_topology().get_pseudo_factory());
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
        interp.compute_point(result.pt, p1.pt, p2.pt, static_cast<const RateLimitedJointSpace&>(m_space.get_space_topology()), m_space.get_time_topology(), dt, dt_total, m_space.get_space_topology().get_pseudo_factory());
        return result;
      };
    };
    
    std::pair<point_type, bool> random_walk(const point_type& p_u) const {
      point_type p_rnd, p_v;
      unsigned int i = 0;
      do {
        p_rnd = point_type(0.0, m_space.get_space_topology().random_point());
        double reach_time = m_space.get_space_topology().distance(p_u.pt, p_rnd.pt);
        p_rnd.time = m_space.get_time_topology().random_point() + reach_time + p_u.time;
        p_v = move_position_toward(p_u, max_edge_length / reach_time, p_rnd);
        ++i;
      } while((p_v.time - p_u.time < min_interval) && (i <= 20));
      if(i > 20) {
        //could not expand vertex u, then just output a random C-free point.
        return std::make_pair(p_v, false);
      };
      return std::make_pair(p_v, true);
    };
    
    
    manip_dynamic_env(const RateLimitedJointSpace& aSpace = RateLimitedJointSpace(),
                      const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr< kte::direct_kinematics_model >(),
                      const shared_ptr< joint_limits_collection<double> >& aJointLimitsMap = shared_ptr< joint_limits_collection<double> >(),
                      double aMinInterval = 0.1, 
                      double aMaxEdgeLength = 1.0,
                      double aInterpTolerance = 1e-6, 
                      unsigned int aInterpMaxIter = 60) : 
                      min_interval(aMinInterval),
                      max_edge_length(aMaxEdgeLength),
                      m_space("manip_dynamic_env_underlying_space", 
                              space_topology(aSpace, aInterpTolerance, aInterpMaxIter), 
                              time_poisson_topology("time-poisson topology", aMinInterval, aMaxEdgeLength)),
                      m_prox_env(aModel, aJointLimitsMap) { };
    
    virtual ~manip_dynamic_env() { };
    
    self& operator<<(const shared_ptr< geom::proxy_query_pair_2D >& aProxy) {
      m_prox_env.m_proxy_env_2D.push_back(aProxy);
      return *this;
    };
    
    self& operator<<(const shared_ptr< geom::proxy_query_pair_3D >& aProxy) {
      m_prox_env.m_proxy_env_3D.push_back(aProxy);
      return *this;
    };
    
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
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_model)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_joint_limits_map)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_2D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_env.m_proxy_env_3D)
        & RK_SERIAL_LOAD_WITH_NAME(m_prox_updaters);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400028,1,"manip_dynamic_env",named_object)
    
    
};


#endif






