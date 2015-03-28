/**
 * \file fehlberg45_integrator_sys.hpp
 *
 * This library implements an integrator that uses the 4-5th Order Fehlberg Method. This
 * method is of the Runge-Kutta-Fehlberg type which are single-step variable-step algorithm for numerical
 * integration in which the error is estimated by the difference between two integration rules (order 4 and 5)
 * and the time-step is adapted on-the-fly to control the norm of that error estimate.
 * These variable-step integrators can be very accurate due to active control of errors, but they can also
 * lead to very long computational times for the integration when applied to stiff systems.
 * This method is the original Runge-Kutta-Fehlberg scheme, however, today, the Dormand-Prince variant is more popular.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_FEHLBERG45_INTEGRATOR_SYS_HPP
#define REAK_FEHLBERG45_INTEGRATOR_SYS_HPP

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/control/systems/state_space_sys_concept.hpp>
#include <ReaK/topologies/spaces/metric_space_concept.hpp>
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>
#include <ReaK/topologies/interpolation/spatial_trajectory_concept.hpp>

#include <ReaK/math/integrators/integration_exceptions.hpp>
#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

namespace ReaK {

namespace ctrl {


namespace detail {


//----------CFEHLBERG45-----------------------------------------------------

/*
     0     |
     1/4   |  1/4
     3/8   |  3/32         9/32
     12/13 |  1932/2197   -7200/2197    7296/2197
     1     |  439/216     -8            3680/513    -845/4104
     1/2   | -8/27         2           -3544/2565    1859/4104   -11/40
     _____________________________________________________________________________________
     O(h^5)|  25/216       0            1408/2565    2197/4104   -1/5          0
     O(h^6)|  16/135       0            6656/12825   28561/56430 -9/50         2/55
*/

template < typename StateSpace, typename StateSpaceSystem, typename InputTrajectory >
void fehlberg45_integrate_impl( const StateSpace& space, const StateSpaceSystem& sys,
                                const typename pp::topology_traits< StateSpace >::point_type& start_point,
                                typename pp::topology_traits< StateSpace >::point_type& end_point,
                                const InputTrajectory& u_traj, double start_time, double end_time, double time_step,
                                double tolerance, double min_step, double max_step ) {
  using std::fabs;
  using std::pow;
  using ReaK::to_vect;
  typedef typename pp::topology_traits< StateSpace >::point_type PointType;
  typedef typename pp::topology_traits< StateSpace >::point_difference_type PointDiffType;

  if( ( time_step == 0.0 ) || ( ( time_step > 0.0 ) && ( start_time > end_time ) )
      || ( ( time_step < 0.0 ) && ( end_time > start_time ) ) || ( tolerance <= 0.0 ) || ( min_step > max_step ) )
    throw impossible_integration( start_time, end_time, time_step );

  typedef typename pp::spatial_trajectory_traits< InputTrajectory >::const_waypoint_descriptor InputWaypoint;
  typedef typename pp::spatial_trajectory_traits< InputTrajectory >::point_type InputType;
  std::pair< InputWaypoint, InputType > u_wp = u_traj.get_waypoint_at_time( start_time );

  PointDiffType dp = sys.get_state_derivative( space, start_point, u_wp.second.pt, start_time );
  double t = start_time;
  end_point = start_point;

  while( ( ( time_step > 0.0 ) && ( t < end_time ) ) || ( ( time_step < 0.0 ) && ( t > end_time ) ) ) {

    PointType prevY = end_point;
    PointDiffType k1 = time_step * dp;
    end_point = space.adjust( end_point, 0.25 * k1 );

    t += 0.25 * time_step;
    u_wp = u_traj.move_time_diff_from( u_wp, 0.25 * time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    PointDiffType k2 = time_step * dp;
    end_point = space.adjust( prevY, ( 3.0 / 32.0 ) * k1 + ( 9.0 / 32.0 ) * k2 );

    t += time_step * 0.125;
    u_wp = u_traj.move_time_diff_from( u_wp, 0.125 * time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    PointDiffType k3 = time_step * dp;
    end_point = space.adjust( prevY, ( 1932.0 / 2197.0 ) * k1 - ( 7200.0 / 2197.0 ) * k2 + ( 7296.0 / 2197.0 ) * k3 );

    t += 57.0 * time_step / 104.0;
    u_wp = u_traj.move_time_diff_from( u_wp, ( 57.0 / 104.0 ) * time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    PointDiffType k4 = time_step * dp;
    end_point
      = space.adjust( prevY, ( 439.0 / 216.0 ) * k1 - 8.0 * k2 + ( 3680.0 / 513.0 ) * k3 - ( 845.0 / 4104.0 ) * k4 );

    t += time_step / 13.0;
    u_wp = u_traj.move_time_diff_from( u_wp, ( 1.0 / 13.0 ) * time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    PointDiffType k5 = time_step * dp;
    end_point = space.adjust( prevY, 2.0 * k2 - ( 8.0 / 27.0 ) * k1 - ( 3544.0 / 2565.0 ) * k3
                                     + ( 1859.0 / 4104.0 ) * k4 - ( 11.0 / 40.0 ) * k5 );

    t -= time_step * 0.5;
    u_wp = u_traj.move_time_diff_from( u_wp, -0.5 * time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    PointDiffType k6 = time_step * dp;

    vect_n< double > err_vect = to_vect< double >( 360.0 * k1 - ( 128.0 / 4275.0 ) * k3 - ( 2197.0 / 75240.0 ) * k4
                                                   + ( 1.0 / 50.0 ) * k5 + ( 2.0 / 55.0 ) * k6 );
    double Rmax = 0.0;
    std::size_t worst_DOF = 0;
    for( std::size_t i = 0; i < err_vect.size(); ++i ) {
      double R = fabs( err_vect[i] / time_step );
      if( R > Rmax ) {
        Rmax = R;
        worst_DOF = i;
      };
    };

    if( Rmax > tolerance ) {
      if( fabs( time_step ) <= min_step )
        throw untolerable_integration( tolerance, Rmax, worst_DOF, time_step, t );

      t -= time_step * 0.5;
      u_wp = u_traj.move_time_diff_from( u_wp, -0.5 * time_step );
      end_point = prevY;
      dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
      double R = 0.84 * pow( tolerance / Rmax, 0.25 );
      if( R < 0.1 )
        time_step *= 0.1;
      else
        time_step *= R;
    } else {
      t += time_step * 0.5;
      u_wp = u_traj.move_time_diff_from( u_wp, 0.5 * time_step );
      end_point
        = space.adjust( prevY, ( 25.0 / 216.0 ) * k1 + ( 1408.0 / 2565.0 ) * k3 + ( 2197.0 / 4104.0 ) * k4 - 0.2 * k5 );
      dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );

      double R = 0.84 * pow( tolerance / Rmax, 0.25 );
      if( R >= 4.0 )
        time_step *= 4.0;
      else if( R > 1.0 )
        time_step *= R;
      // if(((time_step > 0.0) && (t + time_step > aEndTime)) || ((time_step < 0.0) && (t + time_step < aEndTime)))
      // time_step = aEndTime - t;
    };

    if( fabs( time_step ) < min_step )
      time_step *= fabs( min_step / time_step );
    if( fabs( time_step ) > max_step )
      time_step *= fabs( max_step / time_step );
  };
};
};

/**
 * This class is a factory for integrators that use the 4-5th Order Fehlberg Method. This
 * method is of the Runge-Kutta-Fehlberg type which are single-step variable-step algorithm for numerical
 * integration in which the error is estimated by the difference between two integration rules (order 4 and 5)
 * and the time-step is adapted on-the-fly to control the norm of that error estimate.
 * These variable-step integrators can be very accurate due to active control of errors, but they can also
 * lead to very long computational times for the integration when applied to stiff systems.
 * This method is the original Runge-Kutta-Fehlberg scheme, however, today, the Dormand-Prince variant is more popular.
 * \tparam TemporalSpace The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations), see
 * SSSystemConcept.
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times, see
 * pp::SpatialTrajectoryConcept.
 */
template < typename TemporalSpace, typename StateSpaceSystem, typename InputTrajectory >
class fehlberg45_integrator_factory : public named_object {
public:
  typedef fehlberg45_integrator_factory< TemporalSpace, StateSpaceSystem, InputTrajectory > self;
  typedef TemporalSpace topology;
  typedef typename pp::topology_traits< TemporalSpace >::point_type point_type;

  typedef typename pp::temporal_space_traits< TemporalSpace >::space_topology space_topology;
  typedef typename pp::temporal_space_traits< TemporalSpace >::time_topology time_topology;

  BOOST_CONCEPT_ASSERT( ( pp::TemporalSpaceConcept< TemporalSpace > ) );
  BOOST_CONCEPT_ASSERT( ( SSSystemConcept< StateSpaceSystem, space_topology > ) );

private:
  shared_ptr< const TemporalSpace > m_t_space;
  shared_ptr< const StateSpaceSystem > m_sys;
  shared_ptr< const InputTrajectory > m_input_traj;
  double m_time_step;
  double m_tolerance;
  double m_min_step;
  double m_max_step;

public:
  /**
   * This class represents an integration task, from a starting temporal point.
   */
  class extrapolator_type {
  private:
    const self* m_parent;
    const point_type* m_start_point;

  public:
    extrapolator_type( const self* aParent, const point_type* aStartPoint )
        : m_parent( aParent ), m_start_point( aStartPoint ){};

    /**
     * Sets the starting point (by a raw-pointer) of the interpolation.
     * \param aStartPoint A raw-pointer to the starting point of the interpolation.
     */
    void set_start_point( const point_type* aStartPoint ) { m_start_point = aStartPoint; };

    /**
     * Returns the pointer to the starting point of the interpolation.
     * \return The pointer to the starting point of the interpolation.
     */
    const point_type* get_start_point() const { return m_start_point; };

    /**
     * Returns the point integrated up to the given time (within the time-step precision).
     * \param end_time The end of the integration period.
     * \return The point resulting from integration from the starting point up to the given end-time.
     */
    point_type get_point_at_time( double end_time ) const {
      point_type end_point;
      end_point.time = end_time;
      detail::fehlberg45_integrate_impl( m_parent->m_t_space->get_space_topology(), *( m_parent->m_sys ),
                                         m_start_point->pt, end_point.pt, *( m_parent->m_input_traj ),
                                         m_start_point->time, end_point.time, m_parent->m_time_step,
                                         m_parent->m_tolerance, m_parent->m_min_step, m_parent->m_max_step );
      return end_point;
    };
  };

  /**
   * Parametrized Constructor.
   * \param aName The name of the integrator factory object.
   * \param aTSpace A pointer to the temporal space to use.
   * \param aSystem A pointer to the state-space system to be integrated.
   * \param aInputTraj A pointer to the trajectory object which can deliver input vectors at any given time point.
   * \param aTimeStep The integration time-step to use.
   * \param aTolerance The tolerance on the norm of the estimated integration error.
   * \param aMinStep The minimum integration time-step to use (will throw an intolerable_integration exception if
   * reached).
   * \param aMaxStep The maximum integration time-step to use.
   */
  fehlberg45_integrator_factory(
    const std::string& aName = "",
    const shared_ptr< const TemporalSpace >& aTSpace = shared_ptr< const TemporalSpace >(),
    const shared_ptr< const StateSpaceSystem >& aSystem = shared_ptr< const StateSpaceSystem >(),
    const shared_ptr< const InputTrajectory >& aInputTraj = shared_ptr< const InputTrajectory >(),
    double aTimeStep = 1e-3, double aTolerance = 1e-6, double aMinStep = 1e-8, double aMaxStep = 1e-1 )
      : named_object(), m_t_space( aTSpace ), m_sys( aSystem ), m_input_traj( aInputTraj ), m_time_step( aTimeStep ),
        m_tolerance( aTolerance ), m_min_step( aMinStep ), m_max_step( aMaxStep ){};

  /**
   * Sets the pointer to the temporal space used by this integrator factory.
   * \param aTSpace A pointer to the temporal space used by this integrator factory.
   */
  void set_temporal_space( const shared_ptr< const TemporalSpace >& aTSpace ) { m_t_space = aTSpace; };
  /**
   * Returns a pointer to the temporal space used by this integrator factory.
   * \return A pointer to the temporal space used by this integrator factory.
   */
  const shared_ptr< const TemporalSpace >& get_temporal_space() const { return m_t_space; };

  /**
   * Sets the pointer to the state-space system used by this integrator factory.
   * \param aSystem A pointer to the state-space system used by this integrator factory.
   */
  void set_system( const shared_ptr< const StateSpaceSystem >& aSystem ) { m_sys = aSystem; };
  /**
   * Returns a pointer to the state-space system used by this integrator factory.
   * \return A pointer to the state-space system used by this integrator factory.
   */
  const shared_ptr< const StateSpaceSystem >& get_system() const { return m_sys; };

  /**
   * Sets the pointer to the input-vector trajectory used by this integrator factory.
   * \param aInputTraj A pointer to the input-vector trajectory used by this integrator factory.
   */
  void set_input_trajectory( const shared_ptr< const InputTrajectory >& aInputTraj ) { m_input_traj = aInputTraj; };
  /**
   * Returns a pointer to the input-vector trajectory used by this integrator factory.
   * \return A pointer to the input-vector trajectory used by this integrator factory.
   */
  const shared_ptr< const InputTrajectory >& get_input_trajectory() const { return m_input_traj; };

  /**
   * Creates an "extrapolator" from a given starting point (by pointer). This creates
   * an IVP integrator (which is a special kind of extrapolator).
   * \param aStartPoint A raw-pointer to a starting point (raw pointers are used as starting points are usually stored
   * in a waypoint trajectory).
   * \return An IVP integrator from the given starting point.
   */
  extrapolator_type create_extrapolator( const point_type* aStartPoint ) const {
    return extrapolator_type( this, aStartPoint );
  };


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    ReaK::named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( m_t_space ) & RK_SERIAL_SAVE_WITH_NAME( m_sys )
      & RK_SERIAL_SAVE_WITH_NAME( m_input_traj ) & RK_SERIAL_SAVE_WITH_NAME( m_time_step )
      & RK_SERIAL_SAVE_WITH_NAME( m_tolerance ) & RK_SERIAL_SAVE_WITH_NAME( m_min_step )
      & RK_SERIAL_SAVE_WITH_NAME( m_max_step );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    ReaK::named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( m_t_space ) & RK_SERIAL_LOAD_WITH_NAME( m_sys )
      & RK_SERIAL_LOAD_WITH_NAME( m_input_traj ) & RK_SERIAL_LOAD_WITH_NAME( m_time_step )
      & RK_SERIAL_LOAD_WITH_NAME( m_tolerance ) & RK_SERIAL_LOAD_WITH_NAME( m_min_step )
      & RK_SERIAL_LOAD_WITH_NAME( m_max_step );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2431005, 1, "fehlberg45_integrator_factory", named_object )
};
};
};

#endif
