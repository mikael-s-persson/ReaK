/**
 * \file adams_BM3_integrator_sys.hpp
 *
 * This library implements an integrator that uses the 3rd Order Adams-Bashforth-Moulton Method. Adams
 * methods are a type of multi-step predictor-corrector algorithm for numerical integration in which
 * the prediction step is first performed and then a fixed number of correction steps.
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

#ifndef REAK_ADAMS_BM3_INTEGRATOR_SYS_HPP
#define REAK_ADAMS_BM3_INTEGRATOR_SYS_HPP

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


//----------adamsBM_integrator------------------------------------------------------

/*
  Bashforth Coefficients
  p | k | j->    |  1     2     3     4     5     6      |  C
  -------------------------------------------------------------------
  1 | 1 | Bj     |  1                                    | 1/2
  2 | 2 | 2Bj    |  3    -1                              | 5/12
  3 | 3 | 12Bj   |  23   -16    5                        | 3/8
  4 | 4 | 24Bj   |  55   -59    37   -9                  | 251/720
  5 | 5 | 720Bj  |  1901 -2774  2616 -1274  251          | 95/288
  6 | 6 | 1440Bj |  4277 -7923  9982 -7298  2877 -475    | 19087/60480

  Yn = Yn-1 + h*sum(j=1->k, Bj*Fn-j);

  Moulton Coefficients
  p | k | j->    |  0     1     2     3     4     5      |  C
  -------------------------------------------------------------------
  1 | 1 | Bj     |  1                                    | -1/2
  2 | 1 | 2Bj    |  1     1                              | -1/12
  3 | 2 | 12Bj   |  5     8    -1                        | -1/24
  4 | 3 | 24Bj   |  9     19   -5     1                  | -19/720
  5 | 4 | 720Bj  |  251   646  -264   106  -19           | -3/160
  6 | 5 | 1440Bj |  475   1427 -798   482  -173   27     | -863/60480

  Yn = Yn-1 + h*sum(j=0->k, Bj*Fn-j); //note: implicit, hence Yn or Fn is not known but predicted using Bashforth

  Predictor-Corrector Steps:
  Predictor : compute (Yn)0 using Bashforth
  Evaluation : evaluate (Fn)0 using (Yn)0
  Corrector : compute (Yn)1 using Moulton  -> iterate again Eval/Corr by a fixed number of times "Corrections"
  Evaluation : evaluate (Fn)1 using (Yn)1 for the value of Fn-1 of the next time step

*/

template < typename StateSpace, typename StateSpaceSystem, typename InputTrajectory >
void adams_BM3_integrate_impl( const StateSpace& space, const StateSpaceSystem& sys,
                               const typename pp::topology_traits< StateSpace >::point_type& start_point,
                               typename pp::topology_traits< StateSpace >::point_type& end_point,
                               const InputTrajectory& u_traj, double start_time, double end_time, double time_step,
                               std::size_t correction_count ) {
  if( ( time_step == 0.0 ) || ( ( time_step > 0.0 ) && ( t > end_time ) )
      || ( ( time_step < 0.0 ) && ( t < end_time ) ) )
    throw impossible_integration( t, end_time, time_step );

  if( correction_count == 0 )
    correction_count = 1;

  using ReaK::to_vect;
  typedef typename pp::topology_traits< StateSpace >::point_type PointType;
  typedef typename pp::topology_traits< StateSpace >::point_difference_type PointDiffType;

  typedef typename pp::spatial_trajectory_traits< InputTrajectory >::const_waypoint_descriptor InputWaypoint;
  typedef typename pp::spatial_trajectory_traits< InputTrajectory >::point_type InputType;
  std::pair< InputWaypoint, InputType > u_wp = u_traj.get_waypoint_at_time( start_time );


  double back_time = start_time;
  end_point = start_point;
  PointDiffType dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, back_time );

  // Order 1
  back_time -= time_step;
  PointType prevY1 = space.adjust( end_point, -time_step * dp );
  u_wp = u_traj.move_time_diff_from( u_wp, -time_step );
  PointDiffType prevF1 = sys.get_state_derivative( space, prevY1, u_wp.second.pt, back_time );

  for( std::size_t j = 0; j < correction_count; ++j ) {
    prevY1 = space.adjust( end_point, -time_step * prevF1 );
    prevF1 = sys.get_state_derivative( space, prevY1, u_wp.second.pt, back_time );
  };

  // Order 2
  back_time -= time_step;
  PointType prevY2 = space.adjust( prevY1, ( 0.5 * time_step ) * dp - ( 1.5 * time_step ) * prevF1 );
  u_wp = u_traj.move_time_diff_from( u_wp, -time_step );
  PointDiffType prevF2 = sys.get_state_derivative( space, prevY2, u_wp.second.pt, back_time );

  for( std::size_t j = 0; j < correction_count; ++j ) {
    prevY2 = space.adjust( prevY1, ( -0.5 * time_step ) * ( prevF1 + prevF2 ) );
    prevF2 = sys.get_state_derivative( space, prevY2, u_wp.second.pt, back_time );
  };

  // Order 3
  back_time -= time_step;
  PointType prevY3 = space.adjust( prevY2, ( 16.0 * time_step / 12.0 ) * prevF1 - ( 23.0 * time_step / 12.0 ) * prevF2
                                           - ( 5.0 * time_step / 12.0 ) * dp );
  u_wp = u_traj.move_time_diff_from( u_wp, -time_step );
  PointDiffType prevF3 = sys.get_state_derivative( space, prevY3, u_wp.second.pt, back_time );

  for( std::size_t j = 0; j < correction_count; ++j ) {
    prevY3 = space.adjust( prevY2, ( time_step / 12.0 ) * prevF1 - ( 8.0 * time_step / 12.0 ) * prevF2
                                   - ( 5.0 * time_step / 12.0 ) * prevF3 );
    prevF3 = sys.get_state_derivative( space, prevY3, u_wp.second.pt, back_time );
  };


  double t = start_time;
  u_wp = u_traj.get_waypoint_at_time( t );

  while( ( ( time_step > 0.0 ) && ( t < end_time ) ) || ( ( time_step < 0.0 ) && ( t > end_time ) ) ) {

    prevY1 = end_point;
    prevF3 = prevF2;
    prevF2 = prevF1;
    prevF1 = dp;
    end_point = space.adjust( end_point, ( 23.0 * time_step / 12.0 ) * prevF1 - ( 16.0 * time_step / 12.0 ) * prevF2
                                         + ( 5.0 * time_step / 12.0 ) * prevF3 );
    t += time_step;
    u_wp = u_traj.move_time_diff_from( u_wp, time_step );
    dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );

    for( std::size_t j = 0; j < correction_count; ++j ) {
      end_point = space.adjust( prevY1, ( 5.0 * time_step / 12.0 ) * dp + ( 8.0 * time_step / 12.0 ) * prevF1
                                        - ( time_step / 12.0 ) * prevF2 );
      dp = sys.get_state_derivative( space, end_point, u_wp.second.pt, t );
    };
  };
};
};


/**
 * This class is a factory for integrators that use the 3rd Order Adams-Bashforth-Moulton Method. Adams
 * methods are a type of multi-step predictor-corrector algorithm for numerical integration in which
 * the prediction step is first performed and then a fixed number of correction steps.
 * \tparam TemporalSpace The temporal space type (space-time topology) on which the computed trajectories lay.
 * \tparam StateSpaceSystem The continuous-time state-space system type to integrate (governing equations), see
 * SSSystemConcept.
 * \tparam InputTrajectory The trajectory type which can deliver input vectors at given times, see
 * pp::SpatialTrajectoryConcept.
 */
template < typename TemporalSpace, typename StateSpaceSystem, typename InputTrajectory >
class adams_BM3_integrator_factory : public named_object {
public:
  typedef adams_BM3_integrator_factory< TemporalSpace, StateSpaceSystem, InputTrajectory > self;
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
  std::size_t m_correction_count;

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
      detail::adams_BM3_integrate_impl( m_parent->m_t_space->get_space_topology(), *( m_parent->m_sys ),
                                        m_start_point->pt, end_point.pt, *( m_parent->m_input_traj ),
                                        m_start_point->time, end_point.time, m_parent->m_time_step,
                                        m_parent->m_correction_count );
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
   * \param aCorrectionCount The number of correction iterations (should be between 1 and 5 or so).
   */
  adams_BM3_integrator_factory( const std::string& aName = "",
                                const shared_ptr< const TemporalSpace >& aTSpace = shared_ptr< const TemporalSpace >(),
                                const shared_ptr< const StateSpaceSystem >& aSystem
                                = shared_ptr< const StateSpaceSystem >(),
                                const shared_ptr< const InputTrajectory >& aInputTraj
                                = shared_ptr< const InputTrajectory >(),
                                double aTimeStep = 1e-3, std::size_t aCorrectionCount = 3 )
      : named_object(), m_t_space( aTSpace ), m_sys( aSystem ), m_input_traj( aInputTraj ), m_time_step( aTimeStep ),
        m_correction_count( aCorrectionCount ){};

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
      & RK_SERIAL_SAVE_WITH_NAME( m_correction_count );
  };

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    ReaK::named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( m_t_space ) & RK_SERIAL_LOAD_WITH_NAME( m_sys )
      & RK_SERIAL_LOAD_WITH_NAME( m_input_traj ) & RK_SERIAL_LOAD_WITH_NAME( m_time_step )
      & RK_SERIAL_LOAD_WITH_NAME( m_correction_count );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( self, 0xC2431007, 1, "adams_BM3_integrator_factory", named_object )
};
};
};

#endif
