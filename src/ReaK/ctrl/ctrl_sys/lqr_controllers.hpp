/**
 * \file lqr_controllers.hpp
 * 
 * This library provides a class template which can be used to create LQR 
 * control systems, as used in ReaK::ctrl. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2012
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

#ifndef REAK_LQR_CONTROLLERS_HPP
#define REAK_LQR_CONTROLLERS_HPP

#include "linear_ss_system_concept.hpp"
#include "discrete_linear_sss_concept.hpp"

#include "base/named_object.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/vect_alg.hpp"

#include "lin_alg/mat_are_solver.hpp"

namespace ReaK {

namespace ctrl {

/**
 * This class template can be used to create a simple discrete-time infinite-horizon LQR 
 * control system, as used in ReaK::ctrl. A discrete-time infinite-horizon LQR state-space 
 * control system is described by two matrices (Q,R) which make the scaling 
 * on the state-error and the input-vector, respectively.
 * \tparam StateSpaceSystem The system under control.
 * \tparam StateSpaceType The state-space type on which the state-vectors of the controlled system lie.
 * \tparam TrajectoryType The trajectory to be followed.
 */
template <typename StateSpaceSystem, typename StateSpaceType, typename TrajectoryType>
class dt_ih_lqr_controller : public named_object {
  public:
    BOOST_CONCEPT_ASSERT((DiscreteLinearSSSConcept<StateSpaceSystem, StateSpaceType, DiscreteLinearizedSystemType>));
    BOOST_CONCEPT_ASSERT((pp::SpatialTrajectoryConcept<TrajectoryType, pp::temporal_space<StateSpaceType, pp::time_topology> >));
    
    typedef double point_type;
    typedef double point_difference_type;
    typedef double point_derivative_type;
    
    typedef typename discrete_sss_traits<StateSpaceSystem>::point_type input_type;
    typedef typename discrete_sss_traits<StateSpaceSystem>::input_type output_type;
    
    typedef typename mat_traits<typename discrete_linear_sss_traits<StateSpaceSystem>::matrixA_type>::value_type value_type;
    
    typedef mat< value_type, mat_structure::identity > matrixA_type;
    typedef mat< value_type, mat_structure::nil > matrixB_type;
    typedef mat< value_type, mat_structure::nil > matrixC_type;
    typedef mat< value_type, mat_structure::rectangular > matrixD_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = discrete_sss_traits<StateSpaceSystem>::dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = discrete_sss_traits<StateSpaceSystem>::input_dimensions);
    
  private:
    shared_ptr< StateSpaceSystem > m_sys;
    shared_ptr< StateSpaceType > m_state_space;
    shared_ptr< TrajectoryType > m_traj;
    
    
    mat<value_type,mat_structure::square> Q;
    mat<value_type,mat_structure::square> R;
    
  public:
        
    /**
     * Parametrized and default constructor.
     * \param aName The name of the controller.
     * \param aQ The positive-definite matrix to scale the state-errors.
     * \param aR The positive-semi-definite matrix to scale the input-vectors.
     * \param aSys The state-space system under control.
     * \param aSpace The state-space topology.
     * \param aTraj The trajectory to be followed.
     */
    dt_ih_lqr_controller(
      const std::string& aName = "",
      const mat<value_type, mat_structure::square>& aQ = mat<value_type,mat_structure::square>(mat<value_type,mat_structure::identity>(input_dimensions)),
      const mat<value_type, mat_structure::square>& aR = mat<value_type,mat_structure::square>(mat<value_type,mat_structure::identity>(output_dimensions)),
      const shared_ptr< StateSpaceSystem >& aSys = shared_ptr< StateSpaceSystem >(),
      const shared_ptr< StateSpaceType >& aSpace = shared_ptr< StateSpaceType >(),
      const shared_ptr< TrajectoryType >& aTraj = shared_ptr< TrajectoryType >()) :
      m_sys(aSys), m_state_space(aSpace), m_traj(aTraj), 
      Q(aQ), 
      R(aR) {
      setName(aName);
    };
    
    /**
     * Sets the state-space system under control.
     * \param aSys The state-space system under control.
     */
    void set_system(const shared_ptr< StateSpaceSystem >& aSys) {
      m_sys = aSys;
    };
    
    /**
     * Sets the state-space topology.
     * \param aSpace The new state-space topology.
     */
    void set_space(const shared_ptr< StateSpaceType >& aSpace) {
      m_state_space = aSpace;
    };
    
    /**
     * Sets the trajectory to be followed.
     * \param aTraj The trajectory to be followed.
     */
    void set_trajectory(const shared_ptr< TrajectoryType >& aTraj) {
      m_traj = aTraj;
    };
    
    /**
     * Returns the state-space system under control.
     * \return The state-space system under control.
     */
    const shared_ptr< StateSpaceSystem >& get_system() const { return m_sys; };
    /**
     * Returns the state-space topology.
     * \return The state-space topology.
     */
    const shared_ptr< StateSpaceType >& get_space() const { return m_state_space; };
    /**
     * Returns the trajectory to be followed.
     * \return The trajectory to be followed.
     */
    const shared_ptr< TrajectoryType >& get_trajectory() const { return m_traj; };
    
    /**
     * Returns the state-vector's dimension.
     * \return The state-vector's dimension.
     */
    size_type get_state_count() const { return 0; };
    /**
     * Returns the input-vector's dimension.
     * \return The input-vector's dimension.
     */
    size_type get_input_count() const { return Q.get_row_count(); };
    /**
     * Returns the output-vector's dimension.
     * \return The output-vector's dimension.
     */
    size_type get_output_count() const { return R.get_row_count(); };
    
    /**
     * Returns output of the system given the current state, input and time.
     * \param u The current input.
     * \param t The current time.
     * \return The current output.
     */
    output_type get_output(const input_type& u, const double& t = 0) const {
      typedef typename discrete_linear_sss_traits<StateSpaceSystem>::matrixA_type Atype;
      typedef typename discrete_linear_sss_traits<StateSpaceSystem>::matrixB_type Btype;
      
      Atype A;
      Btype B;
      m_sys->get_state_transition_blocks(A, B, *m_state_space, t, t, u, u, from_vect<output_type>(vect_n<value_type>(get_output_count(),value_type(0.0))), from_vect<output_type>(vect_n<value_type>(get_output_count(),value_type(0.0))));
      
      mat< value_type, mat_structure::rectangular > K(get_output_count(), get_input_count());
      solve_IHDT_LQR(A, B, Q, R, K, value_type(1e-4));
      
      return from_vect<output_type>( K * to_vect<value_type>( m_state_space->difference(m_traj->get_point_at_time(t).pt, u) ) );
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(m_sys)
         & RK_SERIAL_SAVE_WITH_NAME(m_state_space)
         & RK_SERIAL_SAVE_WITH_NAME(m_traj)
         & RK_SERIAL_SAVE_WITH_NAME(Q)
         & RK_SERIAL_SAVE_WITH_NAME(R);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(m_sys)
         & RK_SERIAL_LOAD_WITH_NAME(m_state_space)
         & RK_SERIAL_LOAD_WITH_NAME(m_traj)
         & RK_SERIAL_LOAD_WITH_NAME(Q)
         & RK_SERIAL_LOAD_WITH_NAME(R);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300012,1,"dt_ih_lqr_controller",named_object)
    
};




/**
 * This class template can be used to create a simple continuous-time infinite-horizon LQR 
 * control system, as used in ReaK::ctrl. A continuous-time infinite-horizon LQR state-space 
 * control system is described by two matrices (Q,R) which make the scaling 
 * on the state-error and the input-vector, respectively.
 * \tparam StateSpaceSystem The system under control.
 * \tparam StateSpaceType The state-space type on which the state-vectors of the controlled system lie.
 * \tparam TrajectoryType The trajectory to be followed.
 */
template <typename StateSpaceSystem, typename StateSpaceType, typename TrajectoryType>
class ct_ih_lqr_controller : public named_object {
  public:
    BOOST_CONCEPT_ASSERT((LinearSSSystemConcept<StateSpaceSystem, StateSpaceType, LinearizedSystemType>));
    BOOST_CONCEPT_ASSERT((pp::SpatialTrajectoryConcept<TrajectoryType, pp::temporal_space<StateSpaceType, pp::time_topology> >));
    
    typedef double point_type;
    typedef double point_difference_type;
    typedef double point_derivative_type;
    
    typedef typename ss_system_traits<StateSpaceSystem>::point_type input_type;
    typedef typename ss_system_traits<StateSpaceSystem>::input_type output_type;
    
    typedef typename mat_traits<typename linear_ss_system_traits<StateSpaceSystem>::matrixA_type>::value_type value_type;
    
    typedef mat< value_type, mat_structure::identity > matrixA_type;
    typedef mat< value_type, mat_structure::nil > matrixB_type;
    typedef mat< value_type, mat_structure::nil > matrixC_type;
    typedef mat< value_type, mat_structure::rectangular > matrixD_type;
    
    typedef double time_type;
    typedef double time_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = ss_system_traits<StateSpaceSystem>::dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = ss_system_traits<StateSpaceSystem>::input_dimensions);
    
  private:
    shared_ptr< StateSpaceSystem > m_sys;
    shared_ptr< StateSpaceType > m_state_space;
    shared_ptr< TrajectoryType > m_traj;
    
    
    mat<value_type,mat_structure::square> Q;
    mat<value_type,mat_structure::square> R;
    
  public:
        
    /**
     * Parametrized and default constructor.
     * \param aName The name of the controller.
     * \param aQ The positive-definite matrix to scale the state-errors.
     * \param aR The positive-semi-definite matrix to scale the input-vectors.
     * \param aSys The state-space system under control.
     * \param aSpace The state-space topology.
     * \param aTraj The trajectory to be followed.
     */
    ct_ih_lqr_controller(
      const std::string& aName = "",
      const mat<value_type, mat_structure::square>& aQ = mat<value_type,mat_structure::square>(mat<value_type,mat_structure::identity>(input_dimensions)),
      const mat<value_type, mat_structure::square>& aR = mat<value_type,mat_structure::square>(mat<value_type,mat_structure::identity>(output_dimensions)),
      const shared_ptr< StateSpaceSystem >& aSys = shared_ptr< StateSpaceSystem >(),
      const shared_ptr< StateSpaceType >& aSpace = shared_ptr< StateSpaceType >(),
      const shared_ptr< TrajectoryType >& aTraj = shared_ptr< TrajectoryType >()) :
      m_sys(aSys), m_state_space(aSpace), m_traj(aTraj), 
      Q(aQ), 
      R(aR) {
      setName(aName);
    };
    
    /**
     * Sets the state-space system under control.
     * \param aSys The state-space system under control.
     */
    void set_system(const shared_ptr< StateSpaceSystem >& aSys) {
      m_sys = aSys;
    };
    
    /**
     * Sets the state-space topology.
     * \param aSpace The new state-space topology.
     */
    void set_space(const shared_ptr< StateSpaceType >& aSpace) {
      m_state_space = aSpace;
    };
    
    /**
     * Sets the trajectory to be followed.
     * \param aTraj The trajectory to be followed.
     */
    void set_trajectory(const shared_ptr< TrajectoryType >& aTraj) {
      m_traj = aTraj;
    };
    
    /**
     * Returns the state-space system under control.
     * \return The state-space system under control.
     */
    const shared_ptr< StateSpaceSystem >& get_system() const { return m_sys; };
    /**
     * Returns the state-space topology.
     * \return The state-space topology.
     */
    const shared_ptr< StateSpaceType >& get_space() const { return m_state_space; };
    /**
     * Returns the trajectory to be followed.
     * \return The trajectory to be followed.
     */
    const shared_ptr< TrajectoryType >& get_trajectory() const { return m_traj; };
    
    /**
     * Returns the state-vector's dimension.
     * \return The state-vector's dimension.
     */
    size_type get_state_count() const { return 0; };
    /**
     * Returns the input-vector's dimension.
     * \return The input-vector's dimension.
     */
    size_type get_input_count() const { return Q.get_row_count(); };
    /**
     * Returns the output-vector's dimension.
     * \return The output-vector's dimension.
     */
    size_type get_output_count() const { return R.get_row_count(); };
    
    /**
     * Returns output of the system given the current state, input and time.
     * \param u The current input.
     * \param t The current time.
     * \return The current output.
     */
    output_type get_output(const input_type& u, const double& t = 0) const {
      typedef typename linear_ss_system_traits<StateSpaceSystem>::matrixA_type Atype;
      typedef typename linear_ss_system_traits<StateSpaceSystem>::matrixB_type Btype;
      typedef typename linear_ss_system_traits<StateSpaceSystem>::matrixC_type Ctype;
      typedef typename linear_ss_system_traits<StateSpaceSystem>::matrixD_type Dtype;
      
      Atype A;
      Btype B;
      Ctype C;
      Dtype D;
      m_sys->get_linear_blocks(A, B, C, D, *m_state_space, t, u, from_vect<output_type>(vect_n<value_type>(get_output_count(),value_type(0.0))));
      
      mat< value_type, mat_structure::rectangular > K(get_output_count(), get_input_count());
      solve_IHCT_LQR(A, B, Q, R, K, value_type(1e-4));
      
      return from_vect<output_type>( K * to_vect<value_type>( m_state_space->difference(m_traj->get_point_at_time(t).pt, u) ) );
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(m_sys)
         & RK_SERIAL_SAVE_WITH_NAME(m_state_space)
         & RK_SERIAL_SAVE_WITH_NAME(m_traj)
         & RK_SERIAL_SAVE_WITH_NAME(Q)
         & RK_SERIAL_SAVE_WITH_NAME(R);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(m_sys)
         & RK_SERIAL_LOAD_WITH_NAME(m_state_space)
         & RK_SERIAL_LOAD_WITH_NAME(m_traj)
         & RK_SERIAL_LOAD_WITH_NAME(Q)
         & RK_SERIAL_LOAD_WITH_NAME(R);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300013,1,"ct_ih_lqr_controller",named_object)
    
};

  


};

};

#endif










