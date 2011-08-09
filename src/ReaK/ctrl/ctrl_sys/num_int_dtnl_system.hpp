/**
 * \file num_int_dtnl_system.hpp
 * 
 * This library defines a class template which can be used to integrate, numerically, a continuous-time
 * state-space system (SSSystemConcept) to produce a discrete-time state-space system (DiscreteSSSConcept).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef NUM_INT_DTNL_SYSTEM_HPP
#define NUM_INT_DTNL_SYSTEM_HPP

#include "discrete_sss_concept.hpp"
#include "state_space_sys_concept.hpp"
#include "integrators/integrator.hpp"
#include "integrators/fixed_step_integrators.hpp"

#include "base/named_object.hpp"

namespace ReaK {

namespace ctrl { 

/**
 * This class template can be used to integrate, numerically, a continuous-time
 * state-space system (SSSystemConcept) to produce a discrete-time state-space system (DiscreteSSSConcept).
 * \tparam CTSystem The continuous-time state-space system type, should model SSSystemConcept.
 * \tparam NumIntegrator The numerical integrator type to be used.
 */
template <typename CTSystem, typename NumIntegrator = euler_integrator<typename CTSystem::value_type> >
class num_int_dtnl_sys : public named_object, public state_rate_function<typename CTSystem::value_type> {
  public:
    typedef num_int_dtnl_sys<CTSystem> self;
    typedef typename CTSystem::value_type value_type;
    typedef typename CTSystem::size_type size_type;
    
    typedef typename ss_system_traits<CTSystem>::point_type point_type;
    typedef typename ss_system_traits<CTSystem>::point_difference_type point_difference_type;
    typedef typename ss_system_traits<CTSystem>::point_derivative_type point_derivative_type;
  
    typedef typename ss_system_traits<CTSystem>::time_type time_type;
    typedef typename ss_system_traits<CTSystem>::time_difference_type time_difference_type;
  
    typedef typename ss_system_traits<CTSystem>::input_type input_type;
    typedef typename ss_system_traits<CTSystem>::output_type output_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = ss_system_traits<CTSystem>::dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = ss_system_traits<CTSystem>::input_dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = ss_system_traits<CTSystem>::output_dimensions);
  
  private:
    
    CTSystem sys;
    NumIntegrator integ;
    time_difference_type dt;
    
    input_type current_u;
        
  public:
    
    /**
     * Default constructor.
     */
    num_int_dtnl_sys(const std::string& aName = "") : sys(), integ(), dt(), current_u() { 
      setName(aName);
    };
    
    /**
     * Parametrized constructor.
     * \param aSys The continuous-time state-space system to integrate.
     * \param aInteg The numerical integrator to use to compute the state transitions.
     * \param aDt The time-step of this discrete-time system (not the integration time-step).
     */
    num_int_dtnl_sys(const CTSystem& aSys, 
		     const NumIntegrator& aInteg, 
		     const time_difference_type& aDt, 
		     const std::string& aName = "") : sys(aSys), integ(aInteg), dt(aDt), current_u() { 
      setName(aName);
    };
    
    
    virtual void RK_CALL computeStateRate(double aTime,const ReaK::vect_n<value_type>& aState, ReaK::vect_n<value_type>& aStateRate) {
      point_type p = aState;
      aStateRate = sys.get_state_derivative(p,current_u,aTime);
    };
    
    /**
     * Returns the time-step of this discrete-time system.
     * \return The time-step of this discrete-time system.
     */
    time_difference_type get_time_step() const { return dt; };
    
    /**
     * Returns next state of the system given the current state, input and time.
     * \param p The current state.
     * \param u The current input.
     * \param t The current time.
     * \return The next state, at t + get_time_step().
     */
    point_type get_next_state(const point_type& p, const input_type& u, const time_type& t = 0) { RK_UNUSED(t);
      integ.setTime(t);
      integ.clearStateVector();
      integ.addStateElements(p);
      current_u = u;
      boost::shared_ptr< state_rate_function<value_type> > aThis(this,null_deleter());
      integ.setStateRateFunc(aThis);
      integ.integrate(t + dt);
      point_type result(p); size_type i = 0;
      for(typename std::vector<value_type>::const_iterator it = integ.getStateBegin(); it != integ.getStateEnd(); ++it, ++i)
	result[i] = *it;
      return result;
    };
    
    /**
     * Returns output of the system given the current state, input and time.
     * \param p The current state.
     * \param u The current input.
     * \param t The current time.
     * \return The current output.
     */
    output_type get_output(const point_type& p, const input_type& u, const time_type& t = 0) { RK_UNUSED(t);
      return sys.get_output(p,u,t);
    };
    
    /**
     * Adjusts the state by adding a state difference to it.
     */
    point_type adjust(const point_type& p, const point_difference_type& dp) {
      return sys.adjust(p,dp);
    };
        
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(sys)
         & RK_SERIAL_SAVE_WITH_NAME(integ)
         & RK_SERIAL_SAVE_WITH_NAME(dt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(sys)
         & RK_SERIAL_LOAD_WITH_NAME(integ)
         & RK_SERIAL_LOAD_WITH_NAME(dt);
    };
    
    RK_RTTI_MAKE_CONCRETE_2BASE(self,0xC2300006,1,"num_int_dtnl_sys",named_object,state_rate_function<value_type>)
    
    
};






};

};


#endif











