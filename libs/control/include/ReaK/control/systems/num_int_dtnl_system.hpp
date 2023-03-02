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

#ifndef REAK_NUM_INT_DTNL_SYSTEM_HPP
#define REAK_NUM_INT_DTNL_SYSTEM_HPP

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/math/integrators/fixed_step_integrators.hpp>
#include <ReaK/math/integrators/integrator.hpp>
#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

#include "discrete_sss_concept.hpp"
#include "state_space_sys_concept.hpp"

namespace ReaK::ctrl {

/**
 * This class template can be used to integrate, numerically, a continuous-time
 * state-space system (SSSystemConcept) to produce a discrete-time state-space system (DiscreteSSSConcept).
 * \tparam CTSystem The continuous-time state-space system type, should model SSSystemConcept.
 * \tparam NumIntegrator The numerical integrator type to be used.
 */
template <typename CTSystem,
          typename NumIntegrator = euler_integrator<typename vect_traits<
              typename ss_system_traits<CTSystem>::point_type>::value_type>>
class num_int_dtnl_sys : public named_object {
 public:
  using self = num_int_dtnl_sys<CTSystem, NumIntegrator>;

  using ct_system_ptr = std::shared_ptr<CTSystem>;

  using point_type = typename ss_system_traits<CTSystem>::point_type;
  using point_difference_type =
      typename ss_system_traits<CTSystem>::point_difference_type;
  using point_derivative_type =
      typename ss_system_traits<CTSystem>::point_derivative_type;

  using value_type = vect_value_type_t<point_type>;
  using size_type = typename vect_traits<point_type>::size_type;

  using time_type = typename ss_system_traits<CTSystem>::time_type;
  using time_difference_type =
      typename ss_system_traits<CTSystem>::time_difference_type;

  using input_type = typename ss_system_traits<CTSystem>::input_type;
  using output_type = typename ss_system_traits<CTSystem>::output_type;

  static constexpr std::size_t dimensions =
      ss_system_traits<CTSystem>::dimensions;
  static constexpr std::size_t input_dimensions =
      ss_system_traits<CTSystem>::input_dimensions;
  static constexpr std::size_t output_dimensions =
      ss_system_traits<CTSystem>::output_dimensions;

 private:
  ct_system_ptr sys;
  NumIntegrator integ;
  time_difference_type dt;

  template <typename StateSpaceType>
  class rate_function_impl : public state_rate_function<value_type> {
   public:
    const StateSpaceType& state_space;
    ct_system_ptr sys;
    input_type current_u;

    rate_function_impl(const StateSpaceType& aStateSpace,
                       const ct_system_ptr& aSys,
                       const input_type& aCurrentInput)
        : state_space(aStateSpace), sys(aSys), current_u(aCurrentInput) {}

    virtual void computeStateRate(double aTime,
                                  const ReaK::vect_n<value_type>& aState,
                                  ReaK::vect_n<value_type>& aStateRate) {
      using ReaK::from_vect;
      using ReaK::to_vect;
      point_type p = from_vect<point_type>(aState);
      aStateRate = to_vect<value_type>(
          sys->get_state_derivative(state_space, p, current_u, aTime));
    }
  };

 public:
  /**
   * Returns the dimensions of the state vectors.
   * \return The dimensions of the state vectors.
   */
  size_type get_state_dimensions() const { return sys->get_state_dimensions(); }

  /**
   * Returns the dimensions of the input vectors.
   * \return The dimensions of the input vectors.
   */
  size_type get_input_dimensions() const { return sys->get_input_dimensions(); }

  /**
   * Returns the dimensions of the output vectors.
   * \return The dimensions of the output vectors.
   */
  size_type get_output_dimensions() const {
    return sys->get_output_dimensions();
  }

  /**
   * Default constructor.
   */
  num_int_dtnl_sys(const std::string& aName = "")
      : sys(std::make_shared<CTSystem>()), integ(), dt() {
    setName(aName);
  }

  /**
   * Parametrized constructor.
   * \param aSys The continuous-time state-space system to integrate.
   * \param aInteg The numerical integrator to use to compute the state transitions.
   * \param aDt The time-step of this discrete-time system (not the integration time-step).
   */
  num_int_dtnl_sys(const ct_system_ptr& aSys, const NumIntegrator& aInteg,
                   const time_difference_type& aDt,
                   const std::string& aName = "")
      : sys(aSys), integ(aInteg), dt(aDt) {
    setName(aName);
  }

  /**
   * Returns the time-step of this discrete-time system.
   * \return The time-step of this discrete-time system.
   */
  time_difference_type get_time_step() const { return dt; }

  /**
   * Sets the time-step of this discrete-time system.
   * \param aDt The time-step of this discrete-time system.
   */
  void set_time_step(time_difference_type aDt) { dt = aDt; }

  /**
   * Returns next state of the system given the current state, input and time.
   * \tparam StateSpaceType The state-space topology type on which the underlying system operates.
   * \param state_space The state-space topology on which the underlying system operates.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The next state, at t + get_time_step().
   */
  template <typename StateSpaceType>
  point_type get_next_state(const StateSpaceType& state_space,
                            const point_type& p, const input_type& u,
                            const time_type& t = 0) {
    integ.setTime(t);
    integ.clearStateVector();
    vect_n<value_type> result = to_vect<value_type>(p);
    integ.addStateElements(result);
    auto temp = std::make_shared<rate_function_impl<StateSpaceType>>(
        state_space, sys, u);
    integ.setStateRateFunc(temp);
    integ.integrate(t + dt);
    size_type i = 0;
    for (auto it = integ.getStateBegin(); it != integ.getStateEnd();
         ++it, ++i) {
      result[i] = *it;
    }
    integ.setStateRateFunc(std::shared_ptr<state_rate_function<value_type>>());
    return from_vect<point_type>(result);
  }

  /**
   * Returns output of the system given the current state, input and time.
   * \tparam StateSpaceType The state-space topology type on which the underlying system operates.
   * \param state_space The state-space topology on which the underlying system operates.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The current output.
   */
  template <typename StateSpaceType>
  output_type get_output(const StateSpaceType& state_space, const point_type& p,
                         const input_type& u, const time_type& t = 0) {
    RK_UNUSED(t);
    return sys->get_output(state_space, p, u, t);
  }

  virtual void save(ReaK::serialization::oarchive& aA, unsigned int) const {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(sys) & RK_SERIAL_SAVE_WITH_NAME(integ) &
        RK_SERIAL_SAVE_WITH_NAME(dt);
  }
  virtual void load(ReaK::serialization::iarchive& aA, unsigned int) {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(sys) & RK_SERIAL_LOAD_WITH_NAME(integ) &
        RK_SERIAL_LOAD_WITH_NAME(dt);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300006, 1, "num_int_dtnl_sys",
                              named_object)
};

/**
 * This class template specialization can be used to integrate, numerically, a continuous-time
 * system represented as a ReaK::state_rate_function OOP interface to produce a
 * discrete-time state-space system (DiscreteSSSConcept).
 *
 * \tparam T The value-type of the state_rate_function class template.
 * \tparam NumIntegrator The numerical integrator type to be used.
 */
template <typename T, template <typename> class NumIntegrator>
class num_int_dtnl_sys<state_rate_function_with_io<T>, NumIntegrator<T>>
    : public named_object {
 public:
  using self =
      num_int_dtnl_sys<state_rate_function_with_io<T>, NumIntegrator<T>>;
  using value_type = T;
  using size_type = typename vect_traits<vect_n<value_type>>::size_type;

  using ct_system_ptr =
      std::shared_ptr<state_rate_function_with_io<value_type>>;

  using point_type = vect_n<value_type>;
  using point_difference_type = vect_n<value_type>;
  using point_derivative_type = vect_n<value_type>;

  using time_type = double;
  using time_difference_type = double;

  using input_type = vect_n<value_type>;
  using output_type = vect_n<value_type>;

  static constexpr std::size_t dimensions = 0;
  static constexpr std::size_t input_dimensions = 0;
  static constexpr std::size_t output_dimensions = 0;

 private:
  ct_system_ptr sys;
  NumIntegrator<value_type> integ;
  time_difference_type dt;

 public:
  /**
   * Returns the dimensions of the state vectors.
   * \return The dimensions of the state vectors.
   */
  size_type get_state_dimensions() const { return sys->get_state_dimensions(); }

  /**
   * Returns the dimensions of the input vectors.
   * \return The dimensions of the input vectors.
   */
  size_type get_input_dimensions() const { return sys->get_input_dimensions(); }

  /**
   * Returns the dimensions of the output vectors.
   * \return The dimensions of the output vectors.
   */
  size_type get_output_dimensions() const {
    return sys->get_output_dimensions();
  }

  /**
   * Default constructor.
   */
  num_int_dtnl_sys(const std::string& aName = "") : sys(), integ(), dt() {
    setName(aName);
  }

  /**
   * Parametrized constructor.
   * \param aSys The continuous-time state-space system to integrate.
   * \param aInteg The numerical integrator to use to compute the state transitions.
   * \param aDt The time-step of this discrete-time system (not the integration time-step).
   */
  num_int_dtnl_sys(const ct_system_ptr& aSys,
                   const NumIntegrator<value_type>& aInteg,
                   const time_difference_type& aDt,
                   const std::string& aName = "")
      : sys(aSys), integ(aInteg), dt(aDt) {
    setName(aName);
  }

  /**
   * Returns the time-step of this discrete-time system.
   * \return The time-step of this discrete-time system.
   */
  time_difference_type get_time_step() const { return dt; }

  /**
   * Sets the time-step of this discrete-time system.
   * \param aDt The time-step of this discrete-time system.
   */
  void set_time_step(time_difference_type aDt) { dt = aDt; }

  /**
   * Returns next state of the system given the current state, input and time.
   * \tparam StateSpaceType The state-space topology type on which the underlying system operates.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The next state, at t + get_time_step().
   */
  template <typename StateSpaceType>
  point_type get_next_state(const StateSpaceType&, const point_type& p,
                            const input_type& u, const time_type& t = 0) {
    if (!sys) {
      return p;
    }
    integ.setTime(t);
    integ.clearStateVector();
    vect_n<value_type> result = to_vect<value_type>(p);
    integ.addStateElements(p);
    sys->setInput(u);
    integ.setStateRateFunc(sys);
    integ.integrate(t + dt);
    size_type i = 0;
    for (typename std::vector<value_type>::const_iterator it =
             integ.getStateBegin();
         it != integ.getStateEnd(); ++it, ++i) {
      result[i] = *it;
    }
    integ.setStateRateFunc(
        std::shared_ptr<state_rate_function_with_io<value_type>>());
    return from_vect<point_type>(result);
  }

  /**
   * Returns output of the system given the current state, input and time.
   * \tparam StateSpaceType The state-space topology type on which the underlying system operates.
   * \param p The current state.
   * \param u The current input.
   * \param t The current time.
   * \return The current output.
   */
  template <typename StateSpaceType>
  output_type get_output(const StateSpaceType&, const point_type& p,
                         const input_type& u, const time_type& t = 0) {
    RK_UNUSED(t);
    output_type y;

    if (sys) {
      sys->setInput(u);
      sys->computeOutput(t, p, y);
    }
    return y;
  }

  virtual void save(ReaK::serialization::oarchive& aA, unsigned int) const {
    ReaK::named_object::save(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_SAVE_WITH_NAME(sys) & RK_SERIAL_SAVE_WITH_NAME(integ) &
        RK_SERIAL_SAVE_WITH_NAME(dt);
  }
  virtual void load(ReaK::serialization::iarchive& aA, unsigned int) {
    ReaK::named_object::load(
        aA, ReaK::named_object::getStaticObjectType()->TypeVersion());
    aA& RK_SERIAL_LOAD_WITH_NAME(sys) & RK_SERIAL_LOAD_WITH_NAME(integ) &
        RK_SERIAL_LOAD_WITH_NAME(dt);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2300006, 1, "num_int_dtnl_sys",
                              named_object)
};

}  // namespace ReaK::ctrl

#endif
