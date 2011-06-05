/**
 * \file fixed_step_integrators.hpp
 *
 * The following library implements numerical methods for integration of systems
 * of ordinary differential equations using fixed time steps. The implementations
 * are done as described in the following books:\n\n
 *
 * Burden R.L. and Faires J.D., "Numerical Analysis", 8th Edition, Thomson, 2005.\n\n
 *
 * Ascher U.M. and Petzold L.R., "Computer Methods for Ordinary Differential Equations
 * Differential-Algebraic Equations", Society for Industrial and Applied Mathematics, 1998.\n\n
 *
 * The methods implemented are:\n\n
 *
 *   - Backward Euler (euler_integrator)\n
 *   - Midpoint (midpoint_integrator)\n
 *   - Runge-Kutta order 4 (runge_kutta4_integrator)\n
 *   - Runge-Kutta order 5 (as of Fehlberg) (runge_kutta5_integrator)\n
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date july 2010
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

#ifndef FIXED_STEP_INTEGRATORS_HPP
#define FIXED_STEP_INTEGRATORS_HPP

#include "integrator.hpp"

namespace ReaK {


template <class T>
class euler_integrator : public integrator<T> {
  public:

    virtual void RK_CALL integrate(double aEndTime) {
      if ((integrator<T>::mGetStateRate.expired()) || 
          (integrator<T>::mState.q.size() == 0) || 
          (integrator<T>::mStepSize == 0.0) || 
          ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) || 
          ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) 
        throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

      boost::shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
      if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);
  
      func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  
      while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

        integrator<T>::mState += integrator<T>::mStateRate * T(integrator<T>::mStepSize);
    
        integrator<T>::mTime += integrator<T>::mStepSize;

        func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
      };
    };

    euler_integrator(const std::string& aName = "") : integrator<T>(aName) { };
    euler_integrator(const std::string& aName,
                     const ReaK::vect_n<T>& aState,
                     double aStartTime,
                     double aStepSize,
                     boost::weak_ptr< state_rate_function<T> > aGetStateRate) :
                     integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate) { };
    virtual ~euler_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      integrator<T>::save(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      integrator<T>::load(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };

    typedef euler_integrator<T> self;
    typedef integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2210001,1,"euler_integrator",base)
};



template <class T>
class midpoint_integrator : public integrator<T> {
  public:

    virtual void RK_CALL integrate(double aEndTime);

    midpoint_integrator(const std::string& aName = "") : integrator<T>(aName) { };
    midpoint_integrator(const std::string& aName,
                        const ReaK::vect_n<T>& aState,
                        double aStartTime,
                        double aStepSize,
                        boost::weak_ptr< state_rate_function<T> > aGetStateRate) :
                        integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate) { };
    virtual ~midpoint_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      integrator<T>::save(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      integrator<T>::load(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };

    typedef midpoint_integrator<T> self;
    typedef integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2210002,1,"midpoint_integrator",base)
    
};

template <class T>
void RK_CALL midpoint_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) || 
      (integrator<T>::mState.q.size() == 0) || 
      (integrator<T>::mStepSize == 0.0) || 
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) || 
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) 
    throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  boost::shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
  if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);
  
  vect_n<T> w(integrator<T>::mState.q.size());

  func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  
  while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    w = integrator<T>::mState + integrator<T>::mStateRate * T(integrator<T>::mStepSize * 0.5);
    integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,w,integrator<T>::mStateRate);
    
    integrator<T>::mState += integrator<T>::mStateRate * T(integrator<T>::mStepSize);
    integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  };
};




template <class T>
class runge_kutta4_integrator : public integrator<T> {
  public:

    virtual void RK_CALL integrate(double aEndTime);

    runge_kutta4_integrator(const std::string& aName = "") : integrator<T>(aName) { };
    runge_kutta4_integrator(const std::string& aName,
                            const ReaK::vect_n<T>& aState,
                            double aStartTime,
                            double aStepSize,
                            boost::weak_ptr< state_rate_function<T> > aGetStateRate) :
                            integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate) { };
    virtual ~runge_kutta4_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      integrator<T>::save(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      integrator<T>::load(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };

    typedef runge_kutta4_integrator<T> self;
    typedef integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2210003,1,"runge_kutta4_integrator",base)
};


template <class T>
void RK_CALL runge_kutta4_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) || 
      (integrator<T>::mState.q.size() == 0) || 
      (integrator<T>::mStepSize == 0.0) || 
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) || 
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) 
    throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  boost::shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
  if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);
  
  vect_n<T> w(integrator<T>::mState.q.size());
  vect_n<T> k1(integrator<T>::mState.q.size());
  vect_n<T> k2(integrator<T>::mState.q.size());
  vect_n<T> k3(integrator<T>::mState.q.size());

  func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  
  while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    w = integrator<T>::mState;
    integrator<T>::mState += (k1 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(0.5);

    integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = w + (k2 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(0.5);

    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = w + (k3 = integrator<T>::mStateRate * T(integrator<T>::mStepSize));
    
    integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState += (k1 + k2 * T(2.0) + integrator<T>::mStateRate * T(integrator<T>::mStepSize)) / T(6.0) - k3 * T(2.0/3.0);
    
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  };
};







template <class T>
class runge_kutta5_integrator : public integrator<T> {
  public:

    virtual void RK_CALL integrate(double aEndTime);

    runge_kutta5_integrator(const std::string& aName = "") : integrator<T>(aName) { };
    runge_kutta5_integrator(const std::string& aName,
                            const ReaK::vect_n<T>& aState,
                            double aStartTime,
                            double aStepSize,
                            boost::weak_ptr< state_rate_function<T> > aGetStateRate) :
                            integrator<T>(aName, aState, aStartTime, aStepSize, aGetStateRate) { };
    virtual ~runge_kutta5_integrator() { };

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      integrator<T>::save(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      integrator<T>::load(A,integrator<T>::getStaticObjectType()->TypeVersion());
    };

    typedef runge_kutta5_integrator<T> self;
    typedef integrator<T> base;
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2210004,1,"runge_kutta5_integrator",base)
};


template <class T>
void RK_CALL runge_kutta5_integrator<T>::integrate(double aEndTime) {
  if ((integrator<T>::mGetStateRate.expired()) || 
      (integrator<T>::mState.q.size() == 0) || 
      (integrator<T>::mStepSize == 0.0) || 
      ((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime > aEndTime)) || 
      ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime < aEndTime))) 
    throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  boost::shared_ptr< state_rate_function<T> > func_ptr = integrator<T>::mGetStateRate.lock();
  if (!func_ptr) throw impossible_integration(integrator<T>::mTime,aEndTime,integrator<T>::mStepSize);

  vect_n<T> w(integrator<T>::mState.q.size());
  vect_n<T> k1(integrator<T>::mState.q.size());
  vect_n<T> k2(integrator<T>::mState.q.size());
  vect_n<T> k3(integrator<T>::mState.q.size());
  vect_n<T> k4(integrator<T>::mState.q.size());
  vect_n<T> k5(integrator<T>::mState.q.size());
  
  func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);

  while(((integrator<T>::mStepSize > 0.0) && (integrator<T>::mTime < aEndTime)) || ((integrator<T>::mStepSize < 0.0) && (integrator<T>::mTime > aEndTime))) {

    w = integrator<T>::mState;
    integrator<T>::mState += (k1 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(0.25);

    integrator<T>::mTime += integrator<T>::mStepSize * 0.25;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState += ((k2 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(9.0) - k1 * T(5.0)) / T(32.0);
    
    integrator<T>::mTime += integrator<T>::mStepSize * 0.125;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState += (k1 * T(276165.0) - k2 * T(1250865.0) + (k3 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(1167360.0)) / T(351520.0);
    
    integrator<T>::mTime += 57.0 * integrator<T>::mStepSize / 104.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = w + k1 * T(439.0 / 216.0) - k2 * T(8.0) + k3 * T(3680.0 / 513.0) - (k4 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(845.0 / 4104.0);

    integrator<T>::mTime += integrator<T>::mStepSize / 13.0;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = w - k1 * T(8.0 / 27.0) + k2 * T(2.0) - k3 * T(3544.0 / 2565.0) + k4 * T(1859.0 / 4104.0) - (k5 = integrator<T>::mStateRate * T(integrator<T>::mStepSize)) * T(11.0 / 40.0);

    integrator<T>::mTime -= integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
    integrator<T>::mState = w + k1 * T(16.0 / 135.0) + k3 * T(6656.0 / 12825.0) + k4 * T(28561.0 / 56430.0) - k5 * T(9.0 / 50.0) + integrator<T>::mStateRate * T(2.0 * integrator<T>::mStepSize / 55.0);

    integrator<T>::mTime += integrator<T>::mStepSize * 0.5;
    func_ptr->computeStateRate(integrator<T>::mTime,integrator<T>::mState,integrator<T>::mStateRate);
  };
};







};

#endif



