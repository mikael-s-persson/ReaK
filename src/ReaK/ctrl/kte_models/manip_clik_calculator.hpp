/**
 * \file manip_clik_calculator.hpp
 * 
 * This library declares a helper class that can be used to solve the inverse kinematics problem using a 
 * closed-loop formulation (based on a general non-linear optimization method (trust-region Newton sequential 
 * quadratic programming algorithm)). The main class (manip_clik_calculator) takes a cost function (for redundancy 
 * resolution, if necessary), a direct kinematics model, bounds on the joints, and a desired set of dependent poses.
 * The clik calculator tries to solve the general optimization problem of keeping the joint values (incl. velocities) 
 * within the given bounds (unless they are set to +- infinity) using inequality constraints, then tries to match 
 * the desired dependent poses using equality constraints, and finally, tries to minimize the given cost function,
 * which will only have an effect if there are degrees of freedom available to do so.
 * 
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_MANIP_CLIK_CALCULATOR_HPP
#define REAK_MANIP_CLIK_CALCULATOR_HPP

#include "direct_kinematics_model.hpp"
#include "manip_kinematics_helper.hpp"

#include "optimization/optim_exceptions.hpp"
#include "optimization/function_types.hpp"

#include <cmath>

namespace ReaK {

namespace kte {


/**
 * This function manufactures a basic quadratic cost function for use in the CLIK calculator.
 * This quadratic cost function normalizes the range of joint values (between bounds) and 
 * centered around the given preferred posture.
 */
shared_ptr< optim::quadratic_cost_evaluator > create_clik_quad_cost(const vect_n<double>& aPreferredPosture, 
                                                                    const vect_n<double>& aLowerBounds,
                                                                    const vect_n<double>& aUpperBounds,
                                                                    const direct_kinematics_model& aModel);


/**
 * This cost evaluator attempts to keep joints from being too straight, i.e., tries 
 * to keep them at a +/- 90 degrees angle. The cost function itself is a sin-squared 
 * function applied on each joint in its internal list. This is a useful cost function 
 * which is a cheap approximation of a singularity avoidance cost function (which would 
 * normally involve evaluating a manipulability index (or performance / dexterity index) 
 * for the whole manipulator, and that is usually very expensive, involving singular 
 * value decompositions of the Jacobian matrix, or similarly expensive metrics).
 */
class clik_bent_joints_cost_eval : public optim::cost_evaluator {
  public:
    std::vector<int> joint_ids; ///< Holds the indices of the joints to keep bent.
    
    /**
     * Parametrized Constructor.
     * \param aJointIDs The indices of the joints to keep bent.
     */
    clik_bent_joints_cost_eval(const std::vector<int>& aJointIDs) : joint_ids(aJointIDs) { };
    
    explicit clik_bent_joints_cost_eval(int J1 = -1, int J2 = -1, int J3 = -1, int J4 = -1) {
      if(J1 >= 0)
        joint_ids.push_back(J1);
      if(J2 >= 0)
        joint_ids.push_back(J2);
      if(J3 >= 0)
        joint_ids.push_back(J3);
      if(J4 >= 0)
        joint_ids.push_back(J4);
    };
    
    clik_bent_joints_cost_eval(int J1, int J2, int J3, int J4, int J5, int J6 = -1, int J7 = -1, int J8 = -1) {
      if(J1 >= 0)
        joint_ids.push_back(J1);
      if(J2 >= 0)
        joint_ids.push_back(J2);
      if(J3 >= 0)
        joint_ids.push_back(J3);
      if(J4 >= 0)
        joint_ids.push_back(J4);
      if(J5 >= 0)
        joint_ids.push_back(J5);
      if(J6 >= 0)
        joint_ids.push_back(J6);
      if(J7 >= 0)
        joint_ids.push_back(J7);
      if(J8 >= 0)
        joint_ids.push_back(J8);
    };
    
    virtual double compute_cost(const vect_n<double>& x) const {
      double sum = 0.0;
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
        double s = std::sin( x[joint_ids[i]] );
        sum += 1.0 - s * s;
      };
      return sum;
    };
    
    virtual vect_n<double> compute_cost_grad(const vect_n<double>& x) const {
      vect_n<double> result(x.size(),0.0);
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
        result[ joint_ids[i] ] = -2.0 * std::sin( x[joint_ids[i]] ) * std::cos( x[joint_ids[i]] );
      };
      return result;
    };
    
    virtual void compute_cost_hessian(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double, const vect_n<double>&) const {
      H.set_row_count(x.size());
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
        double s = std::sin( x[joint_ids[i]] ); s *= s;
        double c = std::cos( x[joint_ids[i]] ); c *= c;
        H(joint_ids[i], joint_ids[i]) = 2.0 * (s - c);
      };
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      optim::cost_evaluator::save(A,optim::cost_evaluator::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(joint_ids);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      optim::cost_evaluator::load(A,optim::cost_evaluator::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(joint_ids);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(clik_bent_joints_cost_eval,0xC1500004,1,"clik_bent_joints_cost_eval",optim::cost_evaluator)
    
};


/**
 * This function manufactures a bent-joint cost function.
 */
inline
shared_ptr< clik_bent_joints_cost_eval > create_clik_bent_joints(const std::vector<int>& aJointIDs) {
  return shared_ptr< clik_bent_joints_cost_eval >(new clik_bent_joints_cost_eval(aJointIDs));
};

/**
 * This function manufactures a bent-joint cost function.
 */
inline
shared_ptr< clik_bent_joints_cost_eval > create_clik_bent_joints(int J1 = -1, int J2 = -1, int J3 = -1, int J4 = -1) {
  return shared_ptr< clik_bent_joints_cost_eval >(new clik_bent_joints_cost_eval(J1,J2,J3,J4));
};

/**
 * This function manufactures a bent-joint cost function.
 */
inline
shared_ptr< clik_bent_joints_cost_eval > create_clik_bent_joints(int J1, int J2, int J3, int J4, int J5, int J6 = -1, int J7 = -1, int J8 = -1) {
  return shared_ptr< clik_bent_joints_cost_eval >(new clik_bent_joints_cost_eval(J1,J2,J3,J4,J5,J6,J7,J8));
};



/**
 * This function manufactures a mixed cost function which adds the cost of any two cost evaluators.
 */
inline
shared_ptr< optim::added_cost_evaluator > create_clik_mixed_cost(const shared_ptr< optim::cost_evaluator >& aCostFunc1,
                                                                 const shared_ptr< optim::cost_evaluator >& aCostFunc2) {
  return shared_ptr< optim::added_cost_evaluator >(new optim::added_cost_evaluator(aCostFunc1, aCostFunc2));
};



/**
 * This class is a helper class which is used to perform 
 * a closed-loop inverse kinematics (CLIK) calculation to compute the joint positions and 
 * velocities necessary for a given set of dependent positions and velocities (end-effector).
 * The inverse kinematics calculation is done using a non-linear constrained optimization 
 * method (the nlip_newton_tr_factory method, which is a non-linear interior-point Newton 
 * method based on a trust-region search strategy, this method has shown rapid convergence for 
 * inverse kinematics problems). This class automatically builds inequality constraints based on
 * the joint limits provided, to keep joints within a valid range. Then, the class builds equality 
 * constraints such that the end-effector pose and twist matches the desired values. Finally, one 
 * can specify a cost function to be minimized if there is sufficient redundancy to permit such 
 * optimization of the resulting configuration.
 */
class manip_clik_calculator : public shared_object {
  public:
    
    shared_ptr<direct_kinematics_model> model;
    
    shared_ptr<optim::cost_evaluator> cost_eval;
    
    vect_n<double> lower_bounds;
    vect_n<double> upper_bounds;
    
    std::vector< gen_coord<double> > desired_gen_coords;
    std::vector< frame_2D<double> > desired_frame_2D;
    std::vector< frame_3D<double> > desired_frame_3D;
    
    double max_radius; 
    double mu;
    unsigned int max_iter;
    double tol;
    double eta;
    double tau;
    
    /**
     * Default constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     */
    manip_clik_calculator(const shared_ptr<direct_kinematics_model>& aModel = shared_ptr<direct_kinematics_model>(), 
                          const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                          double aMaxRadius = 1.0, 
                          double aMu = 0.1, 
                          unsigned int aMaxIter = 300, 
                          double aTol = 1e-6, 
                          double aEta = 1e-3, 
                          double aTau = 0.99) : 
                          model(aModel), cost_eval(aCostEvaluator),
                          max_radius(aMaxRadius), mu(aMu), max_iter(aMaxIter),
                          tol(aTol), eta(aEta), tau(aTau) { 
      if(model) {
        lower_bounds.resize(model->getJointPositionsCount() + model->getJointVelocitiesCount());
        upper_bounds.resize(model->getJointPositionsCount() + model->getJointVelocitiesCount());
        for(std::size_t i = 0; i < lower_bounds.size(); ++i) {
          lower_bounds[i] = -std::numeric_limits<double>::infinity();
          upper_bounds[i] =  std::numeric_limits<double>::infinity();
        };
      };
    };
    
    /**
     * Default destructor.
     */
    ~manip_clik_calculator() { };
    
    void readDesiredFromModel();
    
    vect_n<double> readJointStatesFromModel();
    
    void writeJointStatesToModel(const vect_n<double>& x) const;
    
    vect_n<double> computeStatesError(const vect_n<double>& x) const;
    
    void runOptimizer( vect_n<double>& x );
    
    /**
     * This function takes its given desired 'end-effector' states and solves the inverse 
     * kinematics problem to find the corresponding joint states. 
     * \pre The joint states at which the manipulator model is before the function call 
     *      will act as the initial guess for the inverse kinematics solution.
     *      Before calling this function, it is assumed that all the parameters of the 
     *      optimization have been properly set.
     * \post The joint states which solve the inverse kinematics problem will be set in 
     *       the manipulator model.
     */
    void solveInverseKinematics();
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(cost_eval)
        & RK_SERIAL_SAVE_WITH_NAME(lower_bounds)
        & RK_SERIAL_SAVE_WITH_NAME(upper_bounds)
        & RK_SERIAL_SAVE_WITH_NAME(max_radius)
        & RK_SERIAL_SAVE_WITH_NAME(mu)
        & RK_SERIAL_SAVE_WITH_NAME(max_iter)
        & RK_SERIAL_SAVE_WITH_NAME(tol)
        & RK_SERIAL_SAVE_WITH_NAME(eta)
        & RK_SERIAL_SAVE_WITH_NAME(tau);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(cost_eval)
        & RK_SERIAL_LOAD_WITH_NAME(lower_bounds)
        & RK_SERIAL_LOAD_WITH_NAME(upper_bounds)
        & RK_SERIAL_LOAD_WITH_NAME(max_radius)
        & RK_SERIAL_LOAD_WITH_NAME(mu)
        & RK_SERIAL_LOAD_WITH_NAME(max_iter)
        & RK_SERIAL_LOAD_WITH_NAME(tol)
        & RK_SERIAL_LOAD_WITH_NAME(eta)
        & RK_SERIAL_LOAD_WITH_NAME(tau);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(manip_clik_calculator,0xC2100051,1,"manip_clik_calculator",shared_object)
    
    
};




};

};

#endif











