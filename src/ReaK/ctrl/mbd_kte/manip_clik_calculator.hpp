/**
 * \file manip_clik_calculator.hpp
 * 
 * 
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

#include "optimization/nl_interior_points_methods.hpp"
#include "optimization/function_types.hpp"

namespace ReaK {

namespace kte {



/**
 * This class is a helper of the manipulator_kinematics_model class which is used to perform 
 * a closed-loop inverse kinematics (CLIK) calculation to compute the joint positions and 
 * velocities necessary for a given set of dependent positions and velocities (end-effector).
 * The inverse kinematics calculation is done using a non-linear constrained optimization 
 * method (the nlip_newton_tr_factory method, which is a non-linear interior-point Newton 
 * method based on a trust-region search strategy, this method has shown rapid convergence for 
 * inverse kinematics problems). This class automatically builds inequality constraints based on
 * the joint limits provided. Then, the class builds equality constraints such that the end-effector
 * pose and twist matches the desired values. Finally, one can specify a cost function to be minimized
 * if there is sufficient redundancy to permit such optimization of the resulting configuration.
 */
class manip_clik_calculator {
  public:
    
    direct_kinematics_model* model;
    
    shared_ptr<const optim::cost_evaluator> cost_eval;
    
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
    manip_clik_calculator(direct_kinematics_model* aModel, 
                          const shared_ptr<const optim::cost_evaluator>& aCostEvaluator = shared_ptr<const optim::cost_evaluator>(),
                          double aMaxRadius = 1.0, 
                          double aMu = 0.1, 
                          unsigned int aMaxIter = 300, 
                          double aTol = 1e-6, 
                          double aEta = 1e-3, 
                          double aTau = 0.99) : 
                          model(aModel), cost_eval(aCostEvaluator),
                          max_radius(aMaxRadius), mu(aMu), max_iter(aMaxIter),
                          tol(aTol), eta(aEta), tau(aTau) { 
      if(!model)
        throw optim::improper_problem("CLIK error: The model pointer cannot be null!");
      
      lower_bounds.resize(model->getJointPositionsCount() + model->getJointVelocitiesCount());
      upper_bounds.resize(model->getJointPositionsCount() + model->getJointVelocitiesCount());
      for(std::size_t i = 0; i < lower_bounds.size(); ++i) {
        lower_bounds[i] = -std::numeric_limits<double>::infinity();
        upper_bounds[i] =  std::numeric_limits<double>::infinity();
      };
      
      
    };
    
    
    /**
     * Default destructor.
     */
    ~manip_clik_calculator() { };
    
    void readDesiredFromModel() {
      
      desired_gen_coords.resize( model->DependentCoords().size() );
      for(std::size_t i = 0; i < desired_gen_coords.size(); ++i)
        desired_gen_coords[i] = *(model->DependentCoords()[i]->mFrame);
      
      desired_frame_2D.resize( model->DependentFrames2D().size() );
      for(std::size_t i = 0; i < desired_frame_2D.size(); ++i)
        desired_frame_2D[i] = *(model->DependentFrames2D()[i]->mFrame);
      
      desired_frame_3D.resize( model->DependentFrames3D().size() );
      for(std::size_t i = 0; i < desired_frame_3D.size(); ++i)
        desired_frame_3D[i] = *(model->DependentFrames3D()[i]->mFrame);
      
    };
    
    
    vect_n<double> readJointStatesFromModel() {
      vect_n<double> x(model->getJointPositionsCount() + model->getJointVelocitiesCount());
      
      vect_ref_view< vect_n<double> > pos_x = x[range(0,model->getJointPositionsCount()-1)];
      manip_kin_mdl_joint_io(model).getJointPositions(pos_x);
      
      vect_ref_view< vect_n<double> > vel_x = x[range(model->getJointPositionsCount(),model->getJointPositionsCount() + model->getJointVelocitiesCount() - 1)];
      manip_kin_mdl_joint_io(model).getJointVelocities(vel_x);
      
      return x;
    };
    
    void writeJointStatesToModel(const vect_n<double>& x) const {
      manip_kin_mdl_joint_io(model).setJointPositions(x[range(0,model->getJointPositionsCount()-1)]);
      manip_kin_mdl_joint_io(model).setJointVelocities(x[range(model->getJointPositionsCount(),model->getJointPositionsCount() + model->getJointVelocitiesCount() - 1)]);
    };
    
    void runOptimizer( vect_n<double>& x ) {
      shared_ptr<const optim::cost_evaluator> tmp_cost_eval = cost_eval;
      if(!cost_eval) {
        tmp_cost_eval = shared_ptr<const optim::cost_evaluator>(
          new optim::quadratic_cost_evaluator(
            vect_n<double>(x.size(), double(0.0)),
            mat<double,mat_structure::symmetric>(
              ( mat<double,mat_structure::nil>(model->getJointPositionsCount(), model->getJointPositionsCount() + model->getJointVelocitiesCount()) |
              ( mat<double,mat_structure::nil>(model->getJointVelocitiesCount(), model->getJointPositionsCount()) & mat<double,mat_structure::identity>(model->getJointVelocitiesCount()) ) )
            ) 
          )
        );
      };
      
    
    
    typedef optim::nlip_newton_tr_factory<optim::oop_cost_function,
                                          optim::oop_cost_grad,
                                          optim::oop_cost_hess,
                                          double,
                                          eq_function,
                                          eq_jac_filler,
                                          ineq_function,
                                          ineq_jac_filler> optim_factory_type;
//     typedef optim::nlip_quasi_newton_tr_factory<optim::oop_cost_function,
//                                                 optim::oop_cost_grad,
//                                              double,
//                                              eq_function,
//                                              eq_jac_filler,
//                                              ineq_function,
//                                              ineq_jac_filler> optim_factory_type;
    
    optim_factory_type optimizer(
                            optim::oop_cost_function(tmp_cost_eval),
                            optim::oop_cost_grad(tmp_cost_eval),
                            optim::oop_cost_hess(tmp_cost_eval),
                            max_radius, mu, max_iter,
                            eq_function(this), eq_jac_filler(this),
                            ineq_function(this), ineq_jac_filler(this),
                            tol, eta, tau);

      
      
      optimizer( x );
    };
    
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
    void solveInverseKinematics() {
      
      readDesiredFromModel();
      
      vect_n<double> x = readJointStatesFromModel();
      
      runOptimizer(x);
      
      writeJointStatesToModel(x);
    };
    
    
    
    
    
};




};

};

#endif











