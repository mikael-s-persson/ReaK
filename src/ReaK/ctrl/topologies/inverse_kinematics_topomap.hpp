/**
 * \file inverse_kinematics_topomap.hpp
 * 
 * This library provides classes that define topological mappings between a joint-space (generalized 
 * coordinates) and the end-effector frame of a serial manipulator. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
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

#ifndef REAK_INVERSE_KINEMATICS_TOPOMAP_HPP
#define REAK_INVERSE_KINEMATICS_TOPOMAP_HPP


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "kte_models/inverse_kinematics_model.hpp"
#include "kte_models/manip_clik_calculator.hpp"

#include "joint_space_topologies.hpp"
#include "joint_space_limits.hpp"
#include "se3_topologies.hpp"
#include "se2_topologies.hpp"
#include "Ndof_spaces.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"
#include "kinetostatics/gen_coord.hpp"

// needed for the KNN IK mapping. NOTE: the KNN-IK should be detached from this header file.
#include <boost/property_map/vector_property_map.hpp>
#include "path_planning/metric_space_search.hpp"

#include "direct_kinematics_topomap.hpp"
#include "inverse_kinematics_topomap_detail.hpp"

namespace ReaK {

namespace pp {



/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_inverse_kin_map : public shared_object {
  public:
    
    typedef manip_inverse_kin_map self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::inverse_kinematics_model > model; 
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model which can do the inverse kinematics calculation.
     */
    manip_inverse_kin_map(const shared_ptr< kte::inverse_kinematics_model >& aModel = shared_ptr<kte::inverse_kinematics_model>()) :
                          model(aModel) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt, space_in, model);
      
      model->doInverseMotion();
      
      detail::read_joint_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400014,1,"manip_inverse_kin_map",shared_object)
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_collection<double> >
class manip_rl_inverse_kin_map : public shared_object {
  public:
    
    typedef manip_rl_inverse_kin_map<RateLimitMap> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::inverse_kinematics_model > model;
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
     */
    manip_rl_inverse_kin_map(const shared_ptr< kte::inverse_kinematics_model >& aModel = shared_ptr<kte::inverse_kinematics_model>(),
                             const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >()) :
                             model(aModel), joint_limits_map(aJointLimitMap) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      model->doInverseMotion();
    
      typedef typename get_rate_illimited_space< OutSpace >::type NormalJointSpace;
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(), model);
      result = joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400015,1,"manip_rl_inverse_kin_map",shared_object)
    
};



/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_clik_kin_map : public shared_object {
  public:
    
    typedef manip_clik_kin_map self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::direct_kinematics_model > model; 
    /** This holds the inverse kinematics calculator factory. */
    mutable kte::manip_clik_calculator clik_calc;
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     */
    manip_clik_kin_map(const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                       const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                       double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                       double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99) :
                       model(aModel),
                       clik_calc(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt, space_in, model);
      
      clik_calc.solveInverseKinematics();
      
      detail::read_joint_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400016,1,"manip_clik_kin_map",shared_object)
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_collection<double> >
class manip_rl_clik_kin_map : public shared_object {
  public:
    
    typedef manip_rl_clik_kin_map<RateLimitMap> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::direct_kinematics_model > model;
    /** This holds the inverse kinematics calculator factory. */
    mutable kte::manip_clik_calculator clik_calc;
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     */
    manip_rl_clik_kin_map(const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                          const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >(),
                          const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                          double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                          double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99) :
                          model(aModel),
                          clik_calc(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau),
                          joint_limits_map(aJointLimitMap) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      clik_calc.solveInverseKinematics();
    
      typedef typename get_rate_illimited_space< OutSpace >::type NormalJointSpace;
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(), model);
      result = joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc)
        & RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc)
        & RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400017,1,"manip_rl_clik_kin_map",shared_object)
    
};







/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename JointStateType>
class manip_clik_fig_kin_map : public manip_clik_kin_map {
  public:
    
    typedef manip_clik_fig_kin_map<JointStateType> self;
    
    JointStateType initial_guess;
    
    /**
     * Parametrized Constructor.
     * \param aInitGuess The initial guess to act as the start for the inverse kinematics search.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     */
    manip_clik_fig_kin_map(const JointStateType& aInitGuess = JointStateType(),
                           const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                           const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                           double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                           double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99) :
                           manip_clik_kin_map(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau),
                           initial_guess(aInitGuess) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_joint_coordinates_impl(initial_guess, space_out, model);
      
      detail::write_dependent_coordinates_impl(pt, space_in, model);
      
      clik_calc.solveInverseKinematics();
      
      detail::read_joint_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      manip_clik_kin_map::save(A,manip_clik_kin_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(initial_guess);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      manip_clik_kin_map::load(A,manip_clik_kin_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(initial_guess);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001A,1,"manip_clik_fig_kin_map",manip_clik_kin_map)
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename JointStateType, typename RateLimitMap = joint_limits_collection<double> >
class manip_rl_clik_fig_kin_map : public manip_rl_clik_kin_map< RateLimitMap > {
  public:
    typedef manip_rl_clik_kin_map< RateLimitMap > base_type;
    typedef manip_rl_clik_fig_kin_map<JointStateType,RateLimitMap> self;
    
    JointStateType initial_guess;
    
    /**
     * Parametrized Constructor.
     * \param aInitGuess The initial guess to act as the start for the inverse kinematics search.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     */
    manip_rl_clik_fig_kin_map(const JointStateType& aInitGuess = JointStateType(),
                              const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                              const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >(),
                              const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                              double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                              double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99) :
                              base_type(aModel, aJointLimitMap, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau),
                              initial_guess(aInitGuess) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      typedef typename get_rate_illimited_space< OutSpace >::type NormalJointSpace;
      
      typename topology_traits<NormalJointSpace>::point_type ip_inter = 
        this->joint_limits_map->map_to_space(initial_guess, space_out, NormalJointSpace());
      detail::write_joint_coordinates_impl(ip_inter, NormalJointSpace(), this->model);
      
      detail::write_dependent_coordinates_impl(pt, space_in, this->model);
      
      this->clik_calc.solveInverseKinematics();
      
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(), this->model);
      result = this->joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(initial_guess);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(initial_guess);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001B,1,"manip_rl_clik_fig_kin_map",base_type)
    
    
};



/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_clik_rnd_restart_map : public manip_clik_kin_map {
  public:
    
    typedef manip_clik_rnd_restart_map self;
    
    std::size_t restart_count;
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     * \param aRestartCount The number of restarts to try before giving up on the inverse kinematics search.
     */
    manip_clik_rnd_restart_map(const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                               const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                               double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                               double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99,
                               std::size_t aRestartCount = 10) :
                               manip_clik_kin_map(aModel, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau),
                               restart_count(aRestartCount) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      clik_calc.readDesiredFromModel();
      
      for(std::size_t i = 0; i < restart_count; ++i) {
      
        detail::write_joint_coordinates_impl(get(random_sampler,space_out)(space_out), space_out, model);
      
        vect_n<double> x = clik_calc.readJointStatesFromModel();
        
        try {
          clik_calc.runOptimizer(x);
        } catch(singularity_error& e) { 
        } catch(maximum_iteration& e) {
        } catch(optim::infeasible_problem& e) {
        };
        
        if( norm_2( clik_calc.computeStatesError(x) ) < clik_calc.tol * 10.0 ) {
          clik_calc.writeJointStatesToModel(x);
          break;
        } else if( i == restart_count - 1 )
          throw optim::infeasible_problem("The inverse kinematics problem cannot be solved!");
      };
      
      detail::read_joint_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      manip_clik_kin_map::save(A,manip_clik_kin_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(restart_count);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      manip_clik_kin_map::load(A,manip_clik_kin_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(restart_count);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001C,1,"manip_clik_rnd_restart_map",manip_clik_kin_map)
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
template <typename RateLimitMap = joint_limits_collection<double> >
class manip_rl_clik_rnd_restart_map : public manip_rl_clik_kin_map<RateLimitMap> {
  public:
    typedef manip_rl_clik_kin_map<RateLimitMap> base_type;
    typedef manip_rl_clik_rnd_restart_map<RateLimitMap> self;
    
    std::size_t restart_count;
    
    /**
     * Parametrized Constructor.
     * \param aModel A pointer to the manipulator model on which the inverse kinematics search is applied.
     * \param aJointLimitMap A pointer to the mapping used to map rate-limited points to the normal joint-space.
     * \param aCostEvaluator The cost-evaluator to use as the objective (minimization) for the CLIK algorithm.
     * \param aMaxRadius The maximum trust-region radius to use (i.e. maximum optimization step).
     * \param aMu The initial strength of the barrier on the inequalities (initial "barrier parameter"), this parameter is positive and should start with a rather large value (relative to the scale of the function) and will be progressively decreased by the algorithm as it progresses).
     * \param aMaxIter The maximum number of iterations to perform.
     * \param aTol The tolerance on the norm of the gradient (and thus the step size).
     * \param aEta The tolerance on the decrease in order to accept a step in the trust region.
     * \param aTau The portion (close to 1.0) of a total step to do without coming too close to the inequality constraint (barrier).
     * \param aRestartCount The number of restarts to try before giving up on the inverse kinematics search.
     */
    manip_rl_clik_rnd_restart_map(const shared_ptr< kte::direct_kinematics_model >& aModel = shared_ptr<kte::direct_kinematics_model>(),
                                  const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >(),
                                  const shared_ptr<optim::cost_evaluator>& aCostEvaluator = shared_ptr<optim::cost_evaluator>(),
                                  double aMaxRadius = 1.0, double aMu = 0.1, unsigned int aMaxIter = 300, 
                                  double aTol = 1e-6, double aEta = 1e-3, double aTau = 0.99,
                                  std::size_t aRestartCount = 10) :
                                  base_type(aModel, aJointLimitMap, aCostEvaluator, aMaxRadius, aMu, aMaxIter, aTol, aEta, aTau),
                                  restart_count(aRestartCount) { };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      typedef typename get_rate_illimited_space< OutSpace >::type NormalJointSpace;
      
      detail::write_dependent_coordinates_impl(pt, space_in, this->model);
      
      this->clik_calc.readDesiredFromModel();
      
      for(std::size_t i = 0; i < restart_count; ++i) {
        
        typename topology_traits<NormalJointSpace>::point_type ip_inter = 
          this->joint_limits_map->map_to_space(get(random_sampler, space_out)(space_out), space_out, NormalJointSpace());
        detail::write_joint_coordinates_impl(ip_inter, NormalJointSpace(), this->model);
        
        vect_n<double> x = this->clik_calc.readJointStatesFromModel();
        
        try {
          this->clik_calc.runOptimizer(x);
        } catch(singularity_error& e) { 
        } catch(maximum_iteration& e) {
        } catch(optim::infeasible_problem& e) {
        };
        
        if( norm_2( this->clik_calc.computeStatesError(x) ) < this->clik_calc.tol * 10.0 ) {
          this->clik_calc.writeJointStatesToModel(x);
          break;
        } else if( i == restart_count - 1 )
          throw optim::infeasible_problem("The inverse kinematics problem cannot be solved!");
      };
      
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(), this->model);
      result = this->joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(restart_count);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(restart_count);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001D,1,"manip_rl_clik_rnd_restart_map",base_type)
    
    
};






#if 0




#if 1   
// NOTE: TODO: this is all wrong.

/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam CLIKCalcFactory The factory type for creating a CLIK calculator.
 */
template <typename CLIKCalcFactory, typename JointSpaceType, typename EESpaceType>
class manip_ik_knn_starts_map : public shared_object {
  public:
    
    typedef manip_ik_knn_starts_map<CLIKCalcFactory,JointSpaceType,EESpaceType> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    /** This holds the inverse kinematics calculator factory. */
    CLIKCalcFactory clik_calc_factory;
    
    mutable shared_ptr< kte::manip_clik_calculator > ik_calc;
    
    shared_ptr< JointSpaceType > joint_space;
    shared_ptr< EESpaceType > ee_space;
    
    std::size_t sample_count;
    std::size_t neighbor_count;
    
    typedef typename topology_traits< EESpaceType >::point_type ee_state_type;
    typedef typename topology_traits< JointSpaceType >::point_type joint_state_type;
    
    mutable boost::vector_property_map< ee_state_type > sample_points;
    mutable std::vector< joint_state_type > sample_jt_points;
    
    typedef dvp_tree< std::size_t, 
                      EESpaceType,
                      boost::vector_property_map< ee_state_type > > ee_bsp_tree_type;
    
    mutable shared_ptr< ee_bsp_tree_type > sample_tree;
    
    void initialize_ik_data(bool force_init = false) const {
      if((force_init) || (model))
        ik_calc = clik_calc_factory.create_calculator(model);
      else
        ik_calc = shared_ptr< kte::manip_clik_calculator >();
      
      if((joint_space) && (ee_space)) {
        std::vector< std::size_t > indices(sample_count);
        manip_direct_kin_map dk_map = manip_direct_kin_map(model);
        for(std::size_t i = 0; i < sample_count; ++i) {
          indices[i] = i;
          sample_jt_points[i] = get(random_sampler,*joint_space)(*joint_space);
          put(sample_points, i, dk_map.map_to_space(sample_jt_points[i],*joint_space,*ee_space));
        };
        sample_tree = shared_ptr< ee_bsp_tree_type >(new ee_bsp_tree_type(
          indices.begin(),
          indices.end(),
          ee_space, 
          sample_points));
      } else
        sample_tree = shared_ptr< ee_bsp_tree_type >();
    };
    
    manip_ik_knn_starts_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
                            const CLIKCalcFactory& aCLIKCalcFactory = CLIKCalcFactory(),
                            const shared_ptr< JointSpaceType >& aJointSpace = shared_ptr< JointSpaceType >(),
                            const shared_ptr< EESpaceType >& aEESpace = shared_ptr< EESpaceType >(),
                            const std::size_t& aSampleCount = 1000,
                            const std::size_t& aNeighborCount = 10) :
                            model(aModel),
                            clik_calc_factory(aCLIKCalcFactory), 
                            ik_calc(),
                            joint_space(aJointSpace),
                            ee_space(aEESpace),
                            sample_count(aSampleCount),
                            neighbor_count(aNeighborCount),
                            sample_points(aSampleCount),
                            sample_jt_points(aSampleCount) {
      initialize_ik_data();
    };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      if(!ik_calc)
        initialize_ik_data(true);
            
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      ik_calc->readDesiredFromModel();
      
      std::vector< joint_state_type > neighbors; neighbors.reserve(neighbor_count);
      
      if(sample_tree) {
        std::vector< std::size_t > neighbors_ids(neighbor_count);
        std::vector< std::size_t >::iterator it_end = sample_tree->find_nearest(pt, neighbors_ids.begin(), neighbor_count);
        neighbors_ids.erase(it_end,neighbors_ids.end());
        for(std::size_t i = 0; i < neighbors_ids.size(); ++i)
          neighbors.push_back(sample_jt_points[neighbors_ids[i]]);
      } else {
        for(std::size_t i = 0; i < neighbor_count; ++i)
          neighbors.push_back(get(random_sampler,space_out)(space_out));
      };
      
      for(std::size_t i = 0; i < neighbors.size(); ++i) {
      
        detail::write_joint_coordinates_impl(neighbors[i], space_out, model);
      
        vect_n<double> x = ik_calc->readJointStatesFromModel();
        
        try {
          ik_calc->runOptimizer(x);
        } catch(singularity_error& e) { RK_UNUSED(e);
        } catch(maximum_iteration& e) { RK_UNUSED(e);
        } catch(optim::infeasible_problem& e) { RK_UNUSED(e);
        };
        
        if( norm_2( kte::manip_clik_calculator::eq_function(ik_calc.get())(x) ) < ik_calc->optimizer.tol * 10.0 ) {
          ik_calc->writeJointStatesToModel(x);
          break;
        } else if( i == neighbors.size() - 1 )
          throw optim::infeasible_problem("The inverse kinematics problem cannot be solved!");
      };
      
      detail::read_joint_coordinates_impl(result,space_out,model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_SAVE_WITH_NAME(joint_space)
        & RK_SERIAL_SAVE_WITH_NAME(ee_space)
        & RK_SERIAL_SAVE_WITH_NAME(sample_count)
        & RK_SERIAL_SAVE_WITH_NAME(neighbor_count);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_LOAD_WITH_NAME(joint_space)
        & RK_SERIAL_LOAD_WITH_NAME(ee_space)
        & RK_SERIAL_LOAD_WITH_NAME(sample_count)
        & RK_SERIAL_LOAD_WITH_NAME(neighbor_count);
      
      sample_points = boost::vector_property_map< ee_state_type >(sample_count);
      sample_jt_points = std::vector< joint_state_type >(sample_count);
      
      initialize_ik_data();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001E,1,"manip_ik_knn_starts_map",shared_object)
    
};

#if 0
/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam CLIKCalcFactory The factory type for creating a CLIK calculator.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename CLIKCalcFactory, typename RateLimitMap, typename JointSpaceType, typename EESpaceType>
class manip_rl_ik_knn_starts_map : public shared_object {
  public:
    
    typedef manip_rl_ik_knn_starts_map<CLIKCalcFactory,RateLimitMap,JointSpaceType,EESpaceType> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model;
    /** This holds the inverse kinematics calculator factory. */
    CLIKCalcFactory clik_calc_factory;
    mutable shared_ptr< kte::manip_clik_calculator > ik_calc;
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    shared_ptr< JointSpaceType > joint_space;
    shared_ptr< EESpaceType > ee_space;
    
    std::size_t sample_count;
    std::size_t neighbor_count;
    
    typedef typename topology_traits< EESpaceType >::point_type ee_state_type;
    typedef typename topology_traits< JointSpaceType >::point_type joint_state_type;
    
    mutable boost::vector_property_map< ee_state_type > sample_points;
    mutable std::vector< joint_state_type > sample_jt_points;
    
    typedef dvp_tree< std::size_t, 
                      EESpaceType,
                      boost::vector_property_map< ee_state_type > > ee_bsp_tree_type;
    
    mutable shared_ptr< ee_bsp_tree_type > sample_tree;
    
    void initialize_ik_data(bool force_init = false) const {
      if((force_init) || ((joint_limits_map) && (model)))
        ik_calc = clik_calc_factory.create_calculator(*joint_limits_map,model);
      else
        ik_calc = shared_ptr< kte::manip_clik_calculator >();
      
      if((joint_space) && (ee_space)) {
        std::vector< std::size_t > indices(sample_count);
        manip_rl_direct_kin_map<RateLimitMap> dk_map = manip_rl_direct_kin_map<RateLimitMap>(joint_limits_map,model);
        for(std::size_t i = 0; i < sample_count; ++i) {
          indices[i] = i;
          sample_jt_points[i] = get(random_sampler,*joint_space)(*joint_space);
          put(sample_points, i, dk_map.map_to_space(sample_jt_points[i],*joint_space,*ee_space));
        };
      
        sample_tree = shared_ptr< ee_bsp_tree_type >(new ee_bsp_tree_type(
          indices.begin(),
          indices.end(),
          *ee_space, 
          sample_points));
      } else
        sample_tree = shared_ptr< ee_bsp_tree_type >();
    };
    
    manip_rl_ik_knn_starts_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
                               const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >(),
                               const CLIKCalcFactory& aCLIKCalcFactory = CLIKCalcFactory(),
                               const shared_ptr< JointSpaceType >& aJointSpace = shared_ptr< JointSpaceType >(),
                               const shared_ptr< EESpaceType >& aEESpace = shared_ptr< EESpaceType >(),
                               const std::size_t& aSampleCount = 1000,
                               const std::size_t& aNeighborCount = 10) :
                               model(aModel),
                               clik_calc_factory(aCLIKCalcFactory), ik_calc(),
                               joint_limits_map(aJointLimitMap),
                               joint_space(aJointSpace),
                               ee_space(aEESpace),
                               sample_count(aSampleCount),
                               neighbor_count(aNeighborCount),
                               sample_points(aSampleCount),
                               sample_jt_points(aSampleCount) {
      initialize_ik_data();
    };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      typedef typename RateLimitMap::normal_space_type NormalJointSpace;
      
      if(!ik_calc)
        initialize_ik_data(true);
      
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      ik_calc->readDesiredFromModel();
      
      std::vector< joint_state_type > neighbors; neighbors.reserve(neighbor_count);
      
      if(sample_tree) {
        std::vector< std::size_t > neighbors_ids(neighbor_count);
        std::vector< std::size_t >::iterator it_end = sample_tree->find_nearest(pt, neighbors_ids.begin(), neighbor_count);
        neighbors_ids.erase(it_end,neighbors_ids.end());
        for(std::size_t i = 0; i < neighbors_ids.size(); ++i)
          neighbors.push_back(sample_jt_points[neighbors_ids[i]]);
      } else {
        for(std::size_t i = 0; i < neighbor_count; ++i)
          neighbors.push_back(get(random_sampler,space_out)(space_out));
      };
      
      for(std::size_t i = 0; i < neighbors.size(); ++i) {
        
        typename topology_traits<NormalJointSpace>::point_type ip_inter = 
          joint_limits_map->map_to_space(neighbors[i], space_out, NormalJointSpace());
        detail::write_joint_coordinates_impl(ip_inter, NormalJointSpace(), model);
        
        vect_n<double> x = ik_calc->readJointStatesFromModel();
        
        try {
          ik_calc->runOptimizer(x);
        } catch(singularity_error& e) { 
        } catch(maximum_iteration& e) {
        } catch(optim::infeasible_problem& e) {
        };
        
        if( norm_2( kte::manip_clik_calculator::eq_function(ik_calc.get())(x) ) < ik_calc->optimizer.tol * 10.0 ) {
          ik_calc->writeJointStatesToModel(x);
          break;
        } else if( i == neighbors.size() - 1 )
          throw optim::infeasible_problem("The inverse kinematics problem cannot be solved!");
      };
      
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl(result_inter, NormalJointSpace(), model);
      result = joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_SAVE_WITH_NAME(joint_limits_map)
        & RK_SERIAL_SAVE_WITH_NAME(joint_space)
        & RK_SERIAL_SAVE_WITH_NAME(ee_space)
        & RK_SERIAL_SAVE_WITH_NAME(sample_count)
        & RK_SERIAL_SAVE_WITH_NAME(neighbor_count);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_LOAD_WITH_NAME(joint_limits_map)
        & RK_SERIAL_LOAD_WITH_NAME(joint_space)
        & RK_SERIAL_LOAD_WITH_NAME(ee_space)
        & RK_SERIAL_LOAD_WITH_NAME(sample_count)
        & RK_SERIAL_LOAD_WITH_NAME(neighbor_count);
      initialize_ik_data();
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC240001F,1,"manip_rl_ik_knn_starts_map",shared_object)
    
    
};
#endif
#endif

#endif

};


};

#endif








