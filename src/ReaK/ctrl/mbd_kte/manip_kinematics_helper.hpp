/**
 * \file manip_kinematics_helper.hpp
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

#ifndef REAK_MANIP_KINEMATICS_HELPER_HPP
#define REAK_MANIP_KINEMATICS_HELPER_HPP

#include "manip_kinematics_model.hpp"

namespace ReaK {

namespace kte {


class manip_kin_mdl_joint_io {
  private:
    const direct_kinematics_model* model;
    
  public:
    manip_kin_mdl_joint_io(const direct_kinematics_model* aModel) : model(aModel) { };
    ~manip_kin_mdl_joint_io() { };
    
    void getJointPositions(double* result) const;
    
    void setJointPositions(const double* aJointPositions);
    
    void getJointVelocities(double* result) const;
    
    void setJointVelocities(const double* aJointVelocities);
    
    void getJointAccelerations(double* result) const;
    
    void setJointAccelerations(const double* aJointAccelerations);
    
    void getDependentPositions(double* result) const;
    
    void getDependentVelocities(double* result) const;
    
    void getDependentAccelerations(double* result) const;
    
};



/**
 * This class is a helper of the direct_kinematics_model class which is used to fill in 
 * the Jacobian matrix (and its time-derivative). This class is useful because it is a friend 
 * of the direct_kinematics_model class (thus, has access to its data members), but also 
 * provides member function templates for filling in the matrices, meaning it can be used to 
 * fill in any kind of matrix type (e.g. enabling the filling of matrix sub-blocks for example).
 */
class manip_kin_mdl_jac_calculator {
  private:
    const direct_kinematics_model* model;
    
    void getJacobianMatrixAndDerivativeImpl(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>* JacDot) const;
    
  public:
    
    /**
     * Default constructor.
     */
    manip_kin_mdl_jac_calculator(const direct_kinematics_model* aModel) : model(aModel) { };
    
    /**
     * Default destructor.
     */
    ~manip_kin_mdl_jac_calculator() { };
    
    /**
     * Get the Jacobian matrix for the system (or twist-shaping matrix). The Jacobian takes the velocity 
     * information of the system coordinates and frames, and maps them to velocity information 
     * of the system's dependent coordinates and frames.
     * \param Jac stores, as output, the calculated system's Jacobian matrix.
     */
    void getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const {
      getJacobianMatrixAndDerivativeImpl(Jac, static_cast<mat<double,mat_structure::rectangular>*>(NULL));
    };
    
    /**
     * Get the Jacobian matrix for the system (or twist-shaping matrix), and its time-derivative. 
     * The Jacobian takes the velocity information of the system coordinates and frames, and maps 
     * them to velocity information of the system's dependent coordinates and frames. The time-derivative
     * of the Jacobian matrix will map the velocity information of the system coordinates and frames 
     * to the acceleration information of the system's dependent coordinates and frames.
     * \param Jac stores, as output, the calculated system's Jacobian matrix.
     * \param JacDot stores, as output, the calculated time-derivative of the system's Jacobian matrix.
     */
    void getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, mat<double,mat_structure::rectangular>& JacDot) const {
      getJacobianMatrixAndDerivativeImpl(Jac,&JacDot);
    };
    
};





};

};

#endif











