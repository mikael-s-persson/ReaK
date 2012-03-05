/**
 * \file manipulator_model_helper.hpp
 *
 * This library declares helper classes to use manipulator models, both kinematic only or 
 * dynamic as well. Essentially, the model of the manipulator is only a KTE chain provided 
 * by the user, but these manipulator-model classes take care of grouping the joints, their
 * limits, and their jacobian matrices.
 * The helper classes provided here are friends of the manipulator model classes and implement
 * some details of their functionalities. The helper classes can also be used elsewhere to access 
 * lower-level functions.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
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

#ifndef REAK_MANIPULATOR_MODEL_HELPER_HPP
#define REAK_MANIPULATOR_MODEL_HELPER_HPP

#include "manipulator_model.hpp"

namespace ReaK {

namespace kte {


/**
 */
class manip_kin_mdl_jac_calculator {
  public:

    /**
     * Default constructor.
     */
    manip_kin_mdl_jac_calculator(const manipulator_kinematics_model* aModel) : model(aModel) { };
    
    /**
     * Default destructor.
     */
    virtual ~manip_kin_mdl_jac_calculator() { };
    
    /**
     * Get the Jacobian matrix for the system (or twist-shaping matrix). The Jacobian takes the velocity 
     * information of the system coordinates and frames, and maps them to velocity information 
     * of the system's dependent coordinates and frames.
     * \param Jac stores, as output, the calculated system's Jacobian matrix.
     */
    template <typename Matrix1>
    void getJacobianMatrix(Matrix1& Jac) const {
      getJacobianMatrixAndDerivativeImpl(Jac,static_cast<mat<double,mat_structure::rectangular>*>(NULL));
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
    template <typename Matrix1, typename Matrix2>
    void getJacobianMatrixAndDerivative(Matrix1& Jac, Matrix2& JacDot) const {
      getJacobianMatrixAndDerivativeImpl(Jac,&JacDot);
    };
    
  private:
    const manipulator_kinematics_model* model;
    
    template <typename Matrix1, typename Matrix2>
    void getJacobianMatrixAndDerivativeImpl(Matrix1& Jac, Matrix2* JacDot) const {
      unsigned int m = model->getDependentVelocitiesCount();
      unsigned int n = model->getJointVelocitiesCount();
      Jac    = mat<double,mat_structure::nil>(m,n);
      if(JacDot)
        *JacDot = mat<double,mat_structure::nil>(m,n);

      unsigned int RowInd = 0;
  
    /****************************************************************************************
     *                             Gen Coords 
     * *************************************************************************************/
  
  
      for(unsigned int i=0; i < model->mCoords.size(); ++i) {
        RowInd = 0;
    
        for(unsigned int j=0; j < model->mDependentGenCoords.size(); ++j) {
          if(model->mDependentGenCoords[j]->mUpStreamJoints.find(model->mCoords[i]) != model->mDependentGenCoords[j]->mUpStreamJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd),range(i,i));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd),range(i,i));
	      model->mDependentGenCoords[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->write_to_matrices(subJac, subJacDot);
	    } else {
	      model->mDependentGenCoords[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->write_to_matrices(subJac);
	    };
          };
          RowInd++;
        };
    
        for(unsigned int j=0; j < model->mDependent2DFrames.size(); ++j) {
          if(model->mDependent2DFrames[j]->mUpStreamJoints.find(model->mCoords[i]) != model->mDependent2DFrames[j]->mUpStreamJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+2),range(i,i));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+2),range(i,i));
	      model->mDependent2DFrames[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac, subJacDot);
            } else {
	      model->mDependent2DFrames[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 3;
        };
    
        for(unsigned int j=0; j < model->mDependent3DFrames.size(); ++j) {
          if(model->mDependent3DFrames[j]->mUpStreamJoints.find(model->mCoords[i]) != model->mDependent3DFrames[j]->mUpStreamJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+5),range(i,i));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+5),range(i,i));
  	      model->mDependent3DFrames[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependent3DFrames[j]
	        ->mUpStreamJoints[model->mCoords[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 6;
        };
      };

  
    /****************************************************************************************
     *                             2D Frames 
     * *************************************************************************************/
  
      unsigned int base_i = model->mCoords.size();
      for(unsigned int i=0; i < model->mFrames2D.size(); ++i) {
        RowInd = 0;

        for(unsigned int j=0; j < model->mDependentGenCoords.size(); ++j) {
          if(model->mDependentGenCoords[j]->mUpStream2DJoints.find(model->mFrames2D[i]) != model->mDependentGenCoords[j]->mUpStream2DJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd),range(3 * i + base_i, 3 * i + base_i + 2));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd),range(3 * i + base_i, 3 * i + base_i + 2));
	      model->mDependentGenCoords[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependentGenCoords[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->write_to_matrices(subJac);
	    };
          };
          RowInd++;
        };

        for(unsigned int j=0; j < model->mDependent2DFrames.size(); ++j) {
          if(model->mDependent2DFrames[j]->mUpStream2DJoints.find(model->mFrames2D[i]) != model->mDependent2DFrames[j]->mUpStream2DJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+2),range(3 * i + base_i, 3 * i + base_i + 2));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+2),range(3 * i + base_i, 3 * i + base_i + 2));
  	      model->mDependent2DFrames[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependent2DFrames[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 3;
        };

        for(unsigned int j=0; j < model->mDependent3DFrames.size(); ++j) {
          if(model->mDependent3DFrames[j]->mUpStream2DJoints.find(model->mFrames2D[i]) != model->mDependent3DFrames[j]->mUpStream2DJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+5),range(3 * i + base_i, 3 * i + base_i + 2));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+5),range(3 * i + base_i, 3 * i + base_i + 2));
	      model->mDependent3DFrames[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependent3DFrames[j]
	        ->mUpStream2DJoints[model->mFrames2D[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 6;
        };
      };

  
  
    /****************************************************************************************
     *                             3D Frames 
     * *************************************************************************************/
  
      base_i = model->mCoords.size() + 3 * model->mFrames2D.size();
      for(unsigned int i=0; i < model->mFrames3D.size(); ++i) {
        RowInd = 0; 

        for(unsigned int j=0; j < model->mDependentGenCoords.size(); ++j) {
          if(model->mDependentGenCoords[j]->mUpStreamJoints.find(model->mCoords[i]) != model->mDependentGenCoords[j]->mUpStreamJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd),range(6 * i + base_i, 6 * i + base_i + 5));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd),range(6 * i + base_i, 6 * i + base_i + 5));
	      model->mDependentGenCoords[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependentGenCoords[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->write_to_matrices(subJac);
	    };
          };
          RowInd++;
        };

        for(unsigned int j=0; j < model->mDependent2DFrames.size(); ++j) {
          if(model->mDependent2DFrames[j]->mUpStream3DJoints.find(model->mFrames3D[i]) != model->mDependent2DFrames[j]->mUpStream3DJoints.end()) {
	    mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+2),range(6 * i + base_i, 6 * i + base_i + 5));
	    if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+2),range(6 * i + base_i, 6 * i + base_i + 5));
	      model->mDependent2DFrames[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependent2DFrames[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->get_jac_relative_to(model->mDependent2DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 3;
        };

        for(unsigned int j=0; j < model->mDependent3DFrames.size(); ++j) {
          if(model->mDependent3DFrames[j]->mUpStream3DJoints.find(model->mFrames3D[i]) != model->mDependent3DFrames[j]->mUpStream3DJoints.end()) {
            mat_sub_block< Matrix1 > subJac = sub(Jac)(range(RowInd,RowInd+5),range(6 * i + base_i, 6 * i + base_i + 5));
            if(JacDot) {
	      mat_sub_block< Matrix2 > subJacDot = sub(*JacDot)(range(RowInd,RowInd+5),range(6 * i + base_i, 6 * i + base_i + 5));
	      model->mDependent3DFrames[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac,subJacDot);
	    } else {
	      model->mDependent3DFrames[j]
	        ->mUpStream3DJoints[model->mFrames3D[i]]
	          ->get_jac_relative_to(model->mDependent3DFrames[j]->mFrame)
	            .write_to_matrices(subJac);
	    };
          };
          RowInd += 6;
        };
      };

    };


    
    
  
};







};

};

#endif











