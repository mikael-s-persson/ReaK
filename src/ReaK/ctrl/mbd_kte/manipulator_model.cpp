
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

#include "manipulator_model.hpp"

#include "lin_alg/mat_num_exceptions.hpp"
#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace kte {


  
manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< gen_coord<double> >& aCoord) {
  if(aCoord)
    mCoords.push_back(aCoord);
  return *this;
};

manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< frame_2D<double> >& aFrame2D) {
  if(aFrame2D)
    mFrames2D.push_back(aFrame2D);
  return *this;
};

manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< frame_3D<double> >& aFrame3D) {
  if(aFrame3D)
    mFrames3D.push_back(aFrame3D);
  return *this;
};

manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< joint_dependent_gen_coord >& aDependentGenCoord) {
  if(aDependentGenCoord)
    mDependentGenCoords.push_back(aDependentGenCoord);
  return *this;
};

manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< joint_dependent_frame_2D >& aDependent2DFrame) {
  if(aDependent2DFrame)
    mDependent2DFrames.push_back(aDependent2DFrame);
  return *this;
};

manipulator_kinematics_model& manipulator_kinematics_model::operator <<(const shared_ptr< joint_dependent_frame_3D >& aDependent3DFrame) {
  if(aDependent3DFrame)
    mDependent3DFrames.push_back(aDependent3DFrame);
  return *this;
};


void manipulator_kinematics_model::getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) {
  getJacobianMatrixAndDerivativeImpl(&Jac,NULL);
};

void manipulator_kinematics_model::getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, 
								  mat<double,mat_structure::rectangular>& JacDot) {
  getJacobianMatrixAndDerivativeImpl(&Jac,&JacDot);
};
  
void manipulator_kinematics_model::getJacobianMatrixAndDerivativeImpl(mat<double,mat_structure::rectangular>* Jac, 
								      mat<double,mat_structure::rectangular>* JacDot) {
  unsigned int m = getDependentVelocitiesCount();
  unsigned int n = getJointVelocitiesCount();
  *Jac    = mat<double,mat_structure::nil>(m,n);
  if(JacDot)
    *JacDot = mat<double,mat_structure::nil>(m,n);

  unsigned int RowInd = 0;
  
/****************************************************************************************
 *                             Gen Coords 
 * *************************************************************************************/
  
  
  for(unsigned int i=0; i<mCoords.size(); ++i) {
    RowInd = 0;
    
    for(unsigned int j=0; j<mDependentGenCoords.size(); ++j) {
      if(mDependentGenCoords[j]->mUpStreamJoints.find(mCoords[i]) != mDependentGenCoords[j]->mUpStreamJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,1,1,RowInd,i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,1,1,RowInd,i);
	  mDependentGenCoords[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->write_to_matrices(subJac,subJacDot);
	} else {
	  mDependentGenCoords[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->write_to_matrices(subJac);
	};
      };
      RowInd++;
    };
    
    for(unsigned int j=0; j<mDependent2DFrames.size(); ++j) {
      if(mDependent2DFrames[j]->mUpStreamJoints.find(mCoords[i]) != mDependent2DFrames[j]->mUpStreamJoints.end()) {

	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,3,1,RowInd,i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,3,1,RowInd,i);
	  mDependent2DFrames[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent2DFrames[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 3;
    };
    
    for(unsigned int j=0; j<mDependent3DFrames.size(); ++j) {
      if(mDependent3DFrames[j]->mUpStreamJoints.find(mCoords[i]) != mDependent3DFrames[j]->mUpStreamJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,6,1,RowInd,i);
	if(JacDot) {
  	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,6,1,RowInd,i);
	  mDependent3DFrames[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent3DFrames[j]
	    ->mUpStreamJoints[mCoords[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 6;
    };
  };

  
/****************************************************************************************
 *                             2D Frames 
 * *************************************************************************************/
  
  unsigned int base_i = mCoords.size();
  for(unsigned int i=0; i<mFrames2D.size(); ++i) {
    RowInd = 0;

    for(unsigned int j=0; j<mDependentGenCoords.size(); ++j) {
      if(mDependentGenCoords[j]->mUpStream2DJoints.find(mFrames2D[i]) != mDependentGenCoords[j]->mUpStream2DJoints.end()) {
		
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,1,3,RowInd,3*i+base_i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,1,3,RowInd,3*i+base_i);
	  mDependentGenCoords[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->write_to_matrices(subJac,subJacDot);
	} else {
	  mDependentGenCoords[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->write_to_matrices(subJac);
	};
      };
      RowInd++;
    };

    for(unsigned int j=0; j<mDependent2DFrames.size(); ++j) {
      if(mDependent2DFrames[j]->mUpStream2DJoints.find(mFrames2D[i]) != mDependent2DFrames[j]->mUpStream2DJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,3,3,RowInd,3*i+base_i);
	if(JacDot) {
  	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,3,3,RowInd,3*i+base_i);
	  mDependent2DFrames[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent2DFrames[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<mDependent3DFrames.size(); ++j) {
      if(mDependent3DFrames[j]->mUpStream2DJoints.find(mFrames2D[i]) != mDependent3DFrames[j]->mUpStream2DJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,6,3,RowInd,3*i+base_i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,6,3,RowInd,3*i+base_i);
	  mDependent3DFrames[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent3DFrames[j]
	    ->mUpStream2DJoints[mFrames2D[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 6;
    };
  };

  
  
/****************************************************************************************
 *                             3D Frames 
 * *************************************************************************************/
  
  base_i = mCoords.size() + 3*mFrames2D.size();
  for(unsigned int i=0; i<mFrames3D.size(); ++i) {
    RowInd = 0; 

    for(unsigned int j=0; j<mDependentGenCoords.size(); ++j) {
      if(mDependentGenCoords[j]->mUpStreamJoints.find(mCoords[i]) != mDependentGenCoords[j]->mUpStreamJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,1,6,RowInd,6*i+base_i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,1,6,RowInd,6*i+base_i);
	  mDependentGenCoords[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->write_to_matrices(subJac,subJacDot);
	} else {
	  mDependentGenCoords[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->write_to_matrices(subJac);
	};
      };
      RowInd++;
    };

    for(unsigned int j=0; j<mDependent2DFrames.size(); ++j) {
      if(mDependent2DFrames[j]->mUpStream3DJoints.find(mFrames3D[i]) != mDependent2DFrames[j]->mUpStream3DJoints.end()) {
	
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,3,6,RowInd,6*i+base_i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,3,6,RowInd,6*i+base_i);
	  mDependent2DFrames[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent2DFrames[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->get_jac_relative_to(mDependent2DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<mDependent3DFrames.size(); ++j) {
      if(mDependent3DFrames[j]->mUpStream3DJoints.find(mFrames3D[i]) != mDependent3DFrames[j]->mUpStream3DJoints.end()) {
		
	mat_sub_block< mat<double,mat_structure::rectangular> > subJac(*Jac,6,6,RowInd,6*i+base_i);
	if(JacDot) {
	  mat_sub_block< mat<double,mat_structure::rectangular> > subJacDot(*JacDot,6,6,RowInd,6*i+base_i);
	  mDependent3DFrames[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac,subJacDot);
	} else {
	  mDependent3DFrames[j]
	    ->mUpStream3DJoints[mFrames3D[i]]
	      ->get_jac_relative_to(mDependent3DFrames[j]->mFrame)
	        .write_to_matrices(subJac);
	};
      };
      RowInd += 6;
    };
  };

};








vect_n<double> manipulator_kinematics_model::getJointPositions() const {
  vect_n<double> result(getJointPositionsCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    result[j] = (*it)->q;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    result[j] = (*it)->Position[0]; ++j;
    result[j] = (*it)->Position[1]; ++j;
    result[j] = (*it)->Rotation[0]; ++j;
    result[j] = (*it)->Rotation[1]; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    result[j] = (*it)->Position[0]; ++j;
    result[j] = (*it)->Position[1]; ++j;
    result[j] = (*it)->Position[2]; ++j;
    result[j] = (*it)->Quat[0]; ++j;
    result[j] = (*it)->Quat[1]; ++j;
    result[j] = (*it)->Quat[2]; ++j;
    result[j] = (*it)->Quat[3]; ++j;
  };
  
  return result;
};
    
void manipulator_kinematics_model::setJointPositions(const vect_n<double>& aJointPositions) {
  if(aJointPositions.size() != getJointPositionsCount())
    throw std::range_error("Joint-position vector has incorrect dimensions!");
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    (*it)->q = aJointPositions[j];

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    (*it)->Position[0] = aJointPositions[j]; ++j;
    (*it)->Position[1] = aJointPositions[j]; ++j;
    (*it)->Rotation = rot_mat_2D<double>(vect<double,2>(aJointPositions[j],aJointPositions[j+1])); j += 2;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    (*it)->Position[0] = aJointPositions[j]; ++j;
    (*it)->Position[1] = aJointPositions[j]; ++j;
    (*it)->Position[2] = aJointPositions[j]; ++j;
    (*it)->Quat = quaternion<double>(vect<double,4>(aJointPositions[j],
                                                    aJointPositions[j+1],
						    aJointPositions[j+2],
						    aJointPositions[j+3])); j += 4;
  };
};
    
vect_n<double> manipulator_kinematics_model::getJointVelocities() const {
  vect_n<double> result(getJointVelocitiesCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    result[j] = (*it)->q_dot;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    result[j] = (*it)->Velocity[0]; ++j;
    result[j] = (*it)->Velocity[1]; ++j;
    result[j] = (*it)->AngVelocity; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    result[j] = (*it)->Velocity[0]; ++j;
    result[j] = (*it)->Velocity[1]; ++j;
    result[j] = (*it)->Velocity[2]; ++j;
    result[j] = (*it)->AngVelocity[0]; ++j;
    result[j] = (*it)->AngVelocity[1]; ++j;
    result[j] = (*it)->AngVelocity[2]; ++j;
  };
  
  return result;  
};
    
void manipulator_kinematics_model::setJointVelocities(const vect_n<double>& aJointVelocities) {
  if(aJointVelocities.size() != getJointVelocitiesCount())
    throw std::range_error("Joint-velocity vector has incorrect dimensions!");
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    (*it)->q_dot = aJointVelocities[j];

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    (*it)->Velocity[0] = aJointVelocities[j]; ++j;
    (*it)->Velocity[1] = aJointVelocities[j]; ++j;
    (*it)->AngVelocity = aJointVelocities[j]; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    (*it)->Velocity[0] = aJointVelocities[j]; ++j;
    (*it)->Velocity[1] = aJointVelocities[j]; ++j;
    (*it)->Velocity[2] = aJointVelocities[j]; ++j;
    (*it)->AngVelocity[0] = aJointVelocities[j]; ++j;
    (*it)->AngVelocity[1] = aJointVelocities[j]; ++j;
    (*it)->AngVelocity[2] = aJointVelocities[j]; ++j;
  };
};
    
vect_n<double> manipulator_kinematics_model::getJointAccelerations() const {
  vect_n<double> result(getJointAccelerationsCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    result[j] = (*it)->q_ddot;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    result[j] = (*it)->Acceleration[0]; ++j;
    result[j] = (*it)->Acceleration[1]; ++j;
    result[j] = (*it)->AngAcceleration; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    result[j] = (*it)->Acceleration[0]; ++j;
    result[j] = (*it)->Acceleration[1]; ++j;
    result[j] = (*it)->Acceleration[2]; ++j;
    result[j] = (*it)->AngAcceleration[0]; ++j;
    result[j] = (*it)->AngAcceleration[1]; ++j;
    result[j] = (*it)->AngAcceleration[2]; ++j;
  };
  
  return result;  
};
    
void manipulator_kinematics_model::setJointAccelerations(const vect_n<double>& aJointAccelerations) {
  if(aJointAccelerations.size() != getJointAccelerationsCount())
    throw std::range_error("Joint-acceleration vector has incorrect dimensions!");
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    (*it)->q_ddot = aJointAccelerations[j];

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    (*it)->Acceleration[0] = aJointAccelerations[j]; ++j;
    (*it)->Acceleration[1] = aJointAccelerations[j]; ++j;
    (*it)->AngAcceleration = aJointAccelerations[j]; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    (*it)->Acceleration[0] = aJointAccelerations[j]; ++j;
    (*it)->Acceleration[1] = aJointAccelerations[j]; ++j;
    (*it)->Acceleration[2] = aJointAccelerations[j]; ++j;
    (*it)->AngAcceleration[0] = aJointAccelerations[j]; ++j;
    (*it)->AngAcceleration[1] = aJointAccelerations[j]; ++j;
    (*it)->AngAcceleration[2] = aJointAccelerations[j]; ++j;
  };
  
};
    
vect_n<double> manipulator_kinematics_model::getDependentPositions() const {
  vect_n<double> result(getDependentPositionsCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = mDependentGenCoords.begin(); 
      it < mDependentGenCoords.end(); ++it, ++j)
    result[j] = (*it)->mFrame->q;

  for(std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = mDependent2DFrames.begin(); 
      it < mDependent2DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Position[0]; ++j;
    result[j] = (*it)->mFrame->Position[1]; ++j;
    result[j] = (*it)->mFrame->Rotation[0]; ++j;
    result[j] = (*it)->mFrame->Rotation[1]; ++j;
  };

  for(std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = mDependent3DFrames.begin(); 
      it < mDependent3DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Position[0]; ++j;
    result[j] = (*it)->mFrame->Position[1]; ++j;
    result[j] = (*it)->mFrame->Position[2]; ++j;
    result[j] = (*it)->mFrame->Quat[0]; ++j;
    result[j] = (*it)->mFrame->Quat[1]; ++j;
    result[j] = (*it)->mFrame->Quat[2]; ++j;
    result[j] = (*it)->mFrame->Quat[3]; ++j;
  };
  
  return result;
};
    
vect_n<double> manipulator_kinematics_model::getDependentVelocities() const {
  vect_n<double> result(getDependentVelocitiesCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = mDependentGenCoords.begin(); 
      it < mDependentGenCoords.end(); ++it, ++j)
    result[j] = (*it)->mFrame->q_dot;

  for(std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = mDependent2DFrames.begin(); 
      it < mDependent2DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Velocity[0]; ++j;
    result[j] = (*it)->mFrame->Velocity[1]; ++j;
    result[j] = (*it)->mFrame->AngVelocity; ++j;
  };

  for(std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = mDependent3DFrames.begin(); 
      it < mDependent3DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Velocity[0]; ++j;
    result[j] = (*it)->mFrame->Velocity[1]; ++j;
    result[j] = (*it)->mFrame->Velocity[2]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[0]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[1]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[2]; ++j;
  };
  
  return result;  
};
    
vect_n<double> manipulator_kinematics_model::getDependentAccelerations() const {
  vect_n<double> result(getDependentAccelerationsCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = mDependentGenCoords.begin(); 
      it < mDependentGenCoords.end(); ++it, ++j)
    result[j] = (*it)->mFrame->q_ddot;

  for(std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = mDependent2DFrames.begin(); 
      it < mDependent2DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Acceleration[0]; ++j;
    result[j] = (*it)->mFrame->Acceleration[1]; ++j;
    result[j] = (*it)->mFrame->AngAcceleration; ++j;
  };

  for(std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = mDependent3DFrames.begin(); 
      it < mDependent3DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Acceleration[0]; ++j;
    result[j] = (*it)->mFrame->Acceleration[1]; ++j;
    result[j] = (*it)->mFrame->Acceleration[2]; ++j;
    result[j] = (*it)->mFrame->AngAcceleration[0]; ++j;
    result[j] = (*it)->mFrame->AngAcceleration[1]; ++j;
    result[j] = (*it)->mFrame->AngAcceleration[2]; ++j;
  };
  
  return result;  
};










  
  
  
manipulator_kinematics_model& manipulator_dynamics_model::operator <<(const shared_ptr< gen_coord<double> >& aCoord) {
  if(aCoord) {
    mMassCalc << aCoord;
    manipulator_kinematics_model::operator<<(aCoord);
  };
  return *this;
};

manipulator_kinematics_model& manipulator_dynamics_model::operator <<(const shared_ptr< frame_2D<double> >& aFrame2D) {
  if(aFrame2D) {
    mMassCalc << aFrame2D;
    manipulator_kinematics_model::operator<<(aFrame2D);
  };
  return *this;
};

manipulator_kinematics_model& manipulator_dynamics_model::operator <<(const shared_ptr< frame_3D<double> >& aFrame3D) {
  if(aFrame3D) {
    mMassCalc << aFrame3D;
    manipulator_kinematics_model::operator<<(aFrame3D);
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator <<(const shared_ptr< inertia_gen >& aInertiaGen) {
  if(aInertiaGen) {
    mMassCalc << aInertiaGen;
    *this << aInertiaGen->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator <<(const shared_ptr< inertia_2D >& aInertia2D) {
  if(aInertia2D) {
    mMassCalc << aInertia2D;
    *this << aInertia2D->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator <<(const shared_ptr< inertia_3D >& aInertia3D) {
  if(aInertia3D) {
    mMassCalc << aInertia3D;
    *this << aInertia3D->CenterOfMass();
  };
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator <<(const shared_ptr< system_input >& aInput) {
  mInputs.push_back(aInput);
  return *this;
};

manipulator_dynamics_model& manipulator_dynamics_model::operator <<(const shared_ptr< system_output >& aOutput) {
  mOutputs.push_back(aOutput);
  return *this;
};
    




vect_n<double> manipulator_dynamics_model::getJointStates() const {
  vect_n<double> result(getJointStatesCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    result[j] = (*it)->q;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    result[j] = (*it)->Position[0]; ++j;
    result[j] = (*it)->Position[1]; ++j;
    result[j] = (*it)->Rotation[0]; ++j;
    result[j] = (*it)->Rotation[1]; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    result[j] = (*it)->Position[0]; ++j;
    result[j] = (*it)->Position[1]; ++j;
    result[j] = (*it)->Position[2]; ++j;
    result[j] = (*it)->Quat[0]; ++j;
    result[j] = (*it)->Quat[1]; ++j;
    result[j] = (*it)->Quat[2]; ++j;
    result[j] = (*it)->Quat[3]; ++j;
  };
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    result[j] = (*it)->q_dot;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    result[j] = (*it)->Velocity[0]; ++j;
    result[j] = (*it)->Velocity[1]; ++j;
    result[j] = (*it)->AngVelocity; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    result[j] = (*it)->Velocity[0]; ++j;
    result[j] = (*it)->Velocity[1]; ++j;
    result[j] = (*it)->Velocity[2]; ++j;
    result[j] = (*it)->AngVelocity[0]; ++j;
    result[j] = (*it)->AngVelocity[1]; ++j;
    result[j] = (*it)->AngVelocity[2]; ++j;
  };
  
  return result;
};
    
void manipulator_dynamics_model::setJointStates(const vect_n<double>& aJointStates) {
  if(aJointStates.size() != getJointStatesCount())
    throw std::range_error("Joint-state vector has incorrect dimensions!");
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    (*it)->q = aJointStates[j];

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    (*it)->Position[0] = aJointStates[j]; ++j;
    (*it)->Position[1] = aJointStates[j]; ++j;
    (*it)->Rotation = rot_mat_2D<double>(vect<double,2>(aJointStates[j],aJointStates[j+1])); j += 2;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    (*it)->Position[0] = aJointStates[j]; ++j;
    (*it)->Position[1] = aJointStates[j]; ++j;
    (*it)->Position[2] = aJointStates[j]; ++j;
    (*it)->Quat = quaternion<double>(vect<double,4>(aJointStates[j],
                                                    aJointStates[j+1],
						    aJointStates[j+2],
						    aJointStates[j+3])); j += 4;
  };
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    (*it)->q_dot = aJointStates[j];

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    (*it)->Velocity[0] = aJointStates[j]; ++j;
    (*it)->Velocity[1] = aJointStates[j]; ++j;
    (*it)->AngVelocity = aJointStates[j]; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    (*it)->Velocity[0] = aJointStates[j]; ++j;
    (*it)->Velocity[1] = aJointStates[j]; ++j;
    (*it)->Velocity[2] = aJointStates[j]; ++j;
    (*it)->AngVelocity[0] = aJointStates[j]; ++j;
    (*it)->AngVelocity[1] = aJointStates[j]; ++j;
    (*it)->AngVelocity[2] = aJointStates[j]; ++j;
  };
};


void RK_CALL manipulator_dynamics_model::computeOutput(double aTime,const ReaK::vect_n<double>& aState, ReaK::vect_n<double>& aOutput) {
  setJointStates(aState);
  
  doMotion();
  clearForce();
  doForce();
  
  aOutput.resize(getOutputsCount());
  
  unsigned int i = 0;
  for(std::vector< shared_ptr<system_output> >::iterator it = mOutputs.begin(); it != mOutputs.end(); ++it)
    for(unsigned int k = 0; k < (*it)->getOutputCount(); ++k)
      aOutput[i++] = (*it)->getOutput(k);
};
    
void RK_CALL manipulator_dynamics_model::setInput(const ReaK::vect_n<double>& aInput) {
  if(aInput.size() != getInputsCount())
    throw std::range_error("The size of the input-vector to the manipulator model is not correct!");
  
  unsigned int i = 0;
  for(std::vector< shared_ptr<system_input> >::iterator it = mInputs.begin(); it != mInputs.end(); ++it)
    for(unsigned int k = 0; k < (*it)->getInputCount(); ++k)
      (*it)->getInput(k) = aInput[i++];
};

vect_n<double> manipulator_dynamics_model::getInput() const {
  vect_n<double> result(getInputsCount());
  
  unsigned int i = 0;
  for(std::vector< shared_ptr<system_input> >::const_iterator it = mInputs.begin(); it != mInputs.end(); ++it)
    for(unsigned int k = 0; k < (*it)->getInputCount(); ++k)
      result[i++] = (*it)->getInput(k);
  return result;
};
    
void RK_CALL manipulator_dynamics_model::computeStateRate(double aTime,const vect_n<double>& aState, vect_n<double>& aStateRate) {
  setJointStates(aState);
  
  doMotion();
  clearForce();
  doForce();
  
  aStateRate.resize(getJointStatesCount());
  
  unsigned int j = 0;
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    aStateRate[j] = (*it)->q_dot;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    aStateRate[j] = (*it)->Velocity[0]; ++j;
    aStateRate[j] = (*it)->Velocity[1]; ++j;
    aStateRate[j] = -(*it)->Rotation[1] * (*it)->AngVelocity; ++j;
    aStateRate[j] =  (*it)->Rotation[0] * (*it)->AngVelocity; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    aStateRate[j] = (*it)->Velocity[0]; ++j;
    aStateRate[j] = (*it)->Velocity[1]; ++j;
    aStateRate[j] = (*it)->Velocity[2]; ++j;
    aStateRate[j] = (*it)->QuatDot[0]; ++j;
    aStateRate[j] = (*it)->QuatDot[1]; ++j; 
    aStateRate[j] = (*it)->QuatDot[2]; ++j; 
    aStateRate[j] = (*it)->QuatDot[3]; ++j; 
  };
  
  for(std::vector< shared_ptr< gen_coord<double> > >::const_iterator it = mCoords.begin(); 
      it < mCoords.end(); ++it, ++j)
    aStateRate[j] = (*it)->f;

  for(std::vector< shared_ptr< frame_2D<double> > >::const_iterator it = mFrames2D.begin(); 
      it < mFrames2D.end(); ++it) {
    aStateRate[j] = (*it)->Force[0]; ++j;
    aStateRate[j] = (*it)->Force[1]; ++j;
    aStateRate[j] = (*it)->Torque; ++j;
  };

  for(std::vector< shared_ptr< frame_3D<double> > >::const_iterator it = mFrames3D.begin(); 
      it < mFrames3D.end(); ++it) {
    aStateRate[j] = (*it)->Force[0]; ++j;
    aStateRate[j] = (*it)->Force[1]; ++j;
    aStateRate[j] = (*it)->Force[2]; ++j;
    aStateRate[j] = (*it)->Torque[0]; ++j;
    aStateRate[j] = (*it)->Torque[1]; ++j; 
    aStateRate[j] = (*it)->Torque[2]; ++j; 
  };
  
  mat<double,mat_structure::symmetric> Msys(getJointAccelerationsCount());
  getMassMatrix(Msys);
  try {
    mat_vect_adaptor< vect_n<double> > acc_as_mat(aStateRate, getJointAccelerationsCount(), 1, getJointPositionsCount());
    linsolve_Cholesky(Msys,acc_as_mat);
  } catch(singularity_error& e) {
    std::stringstream ss; ss << "Mass matrix is singular in the manipulator model '" << getName() << "' at time " << aTime << " seconds.";
    throw singularity_error(ss.str());
  };
};
    
vect_n<double> manipulator_dynamics_model::getDependentStates() const {
  vect_n<double> result(getDependentStatesCount());
  
  unsigned int j = 0;
  
  for(std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = mDependentGenCoords.begin(); 
      it < mDependentGenCoords.end(); ++it, ++j)
    result[j] = (*it)->mFrame->q;

  for(std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = mDependent2DFrames.begin(); 
      it < mDependent2DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Position[0]; ++j;
    result[j] = (*it)->mFrame->Position[1]; ++j;
    result[j] = (*it)->mFrame->Rotation[0]; ++j;
    result[j] = (*it)->mFrame->Rotation[1]; ++j;
  };

  for(std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = mDependent3DFrames.begin(); 
      it < mDependent3DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Position[0]; ++j;
    result[j] = (*it)->mFrame->Position[1]; ++j;
    result[j] = (*it)->mFrame->Position[2]; ++j;
    result[j] = (*it)->mFrame->Quat[0]; ++j;
    result[j] = (*it)->mFrame->Quat[1]; ++j;
    result[j] = (*it)->mFrame->Quat[2]; ++j;
    result[j] = (*it)->mFrame->Quat[3]; ++j;
  };
  
  for(std::vector< shared_ptr< joint_dependent_gen_coord > >::const_iterator it = mDependentGenCoords.begin(); 
      it < mDependentGenCoords.end(); ++it, ++j)
    result[j] = (*it)->mFrame->q_dot;

  for(std::vector< shared_ptr< joint_dependent_frame_2D > >::const_iterator it = mDependent2DFrames.begin(); 
      it < mDependent2DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Velocity[0]; ++j;
    result[j] = (*it)->mFrame->Velocity[1]; ++j;
    result[j] = (*it)->mFrame->AngVelocity; ++j;
  };

  for(std::vector< shared_ptr< joint_dependent_frame_3D > >::const_iterator it = mDependent3DFrames.begin(); 
      it < mDependent3DFrames.end(); ++it) {
    result[j] = (*it)->mFrame->Velocity[0]; ++j;
    result[j] = (*it)->mFrame->Velocity[1]; ++j;
    result[j] = (*it)->mFrame->Velocity[2]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[0]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[1]; ++j;
    result[j] = (*it)->mFrame->AngVelocity[2]; ++j;
  };
  
  return result;
};






void manipulator_dynamics_model::getMassMatrix(mat<double,mat_structure::symmetric>& M) {
  mMassCalc.getMassMatrix(M);
};

void manipulator_dynamics_model::getMassMatrixAndDerivative(mat<double,mat_structure::symmetric>& M, mat<double,mat_structure::square>& M_dot) {
  mMassCalc.getMassMatrixAndDerivative(M,M_dot);
};

void manipulator_dynamics_model::get_TMT_TdMT(mat<double,mat_structure::rectangular>& Tcm, mat<double,mat_structure::symmetric>& Mcm, mat<double,mat_structure::rectangular>& Tcm_dot) {
  mMassCalc.get_TMT_TdMT(Tcm,Mcm,Tcm_dot);
};




};

};








