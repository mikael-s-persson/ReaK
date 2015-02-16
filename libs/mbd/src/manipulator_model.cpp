
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

#include <ReaK/mbd/kte/manipulator_model.hpp>

#include <ReaK/math/lin_alg/mat_num_exceptions.hpp>
#include <ReaK/math/lin_alg/mat_cholesky.hpp>

#include <ReaK/mbd/kte/manipulator_model_helper.hpp>

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


void manipulator_kinematics_model::getJacobianMatrix(mat<double,mat_structure::rectangular>& Jac) const {
  manip_kin_mdl_jac_calculator(this).getJacobianMatrix(Jac);
};

void manipulator_kinematics_model::getJacobianMatrixAndDerivative(mat<double,mat_structure::rectangular>& Jac, 
                                                                  mat<double,mat_structure::rectangular>& JacDot) const {
  manip_kin_mdl_jac_calculator(this).getJacobianMatrixAndDerivative(Jac,JacDot);
};
  






vect_n<double> manipulator_kinematics_model::getJointPositions() const {
  vect_n<double> result(getJointPositionsCount());
  
  manip_kin_mdl_joint_io(this).getJointPositions(result);
  
  return result;
};
    
void manipulator_kinematics_model::setJointPositions(const vect_n<double>& aJointPositions) {
  if(aJointPositions.size() != getJointPositionsCount())
    throw std::range_error("Joint-position vector has incorrect dimensions!");
  
  manip_kin_mdl_joint_io(this).setJointPositions(aJointPositions);
};
    
vect_n<double> manipulator_kinematics_model::getJointVelocities() const {
  vect_n<double> result(getJointVelocitiesCount());
  
  manip_kin_mdl_joint_io(this).getJointVelocities(result);
  
  return result;  
};
    
void manipulator_kinematics_model::setJointVelocities(const vect_n<double>& aJointVelocities) {
  if(aJointVelocities.size() != getJointVelocitiesCount())
    throw std::range_error("Joint-velocity vector has incorrect dimensions!");
  
  manip_kin_mdl_joint_io(this).setJointVelocities(aJointVelocities);
};
    
vect_n<double> manipulator_kinematics_model::getJointAccelerations() const {
  vect_n<double> result(getJointAccelerationsCount());
  
  manip_kin_mdl_joint_io(this).getJointAccelerations(result);
  
  return result;  
};
    
void manipulator_kinematics_model::setJointAccelerations(const vect_n<double>& aJointAccelerations) {
  if(aJointAccelerations.size() != getJointAccelerationsCount())
    throw std::range_error("Joint-acceleration vector has incorrect dimensions!");
  
  manip_kin_mdl_joint_io(this).setJointAccelerations(aJointAccelerations);
};
    
vect_n<double> manipulator_kinematics_model::getDependentPositions() const {
  vect_n<double> result(getDependentPositionsCount());
  
  manip_kin_mdl_joint_io(this).getDependentPositions(result);
  
  return result;
};
    
vect_n<double> manipulator_kinematics_model::getDependentVelocities() const {
  vect_n<double> result(getDependentVelocitiesCount());
  
  manip_kin_mdl_joint_io(this).getDependentVelocities(result);
  
  return result;  
};
    
vect_n<double> manipulator_kinematics_model::getDependentAccelerations() const {
  vect_n<double> result(getDependentAccelerationsCount());
  
  manip_kin_mdl_joint_io(this).getDependentAccelerations(result);
  
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
  
  vect_ref_view< vect_n<double> > pos_result = result[range(0,getJointPositionsCount())];
  manip_kin_mdl_joint_io(this).getJointPositions(pos_result);
  
  vect_ref_view< vect_n<double> > vel_result = result[range(getJointPositionsCount(),getJointPositionsCount() + getJointVelocitiesCount())];
  manip_kin_mdl_joint_io(this).getJointVelocities(vel_result);
  
  return result;
};
    
void manipulator_dynamics_model::setJointStates(const vect_n<double>& aJointStates) {
  if(aJointStates.size() != getJointStatesCount())
    throw std::range_error("Joint-state vector has incorrect dimensions!");
  
  vect_const_ref_view< vect_n<double> > positions = aJointStates[range(0,getJointPositionsCount())];
  manip_kin_mdl_joint_io(this).setJointPositions(positions);
  
  vect_const_ref_view< vect_n<double> > velocities = aJointStates[range(getJointPositionsCount(),getJointPositionsCount() + getJointVelocitiesCount())];
  manip_kin_mdl_joint_io(this).setJointVelocities(velocities);
  
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
      (*it)->setInput(k,aInput[i++]);
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
  } catch(singularity_error& e) { RK_UNUSED(e);
    std::stringstream ss; ss << "Mass matrix is singular in the manipulator model '" << getName() << "' at time " << aTime << " seconds.";
    throw singularity_error(ss.str());
  };
};
    
vect_n<double> manipulator_dynamics_model::getDependentStates() const {
  vect_n<double> result(getDependentStatesCount());
  
  vect_ref_view< vect_n<double> > pos_result = result[range(0,getDependentPositionsCount())];
  manip_kin_mdl_joint_io(this).getDependentPositions(pos_result);
  
  vect_ref_view< vect_n<double> > vel_result = result[range(getDependentPositionsCount(),getDependentPositionsCount() + getDependentVelocitiesCount())];
  manip_kin_mdl_joint_io(this).getDependentVelocities(vel_result);
  
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








