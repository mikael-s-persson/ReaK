
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

#include "mass_matrix_calculator.hpp"

namespace ReaK {

namespace kte {


mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< inertia_gen >::type& aGenInertia) {

  if(aGenInertia)
    mGenInertias.push_back(aGenInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< inertia_2D >::type& a2DInertia) {

  if(a2DInertia)
    m2DInertias.push_back(a2DInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< inertia_3D >::type& a3DInertia) {

  if(a3DInertia)
    m3DInertias.push_back(a3DInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< gen_coord<double> >::type& aCoord) {

  if(aCoord)
    mCoords.push_back(aCoord);

  return *this;
};


mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< frame_2D<double> >::type& aFrame2D) {
  
  if(aFrame2D)
    mFrames2D.push_back(aFrame2D);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_pointer< frame_3D<double> >::type& aFrame3D) {

  if(aFrame3D)
    mFrames3D.push_back(aFrame3D);
  
  return *this;
};

void mass_matrix_calc::getMassMatrix(mat<double,mat_structure::symmetric>& M) {
  mat<double,mat_structure::symmetric> Mcm(0);
  mat<double,mat_structure::rectangular> Tcm(0,0);
  mat<double,mat_structure::rectangular> Tcm_dot(0,0);
  get_TMT_TdMT(Tcm,Mcm,Tcm_dot);

  M = transpose(Tcm) * (Mcm * Tcm);
};

void mass_matrix_calc::getMassMatrixAndDerivative(mat<double,mat_structure::symmetric>& M, mat<double,mat_structure::square>& M_dot) {
  mat<double,mat_structure::symmetric> Mcm(0);
  mat<double,mat_structure::rectangular> Tcm(0,0);
  mat<double,mat_structure::rectangular> Tcm_dot(0,0);
  get_TMT_TdMT(Tcm,Mcm,Tcm_dot);

  M = transpose(Tcm) * (Mcm * Tcm);
  M_dot = transpose(Tcm_dot) * (Mcm * Tcm);
  M_dot += transpose(M_dot);
};

void mass_matrix_calc::get_TMT_TdMT(mat<double,mat_structure::rectangular>& Tcm, mat<double,mat_structure::symmetric>& Mcm, mat<double,mat_structure::rectangular>& Tcm_dot) {
  unsigned int m = 6*m3DInertias.size() + 3*m2DInertias.size() + mGenInertias.size();
  unsigned int n = mCoords.size() + 3 * mFrames2D.size() + 6 * mFrames3D.size();
  Mcm.set_col_count(m);
  Tcm_dot = Tcm = mat<double,mat_structure::rectangular>(m,n);

  unsigned int RowInd = 0;
  
/****************************************************************************************
 *                             Gen Coords 
 * *************************************************************************************/
  
  
  for(unsigned int i=0; i<mCoords.size(); ++i) {
    RowInd = 0;
    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->mUpStreamJoints.find(mCoords[i]) != m3DInertias[j]->mUpStreamJoints.end()) {

        shared_pointer< jacobian_gen_3D<double> >::type Jac = m3DInertias[j]->mUpStreamJoints[mCoords[i]];
        frame_3D<double> f2 = m3DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock());

        rot_mat_3D<double> R(f2.Quat.getRotMat());
        vect<double,3> w_tmp = Jac->qd_avel * R;
        Tcm(RowInd,i) = w_tmp[0];
        Tcm(RowInd+1,i) = w_tmp[1];
        Tcm(RowInd+2,i) = w_tmp[2];
        vect<double,3> v_tmp = (Jac->qd_avel % f2.Position
                                    + Jac->qd_vel) * R;
        Tcm(RowInd+3,i) = v_tmp[0];
        Tcm(RowInd+4,i) = v_tmp[1];
        Tcm(RowInd+5,i) = v_tmp[2];
    
        w_tmp = Jac->qd_aacc * R
               - f2.AngVelocity % w_tmp;
        Tcm_dot(RowInd,i) = w_tmp[0];
        Tcm_dot(RowInd+1,i) = w_tmp[1];
        Tcm_dot(RowInd+2,i) = w_tmp[2];
        v_tmp = (Jac->qd_avel % f2.Velocity + Jac->qd_aacc % f2.Position + Jac->qd_acc) * R
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+3,i) = v_tmp[0];
        Tcm_dot(RowInd+4,i) = v_tmp[1];
        Tcm_dot(RowInd+5,i) = v_tmp[2];

      };
      RowInd += 6;
    };
    
    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->mUpStreamJoints.find(mCoords[i]) != m2DInertias[j]->mUpStreamJoints.end()) {

        shared_pointer< jacobian_gen_2D<double> >::type Jac = m2DInertias[j]->mUpStreamJoints[mCoords[i]];
        frame_2D<double> f2 = m2DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock());
        Tcm(RowInd,i) = Jac->qd_avel;

        vect<double,2> v_tmp = (Jac->qd_avel % f2.Position + Jac->qd_vel) * f2.Rotation;
        Tcm(RowInd+1,i) = v_tmp[0];
        Tcm(RowInd+2,i) = v_tmp[1];
        Tcm_dot(RowInd, i) = Jac->qd_aacc;
        v_tmp = (Jac->qd_avel % f2.Velocity + Jac->qd_aacc % f2.Position + Jac->qd_acc) * f2.Rotation
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+1,i) = v_tmp[0];
        Tcm_dot(RowInd+2,i) = v_tmp[1];

      };
      RowInd += 3;
    };
    
    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->mUpStreamJoints.find(mCoords[i]) != mGenInertias[j]->mUpStreamJoints.end()) {
        Tcm(RowInd,i) = mGenInertias[j]->mUpStreamJoints[mCoords[i]]->qd_qd;
        Tcm_dot(RowInd,i) = mGenInertias[j]->mUpStreamJoints[mCoords[i]]->qd_qdd;
      };
      RowInd++;
    };
  };

  
/****************************************************************************************
 *                             2D Frames 
 * *************************************************************************************/
  
  unsigned int base_i = mCoords.size();
  for(unsigned int i=0; i<mFrames2D.size(); ++i) {
    RowInd = 0;

    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->mUpStream2DJoints.find(mFrames2D[i]) != m3DInertias[j]->mUpStream2DJoints.end()) {

        shared_pointer< jacobian_2D_3D<double> >::type Jac = m3DInertias[j]->mUpStream2DJoints[mFrames2D[i]];
        frame_3D<double> f2 = m3DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock());

        //std::cout << "Joint " << i << " at angle " << mCoords[i]->q << " CM " << m3DInertias[j]->getName() << " gives:" << std::endl;
        //std::cout << "\tPosition = " << f2.Position << "\tVelocity = " << f2.Velocity << "\tAccel. = " << f2.Acceleration << std::endl;
        //std::cout << "\tOrientation = " << f2.Quat << "\tAngVelocity = " << f2.AngVelocity << "\tAngAccel. = " << f2.AngAcceleration << std::endl;

        rot_mat_3D<double> R(f2.Quat.getRotMat());
        vect<double,3> w_tmp = Jac->vel_avel[0] * R;
        Tcm(RowInd,  3*i+base_i) = w_tmp[0];
        Tcm(RowInd+1,3*i+base_i) = w_tmp[1];
        Tcm(RowInd+2,3*i+base_i) = w_tmp[2];
        w_tmp = Jac->vel_aacc[0] * R
               - f2.AngVelocity % w_tmp;
        Tcm_dot(RowInd,  3*i+base_i) = w_tmp[0];
        Tcm_dot(RowInd+1,3*i+base_i) = w_tmp[1];
        Tcm_dot(RowInd+2,3*i+base_i) = w_tmp[2];
	
        w_tmp = Jac->vel_avel[1] * R;
        Tcm(RowInd,  3*i+base_i+1) = w_tmp[0];
        Tcm(RowInd+1,3*i+base_i+1) = w_tmp[1];
        Tcm(RowInd+2,3*i+base_i+1) = w_tmp[2];
        w_tmp = Jac->vel_aacc[1] * R
               - f2.AngVelocity % w_tmp;
        Tcm_dot(RowInd,  3*i+base_i+1) = w_tmp[0];
        Tcm_dot(RowInd+1,3*i+base_i+1) = w_tmp[1];
        Tcm_dot(RowInd+2,3*i+base_i+1) = w_tmp[2];
	
	w_tmp = Jac->avel_avel * R;
        Tcm(RowInd,  3*i+base_i+2) = w_tmp[0];
        Tcm(RowInd+1,3*i+base_i+2) = w_tmp[1];
        Tcm(RowInd+2,3*i+base_i+2) = w_tmp[2];
        w_tmp = Jac->avel_aacc * R
               - f2.AngVelocity % w_tmp;
        Tcm_dot(RowInd,  3*i+base_i+2) = w_tmp[0];
        Tcm_dot(RowInd+1,3*i+base_i+2) = w_tmp[1];
        Tcm_dot(RowInd+2,3*i+base_i+2) = w_tmp[2];
	
        vect<double,3> v_tmp = (Jac->vel_avel[0] % f2.Position
                                    + Jac->vel_vel[0]) * R;
        Tcm(RowInd+3,3*i+base_i) = v_tmp[0];
        Tcm(RowInd+4,3*i+base_i) = v_tmp[1];
        Tcm(RowInd+5,3*i+base_i) = v_tmp[2];
        v_tmp = (Jac->vel_avel[0] % f2.Velocity + Jac->vel_aacc[0] % f2.Position + Jac->vel_acc[0]) * R
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+3,3*i+base_i) = v_tmp[0];
        Tcm_dot(RowInd+4,3*i+base_i) = v_tmp[1];
        Tcm_dot(RowInd+5,3*i+base_i) = v_tmp[2];

        v_tmp = (Jac->vel_avel[1] % f2.Position
                                    + Jac->vel_vel[1]) * R;
        Tcm(RowInd+3,3*i+base_i+1) = v_tmp[0];
        Tcm(RowInd+4,3*i+base_i+1) = v_tmp[1];
        Tcm(RowInd+5,3*i+base_i+1) = v_tmp[2];
        v_tmp = (Jac->vel_avel[1] % f2.Velocity + Jac->vel_aacc[1] % f2.Position + Jac->vel_acc[1]) * R
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+3,3*i+base_i+1) = v_tmp[0];
        Tcm_dot(RowInd+4,3*i+base_i+1) = v_tmp[1];
        Tcm_dot(RowInd+5,3*i+base_i+1) = v_tmp[2];
	
        v_tmp = (Jac->avel_avel % f2.Position
                                    + Jac->avel_vel) * R;
        Tcm(RowInd+3,3*i+base_i+2) = v_tmp[0];
        Tcm(RowInd+4,3*i+base_i+2) = v_tmp[1];
        Tcm(RowInd+5,3*i+base_i+2) = v_tmp[2];	
        v_tmp = (Jac->avel_avel % f2.Velocity + Jac->avel_aacc % f2.Position + Jac->avel_acc) * R
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+3,3*i+base_i+2) = v_tmp[0];
        Tcm_dot(RowInd+4,3*i+base_i+2) = v_tmp[1];
        Tcm_dot(RowInd+5,3*i+base_i+2) = v_tmp[2];

      };
      RowInd += 6;
    };

    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->mUpStream2DJoints.find(mFrames2D[i]) != m2DInertias[j]->mUpStream2DJoints.end()) {

        shared_pointer< jacobian_2D_2D<double> >::type Jac = m2DInertias[j]->mUpStream2DJoints[mFrames2D[i]];
        frame_2D<double> f2 = m2DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock());
        Tcm(RowInd,3*i+base_i) = Jac->vel_avel[0];
	Tcm(RowInd,3*i+base_i+1) = Jac->vel_avel[1];
	Tcm(RowInd,3*i+base_i+2) = Jac->avel_avel;
        Tcm_dot(RowInd, 3*i+base_i) = Jac->vel_aacc[0];
        Tcm_dot(RowInd, 3*i+base_i+1) = Jac->vel_aacc[1];
        Tcm_dot(RowInd, 3*i+base_i+2) = Jac->avel_aacc;

        vect<double,2> v_tmp = (Jac->vel_avel[0] % f2.Position + Jac->vel_vel[0]) * f2.Rotation;
        Tcm(RowInd+1,3*i+base_i) = v_tmp[0];
        Tcm(RowInd+2,3*i+base_i) = v_tmp[1];
        v_tmp = (Jac->vel_avel[0] % f2.Velocity + Jac->vel_aacc[0] % f2.Position + Jac->vel_acc[0]) * f2.Rotation
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+1,3*i+base_i) = v_tmp[0];
        Tcm_dot(RowInd+2,3*i+base_i) = v_tmp[1];
	
	v_tmp = (Jac->vel_avel[1] % f2.Position + Jac->vel_vel[1]) * f2.Rotation;
        Tcm(RowInd+1,3*i+base_i+1) = v_tmp[0];
        Tcm(RowInd+2,3*i+base_i+1) = v_tmp[1];
        v_tmp = (Jac->vel_avel[1] % f2.Velocity + Jac->vel_aacc[1] % f2.Position + Jac->vel_acc[1]) * f2.Rotation
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+1,3*i+base_i+1) = v_tmp[0];
        Tcm_dot(RowInd+2,3*i+base_i+1) = v_tmp[1];

	v_tmp = (Jac->avel_avel % f2.Position + Jac->avel_vel) * f2.Rotation;
        Tcm(RowInd+1,3*i+base_i+2) = v_tmp[0];
        Tcm(RowInd+2,3*i+base_i+2) = v_tmp[1];
        v_tmp = (Jac->avel_avel % f2.Velocity + Jac->avel_aacc % f2.Position + Jac->avel_acc) * f2.Rotation
               - f2.AngVelocity % v_tmp;
        Tcm_dot(RowInd+1,3*i+base_i+2) = v_tmp[0];
        Tcm_dot(RowInd+2,3*i+base_i+2) = v_tmp[1];
	
      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->mUpStream2DJoints.find(mFrames2D[i]) != mGenInertias[j]->mUpStream2DJoints.end()) {
        shared_pointer< jacobian_2D_gen<double> >::type Jac = mGenInertias[j]->mUpStream2DJoints[mFrames2D[i]];
	Tcm(RowInd,3*i+base_i) = Jac->vel_qd[0];
        Tcm(RowInd,3*i+base_i+1) = Jac->vel_qd[1];
        Tcm(RowInd,3*i+base_i+2) = Jac->avel_qd;
        Tcm_dot(RowInd,3*i+base_i) = Jac->vel_qdd[0];
        Tcm_dot(RowInd,3*i+base_i+1) = Jac->vel_qdd[1];
        Tcm_dot(RowInd,3*i+base_i+2) = Jac->avel_qdd;
      };
      RowInd++;
    };
  };

  
  
/****************************************************************************************
 *                             3D Frames 
 * *************************************************************************************/
  
  base_i = mCoords.size() + 3*mFrames2D.size();
  for(unsigned int i=0; i<mFrames3D.size(); ++i) {
    RowInd = 0; 

    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->mUpStream3DJoints.find(mFrames3D[i]) != m3DInertias[j]->mUpStream3DJoints.end()) {
        shared_pointer< jacobian_3D_3D<double> >::type Jac = m3DInertias[j]->mUpStream3DJoints[mFrames3D[i]]; 
        frame_3D<double> f2 = m3DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock()); 
        
        //std::cout << "Joint " << i << " at angle " << mCoords[i]->q << " CM " << m3DInertias[j]->getName() << " gives:" << std::endl;
        //std::cout << "\tPosition = " << f2.Position << "\tVelocity = " << f2.Velocity << "\tAccel. = " << f2.Acceleration << std::endl;
        //std::cout << "\tOrientation = " << f2.Quat << "\tAngVelocity = " << f2.AngVelocity << "\tAngAccel. = " << f2.AngAcceleration << std::endl;

        rot_mat_3D<double> R(f2.Quat.getRotMat()); 
        for(unsigned int k = 0; k < 3; ++k) {
	  vect<double,3> w_tmp = Jac->vel_avel[k] * R;
          Tcm(RowInd,6*i+base_i+k) = w_tmp[0];
          Tcm(RowInd+1,6*i+base_i+k) = w_tmp[1]; 
          Tcm(RowInd+2,6*i+base_i+k) = w_tmp[2];
          w_tmp = Jac->vel_aacc[k] * R
                 - f2.AngVelocity % w_tmp;
          Tcm_dot(RowInd,6*i+base_i+k) = w_tmp[0];
          Tcm_dot(RowInd+1,6*i+base_i+k) = w_tmp[1]; 
          Tcm_dot(RowInd+2,6*i+base_i+k) = w_tmp[2];
	
	  vect<double,3> v_tmp = (Jac->vel_avel[k] % f2.Position
                                    + Jac->vel_vel[k]) * R;
          Tcm(RowInd+3,6*i+base_i+k) = v_tmp[0];
          Tcm(RowInd+4,6*i+base_i+k) = v_tmp[1]; 
          Tcm(RowInd+5,6*i+base_i+k) = v_tmp[2];

          v_tmp = (Jac->vel_avel[k] % f2.Velocity + Jac->vel_aacc[k] % f2.Position + Jac->vel_acc[k]) * R
                 - f2.AngVelocity % v_tmp;
          Tcm_dot(RowInd+3,6*i+base_i+k) = v_tmp[0];
          Tcm_dot(RowInd+4,6*i+base_i+k) = v_tmp[1]; 
          Tcm_dot(RowInd+5,6*i+base_i+k) = v_tmp[2];
	};
        for(unsigned int k = 0; k < 3; ++k) {
	  vect<double,3> w_tmp = Jac->avel_avel[k] * R;
          Tcm(RowInd,6*i+base_i+3+k) = w_tmp[0];
          Tcm(RowInd+1,6*i+base_i+3+k) = w_tmp[1]; 
          Tcm(RowInd+2,6*i+base_i+3+k) = w_tmp[2];
          w_tmp = Jac->avel_aacc[k] * R
                 - f2.AngVelocity % w_tmp;
          Tcm_dot(RowInd,6*i+base_i+3+k) = w_tmp[0];
          Tcm_dot(RowInd+1,6*i+base_i+3+k) = w_tmp[1]; 
          Tcm_dot(RowInd+2,6*i+base_i+3+k) = w_tmp[2];
	
	  vect<double,3> v_tmp = (Jac->avel_avel[k] % f2.Position
                                    + Jac->avel_vel[k]) * R;
          Tcm(RowInd+3,6*i+base_i+3+k) = v_tmp[0];
          Tcm(RowInd+4,6*i+base_i+3+k) = v_tmp[1]; 
          Tcm(RowInd+5,6*i+base_i+3+k) = v_tmp[2];

          v_tmp = (Jac->avel_avel[k] % f2.Velocity + Jac->avel_aacc[k] % f2.Position + Jac->avel_acc[k]) * R
                 - f2.AngVelocity % v_tmp;
          Tcm_dot(RowInd+3,6*i+base_i+3+k) = v_tmp[0];
          Tcm_dot(RowInd+4,6*i+base_i+3+k) = v_tmp[1]; 
          Tcm_dot(RowInd+5,6*i+base_i+3+k) = v_tmp[2];
	};
        
      };
      RowInd += 6;
    };

    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->mUpStream3DJoints.find(mFrames3D[i]) != m2DInertias[j]->mUpStream3DJoints.end()) {

        shared_pointer< jacobian_3D_2D<double> >::type Jac = m2DInertias[j]->mUpStream3DJoints[mFrames3D[i]];
        frame_2D<double> f2 = m2DInertias[j]->mCenterOfMass->getFrameRelativeTo(Jac->Parent.lock());
       
	for(unsigned int k=0; k < 3; ++k) {
	  Tcm(RowInd,6*i+base_i+k) = Jac->vel_avel[k];
          Tcm_dot(RowInd, 6*i+base_i+k) = Jac->vel_aacc[k];
        
          vect<double,2> v_tmp = (Jac->vel_avel[k] % f2.Position + Jac->vel_vel[k]) * f2.Rotation;
          Tcm(RowInd+1,6*i+base_i+k) = v_tmp[0];
          Tcm(RowInd+2,6*i+base_i+k) = v_tmp[1];
          v_tmp = (Jac->vel_avel[k] % f2.Velocity + Jac->vel_aacc[k] % f2.Position + Jac->vel_acc[k]) * f2.Rotation
                 - f2.AngVelocity % v_tmp;
          Tcm_dot(RowInd+1,6*i+base_i+k) = v_tmp[0];
          Tcm_dot(RowInd+2,6*i+base_i+k) = v_tmp[1];
	};
	for(unsigned int k=0; k < 3; ++k) {
	  Tcm(RowInd,6*i+base_i+3+k) = Jac->avel_avel[k];
          Tcm_dot(RowInd, 6*i+base_i+3+k) = Jac->avel_aacc[k];
        
          vect<double,2> v_tmp = (Jac->avel_avel[k] % f2.Position + Jac->avel_vel[k]) * f2.Rotation;
          Tcm(RowInd+1,6*i+base_i+3+k) = v_tmp[0];
          Tcm(RowInd+2,6*i+base_i+3+k) = v_tmp[1];
          v_tmp = (Jac->avel_avel[k] % f2.Velocity + Jac->avel_aacc[k] % f2.Position + Jac->avel_acc[k]) * f2.Rotation
                 - f2.AngVelocity % v_tmp;
          Tcm_dot(RowInd+1,6*i+base_i+3+k) = v_tmp[0];
          Tcm_dot(RowInd+2,6*i+base_i+3+k) = v_tmp[1];
	};
	
      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->mUpStreamJoints.find(mCoords[i]) != mGenInertias[j]->mUpStreamJoints.end()) {
        shared_pointer< jacobian_3D_gen<double> >::type Jac = mGenInertias[j]->mUpStream3DJoints[mFrames3D[i]];
	
	for(unsigned int k=0; k < 3; ++k) {
	  Tcm(RowInd,6*i+base_i+k) = Jac->vel_qd[k];
          Tcm_dot(RowInd,6*i+base_i+k) = Jac->vel_qdd[k];
	};
	for(unsigned int k=0; k < 3; ++k) {
	  Tcm(RowInd,6*i+base_i+3+k) = Jac->avel_qd[k];
          Tcm_dot(RowInd,6*i+base_i+3+k) = Jac->avel_qdd[k];
	};
	
      };
      RowInd++;
    };
  };

  
  RowInd = 0;
  for(unsigned int j=0; j<m3DInertias.size();++j) {
    set_block(Mcm,m3DInertias[j]->mInertiaTensor,RowInd); RowInd += 3;
    Mcm(RowInd,RowInd) = m3DInertias[j]->mMass; RowInd++;
    Mcm(RowInd,RowInd) = m3DInertias[j]->mMass; RowInd++;
    Mcm(RowInd,RowInd) = m3DInertias[j]->mMass; RowInd++;
  };
  for(unsigned int j=0; j<m2DInertias.size(); ++j) {
    Mcm(RowInd,RowInd) = m2DInertias[j]->mMomentOfInertia; RowInd++;
    Mcm(RowInd,RowInd) = m2DInertias[j]->mMass; RowInd++;
    Mcm(RowInd,RowInd) = m2DInertias[j]->mMass; RowInd++;
  };
  for(unsigned int j=0; j<mGenInertias.size(); ++j) {
    Mcm(RowInd,RowInd) = mGenInertias[j]->mMass; RowInd++;
  };

};




};

};








