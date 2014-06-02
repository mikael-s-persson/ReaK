
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

#include <ReaK/ctrl/mbd_kte/mass_matrix_calculator.hpp>

namespace ReaK {

namespace kte {


mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< inertia_gen >& aGenInertia) {

  if(aGenInertia)
    mGenInertias.push_back(aGenInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< inertia_2D >& a2DInertia) {

  if(a2DInertia)
    m2DInertias.push_back(a2DInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< inertia_3D >& a3DInertia) {

  if(a3DInertia)
    m3DInertias.push_back(a3DInertia);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< gen_coord<double> >& aCoord) {

  if(aCoord)
    mCoords.push_back(aCoord);

  return *this;
};


mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< frame_2D<double> >& aFrame2D) {
  
  if(aFrame2D)
    mFrames2D.push_back(aFrame2D);

  return *this;
};

mass_matrix_calc& mass_matrix_calc::operator <<(const shared_ptr< frame_3D<double> >& aFrame3D) {

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
  Mcm = mat<double,mat_structure::symmetric>(m,0.0);
  Tcm_dot = Tcm = mat<double,mat_structure::nil>(m,n);

  unsigned int RowInd = 0;
  
/****************************************************************************************
 *                             Gen Coords 
 * *************************************************************************************/
  
  
  for(unsigned int i=0; i<mCoords.size(); ++i) {
    RowInd = 0;
    
    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->CenterOfMass()->mUpStreamJoints.find(mCoords[i]) != mGenInertias[j]->CenterOfMass()->mUpStreamJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,1,1,RowInd,i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,1,1,RowInd,i);
        mGenInertias[j]
          ->CenterOfMass()
            ->mUpStreamJoints[mCoords[i]]
              ->write_to_matrices(subTcm,subTcm_dot);
        
      };
      RowInd++;
    };
    
    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->CenterOfMass()->mUpStreamJoints.find(mCoords[i]) != m2DInertias[j]->CenterOfMass()->mUpStreamJoints.end()) {

        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,3,1,RowInd,i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,3,1,RowInd,i);
        m2DInertias[j]
          ->CenterOfMass()
            ->mUpStreamJoints[mCoords[i]]
              ->get_jac_relative_to(m2DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);
        
      };
      RowInd += 3;
    };
    
    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->CenterOfMass()->mUpStreamJoints.find(mCoords[i]) != m3DInertias[j]->CenterOfMass()->mUpStreamJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,6,1,RowInd,i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,6,1,RowInd,i);
        m3DInertias[j]
          ->CenterOfMass()
            ->mUpStreamJoints[mCoords[i]]
              ->get_jac_relative_to(m3DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);

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

    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->CenterOfMass()->mUpStream2DJoints.find(mFrames2D[i]) != mGenInertias[j]->CenterOfMass()->mUpStream2DJoints.end()) {
                
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,1,3,RowInd,3*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,1,3,RowInd,3*i+base_i);
        mGenInertias[j]
          ->CenterOfMass()
            ->mUpStream2DJoints[mFrames2D[i]]
              ->write_to_matrices(subTcm,subTcm_dot);
        
      };
      RowInd++;
    };

    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->CenterOfMass()->mUpStream2DJoints.find(mFrames2D[i]) != m2DInertias[j]->CenterOfMass()->mUpStream2DJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,3,3,RowInd,3*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,3,3,RowInd,3*i+base_i);
        m2DInertias[j]
          ->CenterOfMass()
            ->mUpStream2DJoints[mFrames2D[i]]
              ->get_jac_relative_to(m2DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);
              
      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->CenterOfMass()->mUpStream2DJoints.find(mFrames2D[i]) != m3DInertias[j]->CenterOfMass()->mUpStream2DJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,6,3,RowInd,3*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,6,3,RowInd,3*i+base_i);
        m3DInertias[j]
          ->CenterOfMass()
            ->mUpStream2DJoints[mFrames2D[i]]
              ->get_jac_relative_to(m3DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);

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

    for(unsigned int j=0; j<mGenInertias.size(); ++j) {
      if(mGenInertias[j]->CenterOfMass()->mUpStreamJoints.find(mCoords[i]) != mGenInertias[j]->CenterOfMass()->mUpStreamJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,1,6,RowInd,6*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,1,6,RowInd,6*i+base_i);
        mGenInertias[j]
          ->CenterOfMass()
            ->mUpStream3DJoints[mFrames3D[i]]
              ->write_to_matrices(subTcm,subTcm_dot);
        
      };
      RowInd++;
    };

    for(unsigned int j=0; j<m2DInertias.size(); ++j) {
      if(m2DInertias[j]->CenterOfMass()->mUpStream3DJoints.find(mFrames3D[i]) != m2DInertias[j]->CenterOfMass()->mUpStream3DJoints.end()) {
        
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,3,6,RowInd,6*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,3,6,RowInd,6*i+base_i);
        m2DInertias[j]
          ->CenterOfMass()
            ->mUpStream3DJoints[mFrames3D[i]]
              ->get_jac_relative_to(m2DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);

      };
      RowInd += 3;
    };

    for(unsigned int j=0; j<m3DInertias.size(); ++j) {
      if(m3DInertias[j]->CenterOfMass()->mUpStream3DJoints.find(mFrames3D[i]) != m3DInertias[j]->CenterOfMass()->mUpStream3DJoints.end()) {
                
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm(Tcm,6,6,RowInd,6*i+base_i);
        mat_sub_block< mat<double,mat_structure::rectangular> > subTcm_dot(Tcm_dot,6,6,RowInd,6*i+base_i);
        m3DInertias[j]
          ->CenterOfMass()
            ->mUpStream3DJoints[mFrames3D[i]]
              ->get_jac_relative_to(m3DInertias[j]->CenterOfMass()->mFrame)
                .write_to_matrices(subTcm,subTcm_dot);
        
      };
      RowInd += 6;
    };
  };

  
  RowInd = 0;
  for(unsigned int j=0; j<mGenInertias.size(); ++j) {
    Mcm(RowInd,RowInd) = mGenInertias[j]->Mass(); RowInd++;
  };
  for(unsigned int j=0; j<m2DInertias.size(); ++j) {
    Mcm(RowInd,RowInd) = m2DInertias[j]->Mass(); RowInd++;
    Mcm(RowInd,RowInd) = m2DInertias[j]->Mass(); RowInd++;
    Mcm(RowInd,RowInd) = m2DInertias[j]->MomentOfInertia(); RowInd++;
  };
  for(unsigned int j=0; j<m3DInertias.size();++j) {
    Mcm(RowInd,RowInd) = m3DInertias[j]->Mass(); RowInd++;
    Mcm(RowInd,RowInd) = m3DInertias[j]->Mass(); RowInd++;
    Mcm(RowInd,RowInd) = m3DInertias[j]->Mass(); RowInd++;
    set_block(Mcm,m3DInertias[j]->InertiaTensor(),RowInd); RowInd += 3;
  };

};




};

};








