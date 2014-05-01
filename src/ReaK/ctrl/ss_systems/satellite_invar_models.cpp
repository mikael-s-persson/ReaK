
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

#include "satellite_invar_models.hpp"

#include "ctrl_sys/sss_exceptions.hpp"
#include "topologies/se3_topologies.hpp"
#include "topologies/se2_topologies.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


  

satellite2D_imdt_sys::satellite2D_imdt_sys(
  const std::string& aName, double aMass, double aInertiaMoment, double aDt) :
  named_object(), mMass(aMass), mInertiaMoment(aInertiaMoment), mDt(aDt) {
  setName(aName);
  if((mInertiaMoment < std::numeric_limits< double >::epsilon()) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information are singular in airship2D_lin_system's definition");
};  



satellite2D_imdt_sys::point_type satellite2D_imdt_sys::get_next_state(
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::point_type& x, 
  const satellite2D_imdt_sys::input_type& u, 
  const satellite2D_imdt_sys::time_type&) const {
  vect<double,2> dv(mDt / mMass * u[0], mDt / mMass * u[1]);
  double dw = mDt / mInertiaMoment * u[2];
  double delta_theta = mDt * (get_ang_velocity(x) + 0.5 * dw); 
  rot_mat_2D<double> r_new = rot_mat_2D<double>(get_rotation(x)) * rot_mat_2D<double>(delta_theta);
  return make_arithmetic_tuple(
    make_arithmetic_tuple(
      get_position(x) + mDt * (get_velocity(x) + 0.5 * dv),
      get_velocity(x) + dv
    ),
    make_arithmetic_tuple(
      r_new.getAngle(),
      get_ang_velocity(x) + dw
    )
  );
};



satellite2D_imdt_sys::output_type satellite2D_imdt_sys::get_output(
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::point_type& x, 
  const satellite2D_imdt_sys::input_type&, 
  const satellite2D_imdt_sys::time_type&) const {
  return satellite2D_imdt_sys::output_type(get_position(x)[0], get_position(x)[1], get_rotation(x));
};



void satellite2D_imdt_sys::get_state_transition_blocks(
  satellite2D_imdt_sys::matrixA_type& A, 
  satellite2D_imdt_sys::matrixB_type& B, 
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::time_type&, 
  const satellite2D_imdt_sys::time_type&, 
  const satellite2D_imdt_sys::point_type&, 
  const satellite2D_imdt_sys::point_type&, 
  const satellite2D_imdt_sys::input_type&, 
  const satellite2D_imdt_sys::input_type&) const {
  A = mat<double,mat_structure::identity>(6);
  A(0,2) = mDt / mMass;
  A(1,3) = mDt / mMass;  
  A(4,5) = mDt / mInertiaMoment;
  
  B = mat<double,mat_structure::nil>(6,3);
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,0) = mDt;
  B(3,1) = mDt;
  B(4,2) = 0.5 * mDt * mDt / mInertiaMoment;
  B(5,2) = mDt;
};



void satellite2D_imdt_sys::get_output_function_blocks(
  satellite2D_imdt_sys::matrixC_type& C, 
  satellite2D_imdt_sys::matrixD_type& D, 
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::time_type&, 
  const satellite2D_imdt_sys::point_type&, 
  const satellite2D_imdt_sys::input_type&) const {
  C = mat<double,mat_structure::nil>(3,6);
  set_block(C,mat_ident<double>(2),0,0);
  set_block(C,mat_ident<double>(1),2,4);
  
  D = mat<double,mat_structure::nil>(3,3);
};



satellite2D_imdt_sys::invariant_error_type satellite2D_imdt_sys::get_invariant_error(
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::point_type& x, 
  const satellite2D_imdt_sys::input_type&, 
  const satellite2D_imdt_sys::output_type& y, 
  const satellite2D_imdt_sys::time_type&) const {
  rot_mat_2D<double> r_y( y[2] );
      
  return satellite2D_imdt_sys::invariant_error_type(
    y[0] - get_position(x)[0],
    y[1] - get_position(x)[1],
    (invert(rot_mat_2D<double>(get_rotation(x))) * r_y).getAngle()); // s_y * c_x - c_y * s_x
};



satellite2D_imdt_sys::point_type satellite2D_imdt_sys::apply_correction(
  const satellite2D_imdt_sys::state_space_type&, 
  const satellite2D_imdt_sys::point_type& x, 
  const satellite2D_imdt_sys::invariant_correction_type& c, 
  const satellite2D_imdt_sys::input_type&, 
  const satellite2D_imdt_sys::time_type&) const {
  return make_arithmetic_tuple(
    make_arithmetic_tuple(
      get_position(x) + vect<double,2>(c[0],c[1]),
      get_velocity(x) + vect<double,2>(c[2] / mMass, c[3] / mMass)
    ),
    make_arithmetic_tuple(
      (rot_mat_2D<double>(get_rotation(x)) * rot_mat_2D<double>(c[4])).getAngle(),
      get_ang_velocity(x) + c[5]
    )
  );
};







satellite3D_imdt_sys::satellite3D_imdt_sys(
  const std::string& aName,
  double aMass,
  const mat<double,mat_structure::symmetric>& aInertiaMoment,
  double aDt, std::size_t aApproxOrder) :
  satellite3D_inv_dt_system(aName, aMass, aInertiaMoment, aDt), approx_order(aApproxOrder) { }; 


void satellite3D_imdt_sys::get_state_transition_blocks(
  satellite3D_imdt_sys::matrixA_type& A, 
  satellite3D_imdt_sys::matrixB_type& B, 
  const satellite3D_imdt_sys::state_space_type&, 
  const satellite3D_imdt_sys::time_type&, 
  const satellite3D_imdt_sys::time_type&,
  const satellite3D_imdt_sys::point_type& p_0, 
  const satellite3D_imdt_sys::point_type& p_1,
  const satellite3D_imdt_sys::input_type&, 
  const satellite3D_imdt_sys::input_type&) const {
  
  
  A = mat_ident<double>(12);
  A(0,3) = mDt;
  A(1,4) = mDt;
  A(2,5) = mDt;
  
  if(approx_order > 1) { // second-order approx:
    mat<double, mat_structure::square> R = (invert( get_quaternion(p_1) ) * get_quaternion(p_0)).as_rotation().getMat();
    mat<double, mat_structure::square> RJRJ = mat<double, mat_structure::square>(transpose_view(R) * mInertiaMomentInv * R * mInertiaMoment);
    set_block(A, (0.5 * mDt) * (mat_ident<double>(3) + RJRJ), 6, 9);
    set_block(A, RJRJ, 9, 9);
  } else { // first-order approx:
    A(6, 9) = mDt;  //
    A(7,10) = mDt;  // set_block(A, mDt * mat_ident<double>(3), 6, 9);
    A(8,11) = mDt;  // 
  };
  
  B = mat<double,mat_structure::nil>(12,6);
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,2) = 0.5 * mDt * mDt / mMass;
  B(3,0) = mDt / mMass;
  B(4,1) = mDt / mMass;
  B(5,2) = mDt / mMass;
  set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 6, 3);
  set_block(B, mDt * mInertiaMomentInv, 9, 3);
  
};



satellite3D_imdt_sys::point_type satellite3D_imdt_sys::apply_correction(
  const satellite3D_imdt_sys::state_space_type&, 
  const satellite3D_imdt_sys::point_type& x, 
  const satellite3D_imdt_sys::invariant_correction_type& c, 
  const satellite3D_imdt_sys::input_type&, 
  const satellite3D_imdt_sys::time_type&) const {
  
//   std::cout << " (Apply Correction) Got quaternion: " << get_quaternion(x) << std::endl;
//   std::cout << " (Apply Correction) Got position: " << get_position(x) << std::endl;
//   std::cout << " (Apply Correction) Got correction: " << c << std::endl;
  
  unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[6],c[7],c[8]) );
//   std::cout << " (Apply Correction) Got quat-diff: " << q_diff << std::endl;
  unit_quat<double> q_new = get_quaternion(x) * q_diff;
//   std::cout << " (Apply Correction) Got quat-new: " << q_new << std::endl;
  
  vect<double,3> w_new = mInertiaMomentInv * (invert(q_diff).as_rotation() * (mInertiaMoment * get_ang_velocity(x) + vect<double,3>(c[9],c[10],c[11])));
//   std::cout << " (Apply Correction) Got w-new: " << w_new << std::endl;
  
//   std::cout << " (Apply Correction) Got corrected state: " << 
//     to_vect<double>(satellite3D_imdt_sys::point_type(
//       make_arithmetic_tuple(
//         get_position(x) + vect<double,3>(c[0],c[1],c[2]),
//         get_velocity(x) + vect<double,3>(c[3],c[4],c[5])
//       ),
//       make_arithmetic_tuple(
//         q_new,
//         w_new
//       )
//     )) << std::endl;
  
  return satellite3D_imdt_sys::point_type(
    make_arithmetic_tuple(
      get_position(x) + vect<double,3>(c[0],c[1],c[2]),
      get_velocity(x) + vect<double,3>(c[3],c[4],c[5])
    ),
    make_arithmetic_tuple(
      q_new,
      w_new
    )
  );
};



satellite3D_imdt_sys::invariant_frame_type satellite3D_imdt_sys::get_invariant_prior_frame(
  const satellite3D_imdt_sys::state_space_type&, 
  const satellite3D_imdt_sys::point_type& x_prev, 
  const satellite3D_imdt_sys::point_type& x_prior, 
  const satellite3D_imdt_sys::input_type&, 
  const satellite3D_imdt_sys::time_type&) const {
  
  satellite3D_imdt_sys::invariant_frame_type result(mat<double,mat_structure::identity>(12));
  mat<double,mat_structure::square> R_diff((invert(get_quaternion(x_prior)) * get_quaternion(x_prev)).as_rotation().getMat());
  set_block(result, R_diff, 6, 6);
  set_block(result, R_diff, 9, 9);
  return result;
};





satellite3D_gyro_imdt_sys::satellite3D_gyro_imdt_sys(
  const std::string& aName, 
  double aMass, 
  const mat<double,mat_structure::symmetric>& aInertiaMoment,
  double aDt, std::size_t aApproxOrder) :
  satellite3D_gyro_inv_dt_system(aName, aMass, aInertiaMoment, aDt), approx_order(aApproxOrder) { }; 



void satellite3D_gyro_imdt_sys::get_state_transition_blocks(
  satellite3D_gyro_imdt_sys::matrixA_type& A, 
  satellite3D_gyro_imdt_sys::matrixB_type& B, 
  const satellite3D_gyro_imdt_sys::state_space_type&, 
  const satellite3D_gyro_imdt_sys::time_type&, 
  const satellite3D_gyro_imdt_sys::time_type&,
  const satellite3D_gyro_imdt_sys::point_type& p_0, 
  const satellite3D_gyro_imdt_sys::point_type& p_1,
  const satellite3D_gyro_imdt_sys::input_type&, 
  const satellite3D_gyro_imdt_sys::input_type&) const {
  
  
  A = mat_ident<double>(12);
  A(0,3) = mDt;
  A(1,4) = mDt;
  A(2,5) = mDt;
  
  if(approx_order > 1) { // second-order approx:
    mat<double, mat_structure::square> R = (invert( get_quaternion(p_1) ) * get_quaternion(p_0)).as_rotation().getMat();
    mat<double, mat_structure::square> RJRJ = mat<double, mat_structure::square>(transpose_view(R) * mInertiaMomentInv * R * mInertiaMoment);
    set_block(A, (0.5 * mDt) * (mat_ident<double>(3) + RJRJ), 6, 9);
    set_block(A, RJRJ, 9, 9);
  } else { // first-order approx:
    A(6, 9) = mDt;  //
    A(7,10) = mDt;  // set_block(A, mDt * mat_ident<double>(3), 6, 9);
    A(8,11) = mDt;  // 
  };
  
  B = mat<double,mat_structure::nil>(12,6);
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,2) = 0.5 * mDt * mDt / mMass;
  B(3,0) = mDt / mMass;
  B(4,1) = mDt / mMass;
  B(5,2) = mDt / mMass;
  set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 6, 3);
  set_block(B, mDt * mInertiaMomentInv, 9, 3);
  
};


  
satellite3D_gyro_imdt_sys::output_type satellite3D_gyro_imdt_sys::get_output(
  const satellite3D_gyro_imdt_sys::state_space_type&, 
  const satellite3D_gyro_imdt_sys::point_type& x, 
  const satellite3D_gyro_imdt_sys::input_type&, 
  const satellite3D_gyro_imdt_sys::time_type&) const {
  
  const vect<double,3>& pos = get_position(x);
  const unit_quat<double>& q = get_quaternion(x);
  const vect<double,3>& w = get_ang_velocity(x);
  return satellite3D_gyro_imdt_sys::output_type(pos[0], pos[1], pos[2], q[0], q[1], q[2], q[3], w[0], w[1], w[2]);
};



void satellite3D_gyro_imdt_sys::get_output_function_blocks(
  satellite3D_gyro_imdt_sys::matrixC_type& C, 
  satellite3D_gyro_imdt_sys::matrixD_type& D, 
  const satellite3D_gyro_imdt_sys::state_space_type&, 
  const satellite3D_gyro_imdt_sys::time_type&, 
  const satellite3D_gyro_imdt_sys::point_type&, 
  const satellite3D_gyro_imdt_sys::input_type&) const {
  
  C = mat<double,mat_structure::nil>(9,12);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(6),3,6);
  
  D = mat<double,mat_structure::nil>(9,6);
};



satellite3D_gyro_imdt_sys::invariant_error_type satellite3D_gyro_imdt_sys::get_invariant_error(
  const satellite3D_gyro_imdt_sys::state_space_type&, 
  const satellite3D_gyro_imdt_sys::point_type& x, 
  const satellite3D_gyro_imdt_sys::input_type&, 
  const satellite3D_gyro_imdt_sys::output_type& y, 
  const satellite3D_gyro_imdt_sys::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(x)) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(x);
  vect<double,3> dw_IMU = q_diff.as_rotation() * vect<double,3>(y[7],y[8],y[9]) - get_ang_velocity(x);
  return satellite3D_gyro_imdt_sys::invariant_error_type(
    y[0] - pos[0],  y[1] - pos[1],  y[2] - pos[2],
    2.0 * a[0],     2.0 * a[1],     2.0 * a[2],
    dw_IMU[0],      dw_IMU[1],      dw_IMU[2]); 
};




satellite3D_IMU_imdt_sys::state_belief_type satellite3D_IMU_imdt_sys::get_zero_state_belief(double aCovValue = 10.0) const {
  
};

satellite3D_IMU_imdt_sys::output_belief_type satellite3D_IMU_imdt_sys::get_zero_output_belief(double aCovValue) const {
  return output_belief_type(output_type(vect_n<double>(0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)), 
                            covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(15,aCovValue))));
};


satellite3D_IMU_imdt_sys::satellite3D_IMU_imdt_sys(
  const std::string& aName, 
  double aMass,
  const mat<double,mat_structure::symmetric>& aInertiaMoment,
  double aDt,
  const unit_quat<double>& aIMUOrientation,
  const vect<double,3>& aIMULocation,
  const unit_quat<double>& aRoomOrientation,
  const vect<double,3>& aMagFieldVector, 
  std::size_t aApproxOrder) :
  satellite3D_imdt_sys(aName, aMass, aInertiaMoment, aDt, aApproxOrder),
  IMU_orientation(aIMUOrientation), IMU_location(aIMULocation),
  room_orientation(aRoomOrientation), mag_field_vector(aMagFieldVector) { }; 



satellite3D_IMU_imdt_sys::output_type satellite3D_IMU_imdt_sys::get_output(
  const satellite3D_IMU_imdt_sys::state_space_type&, 
  const satellite3D_IMU_imdt_sys::point_type& x, 
  const satellite3D_IMU_imdt_sys::input_type& u,
  const satellite3D_IMU_imdt_sys::time_type&) const {
  
  frame_3D<double> earth(
    shared_ptr< pose_3D<double> >(), 
    vect<double,3>(0.0,0.0,0.0),
    room_orientation.as_rotation(),
    vect<double,3>(0.0,0.0,0.0),vect<double,3>(0.0,0.0,0.0),
    vect<double,3>(0.0,0.0,9.81),vect<double,3>(0.0,0.0,0.0),
    vect<double,3>(0.0,0.0,0.0),vect<double,3>(0.0,0.0,0.0));
  
  shared_ptr< frame_3D<double> > earth_ptr(&earth, null_deleter());
      
  frame_3D<double> sat = get_frame_3D(x);
  sat.Parent = earth_ptr;
  shared_ptr< frame_3D<double> > sat_ptr(&sat, null_deleter());
  sat.Acceleration = vect<double,3>(u[0],u[1],u[2]) * (1.0 / mMass);
  sat.AngAcceleration = mInertiaMomentInv * vect<double,3>(u[3],u[4],u[5])
                        - sat.AngVelocity % (mInertiaMoment * sat.AngVelocity);
      
  frame_3D<double> IMU(
    sat_ptr, 
    IMU_location,
    IMU_orientation.as_rotation(),
    vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0),
    vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0),
    vect<double,3>(0.0,0.0,0.0), vect<double,3>(0.0,0.0,0.0));
  
  frame_3D<double> IMU_gbl = IMU.getGlobalFrame();
  
  vect<double,3> a_IMU = IMU_gbl.rotateFromParent(IMU_gbl.Acceleration);
  vect<double,3> m_IMU = IMU_gbl.rotateFromParent(mag_field_vector);
  
  return satellite3D_IMU_imdt_sys::output_type(
    sat.Position[0], sat.Position[1], sat.Position[2], 
    sat.Quat[0],     sat.Quat[1],     sat.Quat[2],     sat.Quat[3], 
    IMU_gbl.AngVelocity[0], IMU_gbl.AngVelocity[1], IMU_gbl.AngVelocity[2],
    a_IMU[0], a_IMU[1], a_IMU[2],
    m_IMU[0], m_IMU[1], m_IMU[2]);
};
    
void satellite3D_IMU_imdt_sys::get_output_function_blocks(
  satellite3D_IMU_imdt_sys::matrixC_type& C, 
  satellite3D_IMU_imdt_sys::matrixD_type& D, 
  const satellite3D_IMU_imdt_sys::state_space_type&, 
  const satellite3D_IMU_imdt_sys::time_type&, 
  const satellite3D_IMU_imdt_sys::point_type&, 
  const satellite3D_IMU_imdt_sys::input_type&) const {
  
  C = mat<double,mat_structure::nil>(15,12);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(6),3,6);
//   set_block(C,mat_ident<double>(3),6,3);
//   set_block(C,mat_ident<double>(3),9,9); //TODO
  
  D = mat<double,mat_structure::nil>(15,6); //TODO
};
        
satellite3D_IMU_imdt_sys::invariant_error_type satellite3D_IMU_imdt_sys::get_invariant_error(
  const satellite3D_IMU_imdt_sys::state_space_type&, 
  const satellite3D_IMU_imdt_sys::point_type& x, 
  const satellite3D_IMU_imdt_sys::input_type& u, 
  const satellite3D_IMU_imdt_sys::output_type& y, 
  const satellite3D_IMU_imdt_sys::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(x)) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(x);
  vect<double,3> dw_IMU = q_diff.as_rotation() * vect<double,3>(y[7],y[8],y[9]) - get_ang_velocity(x);
  return satellite3D_IMU_imdt_sys::invariant_error_type(
    y[0] - pos[0],  y[1] - pos[1],  y[2] - pos[2],
    2.0 * a[0],     2.0 * a[1],     2.0 * a[2],
    dw_IMU[0],      dw_IMU[1],      dw_IMU[2]);  //TODO
};





};

};

















