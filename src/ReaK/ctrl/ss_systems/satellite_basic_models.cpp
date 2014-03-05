
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

#include "satellite_basic_models.hpp"

#include "ctrl_sys/sss_exceptions.hpp"

#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {



satellite3D_lin_dt_system::satellite3D_lin_dt_system(
  const std::string& aName, double aMass, const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt) :
  named_object(), mMass(aMass), mInertiaMoment(aInertiaMoment), mDt(aDt) {
  setName(aName);
  if(mDt < std::numeric_limits< double >::epsilon())
    throw system_incoherency("The time step is below numerical tolerance in satellite3D_lin_dt_system's definition");
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in satellite3D_lin_dt_system's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in satellite3D_lin_dt_system's definition");
  };
}; 
  

satellite3D_lin_dt_system::point_type satellite3D_lin_dt_system::get_next_state(
  const satellite3D_lin_dt_system::state_space_type&, 
  const satellite3D_lin_dt_system::point_type& x, 
  const satellite3D_lin_dt_system::input_type& u, 
  const satellite3D_lin_dt_system::time_type&) const {
  //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
  
  vect<double,3> half_dp(0.005 * mDt * u[3], 0.005 * mDt * u[4], 0.005 * mDt * u[5]);
  vect<double,3> w0 = get_ang_velocity(x);
  unit_quat<double> half_w0_rot = exp( (0.0025 * mDt) * w0 );
  vect<double,3> dp0 = invert(half_w0_rot).as_rotation() * (mInertiaMoment * w0 + half_dp);
  
  unit_quat<double> q_new = get_quaternion(x);
  
  for(unsigned int i = 0; i < 100; ++i) {
    
    vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dp - (0.01 * mDt) * w0 % (mInertiaMoment * w0)));
    unit_quat<double> half_w1_prev_rot = exp( (0.0025 * mDt) * w1_prev );
    
    for(int i = 0; i < 20; ++i) {
      vect<double,3> w1_next = mInertiaMomentInv * (half_dp 
                                + invert(half_w1_prev_rot).as_rotation() * dp0);
      if(norm_2(w1_next - w1_prev) < 1E-6 * norm_2(w1_next + w1_prev)) {
        w1_prev = w1_next;
        break;
      } else
        w1_prev = w1_next;
    };
      
    q_new = q_new * half_w0_rot * half_w1_prev_rot;
    w0 = w1_prev;
    half_w0_rot = half_w1_prev_rot;
  };
      
  vect<double,3> dv(mDt * u[0] / mMass, mDt * u[1] / mMass, mDt * u[2] / mMass);
  return satellite3D_lin_dt_system::point_type(
    make_arithmetic_tuple(
      get_position(x) + mDt * (get_velocity(x) + 0.5 * dv),
      get_velocity(x) + dv),
    make_arithmetic_tuple(
      q_new, 
      w0));
};



void satellite3D_lin_dt_system::get_state_transition_blocks(
  satellite3D_lin_dt_system::matrixA_type& A, 
  satellite3D_lin_dt_system::matrixB_type& B, 
  const satellite3D_lin_dt_system::state_space_type&, 
  const satellite3D_lin_dt_system::time_type&, 
  const satellite3D_lin_dt_system::time_type&,
  const satellite3D_lin_dt_system::point_type& p_0, 
  const satellite3D_lin_dt_system::point_type&,
  const satellite3D_lin_dt_system::input_type& u_0, 
  const satellite3D_lin_dt_system::input_type&) const {
  
  vect<double,3> w = -mDt * get_ang_velocity(p_0);
  
  A = mat_ident<double>(12);
  A(0,3) = mDt;
  A(1,4) = mDt;  
  A(2,5) = mDt;
  
  mat<double,mat_structure::square> JinvWJ(mInertiaMomentInv * mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment);
  set_block(A, mat_ident<double>(3) + JinvWJ, 9, 9);
  
  w -= 0.5 * (mInertiaMomentInv * ( (mDt * mDt) * vect<double,3>(u_0[3], u_0[4], u_0[5]) - w % (mInertiaMoment * w) ) );
  
  set_block(A, mat_ident<double>(3) + mat<double,mat_structure::skew_symmetric>(w), 6, 6);
  /* This part is wrong (for q-q, not a-a). */
//   w *= 0.5;
//   set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::row_major>(w), 3, 4);
//   set_block(A, mat_ident<double>(3) + mat<double,mat_structure::skew_symmetric>(w), 4, 4);
//   w *= -1.0;
//   set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::column_major>(w), 4, 3);
  /* . */
  
  JinvWJ = mat_ident<double>(3) + 0.5 * JinvWJ;
//   const unit_quat<double>& q = get_quaternion(p_0);
  
  set_block(A, mDt * JinvWJ, 6, 9);
  /* This part is wrong (for q-w, not a-w). */
//   w[0] = -0.5 * mDt * q[1];
//   w[1] = -0.5 * mDt * q[2];
//   w[2] = -0.5 * mDt * q[3];
//   vect<double,3> w_JinvWJ = w * JinvWJ;
//   set_block(A, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >(w_JinvWJ),3,10);
//   set_block(A, ((0.5 * mDt * q[0]) * mat_ident<double>(3) - mat<double,mat_structure::skew_symmetric>(w)) * JinvWJ,4,10);
  /* . */
  
  
  B = mat<double,mat_structure::nil>(12,6);
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,2) = 0.5 * mDt * mDt / mMass;
  B(3,0) = mDt / mMass;
  B(4,1) = mDt / mMass;
  B(5,2) = mDt / mMass;
  
  set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 6, 3);
  /* This part is wrong (for q-t, not a-t). */
//   w = -0.25 * mDt * mDt * vect<double,3>(q[1], q[2], q[3]);
//   vect<double,3> w_Jinv = w * mInertiaMomentInv;
//   set_block(B, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >( w_Jinv ), 3, 3);
//   set_block(B, (0.25 * mDt * mDt * q[0]) * mInertiaMomentInv - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMomentInv, 4, 3);
  /* . */
  
  set_block(B, mDt * mInertiaMomentInv, 9, 3);
  
};



satellite3D_lin_dt_system::output_type satellite3D_lin_dt_system::get_output(
  const satellite3D_lin_dt_system::state_space_type&, 
  const satellite3D_lin_dt_system::point_type& x, 
  const satellite3D_lin_dt_system::input_type&, 
  const satellite3D_lin_dt_system::time_type&) const {
  const vect<double,3>& pos = get_position(x);
  const unit_quat<double>& q = get_quaternion(x);
  return satellite3D_lin_dt_system::output_type(pos[0], pos[1], pos[2], q[0], q[1], q[2], q[3]);
};



void satellite3D_lin_dt_system::get_output_function_blocks(
  satellite3D_lin_dt_system::matrixC_type& C, 
  satellite3D_lin_dt_system::matrixD_type& D, 
  const satellite3D_lin_dt_system::state_space_type&,
  const satellite3D_lin_dt_system::time_type&, 
  const satellite3D_lin_dt_system::point_type&, 
  const satellite3D_lin_dt_system::input_type&) const {
  C = mat<double,mat_structure::nil>(7,13);
  set_block(C,mat<double,mat_structure::identity>(7),0,0);
  
  D = mat<double,mat_structure::nil>(7,6);
};

void RK_CALL satellite3D_lin_dt_system::save(ReaK::serialization::oarchive& A, unsigned int) const {
  named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mMass)
    & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_SAVE_WITH_NAME(mDt);
};

void RK_CALL satellite3D_lin_dt_system::load(ReaK::serialization::iarchive& A, unsigned int) {
  named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mMass)
    & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_LOAD_WITH_NAME(mDt);
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in satellite3D_lin_dt_system's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in satellite3D_lin_dt_system's definition");
  };
};







satellite3D_gyro_lin_dt_system::satellite3D_gyro_lin_dt_system(
  const std::string& aName, double aMass, const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt) :
  satellite3D_lin_dt_system(aName, aMass, aInertiaMoment, aDt) { }; 
  
    
satellite3D_gyro_lin_dt_system::output_type satellite3D_gyro_lin_dt_system::get_output(
  const satellite3D_gyro_lin_dt_system::state_space_type&, 
  const satellite3D_gyro_lin_dt_system::point_type& x, 
  const satellite3D_gyro_lin_dt_system::input_type&, 
  const satellite3D_gyro_lin_dt_system::time_type&) const {
  const vect<double,3>& pos = get_position(x);
  const unit_quat<double>& q = get_quaternion(x);
  const vect<double,3>& w = get_ang_velocity(x);
  return satellite3D_gyro_lin_dt_system::output_type(
    pos[0], pos[1], pos[2], 
    q[0], q[1], q[2], q[3], 
    w[0], w[1], w[2]);
};
    
void satellite3D_gyro_lin_dt_system::get_output_function_blocks(
  satellite3D_gyro_lin_dt_system::matrixC_type& C, 
  satellite3D_gyro_lin_dt_system::matrixD_type& D, 
  const satellite3D_gyro_lin_dt_system::state_space_type&, 
  const satellite3D_gyro_lin_dt_system::time_type&, 
  const satellite3D_gyro_lin_dt_system::point_type&, 
  const satellite3D_gyro_lin_dt_system::input_type&) const {
  C = mat<double,mat_structure::nil>(10,13);
  set_block(C,mat<double,mat_structure::identity>(7),0,0);
  set_block(C,mat<double,mat_structure::identity>(3),7,3);
  
  D = mat<double,mat_structure::nil>(14,6);
};











satellite3D_inv_dt_system::satellite3D_inv_dt_system(
  const std::string& aName, double aMass, const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt) :
  satellite3D_lin_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 



void satellite3D_inv_dt_system::get_state_transition_blocks(
  satellite3D_inv_dt_system::matrixA_type& A, 
  satellite3D_inv_dt_system::matrixB_type& B, 
  const satellite3D_inv_dt_system::state_space_type&, 
  const satellite3D_inv_dt_system::time_type&, 
  const satellite3D_inv_dt_system::time_type&,
  const satellite3D_inv_dt_system::point_type& p_0, 
  const satellite3D_inv_dt_system::point_type&,
  const satellite3D_inv_dt_system::input_type&, 
  const satellite3D_inv_dt_system::input_type&) const {
  
  vect<double,3> w = 0.5 * get_ang_velocity(p_0);
  
  A = mat_ident<double>(12);
  A(0,3) = mDt;
  A(1,4) = mDt;
  A(2,5) = mDt;
  
  mat<double,mat_structure::square> T(mInertiaMomentInv * (mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment
                                                        - mat<double,mat_structure::skew_symmetric>(mInertiaMoment * w)));
  
  set_block(A, mat<double,mat_structure::diagonal>(3,mDt) - (0.5 * mDt * mDt) * T, 6, 9);
  set_block(A, mat_ident<double>(3) - mDt * T, 9, 9);
  
  w *= 2.0 * mDt;
  set_block(A, mat_ident<double>(3) - mat<double,mat_structure::skew_symmetric>(w), 6, 6);
  
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



void satellite3D_inv_dt_system::get_output_function_blocks(
  satellite3D_inv_dt_system::matrixC_type& C, 
  satellite3D_inv_dt_system::matrixD_type& D, 
  const satellite3D_inv_dt_system::state_space_type&, 
  const satellite3D_inv_dt_system::time_type&, 
  const satellite3D_inv_dt_system::point_type&, 
  const satellite3D_inv_dt_system::input_type&) const {
  C = mat<double,mat_structure::nil>(6,12);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(3),3,6);
  
  D = mat<double,mat_structure::nil>(6,6);
};



satellite3D_inv_dt_system::invariant_error_type satellite3D_inv_dt_system::get_invariant_error(
  const satellite3D_inv_dt_system::state_space_type&, 
  const satellite3D_inv_dt_system::point_type& x, 
  const satellite3D_inv_dt_system::input_type&, 
  const satellite3D_inv_dt_system::output_type& y, 
  const satellite3D_inv_dt_system::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(x)) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(x);
  return satellite3D_inv_dt_system::invariant_error_type(
    y[0] - pos[0], y[1] - pos[1], y[2] - pos[2],
    2.0 * a[0], 2.0 * a[1], 2.0 * a[2]); 
};



satellite3D_inv_dt_system::point_type satellite3D_inv_dt_system::apply_correction(
  const satellite3D_inv_dt_system::state_space_type&, 
  const satellite3D_inv_dt_system::point_type& x, 
  const satellite3D_inv_dt_system::invariant_correction_type& c, 
  const satellite3D_inv_dt_system::input_type&, 
  const satellite3D_inv_dt_system::time_type&) const {
  
  unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[6],c[7],c[8]) );
  unit_quat<double> q_new = get_quaternion(x) * q_diff;
  
  return satellite3D_inv_dt_system::point_type(
    make_arithmetic_tuple(
      get_position(x) + vect<double,3>(c[0],c[1],c[2]),
      get_velocity(x) + vect<double,3>(c[3],c[4],c[5])
    ),
    make_arithmetic_tuple(
      q_new,
      get_ang_velocity(x) + vect<double,3>(c[9],c[10],c[11])
    )
  );
};





satellite3D_gyro_inv_dt_system::satellite3D_gyro_inv_dt_system(
  const std::string& aName, double aMass, const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt) :
  satellite3D_inv_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 
  
satellite3D_gyro_inv_dt_system::output_type satellite3D_gyro_inv_dt_system::get_output(
  const satellite3D_gyro_inv_dt_system::state_space_type&, 
  const satellite3D_gyro_inv_dt_system::point_type& x, 
  const satellite3D_gyro_inv_dt_system::input_type&, 
  const satellite3D_gyro_inv_dt_system::time_type&) const {
  
  const vect<double,3>& pos = get_position(x);
  const unit_quat<double>& q = get_quaternion(x);
  const vect<double,3>& w = get_ang_velocity(x);
  return satellite3D_gyro_inv_dt_system::output_type(
    pos[0], pos[1], pos[2], 
    q[0], q[1], q[2], q[3], 
    w[0], w[1], w[2]);
};

void satellite3D_gyro_inv_dt_system::get_output_function_blocks(
  satellite3D_gyro_inv_dt_system::matrixC_type& C, 
  satellite3D_gyro_inv_dt_system::matrixD_type& D, 
  const satellite3D_gyro_inv_dt_system::state_space_type&, 
  const satellite3D_gyro_inv_dt_system::time_type&, 
  const satellite3D_gyro_inv_dt_system::point_type&, 
  const satellite3D_gyro_inv_dt_system::input_type&) const {
  
  C = mat<double,mat_structure::nil>(9,12);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(3),3,6);
  set_block(C,mat_ident<double>(3),6,9);
  
  D = mat<double,mat_structure::nil>(12,6);
};
    
satellite3D_gyro_inv_dt_system::invariant_error_type satellite3D_gyro_inv_dt_system::get_invariant_error(
  const satellite3D_gyro_inv_dt_system::state_space_type&, 
  const satellite3D_gyro_inv_dt_system::point_type& x, 
  const satellite3D_gyro_inv_dt_system::input_type&, 
  const satellite3D_gyro_inv_dt_system::output_type& y, 
  const satellite3D_gyro_inv_dt_system::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(x)) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(x);
  vect<double,3> dw_IMU = q_diff.as_rotation() * vect<double,3>(y[7],y[8],y[9]) - get_ang_velocity(x);
  return satellite3D_gyro_inv_dt_system::invariant_error_type(
    y[0] - pos[0], y[1] - pos[1], y[2] - pos[2],
    2.0 * a[0], 2.0 * a[1], 2.0 * a[2],
    dw_IMU[0], dw_IMU[1], dw_IMU[2]); 
};










};

};









