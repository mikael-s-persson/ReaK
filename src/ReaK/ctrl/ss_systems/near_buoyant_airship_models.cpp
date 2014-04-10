
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

#include "near_buoyant_airship_models.hpp"

#include "ctrl_sys/sss_exceptions.hpp"
#include "topologies/se3_topologies.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


#define RK_D_INF std::numeric_limits<double>::infinity()

shared_ptr< airship3D_imdt_em_sys::temporal_state_space_type > airship3D_imdt_em_sys::get_temporal_state_space(double aStartTime, double aEndTime) const {
  return shared_ptr< temporal_state_space_type >(new temporal_state_space_type(
    "airship3D_em_temporal_space", 
    state_space_type(make_arithmetic_tuple(
      pp::make_se3_space(
        "satellite3D_state_space",
        vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
        vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
        RK_D_INF, RK_D_INF),
      pp::line_segment_topology<double>("mass_imbal_param_space", 0.0, RK_D_INF),
      pp::hyperball_topology< vect<double,3> >("eccentricity_param_space", vect<double,3>(0.0,0.0,0.0), RK_D_INF)
    )),
    pp::time_poisson_topology("airship3D_em_time_space", mDt, (aEndTime - aStartTime) * 0.5)));
};

shared_ptr< airship3D_imdt_em_sys::state_space_type > airship3D_imdt_em_sys::get_state_space() const {
  return shared_ptr< state_space_type >(new state_space_type(make_arithmetic_tuple(
    pp::make_se3_space(
      "satellite3D_state_space",
      vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
      vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
      RK_D_INF, RK_D_INF),
    pp::line_segment_topology<double>("mass_imbal_param_space", 0.0, RK_D_INF),
    pp::hyperball_topology< vect<double,3> >("eccentricity_param_space", vect<double,3>(0.0,0.0,0.0), RK_D_INF)
  )));
};

#undef RK_D_INF


airship3D_imdt_em_sys::state_belief_type airship3D_imdt_em_sys::get_zero_state_belief(double aCovValue) const {
  point_type x_init;
  set_frame_3D(get<0>(x_init), frame_3D<double>());
  get<1>(x_init) = 0.0;
  get<2>(x_init) = vect<double,3>(0.0, 0.0, 0.0);
  return state_belief_type(x_init, covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(16, aCovValue))));
};

airship3D_imdt_em_sys::input_belief_type airship3D_imdt_em_sys::get_zero_input_belief(double aCovValue) const {
  return input_belief_type(input_type(vect_n<double>(6, 0.0)), 
                           covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(6,aCovValue))));
};

airship3D_imdt_em_sys::output_belief_type airship3D_imdt_em_sys::get_zero_output_belief(double aCovValue) const {
  return output_belief_type(output_type(vect_n<double>(0.0,0.0,0.0,1.0,0.0,0.0,0.0)), 
                            covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(6,aCovValue))));
};



airship3D_imdt_em_sys::airship3D_imdt_em_sys(
  const std::string& aName, double aMass, 
  const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt, const vect<double,3>& aGravityAcc) :
  named_object(), mMass(aMass), mInertiaMoment(aInertiaMoment), mDt(aDt), mGravityAcc(aGravityAcc) { 
  setName(aName);
  if(mDt < std::numeric_limits< double >::epsilon())
    throw system_incoherency("The time step is below numerical tolerance in airship3D_imdt_em_sys's definition");
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in airship3D_imdt_em_sys's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in airship3D_imdt_em_sys's definition");
  };
}; 



airship3D_imdt_em_sys::point_type airship3D_imdt_em_sys::get_next_state(
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::point_type& x, 
  const airship3D_imdt_em_sys::input_type& u, 
  const airship3D_imdt_em_sys::time_type&) const {
  //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
  
  typedef pp::se3_1st_order_topology<double>::type SE3StateSpace;
  typedef pp::topology_traits< SE3StateSpace >::point_type SE3StateType;
  
  const SE3StateType& x_se3 = get<0>(x);
  
  vect<double,3> half_dp(0.05 * mDt * u[3], 0.05 * mDt * u[4], 0.05 * mDt * u[5]);
  vect<double,3> w0 = get_ang_velocity(x_se3);
  unit_quat<double> half_w0_rot = exp( (0.025 * mDt) * w0 );
  vect<double,3> dp0 = invert(half_w0_rot).as_rotation() * (mInertiaMoment * w0 + half_dp);
  
  unit_quat<double> q_new = get_quaternion(x_se3);
  vect<double,3> g_torque = (0.01 * mDt * mMass) * (get<2>(x) % (invert(q_new).as_rotation() * mGravityAcc));
  
  for(unsigned int i = 0; i < 10; ++i) {
    
    vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dp + g_torque - (0.1 * mDt) * w0 % (mInertiaMoment * w0)));
    unit_quat<double> half_w1_prev_rot = exp( (0.025 * mDt) * w1_prev );
    
    for(int i = 0; i < 20; ++i) {
      vect<double,3> w1_next = mInertiaMomentInv * (half_dp + g_torque 
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
    
    g_torque = (0.01 * mDt * mMass) * (get<2>(x) % (invert(q_new).as_rotation() * mGravityAcc));
  };
  
  vect<double,3> dv = (mDt / mMass) * (get<1>(x) * mGravityAcc + u[range(0,2)]);
  return airship3D_imdt_em_sys::point_type(
    SE3StateType(
      make_arithmetic_tuple(
        get_position(x_se3) + mDt * (get_velocity(x_se3) + 0.5 * dv),
        get_velocity(x_se3) + dv),
      make_arithmetic_tuple(
        q_new, 
        w0)
    ),
    get<1>(x),
    get<2>(x)
  );
};


airship3D_imdt_em_sys::output_type airship3D_imdt_em_sys::get_output(
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::point_type& x, 
  const airship3D_imdt_em_sys::input_type&, 
  const airship3D_imdt_em_sys::time_type&) const {
  const vect<double,3>& pos = get_position(get<0>(x));
  const unit_quat<double>& q = get_quaternion(get<0>(x));
  return airship3D_imdt_em_sys::output_type(pos[0], pos[1], pos[2], q[0], q[1], q[2], q[3]);
};


void airship3D_imdt_em_sys::get_state_transition_blocks(
  airship3D_imdt_em_sys::matrixA_type& A, 
  airship3D_imdt_em_sys::matrixB_type& B, 
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::time_type&, 
  const airship3D_imdt_em_sys::time_type&,
  const airship3D_imdt_em_sys::point_type& p_0, 
  const airship3D_imdt_em_sys::point_type& p_1,
  const airship3D_imdt_em_sys::input_type&, 
  const airship3D_imdt_em_sys::input_type&) const {
  
  // all states conserved:
  A = mat_ident<double>(16);
  
  // p-v block:
  A(0,3) = mDt;
  A(1,4) = mDt;
  A(2,5) = mDt;
  
  // q-w block:
  A(6, 9) = mDt; 
  A(7,10) = mDt; 
  A(8,11) = mDt; 
  
  // p-m block:
  A(0,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[0];
  A(1,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[1];
  A(2,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[2];
  // v-m block:
  A(3,12) = mDt / mMass * mGravityAcc[0];
  A(4,12) = mDt / mMass * mGravityAcc[1];
  A(5,12) = mDt / mMass * mGravityAcc[2];
  
  // q-r and w-r blocks:
  mat<double, mat_structure::square> R = get_quaternion(get<0>(p_0)).as_rotation().getMat();
  mat<double, mat_structure::square> JRgR(
    mInertiaMomentInv * transpose_view(R) * mat<double,mat_structure::skew_symmetric>(mGravityAcc) * R);
  set_block(A, (-0.5 * mDt * mDt * mMass) * JRgR, 6, 13);
  set_block(A, (-mDt * mMass) * JRgR, 9, 13);
  
  
  B = mat<double,mat_structure::nil>(16,6);
  // (p,v)-f block:
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,2) = 0.5 * mDt * mDt / mMass;
  B(3,0) = mDt / mMass;
  B(4,1) = mDt / mMass;
  B(5,2) = mDt / mMass;
  // (q,w)-t block:
  set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 6, 3);
  set_block(B, mDt * mInertiaMomentInv, 9, 3);
  
};


void airship3D_imdt_em_sys::get_output_function_blocks(
  airship3D_imdt_em_sys::matrixC_type& C, 
  airship3D_imdt_em_sys::matrixD_type& D, 
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::time_type&, 
  const airship3D_imdt_em_sys::point_type&, 
  const airship3D_imdt_em_sys::input_type&) const {
  C = mat<double,mat_structure::nil>(6,16);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(3),3,6);
  
  D = mat<double,mat_structure::nil>(6,6);
};


airship3D_imdt_em_sys::invariant_error_type airship3D_imdt_em_sys::get_invariant_error(
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::point_type& x, 
  const airship3D_imdt_em_sys::input_type&, 
  const airship3D_imdt_em_sys::output_type& y, 
  const airship3D_imdt_em_sys::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(get<0>(x))) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(get<0>(x));
  
  return airship3D_imdt_em_sys::invariant_error_type(
    y[0] - pos[0], y[1] - pos[1], y[2] - pos[2],
    2.0 * a[0], 2.0 * a[1], 2.0 * a[2]); 
};


airship3D_imdt_em_sys::point_type airship3D_imdt_em_sys::apply_correction(
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::point_type& x, 
  const airship3D_imdt_em_sys::invariant_correction_type& c, 
  const airship3D_imdt_em_sys::input_type&, 
  const airship3D_imdt_em_sys::time_type&) const {
  
  typedef pp::se3_1st_order_topology<double>::type SE3StateSpace;
  typedef pp::topology_traits< SE3StateSpace >::point_type SE3StateType;
  
  const SE3StateType& x_se3 = get<0>(x);
  
  unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[6],c[7],c[8]) );
  unit_quat<double> q_new = get_quaternion(x_se3) * q_diff;
  
  vect<double,3> w_new = mInertiaMomentInv * (invert(q_diff).as_rotation() * (mInertiaMoment * get_ang_velocity(x_se3) + vect<double,3>(c[9],c[10],c[11])));
  return airship3D_imdt_em_sys::point_type(
    SE3StateType(
      make_arithmetic_tuple(
        get_position(x_se3) + vect<double,3>(c[0],c[1],c[2]),
        get_velocity(x_se3) + vect<double,3>(c[3],c[4],c[5])
      ),
      make_arithmetic_tuple(
        q_new, 
        w_new
      )
    ),
    get<1>(x) + c[12],
    get<2>(x) + vect<double,3>(c[13],c[14],c[15])
  );
};



airship3D_imdt_em_sys::invariant_frame_type airship3D_imdt_em_sys::get_invariant_prior_frame(
  const airship3D_imdt_em_sys::state_space_type&, 
  const airship3D_imdt_em_sys::point_type& x_prev, 
  const airship3D_imdt_em_sys::point_type& x_prior, 
  const airship3D_imdt_em_sys::input_type&, 
  const airship3D_imdt_em_sys::time_type&) const {
  
  airship3D_imdt_em_sys::invariant_frame_type result(mat<double,mat_structure::identity>(16));
  mat<double,mat_structure::square> R_diff((invert(get_quaternion(get<0>(x_prior))) * get_quaternion(get<0>(x_prev))).as_rotation().getMat());
  set_block(result, R_diff, 6, 6);
  set_block(result, R_diff, 9, 9);
  return result;
};



void RK_CALL airship3D_imdt_em_sys::save(ReaK::serialization::oarchive& A, unsigned int) const {
  named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mMass)
    & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_SAVE_WITH_NAME(mDt)
    & RK_SERIAL_SAVE_WITH_NAME(mGravityAcc);
};

void RK_CALL airship3D_imdt_em_sys::load(ReaK::serialization::iarchive& A, unsigned int) {
  named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mMass)
    & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_LOAD_WITH_NAME(mDt)
    & RK_SERIAL_LOAD_WITH_NAME(mGravityAcc);
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in airship3D_imdt_em_sys's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in airship3D_imdt_em_sys's definition");
  };
};






#define RK_D_INF std::numeric_limits<double>::infinity()

shared_ptr< airship3D_imdt_emd_sys::temporal_state_space_type > airship3D_imdt_emd_sys::get_temporal_state_space(double aStartTime, double aEndTime) const {
  return shared_ptr< temporal_state_space_type >(new temporal_state_space_type(
    "airship3D_emd_temporal_space", 
    state_space_type(make_arithmetic_tuple(
      pp::make_se3_space(
        "satellite3D_state_space",
        vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
        vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
        RK_D_INF, RK_D_INF),
      pp::line_segment_topology<double>("mass_imbal_param_space", 0.0, RK_D_INF),
      pp::hyperball_topology< vect<double,3> >("eccentricity_param_space", vect<double,3>(0.0,0.0,0.0), RK_D_INF),
      pp::line_segment_topology<double>("tr_drag_param_space", 0.0, RK_D_INF),
      pp::line_segment_topology<double>("rot_drag_param_space", 0.0, RK_D_INF)
    )),
    pp::time_poisson_topology("airship3D_emd_time_space", mDt, (aEndTime - aStartTime) * 0.5)));
};

shared_ptr< airship3D_imdt_emd_sys::state_space_type > airship3D_imdt_emd_sys::get_state_space() const {
  return shared_ptr< state_space_type >(new state_space_type(make_arithmetic_tuple(
    pp::make_se3_space(
      "satellite3D_state_space",
      vect<double,3>(-RK_D_INF, -RK_D_INF, -RK_D_INF),
      vect<double,3>( RK_D_INF,  RK_D_INF,  RK_D_INF),
      RK_D_INF, RK_D_INF),
    pp::line_segment_topology<double>("mass_imbal_param_space", 0.0, RK_D_INF),
    pp::hyperball_topology< vect<double,3> >("eccentricity_param_space", vect<double,3>(0.0,0.0,0.0), RK_D_INF),
    pp::line_segment_topology<double>("tr_drag_param_space", 0.0, RK_D_INF),
    pp::line_segment_topology<double>("rot_drag_param_space", 0.0, RK_D_INF)
  )));
};

#undef RK_D_INF


airship3D_imdt_emd_sys::state_belief_type airship3D_imdt_emd_sys::get_zero_state_belief(double aCovValue) const {
  point_type x_init;
  set_frame_3D(get<0>(x_init), frame_3D<double>());
  get<1>(x_init) = 0.0;
  get<2>(x_init) = vect<double,3>(0.0, 0.0, 0.0);
  get<3>(x_init) = 0.0;
  get<4>(x_init) = 0.0;
  return state_belief_type(x_init, covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(18, aCovValue))));
};

airship3D_imdt_emd_sys::input_belief_type airship3D_imdt_emd_sys::get_zero_input_belief(double aCovValue) const {
  return input_belief_type(input_type(vect_n<double>(6, 0.0)), 
                           covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(6,aCovValue))));
};

airship3D_imdt_emd_sys::output_belief_type airship3D_imdt_emd_sys::get_zero_output_belief(double aCovValue) const {
  return output_belief_type(output_type(vect_n<double>(0.0,0.0,0.0,1.0,0.0,0.0,0.0)), 
                            covar_type(covar_type::matrix_type(mat<double,mat_structure::diagonal>(6,aCovValue))));
};



airship3D_imdt_emd_sys::airship3D_imdt_emd_sys(
  const std::string& aName, double aMass, 
  const mat<double,mat_structure::symmetric>& aInertiaMoment, double aDt, const vect<double,3>& aGravityAcc) :
  named_object(), mMass(aMass), mInertiaMoment(aInertiaMoment), mDt(aDt), mGravityAcc(aGravityAcc) { 
  setName(aName);
  if(mDt < std::numeric_limits< double >::epsilon())
    throw system_incoherency("The time step is below numerical tolerance in airship3D_imdt_emd_sys's definition");
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in airship3D_imdt_emd_sys's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in airship3D_imdt_emd_sys's definition");
  };
}; 



airship3D_imdt_emd_sys::point_type airship3D_imdt_emd_sys::get_next_state(
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::point_type& x, 
  const airship3D_imdt_emd_sys::input_type& u, 
  const airship3D_imdt_emd_sys::time_type&) const {
  //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
  
  typedef pp::se3_1st_order_topology<double>::type SE3StateSpace;
  typedef pp::topology_traits< SE3StateSpace >::point_type SE3StateType;
  
  const SE3StateType& x_se3 = get<0>(x);
  
  vect<double,3> half_dp(0.05 * mDt * u[3], 0.05 * mDt * u[4], 0.05 * mDt * u[5]);
  vect<double,3> w0 = get_ang_velocity(x_se3);
  unit_quat<double> half_w0_rot = exp( (0.025 * mDt) * w0 );
  vect<double,3> half_dw0_drag = (-0.05 * mDt * get<4>(x) * norm_2(w0)) * w0;
  vect<double,3> dp0 = invert(half_w0_rot).as_rotation() * (mInertiaMoment * w0 + half_dp + half_dw0_drag);
  
  unit_quat<double> q_new = get_quaternion(x_se3);
  vect<double,3> g_torque = (0.01 * mDt * mMass) * (get<2>(x) % (invert(q_new).as_rotation() * mGravityAcc));
  
  for(unsigned int i = 0; i < 10; ++i) {
    
    vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * (half_dp + half_dw0_drag) + g_torque - (0.1 * mDt) * w0 % (mInertiaMoment * w0)));
    unit_quat<double> half_w1_prev_rot = exp( (0.025 * mDt) * w1_prev );
    vect<double,3> half_dw1_drag = (-0.05 * mDt * get<4>(x) * norm_2(w1_prev)) * w1_prev;
    
    for(int i = 0; i < 20; ++i) {
      vect<double,3> w1_next = mInertiaMomentInv * (half_dp + half_dw1_drag + g_torque 
                                + invert(half_w1_prev_rot).as_rotation() * dp0);
      half_dw1_drag = (-0.05 * mDt * get<4>(x) * norm_2(w1_next)) * w1_next;
      if(norm_2(w1_next - w1_prev) < 1E-6 * norm_2(w1_next + w1_prev)) {
        w1_prev = w1_next;
        break;
      } else
        w1_prev = w1_next;
    };
    
    q_new = q_new * half_w0_rot * half_w1_prev_rot;
    w0 = w1_prev;
    half_w0_rot = half_w1_prev_rot;
    half_dw0_drag = half_dw1_drag;
    dp0 = invert(half_w0_rot).as_rotation() * (mInertiaMoment * w0 + half_dp + half_dw0_drag);
    
    g_torque = (0.01 * mDt * mMass) * (get<2>(x) % (invert(q_new).as_rotation() * mGravityAcc));
  };
  
  const vect<double,3>& v = get_velocity(x_se3);
  vect<double,3> dv = (mDt / mMass) * (get<1>(x) * mGravityAcc - get<3>(x) * norm_2(v) * v + u[range(0,2)]);
  return airship3D_imdt_emd_sys::point_type(
    SE3StateType(
      make_arithmetic_tuple(
        get_position(x_se3) + mDt * (v + 0.5 * dv),
        v + dv),
      make_arithmetic_tuple(
        q_new, 
        w0)
    ),
    get<1>(x), get<2>(x), get<3>(x), get<4>(x)
  );
};


airship3D_imdt_emd_sys::output_type airship3D_imdt_emd_sys::get_output(
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::point_type& x, 
  const airship3D_imdt_emd_sys::input_type&, 
  const airship3D_imdt_emd_sys::time_type&) const {
  const vect<double,3>& pos = get_position(get<0>(x));
  const unit_quat<double>& q = get_quaternion(get<0>(x));
  return airship3D_imdt_emd_sys::output_type(pos[0], pos[1], pos[2], q[0], q[1], q[2], q[3]);
};


void airship3D_imdt_emd_sys::get_state_transition_blocks(
  airship3D_imdt_emd_sys::matrixA_type& A, 
  airship3D_imdt_emd_sys::matrixB_type& B, 
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::time_type&, 
  const airship3D_imdt_emd_sys::time_type&,
  const airship3D_imdt_emd_sys::point_type& p_0, 
  const airship3D_imdt_emd_sys::point_type& p_1,
  const airship3D_imdt_emd_sys::input_type&, 
  const airship3D_imdt_emd_sys::input_type&) const {
  
  // all states conserved:
  A = mat_ident<double>(18);
  
  // p-v block:
  A(0,3) = mDt;
  A(1,4) = mDt;
  A(2,5) = mDt;
  
  // q-w block:
  A(6, 9) = mDt; 
  A(7,10) = mDt; 
  A(8,11) = mDt; 
  
  // p-m block:
  A(0,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[0];
  A(1,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[1];
  A(2,12) = 0.5 * mDt * mDt / mMass * mGravityAcc[2];
  // v-m block:
  A(3,12) = mDt / mMass * mGravityAcc[0];
  A(4,12) = mDt / mMass * mGravityAcc[1];
  A(5,12) = mDt / mMass * mGravityAcc[2];
  
  // (p,v)-d block:
  const vect<double,3>& v = get_velocity(get<0>(p_0));
  double v_mag = norm_2(v);
  A(0,16) = -0.5 * mDt * mDt / mMass * v_mag * v[0];
  A(1,16) = -0.5 * mDt * mDt / mMass * v_mag * v[1];
  A(2,16) = -0.5 * mDt * mDt / mMass * v_mag * v[2];
  A(3,16) = -mDt / mMass * v_mag * v[0];
  A(4,16) = -mDt / mMass * v_mag * v[1];
  A(5,16) = -mDt / mMass * v_mag * v[2];
  
  
  // q-r and w-r blocks:
  mat<double, mat_structure::square> R = get_quaternion(get<0>(p_0)).as_rotation().getMat();
  mat<double, mat_structure::square> JRgR(
    mInertiaMomentInv * transpose_view(R) * mat<double,mat_structure::skew_symmetric>(mGravityAcc) * R);
  set_block(A, (-0.5 * mDt * mDt * mMass) * JRgR, 6, 13);
  set_block(A, (-mDt * mMass) * JRgR, 9, 13);
  
  // (q,w)-d block:
  const vect<double,3>& w = get_ang_velocity(get<0>(p_0));
  double w_mag = norm_2(w);
  vect<double,3> Jw = mInertiaMomentInv * w;
  A(6, 17) = -0.5 * mDt * mDt * w_mag * Jw[0];
  A(7, 17) = -0.5 * mDt * mDt * w_mag * Jw[1];
  A(8, 17) = -0.5 * mDt * mDt * w_mag * Jw[2];
  A(9, 17) = -mDt * w_mag * Jw[0];
  A(10,17) = -mDt * w_mag * Jw[1];
  A(11,17) = -mDt * w_mag * Jw[2];
  
  
  B = mat<double,mat_structure::nil>(18,6);
  // (p,v)-f block:
  B(0,0) = 0.5 * mDt * mDt / mMass;
  B(1,1) = 0.5 * mDt * mDt / mMass;
  B(2,2) = 0.5 * mDt * mDt / mMass;
  B(3,0) = mDt / mMass;
  B(4,1) = mDt / mMass;
  B(5,2) = mDt / mMass;
  // (q,w)-t block:
  set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 6, 3);
  set_block(B, mDt * mInertiaMomentInv, 9, 3);
  
};


void airship3D_imdt_emd_sys::get_output_function_blocks(
  airship3D_imdt_emd_sys::matrixC_type& C, 
  airship3D_imdt_emd_sys::matrixD_type& D, 
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::time_type&, 
  const airship3D_imdt_emd_sys::point_type&, 
  const airship3D_imdt_emd_sys::input_type&) const {
  C = mat<double,mat_structure::nil>(6,18);
  set_block(C,mat_ident<double>(3),0,0);
  set_block(C,mat_ident<double>(3),3,6);
  
  D = mat<double,mat_structure::nil>(6,6);
};


airship3D_imdt_emd_sys::invariant_error_type airship3D_imdt_emd_sys::get_invariant_error(
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::point_type& x, 
  const airship3D_imdt_emd_sys::input_type&, 
  const airship3D_imdt_emd_sys::output_type& y, 
  const airship3D_imdt_emd_sys::time_type&) const {
  
  unit_quat<double> q_diff = invert(get_quaternion(get<0>(x))) * unit_quat<double>(y[3],y[4],y[5],y[6]);
  vect<double,3> a = log(q_diff);
  const vect<double,3>& pos = get_position(get<0>(x));
  
  return airship3D_imdt_emd_sys::invariant_error_type(
    y[0] - pos[0], y[1] - pos[1], y[2] - pos[2],
    2.0 * a[0], 2.0 * a[1], 2.0 * a[2]); 
};


airship3D_imdt_emd_sys::point_type airship3D_imdt_emd_sys::apply_correction(
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::point_type& x, 
  const airship3D_imdt_emd_sys::invariant_correction_type& c, 
  const airship3D_imdt_emd_sys::input_type&, 
  const airship3D_imdt_emd_sys::time_type&) const {
  
  typedef pp::se3_1st_order_topology<double>::type SE3StateSpace;
  typedef pp::topology_traits< SE3StateSpace >::point_type SE3StateType;
  
  const SE3StateType& x_se3 = get<0>(x);
  
  unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[6],c[7],c[8]) );
  unit_quat<double> q_new = get_quaternion(x_se3) * q_diff;
  
  vect<double,3> w_new = mInertiaMomentInv * (invert(q_diff).as_rotation() * (mInertiaMoment * get_ang_velocity(x_se3) + vect<double,3>(c[9],c[10],c[11])));
  return airship3D_imdt_emd_sys::point_type(
    SE3StateType(
      make_arithmetic_tuple(
        get_position(x_se3) + vect<double,3>(c[0],c[1],c[2]),
        get_velocity(x_se3) + vect<double,3>(c[3],c[4],c[5])
      ),
      make_arithmetic_tuple(
        q_new, 
        w_new
      )
    ),
    get<1>(x) + c[12],
    get<2>(x) + vect<double,3>(c[13],c[14],c[15]),
    get<3>(x) + c[16], get<4>(x) + c[17]
  );
};



airship3D_imdt_emd_sys::invariant_frame_type airship3D_imdt_emd_sys::get_invariant_prior_frame(
  const airship3D_imdt_emd_sys::state_space_type&, 
  const airship3D_imdt_emd_sys::point_type& x_prev, 
  const airship3D_imdt_emd_sys::point_type& x_prior, 
  const airship3D_imdt_emd_sys::input_type&, 
  const airship3D_imdt_emd_sys::time_type&) const {
  
  airship3D_imdt_emd_sys::invariant_frame_type result(mat<double,mat_structure::identity>(18));
  mat<double,mat_structure::square> R_diff((invert(get_quaternion(get<0>(x_prior))) * get_quaternion(get<0>(x_prev))).as_rotation().getMat());
  set_block(result, R_diff, 6, 6);
  set_block(result, R_diff, 9, 9);
  return result;
};



void RK_CALL airship3D_imdt_emd_sys::save(ReaK::serialization::oarchive& A, unsigned int) const {
  named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mMass)
    & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_SAVE_WITH_NAME(mDt)
    & RK_SERIAL_SAVE_WITH_NAME(mGravityAcc);
};

void RK_CALL airship3D_imdt_emd_sys::load(ReaK::serialization::iarchive& A, unsigned int) {
  named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mMass)
    & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
    & RK_SERIAL_LOAD_WITH_NAME(mDt)
    & RK_SERIAL_LOAD_WITH_NAME(mGravityAcc);
  if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
    throw system_incoherency("Inertial information is improper in airship3D_imdt_emd_sys's definition");
  try {
    invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
  } catch(singularity_error&) {
    throw system_incoherency("Inertial tensor is singular in airship3D_imdt_emd_sys's definition");
  };
};






};

};

















