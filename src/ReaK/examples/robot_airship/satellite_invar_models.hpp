
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

#ifndef SATELLITE_INVAR_MODELS_HPP
#define SATELLITE_INVAR_MODELS_HPP

#include "base/named_object.hpp"

#include "lin_alg/mat_alg.hpp"
#include "ctrl_sys/sss_exceptions.hpp"
#include "topologies/se3_topologies.hpp"
#include "topologies/se2_topologies.hpp"

#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


  
  
class satellite2D_imdt_sys : public named_object {
  public:
    
    typedef pp::se2_1st_order_topology<double>::type state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::identity> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 6);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
    struct zero_input_trajectory {
      input_type get_point(time_type) const {
	return input_type(0.0,0.0,0.0);
      };
    };
    
  protected:
    double mMass;
    double mInertiaMoment;
    double mDt;
  
  public:
    
    satellite2D_imdt_sys(const std::string& aName = "", 
			 double aMass = 1.0, 
			 double aInertiaMoment = 1.0,
			 double aDt = 0.001) :
			 named_object(),
			 mMass(aMass),
			 mInertiaMoment(aInertiaMoment),
			 mDt(aDt) {
      setName(aName);
      if((mInertiaMoment < std::numeric_limits< double >::epsilon()) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information are singular in airship2D_lin_system's definition");
    };  
  
    virtual ~satellite2D_imdt_sys() { };
        
    time_difference_type get_time_step() const { return mDt; };
    
    point_type get_next_state(const state_space_type&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      vect<double,2> dv(mDt / mMass * u[0], mDt / mMass * u[1]);
      double dw = mDt / mInertiaMoment * u[2];
      double delta_theta = mDt * (x[6] + 0.5 * dw); 
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
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(get_position(x)[0], get_position(x)[1], get_rotation(x));
    };
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const state_space_type&, const time_type&, const point_type&, const input_type&) const {
      A = mat<double,mat_structure::identity>(6);
      A(0,3) = mDt / mMass;
      A(1,4) = mDt / mMass;  
      A(2,5) = mDt / mInertiaMoment;
      
      B = mat<double,mat_structure::nil>(6,3);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mInertiaMoment;
      B(3,0) = mDt;
      B(4,1) = mDt;
      B(5,2) = mDt;
    };
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, const time_type&, const point_type&, const input_type&) const {
      C = mat<double,mat_structure::nil>(3,6);
      set_block(C,mat<double,mat_structure::identity>(3),0,0);
      
      D = mat<double,mat_structure::nil>(3,3);
    };
        
    invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      rot_mat_2D<double> r_y( y[2] );
      
      return invariant_error_type(y[0] - get_position(x)[0],
			          y[1] - get_position(x)[1],
			          (invert(rot_mat_2D<double>(get_rotation(x))) * r_y).getAngle()); // s_y * c_x - c_y * s_x
    };
    
    point_type apply_correction(const state_space_type&, const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const {
      return make_arithmetic_tuple(
	make_arithmetic_tuple(
	  get_position(x) + vect<double,2>(c[0],c[1]),
	  get_velocity(x) + vect<double,2>(c[3] / mMass, c[4] / mMass)
	),
	make_arithmetic_tuple(
	  (rot_mat_2D<double>(get_rotation(x)) * rot_mat_2D<double>(c[2])).getAngle(),
	  get_ang_velocity(x) + c[5]
	)
      );
    };
    
    invariant_frame_type get_invariant_prior_frame(const state_space_type&, const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(6);
    };
    
    invariant_frame_type get_invariant_posterior_frame(const state_space_type&, const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(6);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
        & RK_SERIAL_SAVE_WITH_NAME(mDt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
        & RK_SERIAL_LOAD_WITH_NAME(mDt);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite2D_imdt_sys,0xC231000A,1,"satellite2D_imdt_sys",named_object)
    
};





class satellite3D_imdt_sys : public named_object {
  public:
    
    typedef pp::se3_1st_order_topology<double>::type state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::square> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
    struct zero_input_trajectory {
      input_type get_point(time_type) const {
	return input_type(0.0,0.0,0.0,0.0,0.0,0.0);
      };
    };
    
  private:
    double mMass;
    mat<double,mat_structure::symmetric> mInertiaMoment;
    mat<double,mat_structure::symmetric> mInertiaMomentInv;
    time_difference_type mDt;
    
  public:  
    satellite3D_imdt_sys(const std::string& aName = "", 
			 double aMass = 1.0, 
			 const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			 double aDt = 0.001) :
			 named_object(),
			 mMass(aMass),
			 mInertiaMoment(aInertiaMoment),
			 mDt(aDt) { 
      setName(aName);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information is improper in airship3D_lin_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
	throw system_incoherency("Inertial tensor is singular in airship3D_lin_system's definition");
      };
    }; 
  
    virtual ~satellite3D_imdt_sys() { };
    
    time_difference_type get_time_step() const { return mDt; };
    
    point_type get_next_state(const state_space_type&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
      
      vect<double,3> half_dp(0.005 * mDt * u[3], 0.005 * mDt * u[4], 0.005 * mDt * u[5]);
      vect<double,3> w0 = get_ang_velocity(x);
      unit_quat<double> half_w0_rot = exp( (0.0025 * mDt) * w0 );
      vect<double,3> dp0 = invert(half_w0_rot) * (mInertiaMoment * w0 + half_dp);
      
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
      return make_arithmetic_tuple(
	       make_arithmetic_tuple(get_position(x) + mDt * (get_velocity(x) + 0.5 * dv),
		                     get_velocity(x) + dv),
	       make_arithmetic_tuple(q_new, 
		                     w0));
    };
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      const vect<double,3>& pos = get_position(x);
      const unit_quat<double>& q = get_quaternion(x);
      return output_type(pos[0], pos[1], pos[2], q[0], q[1], q[2], q[3]);
    };
    
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const state_space_type&, 
				     const time_type&, const time_type&,
				     const point_type& p_0, const point_type& p_1,
				     const input_type&, const input_type&) const {
      
      mat<double, mat_structure::square> R = (invert( get_quaternion(p_1) ) * get_quaternion(p_0)).as_rotation().getMat();
      
      A = mat<double,mat_structure::identity>(12);
      A(0,6) = mDt;
      A(1,7) = mDt;  
      A(2,8) = mDt;
      set_block(A, (0.5 * mDt) * (mInertiaMomentInv
                                  + transpose_view(R) * mInertiaMomentInv * R), 3, 9);
            
      B = mat<double,mat_structure::nil>(12,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;
      set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 3, 3);
      B(6,0) = mDt / mMass;
      B(7,1) = mDt / mMass;
      B(8,2) = mDt / mMass;
      set_block(B, mDt * mat<double,mat_structure::identity>(3), 9, 3);
      
    };
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const {
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
    };
        
    invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      unit_quat<double> q_diff = invert(get_quaternion(x)) * unit_quat<double>(y[3],y[4],y[5],y[6]);
      vect<double,3> a = log(q_diff);
      const vect<double,3>& pos = get_position(x);
      return invariant_error_type(y[0] - pos[0],
			          y[1] - pos[1],
			          y[2] - pos[2],
	                          2.0 * a[0],
	                          2.0 * a[1],
	                          2.0 * a[2]); 
    };
    
    point_type apply_correction(const state_space_type&, const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const {
      unit_quat<double> q_diff = exp( 0.5 * vect<double,3>(c[3],c[4],c[5]) );
      unit_quat<double> q_new = get_quaternion(x) * q_diff;
      
      vect<double,3> w_new = mInertiaMomentInv * (invert(q_diff).as_rotation() * (mInertiaMoment * get_ang_velocity(x) + vect<double,3>(c[9],c[10],c[11])));
      return make_arithmetic_tuple(
	make_arithmetic_tuple(
	  get_position(x) + vect<double,3>(c[0],c[1],c[2]),
	  get_velocity(x) + vect<double,3>(c[6],c[7],c[8])
	),
	make_arithmetic_tuple(
	  q_new,
	  w_new
	)
      );
    };
    
    invariant_frame_type get_invariant_prior_frame(const state_space_type&, const point_type& x_prev, const point_type& x_prior, const input_type&, const time_type&) const {
      invariant_frame_type result(mat<double,mat_structure::identity>(12));
      mat<double,mat_structure::square> R_diff((invert(get_quaternion(x_prior)) * get_quaternion(x_prev)).as_rotation().getMat());
      set_block(result, R_diff, 3, 3);
      set_block(result, R_diff, 9, 9);
      return result;
    };
    
    invariant_frame_type get_invariant_posterior_frame(const state_space_type& state_space, const point_type& x_prior, const point_type& x_post, const input_type& u, const time_type& t) const {
      return get_invariant_prior_frame(state_space,x_prior,x_post,u,t);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
	& RK_SERIAL_SAVE_WITH_NAME(mDt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
	& RK_SERIAL_LOAD_WITH_NAME(mDt);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information is improper in airship3D_lin_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
	throw system_incoherency("Inertial tensor is singular in airship3D_lin_system's definition");
      };
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_imdt_sys,0xC231000B,1,"satellite3D_imdt_sys",named_object)
    
};



};

};

#endif

















