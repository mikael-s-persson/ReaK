
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

#include "math/mat_alg.hpp"
#include "ctrl_sys/sss_exceptions.hpp"
#include "rigid_body_states.hpp"

#include "math/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


  
  
class satellite2D_imdt_sys : public named_object {
  public:
    
    typedef rigid_body2D_state<double> point_type;
    typedef typename state_vector_traits< point_type >::state_difference_type point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::identity> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 6);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
    struct zero_input_trajectory {
      vect_n<double> get_point(time_type) const {
	return vect_n<double>(3,0.0);
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
    
    point_type get_next_state(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      input_type dp(mDt * u[0], mDt * u[1], mDt * u[2]);
      double delta_theta = (mDt / mInertiaMoment) * (x[6] + 0.5 * dp[2]); double cdt = cos(delta_theta); double sdt = sin(delta_theta);
      return point_type( x[0] + mDt / mMass * (x[4] + 0.5 * dp[0]),
			 x[1] + mDt / mMass * (x[5] + 0.5 * dp[1]),
			 x[2] * cdt - x[3] * sdt,
			 x[2] * sdt + x[3] * cdt,
			 x[4] + dp[0],
			 x[5] + dp[1],
			 x[6] + dp[2]);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x.position[0], x.position[1], x.rotation[0], x.rotation[1]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
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
      
      C = mat<double,mat_structure::nil>(3,6);
      set_block(C,mat<double,mat_structure::identity>(3),0,0);
      
      D = mat<double,mat_structure::nil>(3,3);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      rot_mat_2D<double> r_y( vect<double,2>(y[2],y[3]) );
      
      return invariant_error_type(y[0] - x.position[0],
			          y[1] - x.position[1],
			          (invert(x.rotation) * r_y).getAngle()); // s_y * c_x - c_y * s_x
    };
    
    point_type apply_correction(const point_type& x, const invariant_correction_type& c, const input_type& u, const time_type& t) const {
      
      return point_type(x.position + vect<double,2>(c[0],c[1]),
			x.rotation * rot_mat_2D<double>(c[2]),
			x.momentum + vect<double,2>(c[3],c[4]),
			x.angular_momentum + c[5]);
    };
    
    invariant_frame_type get_invariant_prior_frame(const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(6);
    };
    
    invariant_frame_type get_invariant_posterior_frame(const point_type&, const point_type&, const input_type&, const time_type&) const {
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
    
    typedef rigid_body3D_state<double> point_type;
    typedef typename state_vector_traits< point_type >::state_difference_type point_difference_type;
  
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
      vect_n<double> get_point(time_type) const {
	return vect_n<double>(6,0.0);
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
    
    point_type get_next_state(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
      
      vect<double,3> half_dh(0.5 * mDt * u[3], 0.5 * mDt * u[4], 0.5 * mDt * u[5]);
      vect<double,3> w0 = mInertiaMomentInv * x.angular_momentum;
      quaternion<double> half_w0_rot = quaternion<double>( exp( quat<double>( (0.25 * mDt) * w0) ) );
      vect<double,3> dh0 = invert(half_w0_rot) * (x.angular_momentum + half_dh);
      
      vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dh - mDt * w0 % x.angular_momentum));
      
      for(int i = 0; i < 20; ++i) {
	vect<double,3> w1_next = mInertiaMomentInv * (half_dh 
	                          + quaternion<double>( exp( quat<double>( (-0.25 * mDt) * w1_prev) ) ) * dh0);
	if(norm(w1_next - w1_prev) < 1E-6 * norm(w1_next + w1_prev)) {
	  w1_prev = w1_next;
	  break;
	} else
	  w1_prev = w1_next;
      };
      
      quaternion<double> q_new = x.rotation * 
                                 half_w0_rot * 
                                 quaternion<double>( exp( quat<double>( (0.25 * mDt) * w1_prev) ) );
				 
      vect<double,3> dp(mDt * u[0], mDt * u[1], mDt * u[2]);
      return point_type( x.position + (mDt / mMass) * (x.momentum + 0.5 * dp),
			 q_new,
			 x.momentum + dp,
			 mInertiaMoment * w1_prev);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
    };
    
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::identity>(12);
      A(0,6) = mDt / mMass;
      A(1,7) = mDt / mMass;  
      A(2,8) = mDt / mMass;
      set_block(A, mDt * mInertiaMomentInv, 3, 9);
            
      B = mat<double,mat_structure::nil>(12,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;
      set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 3, 3);
      B(6,0) = mDt;
      B(7,1) = mDt;
      B(8,2) = mDt;
      B(9,3) = mDt;
      B(10,4) = mDt;
      B(11,5) = mDt;
      
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      quaternion<double> q_diff = invert(x.rotation) 
                                * quaternion<double>(vect<double,4>(y[3],y[4],y[5],y[6]));
      quat<double> a = log(quat<double>(q_diff[0],q_diff[1],q_diff[2],q_diff[3]));
      return invariant_error_type(y[0] - x[0],
			          y[1] - x[1],
			          y[2] - x[2],
	                          2.0 * a[1],
	                          2.0 * a[2],
	                          2.0 * a[3]); 
    };
    
    point_type apply_correction(const point_type& x, const invariant_correction_type& c, const input_type& u, const time_type& t) const {
      quaternion<double> q_diff(exp(quat<double>(0.0, 0.5 * c[3], 0.5 * c[4], 0.5 * c[5])));
      quaternion<double> q_new = x.rotation * q_diff;
      vect<double,3> h_new = mInertiaMomentInv * (invert(q_diff) * (x.angular_momentum + vect<double,3>(c[9],c[10],c[11])));
      return point_type( x.position + vect<double,3>(c[0],c[1],c[2]),
			 q_new,
			 x.momentum + vect<double,3>(c[6],c[7],c[8]),
			 h_new);
    };
    
    invariant_frame_type get_invariant_prior_frame(const point_type& x_prev, const point_type& x_prior, const input_type&, const time_type&) const {
      invariant_frame_type result(mat<double,mat_structure::identity>(12));
      mat<double,mat_structure::square> R_diff((invert(x_prior.rotation) * x_prev.rotation).getMat());
      set_block(result, R_diff, 3, 3);
      set_block(result, R_diff, 9, 9);
      return result;
    };
    
    invariant_frame_type get_invariant_posterior_frame(const point_type& x_prior, const point_type& x_post, const input_type& u, const time_type& t) const {
      return get_invariant_prior_frame(x_prior,x_post,u,t);
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

















