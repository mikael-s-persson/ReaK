
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

#ifndef AIRSHIP3D_LIN_MODEL_HPP
#define AIRSHIP3D_LIN_MODEL_HPP

#include "lin_alg/vect_alg.hpp"
#include "base/named_object.hpp"

#include "lin_alg/mat_alg.hpp"
#include "ctrl_sys/sss_exceptions.hpp"

#include "lin_alg/mat_cholesky.hpp"

#include "kinetostatics/quat_alg.hpp"
#include "kinetostatics/rotations_3D.hpp"

namespace ReaK {


namespace ctrl {
  


class airship3D_lin_system : public named_object {
  public:
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
    typedef vect_n<double> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  protected:
    double mMass;
    mat<double,mat_structure::symmetric> mInertiaMoment;
    mat<double,mat_structure::symmetric> mInertiaMomentInv;
    
  public:
    
    airship3D_lin_system(const std::string& aName = "", 
			 double aMass = 1.0, 
			 const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))) :
			 named_object(),
			 mMass(aMass),
			 mInertiaMoment(aInertiaMoment) {
      setName(aName);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information is improper in airship3D_lin_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
	throw system_incoherency("Inertial tensor is singular in airship3D_lin_system's definition");
      };
    }; 
  
    virtual ~airship3D_lin_system() { };
    
    point_derivative_type get_state_derivative(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      quaternion<double> q(vect<double,4>(x[3],x[4],x[5],x[6]));
      vect<double,3> w(x[10],x[11],x[12]);
      vect<double,4> qd = q.getQuaternionDot(w);
      vect<double,3> aacc = mInertiaMomentInv * ( vect<double,3>(u[3],u[4],u[5]) - w % (mInertiaMoment * w) );
      return point_derivative_type(x[7],
	                           x[8],
				   x[9],
				   qd[0],
				   qd[1],
				   qd[2],
				   qd[3],
				   u[0] / mMass,
				   u[1] / mMass,
				   u[2] / mMass,
				   aacc[0],
				   aacc[1],
				   aacc[2]);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      vect<double,3> w(-x[10],-x[11],-x[12]);
      
      A = mat<double,mat_structure::nil>(13,13);
      A(0,7) = 1.0;
      A(1,8) = 1.0;
      A(2,9) = 1.0;

      w[0] = -0.5 * x[4];
      w[1] = -0.5 * x[5];
      w[2] = -0.5 * x[6];
      set_block(A, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >(w),3,10);
      set_block(A, (0.5 * x[3]) * mat<double,mat_structure::identity>(3) - mat<double,mat_structure::skew_symmetric>(w),4,10);

      B = mat<double,mat_structure::nil>(13,6);
      B(7,0) = 1.0 / mMass;
      B(8,1) = 1.0 / mMass;
      B(9,2) = 1.0 / mMass;
      set_block(B,mInertiaMomentInv,10,3);
      
      C = mat<double,mat_structure::nil>(7,13);
      set_block(C,mat<double,mat_structure::identity>(7),0,0);
      
      D = mat<double,mat_structure::nil>(7,6);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information is improper in airship3D_lin_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
	throw system_incoherency("Inertial tensor is singular in airship3D_lin_system's definition");
      };
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_lin_system,0xC2310005,1,"airship3D_lin_system",named_object)
    
};




class airship3D_inv_system : public airship3D_lin_system {
  public:
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
    typedef vect_n<double> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
    
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::identity> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
   
    airship3D_inv_system(const std::string& aName = "", 
			 double aMass = 1.0, 
			 const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))) :
			 airship3D_lin_system(aName,aMass,aInertiaMoment) { }; 
  
    virtual ~airship3D_inv_system() { };
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::nil>(12,12);
      A(0,6) = 1.0;
      A(1,7) = 1.0;
      A(2,8) = 1.0;
      A(3,9) = 1.0;
      A(4,10) = 1.0;
      A(5,11) = 1.0;
      
      B = mat<double,mat_structure::nil>(12,6);
      B(6,0) = 1.0 / mMass;
      B(7,1) = 1.0 / mMass;
      B(8,2) = 1.0 / mMass;
      set_block(B,mInertiaMomentInv,9,3);
      
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      quaternion<double> q_diff = quaternion<double>(vect<double,4>(y[3],y[4],y[5],y[6])) 
                                * invert(quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6])));
      return invariant_error_type(y[0] - x[0],
			          y[1] - x[1],
			          y[2] - x[2],
	                          2.0 * q_diff[0] * q_diff[1],
	                          2.0 * q_diff[0] * q_diff[2],
	                          2.0 * q_diff[0] * q_diff[3]); 
    };
    
    point_derivative_type apply_correction(const point_type& x, const point_derivative_type& xd, const invariant_correction_type& c, const input_type& u, const time_type& t) const {
      quaternion<double> q(vect<double,4>(x[3],x[4],x[5],x[6]));
      vect<double,4> dq = q.getQuaternionDot(vect<double,3>(c[3],c[4],c[5]));
      return point_type(xd[0] + c[0],
	                xd[1] + c[1],
			xd[2] + c[2],
			xd[3] + dq[0],
			xd[4] + dq[1],
			xd[5] + dq[2],
			xd[6] + dq[3],
			xd[7] + c[6],
			xd[8] + c[7],
			xd[9] + c[8],
			xd[10] + c[9],
			xd[11] + c[10],
			xd[12] + c[11]);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      airship3D_lin_system::save(A,airship3D_lin_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship3D_lin_system::load(A,airship3D_lin_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_inv_system,0xC2310006,1,"airship3D_inv_system",airship3D_lin_system)
    
};





class airship3D_lin_dt_system : public airship3D_lin_system {
  public:
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  protected:
    time_difference_type mDt;
    
  public:
    
    airship3D_lin_dt_system(const std::string& aName = "", 
			    double aMass = 1.0, 
			    const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			    double aDt = 0.001) :
			    airship3D_lin_system(aName,aMass,aInertiaMoment),
			    mDt(aDt) { 
      if(mDt < std::numeric_limits< double >::epsilon())
	throw system_incoherency("The time step is below numerical tolerance in airship2D_lin_dt_system's definition");
    }; 
  
    virtual ~airship3D_lin_dt_system() { };
    
    time_difference_type get_time_step() const { return mDt; };
    
    point_type get_next_state(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
      
      vect<double,3> half_dp(0.5 * mDt * u[3], 0.5 * mDt * u[4], 0.5 * mDt * u[5]);
      vect<double,3> w0(x[10],x[11],x[12]);
      quaternion<double> half_w0_rot = quaternion<double>( exp( quat<double>( (0.25 * mDt) * w0) ) );
      vect<double,3> dp0 = invert(half_w0_rot) * (mInertiaMoment * w0 + half_dp);
      
      vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dp - mDt * w0 % (mInertiaMoment * w0)));
      
      for(int i = 0; i < 20; ++i) {
	vect<double,3> w1_next = mInertiaMomentInv * (half_dp 
	                          + quaternion<double>( exp( quat<double>( (-0.25 * mDt) * w1_prev) ) ) * dp0);
	if(norm(w1_next - w1_prev) < 1E-6 * norm(w1_next + w1_prev)) {
	  w1_prev = w1_next;
	  break;
	} else
	  w1_prev = w1_next;
      };
      
      quaternion<double> q_new = quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6])) * 
                                 half_w0_rot * 
                                 quaternion<double>( exp( quat<double>( (0.25 * mDt) * w1_prev) ) );
				 
      vect<double,3> dv(mDt * u[0] / mMass, mDt * u[1] / mMass, mDt * u[2] / mMass);
      return point_type( x[0] + mDt * (x[7] + 0.5 * dv[0]),
			 x[1] + mDt * (x[8] + 0.5 * dv[1]),
			 x[2] + mDt * (x[9] + 0.5 * dv[2]),
			 q_new[0],
			 q_new[1],
			 q_new[2],
			 q_new[3],
			 x[7] + dv[0],
			 x[8] + dv[1],
			 x[9] + dv[2],
			 w1_prev[0],
			 w1_prev[1],
			 w1_prev[2]);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      vect<double,3> w(-mDt * x[10],-mDt * x[11],-mDt * x[12]);
      
      A = mat<double,mat_structure::identity>(13);
      A(0,7) = mDt;
      A(1,8) = mDt;  
      A(2,9) = mDt;
      mat<double,mat_structure::square> JinvWJ(mInertiaMomentInv * mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment);
      set_block(A, mat<double,mat_structure::identity>(3) + JinvWJ, 10, 10);
      
      w[0] = -0.5 * mDt * x[4];
      w[1] = -0.5 * mDt * x[5];
      w[2] = -0.5 * mDt * x[6];
      JinvWJ = mat<double,mat_structure::identity>(3) + 0.5 * JinvWJ;
      vect<double,3> w_JinvWJ = w * JinvWJ;
      set_block(A, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >(w_JinvWJ),3,10);
      set_block(A, ((0.5 * mDt * x[3]) * mat<double,mat_structure::identity>(3) - mat<double,mat_structure::skew_symmetric>(w)) * JinvWJ,4,10);
      
      B = mat<double,mat_structure::nil>(13,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;

      w *= mDt * 0.5;
      vect<double,3> w_Jinv = w * mInertiaMomentInv;
      set_block(B, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >( w_Jinv ), 3, 3);
      set_block(B, (0.25 * mDt * mDt * x[3]) * mInertiaMomentInv - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMomentInv, 4, 3);
      
      B(7,0) = mDt / mMass;
      B(8,1) = mDt / mMass;
      B(9,2) = mDt / mMass;
      set_block(B, mDt * mInertiaMomentInv, 10, 3);
      
      C = mat<double,mat_structure::nil>(7,13);
      set_block(C,mat<double,mat_structure::identity>(7),0,0);
      
      D = mat<double,mat_structure::nil>(7,6);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      airship3D_lin_system::save(A,airship3D_lin_system::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mDt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship3D_lin_system::load(A,airship3D_lin_system::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mDt);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_lin_dt_system,0xC2310007,1,"airship3D_lin_dt_system",airship3D_lin_system)
    
};





class airship3D_inv_dt_system : public airship3D_lin_dt_system {
  public:
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::identity> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
    
    airship3D_inv_dt_system(const std::string& aName = "", 
			    double aMass = 1.0, 
			    const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			    double aDt = 0.001) :
			    airship3D_lin_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 
  
    virtual ~airship3D_inv_dt_system() { };
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      vect<double,3> w(x[10],x[11],x[12]);
      
      A = mat<double,mat_structure::identity>(12);
      A(0,6) = mDt;
      A(1,7) = mDt;  
      A(2,8) = mDt;
      
      mat<double,mat_structure::square> T(mInertiaMomentInv * (mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment
                                                            - mat<double,mat_structure::skew_symmetric>(mInertiaMoment * w)));
      
      set_block(A, mat<double,mat_structure::diagonal>(3,mDt) - (0.5 * mDt * mDt) * T, 3, 9);
      set_block(A, mat<double,mat_structure::identity>(3) - mDt * T, 9, 9);
            
      B = mat<double,mat_structure::nil>(12,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;
      set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 3, 3);
      B(6,0) = mDt / mMass;
      B(7,1) = mDt / mMass;
      B(8,2) = mDt / mMass;
      set_block(B, mDt * mInertiaMomentInv, 9, 3);
      
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      quaternion<double> q_diff = invert(quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6]))) 
                                * quaternion<double>(vect<double,4>(y[3],y[4],y[5],y[6]));
      quat<double> a = log(quat<double>(q_diff[0],q_diff[1],q_diff[2],q_diff[3]));
      return invariant_error_type(y[0] - x[0],
			          y[1] - x[1],
			          y[2] - x[2],
	                          2.0 * a[1],
	                          2.0 * a[2],
	                          2.0 * a[3]); 
    };
    
    point_type apply_correction(const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const {
      quaternion<double> q_new = quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6])) * 
                                 quaternion<double>(exp(quat<double>(0.0, 0.5 * c[3], 0.5 * c[4], 0.5 * c[5])));
      return point_type(x[0] + c[0],
	                x[1] + c[1],
			x[2] + c[2],
			q_new[0],
			q_new[1],
			q_new[2],
			q_new[3],
			x[7] + c[6],
			x[8] + c[7],
			x[9] + c[8],
			x[10] + c[9],
			x[11] + c[10],
			x[12] + c[11]);
    };
    
    invariant_frame_type get_invariant_prior_frame(const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(invariant_correction_dimensions);
    };
    
    invariant_frame_type get_invariant_posterior_frame(const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(invariant_correction_dimensions);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      airship3D_lin_dt_system::save(A,airship3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship3D_lin_dt_system::load(A,airship3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_inv_dt_system,0xC2310008,1,"airship3D_inv_dt_system",airship3D_lin_dt_system)
    
};





class airship3D_inv_mom_dt_system : public airship3D_lin_dt_system {
  public:
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
  
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
    
    
    airship3D_inv_mom_dt_system(const std::string& aName = "", 
			        double aMass = 1.0, 
			        const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			        double aDt = 0.001) :
			        airship3D_lin_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 
  
    virtual ~airship3D_inv_mom_dt_system() { };
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      vect<double,3> w(x[10],x[11],x[12]);
      
      A = mat<double,mat_structure::identity>(12);
      A(0,6) = mDt;
      A(1,7) = mDt;  
      A(2,8) = mDt;
      set_block(A, mDt * mInertiaMomentInv, 3, 9);
            
      B = mat<double,mat_structure::nil>(12,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;
      set_block(B, (0.5 * mDt * mDt) * mInertiaMomentInv, 3, 3);
      B(6,0) = mDt / mMass;
      B(7,1) = mDt / mMass;
      B(8,2) = mDt / mMass;
      B(9,3) = mDt;
      B(10,4) = mDt;
      B(11,5) = mDt;
      
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      quaternion<double> q_diff = invert(quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6]))) 
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
      quaternion<double> q_new = quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6])) * 
                                 q_diff;
      vect<double,3> w_new = mInertiaMomentInv * (invert(q_diff) * (mInertiaMoment * vect<double,3>(x[10],x[11],x[12]) + vect<double,3>(c[9],c[10],c[11])));
      return point_type(x[0] + c[0],
	                x[1] + c[1],
			x[2] + c[2],
			q_new[0],
			q_new[1],
			q_new[2],
			q_new[3],
			x[7] + c[6],
			x[8] + c[7],
			x[9] + c[8],
			w_new[0],
			w_new[1],
			w_new[2]);
    };
    
    invariant_frame_type get_invariant_prior_frame(const point_type& x_prev, const point_type& x_prior, const input_type&, const time_type&) const {
      invariant_frame_type result(mat<double,mat_structure::identity>(12));
      mat<double,mat_structure::square> R_diff((invert(quaternion<double>(vect<double,4>(x_prior[3],x_prior[4],x_prior[5],x_prior[6])))
                                                 * quaternion<double>(vect<double,4>(x_prev[3],x_prev[4],x_prev[5],x_prev[6]))).getMat());
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
      airship3D_lin_dt_system::save(A,airship3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship3D_lin_dt_system::load(A,airship3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_inv_mom_dt_system,0xC2310009,1,"airship3D_inv_mom_dt_system",airship3D_lin_dt_system)
    
};






};

};

#endif









