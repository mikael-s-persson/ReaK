/**
 * \file quadrotor_system.hpp
 * 
 * This library implements a state-space system for a quad-rotor aircraft.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */


/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef RK_QUADROTOR_SYSTEM_HPP
#define RK_QUADROTOR_SYSTEM_HPP

#include "lin_alg/vect_alg.hpp"
#include "base/named_object.hpp"

#include "lin_alg/mat_alg.hpp"
#include "ctrl_sys/sss_exceptions.hpp"
#include "topologies/se3_topologies.hpp"

#include "lin_alg/mat_cholesky.hpp"

#include "kinetostatics/quat_alg.hpp"
#include "kinetostatics/rotations_3D.hpp"

#include "topologies/vector_topology.hpp"

namespace ReaK {


namespace ctrl {
  


class quadrotor_system : public named_object {
  public:
    
    typedef pp::se3_1st_order_topology<double>::type state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    typedef point_difference_type point_derivative_type;
    
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect<double,4> input_type;
    typedef point_type output_type;
    
    typedef point_difference_type invariant_error_type;
    typedef point_difference_type invariant_correction_type;
    typedef mat<double,mat_structure::identity> invariant_frame_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 4);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 12);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  protected:
    double mMass;
    mat<double,mat_structure::symmetric> mInertiaMoment;
    mat<double,mat_structure::symmetric> mInertiaMomentInv;
    
    mat<double,mat_structure::diagonal> mTransDragCoefs;
    mat<double,mat_structure::diagonal> mRotDragCoefs;
    
  public:
    
    quadrotor_system(const std::string& aName = "", 
                     double aMass = 1.0, 
                     const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
                     const mat<double,mat_structure::diagonal>& aTransDragCoefs = (mat<double,mat_structure::diagonal>(3, 0.5)),
                     const mat<double,mat_structure::diagonal>& aRotDragCoefs = (mat<double,mat_structure::diagonal>(3, 0.5))
                    ) :
                     named_object(),
                     mMass(aMass),
                     mInertiaMoment(aInertiaMoment),
                     mTransDragCoefs(aTransDragCoefs),
                     mRotDragCoefs(aRotDragCoefs) {
      setName(aName);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
        throw system_incoherency("Inertial information is improper in airship3D_lin_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
        throw system_incoherency("Inertial tensor is singular in airship3D_lin_system's definition");
      };
    }; 
  
    virtual ~quadrotor_system() { };
    
    point_derivative_type get_state_derivative(const state_space_type&, const point_type& x, const input_type& u, time_type = 0.0) const {
      using std::fabs;
      
      quaternion<double> q = get_quaternion(x).as_rotation();
      vect<double,3> w = get_ang_velocity(x);
      vect<double,3> w_sqr( w[0] * fabs(w[0]), w[1] * fabs(w[1]), w[2] * fabs(w[2]) );
      vect<double,3> aacc = mInertiaMomentInv * ( vect<double,3>(u[1],u[2],u[3]) - w % (mInertiaMoment * w) - mRotDragCoefs * w_sqr );
      
      vect<double,3> local_v = invert(q) * get_velocity(x);
      vect<double,3> v = get_velocity(x);
      local_v[0] *= fabs(local_v[0]) / mMass;
      local_v[1] *= fabs(local_v[1]) / mMass;
      local_v[2] *= fabs(local_v[2]) / mMass;
      return point_derivative_type(
        make_arithmetic_tuple(
          v,  // velocity -> derivative of position
          vect<double,3>(0.0, 0.0, 9.81) - q * (mTransDragCoefs * local_v + vect<double,3>(0.0, 0.0, u[0] / mMass))
        ),  
        make_arithmetic_tuple(
          w,  // angular velocity -> invariant derivative of rotation
          aacc
        )
      );
    };
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type&, const time_type t = 0.0) const {
      return x;
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const state_space_type&, const time_type& t, const point_type& x, const input_type& u) const {
      using std::fabs;
      
      quaternion<double> q = get_quaternion(x).as_rotation();
      mat<double, mat_structure::square> R = q.getMat();
      
      A = mat<double,mat_structure::nil>(12,12);
      
      // velocity to position partial derivative:
      A(0,3) = 1.0;
      A(1,4) = 1.0;
      A(2,5) = 1.0;
      
      vect<double,3> local_v = invert(q) * get_velocity(x);
      mat<double,mat_structure::diagonal> dV(vect<double,3>(
        -2.0 * fabs(local_v[0]) / mMass, 
        -2.0 * fabs(local_v[1]) / mMass, 
        -2.0 * fabs(local_v[2]) / mMass));
      
      // velocity - velocity partial derivative:
      set_block(A, R * mTransDragCoefs * dV, 3, 3);
      
      // velocity - quaternion partial derivative:
      local_v[0] *= fabs(local_v[0]) / mMass;
      local_v[1] *= fabs(local_v[1]) / mMass;
      local_v[2] *= fabs(local_v[2]) / mMass;
      set_block(A, R * (mat<double,mat_structure::skew_symmetric>(mTransDragCoefs * local_v) 
                      - mTransDragCoefs * mat<double,mat_structure::skew_symmetric>(local_v) 
                      + mat<double,mat_structure::skew_symmetric>(vect<double,3>(0.0, 0.0, u[0] / mMass))), 
                3, 6);
//       set_block(A, R * (mat<double,mat_structure::skew_symmetric>(vect<double,3>(0.0, 0.0, u[0] / mMass))), 
//                 3, 6);
//       set_block(A, mat<double,mat_structure::skew_symmetric>(vect<double,3>(0.0, 0.0, 9.81)), 
//                 3, 6);
      
      // angular velocity to quaternion partial derivative:
      A(6,9)  = 1.0;
      A(7,10) = 1.0;  // identity
      A(8,11) = 1.0;
      
      vect<double,3> w = get_ang_velocity(x);
      mat<double,mat_structure::diagonal> dW(vect<double,3>(
        -2.0 * fabs(w[0]), 
        -2.0 * fabs(w[1]), 
        -2.0 * fabs(w[2])));
//       set_block(A, mInertiaMomentInv * ( mat<double,mat_structure::skew_symmetric>(mInertiaMoment * w)
//                                        - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment), 
//                 9, 9);
      set_block(A, mInertiaMomentInv * ( mRotDragCoefs * dW 
                                       - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment
                                       + mat<double,mat_structure::skew_symmetric>(mInertiaMoment * w) ), 
                9, 9);
      
      B = mat<double,mat_structure::nil>(12,4);
      vect<double,3> t_global = R * vect<double,3>(0.0, 0.0, -1.0 / mMass);
//       vect<double,3> t_global = vect<double,3>(0.0, 0.0, -1.0 / mMass);
      set_block(B, mat_vect_adaptor< vect<double,3> >(t_global), 3, 0);
      set_block(B, mInertiaMomentInv, 9, 1);
      
      C = mat<double,mat_structure::nil>(6,12);
      set_block(C,mat<double,mat_structure::identity>(6),0,0);
      
      D = mat<double,mat_structure::nil>(6,6);
      
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mMass)
        & RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment)
        & RK_SERIAL_SAVE_WITH_NAME(mTransDragCoefs)
        & RK_SERIAL_SAVE_WITH_NAME(mRotDragCoefs);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mMass)
        & RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment)
        & RK_SERIAL_LOAD_WITH_NAME(mTransDragCoefs)
        & RK_SERIAL_LOAD_WITH_NAME(mRotDragCoefs);
      if((mInertiaMoment.get_row_count() != 3) || (mMass < std::numeric_limits< double >::epsilon()))
        throw system_incoherency("Inertial information is improper in quadrotor_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
        throw system_incoherency("Inertial tensor is singular in quadrotor_system's definition");
      };
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(quadrotor_system,0xC2310005,1,"quadrotor_system",named_object)
    
};



#if 0

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
                            const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
                            double aDt = 0.001) :
                            airship3D_lin_system(aName,aMass,aInertiaMoment),
                            mDt(aDt) { 
      if(mDt < std::numeric_limits< double >::epsilon())
        throw system_incoherency("The time step is below numerical tolerance in airship2D_lin_dt_system's definition");
    }; 
  
    virtual ~airship3D_lin_dt_system() { };
    
    time_difference_type get_time_step() const { return mDt; };
    
    void set_time_step(time_difference_type aDt) { mDt = aDt; };
    
    point_type get_next_state(const pp::vector_topology< vect_n<double> >&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      //this function implements the momentum-conserving trapezoidal rule (variational integrator). This is very similar to the symplectic variational midpoint integrator over Lie Groups.
      
      /*
      vect<double,3> half_dp(0.5 * mDt * u[3], 0.5 * mDt * u[4], 0.5 * mDt * u[5]);
      vect<double,3> w0(x[10],x[11],x[12]);
      quaternion<double> half_w0_rot = quaternion<double>( exp( quat<double>( (0.25 * mDt) * w0) ) );
      vect<double,3> dp0 = invert(half_w0_rot) * (mInertiaMoment * w0 + half_dp);
      
      vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dp - mDt * w0 % (mInertiaMoment * w0)));
      
      for(int i = 0; i < 20; ++i) {
        vect<double,3> w1_next = mInertiaMomentInv * (half_dp 
                                  + quaternion<double>( exp( quat<double>( (-0.25 * mDt) * w1_prev) ) ) * dp0);
        if(norm_2(w1_next - w1_prev) < 1E-6 * norm_2(w1_next + w1_prev)) {
          w1_prev = w1_next;
          break;
        } else
          w1_prev = w1_next;
      };
      
      quaternion<double> q_new = quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6])) * 
                                 half_w0_rot * 
                                 quaternion<double>( exp( quat<double>( (0.25 * mDt) * w1_prev) ) );
      */
      
      vect<double,3> half_dp(0.005 * mDt * u[3], 0.005 * mDt * u[4], 0.005 * mDt * u[5]);
      vect<double,3> w0(x[10],x[11],x[12]);
      quaternion<double> half_w0_rot = quaternion<double>( exp( quat<double>( (0.0025 * mDt) * w0) ) );
      vect<double,3> dp0 = invert(half_w0_rot) * (mInertiaMoment * w0 + half_dp);
      
      quaternion<double> q_new = quaternion<double>(vect<double,4>(x[3],x[4],x[5],x[6]));
      
      for(unsigned int i = 0; i < 100; ++i) {
      
        vect<double,3> w1_prev = w0 + (mInertiaMomentInv * (2.0 * half_dp - (0.01 * mDt) * w0 % (mInertiaMoment * w0)));
        quaternion<double> half_w1_prev_rot = quaternion<double>( exp( quat<double>( (0.0025 * mDt) * w1_prev) ) );
        
        for(int i = 0; i < 20; ++i) {
          vect<double,3> w1_next = mInertiaMomentInv * (half_dp 
                                    + invert(half_w1_prev_rot) * dp0);
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
                         w0[0],
                         w0[1],
                         w0[2]);
    };
    
    output_type get_output(const pp::vector_topology< vect_n<double> >&, const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
    };
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const pp::vector_topology< vect_n<double> >&, 
                                     const time_type& t_0, const time_type&,
                                     const point_type& p_0, const point_type&,
                                     const input_type& u_0, const input_type&) const {
      vect<double,3> w(-mDt * p_0[10],-mDt * p_0[11],-mDt * p_0[12]);
      
      A = mat<double,mat_structure::identity>(13);
      A(0,7) = mDt;
      A(1,8) = mDt;  
      A(2,9) = mDt;
      mat<double,mat_structure::square> JinvWJ(mInertiaMomentInv * mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment);
      set_block(A, mat<double,mat_structure::identity>(3) + JinvWJ, 10, 10);
      
      w *= 0.5;
      set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::row_major>(w), 3, 4);
      set_block(A, mat<double,mat_structure::identity>(3) + mat<double,mat_structure::skew_symmetric>(w), 4, 4);
      w *= -1.0;
      set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::column_major>(w), 4, 3);
      
      w[0] = -0.5 * mDt * p_0[4];
      w[1] = -0.5 * mDt * p_0[5];
      w[2] = -0.5 * mDt * p_0[6];
      set_block(A, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >(w),3,10);
      set_block(A, (0.5 * mDt * p_0[3]) * mat<double,mat_structure::identity>(3) - mat<double,mat_structure::skew_symmetric>(w),4,10);
      
      B = mat<double,mat_structure::nil>(13,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;

      w *= mDt * 0.5;
      vect<double,3> w_Jinv = w * mInertiaMomentInv;
      set_block(B, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >( w_Jinv ), 3, 3);
      set_block(B, (0.25 * mDt * mDt * p_0[3]) * mInertiaMomentInv - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMomentInv, 4, 3);
      
      B(7,0) = mDt / mMass;
      B(8,1) = mDt / mMass;
      B(9,2) = mDt / mMass;
      set_block(B, mDt * mInertiaMomentInv, 10, 3);
      
    };
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const pp::vector_topology< vect_n<double> >&, 
                                    const time_type&, const point_type&, const input_type&) const {
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



class airship3D_lin2_dt_system : public airship3D_lin_dt_system {
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
    
    
    airship3D_lin2_dt_system(const std::string& aName = "", 
                             double aMass = 1.0, 
                             const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
                             double aDt = 0.001) :
                             airship3D_lin_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 
  
    virtual ~airship3D_lin2_dt_system() { };
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const pp::vector_topology< vect_n<double> >&, 
                                     const time_type& t_0, const time_type&,
                                     const point_type& p_0, const point_type&,
                                     const input_type& u_0, const input_type&) const {
      vect<double,3> w(-mDt * p_0[10],-mDt * p_0[11],-mDt * p_0[12]);
      
      A = mat<double,mat_structure::identity>(13);
      A(0,7) = mDt;
      A(1,8) = mDt;  
      A(2,9) = mDt;
      mat<double,mat_structure::square> JinvWJ(mInertiaMomentInv * mat<double,mat_structure::skew_symmetric>(w) * mInertiaMoment);
      set_block(A, mat<double,mat_structure::identity>(3) + JinvWJ, 10, 10);
      w -= 0.5 * (mInertiaMomentInv * ( (mDt * mDt) * vect<double,3>(u_0[3], u_0[4], u_0[5]) - w % (mInertiaMoment * w) ) );
      w *= 0.5;
      set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::row_major>(w), 3, 4);
      set_block(A, mat<double,mat_structure::identity>(3) + mat<double,mat_structure::skew_symmetric>(w), 4, 4);
      w *= -1.0;
      set_block(A, mat_vect_adaptor<vect<double,3>,mat_alignment::column_major>(w), 4, 3);
      
      
      w[0] = -0.5 * mDt * p_0[4];
      w[1] = -0.5 * mDt * p_0[5];
      w[2] = -0.5 * mDt * p_0[6];
      JinvWJ = mat<double,mat_structure::identity>(3) + 0.5 * JinvWJ;
      vect<double,3> w_JinvWJ = w * JinvWJ;
      set_block(A, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >(w_JinvWJ),3,10);
      set_block(A, ((0.5 * mDt * p_0[3]) * mat<double,mat_structure::identity>(3) - mat<double,mat_structure::skew_symmetric>(w)) * JinvWJ,4,10);
      
      B = mat<double,mat_structure::nil>(13,6);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mMass;

      w *= mDt * 0.5;
      vect<double,3> w_Jinv = w * mInertiaMomentInv;
      set_block(B, mat_vect_adaptor< vect<double,3>, mat_alignment::row_major >( w_Jinv ), 3, 3);
      set_block(B, (0.25 * mDt * mDt * p_0[3]) * mInertiaMomentInv - mat<double,mat_structure::skew_symmetric>(w) * mInertiaMomentInv, 4, 3);
      
      B(7,0) = mDt / mMass;
      B(8,1) = mDt / mMass;
      B(9,2) = mDt / mMass;
      set_block(B, mDt * mInertiaMomentInv, 10, 3);
      
    };
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const pp::vector_topology< vect_n<double> >&, 
                                    const time_type&, const point_type&, const input_type&) const {
      C = mat<double,mat_structure::nil>(7,13);
      set_block(C,mat<double,mat_structure::identity>(7),0,0);
      
      D = mat<double,mat_structure::nil>(7,6);
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

    RK_RTTI_MAKE_CONCRETE_1BASE(airship3D_lin2_dt_system,0xC231000D,1,"airship3D_lin2_dt_system",airship3D_lin_dt_system)
    
};

#endif


};

};

#endif









