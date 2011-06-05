
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

#ifndef AIRSHIP2D_LIN_MODEL_HPP
#define AIRSHIP2D_LIN_MODEL_HPP

#include "math/vect_alg.hpp"
#include "base/named_object.hpp"

#include "math/mat_alg.hpp"
#include <ctrl_sys/sss_exceptions.hpp>


namespace ReaK {


namespace ctrl {
  


class airship2D_lin_system : public named_object {
  public:
    
    typedef vect<double,7> point_type;
    typedef vect<double,7> point_difference_type;
    typedef vect<double,7> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect<double,3> input_type;
    typedef vect<double,4> output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  protected:
    double mMass;
    double mInertiaMoment;
    
  public:
    
    airship2D_lin_system(const std::string& aName = "", 
			 double aMass = 1.0, 
			 double aInertiaMoment = 1.0) :
			 named_object(aName),
			 mMass(aMass),
			 mInertiaMoment(aInertiaMoment) { 
      if((mInertiaMoment < std::numeric_limits< double >::epsilon()) || (mMass < std::numeric_limits< double >::epsilon()))
	throw system_incoherency("Inertial information are singular in airship2D_lin_system's definition");
    }; 
  
    virtual ~airship2D_lin_system() { };
    
    point_derivative_type get_state_derivative(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return point_derivative_type(x[4],
	                           x[5],
				   -x[4] * x[6],
				   x[3] * x[6],
				   u[0] / mMass,
				   u[1] / mMass,
				   u[2] / mInertiaMoment);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::nil>(7,7);
      A(0,4) = 1.0;
      A(1,5) = 1.0;
      A(2,6) = -x[4];
      A(3,6) = x[3];
      
      B = mat<double,mat_structure::nil>(7,3);
      B(4,0) = 1.0 / mMass;
      B(5,1) = 1.0 / mMass;
      B(6,2) = 1.0 / mInertiaMoment;
      
      C = mat<double,mat_structure::nil>(4,7);
      set_block(C,mat<double,mat_structure::identity>(4),0,0);
      
      D = mat<double,mat_structure::nil>(4,3);
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
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship2D_lin_system,0xC2310001,1,"airship2D_lin_system",named_object)
    
};




class airship2D_inv_system : public airship2D_lin_system {
  public:
    
    typedef vect<double,7> point_type;
    typedef vect<double,7> point_difference_type;
    typedef vect<double,7> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect<double,3> input_type;
    typedef vect<double,4> output_type;
    
    typedef vect<double,3> invariant_type;
    typedef vect<double,4> output_error_type;
    typedef mat<double,mat_structure::square> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_error_dimensions = 4);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
   
    airship2D_inv_system(const std::string& aName = "", 
			 double aMass = 1.0, 
			 double aInertiaMoment = 1.0) :
			 airship2D_lin_system(aName,aMass,aInertiaMoment) { }; 
  
    virtual ~airship2D_inv_system() { };
    
    point_derivative_type get_state_derivative(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return point_derivative_type(x[4],
	                           x[5],
				   -x[4] * x[6],
				   x[3] * x[6],
				   u[0] / mMass,
				   u[1] / mMass,
				   u[2] / mInertiaMoment);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::nil>(7,7);
      A(0,4) = 1.0;
      A(1,5) = 1.0;
      A(3,6) = 1.0;
      
      B = mat<double,mat_structure::nil>(7,3);
      B(4,0) = 1.0 / mMass;
      B(5,1) = 1.0 / mMass;
      B(6,2) = 1.0 / mInertiaMoment;
      
      C = mat<double,mat_structure::nil>(4,7);
      set_block(C,mat<double,mat_structure::identity>(4),0,0);
      
      D = mat<double,mat_structure::nil>(4,3);
    };
    
    void get_invariant(invariant_type& i, const time_type& t, const point_type& x, const input_type& u) const {
      i = u;
    };
    
    void get_output_error(output_error_type& e, const time_type& t, const point_type& x, const input_type& u, const output_type& y) const {
      e[0] = y[0] - x[0];
      e[1] = y[1] - x[1];
      e[2] = y[2] * x[2] + y[3] * x[3]; // c_y * c_x + s_y * s_x
      e[3] = y[3] * x[2] - y[2] * x[3]; // s_y * c_x - c_y * s_x
    };
    
    void get_invariant_frame(invariant_frame_type& W, const time_type& t, const point_type& x) const {
      W = mat<double,mat_structure::identity>(7);
      W(2,2) = x[2]; W(2,3) = -x[3];  // R(q)
      W(3,2) = x[3]; W(3,3) = x[2];   //
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      airship2D_lin_system::save(A,airship2D_lin_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship2D_lin_system::load(A,airship2D_lin_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship2D_inv_system,0xC2310002,1,"airship2D_inv_system",airship2D_lin_system)
    
};





class airship2D_lin_dt_system : public airship2D_lin_system {
  public:
    
    typedef vect<double,7> point_type;
    typedef vect<double,7> point_difference_type;
    typedef vect<double,7> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect<double,3> input_type;
    typedef vect<double,4> output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  private:
    time_difference_type mDt;
    
  public:
    
    airship2D_lin_dt_system(const std::string& aName = "", 
			    double aMass = 1.0, 
			    double aInertiaMoment = 1.0,
			    double aDt = 0.001) :
			    airship2D_lin_system(aName,aMass,aInertiaMoment),
			    mDt(aDt) { 
      if(mDt < std::numeric_limits< double >::epsilon())
	throw system_incoherency("The time step is below numerical tolerance in airship2D_lin_dt_system's definition");
    }; 
  
    virtual ~airship2D_lin_dt_system() { };
    
    time_difference_type get_time_step() const { return mDt; };
    
    point_type get_next_state(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      input_type dv(mDt * u[0] / mMass, mDt * u[1] / mMass, mDt * u[2] / mInertiaMoment);
      double delta_theta = mDt * (x[6] + 0.5 * dv[2]); double cdt = cos(delta_theta) - 1.0; double sdt = sin(delta_theta);
      point_type r = x + point_difference_type( mDt * (x[4] + 0.5 * dv[0]),
						mDt * (x[5] + 0.5 * dv[1]),
					        x[3] * cdt - x[4] * sdt,
						x[3] * sdt + x[4] * cdt,
						dv[0],
						dv[1],
						dv[2]);
      using std::sqrt;
      double rcs_mag = sqrt( r[2] * r[2] + r[3] * r[3] );
      r[2] /= rcs_mag;
      r[3] /= rcs_mag;
      return r;
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::identity>(7);
      A(0,4) = mDt;
      A(1,5) = mDt;  
      A(2,6) = -mDt * x[4];
      A(3,6) = mDt * x[3];
      
      B = mat<double,mat_structure::nil>(7,3);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = -0.5 * x[4] * mDt * mDt / mInertiaMoment;
      B(3,2) = 0.5 * x[3] * mDt * mDt / mInertiaMoment;
      B(4,0) = mDt / mMass;
      B(5,1) = mDt / mMass;
      B(6,2) = mDt / mInertiaMoment;
      
      C = mat<double,mat_structure::nil>(4,7);
      set_block(C,mat<double,mat_structure::identity>(4),0,0);
      
      D = mat<double,mat_structure::nil>(4,3);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      airship2D_lin_system::save(A,airship2D_lin_system::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mDt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship2D_lin_system::load(A,airship2D_lin_system::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mDt);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship2D_lin_dt_system,0xC2310003,1,"airship2D_lin_dt_system",airship2D_lin_system)
    
};





};

};

#endif









