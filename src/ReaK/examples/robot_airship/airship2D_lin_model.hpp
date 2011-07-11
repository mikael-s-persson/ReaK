
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
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
    typedef vect_n<double> point_derivative_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
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
			 named_object(),
			 mMass(aMass),
			 mInertiaMoment(aInertiaMoment) {
      setName(aName);
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
      A(2,6) = -x[3];
      A(3,6) = x[2];
      
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
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 6);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
   
    airship2D_inv_system(const std::string& aName = "", 
			 double aMass = 1.0, 
			 double aInertiaMoment = 1.0) :
			 airship2D_lin_system(aName,aMass,aInertiaMoment) { }; 
  
    virtual ~airship2D_inv_system() { };
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::nil>(6,6);
      A(0,3) = 1.0;
      A(1,4) = 1.0;
      A(2,5) = 1.0;
      
      B = mat<double,mat_structure::nil>(6,3);
      B(3,0) = 1.0 / mMass;
      B(4,1) = 1.0 / mMass;
      B(5,2) = 1.0 / mInertiaMoment;
      
      C = mat<double,mat_structure::nil>(3,6);
      set_block(C,mat<double,mat_structure::identity>(3),0,0);
      
      D = mat<double,mat_structure::nil>(3,3);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      vect<double,2> qy(y[2],y[3]); qy = unit(qy);
      vect<double,2> qx(x[2],x[3]); qx = unit(qx);
      return invariant_error_type(y[0] - x[0],
			          y[1] - x[1],
			          //qy[0] * qx[0] + qy[1] * qx[1], // c_y * c_x + s_y * s_x
                                  qy[1] * qx[0] - qy[0] * qx[1]); // s_y * c_x - c_y * s_x
    };
    
    point_derivative_type apply_correction(const point_type& x, const point_derivative_type& xd, const invariant_correction_type& c, const input_type& u, const time_type& t) const {
      return point_type(xd[0] + c[0],
	                xd[1] + c[1],
			xd[2] - x[3] * c[2],
			xd[3] + x[2] * c[2],
			xd[4] + c[3],
			xd[5] + c[4],
			xd[6] + c[5]);
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
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
  protected:
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
      double delta_theta = mDt * (x[6] + 0.5 * dv[2]); double cdt = cos(delta_theta); double sdt = sin(delta_theta);
      return point_type( x[0] + mDt * (x[4] + 0.5 * dv[0]),
			 x[1] + mDt * (x[5] + 0.5 * dv[1]),
			 x[2] * cdt - x[3] * sdt,
			 x[2] * sdt + x[3] * cdt,
			 x[4] + dv[0],
			 x[5] + dv[1],
			 x[6] + dv[2]);
    };
    
    output_type get_output(const point_type& x, const input_type& u, const time_type t = 0.0) const {
      return output_type(x[0], x[1], x[2], x[3]);
    };
    
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::identity>(7);
      A(0,4) = mDt;
      A(1,5) = mDt;  
      A(2,6) = -mDt * x[3];
      A(3,6) = mDt * x[2];
      
      B = mat<double,mat_structure::nil>(7,3);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = -0.5 * x[3] * mDt * mDt / mInertiaMoment;
      B(3,2) = 0.5 * x[2] * mDt * mDt / mInertiaMoment;
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





class airship2D_inv_dt_system : public airship2D_lin_dt_system {
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
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 4);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 3);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 6);
    
    typedef mat<double,mat_structure::square> matrixA_type;
    typedef mat<double,mat_structure::rectangular> matrixB_type;
    typedef mat<double,mat_structure::rectangular> matrixC_type;
    typedef mat<double,mat_structure::nil> matrixD_type;
    
    
    airship2D_inv_dt_system(const std::string& aName = "", 
			    double aMass = 1.0, 
			    double aInertiaMoment = 1.0,
			    double aDt = 0.001) :
			    airship2D_lin_dt_system(aName,aMass,aInertiaMoment,aDt) { }; 
  
    virtual ~airship2D_inv_dt_system() { };
        
    void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C, matrixD_type& D, const time_type& t, const point_type& x, const input_type& u) const {
      A = mat<double,mat_structure::identity>(6);
      A(0,3) = mDt;
      A(1,4) = mDt;  
      A(2,5) = mDt;
      
      B = mat<double,mat_structure::nil>(6,3);
      B(0,0) = 0.5 * mDt * mDt / mMass;
      B(1,1) = 0.5 * mDt * mDt / mMass;
      B(2,2) = 0.5 * mDt * mDt / mInertiaMoment;
      B(3,0) = mDt / mMass;
      B(4,1) = mDt / mMass;
      B(5,2) = mDt / mInertiaMoment;
      
      C = mat<double,mat_structure::nil>(3,6);
      set_block(C,mat<double,mat_structure::identity>(3),0,0);
      
      D = mat<double,mat_structure::nil>(3,3);
    };
        
    invariant_error_type get_invariant_error(const point_type& x, const input_type& u, const output_type& y, const time_type& t) const {
      return invariant_error_type(y[0] - x[0],
			          y[1] - x[1],
			          y[3] * x[2] - y[2] * x[3]); // s_y * c_x - c_y * s_x
    };
    
    point_type apply_correction(const point_type& x, const invariant_correction_type& c, const input_type& u, const time_type& t) const {
      using std::fabs;
      using std::sqrt;
      double sc = c[2];
      if(fabs(sc) > 1) 
	sc /= (fabs(sc) + std::numeric_limits< double >::epsilon());
      double cc = sqrt(1 - sc * sc);
      return point_type(x[0] + c[0],
	                x[1] + c[1],
			x[2] * cc - x[3] * sc,
			x[3] * cc + x[2] * sc,
			x[4] + c[3],
			x[5] + c[4],
			x[6] + c[5]);
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
      airship2D_lin_dt_system::save(A,airship2D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      airship2D_lin_dt_system::load(A,airship2D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(airship2D_inv_dt_system,0xC2310004,1,"airship2D_inv_dt_system",airship2D_lin_dt_system)
    
};








};

};

#endif









