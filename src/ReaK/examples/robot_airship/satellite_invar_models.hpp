
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

#include "satellite_basic_models.hpp"

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
			 double aDt = 0.001);  
  
    virtual ~satellite2D_imdt_sys() { };
        
    time_difference_type get_time_step() const { return mDt; };
    
    point_type get_next_state(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const state_space_type&, const time_type&, const point_type&, const input_type&) const;
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, const time_type&, const point_type&, const input_type&) const;
        
    invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const;
    
    point_type apply_correction(const state_space_type&, const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const;
    
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

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite2D_imdt_sys,0xC2310012,1,"satellite2D_imdt_sys",named_object)
    
};





class satellite3D_imdt_sys : public satellite3D_inv_dt_system {
  public:
    
    typedef satellite3D_inv_dt_system::state_space_type state_space_type;
    
    typedef satellite3D_inv_dt_system::point_type point_type;
    typedef satellite3D_inv_dt_system::point_difference_type point_difference_type;
  
    typedef satellite3D_inv_dt_system::time_type time_type;
    typedef satellite3D_inv_dt_system::time_difference_type time_difference_type;
  
    typedef satellite3D_inv_dt_system::input_type input_type;
    typedef satellite3D_inv_dt_system::output_type output_type;
  
    typedef satellite3D_inv_dt_system::invariant_error_type invariant_error_type;
    typedef satellite3D_inv_dt_system::invariant_correction_type invariant_correction_type;
    typedef satellite3D_inv_dt_system::invariant_frame_type invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef satellite3D_inv_dt_system::matrixA_type matrixA_type;
    typedef satellite3D_inv_dt_system::matrixB_type matrixB_type;
    typedef satellite3D_inv_dt_system::matrixC_type matrixC_type;
    typedef satellite3D_inv_dt_system::matrixD_type matrixD_type;
    
    typedef satellite3D_inv_dt_system::zero_input_trajectory zero_input_trajectory;
    
    
  public:  
    satellite3D_imdt_sys(const std::string& aName = "", 
			 double aMass = 1.0, 
			 const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			 double aDt = 0.001); 
  
    virtual ~satellite3D_imdt_sys() { };
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const state_space_type&, 
				     const time_type&, const time_type&,
				     const point_type& p_0, const point_type& p_1,
				     const input_type&, const input_type&) const;
    
    point_type apply_correction(const state_space_type&, const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const;
    
    invariant_frame_type get_invariant_prior_frame(const state_space_type&, const point_type& x_prev, const point_type& x_prior, const input_type&, const time_type&) const;
    
    invariant_frame_type get_invariant_posterior_frame(const state_space_type& state_space, const point_type& x_prior, const point_type& x_post, const input_type& u, const time_type& t) const {
      return get_invariant_prior_frame(state_space,x_prior,x_post,u,t);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_inv_dt_system::save(A,satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_inv_dt_system::load(A,satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_imdt_sys,0xC2310015,1,"satellite3D_imdt_sys",satellite3D_inv_dt_system)
    
};







class satellite3D_gyro_imdt_sys : public satellite3D_imdt_sys {
  public:
    
    typedef satellite3D_imdt_sys::state_space_type state_space_type;
    
    typedef satellite3D_imdt_sys::point_type point_type;
    typedef satellite3D_imdt_sys::point_difference_type point_difference_type;
  
    typedef satellite3D_imdt_sys::time_type time_type;
    typedef satellite3D_imdt_sys::time_difference_type time_difference_type;
  
    typedef satellite3D_imdt_sys::input_type input_type;
    typedef satellite3D_imdt_sys::output_type output_type;
  
    typedef satellite3D_imdt_sys::invariant_error_type invariant_error_type;
    typedef satellite3D_imdt_sys::invariant_correction_type invariant_correction_type;
    typedef satellite3D_imdt_sys::invariant_frame_type invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 14);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 12);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef satellite3D_imdt_sys::matrixA_type matrixA_type;
    typedef satellite3D_imdt_sys::matrixB_type matrixB_type;
    typedef satellite3D_imdt_sys::matrixC_type matrixC_type;
    typedef satellite3D_imdt_sys::matrixD_type matrixD_type;
    
    typedef satellite3D_imdt_sys::zero_input_trajectory zero_input_trajectory;
    
  public:  
    satellite3D_gyro_imdt_sys(const std::string& aName = "", 
			      double aMass = 1.0, 
			      const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			      double aDt = 0.001); 
  
    virtual ~satellite3D_gyro_imdt_sys() { };
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
        
    invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_imdt_sys::save(A,satellite3D_imdt_sys::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_imdt_sys::load(A,satellite3D_imdt_sys::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_gyro_imdt_sys,0xC2310016,1,"satellite3D_gyro_imdt_sys",satellite3D_imdt_sys)
    
};








class satellite3D_IMU_imdt_sys : public satellite3D_imdt_sys {
  public:
    
    typedef satellite3D_imdt_sys::state_space_type state_space_type;
    
    typedef satellite3D_imdt_sys::point_type point_type;
    typedef satellite3D_imdt_sys::point_difference_type point_difference_type;
  
    typedef satellite3D_imdt_sys::time_type time_type;
    typedef satellite3D_imdt_sys::time_difference_type time_difference_type;
  
    typedef satellite3D_imdt_sys::input_type input_type;
    typedef satellite3D_imdt_sys::output_type output_type;
  
    typedef satellite3D_imdt_sys::invariant_error_type invariant_error_type;
    typedef satellite3D_imdt_sys::invariant_correction_type invariant_correction_type;
    typedef satellite3D_imdt_sys::invariant_frame_type invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 20);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 12);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef satellite3D_imdt_sys::matrixA_type matrixA_type;
    typedef satellite3D_imdt_sys::matrixB_type matrixB_type;
    typedef satellite3D_imdt_sys::matrixC_type matrixC_type;
    typedef satellite3D_imdt_sys::matrixD_type matrixD_type;
    
    typedef satellite3D_imdt_sys::zero_input_trajectory zero_input_trajectory;
    
  protected:
    
    unit_quat<double> IMU_orientation;
    vect<double,3> IMU_location;
    unit_quat<double> room_orientation;
    vect<double,3> mag_field_vector;
    
  public:  
    satellite3D_IMU_imdt_sys(const std::string& aName = "", 
			     double aMass = 1.0, 
			     const mat<double,mat_structure::symmetric>& aInertiaMoment = mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
			     double aDt = 0.001,
			     const unit_quat<double>& aIMUOrientation = unit_quat<double>(),
			     const vect<double,3>& aIMULocation = vect<double,3>(),
			     const unit_quat<double>& aRoomOrientation = unit_quat<double>(),
			     const vect<double,3>& aMagFieldVector = vect<double,3>(1.0,0.0,0.0)); 
  
    virtual ~satellite3D_IMU_imdt_sys() { };
    
    output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
        
    invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_imdt_sys::save(A,satellite3D_imdt_sys::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
        & RK_SERIAL_SAVE_WITH_NAME(IMU_location)
	& RK_SERIAL_SAVE_WITH_NAME(room_orientation)
	& RK_SERIAL_SAVE_WITH_NAME(mag_field_vector);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_imdt_sys::load(A,satellite3D_imdt_sys::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
        & RK_SERIAL_LOAD_WITH_NAME(IMU_location)
	& RK_SERIAL_LOAD_WITH_NAME(room_orientation)
	& RK_SERIAL_LOAD_WITH_NAME(mag_field_vector);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_IMU_imdt_sys,0xC2310017,1,"satellite3D_IMU_imdt_sys",satellite3D_imdt_sys)
    
};







};

};

#endif

















