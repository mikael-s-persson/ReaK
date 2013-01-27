
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

#ifndef SATELLITE_BASIC_MODELS_HPP
#define SATELLITE_BASIC_MODELS_HPP

#include "base/named_object.hpp"

#include "ctrl_sys/sss_exceptions.hpp"
#include "topologies/se3_topologies.hpp"

#include "lin_alg/mat_alg.hpp"
#include "lin_alg/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


  




class satellite3D_lin_dt_system : public named_object {
  public:
  
    typedef pp::se3_1st_order_topology<double>::type state_space_type;
    
    typedef pp::topology_traits< state_space_type >::point_type point_type;
    typedef pp::topology_traits< state_space_type >::point_difference_type point_difference_type;
    
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
    typedef mat<double,mat_structure::rectangular> matrixD_type;
    
    struct zero_input_trajectory {
      input_type get_point(time_type) const {
	return input_type(0.0,0.0,0.0,0.0,0.0,0.0);
      };
    };
    
  protected:
    double mMass;
    mat<double,mat_structure::symmetric> mInertiaMoment;
    mat<double,mat_structure::symmetric> mInertiaMomentInv;
    time_difference_type mDt;
    
  public:
    
    satellite3D_lin_dt_system(const std::string& aName = "", 
			      double aMass = 1.0, 
			      const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
			      double aDt = 0.001); 
  
    virtual ~satellite3D_lin_dt_system() { };
    
    time_difference_type get_time_step() const { return mDt; };
    
    virtual point_type get_next_state(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    virtual void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, 
				     const state_space_type&, 
				     const time_type& t_0, const time_type&,
				     const point_type& p_0, const point_type&,
				     const input_type&, const input_type&) const;
    
    virtual output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    virtual void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
    
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
	throw system_incoherency("Inertial information is improper in satellite3D_lin_dt_system's definition");
      try {
        invert_Cholesky(mInertiaMoment,mInertiaMomentInv);
      } catch(singularity_error&) {
	throw system_incoherency("Inertial tensor is singular in satellite3D_lin_dt_system's definition");
      };
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_lin_dt_system,0xC2310013,1,"satellite3D_lin_dt_system",named_object)
    
};





class satellite3D_gyro_lin_dt_system : public satellite3D_lin_dt_system {
  public:
  
    typedef satellite3D_lin_dt_system::state_space_type state_space_type;
    
    typedef satellite3D_lin_dt_system::point_type point_type;
    typedef satellite3D_lin_dt_system::point_difference_type point_difference_type;
    
    typedef satellite3D_lin_dt_system::time_type time_type;
    typedef satellite3D_lin_dt_system::time_difference_type time_difference_type;
  
    typedef satellite3D_lin_dt_system::input_type input_type;
    typedef satellite3D_lin_dt_system::output_type output_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 14);
    
    typedef satellite3D_lin_dt_system::matrixA_type matrixA_type;
    typedef satellite3D_lin_dt_system::matrixB_type matrixB_type;
    typedef satellite3D_lin_dt_system::matrixC_type matrixC_type;
    typedef satellite3D_lin_dt_system::matrixD_type matrixD_type;

    typedef satellite3D_lin_dt_system::zero_input_trajectory zero_input_trajectory;
    
  public:
    
    satellite3D_gyro_lin_dt_system(const std::string& aName = "", 
			           double aMass = 1.0, 
			           const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
			           double aDt = 0.001); 
  
    virtual ~satellite3D_gyro_lin_dt_system() { };
    
    virtual output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    virtual void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_lin_dt_system::save(A,satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_lin_dt_system::load(A,satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_gyro_lin_dt_system,0xC2310018,1,"satellite3D_gyro_lin_dt_system",satellite3D_lin_dt_system)
    
};











class satellite3D_inv_dt_system : public satellite3D_lin_dt_system {
  public:
    
    typedef satellite3D_lin_dt_system::state_space_type state_space_type;
    
    typedef satellite3D_lin_dt_system::point_type point_type;
    typedef satellite3D_lin_dt_system::point_difference_type point_difference_type;
  
    typedef satellite3D_lin_dt_system::time_type time_type;
    typedef satellite3D_lin_dt_system::time_difference_type time_difference_type;
  
    typedef satellite3D_lin_dt_system::input_type input_type;
    typedef satellite3D_lin_dt_system::output_type output_type;
  
    typedef vect_n<double> invariant_error_type;
    typedef vect_n<double> invariant_correction_type;
    typedef mat<double,mat_structure::square> invariant_frame_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 7);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 6);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef satellite3D_lin_dt_system::matrixA_type matrixA_type;
    typedef satellite3D_lin_dt_system::matrixB_type matrixB_type;
    typedef satellite3D_lin_dt_system::matrixC_type matrixC_type;
    typedef satellite3D_lin_dt_system::matrixD_type matrixD_type;
    
    typedef satellite3D_lin_dt_system::zero_input_trajectory zero_input_trajectory;
    
    satellite3D_inv_dt_system(const std::string& aName = "", 
			      double aMass = 1.0, 
			      const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
			      double aDt = 0.001); 
  
    virtual ~satellite3D_inv_dt_system() { };
    
    void get_state_transition_blocks(matrixA_type& A, matrixB_type& B, const state_space_type&, 
				     const time_type& t_0, const time_type&,
				     const point_type& p_0, const point_type&,
				     const input_type&, const input_type&) const;
    
    void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
    
    virtual invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const;
    
    virtual point_type apply_correction(const state_space_type&, const point_type& x, const invariant_correction_type& c, const input_type&, const time_type&) const;
    
    virtual invariant_frame_type get_invariant_prior_frame(const state_space_type&, const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(mat<double,mat_structure::identity>(invariant_correction_dimensions));
    };
    
    virtual invariant_frame_type get_invariant_posterior_frame(const state_space_type&, const point_type&, const point_type&, const input_type&, const time_type&) const {
      return invariant_frame_type(mat<double,mat_structure::identity>(invariant_correction_dimensions));
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_lin_dt_system::save(A,satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_lin_dt_system::load(A,satellite3D_lin_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_inv_dt_system,0xC2310014,1,"satellite3D_inv_dt_system",satellite3D_lin_dt_system)
    
};




class satellite3D_gyro_inv_dt_system : public satellite3D_inv_dt_system {
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
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 14);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_error_dimensions = 12);
    BOOST_STATIC_CONSTANT(std::size_t, invariant_correction_dimensions = 12);
    
    typedef satellite3D_inv_dt_system::matrixA_type matrixA_type;
    typedef satellite3D_inv_dt_system::matrixB_type matrixB_type;
    typedef satellite3D_inv_dt_system::matrixC_type matrixC_type;
    typedef satellite3D_inv_dt_system::matrixD_type matrixD_type;

    typedef satellite3D_inv_dt_system::zero_input_trajectory zero_input_trajectory;
    
  public:
    
    satellite3D_gyro_inv_dt_system(const std::string& aName = "", 
			           double aMass = 1.0, 
			           const mat<double,mat_structure::symmetric>& aInertiaMoment = (mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3))),
			           double aDt = 0.001); 
  
    virtual ~satellite3D_gyro_inv_dt_system() { };
    
    virtual output_type get_output(const state_space_type&, const point_type& x, const input_type& u, const time_type& t = 0.0) const;
    
    virtual void get_output_function_blocks(matrixC_type& C, matrixD_type& D, const state_space_type&, 
				    const time_type&, const point_type&, const input_type&) const;
    
    virtual invariant_error_type get_invariant_error(const state_space_type&, const point_type& x, const input_type& u, const output_type& y, const time_type& t) const;
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      satellite3D_inv_dt_system::save(A,satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      satellite3D_inv_dt_system::load(A,satellite3D_inv_dt_system::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(satellite3D_gyro_inv_dt_system,0xC2310019,1,"satellite3D_gyro_inv_dt_system",satellite3D_inv_dt_system)
    
};




};

};

#endif











