/**
 * \file discretized_lti_sys.hpp
 * 
 * This library provides a class template which can take a continuous-time linear time-invariant 
 * state-space system and turn it into a discrete-time linear time-invariant state-space system.
 * The class template uses a matrix exponential method to compute the system matrix of the 
 * discrete-time system, and thus, the resulting system is exact, as far as the matrix exponential
 * is correct (see mat_exp_methods.hpp).
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
 */

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

#ifndef DISCRETIZED_LTI_SYS_HPP
#define DISCRETIZED_LTI_SYS_HPP

#include "linear_ss_system_concept.hpp"
#include "discrete_linear_sss_concept.hpp"

#include "math/mat_concepts.hpp"
#include "math/mat_alg.hpp"
#include "math/mat_exp_methods.hpp"
#include "math/mat_qr_decomp.hpp"

#include "base/named_object.hpp"

namespace ReaK {

namespace ctrl {


/**
 * This class template can take a continuous-time linear time-invariant 
 * state-space system and turn it into a discrete-time linear time-invariant state-space system.
 * The class template uses a matrix exponential method to compute the system matrix of the 
 * discrete-time system, and thus, the resulting system is exact, as far as the matrix exponential
 * is correct (see mat_exp_methods.hpp).
 *
 * \tparam LTISystem The linear time-invariance, continuous-time state-space system type, should model LinearSSSystemConcept with LTISystemType.
 */
template <typename LTISystem>
class discretized_lti_sys : public named_object {
  public:
    typedef discretized_lti_sys<LTISystem> self;
    typedef LTISystem::value_type value_type;
    typedef LTISystem::size_type size_type;
    
    typedef typename ss_system_traits<LTISystem>::point_type point_type;
    typedef typename ss_system_traits<LTISystem>::point_difference_type point_difference_type;
  
    typedef typename ss_system_traits<LTISystem>::time_type time_type;
    typedef typename ss_system_traits<LTISystem>::time_difference_type time_difference_type;
  
    typedef typename ss_system_traits<LTISystem>::input_type input_type;
    typedef typename ss_system_traits<LTISystem>::output_type output_type;
    
    typedef typename linear_ss_system_traits<LTISystem>::matrixA_type matrixA_type;
    typedef typename linear_ss_system_traits<LTISystem>::matrixB_type matrixB_type;
    typedef typename linear_ss_system_traits<LTISystem>::matrixC_type matrixC_type;
    typedef typename linear_ss_system_traits<LTISystem>::matrixD_type matrixD_type;
  
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = ss_system_traits<LTISystem>::dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = ss_system_traits<LTISystem>::input_dimensions);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = ss_system_traits<LTISystem>::output_dimensions);
        
  private:
    LTISystem sys;
    time_difference_type dt;
    
    matrixA_type Ad;
    matrixB_type Bd;
    matrixC_type Cd;
    matrixD_type Dd;
    
  public:
    
    /**
     * Parametrized and default constructor.
     * \param aX_size The size of the state-vector.
     * \param aU_size The size of the input-vector.
     * \param aY_size The size of the output-vector.
     * \param aDt The time-step of the discrete-time system.
     */
    discretized_lti_sys(size_type aX_size = 0, size_type aU_size = 0, size_type aY_size = 0, const time_difference_type& aDt = 1, const std::string& aName = "") :
                        sys(aX_size,aU_size,aY_size,aName + "_continuous"), dt(aDt) {
      setName(aName);
      
      sys.get_linear_blocks(Ad,Bd,Cd,Dd);
      
      Ad = mat<value_type, mat_structure::identity>(aX_size);
    };

    /**
     * Standard copy-constructor.
     */
    discretized_lti_sys(const self& rhs) : sys(rhs.sys), dt(rhs.dt), Ad(rhs.Ad), Bd(rhs.Bd), Cd(rhs.Cd), Dd(rhs.Dd) {
      setName(rhs.getName());
    };
    
    /**
     * Parametrized constructor.
     * \param aA The continuous-time system matrix A.
     * \param aB The continuous-time system matrix B.
     * \param aC The continuous-time system matrix C.
     * \param aD The continuous-time system matrix D.
     * \param aDt The time-step of the discrete-time system.
     */
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    discretized_lti_sys(const MatrixA& aA, const MatrixB& aB, const MatrixC& aC, const MatrixD& aD, const time_difference_type& aDt, const std::string& aName = "",
                        typename boost::enable_if_c< is_readable_matrix<MatrixA>::value &&
                                                     is_readable_matrix<MatrixB>::value &&
                                                     is_readable_matrix<MatrixC>::value &&
                                                     is_readable_matrix<MatrixD>::value , void* >::type dummy = NULL) : 
                        sys(aA,aB,aC,aD,aName + "_continuous"), dt(aDt) { 
      setName(aName);
      
      sys.get_linear_blocks(Ad,Bd,Cd,Dd);
      
      mat<value_type, mat_structure::square> A_aug(Ad.get_col_count() + Bd.get_col_count(), value_type(0));
      set_block(A_aug, Ad * dt, 0, 0);
      set_block(A_aug, Bd * dt, 0, Ad.get_col_count());
      mat<value_type, mat_structure::square> A_aug_exp(Ad.get_col_count() + Bd.get_col_count(), value_type(0));
      exp_PadeSAS(A_aug,A_aug_exp,QR_linlsqsolver());
      Ad = get_block(A_aug_exp, 0, 0, Ad.get_row_count(), Ad.get_col_count());
      Bd = get_block(A_aug_exp, 0, Ad.get_col_count(), Bd.get_row_count(), Bd.get_col_count());
    };
  
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.sys,rhs.sys);
      swap(lhs.dt,rhs.dt);
      swap(lhs.Ad,rhs.Ad);
      swap(lhs.Bd,rhs.Bd);
      swap(lhs.Cd,rhs.Cd);
      swap(lhs.Dd,rhs.Dd);
    };
 
    /**
     * Standard assignment operator.
     */
    self& operator =(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    /**
     * Returns the underlying continuous-time system.
     */
    const LTISystem& getSys() const { return sys; };

    /**
     * Fills the given matrices with the discrete-time system matrices.
     * \param aA Stores, as output, the system matrix A.
     * \param aB Stores, as output, the system matrix B.
     * \param aC Stores, as output, the system matrix C.
     * \param aD Stores, as output, the system matrix D.
     */
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD) const {
      aA = Ad;
      aB = Bd;
      aC = Cd;
      aD = Dd;
    };
    
    /**
     * Fills the given matrices with the discrete-time system matrices.
     * \param aA Stores, as output, the system matrix A.
     * \param aB Stores, as output, the system matrix B.
     * \param aC Stores, as output, the system matrix C.
     * \param aD Stores, as output, the system matrix D.
     */
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD,const time_type&) const {
      aA = Ad;
      aB = Bd;
      aC = Cd;
      aD = Dd;
    };
    
    /**
     * Fills the given matrices with the discrete-time system matrices.
     * \param aA Stores, as output, the system matrix A.
     * \param aB Stores, as output, the system matrix B.
     * \param aC Stores, as output, the system matrix C.
     * \param aD Stores, as output, the system matrix D.
     */
    template <typename MatrixA, typename MatrixB, typename MatrixC, typename MatrixD>
    typename boost::enable_if_c< is_writable_matrix<MatrixA>::value &&
                                 is_writable_matrix<MatrixB>::value &&
                                 is_writable_matrix<MatrixC>::value &&
                                 is_writable_matrix<MatrixD>::value,
    void >::type get_linear_blocks(MatrixA& aA, MatrixB& aB, MatrixC& aC, MatrixD& aD, const time_type&, const point_type&, const input_type&) const {
      aA = Ad;
      aB = Bd;
      aC = Cd;
      aD = Dd;
    };
    
    /**
     * Returns the time-step of the discrete-time system.
     */
    time_difference_type get_time_step() const { return dt; };
    
    /**
     * Returns next state of the system given the current state, input and time.
     * \param p The current state.
     * \param u The current input.
     * \param t The current time.
     * \return The next state, at t + get_time_step().
     */
    point_type get_next_state(const point_type& p, const input_type& u, const time_type& t = 0) const { RK_UNUSED(t);
      return Ad * p + Bd * u;
    };
    
    /**
     * Returns output of the system given the current state, input and time.
     * \param p The current state.
     * \param u The current input.
     * \param t The current time.
     * \return The current output.
     */
    output_type get_output(const point_type& p, const input_type& u, const time_type& t = 0) const { RK_UNUSED(t);
      return Cd * p + Dd * u;
    };
    
    /**
     * Adjusts the state by adding a state difference to it.
     */
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return sys.adjust(p,dp);
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(sys)
         & RK_SERIAL_SAVE_WITH_NAME(dt);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(sys)
         & RK_SERIAL_LOAD_WITH_NAME(dt);
      
      sys.get_linear_blocks(Ad,Bd,Cd,Dd);
      
      mat<value_type, mat_structure::square> A_aug(Ad.get_col_count() + Bd.get_col_count(), value_type(0));
      set_block(A_aug, Ad * dt, 0, 0);
      set_block(A_aug, Bd * dt, 0, Ad.get_col_count());
      mat<value_type, mat_structure::square> A_aug_exp(Ad.get_col_count() + Bd.get_col_count(), value_type(0));
      exp_PadeSAS(A_aug,A_aug_exp,QR_linlsqsolver());
      Ad = get_block(A_aug_exp, 0, 0, Ad.get_row_count(), Ad.get_col_count());
      Bd = get_block(A_aug_exp, 0, Ad.get_col_count(), Bd.get_row_count(), Bd.get_col_count());
	 
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300005,1,"discretized_lti_sys",named_object)
  
};
  





};

};

#endif








