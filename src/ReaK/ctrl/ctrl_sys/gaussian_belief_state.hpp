/**
 * \file gaussian_belief_state.hpp
 * 
 * This library provides a number of class templates that can be used to represent and 
 * use a Gaussian belief-state (i.e. a mean-state and a covariance matrix). Those class
 * templates include the belief-state itself, of course, a sampler to generate random
 * states from the belief-state's probability distribution, and a PDF (probability 
 * distribution function) to compute the probability of a given state vector.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2011
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

#ifndef REAK_GAUSSIAN_BELIEF_STATE_HPP
#define REAK_GAUSSIAN_BELIEF_STATE_HPP

#include "covariance_concept.hpp"
#include "belief_state_concept.hpp"

#include "lin_alg/mat_cholesky.hpp"
#include "lin_alg/mat_svd_method.hpp"
#include "lin_alg/mat_qr_decomp.hpp"

#include "base/named_object.hpp"

#include <boost/noncopyable.hpp>

#include "base/global_rng.hpp"

#include <boost/random/normal_distribution.hpp>

namespace ReaK {


namespace ctrl {
  
  
/**
 * This class template is a functor that can compute the probability that a given state is part of 
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * \tparam BeliefState The belief-state type to represent the Gaussian distribution of the state vector, should model the ContinuousBeliefStateConcept.
 * \tparam Storage The storage strategy of the Covariance matrix, this is used for specializing the PDF implementation for the most efficient way to compute the probabilities.
 */
template <typename BeliefState, covariance_storage::tag Storage = covariance_mat_traits< typename continuous_belief_state_traits<BeliefState>::covariance_type >::storage>
struct gaussian_pdf {
  typedef gaussian_pdf<BeliefState, Storage> self;
  
  typedef typename continuous_belief_state_traits<BeliefState>::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
    
  typedef typename covariance_mat_traits<covariance_type>::value_type scalar_type;
  typedef typename covariance_mat_traits<covariance_type>::size_type size_type;
  typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  
  state_type mean_state;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> L; 
  scalar_type factor;
  
  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  gaussian_pdf(const BeliefState& aBelief) : mean_state(aBelief.get_mean_state()), L(aBelief.get_covariance().size()), factor(-1) {
    const matrix_type& E = aBelief.get_covariance().get_matrix();
    try {
      decompose_Cholesky(E,L);
    } catch(singularity_error&) { return; };
    factor = scalar_type(1);
    for(size_type i = 0; i < L.get_row_count(); ++i)
      factor *= scalar_type(6.28318530718) * L(i,i);
  };
  
  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs) : mean_state(rhs.mean_state), L(rhs.L), factor(rhs.factor) { };
  
  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) { 
    using std::swap;
    swap(lhs.mean_state,rhs.mean_state);
    swap(lhs.L,rhs.L);
    swap(lhs.factor,rhs.factor);
  };
  
  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this,rhs);
    return *this;
  };
  
  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <typename Topology>
  scalar_type operator()(const state_type& v, const Topology& space) const {
    using std::sqrt;
    using std::exp;
    using ReaK::to_vect;
    
    if(factor <= scalar_type(0))
      return scalar_type(0);
      
    vect_n<scalar_type> d = to_vect<scalar_type>(space.difference(v, mean_state));
    mat< typename mat_traits<matrix_type>::value_type, mat_structure::rectangular> b(d.size(),1);
    for(size_type i = 0; i < d.size(); ++i) 
      b(i,0) = d[i];
    ::ReaK::detail::backsub_Cholesky_impl(L,b);
    scalar_type sum = scalar_type(0);
    for(size_type i = 0; i < d.size(); ++i) 
      sum += d[i] * b(i,0);
    return exp(scalar_type(-0.5) * sum) / sqrt(factor);
  };
  
  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) * (log( P.factor ) + scalar_type(P.L.get_row_count()));
  };
  
};


/**
 * This class template is a functor that can compute the probability that a given state is part of 
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * This class template specialization uses the fact that the covariance is represented as an information 
 * matrix in order to have a more efficient implementation.
 * \tparam BeliefState The belief-state type to represent the Gaussian distribution of the state vector, should model the ContinuousBeliefStateConcept.
 */
template <typename BeliefState>
struct gaussian_pdf<BeliefState, covariance_storage::information> {
  typedef gaussian_pdf<BeliefState, covariance_storage::information> self;
  
  typedef typename continuous_belief_state_traits<BeliefState>::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
    
  typedef typename covariance_mat_traits<covariance_type>::value_type scalar_type;
  typedef typename covariance_mat_traits<covariance_type>::size_type size_type;
  typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  
  state_type mean_state;
  matrix_type E_inv;
  scalar_type factor;
  
  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  gaussian_pdf(const BeliefState& aBelief) : mean_state(aBelief.get_mean_state()), E_inv(aBelief.get_covariance().get_inverse_matrix()), factor(-1) { 
    factor = determinant_Cholesky(E_inv);
    if(fabs(factor) < std::numeric_limits< scalar_type >::epsilon()) {
      factor = scalar_type(-1);
    } else {
      factor = scalar_type(1) / factor;
      for(size_type i = 0; i < E_inv.get_row_count(); ++i)
        factor *= scalar_type(6.28318530718);
    };
  };
  
  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs) : mean_state(rhs.mean_state), E_inv(rhs.E_inv), factor(rhs.factor) { };
  
  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) { 
    using std::swap;
    swap(lhs.mean_state,rhs.mean_state);
    swap(lhs.E_inv,rhs.E_inv);
    swap(lhs.factor,rhs.factor);
  };
  
  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this,rhs);
    return *this;
  };
  
  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <typename Topology>
  scalar_type operator()(const state_type& v, const Topology& space) const {
    using std::sqrt;
    using std::exp;
    using std::fabs;
    using ReaK::to_vect;
    typedef typename pp::topology_traits<Topology>::point_difference_type state_difference_type;
    BOOST_CONCEPT_ASSERT((ReadableVectorConcept<state_difference_type>));
    BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<covariance_type,state_difference_type>));
      
    if(factor <= scalar_type(0)) 
      return scalar_type(0);
    
    vect_n<scalar_type> d = to_vect<scalar_type>(space.difference(v, mean_state));
    return exp(scalar_type(-0.5) * (d * (E_inv * d))) / sqrt(factor);
  };
  
  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) * (log( P.factor ) + scalar_type(P.E_inv.get_row_count()));
  };
  
};



/**
 * This class template is a functor that can compute the probability that a given state is part of 
 * a gaussian belief-state. In other words, this is a probability distribution functor (PDF).
 * This class template specialization uses the fact that the covariance is represented as a decomposition of
 * the covariance matrix in order to have a more efficient implementation.
 * \tparam BeliefState The belief-state type to represent the Gaussian distribution of the state vector, should model the ContinuousBeliefStateConcept.
 */
template <typename BeliefState>
struct gaussian_pdf<BeliefState, covariance_storage::decomposed> {
  typedef gaussian_pdf<BeliefState, covariance_storage::decomposed> self;
  
  typedef typename continuous_belief_state_traits<BeliefState>::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
    
  typedef typename covariance_mat_traits<covariance_type>::value_type scalar_type;
  typedef typename covariance_mat_traits<covariance_type>::size_type size_type;
  typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  
  state_type mean_state;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> QX;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> RX;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> QY;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> RY;
  scalar_type factor;
  
  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  gaussian_pdf(const BeliefState& aBelief) : mean_state(aBelief.get_mean_state()),
                                             factor(-1) {
    decompose_QR(aBelief.get_covariance().get_covarying_block(),QX,RX);
    decompose_QR(aBelief.get_covariance().get_informing_inv_block(),QY,RY);
    
    factor = scalar_type(1);
    for(size_type i = 0; i < RX.get_row_count(); ++i)
      factor *= scalar_type(6.28318530718) * RX(i,i) / RY(i,i);
  };
  
  /**
   * Standard copy-constructor.
   */
  gaussian_pdf(const self& rhs) : mean_state(rhs.mean_state), QX(rhs.QX), RX(rhs.RX), QY(rhs.QY), RY(rhs.RY), factor(rhs.factor) { };
  
  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) { 
    using std::swap;
    swap(lhs.mean_state,rhs.mean_state);
    swap(lhs.QX,rhs.QX);
    swap(lhs.RX,rhs.RX);
    swap(lhs.QY,rhs.QY);
    swap(lhs.RY,rhs.RY);
    swap(lhs.factor,rhs.factor);
  };
  
  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this,rhs);
    return *this;
  };
  
  /**
   * The call-operator that computes the probability for a given state-vector.
   * \param v The state-vector for which the probability is sought.
   * \param space The state-space on which the state-vectors are distributed.
   * \return The probability of the given state-vector.
   */
  template <typename Topology>
  scalar_type operator()(const state_type& v, const Topology& space) const {
    using std::sqrt;
    using std::exp;
    using ReaK::to_vect;
    typedef typename pp::topology_traits<Topology>::point_difference_type state_difference_type;
    BOOST_CONCEPT_ASSERT((WritableVectorConcept<state_difference_type>));
    BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<covariance_type,state_difference_type>));
    
    if(factor <= scalar_type(0))
      return scalar_type(0);
    
    vect_n<scalar_type> d = to_vect<scalar_type>(space.difference(v, mean_state));
    vect_n<scalar_type> d_tmp = d * QX;  //QX^T d
    mat_vect_adaptor<state_difference_type> d_m(d_tmp);
    backsub_R(RX,d_m);
    return exp( scalar_type(-0.5) * ( d * ( QY * (RY * d_tmp) ) ) ) / sqrt(factor);
  };
  
  /**
   * This function computes the entropy of a gaussian probability distribution functor.
   * \param P The gaussian probability distribution functor.
   * \return The entropy of the pdf.
   */
  friend scalar_type entropy(const self& P) {
    using std::log;
    return scalar_type(0.5) * (log( P.factor ) + scalar_type(P.QX.get_row_count()));
  };
  
};


/*
template <typename Covariance, covariance_storage::tag Storage>
scalar_type KL_divergence(const gaussian_pdf<Covariance,Storage>& N0, 
                          const gaussian_pdf<Covariance,Storage>& N1) {
  using std::log;
  typedef typename covariance_mat_traits<Covariance>::point_type StateType;
  typedef typename state_vector_traits<StateType>::state_difference_type StateDiffType;
  typedef typename covariance_mat_traits<Covariance>::matrix_type MatType;
  StateDiffType d_mu = diff(N1.mean_state, N0.mean_state);
  mat_vect_adaptor<state_difference_type> d_mu_mat(d_mu);
  MatType S0(N0.covar.get_matrix());
  MatType S1inv(N1.covar.get_inverse_matrix());
  return scalar_type(0.5) * ( trace(S1inv * S0) - log(N1(N0.mean_state)) ) - entropy(N0);
};*/
    
/**
 * This function template computes the symmetric KL-divergence between two Gaussian probability 
 * distribution function objects.
 * \tparam BeliefState The belief-state type for the Gaussian PDFs.
 * \tparam Storage The storage strategy for the covariance matrices.
 * \tparam Topology The topology type of the underlying state representations.
 * \param N0 The first PDF.
 * \param N1 The second PDF.
 * \param space The topology of the underlying state representations.
 * \return The symmetric KL-divergence between the two Gaussian PDFs.
 */
template <typename BeliefState, covariance_storage::tag Storage, typename Topology>
typename gaussian_pdf<BeliefState,Storage>::scalar_type symKL_divergence(const gaussian_pdf<BeliefState,Storage>& N0, 
                                                                         const gaussian_pdf<BeliefState,Storage>& N1,
                                                                         const Topology& space) {
  using std::log;
  typedef typename gaussian_pdf<BeliefState,Storage>::scalar_type ScalarType;
  return ScalarType(-0.5) * log(N1(N0.mean_state, space) * N0(N1.mean_state, space)) - entropy(N0) - entropy(N1);
};
    



/**
 * This class template is a callable object (functor) which can generate random samples of 
 * state-vectors taken from a gaussian belief-state.
 * \tparam BeliefState The belief-state type to represent the Gaussian distribution of the state vector, should model the ContinuousBeliefStateConcept.
 */
template <typename BeliefState>
struct gaussian_sampler {
  typedef gaussian_sampler<BeliefState> self;
  
  typedef typename continuous_belief_state_traits<BeliefState>::state_type state_type;
  typedef typename continuous_belief_state_traits<BeliefState>::covariance_type covariance_type;
    
  typedef typename covariance_mat_traits<covariance_type>::value_type scalar_type;
  typedef typename covariance_mat_traits<covariance_type>::size_type size_type;
  typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
  
  BOOST_CONCEPT_ASSERT((ContinuousBeliefStateConcept<BeliefState>));
  
  state_type mean_state;
  mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> L;
  
  /**
   * Parametrized constructor.
   * \param aBelief The belief-state of the gaussian distribution.
   */
  gaussian_sampler(const BeliefState& aBelief) : mean_state(aBelief.get_mean_state()), L(aBelief.get_covariance().size()) {
    using std::sqrt;
    const matrix_type& C = aBelief.get_covariance().get_matrix();
    try {
      decompose_Cholesky(C,L);
    } catch(singularity_error&) { 
      mat< typename mat_traits<matrix_type>::value_type, mat_structure::diagonal> E(aBelief.get_covariance().size());
      mat< typename mat_traits<matrix_type>::value_type, mat_structure::square> U(aBelief.get_covariance().size()), V(aBelief.get_covariance().size());
      decompose_SVD(C,U,E,V);
      for(size_type i = 0; i < aBelief.get_covariance().size(); ++i)
        E(i,i) = sqrt(E(i,i));
      L = U * E;
    };
  };
  
  /**
   * Standard copy-constructor.
   */
  gaussian_sampler(const self& rhs) : mean_state(rhs.mean_state), L(rhs.L) { };
  
  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) { 
    using std::swap;
    swap(lhs.mean_state,rhs.mean_state);
    swap(lhs.L,rhs.L);
  };
  
  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this,rhs);
    return *this;
  };
  
  /**
   * The call-operator which can be used to generate a random state-sample from the gaussian probability distribution.
   * \param space The state-space on which the state-vectors are distributed.
   * \return A random state-sample from the gaussian probability distribution
   */
  template <typename Topology>
  state_type operator()(const Topology& space) const {
    using ReaK::to_vect;
    using ReaK::from_vect;
    
    boost::variate_generator< global_rng_type&, boost::normal_distribution<scalar_type> > var_rnd(get_global_rng(), boost::normal_distribution<scalar_type>());
    
    typedef typename pp::topology_traits<Topology>::point_difference_type state_difference_type;
    BOOST_CONCEPT_ASSERT((WritableVectorConcept<state_difference_type>));
    BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<covariance_type,state_difference_type>));
    
    vect_n<scalar_type> z = to_vect<scalar_type>(space.difference(mean_state, mean_state));
    for(size_type i = 0; i < z.size(); ++i)
      z[i] = var_rnd();
    
    return space.adjust(mean_state, from_vect<state_difference_type>(L * z));
  };
  
};




/**
 * This class template is used to represent a Gaussian belief-state, which is essentially a Gaussian
 * probability distribution which characterizes the estimation of a state-vector.
 * 
 * \tparam StateType The state vector type which represents the mean-state of the belief, should model 
 *         StateVectorConcept, by default it is the state-type associated with the covariance matrix type.
 * \tparam Covariance The covariance matrix type which represents the covariance of the state estimate, should 
 *         model the CovarianceMatrixConcept.
 */
template <typename StateType, typename Covariance>
class gaussian_belief_state : public virtual shared_object {
  public:
    typedef gaussian_belief_state<StateType, Covariance> self;
    
    typedef StateType state_type;
    typedef StateType state_difference_type;
    typedef Covariance covariance_type;
    
    typedef typename covariance_mat_traits<covariance_type>::value_type scalar_type;
    typedef typename covariance_mat_traits<covariance_type>::size_type size_type;
    typedef typename covariance_mat_traits<covariance_type>::matrix_type matrix_type;
    
    typedef gaussian_pdf<self> pdf_type;
    typedef gaussian_sampler<self> random_sampler_type;
    
    BOOST_CONCEPT_ASSERT((CovarianceMatrixConcept<Covariance,state_difference_type>));
    
    BOOST_STATIC_CONSTANT(belief_distribution::tag, distribution = belief_distribution::unimodal);
    BOOST_STATIC_CONSTANT(belief_representation::tag, representation = belief_representation::gaussian);
    
  private:
    state_type mean_state;
    covariance_type covar;
    
  public:
    
    /**
     * Returns the probability distribution functor associated with this belief-state's probability distribution.
     * \return The probability distribution functor associated with this belief-state's probability distribution.
     */
    pdf_type get_pdf() const { return pdf_type(*this); };
    
    /**
     * Returns the most-likely state (i.e. the mean-state for a Gaussian distribution).
     * \return The most-likely state.
     */
    const state_type& get_most_likely_state() const { return mean_state; };
    
    /**
     * Returns the random sampler functor associated with this belief-state's probability distribution.
     * \return The random sampler functor associated with this belief-state's probability distribution.
     */
    random_sampler_type get_random_sampler() const { 
      return random_sampler_type(*this);
    };
    
    /**
     * Returns the mean-state state.
     * \return The mean-state.
     */
    const state_type& get_mean_state() const { return mean_state; };
    /**
     * Returns the covariance.
     * \return The covariance.
     */
    const covariance_type& get_covariance() const { return covar; };
    
    /**
     * Sets the mean-state.
     * \param aMeanState The new mean-state for this gaussian belief-state.
     */
    void set_mean_state(const state_type& aMeanState) { mean_state = aMeanState; };
    /**
     * Sets the covariance.
     * \param aCov The new covariance for this gaussian belief-state.
     */
    void set_covariance(const covariance_type& aCov) { covar = aCov; };
    
    /**
     * Returns the size of the covariance matrix of this gaussian belief-state.
     */
    size_type size() const { return covar.size(); };
    
    /**
     * Parametrized and default constructor.
     * \param aMeanState The mean-state of the gaussian belief-state.
     * \param aCov The covariance of the gaussian belief-state.
     */
    gaussian_belief_state(const state_type& aMeanState = state_type(), 
                          const covariance_type& aCov = covariance_type()) : 
                          mean_state(aMeanState), 
                          covar(aCov) { };

    /**
     * Standard copy-constructor.
     */
    gaussian_belief_state(const self& rhs) : mean_state(rhs.mean_state), covar(rhs.covar) { };
    
    /**
     * Standard swap function.
     */
    friend void swap(self& lhs, self& rhs) {
      using std::swap;
      swap(lhs.mean_state,rhs.mean_state);
      swap(lhs.covar,rhs.covar);
    };
    
    /**
     * Standard assignment operator.
     */
    self& operator=(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::shared_object::save(aA,ReaK::shared_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(mean_state)
         & RK_SERIAL_SAVE_WITH_NAME(covar);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::shared_object::load(aA,ReaK::shared_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(mean_state)
         & RK_SERIAL_LOAD_WITH_NAME(covar);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300010,1,"gaussian_belief_state",shared_object)
    
    
};

template <typename Covariance, typename StateType>
struct is_belief_state< gaussian_belief_state<Covariance,StateType> > {
  BOOST_STATIC_CONSTANT(bool, value = true);
  typedef is_belief_state< gaussian_belief_state<Covariance,StateType> > type;
};

template <typename Covariance, typename StateType>
struct is_continuous_belief_state< gaussian_belief_state<Covariance,StateType> > {
  BOOST_STATIC_CONSTANT(bool, value = true);
  typedef is_continuous_belief_state< gaussian_belief_state<Covariance,StateType> > type;
};


/**
 * This function computes the symmetric KL-divergence between two belief-states.
 * \tparam StateType The state-type for the Gaussian belief-states.
 * \tparam Covariance The covariance matrix type for the Gaussian belief-states, should model the CovarianceMatrixConcept.
 * \tparam Topology The belief-space on which the belief-states lie, should model the CovarianceMatrixConcept.
 * \param P The first Gaussian belief-state.
 * \param Q The second Gaussian belief-state.
 * \return The symmetric KL-divergence between the two Gaussian belief-states.
 */
template <typename StateType, typename Covariance, typename Topology>
typename gaussian_belief_state<StateType,Covariance>::scalar_type 
 symKL_divergence(const gaussian_belief_state<StateType,Covariance>& P, 
                  const gaussian_belief_state<StateType,Covariance>& Q,
                  const Topology& space) {
  return symKL_divergence(P.get_pdf(),Q.get_pdf(),space.get_state_topology());
};

/**
 * This function computes the symmetric KL-divergence between two belief-states.
 * \tparam StateType The state-type for the Gaussian belief-state.
 * \tparam Covariance The covariance matrix type for the Gaussian belief-state, should model the CovarianceMatrixConcept.
 * \param P The Gaussian belief-state.
 * \return The entropy of the Gaussian belief-state.
 */
template <typename StateType, typename Covariance>
typename gaussian_belief_state<StateType,Covariance>::scalar_type 
 entropy(const gaussian_belief_state<StateType,Covariance>& P) {
  return entropy(P.get_pdf());
};








};

};

#endif








