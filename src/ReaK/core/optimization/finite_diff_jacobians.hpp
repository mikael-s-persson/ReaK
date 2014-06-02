/**
 * \file finite_diff_jacobians.hpp
 *
 * The following library provides implementations of finite difference methods for evaluating a Jacobian
 * for a function of multiple variables and multiple outputs (or single).
 *
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2011
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

#ifndef REAK_FINITE_DIFF_JACOBIANS_HPP
#define REAK_FINITE_DIFF_JACOBIANS_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/lin_alg/mat_alg.hpp>

#include "optim_exceptions.hpp"

#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/bind.hpp>


namespace ReaK {
  
  
namespace optim {
  
  
namespace detail {
  
  

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_2pts_forward_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);
  
  Vector2 y1 = y;

  for(SizeType j=0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * fabs(x[j]);
    if(d < delta)
      d = delta;

    ValueType tmp = x[j];
    x[j] += d;

    y1 = f(x);

    x[j] = tmp; /* restore */

    slice(jac)(range(SizeType(0),M),j) = (y1 - y) * (1.0 / d);
    
  };
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_2pts_central_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);

  Vector2 y_prev = y;
  Vector2 y_next = y;

  for(SizeType j=0; j < N; ++j) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d= ValueType(1E-04) * fabs(x[j]);
    if( d < delta )
      d = delta;

    ValueType tmp = x[j];
    x[j] -= d;
    y_prev = f(x);
    
    x[j] = tmp + d;
    y_next = f(x);
    x[j]=tmp; /* restore */
    
    slice(jac)(range(SizeType(0),M),j) = (y_next - y_prev) * (0.5 / d);
  };
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
void compute_jacobian_5pts_central_impl(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta) {
  typedef typename vect_traits<Vector1>::value_type ValueType;
  typedef typename vect_traits<Vector1>::size_type SizeType;
  using std::fabs;
  
  SizeType N = x.size();
  SizeType M = y.size();
  if(jac.get_row_count() != M) 
    jac.set_row_count(M);
  if(jac.get_col_count() != N)
    jac.set_col_count(N);
  
  Vector2 y0 = y;
  Vector2 y1 = y;
  Vector2 y2 = y;
  Vector2 y3 = y;
  
  for(SizeType i=0; i < N; ++i) {
    /* determine d=max(1E-04*|p[j]|, delta), see HZ */
    ValueType d = ValueType(1E-04) * fabs(x[i]);
    if( d < delta )
      d = delta;
    
    ValueType tmp = x[i];
    x[i] -= 2.0 * d;
    y0 = f(x);
    x[i] += d;
    y1 = f(x);
    x[i] += 2.0 * d;
    y2 = f(x);
    x[i] += d;
    y3 = f(x);
    
    x[i] = tmp; // restore
    
    slice(jac)(range(SizeType(0),M),i) = (y0 - 8.0 * (y1 - y2) - y3) * (1.0 / (12.0 * d));
  };
};



template <typename Function, typename Vector, typename Scalar>
vect<Scalar,1> scalar_return_function_to_vect_function(Function f, const Vector& x) {
  vect<Scalar,1> result;
  result[0] = f(x);
  return result;
};

template <typename Function, typename Scalar, typename Vector>
Vector scalar_param_function_to_vect_function(Function f, const vect<Scalar,1>& x) {
  return f(x[0]);
};

template <typename Function, typename Scalar1, typename Scalar2>
vect<Scalar2,1> scalar_param_ret_function_to_vect_function(Function f, const vect<Scalar1,1>& x) {
  vect<Scalar2,1> result;
  result[0] = f(x[0]);
  return result;
};




};




template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    is_fully_writable_matrix<Matrix> 
  >,
void >::type compute_jacobian_2pts_forward(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  detail::compute_jacobian_2pts_forward_impl(f,x,y,jac,delta);
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    boost::mpl::not_< is_fully_writable_matrix<Matrix> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_forward(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  mat< typename vect_traits<Vector1>::value_type, mat_structure::rectangular> jac_tmp(y.size(), x.size());
  detail::compute_jacobian_2pts_forward_impl(f,x,y,jac_tmp,delta);
  jac = jac_tmp;
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_forward(Function f, Vector& x, const Scalar& y, Matrix& jac, typename vect_traits<Vector>::value_type delta = typename vect_traits<Vector>::value_type(1e-6)) {
  vect< Scalar, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_2pts_forward( boost::bind< vect< Scalar, 1> >(detail::scalar_return_function_to_vect_function<Function,Vector,Scalar>,f,_1),x,y_tmp,jac,delta);
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_forward(Function f, Scalar& x, const Vector& y, Matrix& jac, Scalar delta = Scalar(1e-6)) {
  vect< Scalar, 1> x_tmp; x_tmp[0] = x;
  compute_jacobian_2pts_forward( boost::bind< vect< Scalar, 1> >(detail::scalar_param_function_to_vect_function<Function,Scalar,Vector>,f,_1),x_tmp,y,jac,delta);
  x = x_tmp[0];
};


template <typename Function, typename Scalar1, typename Scalar2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    boost::mpl::not_< is_readable_vector<Scalar1> >,
    boost::mpl::not_< is_readable_vector<Scalar2> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_forward(Function f, Scalar1& x, const Scalar2& y, Matrix& jac, Scalar1 delta = Scalar1(1e-6)) {
  vect< Scalar1, 1> x_tmp; x_tmp[0] = x;
  vect< Scalar2, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_2pts_forward( boost::bind< vect< Scalar2, 1> >(detail::scalar_param_ret_function_to_vect_function<Function,Scalar1,Scalar2>,f,_1),x_tmp,y_tmp,jac,delta);
  x = x_tmp[0];
};








template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    is_fully_writable_matrix<Matrix> 
  >,
void >::type compute_jacobian_2pts_central(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  detail::compute_jacobian_2pts_central_impl(f,x,y,jac,delta);
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    boost::mpl::not_< is_fully_writable_matrix<Matrix> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_central(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  mat< typename vect_traits<Vector1>::value_type, mat_structure::rectangular> jac_tmp(y.size(), x.size());
  detail::compute_jacobian_2pts_central_impl(f,x,y,jac_tmp,delta);
  jac = jac_tmp;
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_central(Function f, Vector& x, const Scalar& y, Matrix& jac, typename vect_traits<Vector>::value_type delta = typename vect_traits<Vector>::value_type(1e-6)) {
  vect< Scalar, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_2pts_central( boost::bind< vect< Scalar, 1> >(detail::scalar_return_function_to_vect_function<Function,Vector,Scalar>,f,_1),x,y_tmp,jac,delta);
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_central(Function f, Scalar& x, const Vector& y, Matrix& jac, Scalar delta = Scalar(1e-6)) {
  vect< Scalar, 1> x_tmp; x_tmp[0] = x;
  compute_jacobian_2pts_central( boost::bind< vect< Scalar, 1> >(detail::scalar_param_function_to_vect_function<Function,Scalar,Vector>,f,_1),x_tmp,y,jac,delta);
  x = x_tmp[0];
};


template <typename Function, typename Scalar1, typename Scalar2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    boost::mpl::not_< is_readable_vector<Scalar1> >,
    boost::mpl::not_< is_readable_vector<Scalar2> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_2pts_central(Function f, Scalar1& x, const Scalar2& y, Matrix& jac, Scalar1 delta = Scalar1(1e-6)) {
  vect< Scalar1, 1> x_tmp; x_tmp[0] = x;
  vect< Scalar2, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_2pts_central( boost::bind< vect< Scalar2, 1> >(detail::scalar_param_ret_function_to_vect_function<Function,Scalar1,Scalar2>,f,_1),x_tmp,y_tmp,jac,delta);
  x = x_tmp[0];
};






template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    is_fully_writable_matrix<Matrix> 
  >,
void >::type compute_jacobian_5pts_central(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  detail::compute_jacobian_5pts_central_impl(f,x,y,jac,delta);
};

template <typename Function, typename Vector1, typename Vector2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector1>,
    is_readable_vector<Vector2>,
    boost::mpl::not_< is_fully_writable_matrix<Matrix> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_5pts_central(Function f, Vector1& x, const Vector2& y, Matrix& jac, typename vect_traits<Vector1>::value_type delta = typename vect_traits<Vector1>::value_type(1e-6)) {
  mat< typename vect_traits<Vector1>::value_type, mat_structure::rectangular> jac_tmp(y.size(), x.size());
  detail::compute_jacobian_5pts_central_impl(f,x,y,jac_tmp,delta);
  jac = jac_tmp;
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_5pts_central(Function f, Vector& x, const Scalar& y, Matrix& jac, typename vect_traits<Vector>::value_type delta = typename vect_traits<Vector>::value_type(1e-6)) {
  vect< Scalar, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_5pts_central( boost::bind< vect< Scalar, 1> >(detail::scalar_return_function_to_vect_function<Function,Vector,Scalar>,f,_1),x,y_tmp,jac,delta);
};

template <typename Function, typename Vector, typename Scalar, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    is_readable_vector<Vector>,
    boost::mpl::not_< is_readable_vector<Scalar> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_5pts_central(Function f, Scalar& x, const Vector& y, Matrix& jac, Scalar delta = Scalar(1e-6)) {
  vect< Scalar, 1> x_tmp; x_tmp[0] = x;
  compute_jacobian_5pts_central( boost::bind< vect< Scalar, 1> >(detail::scalar_param_function_to_vect_function<Function,Scalar,Vector>,f,_1),x_tmp,y,jac,delta);
  x = x_tmp[0];
};


template <typename Function, typename Scalar1, typename Scalar2, typename Matrix>
typename boost::enable_if< 
  boost::mpl::and_<
    boost::mpl::not_< is_readable_vector<Scalar1> >,
    boost::mpl::not_< is_readable_vector<Scalar2> >,
    is_writable_matrix<Matrix>
  >,
void >::type compute_jacobian_5pts_central(Function f, Scalar1& x, const Scalar2& y, Matrix& jac, Scalar1 delta = Scalar1(1e-6)) {
  vect< Scalar1, 1> x_tmp; x_tmp[0] = x;
  vect< Scalar2, 1> y_tmp; y_tmp[0] = y;
  compute_jacobian_5pts_central( boost::bind< vect< Scalar2, 1> >(detail::scalar_param_ret_function_to_vect_function<Function,Scalar1,Scalar2>,f,_1),x_tmp,y_tmp,jac,delta);
  x = x_tmp[0];
};







};

};

#endif








