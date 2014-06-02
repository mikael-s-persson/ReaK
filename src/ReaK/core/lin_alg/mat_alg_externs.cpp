
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

#include <ReaK/core/base/defs.hpp>

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/core/lin_alg/mat_alg.hpp>

namespace ReaK {

template class mat<double, mat_structure::rectangular, mat_alignment::column_major>;
template class mat<double, mat_structure::rectangular, mat_alignment::row_major>;

template class mat<float, mat_structure::rectangular, mat_alignment::column_major>;
template class mat<float, mat_structure::rectangular, mat_alignment::row_major>;

template mat<double, mat_structure::rectangular, mat_alignment::column_major>& mat<double, mat_structure::rectangular, mat_alignment::column_major>::operator +=(const mat<double, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<double, mat_structure::rectangular, mat_alignment::column_major>& mat<double, mat_structure::rectangular, mat_alignment::column_major>::operator -=(const mat<double, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<double, mat_structure::rectangular, mat_alignment::column_major>& mat<double, mat_structure::rectangular, mat_alignment::column_major>::operator *=(const mat<double, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<double, mat_structure::rectangular, mat_alignment::row_major>& mat<double, mat_structure::rectangular, mat_alignment::row_major>::operator +=(const mat<double, mat_structure::rectangular, mat_alignment::row_major>& M);
template mat<double, mat_structure::rectangular, mat_alignment::row_major>& mat<double, mat_structure::rectangular, mat_alignment::row_major>::operator -=(const mat<double, mat_structure::rectangular, mat_alignment::row_major>& M);
template mat<double, mat_structure::rectangular, mat_alignment::row_major>& mat<double, mat_structure::rectangular, mat_alignment::row_major>::operator *=(const mat<double, mat_structure::rectangular, mat_alignment::row_major>& M);

template mat<float, mat_structure::rectangular, mat_alignment::column_major>& mat<float, mat_structure::rectangular, mat_alignment::column_major>::operator +=(const mat<float, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<float, mat_structure::rectangular, mat_alignment::column_major>& mat<float, mat_structure::rectangular, mat_alignment::column_major>::operator -=(const mat<float, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<float, mat_structure::rectangular, mat_alignment::column_major>& mat<float, mat_structure::rectangular, mat_alignment::column_major>::operator *=(const mat<float, mat_structure::rectangular, mat_alignment::column_major>& M);
template mat<float, mat_structure::rectangular, mat_alignment::row_major>& mat<float, mat_structure::rectangular, mat_alignment::row_major>::operator +=(const mat<float, mat_structure::rectangular, mat_alignment::row_major>& M);
template mat<float, mat_structure::rectangular, mat_alignment::row_major>& mat<float, mat_structure::rectangular, mat_alignment::row_major>::operator -=(const mat<float, mat_structure::rectangular, mat_alignment::row_major>& M);
template mat<float, mat_structure::rectangular, mat_alignment::row_major>& mat<float, mat_structure::rectangular, mat_alignment::row_major>::operator *=(const mat<float, mat_structure::rectangular, mat_alignment::row_major>& M);



template class mat<double, mat_structure::square, mat_alignment::column_major>;
template class mat<double, mat_structure::square, mat_alignment::row_major>;

template class mat<float, mat_structure::square, mat_alignment::column_major>;
template class mat<float, mat_structure::square, mat_alignment::row_major>;

template mat<double, mat_structure::square, mat_alignment::column_major>& mat<double, mat_structure::square, mat_alignment::column_major>::operator +=(const mat<double, mat_structure::square, mat_alignment::column_major>& M);
template mat<double, mat_structure::square, mat_alignment::column_major>& mat<double, mat_structure::square, mat_alignment::column_major>::operator -=(const mat<double, mat_structure::square, mat_alignment::column_major>& M);
template mat<double, mat_structure::square, mat_alignment::column_major>& mat<double, mat_structure::square, mat_alignment::column_major>::operator *=(const mat<double, mat_structure::square, mat_alignment::column_major>& M);
template mat<double, mat_structure::square, mat_alignment::row_major>& mat<double, mat_structure::square, mat_alignment::row_major>::operator +=(const mat<double, mat_structure::square, mat_alignment::row_major>& M);
template mat<double, mat_structure::square, mat_alignment::row_major>& mat<double, mat_structure::square, mat_alignment::row_major>::operator -=(const mat<double, mat_structure::square, mat_alignment::row_major>& M);
template mat<double, mat_structure::square, mat_alignment::row_major>& mat<double, mat_structure::square, mat_alignment::row_major>::operator *=(const mat<double, mat_structure::square, mat_alignment::row_major>& M);

template mat<float, mat_structure::square, mat_alignment::column_major>& mat<float, mat_structure::square, mat_alignment::column_major>::operator +=(const mat<float, mat_structure::square, mat_alignment::column_major>& M);
template mat<float, mat_structure::square, mat_alignment::column_major>& mat<float, mat_structure::square, mat_alignment::column_major>::operator -=(const mat<float, mat_structure::square, mat_alignment::column_major>& M);
template mat<float, mat_structure::square, mat_alignment::column_major>& mat<float, mat_structure::square, mat_alignment::column_major>::operator *=(const mat<float, mat_structure::square, mat_alignment::column_major>& M);
template mat<float, mat_structure::square, mat_alignment::row_major>& mat<float, mat_structure::square, mat_alignment::row_major>::operator +=(const mat<float, mat_structure::square, mat_alignment::row_major>& M);
template mat<float, mat_structure::square, mat_alignment::row_major>& mat<float, mat_structure::square, mat_alignment::row_major>::operator -=(const mat<float, mat_structure::square, mat_alignment::row_major>& M);
template mat<float, mat_structure::square, mat_alignment::row_major>& mat<float, mat_structure::square, mat_alignment::row_major>::operator *=(const mat<float, mat_structure::square, mat_alignment::row_major>& M);





template class mat<double, mat_structure::symmetric>;
template class mat<float, mat_structure::symmetric>;


template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator +=(const mat<double,mat_structure::symmetric>& M);
template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator +=(const mat<double,mat_structure::diagonal>& M);
template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator -=(const mat<double,mat_structure::symmetric>& M);
template mat<double,mat_structure::symmetric>& mat<double,mat_structure::symmetric>::operator -=(const mat<double,mat_structure::diagonal>& M);

template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator +=(const mat<float,mat_structure::symmetric>& M);
template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator +=(const mat<float,mat_structure::diagonal>& M);
template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator -=(const mat<float,mat_structure::symmetric>& M);
template mat<float,mat_structure::symmetric>& mat<float,mat_structure::symmetric>::operator -=(const mat<float,mat_structure::diagonal>& M);



template class mat<double, mat_structure::nil>;
template class mat<float, mat_structure::nil>;

template mat<double,mat_structure::nil> mat_nil<double>(mat<double,mat_structure::nil>::size_type aRowCount, mat<double,mat_structure::nil>::size_type aColCount);
template mat<float,mat_structure::nil> mat_nil<float>(mat<float,mat_structure::nil>::size_type aRowCount, mat<float,mat_structure::nil>::size_type aColCount);

template vect<double,2> operator *< double, 2 >(const mat<double,mat_structure::nil>& M,const vect<double,2>& V);
template vect<double,3> operator *< double, 3 >(const mat<double,mat_structure::nil>& M,const vect<double,3>& V);
template vect<double,4> operator *< double, 4 >(const mat<double,mat_structure::nil>& M,const vect<double,4>& V);
template vect<double,6> operator *< double, 6 >(const mat<double,mat_structure::nil>& M,const vect<double,6>& V);
template vect_n<double> operator *< double, vect_n<double> >(const mat<double,mat_structure::nil>& M, const vect_n<double>& V);

template vect<double,2> operator *< double, 2 >(const vect<double,2>& V,const mat<double,mat_structure::nil>& M);
template vect<double,3> operator *< double, 3 >(const vect<double,3>& V,const mat<double,mat_structure::nil>& M);
template vect<double,4> operator *< double, 4 >(const vect<double,4>& V,const mat<double,mat_structure::nil>& M);
template vect<double,6> operator *< double, 6 >(const vect<double,6>& V,const mat<double,mat_structure::nil>& M);
template vect_n<double> operator *< double, vect_n<double> >(const vect_n<double>& V,const mat<double,mat_structure::nil>& M);

template vect<float,2> operator *< float, 2 >(const mat<float,mat_structure::nil>& M,const vect<float,2>& V);
template vect<float,3> operator *< float, 3 >(const mat<float,mat_structure::nil>& M,const vect<float,3>& V);
template vect<float,4> operator *< float, 4 >(const mat<float,mat_structure::nil>& M,const vect<float,4>& V);
template vect<float,6> operator *< float, 6 >(const mat<float,mat_structure::nil>& M,const vect<float,6>& V);
template vect_n<float> operator *< float, vect_n<float> >(const mat<float,mat_structure::nil>& M, const vect_n<float>& V);

template vect<float,2> operator *< float, 2 >(const vect<float,2>& V,const mat<float,mat_structure::nil>& M);
template vect<float,3> operator *< float, 3 >(const vect<float,3>& V,const mat<float,mat_structure::nil>& M);
template vect<float,4> operator *< float, 4 >(const vect<float,4>& V,const mat<float,mat_structure::nil>& M);
template vect<float,6> operator *< float, 6 >(const vect<float,6>& V,const mat<float,mat_structure::nil>& M);
template vect_n<float> operator *< float, vect_n<float> >(const vect_n<float>& V,const mat<float,mat_structure::nil>& M);




template class mat<double, mat_structure::identity>;
template class mat<float, mat_structure::identity>;

template mat<double,mat_structure::identity> mat_ident<double>(mat<double,mat_structure::identity>::size_type aRowCount);
template mat<float,mat_structure::identity> mat_ident<float>(mat<float,mat_structure::identity>::size_type aRowCount);

template vect<double,2> operator *<double, vect<double,2>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,2>& V);
template vect<double,3> operator *<double, vect<double,3>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,3>& V);
template vect<double,4> operator *<double, vect<double,4>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,4>& V);
template vect<double,6> operator *<double, vect<double,6>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M,const vect<double,6>& V);
template vect_n<double> operator *<double, vect_n<double>, mat_alignment::column_major, std::allocator<double> >(const mat<double,mat_structure::identity>& M, const vect_n<double>& V);

template vect<double,2> operator *<double, vect<double,2>, mat_alignment::column_major, std::allocator<double> >(const vect<double,2>& V,const mat<double,mat_structure::identity>& M);
template vect<double,3> operator *<double, vect<double,3>, mat_alignment::column_major, std::allocator<double> >(const vect<double,3>& V,const mat<double,mat_structure::identity>& M);
template vect<double,4> operator *<double, vect<double,4>, mat_alignment::column_major, std::allocator<double> >(const vect<double,4>& V,const mat<double,mat_structure::identity>& M);
template vect<double,6> operator *<double, vect<double,6>, mat_alignment::column_major, std::allocator<double> >(const vect<double,6>& V,const mat<double,mat_structure::identity>& M);
template vect_n<double> operator *<double, vect_n<double>, mat_alignment::column_major, std::allocator<double> >(const vect_n<double>& V,const mat<double,mat_structure::identity>& M);

template vect<float,2> operator *<float, vect<float,2>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,2>& V);
template vect<float,3> operator *<float, vect<float,3>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,3>& V);
template vect<float,4> operator *<float, vect<float,4>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,4>& V);
template vect<float,6> operator *<float, vect<float,6>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M,const vect<float,6>& V);
template vect_n<float> operator *<float, vect_n<float>, mat_alignment::column_major, std::allocator<float> >(const mat<float,mat_structure::identity>& M, const vect_n<float>& V);

template vect<float,2> operator *<float, vect<float,2>, mat_alignment::column_major, std::allocator<float> >(const vect<float,2>& V,const mat<float,mat_structure::identity>& M);
template vect<float,3> operator *<float, vect<float,3>, mat_alignment::column_major, std::allocator<float> >(const vect<float,3>& V,const mat<float,mat_structure::identity>& M);
template vect<float,4> operator *<float, vect<float,4>, mat_alignment::column_major, std::allocator<float> >(const vect<float,4>& V,const mat<float,mat_structure::identity>& M);
template vect<float,6> operator *<float, vect<float,6>, mat_alignment::column_major, std::allocator<float> >(const vect<float,6>& V,const mat<float,mat_structure::identity>& M);
template vect_n<float> operator *<float, vect_n<float>, mat_alignment::column_major, std::allocator<float> >(const vect_n<float>& V,const mat<float,mat_structure::identity>& M);


template mat<double,mat_structure::scalar> operator *(const mat<double,mat_structure::identity>& M, const double& S);
template mat<double,mat_structure::scalar> operator *(const double& S, const mat<double,mat_structure::identity>& M);

template mat<float,mat_structure::scalar> operator *(const mat<float,mat_structure::identity>& M, const float& S);
template mat<float,mat_structure::scalar> operator *(const float& S, const mat<float,mat_structure::identity>& M);




template class mat<double, mat_structure::scalar>;
template class mat<float, mat_structure::scalar>;



template class mat<double, mat_structure::skew_symmetric>;
template class mat<float, mat_structure::skew_symmetric>;



template class mat<double, mat_structure::diagonal>;
template class mat<float, mat_structure::diagonal>;


template class mat<double, mat_structure::permutation>;
template class mat<float, mat_structure::permutation>;


template vect<double,2> operator *< double, vect<double,2> >(const mat<double,mat_structure::permutation>& M, const vect<double,2>& V);
template vect<double,3> operator *< double, vect<double,3> >(const mat<double,mat_structure::permutation>& M, const vect<double,3>& V);
template vect<double,4> operator *< double, vect<double,4> >(const mat<double,mat_structure::permutation>& M, const vect<double,4>& V);
template vect<double,6> operator *< double, vect<double,6> >(const mat<double,mat_structure::permutation>& M, const vect<double,6>& V);
template vect_n<double> operator *< double, vect_n<double> >(const mat<double,mat_structure::permutation>& M, const vect_n<double>& V);

template vect<double,2> operator *< double, vect<double,2> >(const vect<double,2>& V, const mat<double,mat_structure::permutation>& M);
template vect<double,3> operator *< double, vect<double,3> >(const vect<double,3>& V, const mat<double,mat_structure::permutation>& M);
template vect<double,4> operator *< double, vect<double,4> >(const vect<double,4>& V, const mat<double,mat_structure::permutation>& M);
template vect<double,6> operator *< double, vect<double,6> >(const vect<double,6>& V, const mat<double,mat_structure::permutation>& M);
template vect_n<double> operator *< double, vect_n<double> >(const vect_n<double>& V, const mat<double,mat_structure::permutation>& M);

template mat<double,mat_structure::square> operator *(const mat<double,mat_structure::permutation>& M, const double& S);
template mat<double,mat_structure::square> operator *(const double& S,const mat<double,mat_structure::permutation>& M);

template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::rectangular>& M1, const mat<double,mat_structure::permutation>& M2);
template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::square>& M1, const mat<double,mat_structure::permutation>& M2);

template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::permutation>& M1, const mat<double,mat_structure::rectangular>& M2);
template mat<double,mat_structure::rectangular> operator *(const mat<double,mat_structure::permutation>& M1, const mat<double,mat_structure::square>& M2);


template vect<float,2> operator *< float, vect<float,2> >(const mat<float,mat_structure::permutation>& M, const vect<float,2>& V);
template vect<float,3> operator *< float, vect<float,3> >(const mat<float,mat_structure::permutation>& M, const vect<float,3>& V);
template vect<float,4> operator *< float, vect<float,4> >(const mat<float,mat_structure::permutation>& M, const vect<float,4>& V);
template vect<float,6> operator *< float, vect<float,6> >(const mat<float,mat_structure::permutation>& M, const vect<float,6>& V);
template vect_n<float> operator *< float, vect_n<float> >(const mat<float,mat_structure::permutation>& M, const vect_n<float>& V);

template vect<float,2> operator *< float, vect<float,2> >(const vect<float,2>& V, const mat<float,mat_structure::permutation>& M);
template vect<float,3> operator *< float, vect<float,3> >(const vect<float,3>& V, const mat<float,mat_structure::permutation>& M);
template vect<float,4> operator *< float, vect<float,4> >(const vect<float,4>& V, const mat<float,mat_structure::permutation>& M);
template vect<float,6> operator *< float, vect<float,6> >(const vect<float,6>& V, const mat<float,mat_structure::permutation>& M);
template vect_n<float> operator *< float, vect_n<float> >(const vect_n<float>& V, const mat<float,mat_structure::permutation>& M);

template mat<float,mat_structure::square> operator *(const mat<float,mat_structure::permutation>& M, const float& S);
template mat<float,mat_structure::square> operator *(const float& S,const mat<float,mat_structure::permutation>& M);

template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::rectangular>& M1, const mat<float,mat_structure::permutation>& M2);
template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::square>& M1, const mat<float,mat_structure::permutation>& M2);

template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::permutation>& M1, const mat<float,mat_structure::rectangular>& M2);
template mat<float,mat_structure::rectangular> operator *(const mat<float,mat_structure::permutation>& M1, const mat<float,mat_structure::square>& M2);















template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::square,mat_alignment::column_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::square,mat_alignment::row_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::symmetric>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::skew_symmetric>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::nil>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::identity>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::scalar>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::diagonal>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<double,mat_structure::permutation>& M);

template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);


//rect and square matrix.
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);

//template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::symmetric>& M2);
//template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::symmetric>& M2);

//template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::skew_symmetric>& M2);
//template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::skew_symmetric>& M2);


//square and rect matrix.
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);

//template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *< mat<double,mat_structure::rectangular,mat_alignment::column_major> >(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
//template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *< mat<double,mat_structure::rectangular,mat_alignment::row_major> >(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);

//template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *< mat<double,mat_structure::rectangular,mat_alignment::column_major> >(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
//template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *< mat<double,mat_structure::rectangular,mat_alignment::row_major> >(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);


//square and square matrix.
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::skew_symmetric>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::symmetric>& M2);

template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::symmetric>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::skew_symmetric>& M2);


//dense and diagonal matrix.
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::diagonal>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::diagonal>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::diagonal>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::diagonal>& M2);
template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::diagonal>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::diagonal>& M2);

//diagonal and dense matrix.
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::symmetric>& M2);
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::skew_symmetric>& M2);

//diagonal and diagonal matrix.
template mat_product_result< mat<double,mat_structure::diagonal>, mat<double,mat_structure::diagonal> >::type operator *(const mat<double,mat_structure::diagonal>& M1, const mat<double,mat_structure::diagonal>& M2);


//dense and nil matrix.
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::nil>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::nil>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::nil>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::nil>& M2);
template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::nil>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::nil>& M2);

//nil and dense matrix.
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::symmetric>& M2);
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::skew_symmetric>& M2);

//nil and nil matrix.
template mat_product_result< mat<double,mat_structure::nil>, mat<double,mat_structure::nil> >::type operator *(const mat<double,mat_structure::nil>& M1, const mat<double,mat_structure::nil>& M2);


//dense and scalar matrix.
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::scalar>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::scalar>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::scalar>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::scalar>& M2);
template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::scalar>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::scalar>& M2);

//scalar and dense matrix.
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::symmetric>& M2);
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::skew_symmetric>& M2);

//scalar and scalar matrix.
template mat_product_result< mat<double,mat_structure::scalar>, mat<double,mat_structure::scalar> >::type operator *(const mat<double,mat_structure::scalar>& M1, const mat<double,mat_structure::scalar>& M2);


//dense and identity matrix.
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::column_major>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<double,mat_structure::identity>& M2);
template mat_product_result< mat<double,mat_structure::rectangular,mat_alignment::row_major>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<double,mat_structure::identity>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::column_major>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::square,mat_alignment::column_major>& M1, const mat<double,mat_structure::identity>& M2);
template mat_product_result< mat<double,mat_structure::square,mat_alignment::row_major>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::square,mat_alignment::row_major>& M1, const mat<double,mat_structure::identity>& M2);
template mat_product_result< mat<double,mat_structure::symmetric>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::symmetric>& M1, const mat<double,mat_structure::identity>& M2);
template mat_product_result< mat<double,mat_structure::skew_symmetric>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::skew_symmetric>& M1, const mat<double,mat_structure::identity>& M2);

//identity and dense matrix.
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::symmetric> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::symmetric>& M2);
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::skew_symmetric> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::skew_symmetric>& M2);

//identity and identity matrix.
template mat_product_result< mat<double,mat_structure::identity>, mat<double,mat_structure::identity> >::type operator *(const mat<double,mat_structure::identity>& M1, const mat<double,mat_structure::identity>& M2);





template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::square,mat_alignment::column_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::square,mat_alignment::row_major>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::symmetric>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::skew_symmetric>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::nil>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::identity>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::scalar>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::diagonal>& M);
template std::ostream& operator <<(std::ostream& out_stream, const mat<float,mat_structure::permutation>& M);

template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);


//rect and square matrix.
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);

//template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::symmetric>& M2);
//template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::symmetric>& M2);

//template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::skew_symmetric>& M2);
//template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::skew_symmetric>& M2);


//square and rect matrix.
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);

//template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *< mat<float,mat_structure::rectangular,mat_alignment::column_major> >(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
//template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *< mat<float,mat_structure::rectangular,mat_alignment::row_major> >(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);

//template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *< mat<float,mat_structure::rectangular,mat_alignment::column_major> >(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
//template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *< mat<float,mat_structure::rectangular,mat_alignment::row_major> >(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);


//square and square matrix.
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);

template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::skew_symmetric>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::symmetric>& M2);

template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::symmetric>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::skew_symmetric>& M2);


//dense and diagonal matrix.
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::diagonal>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::diagonal>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::diagonal>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::diagonal>& M2);
template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::diagonal>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::diagonal>& M2);

//diagonal and dense matrix.
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::symmetric>& M2);
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::skew_symmetric>& M2);

//diagonal and diagonal matrix.
template mat_product_result< mat<float,mat_structure::diagonal>, mat<float,mat_structure::diagonal> >::type operator *(const mat<float,mat_structure::diagonal>& M1, const mat<float,mat_structure::diagonal>& M2);


//dense and nil matrix.
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::nil>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::nil>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::nil>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::nil>& M2);
template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::nil>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::nil>& M2);

//nil and dense matrix.
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::symmetric>& M2);
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::skew_symmetric>& M2);

//nil and nil matrix.
template mat_product_result< mat<float,mat_structure::nil>, mat<float,mat_structure::nil> >::type operator *(const mat<float,mat_structure::nil>& M1, const mat<float,mat_structure::nil>& M2);


//dense and scalar matrix.
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::scalar>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::scalar>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::scalar>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::scalar>& M2);
template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::scalar>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::scalar>& M2);

//scalar and dense matrix.
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::symmetric>& M2);
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::skew_symmetric>& M2);

//scalar and scalar matrix.
template mat_product_result< mat<float,mat_structure::scalar>, mat<float,mat_structure::scalar> >::type operator *(const mat<float,mat_structure::scalar>& M1, const mat<float,mat_structure::scalar>& M2);


//dense and identity matrix.
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::column_major>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M1, const mat<float,mat_structure::identity>& M2);
template mat_product_result< mat<float,mat_structure::rectangular,mat_alignment::row_major>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M1, const mat<float,mat_structure::identity>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::column_major>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::square,mat_alignment::column_major>& M1, const mat<float,mat_structure::identity>& M2);
template mat_product_result< mat<float,mat_structure::square,mat_alignment::row_major>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::square,mat_alignment::row_major>& M1, const mat<float,mat_structure::identity>& M2);
template mat_product_result< mat<float,mat_structure::symmetric>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::symmetric>& M1, const mat<float,mat_structure::identity>& M2);
template mat_product_result< mat<float,mat_structure::skew_symmetric>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::skew_symmetric>& M1, const mat<float,mat_structure::identity>& M2);

//identity and dense matrix.
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::rectangular,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::rectangular,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::square,mat_alignment::column_major> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::square,mat_alignment::column_major>& M2);
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::square,mat_alignment::row_major> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::square,mat_alignment::row_major>& M2);
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::symmetric> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::symmetric>& M2);
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::skew_symmetric> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::skew_symmetric>& M2);

//identity and identity matrix.
template mat_product_result< mat<float,mat_structure::identity>, mat<float,mat_structure::identity> >::type operator *(const mat<float,mat_structure::identity>& M1, const mat<float,mat_structure::identity>& M2);





};

#else

namespace ReaK {

void dummy_mat_alg_externs_symbol() { };

};



#endif














