
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

#include <ReaK/core/lin_alg/mat_num.hpp>

namespace ReaK {



template double norm_1(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M);
template double norm_1(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M);
template double norm_1(const mat<double,mat_structure::square,mat_alignment::column_major>& M);
template double norm_1(const mat<double,mat_structure::square,mat_alignment::row_major>& M);
template double norm_1(const mat<double,mat_structure::symmetric>& M);
template double norm_1(const mat<double,mat_structure::skew_symmetric>& M);
template double norm_1(const mat<double,mat_structure::diagonal>& M);
template double norm_1(const mat<double,mat_structure::scalar>& M);
template double norm_1(const mat<double,mat_structure::identity>& M);
template double norm_1(const mat<double,mat_structure::nil>& M);

template double norm_inf(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M);
template double norm_inf(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M);
template double norm_inf(const mat<double,mat_structure::square,mat_alignment::column_major>& M);
template double norm_inf(const mat<double,mat_structure::square,mat_alignment::row_major>& M);
template double norm_inf(const mat<double,mat_structure::symmetric>& M);
template double norm_inf(const mat<double,mat_structure::skew_symmetric>& M);
template double norm_inf(const mat<double,mat_structure::diagonal>& M);
template double norm_inf(const mat<double,mat_structure::scalar>& M);
template double norm_inf(const mat<double,mat_structure::identity>& M);
template double norm_inf(const mat<double,mat_structure::nil>& M);

template double elem_norm_2(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M);
template double elem_norm_2(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M);
template double elem_norm_2(const mat<double,mat_structure::square,mat_alignment::column_major>& M);
template double elem_norm_2(const mat<double,mat_structure::square,mat_alignment::row_major>& M);
template double elem_norm_2(const mat<double,mat_structure::symmetric>& M);
template double elem_norm_2(const mat<double,mat_structure::skew_symmetric>& M);
template double elem_norm_2(const mat<double,mat_structure::diagonal>& M);
template double elem_norm_2(const mat<double,mat_structure::scalar>& M);
template double elem_norm_2(const mat<double,mat_structure::identity>& M);
template double elem_norm_2(const mat<double,mat_structure::nil>& M);

template double elem_norm_max(const mat<double,mat_structure::rectangular,mat_alignment::column_major>& M);
template double elem_norm_max(const mat<double,mat_structure::rectangular,mat_alignment::row_major>& M);
template double elem_norm_max(const mat<double,mat_structure::square,mat_alignment::column_major>& M);
template double elem_norm_max(const mat<double,mat_structure::square,mat_alignment::row_major>& M);
template double elem_norm_max(const mat<double,mat_structure::symmetric>& M);
template double elem_norm_max(const mat<double,mat_structure::skew_symmetric>& M);
template double elem_norm_max(const mat<double,mat_structure::diagonal>& M);
template double elem_norm_max(const mat<double,mat_structure::scalar>& M);
template double elem_norm_max(const mat<double,mat_structure::identity>& M);
template double elem_norm_max(const mat<double,mat_structure::nil>& M);


template float norm_1(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M);
template float norm_1(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M);
template float norm_1(const mat<float,mat_structure::square,mat_alignment::column_major>& M);
template float norm_1(const mat<float,mat_structure::square,mat_alignment::row_major>& M);
template float norm_1(const mat<float,mat_structure::symmetric>& M);
template float norm_1(const mat<float,mat_structure::skew_symmetric>& M);
template float norm_1(const mat<float,mat_structure::diagonal>& M);
template float norm_1(const mat<float,mat_structure::scalar>& M);
template float norm_1(const mat<float,mat_structure::identity>& M);
template float norm_1(const mat<float,mat_structure::nil>& M);

template float norm_inf(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M);
template float norm_inf(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M);
template float norm_inf(const mat<float,mat_structure::square,mat_alignment::column_major>& M);
template float norm_inf(const mat<float,mat_structure::square,mat_alignment::row_major>& M);
template float norm_inf(const mat<float,mat_structure::symmetric>& M);
template float norm_inf(const mat<float,mat_structure::skew_symmetric>& M);
template float norm_inf(const mat<float,mat_structure::diagonal>& M);
template float norm_inf(const mat<float,mat_structure::scalar>& M);
template float norm_inf(const mat<float,mat_structure::identity>& M);
template float norm_inf(const mat<float,mat_structure::nil>& M);

template float elem_norm_2(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M);
template float elem_norm_2(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M);
template float elem_norm_2(const mat<float,mat_structure::square,mat_alignment::column_major>& M);
template float elem_norm_2(const mat<float,mat_structure::square,mat_alignment::row_major>& M);
template float elem_norm_2(const mat<float,mat_structure::symmetric>& M);
template float elem_norm_2(const mat<float,mat_structure::skew_symmetric>& M);
template float elem_norm_2(const mat<float,mat_structure::diagonal>& M);
template float elem_norm_2(const mat<float,mat_structure::scalar>& M);
template float elem_norm_2(const mat<float,mat_structure::identity>& M);
template float elem_norm_2(const mat<float,mat_structure::nil>& M);

template float elem_norm_max(const mat<float,mat_structure::rectangular,mat_alignment::column_major>& M);
template float elem_norm_max(const mat<float,mat_structure::rectangular,mat_alignment::row_major>& M);
template float elem_norm_max(const mat<float,mat_structure::square,mat_alignment::column_major>& M);
template float elem_norm_max(const mat<float,mat_structure::square,mat_alignment::row_major>& M);
template float elem_norm_max(const mat<float,mat_structure::symmetric>& M);
template float elem_norm_max(const mat<float,mat_structure::skew_symmetric>& M);
template float elem_norm_max(const mat<float,mat_structure::diagonal>& M);
template float elem_norm_max(const mat<float,mat_structure::scalar>& M);
template float elem_norm_max(const mat<float,mat_structure::identity>& M);
template float elem_norm_max(const mat<float,mat_structure::nil>& M);





template void invert_gaussian(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
template void invert_gaussian(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_inv, double NumTol);
template void invert_gaussian(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
template void invert_gaussian(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_inv, double NumTol);

template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& b, vect_n<unsigned int>& P, double NumTol);
template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& b, vect_n<unsigned int>& P, double NumTol);
template void linsolve_PLU(mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& b, vect_n<unsigned int>& P, double NumTol);
template void linsolve_PLU(mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& b, vect_n<unsigned int>& P, double NumTol);

template void linsolve_PLU(mat<double,mat_structure::rectangular>& A, vect_n<double>& b, vect_n<unsigned int>& P, double NumTol);
template void linsolve_PLU(mat<double,mat_structure::square>& A, vect_n<double>& b, vect_n<unsigned int>& P, double NumTol);

template void invert_PLU(mat<double,mat_structure::rectangular> A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
template void invert_PLU(mat<double,mat_structure::rectangular> A, mat<double,mat_structure::square>& A_inv, double NumTol);
template void invert_PLU(mat<double,mat_structure::square> A, mat<double,mat_structure::rectangular>& A_inv, double NumTol);
template void invert_PLU(mat<double,mat_structure::square> A, mat<double,mat_structure::square>& A_inv, double NumTol);


template void invert_gaussian(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
template void invert_gaussian(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_inv, float NumTol);
template void invert_gaussian(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
template void invert_gaussian(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_inv, float NumTol);

template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& b, vect_n<unsigned int>& P, float NumTol);
template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& b, vect_n<unsigned int>& P, float NumTol);
template void linsolve_PLU(mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& b, vect_n<unsigned int>& P, float NumTol);
template void linsolve_PLU(mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& b, vect_n<unsigned int>& P, float NumTol);

template void linsolve_PLU(mat<float,mat_structure::rectangular>& A, vect_n<float>& b, vect_n<unsigned int>& P, float NumTol);
template void linsolve_PLU(mat<float,mat_structure::square>& A, vect_n<float>& b, vect_n<unsigned int>& P, float NumTol);

template void invert_PLU(mat<float,mat_structure::rectangular> A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
template void invert_PLU(mat<float,mat_structure::rectangular> A, mat<float,mat_structure::square>& A_inv, float NumTol);
template void invert_PLU(mat<float,mat_structure::square> A, mat<float,mat_structure::rectangular>& A_inv, float NumTol);
template void invert_PLU(mat<float,mat_structure::square> A, mat<float,mat_structure::square>& A_inv, float NumTol);






template void decompose_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& L, double NumTol);
template void decompose_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& L, double NumTol);

template void linsolve_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& b, double NumTol);
template void linsolve_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& b, double NumTol);

template void linsolve_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::symmetric>& b, double NumTol);
template void linsolve_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& b, double NumTol);

template void invert_Cholesky(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_inv, double NumTol);
template void invert_Cholesky(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& A_inv, double NumTol);


template void decompose_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& L, float NumTol);
template void decompose_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& L, float NumTol);

template void linsolve_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& b, float NumTol);
template void linsolve_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& b, float NumTol);

template void linsolve_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::symmetric>& b, float NumTol);
template void linsolve_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& b, float NumTol);

template void invert_Cholesky(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_inv, float NumTol);
template void invert_Cholesky(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& A_inv, float NumTol);







template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::rectangular>& Q, double NumTol);
template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& Q, double NumTol);

template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::rectangular>& E, mat<double,mat_structure::rectangular>& Q, double NumTol);
template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::rectangular>& E, mat<double,mat_structure::square>& Q, double NumTol);
template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& E, mat<double,mat_structure::rectangular>& Q, double NumTol);
template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& E, mat<double,mat_structure::square>& Q, double NumTol);
template void eigensolve_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& E, mat<double,mat_structure::square>& Q, double NumTol);

template void linlsq_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void linlsq_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

template void pseudoinvert_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::square>& A_inv, double NumTol);
template void pseudoinvert_Jacobi(const mat<double,mat_structure::symmetric>& A, mat<double,mat_structure::symmetric>& A_inv, double NumTol);


template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::rectangular>& Q, float NumTol);
template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& Q, float NumTol);

template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::rectangular>& E, mat<float,mat_structure::rectangular>& Q, float NumTol);
template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::rectangular>& E, mat<float,mat_structure::square>& Q, float NumTol);
template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& E, mat<float,mat_structure::rectangular>& Q, float NumTol);
template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& E, mat<float,mat_structure::square>& Q, float NumTol);
template void eigensolve_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& E, mat<float,mat_structure::square>& Q, float NumTol);

template void linlsq_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void linlsq_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

template void pseudoinvert_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::square>& A_inv, float NumTol);
template void pseudoinvert_Jacobi(const mat<float,mat_structure::symmetric>& A, mat<float,mat_structure::symmetric>& A_inv, float NumTol);






template void decompose_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& Q, mat<double,mat_structure::rectangular>& R, double NumTol);
template void decompose_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::rectangular>& R, double NumTol);
template void decompose_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& R, double NumTol);

template double determinant_QR(const mat<double,mat_structure::rectangular>& A, double NumTol);
template double determinant_QR(const mat<double,mat_structure::square>& A, double NumTol);

template void linlsq_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void linlsq_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void linlsq_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

template void minnorm_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void minnorm_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void minnorm_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

template void pseudoinvert_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_QR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
template void pseudoinvert_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_QR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);

template void linlsq_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void linlsq_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void linlsq_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

template void minnorm_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void minnorm_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& x, const mat<double,mat_structure::rectangular>& b, double NumTol);
template void minnorm_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& x, const mat<double,mat_structure::square>& b, double NumTol);

template void pseudoinvert_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_RRQR(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
template void pseudoinvert_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_RRQR(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);


template void decompose_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& Q, mat<float,mat_structure::rectangular>& R, float NumTol);
template void decompose_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::rectangular>& R, float NumTol);
template void decompose_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& R, float NumTol);

template float determinant_QR(const mat<float,mat_structure::rectangular>& A, float NumTol);
template float determinant_QR(const mat<float,mat_structure::square>& A, float NumTol);

template void linlsq_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void linlsq_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void linlsq_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

template void minnorm_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void minnorm_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void minnorm_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

template void pseudoinvert_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_QR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
template void pseudoinvert_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_QR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);

template void linlsq_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void linlsq_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void linlsq_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

template void minnorm_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void minnorm_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& x, const mat<float,mat_structure::rectangular>& b, float NumTol);
template void minnorm_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& x, const mat<float,mat_structure::square>& b, float NumTol);

template void pseudoinvert_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_RRQR(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
template void pseudoinvert_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_RRQR(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);







template void decompose_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::rectangular>& V, double NumTol);
template void decompose_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);
template void decompose_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);
template void decompose_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& U, mat<double,mat_structure::diagonal>& E, mat<double,mat_structure::square>& V, double NumTol);

template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::rectangular>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::square>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::square>& U, const mat<double,mat_structure::diagonal>& E, const mat<double,mat_structure::square>& V, mat<double,mat_structure::square>& A_pinv, double NumTol);

template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::rectangular>& A_pinv, double NumTol);
template void pseudoinvert_SVD(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& A_pinv, double NumTol);


template void decompose_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::rectangular>& V, float NumTol);
template void decompose_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);
template void decompose_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);
template void decompose_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& U, mat<float,mat_structure::diagonal>& E, mat<float,mat_structure::square>& V, float NumTol);

template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::rectangular>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::square>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::square>& U, const mat<float,mat_structure::diagonal>& E, const mat<float,mat_structure::square>& V, mat<float,mat_structure::square>& A_pinv, float NumTol);

template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::rectangular>& A_pinv, float NumTol);
template void pseudoinvert_SVD(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& A_pinv, float NumTol);
 





template void balance(mat<double,mat_structure::rectangular>& A, vect_n<int>& D);
template void balance(mat<double,mat_structure::square>& A, vect_n<int>& D);

template void balance_pencil(mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::diagonal>& Dl, mat<double,mat_structure::diagonal>& Dr);
template void balance_pencil(mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& B, mat<double,mat_structure::diagonal>& Dl, mat<double,mat_structure::diagonal>& Dr);
 

template void balance(mat<float,mat_structure::rectangular>& A, vect_n<int>& D);
template void balance(mat<float,mat_structure::square>& A, vect_n<int>& D);

template void balance_pencil(mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::diagonal>& Dl, mat<float,mat_structure::diagonal>& Dr);
template void balance_pencil(mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& B, mat<float,mat_structure::diagonal>& Dl, mat<float,mat_structure::diagonal>& Dr);






template void exp_PadeSAS(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& X, QR_linlsqsolver linsolve, double NumTol);
template void exp_PadeSAS(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& X, RRQR_linlsqsolver linsolve, double NumTol);


template void exp_PadeSAS(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& X, QR_linlsqsolver linsolve, float NumTol);
template void exp_PadeSAS(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& X, RRQR_linlsqsolver linsolve, float NumTol);






template void decompose_Hess(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& H, double NumTol);
template void decompose_Hess(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& Q, mat<double,mat_structure::rectangular>& H, double NumTol);
template void decompose_Hess(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::rectangular>& H, double NumTol);

template void decompose_Hess(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& H, double NumTol);
template void decompose_Hess(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& H, double NumTol);
template void decompose_Hess(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& H, double NumTol);

template void reduce_HessTri(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::square>& B, mat<double,mat_structure::square>& H, mat<double,mat_structure::square>& R, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& Z, double NumTol);
template void reduce_HessTri(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::rectangular>& H, mat<double,mat_structure::rectangular>& R, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& Z, double NumTol);
template void reduce_HessTri(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::rectangular>& H, mat<double,mat_structure::rectangular>& R, mat<double,mat_structure::rectangular>& Q, mat<double,mat_structure::rectangular>& Z, double NumTol);

template void reduce_HessTri(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::square>& B, mat<double,mat_structure::square>& H, mat<double,mat_structure::square>& R, double NumTol);
template void reduce_HessTri(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::rectangular>& H, mat<double,mat_structure::rectangular>& R, double NumTol);


template void decompose_Hess(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& H, float NumTol);
template void decompose_Hess(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& Q, mat<float,mat_structure::rectangular>& H, float NumTol);
template void decompose_Hess(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::rectangular>& H, float NumTol);

template void decompose_Hess(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& H, float NumTol);
template void decompose_Hess(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& H, float NumTol);
template void decompose_Hess(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& H, float NumTol);

template void reduce_HessTri(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::square>& B, mat<float,mat_structure::square>& H, mat<float,mat_structure::square>& R, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& Z, float NumTol);
template void reduce_HessTri(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::rectangular>& H, mat<float,mat_structure::rectangular>& R, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& Z, float NumTol);
template void reduce_HessTri(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::rectangular>& H, mat<float,mat_structure::rectangular>& R, mat<float,mat_structure::rectangular>& Q, mat<float,mat_structure::rectangular>& Z, float NumTol);

template void reduce_HessTri(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::square>& B, mat<float,mat_structure::square>& H, mat<float,mat_structure::square>& R, float NumTol);
template void reduce_HessTri(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::rectangular>& H, mat<float,mat_structure::rectangular>& R, float NumTol);






template void decompose_RealSchur(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& T, double NumTol);
template void decompose_RealSchur(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::square>& Q, mat<double,mat_structure::rectangular>& T, double NumTol);

template void decompose_RealSchur(const mat<double,mat_structure::square>& A, mat<double,mat_structure::square>& T, double NumTol);
template void decompose_RealSchur(const mat<double,mat_structure::rectangular>& A, mat<double,mat_structure::rectangular>& T, double NumTol);

template void decompose_GenRealSchur(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::square>& B, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& Z, mat<double,mat_structure::square>& T, mat<double,mat_structure::square>& R, double NumTol);
template void decompose_GenRealSchur(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::square>& Q, mat<double,mat_structure::square>& Z, mat<double,mat_structure::rectangular>& T, mat<double,mat_structure::rectangular>& R, double NumTol);
template void decompose_GenRealSchur(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::rectangular>& Q, mat<double,mat_structure::rectangular>& Z, mat<double,mat_structure::rectangular>& T, mat<double,mat_structure::rectangular>& R, double NumTol);

template void decompose_GenRealSchur(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::square>& B, mat<double,mat_structure::square>& T, mat<double,mat_structure::square>& R, double NumTol);
template void decompose_GenRealSchur(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, mat<double,mat_structure::rectangular>& T, mat<double,mat_structure::rectangular>& R, double NumTol);


template void decompose_RealSchur(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& T, float NumTol);
template void decompose_RealSchur(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::square>& Q, mat<float,mat_structure::rectangular>& T, float NumTol);

template void decompose_RealSchur(const mat<float,mat_structure::square>& A, mat<float,mat_structure::square>& T, float NumTol);
template void decompose_RealSchur(const mat<float,mat_structure::rectangular>& A, mat<float,mat_structure::rectangular>& T, float NumTol);

template void decompose_GenRealSchur(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::square>& B, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& Z, mat<float,mat_structure::square>& T, mat<float,mat_structure::square>& R, float NumTol);
template void decompose_GenRealSchur(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::square>& Q, mat<float,mat_structure::square>& Z, mat<float,mat_structure::rectangular>& T, mat<float,mat_structure::rectangular>& R, float NumTol);
template void decompose_GenRealSchur(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::rectangular>& Q, mat<float,mat_structure::rectangular>& Z, mat<float,mat_structure::rectangular>& T, mat<float,mat_structure::rectangular>& R, float NumTol);

template void decompose_GenRealSchur(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::square>& B, mat<float,mat_structure::square>& T, mat<float,mat_structure::square>& R, float NumTol);
template void decompose_GenRealSchur(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, mat<float,mat_structure::rectangular>& T, mat<float,mat_structure::rectangular>& R, float NumTol);






template void solve_care_problem(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::square>& Q, const mat<double,mat_structure::square>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_care_problem(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::symmetric>& Q, const mat<double,mat_structure::symmetric>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_care_problem(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::diagonal>& Q, const mat<double,mat_structure::diagonal>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_care_problem(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::square>& Q, const mat<double,mat_structure::square>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);
template void solve_care_problem(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::symmetric>& Q, const mat<double,mat_structure::symmetric>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);
template void solve_care_problem(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::diagonal>& Q, const mat<double,mat_structure::diagonal>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);

template void solve_dare_problem(const mat<double,mat_structure::square>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::square>& Q, const mat<double,mat_structure::square>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<double,mat_structure::square>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::symmetric>& Q, const mat<double,mat_structure::symmetric>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<double,mat_structure::square>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::diagonal>& Q, const mat<double,mat_structure::diagonal>& R, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<double,mat_structure::rectangular>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::square>& Q, const mat<double,mat_structure::square>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<double,mat_structure::rectangular>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::symmetric>& Q, const mat<double,mat_structure::symmetric>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<double,mat_structure::rectangular>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::diagonal>& Q, const mat<double,mat_structure::diagonal>& R, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);

template void solve_ctsf_problem(const mat<double,mat_structure::square>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::rectangular>& C, const mat<double,mat_structure::rectangular>& D, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_ctsf_problem(const mat<double,mat_structure::rectangular>& A, const mat<double,mat_structure::rectangular>& B, 
                                        const mat<double,mat_structure::rectangular>& C, const mat<double,mat_structure::rectangular>& D, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);

template void solve_dtsf_problem(const mat<double,mat_structure::square>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::rectangular>& H, const mat<double,mat_structure::rectangular>& J, 
                                        mat<double,mat_structure::square>& P, double NumTol, bool UseBalancing);
template void solve_dtsf_problem(const mat<double,mat_structure::rectangular>& F, const mat<double,mat_structure::rectangular>& G, 
                                        const mat<double,mat_structure::rectangular>& H, const mat<double,mat_structure::rectangular>& J, 
                                        mat<double,mat_structure::rectangular>& P, double NumTol, bool UseBalancing);


template void solve_care_problem(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::square>& Q, const mat<float,mat_structure::square>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_care_problem(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::symmetric>& Q, const mat<float,mat_structure::symmetric>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_care_problem(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::diagonal>& Q, const mat<float,mat_structure::diagonal>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_care_problem(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::square>& Q, const mat<float,mat_structure::square>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);
template void solve_care_problem(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::symmetric>& Q, const mat<float,mat_structure::symmetric>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);
template void solve_care_problem(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::diagonal>& Q, const mat<float,mat_structure::diagonal>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);

template void solve_dare_problem(const mat<float,mat_structure::square>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::square>& Q, const mat<float,mat_structure::square>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<float,mat_structure::square>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::symmetric>& Q, const mat<float,mat_structure::symmetric>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<float,mat_structure::square>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::diagonal>& Q, const mat<float,mat_structure::diagonal>& R, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<float,mat_structure::rectangular>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::square>& Q, const mat<float,mat_structure::square>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<float,mat_structure::rectangular>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::symmetric>& Q, const mat<float,mat_structure::symmetric>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);
template void solve_dare_problem(const mat<float,mat_structure::rectangular>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::diagonal>& Q, const mat<float,mat_structure::diagonal>& R, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);

template void solve_ctsf_problem(const mat<float,mat_structure::square>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::rectangular>& C, const mat<float,mat_structure::rectangular>& D, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_ctsf_problem(const mat<float,mat_structure::rectangular>& A, const mat<float,mat_structure::rectangular>& B, 
                                        const mat<float,mat_structure::rectangular>& C, const mat<float,mat_structure::rectangular>& D, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);

template void solve_dtsf_problem(const mat<float,mat_structure::square>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::rectangular>& H, const mat<float,mat_structure::rectangular>& J, 
                                        mat<float,mat_structure::square>& P, float NumTol, bool UseBalancing);
template void solve_dtsf_problem(const mat<float,mat_structure::rectangular>& F, const mat<float,mat_structure::rectangular>& G, 
                                        const mat<float,mat_structure::rectangular>& H, const mat<float,mat_structure::rectangular>& J, 
                                        mat<float,mat_structure::rectangular>& P, float NumTol, bool UseBalancing);



};

#else

namespace ReaK {

void dummy_mat_num_externs_symbol() { };

};



#endif














