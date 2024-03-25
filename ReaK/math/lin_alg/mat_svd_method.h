/**
 * \file mat_svd_method.h
 *
 * This library provides a number of functions related to performing a Singular Value Decomposition (SVD)
 * on a matrix, e.g., to invert a matrix, to pseudo-invert a matrix, to solve a linear system with
 * least-square error and to find the eigen-values of a symmetric matrix. SVD is not very efficient but
 * very powerful. The SVD implementation used here is based on the implementation from CLARAty,
 * developed by the Jet Propulsion Laboratory.
 *
 * According to performance tests, PLU methods are as good as Cholesky methods in terms of speed.
 * And they are both the best for well-conditioned matrices. For ill-conditioned matrices, QR-decomposition
 * methods are only a little slower then PLU (about 20% slower, same time-complexity) but provide better
 * numerical stability. The Jacobi methods are significantly slower, but this implementation is in need
 * of a revision for performance enhancement. And, of course, SVD is also very slow (slightly faster than
 * Jacobi) but it is based on a LAPACK implementation that is very poorly written, and it has not been
 * updated since.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_MATH_LIN_ALG_MAT_SVD_METHOD_H_
#define REAK_MATH_LIN_ALG_MAT_SVD_METHOD_H_

#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/mat_concepts.h"
#include "ReaK/math/lin_alg/mat_num_exceptions.h"

#include <type_traits>

namespace ReaK {

/*************************************************************************
                    Singular Value Decomposition (SVD)
*************************************************************************/

/**
 * Singular Value Decomposition.
 *
 * For an N-by-M matrix A with N >= M, the singular value decomposition is
 * an N-by-M orthogonal matrix U, an M-by-M diagonal matrix S, and
 * an M-by-M orthogonal matrix V so that A = U*S*V'.
 *
 * The singular values, sigma(k) = S(k, k), are ordered so that
 * sigma(0) >= sigma(1) >= ... >= sigma(M-1).
 *
 * The singular value decompostion always exists, so the constructor will
 * never fail.  The matrix condition number and the effective numerical
 * rank can be computed from this decomposition.
 *
 * NIST Disclaimer:
 * This software was developed at the National Institute of Standards and
 * Technology (NIST) by employees of the Federal Government in the course
 * of their official duties. Pursuant to title 17 Section 105 of the
 * United States Code this software is not subject to copyright
 * protection and is in the public domain. NIST assumes no responsibility
 * whatsoever for its use by other parties, and makes no guarantees,
 * expressed or implied, about its quality, reliability, or any other
 * characteristic.
 *
 * &copy; 2006, Jet Propulsion Laboratory, California Institute of Technology<br>
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param U output unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order (size = min(RowCount,ColCount))
 * \param V output unitary matrix (ColCount x ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author NIST
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2,
          WritableMatrix Matrix3, WritableMatrix Matrix4>
void decompose_SVD(const Matrix1& A, Matrix2& U, Matrix3& E, Matrix4& V,
                   mat_value_type_t<Matrix1> NumTol = 1E-15) requires
    DiagonalMatrix<Matrix3> {
  using ValueType = mat_value_type_t<Matrix1>;
  using std::abs;
  using std::sqrt;
  using std::swap;

  mat<ValueType, mat_structure::rectangular> At;
  mat<ValueType, mat_structure::rectangular> Ut;
  mat<ValueType, mat_structure::rectangular> Vt;
  if (A.get_row_count() < A.get_col_count()) {
    At = transpose_view(A);
  } else {
    At = A;
  }

  int N = At.get_row_count();
  int M = At.get_col_count();
  int nu = M;
  if (nu > N) {
    nu = N;
  }
  if ((N == 0) || (M == 0)) {
    throw std::range_error("SVD-Decomp: Matrix A has 0 rows or columns!");
  }
  int nct = 0;
  int nrt = 0;
  // int nct = min(N-1,M);
  // int nrt = max(0,min(int(M-2),int(N)));
  if (N >= M + 1) {
    nct = M;
  } else {
    nct = N - 1;
  }
  if (M >= N + 2) {
    nrt = N;
  } else {
    if (M >= 2) {
      nrt = M - 2;
    } else {
      nrt = 0;
    }
  }
  int max_iter = nct;
  if (max_iter < nrt) {
    max_iter = nrt;
  }

  if ((Ut.get_row_count() != N) || (Ut.get_col_count() != nu)) {
    Ut.set_row_count(N);
    Ut.set_col_count(nu);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        if (i == j) {
          Ut(i, i) = 1;
        } else {
          Ut(i, j) = 0;
        }
      }
    }
  }
  if ((Vt.get_row_count() != M) || (Vt.get_col_count() != nu)) {
    Vt.set_row_count(M);
    Vt.set_col_count(nu);
    for (int i = 0; i < M; ++i) {
      for (int j = 0; j < nu; ++j) {
        if (i == j) {
          Vt(i, i) = 1;
        } else {
          Vt(i, j) = 0;
        }
      }
    }
  }
  if ((E.get_row_count() != nu) || (E.get_col_count() != nu)) {
    E = mat<ValueType, mat_structure::identity>(nu);
  }
  vect_n<ValueType> e(M);
  vect_n<ValueType> work(N);

  // Reduce A to bidiagonal form, storing the diagonal elements
  // in E and the super-diagonal elements in e.
  for (int k = 0; k < max_iter; ++k) {

    if (k < nct) {

      // Compute the transformation for the k-th column and
      // place the k-th diagonal in E(k).
      // Compute 2-norm of k-th column without under/overflow.
      E(k, k) = 0;
      for (int i = k; i < N; ++i) {
        E(k, k) = sqrt(E(k, k) * E(k, k) + At(i, k) * At(i, k));
      }

      if (abs(E(k, k)) > NumTol) {
        if (At(k, k) < 0.0) {
          E(k, k) = -E(k, k);
        }
        for (int i = k; i < N; ++i) {
          At(i, k) /= E(k, k);
        }
        At(k, k) += 1.0;
      }
      E(k, k) = -E(k, k);
    }
    for (int j = k + 1; j < M; ++j) {
      if ((k < nct) && (abs(E(k, k)) > NumTol)) {

        // Apply the transformation.
        ValueType t(0.0);
        for (int i = k; i < N; ++i) {
          t += At(i, k) * At(i, j);
        }

        t = -t / At(k, k);

        for (int i = k; i < N; ++i) {
          At(i, j) += t * At(i, k);
        }
      }

      // Place the k-th row of A into e for the
      // subsequent calculation of the row transformation.
      e[j] = At(k, j);
    }
    if (k < nct) {
      // Place the transformation in U for subsequent back
      // multiplication.
      for (int i = k; i < N; ++i) {
        Ut(i, k) = At(i, k);
      }
    }
    if (k < nrt) {

      // Compute the k-th row transformation and place the
      // k-th super-diagonal in e(k).
      // Compute 2-norm without under/overflow.
      e[k] = 0;
      for (int i = k + 1; i < M; ++i) {
        e[k] = sqrt(e[k] * e[k] + e[i] * e[i]);
      }
      if (abs(e[k]) > NumTol) {
        if (e[k + 1] < 0.0) {
          e[k] = -e[k];
        }
        for (int i = k + 1; i < M; ++i) {
          e[i] /= e[k];
        }
        e[k + 1] += 1.0;
      }
      e[k] = -e[k];
      if ((k + 1 < N) & (abs(e[k]) > NumTol)) {

        // Apply the transformation.

        for (int i = k + 1; i < N; ++i) {
          work[i] = 0.0;
        }
        for (int j = k + 1; j < M; ++j) {
          for (int i = k + 1; i < N; ++i) {
            work[i] += e[j] * At(i, j);
          }
        }
        for (int j = k + 1; j < M; ++j) {
          ValueType t = -e[j] / e[k + 1];
          for (int i = k + 1; i < N; ++i) {
            At(i, j) += t * work[i];
          }
        }
      }

      // Place the transformation in _V for subsequent
      // back multiplication.

      for (int i = k + 1; i < M; ++i) {
        Vt(i, k) = e[i];
      }
    }
  }

  // Set up the final bidiagonal matrix or order p.

  int p = M;
  if (p > N + 1) {
    p = N + 1;
  }
  if (nct < M) {
    E(nct, nct) = At(nct, nct);
  }
  if (N < p) {
    E(p - 1, p - 1) = 0.0;
  }
  if (nrt + 1 < p) {
    e[nrt] = At(nrt, p - 1);
  }

  e[p - 1] = 0.0;

  // Generate U.

  for (int j = nct; j < nu; ++j) {
    for (int i = 0; i < N; ++i) {
      Ut(i, j) = 0.0;
    }
    Ut(j, j) = 1.0;
  }

  for (int k = static_cast<int>(nct - 1); k >= 0; --k) {
    if (abs(E(k, k)) > NumTol) {
      for (int j = k + 1; j < nu; ++j) {
        ValueType t = 0;
        for (int i = k; i < N; ++i) {
          t += Ut(i, k) * Ut(i, j);
        }
        t = -t / Ut(k, k);
        for (int i = k; i < N; ++i) {
          Ut(i, j) += t * Ut(i, k);
        }
      }
      for (int i = k; i < N; ++i) {
        Ut(i, k) = -Ut(i, k);
      }
      Ut(k, k) = 1.0 + Ut(k, k);
      for (int i = 0; int(i) < k - 1; ++i) {
        Ut(i, k) = 0.0;
      }
    } else {
      for (int i = 0; i < N; ++i) {
        Ut(i, k) = 0.0;
      }
      Ut(k, k) = 1.0;
    }
  }

  // Generate V.

  for (int k = static_cast<int>(nu - 1); k >= 0; --k) {
    if ((int(k) < nrt) & (abs(e[k]) > NumTol)) {
      for (int j = k + 1; j < nu; ++j) {
        ValueType t = 0;
        for (int i = k + 1; i < M; ++i) {
          t += Vt(i, k) * Vt(i, j);
        }
        t = -t / Vt(k + 1, k);
        for (int i = k + 1; i < M; ++i) {
          Vt(i, j) += t * Vt(i, k);
        }
      }
    }
    for (int i = 0; i < M; ++i) {
      Vt(i, k) = 0.0;
    }
    Vt(k, k) = 1.0;
  }

  // Main iteration loop for the singular values.
  int pp = static_cast<int>(p - 1);
  int iter = 0;
  while (p > 0) {
    int k = 0;
    int kase = 0;

    // Here is where a test for too many iterations would go.

    // This section of the program inspects for
    // negligible elements in the E and e arrays.  On
    // completion the variables kase and k are set as follows.

    // kase = 1     if E(p) and e(k-1) are negligible and k<p
    // kase = 2     if E(k) is negligible and k<p
    // kase = 3     if e(k-1) is negligible, k<p, and
    //              E(k), ..., E(p)
    //              are not negligible (qr step).
    // kase = 4     if e(p-1) is negligible (convergence).

    for (k = static_cast<int>(p - 2); k >= -1; --k) {
      if (k == -1) {
        break;
      }
      if (abs(e[k]) <= NumTol * (abs(E(k, k)) + abs(E(k + 1, k + 1)))) {
        e[k] = 0.0;
        break;
      }
    }
    if (k == int(p - 2)) {
      kase = 4;
    } else {
      int ks = 0;
      for (ks = static_cast<int>(p - 1); ks >= k; --ks) {
        if (ks == k) {
          break;
        }
        ValueType t = (ks != int(p) ? abs(e[ks]) : 0.0) +
                      (ks != k + 1 ? abs(e[ks - 1]) : 0.0);
        if (abs(E(ks, ks)) <= NumTol * t) {
          E(ks, ks) = 0.0;
          break;
        }
      }
      if (ks == k) {
        kase = 3;
      } else if (ks == int(p - 1)) {
        kase = 1;
      } else {
        kase = 2;
        k = ks;
      }
    }
    k++;

    // Perform the task indicated by kase.

    switch (kase) {

      // Deflate negligible E(p).
      case 1: {
        ValueType f = e[p - 2];
        e[p - 2] = 0.0;
        for (int j = static_cast<int>(p - 2); j >= k; --j) {
          ValueType t = sqrt(E(j, j) * E(j, j) + f * f);
          ValueType cs = E(j, j) / t;
          ValueType sn = f / t;
          E(j, j) = t;
          if (j != k) {
            f = -sn * e[j - 1];
            e[j - 1] = cs * e[j - 1];
          }
          for (int i = 0; i < M; ++i) {
            t = cs * Vt(i, j) + sn * Vt(i, p - 1);
            Vt(i, p - 1) = -sn * Vt(i, j) + cs * Vt(i, p - 1);
            Vt(i, j) = t;
          }
        }
      } break;

        // Split at negligible E(k).

      case 2: {
        ValueType f = e[k - 1];
        e[k - 1] = 0.0;
        for (int j = k; j < p; ++j) {
          ValueType t = sqrt(E(j, j) * E(j, j) + f * f);
          ValueType cs = E(j, j) / t;
          ValueType sn = f / t;
          E(j, j) = t;
          f = -sn * e[j];
          e[j] = cs * e[j];

          for (int i = 0; i < N; ++i) {
            t = cs * Ut(i, j) + sn * Ut(i, k - 1);
            Ut(i, k - 1) = -sn * Ut(i, j) + cs * Ut(i, k - 1);
            Ut(i, j) = t;
          }
        }
      } break;

        // Perform one qr step.

      case 3: {

        // Calculate the shift.

        // ValueType scale = MAX(MAX(MAX(MAX(abs(E_a[p-1]),abs(E_a[p-2])),abs(e[p-2])),abs(E_a[k])),abs(e[k]));
        ValueType scale = abs(E(p - 1, p - 1));
        if (scale < abs(E(p - 2, p - 2))) {
          scale = abs(E(p - 2, p - 2));
        }
        if (scale < abs(E(k, k))) {
          scale = abs(E(k, k));
        }
        if (scale < abs(e[p - 2])) {
          scale = abs(e[p - 2]);
        }
        if (scale < abs(e[k])) {
          scale = abs(e[k]);
        }
        ValueType sp = E(p - 1, p - 1) / scale;
        ValueType spm1 = E(p - 2, p - 2) / scale;
        ValueType epm1 = e[p - 2] / scale;
        ValueType sk = E(k, k) / scale;
        ValueType ek = e[k] / scale;
        ValueType b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
        ValueType c = (sp * epm1) * (sp * epm1);
        ValueType shift = 0.0;
        if ((b != 0.0) | (c != 0.0)) {
          shift = sqrt(b * b + c);
          if (b < 0.0) {
            shift = -shift;
          }
          shift = c / (b + shift);
        }
        ValueType f = (sk + sp) * (sk - sp) + shift;
        ValueType g = sk * ek;

        // Chase zeros.

        for (int j = k; j < p - 1; ++j) {
          ValueType t = sqrt(f * f + g * g);
          ValueType cs = f / t;
          ValueType sn = g / t;
          if (j != int(k)) {
            e[j - 1] = t;
          }
          f = cs * E(j, j) + sn * e[j];
          e[j] = cs * e[j] - sn * E(j, j);
          g = sn * E(j + 1, j + 1);
          E(j + 1, j + 1) = cs * E(j + 1, j + 1);

          for (int i = 0; i < M; ++i) {
            t = cs * Vt(i, j) + sn * Vt(i, j + 1);
            Vt(i, j + 1) = -sn * Vt(i, j) + cs * Vt(i, j + 1);
            Vt(i, j) = t;
          }
          t = sqrt(f * f + g * g);
          cs = f / t;
          sn = g / t;
          E(j, j) = t;
          f = cs * e[j] + sn * E(j + 1, j + 1);
          E(j + 1, j + 1) = -sn * e[j] + cs * E(j + 1, j + 1);
          g = sn * e[j + 1];
          e[j + 1] = cs * e[j + 1];
          if (j < N - 1) {
            for (int i = 0; i < N; ++i) {
              t = cs * Ut(i, j) + sn * Ut(i, j + 1);
              Ut(i, j + 1) = -sn * Ut(i, j) + cs * Ut(i, j + 1);
              Ut(i, j) = t;
            }
          }
        }
        e[p - 2] = f;
        iter = iter + 1;
      } break;

        // Convergence.

      case 4: {

        // Make the singular values positive.

        if (E(k, k) <= 0.0) {
          E(k, k) = (E(k, k) < 0.0 ? -E(k, k) : 0.0);
          for (int i = 0; i <= int(pp); ++i) {
            Vt(i, k) = -Vt(i, k);
          }
        }

        // Order the singular values.

        while (k < pp) {
          if (E(k, k) >= E(k + 1, k + 1)) {
            break;
          }
          swap(E(k, k), E(k + 1, k + 1));
          if (k < int(M - 1)) {
            for (int i = 0; i < M; ++i) {
              swap(Vt(i, k), Vt(i, k + 1));
            }
          }
          if (k < int(N - 1)) {
            for (int i = 0; i < N; ++i) {
              swap(Ut(i, k), Ut(i, k + 1));
            }
          }
          k++;
        }
        iter = 0;
        p--;
      } break;
    }
  }

  if (A.get_row_count() < A.get_col_count()) {
    U = Vt;
    V = Ut;
  } else {
    U = Ut;
    V = Vt;
  }
}

/**
 * This function returns the two norm of a singular value decomposition, i.e. the
 * highest singular value.
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \return the two-norm.
 *
 * \throws std::range_error if the vector E is empty.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> two_norm_SVD(const Matrix& E) {
  if (E.get_row_count() == 0) {
    throw std::range_error(
        "No singular values available for 2-norm evaluation!");
  }
  return E(0, 0);
}

/**
 * This function returns the two norm (induced norm) of a matrix via the computation of
 * its singular value decomposition.
 *
 * \param M rectangular matrix (RowCount x ColCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the two-norm.
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> norm_2(const Matrix& A,
                                mat_value_type_t<Matrix> NumTol = 1E-15) {
  using ValueType = mat_value_type_t<Matrix>;

  std::size_t nu = (A.get_row_count() > A.get_col_count() ? A.get_col_count()
                                                          : A.get_row_count());
  mat<ValueType, mat_structure::rectangular> U(A.get_row_count(), nu);
  mat<ValueType, mat_structure::diagonal> E(nu);
  mat<ValueType, mat_structure::rectangular> V(A.get_col_count(), nu);

  decompose_SVD(A, U, E, V, NumTol);

  return E(0, 0);
}

/**
 * This function returns the condition number of a singular value decomposition,
 * i.e. max(E)/min(E).
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \return the condition number.
 *
 * \throws std::range_error if the vector E is empty.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix>
mat_value_type_t<Matrix> condition_number_SVD(const Matrix& E) {
  if (E.get_row_count() == 0) {
    throw std::range_error(
        "No singular values available for condition number evaluation!");
  }
  return E(0, 0) / E(E.get_row_count() - 1, E.get_row_count() - 1);
}

/**
 * This function computes the effective numerical rank of a singular value
 * decomposition.
 *
 * \param E vector of singular values, sorted in decreasing order.
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 * \return the numerical rank.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix>
int numrank_SVD(const Matrix& E, mat_value_type_t<Matrix> NumTol = 1E-8) {
  using std::abs;
  int r = 0;
  for (mat_size_type_t<Matrix> i = 0; i < E.get_row_count(); ++i) {
    if (abs(E(i, i)) > abs(E(0, 0)) * NumTol) {
      ++r;
    }
  }
  return r;
}

/**
 * Computes the pseudo-inverse of a RowCount x ColCount matrix already SV-Decomposed.
 *
 * \param U input unitary matrix (RowCount x min(RowCount,ColCount)).
 * \param E vector of singular values, sorted in decreasing order.
 * \param V input unitary matrix (ColCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \throws std::range_error if the dimensions don't match.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, ReadableMatrix Matrix2,
          ReadableMatrix Matrix3, WritableMatrix Matrix4>
void pseudoinvert_SVD(const Matrix1& U, const Matrix2& E, const Matrix3& V,
                      Matrix4& A_pinv,
                      mat_value_type_t<Matrix2> NumTol = 1E-15) {
  if constexpr (!FullyWritableMatrix<Matrix4>) {
    mat<mat_value_type_t<Matrix4>, mat_structure::rectangular> A_pinv_tmp(
        V.get_col_count(), U.get_row_count());
    pseudoinvert_SVD(U, E, V, A_pinv_tmp, NumTol);
    A_pinv = A_pinv_tmp;
  } else {
    if ((U.get_col_count() != E.get_row_count()) ||
        (E.get_row_count() > V.get_row_count())) {
      throw std::range_error(
          "Dimensions of the U E V matrices don't match a singular value "
          "decomposition!");
    }
    using std::abs;

    A_pinv.set_row_count(V.get_row_count());
    A_pinv.set_col_count(U.get_row_count());
    for (int i = 0; i < V.get_row_count(); ++i) {
      for (int j = 0; j < U.get_row_count(); ++j) {
        A_pinv(i, j) = 0.0;
        for (int k = 0; k < E.get_row_count(); ++k) {
          if (abs(E(k, k)) > abs(E(0, 0)) * NumTol) {
            A_pinv(i, j) += V(i, k) * U(j, k) / E(k, k);
          }
        }
      }
    }
  }
}

/**
 * Computes the pseudo-inverse of a NxM matrix A, by performing Singular Value Decomposition.
 *
 * \param A rectangular matrix (RowCount x ColCount).
 * \param A_pinv the pseudo-inverse (ColCount x RowCount).
 * \param NumTol tolerance for considering a value to be zero in avoiding divisions
 *               by zero and singularities.
 *
 * \author Mikael Persson
 */
template <ReadableMatrix Matrix1, WritableMatrix Matrix2>
void pseudoinvert_SVD(const Matrix1& A, Matrix2& A_pinv,
                      mat_value_type_t<Matrix1> NumTol = 1E-15) {
  if constexpr (!FullyWritableMatrix<Matrix2>) {
    mat<mat_value_type_t<Matrix2>, mat_structure::rectangular> A_pinv_tmp(
        A.get_col_count(), A.get_row_count());
    pseudoinvert_SVD(A, A_pinv_tmp, NumTol);
    A_pinv = A_pinv_tmp;
  } else {
    using ValueType = mat_value_type_t<Matrix1>;
    using std::abs;

    int nu = (A.get_row_count() > A.get_col_count() ? A.get_col_count()
                                                    : A.get_row_count());
    mat<ValueType, mat_structure::rectangular> U(A.get_row_count(), nu);
    mat<ValueType, mat_structure::diagonal> E(nu);
    mat<ValueType, mat_structure::rectangular> V(A.get_col_count(), nu);

    decompose_SVD(A, U, E, V, NumTol);

    A_pinv.set_row_count(V.get_row_count());
    A_pinv.set_col_count(U.get_row_count());

    for (int i = 0; i < V.get_row_count(); ++i) {
      for (int j = 0; j < U.get_row_count(); ++j) {
        A_pinv(i, j) = 0.0;
        for (int k = 0; k < E.get_row_count(); ++k) {
          if (abs(E(k, k)) > abs(E(0, 0)) * NumTol) {
            A_pinv(i, j) += V(i, k) * U(j, k) / E(k, k);
          }
        }
      }
    }
  }
}

/**
 * Functor to wrap a call to a SVD-based linear-least-square solver.
 */
struct SVD_linlsqsolver {
  template <typename Matrix1, typename Matrix2, typename Matrix3>
  void operator()(const Matrix1& A, Matrix2& X, const Matrix3& B,
                  mat_value_type_t<Matrix1> NumTol = 1E-8) {
    mat<mat_value_type_t<Matrix1>, mat_structure::rectangular> A_pinv(
        A.get_col_count(), A.get_row_count());
    pseudoinvert_SVD(A, A_pinv, NumTol);
    X = A_pinv * B;
  }
};

}  // namespace ReaK

#endif  // REAK_MATH_LIN_ALG_MAT_SVD_METHOD_H_
