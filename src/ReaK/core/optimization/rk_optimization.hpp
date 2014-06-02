/**
 * \file rk_optimization.hpp
 *
 * This library implements a set of optimization routines for linear and non-linear
 * equations. Some implementations come from textbooks (pseudo-code translated to
 * C++ code) and some come from porting free source codes.
 *
 * \author Mikael Persson, B.Eng. Honours student in Mech. Eng., McGill University.
 * \date december 2007
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

#include "rk_mat_num.hpp"

/// Return code for an infeasible problem, i.e. the constrained space is empty.
#define OPT_ERR_INFEASIBLE -1
/// Return code for an unbounded problem, i.e. the constrained space opens to infinity towards the optimum.
#define OPT_ERR_UNBOUNDED  -2
/// Return code for an improper problem, i.e. the input to the optimization are problematic.
#define OPT_ERR_IMPROPER   -3
/// Return code to indicate that the maximum number of iterations was reached.
#define OPT_MAX_ITERATION  2
/// Returm code to indicate that an optimum has been reached.
#define OPT_OPTIMAL 1

#ifndef RK_OPTIMIZATION_HPP
#define RK_OPTIMIZATION_HPP

///This type is a function pointer to a function which maps an input vector u to an output vector y.
typedef __stdcall void (*FMapFunction) (void* UserData,void* u,unsigned int u_count,void* y,unsigned int y_count);
///This type is a function pointer to a function which maps an input vector u to an output vector y, with the weights (or parameters) w.
typedef __stdcall void (*FBasisFunction) (void* UserData,void* u,unsigned int u_count,void* w, unsigned int w_count,void* y,unsigned int y_count);

///This structure holds the information relative to a basis map.
typedef struct {
  unsigned int input_count; ///< The size of the input vector.
  unsigned int basis_count; ///< The number of weights or parameters.
  unsigned int output_count; ///< The size of the output vector.
  FBasisFunction func; ///< The function pointer for the map.
  void* weights; ///< The current parameter values or estimates (i.e. for parameter identification).
} TBasisMap;

///This structure holds the information relative to a non-linear map.
typedef struct {
  unsigned int input_count; ///< The size of the input vector.
  unsigned int output_count; ///< The size of the output vector.
  FMapFunction func; ///< The function pointer for the map.
} TNonLinMap;

/**
 * Forward finite 2-point difference approximation to the parameter-jacobian of
 * a set of basis functions.
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see finite_diff_jacobians.hpp).
 */
template <class T>
void WeightJacobian2PtsForward(void* UserData, ///< The pointer that is passed uniterpreted to the map function.
                               TBasisMap Map, ///< Basis map to be differentiated.
                               T* u, ///< Input to provide uninterpreted to the function.
                               T* hx, ///< Output of the function with the current weights.
                               T* hxx, ///< Work array for evaluating function(p+delta).
                               T delta, ///< The increment for computing the jacobian.
                               T* jac, ///< The pointer to the memory in which to store the jacobian (size Map.weight_count x Map.output_count).
                               bool rowMajor = false);

/**
 * Central finite 2-point difference approximation to the parameter-jacobian of
 * a set of basis functions.
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see finite_diff_jacobians.hpp).
 */
template <class T>
void WeightJacobian2PtsCentral(void* UserData, ///< The pointer that is passed uniterpreted to the map function.
                               TBasisMap Map, ///< Basis map to be differentiated.
                               T* u, ///< Input to provide uninterpreted to the function.
                               T* hxm, ///< Work array for evaluating function(p-delta).
                               T* hxp, ///< Work array for evaluating function(p+delta).
                               T delta, ///< The increment for computing the jacobian.
                               T* jac, ///< The pointer to the memory in which to store the jacobian (size Map.weight_count x Map.output_count).
                               bool rowMajor = false);

/**
 * This function calculates the Jacobian of a TBasisMapf using a 5 point central
 * difference method. This has the advantage of computing the values only four
 * times (for each partial derivative) and still gaining an order 5 estimation.
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see finite_diff_jacobians.hpp).
 */
template <class T>
void Jacobian5Pts(void* UserData, ///< The pointer that is passed uniterpreted to the map function.
                  TBasisMap Map, ///< The map information.
                  T* u, ///< The input vector at which to evaluate the jacobian (size Map.input_count).
                  T du, ///< The increment to give to the input vector for the finite difference.
                  T* Jac, ///< The pointer to the memory in which to store the jacobian (size Map.input_count x Map.output_count).
                  bool rowMajor = false);

/**
 * This function calculates the Jacobian of a TNonLinMapf using a 5 point central
 * difference method. This has the advantage of computing the values only four
 * times (for each partial derivative) and still gaining an order 5 estimation.
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see finite_diff_jacobians.hpp).
 */
template <class T>
void Jacobian5Pts(void* UserData, ///< The pointer that is passed uniterpreted to the map function.
                  TNonLinMap Map, ///< The map information.
                  T* u, ///< The input vector at which to evaluate the jacobian (size Map.input_count).
                  T du, ///< The increment to give to the input vector for the finite difference.
                  T* Jac, ///< The pointer to the memory in which to store the jacobian (size Map.input_count x Map.output_count).
                  bool rowMajor = false);

/**
 * This function is an implementation of the general two-phase revised simplex
 * method for bounded variables. It solves the following problem: \n
 * \n
 *           max c'x \n
 *               Ax = b \n
 *             l <= x <= u \n
 * \n
 * The implementation was inspired from the algorithm described in the book:\n
 *   Chvatal, Vasek, Linear Programming, W. H. Freeman and Company, 1983.
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see simplex_method.hpp).
 */
template <class T>
int SimplexMethod(T* A, ///< The constraint matrix of dimension M*N.
                  T* b, ///< The b vector of dimension M.
                  unsigned int M, ///< The size of the b vector.
                  T* c, ///< The cost vector of dimension N.
                  unsigned int N, ///< The size of the c vector.
                  T* x0, ///< The initial guess for the optimal vector of dimension N.
                  T* l, ///< The lower bound vector, if there is no lower bound for a variable, set to -RK_F_INF.
                  T* u ///< The upper bound vector, if there is no upper bound for a variable, set to RK_F_INF.
                  );

//int QuadraticMethod(/* params */);   //Mehrotra predictor-corrector method

/**
 * This function solves the linear least-square fit. In other words, it solves
 * the following linear system:\n
 *                                    Ax ~ b \n
 * by computing the least-square approximation:\n\n   A'Ax = A'b \n\n
 *
 * The method used for the solution is the QRDecomposition which is numerically
 * very stable and not too computationally intensive. The QR decomposition yields
 * the following compuation:\n\n   Rx = Q'b (where R is upper triangular)
 *
 * \author Mikael Persson
 * \note This function has been ported to the new framework (see mat_num.hpp).
 */
template <class T>
void LinearLeastSquare(T* A, ///< A matrix of M rows and N columns, where M > N for a solution to exist.
                       unsigned int N, ///< The size of the vector x.
                       unsigned int M, ///< The size of the vector b.
                       T* x, ///< The x vector, to be approximated and store the result of the algorithm.
                       T* b ///< The b vector.
                       );

/**
 * This function solves the linear least-square curve fitting. In other words, it
 * finds the best approximation for the weights of a linear combination of basis
 * functions. This function basically calls the basis map with the one value for
 * each weight consecutively to construct the matrix A of the linear least-square
 * approximation algorithm, then it calls LinearLeastSquare.
 *
 * \author Mikael Persson
 */
template <class T>
void LinearLeastSquare(void* UserData, ///< The pointer that is passed uninterpreted to the basis map.
                       TBasisMap Basis, ///< The basis map, with the weight vector which will hold the result.
                       T* input, ///< The input vector of dimension Basis.input_count.
                       T* output ///< The output vector of dimension Basis.output_count.
                       );

/**
 * Check the weight-jacobian of basis nonlinear functions in input_count variables
 * evaluated at current parameters, for consistency with the function itself.
 *
 * Based on fortran77 subroutine CHKDER by
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
 * Argonne National Laboratory. MINPACK project. March 1980.
 *
 * The function does not perform reliably if cancellation or
 * rounding errors cause a severe loss of significance in the
 * evaluation of a function. therefore, none of the parameters
 * should be unusually small (in particular, zero) or any
 * other value which may cause loss of significance.
 */
template <class T>
void CheckWeightJacobian(void* UserData, ///< The pointer that is passed uniterpreted to the map function.
                         TBasisMap Map, ///< Map which computes the derivated function. Holds the current weights.
                         FBasisFunction Jacf, ///< Function which computes the derivated function. The output is Map.output_count x Map.basis_count.
                         T* u, ///< Input to provide uninterpreted to the function.
                         T* err /**< array of length output_count. On output, err contains measures
                                     * of correctness of the respective gradients. if there is
                                     * no severe loss of significance, then if err[i] is 1.0 the
                                     * i-th gradient is correct, while if err[i] is 0.0 the i-th
                                     * gradient is incorrect. For values of err between 0.0 and 1.0,
                                     * the categorization is less certain. In general, a value of
                                     * err[i] greater than 0.5 indicates that the i-th gradient is
                                     * probably correct, while a value of err[i] less than 0.5
                                     * indicates that the i-th gradient is probably incorrect.
                                     */
                         );

/**
 * This function computes in C the covariance matrix corresponding to a least
 * squares fit. JtJ is the approximate Hessian at the solution (i.e. J^T*J, where
 * J is the jacobian at the solution), sumsq is the sum of squared residuals
 * (i.e. goodnes of fit) at the solution, M is the number of parameters (variables)
 * and N the number of observations. JtJ can coincide with C.
 *
 * if JtJ is of full rank, C is computed as sumsq/(n-m)*(JtJ)^-1
 * otherwise C=sumsq/(n-r)*(JtJ)^+ where r is JtJ's rank and ^+ denotes
 * the pseudoinverse. The diagonal of C is made up from the estimates of
 * the variances of the estimated regression coefficients.
 * See the documentation of routine E04YCF from the NAG fortran lib
 *
 * The function returns the rank of JtJ if successful, 0 on error
 *
 * JtJ and C are MxM
 *
 */
template <class T>
int ComputeCovariance(T* JtJ, T* C, T sumsq, unsigned int N, unsigned int M);

#define LM_INIT_MU       1E-03
#define LM_STOP_THRESH   1E-17
#define LM_DIFF_DELTA    1E-06

#define ONE_THIRD     0.3333333334 /* 1.0/3.0 */



/**
 * This function seeks the parameter vector p that best describes the measurements vector x.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized.
 *
 * This function requires an analytic jacobian.
 *
 * Returns the number of iterations (>=0) if successfull, -1 if failed
 *
 * For more details, see H.B. Nielsen's (http://www.imm.dtu.dk/~hbn) IMM/DTU
 * tutorial at http://www.imm.dtu.dk/courses/02611/nllsq.pdf
 */
template <class T>
int LMNonLinLsq(void* UserData, ///< A user pointer to be provided uninterpreted to the function and jacobian.
                TBasisMap Map, /**< The functional relation describing measurements. A p \in R^Map.basis_count yields a \hat{x} \in  R^Map.output_count */
                FBasisFunction Jacf, /**< The function to evaluate the weight jacobian \part x / \part p */
                T* u, ///< The input to the measurement vector to be provided uninterpreted to Map.
                T* x, /**< The actual measured results (x:outputs). */
                unsigned int itmax, /**< The maximum number of iterations. */
                T tau, ///< The scale factor for the initial mu.
                T epsj, ///< Stopping threshold for the jacobian norm.
                T epsw, ///< Stopping threshold for the weight changes.
                T epsx, ///< Stopping threshold for the error change.
                T info[9],     /**< The information regarding the minimization. Set to NULL if don't care
                                    * info[0]= ||e||_2 at initial weights.
                                    * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated weights.
                                    * info[5]= # iterations,
                                    * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                                    *                                 2 - stopped by small D(weights)
                                    *                                 3 - stopped by itmax
                                    *                                 4 - singular matrix. Restart from current p with increased mu
                                    *                                 5 - no further error reduction is possible. Restart with increased mu
                                    *                                 6 - stopped by small ||e||_2
                                    * info[7]= # function evaluations
                                    * info[8]= # jacobian evaluations
                                    */
                T* covar);    /**< Covariance matrix corresponding to LS solution; Map.basis_count x Map.basis_count. Set to NULL if not needed. */

/** Secant version of the LEVMAR_DER() function above: the jacobian is approximated with
 * the aid of finite differences (forward or central, see the comment for the opts argument)
 */
template <class T>
int LMNonLinLsq(void* UserData, ///< A user pointer to be provided uninterpreted to the function.
                TBasisMap Map, /**< The functional relation describing measurements. A p \in R^Map.basis_count yields a \hat{x} \in  R^Map.output_count */
                T* u, ///< The input to the measurement vector to be provided uninterpreted to Map.
                T* x, /**< The actual measured results (x:outputs). */
                int itmax, /**< The maximum number of iterations. */
                T tau, ///< The scale factor for the initial mu.
                T epsj, ///< Stopping threshold for the jacobian norm.
                T epsw, ///< Stopping threshold for the weight change.
                T epsx, ///< Stopping threshold for the error change.
                T delta, ///< The step to take for finite difference jacobian.
                T info[9], /**< The information regarding the minimization. Set to NULL if don't care
                                * info[0]= ||e||_2 at initial p.
                                * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                                * info[5]= # iterations,
                                * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                                *                                 2 - stopped by small Dp
                                *                                 3 - stopped by itmax
                                *                                 4 - singular matrix. Restart from current p with increased mu
                                *                                 5 - no further error reduction is possible. Restart with increased mu
                                *                                 6 - stopped by small ||e||_2
                                * info[7]= # function evaluations
                                * info[8]= # jacobian evaluations
                                */
                T* covar); /**< Covariance matrix corresponding to LS solution; Map.basis_count x Map.basis_count. Set to NULL if not needed. */

/** This function projects a vector p to a box shaped feasible set. p is a Nx1 vector.
 * Either lb, ub can be NULL. If not NULL, they are Nx1 vectors
 */
template <class T>
void ProjectOnBox(T* p, T* lb, T* ub, unsigned int N);

/** This function checks box constraints for consistency */
template <class T>

int CheckBoxConsistency(T* lb, T* ub, unsigned int N);

/**
 * This function seeks the parameter vector p that best describes the measurements
 * vector x under box constraints.
 * More precisely, given a vector function  func : R^m --> R^n with n>=m,
 * it finds p s.t. func(p) ~= x, i.e. the squared second order (i.e. L2) norm of
 * e=x-func(p) is minimized under the constraints lb[i]<=p[i]<=ub[i].
 * If no lower bound constraint applies for p[i], use -DBL_MAX/-FLT_MAX for lb[i];
 * If no upper bound constraint applies for p[i], use DBL_MAX/FLT_MAX for ub[i].
 *
 * This function requires an analytic jacobian.
 *
 * Returns the number of iterations (>=0) if successfull, -1 if failed
 *
 * For details, see C. Kanzow, N. Yamashita and M. Fukushima: "Levenberg-Marquardt
 * methods for constrained nonlinear equations with strong local convergence properties",
 * Journal of Computational and Applied Mathematics 172, 2004, pp. 375-397.
 * Also, see H.B. Nielsen's (http://www.imm.dtu.dk/~hbn) IMM/DTU tutorial on
 * unconrstrained Levenberg-Marquardt at http://www.imm.dtu.dk/courses/02611/nllsq.pdf
 */
template <class T>
int LMBoxNonLinLsq(void* UserData, ///< A user pointer to be provided uninterpreted to the function.
  TBasisMap Map, /**< The functional relation describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  FBasisFunction Jacf, /**< The function to evaluate the weight jacobian \part x / \part p */
  T* u, ///< The input to the measurement vector to be provided uninterpreted to Map.
  T* x, /**< The actual measured results (x:outputs). */
  T* lb, /**< The vector of lower bounds. If NULL, no lower bounds apply */
  T* ub, /**< The vector of upper bounds. If NULL, no upper bounds apply */
  int itmax, /**< The maximum number of iterations. */
  T tau, ///< The scale factor for the initial mu.
  T epsj, ///< Stopping threshold for the jacobian norm.
  T epsw, ///< Stopping threshold for the weight change.
  T epsf, ///< Stopping threshold for the error change.
  T info[9], /**< The information regarding the minimization. Set to NULL if don't care
                  * info[0]= ||e||_2 at initial p.
                  * info[1-4]=[ ||e||_2, ||J^T e||_inf,  ||Dp||_2, mu/max[J^T J]_ii ], all computed at estimated p.
                  * info[5]= # iterations,
                  * info[6]=reason for terminating: 1 - stopped by small gradient J^T e
                  *                                 2 - stopped by small Dp
                  *                                 3 - stopped by itmax
                  *                                 4 - singular matrix. Restart from current p with increased mu
                  *                                 5 - no further error reduction is possible. Restart with increased mu
                  *                                 6 - stopped by small ||e||_2
                  * info[7]= # function evaluations
                  * info[8]= # jacobian evaluations
                  */
  T* covar);    /**< Covariance matrix corresponding to LS solution; mxm. Set to NULL if not needed. */

/**
 * Limited Memory BFGS Method For Large Scale Optimization\n
 * \n
 * \author Jorge Nocedal
 * \n
 * The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
 * Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
 * of memory.
 *
 * The subroutine generates the approximation of an inverse Hessian matrix by
 * using information about the last M steps of the algorithm  (instead of N).
 * It lessens a required amount of memory from a value  of  order  N^2  to  a
 * value of order 2*N*M.
 */
template <class T>
int LBFGSMinimize(void* UserData, ///< The pointer that is passed uninterpreted to the basis map.
                  TBasisMap Map, ///< Map which computes the cost function.
                  FBasisFunction Grad, ///< Function which will evaluate the gradient vector.
                  unsigned int M, ///< Number of corrections in the BFGS scheme of Hessian approximation update. Recommended value:  3<=M<=7. The smaller value causes worse convergence, the bigger will not cause a considerably better convergence, but will cause a fall in the performance. M<=N.
                  T* x, ///< Initial solution approximation and also stores the final result.
                  T epsg, ///< Tolerance for the euclidian norm of the gradient vector
                  T epsf, ///< Relative (or absolute for less than 1) tolerance for the change in the cost function.
                  T epsx, ///< Tolerance for the step size.
                  unsigned int maxits ///< Maximum number of iterations.
                  );

//int SecantRootFinding(/* params */);
//int NewtonRootFinding(/* params */);

#endif //RK_OptimizationH

#include <RKOptimization.cpp>

//------------------------------------------------------------------------------
