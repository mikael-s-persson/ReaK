# ReaK Library

Version: 0.28.0

Brief description: Software platform and algorithms for multi-body dynamics simulation, 
control, estimation, and path-planning. Intended for robotics software development and testing.


## Table of Contents

- [License Notice](#license-notice)
- [Author](#author)
- [Detailed Description](#detailed-description)
 - [Core Math Utilities](#core-math-utilities)
 - [Multibody Dynamics](#multibody-dynamics)
 - [Control and Motion-planning](#control-and-motion-planning)
- [Installation](#installation)
 - [Dependencies](#dependencies)
 - [Folder Contents](#folder-contents)
- [List of Algorithms and Data-structures](#list-of-algorithms-and-data-structures)


## License Notice

> Copyright 2011 Sven Mikael Persson
> 
> THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
> 
> This program is free software: you can redistribute it and/or modify
> it under the terms of the GNU General Public License as published by
> the Free Software Foundation, either version 3 of the License, or
> (at your option) any later version.
> 
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU General Public License for more details.
> 
> You should have received a copy of the GNU General Public License
> along with this program (as LICENSE in the root folder).  
> If not, see <http://www.gnu.org/licenses/>.

## Author

Main contributor and founder: Sven Mikael Persson, M.Sc.(Tech.) <mikael.s.persson@gmail.com>



## Detailed Description

ReaK (pronounced as the 'reac' in 'reaction' or 'reactor') is a software platform 
which is the result of many years of accumulation of C++ code used by the original 
author (Sven Mikael Persson) for various project in the fields of multibody dynamics 
and control. At the core of the ReaK platform are several general-purpose utilities 
which facilitate serialization / deserialization of objects, memory sharing between 
distributed software modules, run-time type identification, and data input / output. 

### Core Math Utilities

The core math libraries included in ReaK handle basic linear algebra methods (fixed 
and variable size vectors, variable size matrices of various structures and alignments, 
and matrix composition, views and slices), 2D and 3D geometric calculations (rotations 
and kinetostatic frame transformations), matrix numerical methods (LU, Cholesky, QR, 
Jacobi, SVD, PadeSAS matrix exponentials, Redheffer star-product, Algebraic Riccati 
Equations, Schur Decomposition, and matrix norms), 
numerical integration methods (fixed-step, variable-step, and multi-step 
predictor-correctors, both closed-form and iterative), and a set of optimization 
routines. Performance optimization of these libraries is limited to good coding style, 
and thus, do not expect these math libraries to be the fastest available, they were 
designed to be easy to use and interoperable, not for performance-critical applications.

### Multibody Dynamics

The multibody dynamics elements of this library were developed according to the 
Kinetostatic Transmission Elements (KTEs) framework, as originally developed by 
Prof. Andres Kesckemethy at Graz University (now at the University of Duisburg-Essen), 
this is not, however, developed from other existing code that use KTEs, this is an 
original implementation which was done during the course of a Master's degree, by 
the original author, in Space Robotics and Automation at Aalto University, School of 
Science and Technology in Helsinki, Finland. This framework allows serial kinematic 
chains to be modeled in a modular and flexible fashion, to be used for model-based 
control of a robotic system and high-fidelity dynamics simulations [1]. The 
construction of a dynamics model is done via a serial chain of KTEs which model 
simple (or complex) transmission of motion and forces (hence the 'kineto' and 'static' 
in Kinetostatic Transmission Elements). Available KTEs include, but not limited to, 
the following: inertial elements, (torsional) springs, (torsional) dampers, revolute 
joints, prismatic joints, free joints, rigid links, flexible Euler-Bernoulli beams, 
force actuators, driving actuators, geometric constraints (point-on-line and point-on-plane), 
dry and viscous friction, virtual-to-real model interfaces (for Virtual Model Control (VMC)), 
and state measurements and direct controls (no motor/controller model). Additionally, 
some utility classes are available, which are not KTEs but work in parallel with KTEs 
to extract higher-level information about the KTE chain. The main utility class is the 
`mass_matrix_calc` class which can, once given a list of degrees-of-freedom, joint motion 
jacobians and inertial elements, be used to compute the system's mass-matrix as well as 
its composing elements (twist-shaping matrix and aggregate, constant mass matrix), and 
their derivatives (and thus, also the time-derivative of the system mass matrix, which 
is useful in model-based control and estimation).

### Control and Motion-planning

Finally, for the purpose of estimation and motion-planning, the ReaK platform includes 
several generic algorithms and concepts for state representation and estimation, and 
probabilistic motion-planning methods. Note that this part of ReaK is under active 
development (as part of the author's Ph.D. research), so it should be considered as 
very experimental at this stage and incomplete parts are to be expected. First, the 
state estimation algorithms and concepts include several variations of the Kalman 
filtering method, including the '(Extended-)Kalman Filter' ([E]KF), the '(Extended-)Kalman-Bucy 
Filter' ([E]KBF), the 'Hybrid Kalman Filter' (HKF), the 'Unscented Kalman Filter' (UKF), 
the 'Aggregate Kalman Filter' (AKF), the 'Symplectic Kalman Filter' (SKF), and 'Invariant' 
versions of most of the filters (IKF, IKBF, IAKF, and ISKF). The implementations are 
all generic, based on a certain number of concepts (using the 'Boost.Concept-Check' 
library) that define fundamental constructs such as a continuous-time state-space system, 
a continuous-time linear state-space system (including linear-time-invariant (LTI), 
linear-time-varying (LTV), and linearized at the current time, state and input), 
discrete-time versions of those state-space system concepts, an invariant state-space 
system, a covariance matrix representation, and a belief-state representation. Additionally, 
a Gaussian belief-state class is also provided for convenience. Second, the path-planning 
algorithms include several basic implementations and concepts related to classic probabilistic 
path-planning methods and other related utilities. Most algorithms build upon the framework of 
the Boost.Graph library, in the same programming style. Algorithms include the 'Anytime Dynamic A*' (AD*), 
the 'Rapidly-exploring Random Tree' (RRT), the 'Probabilistic Roadmap' (PRM), their optimizing 
variants (RRT* and PRM*), the 'Rapidly-exploring Random Graph' (RRG), the RRT-Connect (or 
bi-directional RRT), the 'Flexible Anytime Dynamic Probabilistic Roadmap' (FADPRM), and 
several variants of a novel algorithm called Sampling-based A* (SBA*) such as an anytime 
version, a simulated annealing version, and a bi-directional version of all of them. 
Several concepts are defined such as a metric space, a temporal space, a reachability space, 
as well as spatial paths and trajectories. Additionally, some of the utilities provided 
include linear-search through topological point-sets (for nearest-neighbor queries, through 
extensive search), a 'Dynamic Vantage-Point Tree' (DVP-tree) implementation for fast 
nearest-neighbor queries via the partitioning of a general metric space (symmetric or 
asymmetric), and a multi-index sorting of points of a reachability 
space. There are also several generic interpolation methods that can operate over any 
Lie Group topology. This part of ReaK is still under active development as of the summer of 2014, 
contributers are welcomed!


## Installation

### Dependencies

External Dependencies:
- CMake build system.
- The Boost Library (www.boost.org), version 1.49.0 or later (note: broken for 1.55 and 1.56).
- BGL "Extras" library (see https://github.com/mikael-s-persson/boost_graph_ext_mp)

Optional External Dependencies (for test-programs):
- OpenCV
- Qt 4.6 or later
- Coin3D and SoQt

### Folder Contents

Folder Descriptions:
- ./src is the top-level folder for all source files and the top-level CMakeLists.txt file.
- ./bin is the folder in which all compiled, executable binary will be put when built with cmake.
- ./lib is the folder in which all compiled libraries will be put when built with cmake.
- ./dox is the folder where the Doxyfile is and is the working directory for generating doxygen documentation for the ReaK platform.
- ./include is the folder where headers for the platform are installed after running "make install"

File Descriptions:
- README.md is this file.
- LICENSE contains the text version of the GNU GPLv3 license agreement.
- TODO_list.txt is the current list of things to be implemented in the near future.


## List of Algorithms and Data-structures

### ReaK.core

- `shared_object`: The stem of a class hierarchy enabling safe object sharing between modules (executables and DLL/.so files).
- RTTI: A complete, template-aware run-time type identification facility with:
  - generic object factories.
  - dynamic casting up, down and across class hierarchies.
  - both intrusive and non-intrusive type-id specifications.
  - cross-module and persistent sharing of type specifications and factories.
- Serialization: A serialization library capable of flattening general object hierarchies:
  - XML output.
  - binary output.
  - Google Protocol Buffer output.
  - both intrusive and non-intrusive interfaces.
- Recorders: A data-streaming library for file and network streaming of floating-point data streams.
  - File formats: binary and text (space-, comma- or tab-separated value files, easily importable in Excel or Matlab).
  - Network protocols: UDP, TCP and "raw" UDP (without meta-data).

### ReaK.math

- Linear Algebra: Complete set of matrix and vector algebra:
  - vector classes:
    - concept-check classes, type-traits and meta-functions for specifying vectors.
    - fixed-size vectors, see `ReaK::vect<T,N>`.
    - variable-size vectors, see `ReaK::vect_n<T>`.
    - scalar vector (vector with all elements equal), see `ReaK::vect_scalar<T,N>`.
    - vector views, as copies or reference wrappers.
  - tuple library extensions:
    - `arithmetic_tuple`: an extension to standard or boost tuple that allows vector-space arithmetic operations.
    - vector-tuple conversions.
  - matrix class template: (`ReaK::mat<T,Struct,Align,Alloc>`):
    - concept-check classes, type-traits and meta-functions for specifying matrices.
    - supports different (sparse) structures:
        - rectangular and square (dense matrices).
        - symmetric and skew-symmetric (stores half the matrix only).
        - diagonal and scalar (stores diagonal only).
        - nil and identity (fixed values, no actual storage, short-cutted operators).
    - supports different data alignments:
        - column-major: stores columns contiguously (default).
        - row-major: stores rows contiguously.
    - all matrix-matrix operations and matrix-vector operations are overloaded for each relevant case (special matrix structures).
  - matrix views and slices:
    - matrix-vector adaptor (make a vector look like a single-row/column matrix), see mat_vector_adaptor.hpp.
    - matrix slices (make a row or column of matrix look like a vector), see mat_slices.hpp.
    - matrix view (make views on sub-matrices of a matrix), see mat_views.hpp.
    - matrix concatenation (assemble large matrices for sub-matrix blocks), see mat_composite_adaptor.hpp.
    - matrix transpose views, see mat_transpose_view.hpp.
- Matrix Numerical Methods:
  - Gaussian elimination-based inversion of matrix, see `ReaK::invert_gaussian()`.
  - PLU Decomposition (square, well-conditioned matrices):
    - solve linear system with PLU-decomposition, see `ReaK::linsolve_PLU()`.
    - invert a well-conditioned matrix with PLU-decomposition, see `ReaK::invert_PLU()`.
  - Cholesky decomposition (symmetric positive-definite matrices):
    - decomposition of symmetric positive-definite matrix, see `ReaK::decompose_Cholesky()`.
    - finding the determinant of symmetric matrix through Cholesky decomp., see `ReaK::determinant_Cholesky()`.
    - solve linear system with Cholesky decomp., see `ReaK::linsolve_Cholesky()`.
    - invert a positive-definite symmetric matrix through Cholesky decomp., see `ReaK::invert_Cholesky()`.
  - Givens rotations:
    - construct stable Givens rotations as 2x2 pseudo-matrices.
    - perform Givens rotations on matrices (or matrix views).
  - Householder reflections:
    - construct stable Householder reflections (and reversed reflections) in any dimensions, as pseudo-matrices.
    - perform Householder reflection products on matrices (or matrix views).
  - Matrix Dampening, see `ReaK::mat_damped_matrix<M1,M2>`.
  - Matrix Balancing:
    - single matrix, see `ReaK::balance()`.
    - matrix pencil, see `ReaK::balance_pencil()`.
  - QR Decomposition (general matrices):
    - decomposition of a general matrix, see `ReaK::decompose_QR()`.
    - finding the determinant of a general matrix through QR decomp., see `ReaK::determinant_QR()`.
    - solve linear least-square with QR decomp., see `ReaK::linlsq_QR()`.
    - solve linear least-square with Rank-Revealing QR decomp., see `ReaK::linlsq_RRQR()`.
    - solve minimum-norm problem with a QR decomp., see `ReaK::minnorm_QR()`.
    - perform a backsubstitution involving a right-triangular matrix, see `ReaK::backsub_R()`.
    - invert a well-conditioned matrix through a QR decomp., see `ReaK::invert_QR()`.
    - pseudo-invert a matrix through QR decomp. (left Moore-Penrose pseudo-inverse), see `ReaK::pseudoinvert_QR()`.
  - Jacobi method (symmetric matrices only):
    - find the determinant of a symmetric matrix through the Jacobi method, see `ReaK::determinant_Jacobi()`.
    - solve linear least-square with the Jacobi method, see `ReaK::linlsq_Jacobi()`.
    - pseudo-invert a matrix through the Jacobi method, see `ReaK::pseudoinvert_Jacobi()`.
    - solve for eigen-values / -vectors through Jacobi method, see `ReaK::eigensolve_Jacobi()`.
  - Singular-value decomposition (SVD) (general matrices):
    - decomposition of a general matrix, see `ReaK::decompose_SVD()`.
    - obtain the numerical condition number of a set of singular-values, see `ReaK::condition_number_SVD()`.
    - obtain the numerical rand of a set of singular-values, see `ReaK::numrank_SVD()`.
    - pseudo-invert a general matrix through SVD, see `ReaK::pseudoinvert_SVD()`.
  - Hessenberg decomposition:
    - decomposition of a general matrix into Hessenberg form, see `ReaK::decompose_Hess()`.
    - reduction of a general matrix pencil into a Hessenberg-Triangular form, see `ReaK::reduce_HessTri()`.
  - Schur decomposition:
    - real-Schur decomposition of a general matrix, see `ReaK::decompose_RealSchur()`.
    - generalized real-Schur decomposition of a general matrix pencil, see `ReaK::decompose_GenRealSchur()`.
  - matrix exponential through Pade Square-And-Sum algorithm, see `ReaK::exp_PadeSAS()`.
  - matrix Redheffer star-product of hamiltonian matrices, see `ReaK::star_product()`.
  - matrix norm calculations (L1, L2, Frobenius, LInf, etc.).
  - Algebraic Riccati Equation solvers:
    - solve the continuous-time algebraic Riccati equation (CARE), see `ReaK::solve_care_problem()`.
    - solve the infinite-horizon, continuous-time LQR problem, see `ReaK::solve_IHCT_LQR()`.
    - solve the infinite-horizon, continuous-time LQG problem, see `ReaK::solve_IHCT_LQG()`.
    - solve the discrete-time algebraic Riccati equation (DARE), see `ReaK::solve_dare_problem()`.
    - solve the infinite-horizon, discrete-time LQR problem, see `ReaK::solve_IHDT_LQR()`.
    - solve the infinite-horizon, discrete-time LQG problem, see `ReaK::solve_IHDT_LQG()`.
    - solve the continuous-time spectral factorization (CTSF), see `ReaK::solve_ctsf_problem()`.
    - solve the discrete-time spectral factorization (DTSF), see `ReaK::solve_dtsf_problem()`.
- Kinetostatics calculations:
  - 2D rotations:
    - 2D rotation matrix representation, see `ReaK::rot_mat_2D<T>`.
    - 2D homogeneous transformation matrix, see `ReaK::trans_mat_2D<T>`.
  - 3D rotations:
    - 3D rotation matrix representation, see `ReaK::rot_mat_3D<T>`.
    - unit-quaternion representation of 3D rotations, see `ReaK::quaternion<T>`.
    - Euler Angles (Tait-Bryant) representation of 3D rotations, see `ReaK::euler_angles_TB<T>`.
    - Axis-angle representation of 3D rotations, see `ReaK::axis_angle<T>`.
    - 3D homogeneous transformation matrix, see `ReaK::trans_mat_3D<T>`.
  - Quaternionic algebra:
    - quaternion class with complete algebraic operations and functions, see `ReaK::quat<T>`.
    - unit-quaternion class with complete algebraic operations and functions, see `ReaK::unit_quat<T>`.
    - inter-operability with 3D rotations and 3D vectors (`ReaK::vect<T,3>`).
  - Kinetostatic frames:
    - 2D/3D pose classes to represent translation and rotation, see `ReaK::pose_2D<T>` and `ReaK::pose_3D<T>`.
    - 2D/3D frame classes to represent full kinematics and statics, see `ReaK::frame_2D<T>` and `ReaK::frame_3D<T>`.
    - hierarchial frame dependencies (relative frames) with kinematic calculations ("rotating frame" formulae).
    - generalized coordinates to represent full kinematics and statics of single-value coordinates, see `ReaK::gen_coord<T>`.
    - complete set of motion-Jacobian representations between any kind of kinetostatic frames (2D, 3D, and generalized).
- Root-finding methods:
  - bisection method (i.e., binary-search for the root), see `ReaK::bisection_method()`.
  - secant methods:
    - secant method (basic variant), see `ReaK::secant_method()`.
    - Illinois method, see `ReaK::illinois_method()`.
    - Ford-3 method, see `ReaK::ford3_method()`.
    - Brent method, see `ReaK::brent_method()`.
    - Ridders method, see `ReaK::ridders_method()`.
  - Broyden methods, see `ReaK::broyden_good_method()` and `ReaK::broyden_fast_method()`.
  - Newton-Raphson method, see `ReaK::newton_raphson_method()`.
- Sorting algorithms:
  - Bubble-sort, see `ReaK::sorting::bubble_sort()`.
  - Insertion-sort, see `ReaK::sorting::insertion_sort()`.
  - Selection-sort, see `ReaK::sorting::selection_sort()`.
  - Comb-sort, see `ReaK::sorting::comb_sort()`.
  - Heap-sort, see `ReaK::sorting::heap_sort()`.
  - Merge-sort, see `ReaK::sorting::merge_sort()`.
  - Shell-sort, see `ReaK::sorting::shell_sort()`.
  - Quick-sort, see `ReaK::sorting::quick_sort()`:
    - first element as pivot, see `ReaK::sorting::first_pivot`.
    - random element as pivot, see `ReaK::sorting::random_pivot`.
    - median-of-3 as pivot, see `ReaK::sorting::median_of_3_pivots`.
    - median-of-3-random as pivot, see `ReaK::sorting::median_of_3_random_pivots`.
  - Intro-sort, see `ReaK::sorting::intro_sort()`.
- Optimization algorithms:
  - Line-search methods (one-dimensional optimization):
    - Dichotomous search, see `ReaK::optim::dichotomous_search()`.
    - Golden-section search, see `ReaK::optim::golden_section_search()`
    - Fibonacci search, see `ReaK::optim::fibonacci_search()`
    - Back-tracking search (used by other general optim. methods), see `ReaK::optim::backtracking_search()`
    - Expand-and-zoom search (used by other general optim. methods), see `ReaK::optim::expand_and_zoom_search()`
  - Linear Programming (LP):
    - primal-dual simplex method (note: not working yet, buggy), see `ReaK::optim::simplex_method()`.
    - Mehrotra's interior-point method (note: not working, unstable), see `ReaK::optim::mehrotra_method()`.
  - Quadratic Programming (QP):
    - Equality-constrained:
        - null-space direct method, see `ReaK::optim::null_space_QP_method()` and `ReaK::optim::null_space_RRQP_method()`.
        - projected conjugate gradient method, see `ReaK::optim::projected_CG_method()`.
        - mehrotra's QP method, see `ReaK::optim::mehrotra_QP_method()`.
    - Inequality-constrained:
        - mehrotra's QP method, see `ReaK::optim::mehrotra_QP_method()`.
  - Non-Linear Least-square
    - Unconstrained (aside from limiters)
        - Gauss-Newton method (performance: reasonable), see `ReaK::optim::gauss_newton_nllsq()` and see `ReaK::optim::limited_gauss_newton_nllsq()`.
        - Jacobian-transpose method (performance: shit), see `ReaK::optim::jacobian_transpose_nllsq()` and see `ReaK::optim::limited_jacobian_transpose_nllsq()`.
        - Levenberg-Marquardt method (DLS with trust-region) (performance: best), see `ReaK::optim::levenberg_marquardt_nllsq()` and see `ReaK::optim::limited_levenberg_marquardt_nllsq()`.
  - Non-Linear Optimization problems
    - Unconstrained (aside from limiters)
        - Nelder-Mead method (performance: bad, expected), see `ReaK::optim::nelder_mead_method()`.
        - Quasi-Newton line-search methods (performance: good, best update is bfgs), see `ReaK::optim::quasi_newton_line_search()` and `ReaK::optim::bfgs_method()`.
        - Quasi-Newton trust-region methods (performance: good, best update is sr1), see `ReaK::optim::quasi_newton_trust_region()` and `ReaK::optim::sr1_tr_method()`.
        - Conjugate-Gradient method (performance: so so), see `ReaK::optim::non_linear_conj_grad_method()`.
        - Newton line-search methods (performance: good, but quasi-newton is more stable), .
        - Newton trust-region methods (performance: good, but quasi-newton is more stable), .
    - Equality-constrained
        - Bound-constrained Newton methods (Augmented Lagrangian methods) (performance: bad, expected), see `ReaK::optim::eq_cnstr_newton_method_tr()` or `ReaK::optim::constraint_newton_method_tr()`.
        - Byrd-Omojokun SQP method (performance: OK), see `ReaK::optim::make_bosqp_newton_tr()` or `ReaK::optim::make_bosqp_quasi_newton_tr()`.
        - Line-search Interior-point method (performance: sucks, must be a bug), see `ReaK::optim::make_nlip_newton_ls()` or `ReaK::optim::make_nlip_quasi_newton_ls()`.
        - Trust-region Interior-point method (performance: good), see `ReaK::optim::make_nlip_newton_tr()` or `ReaK::optim::make_nlip_quasi_newton_tr()`.
    - Inequality-constrained
        - Bound-constrained Newton methods (Augmented Lagrangian methods) (performance: bad, expected), see `ReaK::optim::eq_cnstr_newton_method_tr()` or `ReaK::optim::constraint_newton_method_tr()`.
        - Byrd-Omojokun SQP method (with non-negative limiters) (performance: OK), see `ReaK::optim::make_bosqp_newton_tr()` or `ReaK::optim::make_bosqp_quasi_newton_tr()`.
        - Line-search Interior-point method (performance: sucks, must be a bug), see `ReaK::optim::make_nlip_newton_ls()` or `ReaK::optim::make_nlip_quasi_newton_ls()`.
        - Trust-region Interior-point method (performance: good), see `ReaK::optim::make_nlip_newton_tr()` or `ReaK::optim::make_nlip_quasi_newton_tr()`.
- Numerical integration methods:
  - Fixed-step single-step integration methods:
    - Euler method (forward-euler) (same as matlab ode1), see `ReaK::euler_integrator< T >`.
    - Midpoint method, see `ReaK::midpoint_integrator< T >`.
    - Runge-Kutta 4 method (same as matlab ode4), see `ReaK::runge_kutta4_integrator< T >`.
    - Runge-Kutta 5 method (same as matlab ode5), see `ReaK::runge_kutta5_integrator< T >`.
  - Fixed-step predictor-corrector (multi-step) integration methods:
    - Adams-Bashforth-Moulton method order 3 (same as matlab ode13), see `ReaK::adamsBM3_integrator< T >`.
    - Adams-Bashforth-Moulton method order 5 (same as matlab ode15), see `ReaK::adamsBM5_integrator< T >`.
    - Hamming's modified method, see `ReaK::hamming_mod_integrator< T >`.
    - Hamming's iterated modified method, see `ReaK::hamming_iter_mod_integrator< T >`.
  - Variable-step single-step integration methods:
    - Runge-Kutta-Fehlberg method order 4-5, see `ReaK::fehlberg45_integrator< T >`.
    - Dormand-Prince method order 4-5 (same as matlab ode45), see `ReaK::dormand_prince45_integrator< T >`.

### ReaK.mbd
- Kinetostatic Transmission Elements (KTE) Multibody Dynamics:
  - Elements (in namespace `ReaK::kte`):
    - Linear Spring, see `spring_gen`, `spring_2D` and `spring_3D`.
    - Torsional Spring, see `torsion_spring_2D` and `torsion_spring_3D`.
    - Linear Damper, see `damper_gen`, `damper_2D` and `damper_3D`.
    - Torsional Damper, see `torsion_damper_2D` and `torsion_damper_3D`.
    - Inertia (body mass), see `inertia_gen`, `inertia_2D` and `inertia_3D`.
    - Revolute Joint, see `revolute_joint_2D` and `revolute_joint_3D`.
    - Prismatic Joint, see `prismatic_joint_2D` and `prismatic_joint_3D`.
    - Free Joint, see `free_joint_2D` and `free_joint_3D`.
    - Rigid Link, see `rigid_link_gen`, `rigid_link_2D` and `rigid_link_3D`.
    - Dry Revolute Joint (dry friction), see `dry_revolute_joint_2D` and `dry_revolute_joint_3D`.
    - Flexible Beam (Linearized Euler-Bernoulli), see `flexible_beam_2D` and `flexible_beam_3D`.
    - Inertial Beam (flexible beam with mass), see `inertial_beam_2D` and `inertial_beam_3D`.
    - Driving Actuator (to apply force to joints), see `driving_actuator_gen`, `driving_actuator_2D` and `driving_actuator_3D`.
    - Min-distance to Line, see `line_point_mindist_2D` and `line_point_mindist_3D`.
    - Min-distance to Plane, see `plane_point_mindist_3D`.
    - Sticky Revolute Joint for Virtual Model Control (to compensate stiction friction), see `vmc_revolute_joint_2D` and `vmc_revolute_joint_3D`.
    - Virtual Model Control KTE Interface (traverse the real-to-virtual boundary), see `virtual_kte_interface_gen`, `virtual_kte_interface_2D` and `virtual_kte_interface_3D`.
  - KTE-related Utilities (in namespace `ReaK::kte`):
    - KTE Map Chain (to build a complete model), see `kte_map_chain`.
    - Mass Matrix Calculator, see `mass_matrix_calc`.
    - Direct State Control, see `state_controls.hpp`.
    - Direct State Measures, see `state_measures.hpp`.
    - Manipulator's Kinematics Model (serial chain), see `manipulator_kinematics_model`.
    - Manipulator's Dynamics Model (serial chain), see `manipulator_dynamics_model`.
    - KTE Chain Visitation, see `kte_chain_visitation.hpp`.
    - Jacobian Joint Mappings (to compute jacobians between joints and dependent frames), see `joint_dependent_gen_coord`, `joint_dependent_frame_2D` and `joint_dependent_frame_3D`.
- Dynamic and Kinematic Models (in namespace `ReaK::kte`):
  - Direct/Inverse Kinematics Model (base-classes), see `direct_kinematics_model` and `inverse_kinematics_model`.
  - Inverse Dynamics Model (base-class), see `inverse_dynamics_model`.
  - Free-floating Platform, see `free_floater_2D_kinematics` and `free_floater_3D_kinematics`.
  - Closed-loop Inverse Kinematics Model, see `manip_clik_calculator`.
  - UAV Kinematics, see `UAV_kinematics`.
  - Manipulator Models with Closed-form Inverse Kinematics:
    - 3 Revolute Joint Manipulator (3R, shoulder-elbow-wrist), see `manip_3R_2D_kinematics` and `manip_3R_3D_kinematics`.
    - 6 dof Decoupled Manipulator (3R-3R), see `manip_3R3R_kinematics`.
    - SCARA Manipulator, see `manip_SCARA_kinematics`.
    - SSRMS Space Manipulator (CanadArm-2, 7 dof), see `manip_SSRMS_kinematics`.
    - ERA Space Manipulator (European Robotic Arm, 7 dof), see `manip_ERA_kinematics`.
    - 6 dof Manipulator on a Linear Track (P-3R-3R topology, 7 dof), see `manip_P3R3R_kinematics`.

### ReaK.topologies
- Topologies (in namespace `ReaK::pp`):
  - Concepts:
    - `TopologyConcept` and `LieGroupConcept`.
    - `DistanceMetricConcept`, `MetricSpaceConcept` and `ProperMetricSpaceConcept`.
    - `BoundedSpaceConcept`, `BoxBoundedSpaceConcept` and `SphereBoundedSpaceConcept`.
    - `ProbDistFunctionConcept` and `ProbabilityDistributionConcept`.
    - `RandomSamplerConcept` and `PointDistributionConcept`.
    - `ReachabilitySpaceConcept`.
    - `ReversibleSpaceConcept` (bi-directional or undirected).
    - `SteerableSpaceConcept` (non-trivial point-to-point navigation).
    - `SubSpaceConcept` (sub-space vs. super-spaces).
    - `TangentBundleConcept` (differentiable spaces).
    - `TemporalSpaceConcept` (time-space topology).
    - `BijectionConcept`, `HomeomorphismConcept` and `DiffeomorphismConcept` (all kinds of mappings between topologies).
  - Generic Topologies and Utilities:
    - Vector Topology, see `vector_topology`.
    - Vector Distance Metrics, see `manhattan_distance_metric`, `euclidean_distance_metric`, `inf_norm_distance_metric` and `p_norm_distance_metric`.
    - Temporal Space (and related mappings), see `temporal_space`.
    - Temporal Space Distance Metrics, see `spatial_distance_only` and `time_distance_only`.
    - Metric-space Tuple (glue several topologies together), see `metric_space_tuple`.
    - Metric-space Tuple Distance Metrics, see `manhattan_tuple_distance`, `euclidean_tuple_distance`, `inf_norm_tuple_distance` and `p_norm_tuple_distance`.
    - Differentiable Space and Rate-limited Differentiable Space, see `differentiable_space` and `reach_time_diff_space`.
  - Concrete Topologies and Utilities:
    - Line-segment Topology (1 dof), see `line_segment_topology`.
    - Time Topologies, see `time_topology` and `time_poisson_topology`.
    - Hyper-ball and Hyper-box Topologies, see `hyperbox_topology` and `hyperball_topology`.
    - Direct/Inverse Kinematics Topological Maps (using ReaK.mbd manipulator models), see `manip_direct_kin_map` and `manip_inverse_kin_map`.
    - SE(2) Topologies, see `se2_0th_order_topology`, `se2_1st_order_topology` and `se2_2nd_order_topology`.
    - SO(3) Topologies, see `quaternion_topology`, `so3_0th_order_topology`, `so3_1st_order_topology` and `so3_2nd_order_topology`.
    - SE(3) Topologies, see `se3_0th_order_topology`, `se3_1st_order_topology` and `se3_2nd_order_topology`.
    - N-dof Joint-space Topologies, see `Ndof_space` and `Ndof_rl_space`.
    - N-dof Joint-space Limits, see `Ndof_limits`.
    - N-dof Interpolated Topologies (linear, cubic, quintic, sustained velocity pulse, and sustained acceleration pulse).
    - Manipulator Collision-free Workspace Topologies, see `manip_quasi_static_env` and `manip_dynamic_env`.
    - Proximity-query Applicators, see `proxy_model_applicator` (static) and `proxy_traj_applicator` (dynamic).
    - 2D Point-robot Topology (from environment black-white image, for simple path-planning tests), see `ptrobot2D_test_world`.
- Interpolation (in namespace `ReaK::pp`):
  - Concepts:
    - `InterpolatorConcept`, `LimitedInterpolatorConcept` and `InterpolatorFactoryConcept`.
    - `ExtrapolatorConcept` and `ExtrapolatorFactoryConcept`.
    - Static Paths: `SpatialPathConcept` (random-access) and `SequentialPathConcept` (sequential).
    - Trajectories: `SpatialTrajectoryConcept` (random-access) and `SequentialTrajectoryConcept` (sequential).
    - `PredictedTrajectoryConcept`.
  - Generic Utilities:
    - Constant Trajectory, see `constant_trajectory`.
    - Waypoint Container (base-class for paths and trajectories built on waypoints), see `waypoint_container`.
    - Discrete Point Paths / Trajectories, see `discrete_point_path` and `discrete_point_trajectory`.
    - Point-to-Point Paths / Trajectories, see `point_to_point_path` and `point_to_point_trajectory`.
    - Generic Interpolator Generators (apply interpolators to complex topologies), see `generic_interpolator`.
    - Interpolated Topology / Trajectories, see `interpolated_topology` and `interpolated_trajectory`.
    - Object-Oriented Bindings for Paths and Trajectories, see `[seq_]path_base`, `[seq_]path_wrapper`, `[seq_]trajectory_base` and `[seq_]trajectory_wrapper`.
    - Transformed Trajectory (via topological mapping), see `transformed_trajectory`.
  - Concrete Interpolation Methods (note, can be combined with most generic utilities listed above):
    - Linear Interpolation (1st-order, 1st-degree), see `linear_interpolator`.
    - Cubic Hermite-spline Interpolation (c-spline) (2nd-order, 3rd-degree), see `cubic_hermite_interpolator`.
    - Quintic Hermite-spline Interpolation (q-spline) (3rd-order, 5th-degree), see `quintic_hermite_interpolator`.
    - N-dof Sustained Velocity Pulse (2nd-order continuous, limited velocity / accel), see `svp_Ndof_interpolator`.
    - N-dof Sustained Acceleration Pulse (3rd-order continuous, limited velocity / accel / jerk), see `sap_Ndof_interpolator`.
    - N-dof SVP / SAP Metrics (optimal travel-time) and Samplers (motion constraints respected while sampling points in the topology).

### ReaK.ctrl_est
- State-space system library:
  - state-space system concepts:
    - continuous-time state-space system (input-output system):
        - linear time-invariant (LTI)
        - linear time-varying (LTV)
        - linearized
        - non-linear
    - discrete-time state-space system (input-output system)
        - linear time-invariant (LTI)
        - linear time-varying (LTV)
        - linearized
        - non-linear
    - invariant system concept.
  - generic state-space systems:
    - continuous-time LTI state-space system (vector-space).
    - discrete-time LTI state-space system (vector-space).
    - discretized LTI state-space system (vector-space).
    - numerically-integrated non-linear discrete-time system.
    - KTE-based non-linear system.
- State-estimation library:
  - state estimation concepts:
    - Gaussian belief-state concepts
    - covariance matrix concept
    - state-estimator concept
  - covariance representations:
    - covariance matrix
    - information matrix
    - decomposed covariance matrix
    - Gaussian belief-space (topology of Gaussian belief-states).
  - continuous-time Kalman filters:
    - (Extended-)Kalman-Bucy Filter (note: not tested).
    - Invariant (Extended-)Kalman-Bucy Filter (note: not tested).
    - Hybrid (Extended-)Kalman-Bucy Filter (note: not tested).
  - discrete-time Kalman filters:
    - (Extended-)Kalman Filter (EKF).
    - Aggregate Kalman Filter (or Hamiltonian Kalman Filter).
    - Symplectic Kalman Filter (based on Symplectic covariance mappings).
    - Unscented Kalman Filter (note: not tested).
    - Invariant (Extended-)Kalman Filter (IEKF).
    - Invariant Symplectic Kalman Filter.
    - Invariant Aggregate Kalman Filter (note: not useful).


