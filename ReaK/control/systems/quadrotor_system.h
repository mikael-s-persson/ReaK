/**
 * \file quadrotor_system.h
 *
 * This library implements a state-space system for a quad-rotor aircraft.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_CONTROL_SYSTEMS_QUADROTOR_SYSTEM_H_
#define REAK_CONTROL_SYSTEMS_QUADROTOR_SYSTEM_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/lin_alg/mat_alg.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include "ReaK/topologies/spaces/se3_topologies.h"

namespace ReaK::ctrl {

/**
 * This class implements a state-space system for a quad-rotor aircraft. The modeling
 * includes rotational and translational drag. This class models SSSystemConcept and
 * LinearSSSystemConcept(LinearizedSystemType).
 */
class quadrotor_system : public named_object {
 public:
  using state_space_type = pp::se3_1st_order_topology<double>::type;

  using point_type = pp::topology_traits<state_space_type>::point_type;
  using point_difference_type =
      pp::topology_traits<state_space_type>::point_difference_type;
  using point_derivative_type = point_difference_type;

  using time_type = double;
  using time_difference_type = double;

  using input_type = vect<double, 4>;
  using output_type = point_type;

  using invariant_error_type = point_difference_type;
  using invariant_correction_type = point_difference_type;
  using invariant_frame_type = mat<double, mat_structure::identity>;

  static constexpr std::size_t dimensions = 13;
  static constexpr std::size_t input_dimensions = 4;
  static constexpr std::size_t output_dimensions = 13;
  static constexpr std::size_t invariant_error_dimensions = 12;
  static constexpr std::size_t invariant_correction_dimensions = 12;

  using matrixA_type = mat<double, mat_structure::square>;
  using matrixB_type = mat<double, mat_structure::rectangular>;
  using matrixC_type = mat<double, mat_structure::rectangular>;
  using matrixD_type = mat<double, mat_structure::nil>;

 protected:
  double mMass;
  mat<double, mat_structure::symmetric> mInertiaMoment;
  mat<double, mat_structure::symmetric> mInertiaMomentInv;

  mat<double, mat_structure::diagonal> mTransDragCoefs;
  mat<double, mat_structure::diagonal> mRotDragCoefs;

 public:
  /**
   * Returns the dimensions of the states of the system.
   * \return The dimensions of the states of the system.
   */
  std::size_t get_state_dimensions() const { return 13; }

  /**
   * Returns the dimensions of the input of the system.
   * \return The dimensions of the input of the system.
   */
  std::size_t get_input_dimensions() const { return 4; }

  /**
   * Returns the dimensions of the output of the system.
   * \return The dimensions of the output of the system.
   */
  std::size_t get_output_dimensions() const { return 13; }

  /**
   * Returns the dimensions of the invariant errors of the system.
   * \return The dimensions of the invariant errors of the system.
   */
  std::size_t get_invariant_error_dimensions() const { return 12; }

  /**
   * Returns the dimensions of the corrections to the states of the system.
   * \return The dimensions of the corrections to the states of the system.
   */
  std::size_t get_correction_dimensions() const { return 12; }

  /**
   * Default and parametrized constructor.
   * \param aName The name of this object.
   * \param aMass The mass of the quad-rotor.
   * \param aInertiaMoment The inertia-moment (inertia tensor) of the quad-rotor.
   * \param aTransDragCoefs The body-fixed translational drag coefficients of the quad-rotor.
   * \param aRotDragCoefs The body-fixed rotational drag coefficients of the quad-rotor.
   */
  explicit quadrotor_system(
      const std::string& aName, double aMass = 1.0,
      const mat<double, mat_structure::symmetric>& aInertiaMoment =
          (mat<double, mat_structure::symmetric>(
              mat<double, mat_structure::identity>(3))),
      const mat<double, mat_structure::diagonal>& aTransDragCoefs =
          (mat<double, mat_structure::diagonal>(3, 0.5)),
      const mat<double, mat_structure::diagonal>& aRotDragCoefs =
          (mat<double, mat_structure::diagonal>(3, 0.5)));

  quadrotor_system() : quadrotor_system("") {}

  ~quadrotor_system() override = default;

  /**
   * This function computes the state-vector derivative at a given point in state (x), input (u), and time (t).
   * \param space The state-space object (ignored).
   * \param x The current state of the quad-rotor.
   * \param u The current input to the quad-rotor.
   * \param t The current time (ignored).
   * \return The state-vector derivative at a given point in state (x), input (u), and time (t).
   */
  point_derivative_type get_state_derivative(const state_space_type& space,
                                             const point_type& x,
                                             const input_type& u,
                                             time_type t = 0.0) const;

  /**
   * This function returns the output-vector at a given state (x).
   * \param space The state-space object (ignored).
   * \param x The current state of the quad-rotor.
   * \param u The current input to the quad-rotor (ignored).
   * \param t The current time (ignored).
   * \return The output-vector at a given state (x).
   */
  output_type get_output(const state_space_type& space, const point_type& x,
                         const input_type& u, const time_type t = 0.0) const {
    RK_UNUSED(space);
    RK_UNUSED(u);
    RK_UNUSED(t);
    return x;
  }

  /**
   * This function fills the linearization matrices (A,B,C,D) at a given point in state (x), input (u), and time (t).
   * \param A Holds, as output, the jacobian matrix from the state to the state-derivative.
   * \param B Holds, as output, the jacobian matrix from the input to the state-derivative.
   * \param C Holds, as output, the jacobian matrix from the state to the output.
   * \param D Holds, as output, the jacobian matrix from the input to the output.
   * \param space The state-space object (ignored).
   * \param t The current time (ignored).
   * \param u The current input to the quad-rotor.
   * \param x The current state of the quad-rotor.
   */
  void get_linear_blocks(matrixA_type& A, matrixB_type& B, matrixC_type& C,
                         matrixD_type& D, const state_space_type& space,
                         const time_type& t, const point_type& x,
                         const input_type& u) const;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override;
  void load(ReaK::serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(quadrotor_system, 0xC2310005, 1,
                              "quadrotor_system", named_object)
};

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_QUADROTOR_SYSTEM_H_
