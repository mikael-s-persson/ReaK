

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

#include <ReaK/control/systems/quadrotor_system.hpp>

#include <ReaK/control/systems/sss_exceptions.hpp>

#include <ReaK/math/lin_alg/mat_cholesky.hpp>

#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/kinetostatics/rotations_3D.hpp>

#include <ReaK/topologies/spaces/vector_topology.hpp>

namespace ReaK::ctrl {

quadrotor_system::quadrotor_system(
    const std::string& aName, double aMass,
    const mat<double, mat_structure::symmetric>& aInertiaMoment,
    const mat<double, mat_structure::diagonal>& aTransDragCoefs,
    const mat<double, mat_structure::diagonal>& aRotDragCoefs)
    : mMass(aMass),
      mInertiaMoment(aInertiaMoment),
      mTransDragCoefs(aTransDragCoefs),
      mRotDragCoefs(aRotDragCoefs) {
  setName(aName);
  if ((mInertiaMoment.get_row_count() != 3) ||
      (mMass < std::numeric_limits<double>::epsilon())) {
    throw system_incoherency(
        "Inertial information is improper in airship3D_lin_system's "
        "definition");
  }
  try {
    invert_Cholesky(mInertiaMoment, mInertiaMomentInv);
  } catch (singularity_error&) {
    throw system_incoherency(
        "Inertial tensor is singular in airship3D_lin_system's definition");
  }
}

quadrotor_system::point_derivative_type quadrotor_system::get_state_derivative(
    const state_space_type& /*unused*/, const point_type& x,
    const input_type& u, time_type /*unused*/) const {
  using std::abs;

  quaternion<double> q = get_quaternion(x).as_rotation();
  vect<double, 3> w = get_ang_velocity(x);
  vect<double, 3> w_sqr(w[0] * abs(w[0]), w[1] * abs(w[1]), w[2] * abs(w[2]));
  vect<double, 3> aacc =
      mInertiaMomentInv * (vect<double, 3>(u[1], u[2], u[3]) -
                           w % (mInertiaMoment * w) - mRotDragCoefs * w_sqr);

  vect<double, 3> local_v = invert(q) * get_velocity(x);
  vect<double, 3> v = get_velocity(x);
  local_v[0] *= abs(local_v[0]) / mMass;
  local_v[1] *= abs(local_v[1]) / mMass;
  local_v[2] *= abs(local_v[2]) / mMass;
  return {
      make_arithmetic_tuple(v,  // velocity -> derivative of position
                            vect<double, 3>(0.0, 0.0, 9.81) -
                                q * (mTransDragCoefs * local_v +
                                     vect<double, 3>(0.0, 0.0, u[0] / mMass))),
      make_arithmetic_tuple(
          w,  // angular velocity -> invariant derivative of rotation
          aacc)};
}

void quadrotor_system::get_linear_blocks(matrixA_type& A, matrixB_type& B,
                                         matrixC_type& C, matrixD_type& D,
                                         const state_space_type& /*unused*/,
                                         const time_type& /*unused*/,
                                         const point_type& x,
                                         const input_type& u) const {
  using std::abs;

  quaternion<double> q = get_quaternion(x).as_rotation();
  mat<double, mat_structure::square> R = q.getMat();

  A = mat<double, mat_structure::nil>(12, 12);

  // velocity to position partial derivative:
  A(0, 3) = 1.0;
  A(1, 4) = 1.0;
  A(2, 5) = 1.0;

  vect<double, 3> local_v = invert(q) * get_velocity(x);
  mat<double, mat_structure::diagonal> dV(vect<double, 3>(
      -2.0 * abs(local_v[0]) / mMass, -2.0 * abs(local_v[1]) / mMass,
      -2.0 * abs(local_v[2]) / mMass));

  // velocity - velocity partial derivative:
  set_block(A, R * mTransDragCoefs * dV, 3, 3);

  // velocity - quaternion partial derivative:
  local_v[0] *= abs(local_v[0]) / mMass;
  local_v[1] *= abs(local_v[1]) / mMass;
  local_v[2] *= abs(local_v[2]) / mMass;
  set_block(A,
            R * (mat<double, mat_structure::skew_symmetric>(mTransDragCoefs *
                                                            local_v) -
                 mTransDragCoefs *
                     mat<double, mat_structure::skew_symmetric>(local_v) +
                 mat<double, mat_structure::skew_symmetric>(
                     vect<double, 3>(0.0, 0.0, u[0] / mMass))),
            3, 6);

  // angular velocity to quaternion partial derivative:
  A(6, 9) = 1.0;
  A(7, 10) = 1.0;  // identity
  A(8, 11) = 1.0;

  vect<double, 3> w = get_ang_velocity(x);
  mat<double, mat_structure::diagonal> dW(
      vect<double, 3>(-2.0 * abs(w[0]), -2.0 * abs(w[1]), -2.0 * abs(w[2])));
  set_block(
      A,
      mInertiaMomentInv *
          (mRotDragCoefs * dW -
           mat<double, mat_structure::skew_symmetric>(w) * mInertiaMoment +
           mat<double, mat_structure::skew_symmetric>(mInertiaMoment * w)),
      9, 9);

  B = mat<double, mat_structure::nil>(12, 4);
  vect<double, 3> t_global = R * vect<double, 3>(0.0, 0.0, -1.0 / mMass);
  set_block(B, mat_vect_adaptor<vect<double, 3>>(t_global), 3, 0);
  set_block(B, mInertiaMomentInv, 9, 1);

  C = mat<double, mat_structure::nil>(6, 12);
  set_block(C, mat<double, mat_structure::identity>(6), 0, 0);

  D = mat<double, mat_structure::nil>(6, 6);
}

/*******************************************************************************
                ReaK's RTTI and Serialization interfaces
*******************************************************************************/

void quadrotor_system::save(ReaK::serialization::oarchive& A,
                            unsigned int /*unused*/) const {
  named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(mMass) &
      RK_SERIAL_SAVE_WITH_NAME(mInertiaMoment) &
      RK_SERIAL_SAVE_WITH_NAME(mTransDragCoefs) &
      RK_SERIAL_SAVE_WITH_NAME(mRotDragCoefs);
}

void quadrotor_system::load(ReaK::serialization::iarchive& A,
                            unsigned int /*unused*/) {
  named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_LOAD_WITH_NAME(mMass) &
      RK_SERIAL_LOAD_WITH_NAME(mInertiaMoment) &
      RK_SERIAL_LOAD_WITH_NAME(mTransDragCoefs) &
      RK_SERIAL_LOAD_WITH_NAME(mRotDragCoefs);
  if ((mInertiaMoment.get_row_count() != 3) ||
      (mMass < std::numeric_limits<double>::epsilon())) {
    throw system_incoherency(
        "Inertial information is improper in quadrotor_system's definition");
  }
  try {
    invert_Cholesky(mInertiaMoment, mInertiaMomentInv);
  } catch (singularity_error&) {
    throw system_incoherency(
        "Inertial tensor is singular in quadrotor_system's definition");
  }
}

}  // namespace ReaK::ctrl
