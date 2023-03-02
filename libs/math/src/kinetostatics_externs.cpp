
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

#if 0

#include <ReaK/math/kinetostatics/kinetostatics.hpp>
#include <ReaK/math/kinetostatics/motion_jacobians.hpp>
#include <ReaK/math/kinetostatics/quat_alg.hpp>
#include <ReaK/math/kinetostatics/rotations.hpp>

namespace ReaK {

#if 0
template class rot_mat_2D<double>;
template class trans_mat_2D<double>;

template std::ostream& operator <<(std::ostream& out_stream, const rot_mat_2D<double>& R);
template std::ostream& operator <<(std::ostream& out_stream, const trans_mat_2D<double>& M);


template class rot_mat_2D<float>;
template class trans_mat_2D<float>;

template std::ostream& operator <<(std::ostream& out_stream, const rot_mat_2D<float>& R);
template std::ostream& operator <<(std::ostream& out_stream, const trans_mat_2D<float>& M);



template class rot_mat_3D<double>;
template class quaternion<double>;
template class euler_angles_TB<double>;
template class axis_angle<double>;
template class trans_mat_3D<double>;

template std::ostream& operator <<(std::ostream& out_stream, const rot_mat_3D<double>& R);
template std::ostream& operator <<(std::ostream& out_stream, const quaternion<double>& Q);
template std::ostream& operator <<(std::ostream& out_stream, const euler_angles_TB<double>& E);
template std::ostream& operator <<(std::ostream& out_stream, const axis_angle<double>& A);
template std::ostream& operator <<(std::ostream& out_stream, const trans_mat_3D<double>& M);


template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const quaternion<double>& Q);
template rot_mat_3D<double> operator *(const quaternion<double>& Q, const rot_mat_3D<double>& R);
template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
template rot_mat_3D<double> operator *(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
template quaternion<double> operator *(const quaternion<double>& Q, const euler_angles_TB<double>& E);
template quaternion<double> operator *(const euler_angles_TB<double>& E, const quaternion<double>& Q);
template rot_mat_3D<double> operator *(const rot_mat_3D<double>& R, const axis_angle<double>& A);
template rot_mat_3D<double> operator *(const axis_angle<double>& A, const rot_mat_3D<double>& R);
template quaternion<double> operator *(const quaternion<double>& Q, const axis_angle<double>& A);
template quaternion<double> operator *(const axis_angle<double>& A, const quaternion<double>& Q);
template rot_mat_3D<double> operator *(const euler_angles_TB<double>& E, const axis_angle<double>& A);
template rot_mat_3D<double> operator *(const axis_angle<double>& A, const euler_angles_TB<double>& E);
template trans_mat_3D<double> operator *(const rot_mat_3D<double>& R, const trans_mat_3D<double>& M);
template trans_mat_3D<double> operator *(const quaternion<double>& Q, const trans_mat_3D<double>& M);
template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const quaternion<double>& Q);
template trans_mat_3D<double> operator *(const euler_angles_TB<double>& E, const trans_mat_3D<double>& M);
template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const euler_angles_TB<double>& E);
template trans_mat_3D<double> operator *(const axis_angle<double>& A, const trans_mat_3D<double>& M);
template trans_mat_3D<double> operator *(const trans_mat_3D<double>& M, const axis_angle<double>& A);

template bool operator ==(const rot_mat_3D<double>& R, const quaternion<double>& Q);
template bool operator !=(const rot_mat_3D<double>& R, const quaternion<double>& Q);
template bool operator ==(const quaternion<double>& Q, const rot_mat_3D<double>& R);
template bool operator !=(const quaternion<double>& Q, const rot_mat_3D<double>& R);
template bool operator ==(const quaternion<double>& Q, const euler_angles_TB<double>& E);
template bool operator !=(const quaternion<double>& Q, const euler_angles_TB<double>& E);
template bool operator ==(const euler_angles_TB<double>& E, const quaternion<double>& Q);
template bool operator !=(const euler_angles_TB<double>& E, const quaternion<double>& Q);
template bool operator ==(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
template bool operator !=(const rot_mat_3D<double>& R, const euler_angles_TB<double>& E);
template bool operator ==(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
template bool operator !=(const euler_angles_TB<double>& E, const rot_mat_3D<double>& R);
template bool operator ==(const rot_mat_3D<double>& R, const axis_angle<double>& A);
template bool operator !=(const rot_mat_3D<double>& R, const axis_angle<double>& A);
template bool operator ==(const axis_angle<double>& A, const rot_mat_3D<double>& R);
template bool operator !=(const axis_angle<double>& A, const rot_mat_3D<double>& R);
template bool operator ==(const quaternion<double>& Q, const axis_angle<double>& A);
template bool operator !=(const quaternion<double>& Q, const axis_angle<double>& A);
template bool operator ==(const axis_angle<double>& A, const quaternion<double>& Q);
template bool operator !=(const axis_angle<double>& A, const quaternion<double>& Q);
template bool operator ==(const euler_angles_TB<double>& E, const axis_angle<double>& A);
template bool operator !=(const euler_angles_TB<double>& E, const axis_angle<double>& A);
template bool operator ==(const axis_angle<double>& A, const euler_angles_TB<double>& E);
template bool operator !=(const axis_angle<double>& A, const euler_angles_TB<double>& E);



template class rot_mat_3D<float>;
template class quaternion<float>;
template class euler_angles_TB<float>;
template class axis_angle<float>;
template class trans_mat_3D<float>;

template std::ostream& operator <<(std::ostream& out_stream, const rot_mat_3D<float>& R);
template std::ostream& operator <<(std::ostream& out_stream, const quaternion<float>& Q);
template std::ostream& operator <<(std::ostream& out_stream, const euler_angles_TB<float>& E);
template std::ostream& operator <<(std::ostream& out_stream, const axis_angle<float>& A);
template std::ostream& operator <<(std::ostream& out_stream, const trans_mat_3D<float>& M);


template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const quaternion<float>& Q);
template rot_mat_3D<float> operator *(const quaternion<float>& Q, const rot_mat_3D<float>& R);
template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
template rot_mat_3D<float> operator *(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
template quaternion<float> operator *(const quaternion<float>& Q, const euler_angles_TB<float>& E);
template quaternion<float> operator *(const euler_angles_TB<float>& E, const quaternion<float>& Q);
template rot_mat_3D<float> operator *(const rot_mat_3D<float>& R, const axis_angle<float>& A);
template rot_mat_3D<float> operator *(const axis_angle<float>& A, const rot_mat_3D<float>& R);
template quaternion<float> operator *(const quaternion<float>& Q, const axis_angle<float>& A);
template quaternion<float> operator *(const axis_angle<float>& A, const quaternion<float>& Q);
template rot_mat_3D<float> operator *(const euler_angles_TB<float>& E, const axis_angle<float>& A);
template rot_mat_3D<float> operator *(const axis_angle<float>& A, const euler_angles_TB<float>& E);
template trans_mat_3D<float> operator *(const rot_mat_3D<float>& R, const trans_mat_3D<float>& M);
template trans_mat_3D<float> operator *(const quaternion<float>& Q, const trans_mat_3D<float>& M);
template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const quaternion<float>& Q);
template trans_mat_3D<float> operator *(const euler_angles_TB<float>& E, const trans_mat_3D<float>& M);
template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const euler_angles_TB<float>& E);
template trans_mat_3D<float> operator *(const axis_angle<float>& A, const trans_mat_3D<float>& M);
template trans_mat_3D<float> operator *(const trans_mat_3D<float>& M, const axis_angle<float>& A);

template bool operator ==(const rot_mat_3D<float>& R, const quaternion<float>& Q);
template bool operator !=(const rot_mat_3D<float>& R, const quaternion<float>& Q);
template bool operator ==(const quaternion<float>& Q, const rot_mat_3D<float>& R);
template bool operator !=(const quaternion<float>& Q, const rot_mat_3D<float>& R);
template bool operator ==(const quaternion<float>& Q, const euler_angles_TB<float>& E);
template bool operator !=(const quaternion<float>& Q, const euler_angles_TB<float>& E);
template bool operator ==(const euler_angles_TB<float>& E, const quaternion<float>& Q);
template bool operator !=(const euler_angles_TB<float>& E, const quaternion<float>& Q);
template bool operator ==(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
template bool operator !=(const rot_mat_3D<float>& R, const euler_angles_TB<float>& E);
template bool operator ==(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
template bool operator !=(const euler_angles_TB<float>& E, const rot_mat_3D<float>& R);
template bool operator ==(const rot_mat_3D<float>& R, const axis_angle<float>& A);
template bool operator !=(const rot_mat_3D<float>& R, const axis_angle<float>& A);
template bool operator ==(const axis_angle<float>& A, const rot_mat_3D<float>& R);
template bool operator !=(const axis_angle<float>& A, const rot_mat_3D<float>& R);
template bool operator ==(const quaternion<float>& Q, const axis_angle<float>& A);
template bool operator !=(const quaternion<float>& Q, const axis_angle<float>& A);
template bool operator ==(const axis_angle<float>& A, const quaternion<float>& Q);
template bool operator !=(const axis_angle<float>& A, const quaternion<float>& Q);
template bool operator ==(const euler_angles_TB<float>& E, const axis_angle<float>& A);
template bool operator !=(const euler_angles_TB<float>& E, const axis_angle<float>& A);
template bool operator ==(const axis_angle<float>& A, const euler_angles_TB<float>& E);
template bool operator !=(const axis_angle<float>& A, const euler_angles_TB<float>& E);

#endif


template class quat< double >;
template class unit_quat< double >;

template class quat< float >;
template class unit_quat< float >;


template class gen_coord< double >;
template std::ostream& operator<<( std::ostream& out, const gen_coord< double >& g );

template class gen_coord< float >;
template std::ostream& operator<<( std::ostream& out, const gen_coord< float >& g );


template class pose_2D< double >;
template std::ostream& operator<<( std::ostream& out, const pose_2D< double >& g );

template class pose_2D< float >;
template std::ostream& operator<<( std::ostream& out, const pose_2D< float >& g );


template class pose_3D< double >;
template std::ostream& operator<<( std::ostream& out, const pose_3D< double >& g );

template class pose_3D< float >;
template std::ostream& operator<<( std::ostream& out, const pose_3D< float >& g );


template class frame_2D< double >;
template std::ostream& operator<<( std::ostream& out, const frame_2D< double >& g );

template class frame_2D< float >;
template std::ostream& operator<<( std::ostream& out, const frame_2D< float >& g );


template class frame_3D< double >;
template std::ostream& operator<<( std::ostream& out, const frame_3D< double >& g );

template class frame_3D< float >;
template std::ostream& operator<<( std::ostream& out, const frame_3D< float >& g );


template class jacobian_gen_gen< double >;
template class jacobian_gen_2D< double >;
template class jacobian_gen_3D< double >;
template class jacobian_2D_gen< double >;
template class jacobian_2D_2D< double >;
template class jacobian_2D_3D< double >;
template class jacobian_3D_gen< double >;
template class jacobian_3D_2D< double >;
template class jacobian_3D_3D< double >;

template class jacobian_gen_gen< float >;
template class jacobian_gen_2D< float >;
template class jacobian_gen_3D< float >;
template class jacobian_2D_gen< float >;
template class jacobian_2D_2D< float >;
template class jacobian_2D_3D< float >;
template class jacobian_3D_gen< float >;
template class jacobian_3D_2D< float >;
template class jacobian_3D_3D< float >;
};

#else

namespace ReaK {

void dummy_kinetostatics_externs_symbol(){};
};

#endif
