
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
// #ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/math/lin_alg/vect_alg.hpp>

namespace ReaK {

template class vect<double,2>;
template class vect<double,3>;
template class vect<double,4>;
template class vect<double,6>;
template class vect_n<double>;

template class vect<float,2>;
template class vect<float,3>;
template class vect<float,4>;
template class vect<float,6>;
template class vect_n<float>;


template vect<double,2> diff(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3> diff(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4> diff(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6> diff(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double> diff(const vect_n<double>& v1, const vect_n<double>& v2);

template vect<double,2> add(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3> add(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4> add(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6> add(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double> add(const vect_n<double>& v1, const vect_n<double>& v2);

template double operator %(const vect<double,2>& V1,const vect<double,2>& V2) BOOST_NOEXCEPT;
template vect<double,2> operator %(const double& S, const vect<double,2>& V) BOOST_NOEXCEPT;
template vect<double,2> operator %(const vect<double,2>& V,const double& S) BOOST_NOEXCEPT;
template vect<double,3> operator %(const vect<double,3>& V1 , const vect<double,3>& V2) BOOST_NOEXCEPT;


template vect<double,2>& operator +=(vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3>& operator +=(vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4>& operator +=(vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6>& operator +=(vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double>& operator +=(vect_n<double>& v1, const vect_n<double>& v2);

template vect<double,2>& operator -=(vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3>& operator -=(vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4>& operator -=(vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6>& operator -=(vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double>& operator -=(vect_n<double>& v1, const vect_n<double>& v2);

template vect<double,2>& operator *=(vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,3>& operator *=(vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,4>& operator *=(vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,6>& operator *=(vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
template vect_n<double>& operator *=(vect_n<double>& v1, const double& v2) BOOST_NOEXCEPT;

template vect<double,2>& operator /=(vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,3>& operator /=(vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,4>& operator /=(vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,6>& operator /=(vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
template vect_n<double>& operator /=(vect_n<double>& v1, const double& v2) BOOST_NOEXCEPT;

template vect<double,2> operator +(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3> operator +(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4> operator +(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6> operator +(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double> operator +(const vect_n<double>& v1, const vect_n<double>& v2);

template vect<double,2> operator -(const vect<double,2>& v1) BOOST_NOEXCEPT;
template vect<double,3> operator -(const vect<double,3>& v1) BOOST_NOEXCEPT;
template vect<double,4> operator -(const vect<double,4>& v1) BOOST_NOEXCEPT;
template vect<double,6> operator -(const vect<double,6>& v1) BOOST_NOEXCEPT;
template vect_n<double> operator -(const vect_n<double>& v1);

template vect<double,2> operator -(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3> operator -(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4> operator -(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6> operator -(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double> operator -(const vect_n<double>& v1, const vect_n<double>& v2);

template double operator *(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template double operator *(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template double operator *(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template double operator *(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template double operator *(const vect_n<double>& v1, const vect_n<double>& v2);

template vect<double,2> operator *(const vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,3> operator *(const vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,4> operator *(const vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,6> operator *(const vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
template vect_n<double> operator *(const vect_n<double>& v1, const double& v2);

template vect<double,2> operator *(const double& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template vect<double,3> operator *(const double& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template vect<double,4> operator *(const double& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template vect<double,6> operator *(const double& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template vect_n<double> operator *(const double& v1, const vect_n<double>& v2);

template vect<double,2> operator /(const vect<double,2>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,3> operator /(const vect<double,3>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,4> operator /(const vect<double,4>& v1, const double& v2) BOOST_NOEXCEPT;
template vect<double,6> operator /(const vect<double,6>& v1, const double& v2) BOOST_NOEXCEPT;
template vect_n<double> operator /(const vect_n<double>& v1, const double& v2);


template bool operator ==(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

template bool operator !=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

template bool operator <=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

template bool operator >=(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

template bool operator <(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;

template bool operator >(const vect<double,2>& v1, const vect<double,2>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<double,3>& v1, const vect<double,3>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<double,4>& v1, const vect<double,4>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<double,6>& v1, const vect<double,6>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect_n<double>& v1, const vect_n<double>& v2) BOOST_NOEXCEPT;


template std::ostream& operator <<(std::ostream& out_stream, const vect<double,2>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<double,3>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<double,4>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<double,6>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect_n<double>& V);



template vect<float,2> diff(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3> diff(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4> diff(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6> diff(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float> diff(const vect_n<float>& v1, const vect_n<float>& v2);

template vect<float,2> add(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3> add(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4> add(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6> add(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float> add(const vect_n<float>& v1, const vect_n<float>& v2);

template float operator %(const vect<float,2>& V1,const vect<float,2>& V2) BOOST_NOEXCEPT;
template vect<float,2> operator %(const float& S, const vect<float,2>& V) BOOST_NOEXCEPT;
template vect<float,2> operator %(const vect<float,2>& V,const float& S) BOOST_NOEXCEPT;
template vect<float,3> operator %(const vect<float,3>& V1 , const vect<float,3>& V2) BOOST_NOEXCEPT;


template vect<float,2>& operator +=(vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3>& operator +=(vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4>& operator +=(vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6>& operator +=(vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float>& operator +=(vect_n<float>& v1, const vect_n<float>& v2);

template vect<float,2>& operator -=(vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3>& operator -=(vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4>& operator -=(vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6>& operator -=(vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float>& operator -=(vect_n<float>& v1, const vect_n<float>& v2);

template vect<float,2>& operator *=(vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,3>& operator *=(vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,4>& operator *=(vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,6>& operator *=(vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
template vect_n<float>& operator *=(vect_n<float>& v1, const float& v2) BOOST_NOEXCEPT;

template vect<float,2>& operator /=(vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,3>& operator /=(vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,4>& operator /=(vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,6>& operator /=(vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
template vect_n<float>& operator /=(vect_n<float>& v1, const float& v2) BOOST_NOEXCEPT;

template vect<float,2> operator +(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3> operator +(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4> operator +(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6> operator +(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float> operator +(const vect_n<float>& v1, const vect_n<float>& v2);

template vect<float,2> operator -(const vect<float,2>& v1) BOOST_NOEXCEPT;
template vect<float,3> operator -(const vect<float,3>& v1) BOOST_NOEXCEPT;
template vect<float,4> operator -(const vect<float,4>& v1) BOOST_NOEXCEPT;
template vect<float,6> operator -(const vect<float,6>& v1) BOOST_NOEXCEPT;
template vect_n<float> operator -(const vect_n<float>& v1);

template vect<float,2> operator -(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3> operator -(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4> operator -(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6> operator -(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float> operator -(const vect_n<float>& v1, const vect_n<float>& v2);

template float operator *(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template float operator *(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template float operator *(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template float operator *(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template float operator *(const vect_n<float>& v1, const vect_n<float>& v2);

template vect<float,2> operator *(const vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,3> operator *(const vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,4> operator *(const vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,6> operator *(const vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
template vect_n<float> operator *(const vect_n<float>& v1, const float& v2);

template vect<float,2> operator *(const float& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template vect<float,3> operator *(const float& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template vect<float,4> operator *(const float& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template vect<float,6> operator *(const float& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template vect_n<float> operator *(const float& v1, const vect_n<float>& v2);

template vect<float,2> operator /(const vect<float,2>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,3> operator /(const vect<float,3>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,4> operator /(const vect<float,4>& v1, const float& v2) BOOST_NOEXCEPT;
template vect<float,6> operator /(const vect<float,6>& v1, const float& v2) BOOST_NOEXCEPT;
template vect_n<float> operator /(const vect_n<float>& v1, const float& v2);


template bool operator ==(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator ==(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

template bool operator !=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator !=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

template bool operator <=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator <=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

template bool operator >=(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator >=(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

template bool operator <(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator <(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;

template bool operator >(const vect<float,2>& v1, const vect<float,2>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<float,3>& v1, const vect<float,3>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<float,4>& v1, const vect<float,4>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect<float,6>& v1, const vect<float,6>& v2) BOOST_NOEXCEPT;
template bool operator >(const vect_n<float>& v1, const vect_n<float>& v2) BOOST_NOEXCEPT;


template std::ostream& operator <<(std::ostream& out_stream, const vect<float,2>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<float,3>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<float,4>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect<float,6>& V);
template std::ostream& operator <<(std::ostream& out_stream, const vect_n<float>& V);





};

#else

namespace ReaK {

void dummy_vect_alg_externs_symbol() { };

};


#endif














