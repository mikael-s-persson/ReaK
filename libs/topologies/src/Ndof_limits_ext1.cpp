
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

#include <ReaK/topologies/spaces/Ndof_limits.hpp>

namespace ReaK::pp {

template struct Ndof_limits<double>;

template struct Ndof_limits<double, 1>;
template struct Ndof_limits<double, 2>;
template struct Ndof_limits<double, 3>;
template struct Ndof_limits<double, 4>;
template struct Ndof_limits<double, 5>;
template struct Ndof_limits<double, 6>;
template struct Ndof_limits<double, 7>;
template struct Ndof_limits<double, 8>;
template struct Ndof_limits<double, 9>;
template struct Ndof_limits<double, 10>;

}  // namespace ReaK::pp
