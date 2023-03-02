
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

#include <ReaK/topologies/spaces/Ndof_sap_spaces.hpp>

namespace ReaK::pp {

#define RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(NDOF)                      \
  template class interpolated_topology<Ndof_rl_space_t<double, NDOF, 2>,      \
                                       sap_Ndof_interpolation_tag>;           \
                                                                              \
  template class interpolated_topology<                                       \
      temporal_space<Ndof_rl_space_t<double, NDOF, 2>, time_poisson_topology, \
                     reach_plus_time_metric>,                                 \
      sap_Ndof_interpolation_tag>;

RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(1)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(2)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(3)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(4)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(5)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(6)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(7)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(8)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(9)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DEFS(10)

template class interpolated_topology<Ndof_rl_space_t<double, 0, 2>,
                                     sap_Ndof_interpolation_tag>;

template class interpolated_topology<
    temporal_space<Ndof_rl_space_t<double, 0, 2>, time_poisson_topology,
                   reach_plus_time_metric>,
    sap_Ndof_interpolation_tag>;

}  // namespace ReaK::pp
