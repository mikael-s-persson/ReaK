/**
 * \file Ndof_sap_spaces.h
 *
 * This library provides classes to represent N-dof SAP-interpolated spaces of either static or
 * dynamic (run-time) dimensions, and of differentiation order 0, 1 or 2.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
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

#ifndef REAK_TOPOLOGIES_SPACES_NDOF_SAP_SPACES_H_
#define REAK_TOPOLOGIES_SPACES_NDOF_SAP_SPACES_H_

#include "ReaK/topologies/interpolation/sap_Ndof_reach_topologies.h"
#include "ReaK/topologies/spaces/Ndof_spaces.h"
#include "ReaK/topologies/spaces/reachability_space.h"
#include "ReaK/topologies/spaces/temporal_space.h"
#include "ReaK/topologies/spaces/time_poisson_topology.h"

namespace ReaK::pp {

#define RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(NDOF)                      \
  extern template class interpolated_topology<                                \
      Ndof_rl_space_t<double, NDOF, 2>, sap_Ndof_interpolation_tag>;          \
                                                                              \
  extern template class interpolated_topology<                                \
      temporal_space<Ndof_rl_space_t<double, NDOF, 2>, time_poisson_topology, \
                     reach_plus_time_metric>,                                 \
      sap_Ndof_interpolation_tag>;

RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(1)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(2)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(3)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(4)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(5)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(6)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(7)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(8)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(9)
RK_NDOF_SAP_SPACES_MAKE_NORMAL_EXTERN_DECL(10)

extern template class interpolated_topology<Ndof_rl_space_t<double, 0, 2>,
                                            sap_Ndof_interpolation_tag>;

extern template class interpolated_topology<
    temporal_space<Ndof_rl_space_t<double, 0, 2>, time_poisson_topology,
                   reach_plus_time_metric>,
    sap_Ndof_interpolation_tag>;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_NDOF_SAP_SPACES_H_
