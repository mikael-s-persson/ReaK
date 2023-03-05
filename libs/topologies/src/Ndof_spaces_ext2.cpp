
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

#include "ReaK/core/base/defs.hpp"

#include "ReaK/topologies/spaces/Ndof_spaces.hpp"

namespace ReaK::pp {

#define RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(NDOF)                    \
  template class differentiable_space<                                       \
      time_topology, arithmetic_tuple<hyperbox_topology<                     \
                         vect<double, NDOF>, inf_norm_distance_metric>>>;    \
  template class differentiable_space<                                       \
      time_topology,                                                         \
      arithmetic_tuple<                                                      \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>>>; \
  template class differentiable_space<                                       \
      time_topology,                                                         \
      arithmetic_tuple<                                                      \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>,   \
          hyperbox_topology<vect<double, NDOF>, inf_norm_distance_metric>>>;

RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(1)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(2)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(3)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(4)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(5)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(6)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(7)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(8)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(9)
RK_NDOF_SPACES_MAKE_RATELIMITED_EXTERN_DEFS(10)

template class differentiable_space<
    time_topology, arithmetic_tuple<hyperbox_topology<
                       vect_n<double>, inf_norm_distance_metric>>>;
template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>>>;
template class differentiable_space<
    time_topology,
    arithmetic_tuple<
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>,
        hyperbox_topology<vect_n<double>, inf_norm_distance_metric>>>;

}  // namespace ReaK::pp
