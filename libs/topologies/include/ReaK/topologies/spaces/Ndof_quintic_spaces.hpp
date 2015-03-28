/**
 * \file Ndof_quintic_spaces.hpp
 *
 * This library provides classes to represent N-dof cubic-interpolated spaces of either static or
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

#ifndef REAK_NDOF_QUINTIC_SPACES_HPP
#define REAK_NDOF_QUINTIC_SPACES_HPP

#include "Ndof_spaces.hpp"
#include "temporal_space.hpp"
#include "time_poisson_topology.hpp"
#include "reachability_space.hpp"
#include <ReaK/topologies/interpolation/interpolated_topologies.hpp>
#include <ReaK/topologies/interpolation/quintic_hermite_interp.hpp>


#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

namespace ReaK {

namespace pp {

#define RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( NDOF )                                                  \
  extern template class interpolated_topology< Ndof_rl_space< double, NDOF, 2 >::type,                          \
                                               quintic_hermite_interpolation_tag >;                             \
                                                                                                                \
  extern template class interpolated_topology< temporal_space< Ndof_rl_space< double, NDOF, 2 >::type,          \
                                                               time_poisson_topology, reach_plus_time_metric >, \
                                               quintic_hermite_interpolation_tag >;

RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 1 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 2 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 3 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 4 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 5 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 6 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 7 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 8 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 9 )
RK_NDOF_QUINTIC_SPACES_MAKE_NORMAL_EXTERN_DECL( 10 )

extern template class interpolated_topology< Ndof_rl_space< double, 0, 2 >::type, quintic_hermite_interpolation_tag >;

extern template class interpolated_topology< temporal_space< Ndof_rl_space< double, 0, 2 >::type, time_poisson_topology,
                                                             reach_plus_time_metric >,
                                             quintic_hermite_interpolation_tag >;
};
};

#endif

#endif
