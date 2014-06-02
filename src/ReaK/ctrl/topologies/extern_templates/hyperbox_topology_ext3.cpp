
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

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/ctrl/topologies/hyperbox_topology.hpp>

namespace ReaK {

namespace pp {

template class hyperbox_topology< vect<double,1>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,2>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,3>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,4>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,5>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,6>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,7>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,8>, inf_norm_distance_metric >;
template class hyperbox_topology< vect<double,9>, inf_norm_distance_metric >;
template class hyperbox_topology< vect_n<double>, inf_norm_distance_metric >;

};

};

#else

namespace ReaK {

namespace pp {

void dummy_hyperbox_topology_externs_3_symbol() { };

};

};

#endif














