
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

#include <iostream>

using std::size_t;
using std::ptrdiff_t;


#include "topologies/time_topology.hpp"
#include "topologies/time_poisson_topology.hpp"
#include "topologies/line_topology.hpp"
#include "topologies/differentiable_space.hpp"
#include "topologies/temporal_space.hpp"

#include "linear_interp.hpp"
//#include "cubic_hermite_interp.hpp"
//#include "quintic_hermite_interp.hpp"

#include "recorders/ssv_recorder.hpp"


int main(int argc, char** argv) {

  if(argc < 5) {

    std::cout << "Error: Arguments to the program were incorrect!" << std::endl
              << "Usage:" << std::endl
              << "\t\t./test_interp [time_step] [interp_time_step] [max_time] [amplitude]" << std::endl;

    return 1;
  };
  
  double time_step;
  { std::stringstream ss(argv[1]);
    ss >> time_step; };
    
  double interp_time_step;
  { std::stringstream ss(argv[2]);
    ss >> interp_time_step; };
    
  double max_time;
  { std::stringstream ss(argv[3]);
    ss >> max_time; };
    
  double amplitude;
  { std::stringstream ss(argv[4]);
    ss >> amplitude; };
    
  typedef ReaK::arithmetic_tuple< 
            ReaK::pp::line_segment_topology<double>, 
            ReaK::pp::line_segment_topology<double>, 
            ReaK::pp::line_segment_topology<double>, 
            ReaK::pp::line_segment_topology<double> 
          > SpaceTupleType;

  typedef ReaK::pp::differentiable_space< ReaK::pp::time_poisson_topology, SpaceTupleType > TopoType;
  typedef ReaK::pp::metric_topology_traits<TopoType>::point_type PointType;
  typedef ReaK::pp::temporal_space< TopoType, ReaK::pp::time_poisson_topology> TempTopoType;
  typedef ReaK::pp::metric_topology_traits<TempTopoType>::point_type TempPointType;
    
  TempTopoType topo( "temporal_space",
    SpaceTupleType(ReaK::pp::line_segment_topology<double>("pos_topo",-2.0 * amplitude, 2.0 * amplitude),
                   ReaK::pp::line_segment_topology<double>("vel_topo",-2.0 * amplitude, 2.0 * amplitude),
                   ReaK::pp::line_segment_topology<double>("acc_topo",-2.0 * amplitude, 2.0 * amplitude),
                   ReaK::pp::line_segment_topology<double>("jerk_topo",-2.0 * amplitude, 2.0 * amplitude)));
  
  std::vector< TempPointType > pts;
  
  for(double t = 0.0; t <= max_time; t += interp_time_step) {
    pts.push_back( TempPointType(t,PointType(amplitude * std::sin(t), amplitude * std::cos(t), -amplitude * std::sin(t), -amplitude * std::cos(t))) );
  };
    
  try {
    
    ReaK::recorder::ssv_recorder output_rec("test_interp_results/linear_interp.ssv");
    output_rec << "time" << "pos" << "vel" << "acc" << "jerk" << ReaK::recorder::data_recorder::end_name_row;
    
    ReaK::pp::linear_interp<TempTopoType> interp(pts.begin(), pts.end(), topo);
    
    for(double t = 0.0; t <= max_time; t += time_step) {
      TempPointType p = interp.get_point_at_time(t);
      output_rec << p.time << std::get<0>(p.pt) << std::get<1>(p.pt) << std::get<2>(p.pt) << std::get<3>(p.pt) << ReaK::recorder::data_recorder::end_value_row;
    };
    output_rec << ReaK::recorder::data_recorder::flush;
    
  } catch (...) {
    std::cout << "Error: An exception was thrown during the execution of this test!" << std::endl;
    return 1;
  };

  return 0;
};















