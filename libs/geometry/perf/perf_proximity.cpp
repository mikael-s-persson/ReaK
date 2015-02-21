
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

#include <ReaK/core/base/global_rng.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>
#include <ReaK/geometry/proximity/proximity_finder_3D.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>

#include <ReaK/geometry/shapes/plane.hpp>
#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/sphere.hpp>
#include <ReaK/geometry/shapes/cylinder.hpp>
#include <ReaK/geometry/shapes/capped_cylinder.hpp>

#include <boost/random/uniform_real_distribution.hpp>

#include <ReaK/core/base/chrono_incl.hpp>

#include <iostream>


using namespace ReaK;

struct proxy_query_generator {
  shared_ptr< geom::proxy_query_model_3D > query;
  shared_ptr< pose_3D<double> > anchor;
  
  static pose_3D<double> random_pose() {
    boost::random::uniform_real_distribution<> ud(-1.0, 1.0);
    
    pose_3D<double> rel_pose;
    rel_pose.Position[0] = ud(get_global_rng()) * 5.0;
    rel_pose.Position[1] = ud(get_global_rng()) * 5.0;
    rel_pose.Position[2] = ud(get_global_rng()) * 5.0;
    rel_pose.Quat = quaternion<double>(vect<double,4>(ud(get_global_rng()),
                                                      ud(get_global_rng()),
                                                      ud(get_global_rng()),
                                                      ud(get_global_rng())));
    
    return rel_pose;
  };
  
  proxy_query_generator() : 
    query(new geom::proxy_query_model_3D()), 
    anchor(new pose_3D<double>(random_pose())) { };
  
  void addOneRandomShape() const {
    using std::fabs;
    boost::random::uniform_real_distribution<> ud(0.0, 1.0);
    
    pose_3D<double> rel_pose = random_pose();
    
    double chosen_shape = ud(get_global_rng());
    if( chosen_shape < 0.2 ) {
      // create plane:
      query->addShape(shared_ptr< geom::plane >(new geom::plane(
        "plane", anchor, rel_pose, 
        vect<double,2>(ud(get_global_rng()) * 2.0, ud(get_global_rng()) * 2.0))));
    } else if( chosen_shape < 0.4 ) {
      // create box:
      query->addShape(shared_ptr< geom::box >(new geom::box(
        "box", anchor, rel_pose, 
        vect<double,3>(ud(get_global_rng()) * 2.0, 
                       ud(get_global_rng()) * 2.0,
                       ud(get_global_rng()) * 2.0))));
    } else if( chosen_shape < 0.6 ) {
      // create sphere:
      query->addShape(shared_ptr< geom::sphere >(new geom::sphere(
        "sphere", anchor, rel_pose, ud(get_global_rng()) * 2.0)));
    } else if( chosen_shape < 0.8 ) {
      // create cylinder:
      query->addShape(shared_ptr< geom::cylinder >(new geom::cylinder(
        "cylinder", anchor, rel_pose, 
        ud(get_global_rng()) * 2.0, 
        ud(get_global_rng()) * 1.0)));
    } else {
      // create capped_cylinder:
      query->addShape(shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder(
        "capped_cylinder", anchor, rel_pose, 
        ud(get_global_rng()) * 2.0, 
        ud(get_global_rng()) * 1.0)));
    };
  };
  
};


int main(int argc, const char* argv[]) {
  
  using namespace ReaKaux::chrono;
  
  const int num_runs   = ((argc > 1) ? std::atoi(argv[1]) : 100);
  const int num_passes = ((argc > 2) ? std::atoi(argv[2]) : 100);
  const int num_shapes = ((argc > 3) ? std::atoi(argv[3]) : 50);
  
  high_resolution_clock::duration accum_dt = high_resolution_clock::duration::zero();
  
  for(int i = 0; i < num_runs; ++i) {
    proxy_query_generator pg1, pg2;
    for(int j = 0; j < num_shapes; ++j) {
      pg1.addOneRandomShape();
      pg2.addOneRandomShape();
    };
    
    geom::proxy_query_pair_3D pqp("", pg1.query, pg2.query);
    
    std::vector< geom::proximity_record_3D > col_pts;
    col_pts.reserve(num_shapes * num_shapes);
    
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int j = 0; j < num_passes; ++j) {
      col_pts.clear();
      pqp.gatherCollisionPoints(col_pts);
    };
    accum_dt += high_resolution_clock::now() - t1;
    
  };
  
  std::cout << (duration_cast<nanoseconds>(accum_dt).count() / (num_runs * num_passes)) << std::endl;
  
  return 0;
};








