
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

#include "kte_map_chain.hpp"

#include "inertia.hpp"
#include "revolute_joint.hpp"
#include "rigid_link.hpp"
#include "jacobian_joint_map.hpp"
#include "mass_matrix_calculator.hpp"

#include "recorders/ssv_recorder.hpp"

#include "serialization/xml_archiver.hpp"


using namespace ReaK;

using namespace ReaK::kte;

using namespace serialization;

int main() {
  
#if 1
  shared_ptr<frame_2D<double> > base_frame = rtti::rk_dynamic_ptr_cast< frame_2D<double> >(frame_2D<double>::Create());
  shared_ptr<frame_2D<double> > joint_frame = rtti::rk_dynamic_ptr_cast< frame_2D<double> >(frame_2D<double>::Create());
  shared_ptr<jacobian_gen_2D<double> > joint_jacobian = rtti::rk_dynamic_ptr_cast< jacobian_gen_2D<double> >(jacobian_gen_2D<double>::Create());
  shared_ptr<frame_2D<double> > end_frame = rtti::rk_dynamic_ptr_cast< frame_2D<double> >(frame_2D<double>::Create());
  shared_ptr<gen_coord<double> > joint_coord = rtti::rk_dynamic_ptr_cast< gen_coord<double> >(gen_coord<double>::Create());

  base_frame->Acceleration = vect<double,2>(0,9.81); //add gravity

  //create revolute joint
  shared_ptr<revolute_joint_2D> rev_joint(new revolute_joint_2D("joint1",joint_coord,base_frame,joint_frame,joint_jacobian),scoped_deleter());
  //create link of lenght 0.5 meters
  shared_ptr<rigid_link_2D> link1(new rigid_link_2D("link1",joint_frame,end_frame,pose_2D<double>(weak_ptr<pose_2D<double> >(),vect<double,2>(0.5,0.0),rot_mat_2D<double>(0.0))),scoped_deleter());

  jacobian_joint_map_2D joint_map1;
  joint_map1[joint_coord] = joint_jacobian;
  //create end mass of 1.0 kg (point mass only)
  shared_ptr<inertia_2D> mass1(new inertia_2D("mass1",
                                                        shared_ptr<joint_dependent_frame_2D>(new joint_dependent_frame_2D(end_frame,joint_map1),scoped_deleter()),
                                                        1.0,0.0),scoped_deleter());

  shared_ptr<mass_matrix_calc> mass_mat1(new mass_matrix_calc("mass_mat1"),scoped_deleter());
  
  *mass_mat1 << mass1;
  *mass_mat1 << joint_coord;
  
  kte_map_chain pendulum("pendulum");
  pendulum << rev_joint << link1 << mass1;
  {
    xml_oarchive pendulum_arc("pendulum.xml");
    oarchive& arc_ref = pendulum_arc;
    arc_ref << joint_coord << pendulum << end_frame << mass_mat1;

  };
#else
  shared_ptr< gen_coord<double> > joint_coord;
  shared_ptr< frame_2D<double> > end_frame;
  shared_ptr< mass_matrix_calc > mass_mat1;

  kte_map_chain pendulum;

  {
    xml_iarchive pendulum_arc("pendulum.xml");
    iarchive& arc_ref = pendulum_arc;
    arc_ref >> joint_coord >> pendulum >> end_frame >> mass_mat1;
  };

  std::cout << "Pendulum model loaded.. staring simulation!" << std::endl;

#endif

  recorder::ssv_recorder output_rec("pendulum_results.ssvdat");
  output_rec << "time" << "q" << "qd" << "qdd" << "f" << recorder::data_recorder::end_name_row;

  double sim_time = 0.0;

  for (;sim_time < 50.0;sim_time += 0.001) {

    joint_coord->q_ddot = 0.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl = joint_coord->f;

    mat<double,mat_structure::symmetric> M;
    mass_mat1->getMassMatrix(M);

    joint_coord->q_ddot = 1.0;
    pendulum.doMotion();
    pendulum.clearForce();
    pendulum.doForce();
    double f_nl_in = joint_coord->f;
    
    
    //std::cout << (f_nl - f_nl_in) << " " << M(0,0) << std::endl;

    joint_coord->q_ddot = f_nl / (f_nl - f_nl_in);

    output_rec << sim_time
               << joint_coord->q
	       << joint_coord->q_dot
	       << joint_coord->q_ddot
	       << f_nl << recorder::data_recorder::end_value_row;

    joint_coord->q += joint_coord->q_dot * 0.001;
    joint_coord->q_dot += joint_coord->q_ddot * 0.001;

  };

  output_rec << recorder::data_recorder::close;

};






