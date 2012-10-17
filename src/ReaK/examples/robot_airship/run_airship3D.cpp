
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

#include "mbd_kte/inertia.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "mbd_kte/driving_actuator.hpp"
#include "mbd_kte/state_measures.hpp"
#include "mbd_kte/free_joints.hpp"

#include "mbd_kte/kte_map_chain.hpp"

#include "ctrl_sys/kte_nl_system.hpp"
#include "ctrl_sys/num_int_dtnl_system.hpp"
#include "integrators/variable_step_integrators.hpp"

#include "serialization/xml_archiver.hpp"
#include "recorders/ssv_recorder.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "airship3D_lin_model.hpp"

#include <iomanip>

#include <ctime>

int main(int argc, char** argv) {
  using namespace ReaK;
  
  if(argc < 9) {
    std::cout << "Usage:\n"
	      << "\t./run_airship3D [model_filename.xml] [init_conditions.xml] [inertia_data.xml] [result_filename.ssv] [time_step] [end_time] [Qu.xml] [R.xml]\n"
	      << "\t\t model_filename.xml:\t The filename for the airship model to use.\n"
	      << "\t\t init_conditions.xml:\t The filename for the airship's initial conditions.\n"
	      << "\t\t inertia_data.xml:\t The filename for the airship's inertial data.\n"
	      << "\t\t result_filename.ssv:\t The filename where to record the results as a space-separated values file.\n"
	      << "\t\t time_step:\t\t The time_step of the output points.\n"
	      << "\t\t end_time:\t\t The end-time of the simulation.\n"
	      << "\t\t Qu.xml:\t\t The filename for the airship's input disturbance covariance matrix.\n"
	      << "\t\t R.xml:\t\t The filename for the airship's measurement noise covariance matrix." << std::endl;
    return 0;
  };
  
  std::string model_filename(argv[1]);
  std::string init_filename(argv[2]);
  std::string inertia_filename(argv[3]);
  std::string results_filename(argv[4]);
  
  double time_step = 0.001;
  std::stringstream(argv[5]) >> time_step;
  double end_time = 30.0;
  std::stringstream(argv[6]) >> end_time;
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  std::string Qu_filename(argv[7]);
  std::string R_filename(argv[8]);
  
  shared_ptr< frame_3D<double> > airship3D_frame;
  shared_ptr< kte::position_measure_3D > airship3D_position;
  shared_ptr< kte::rotation_measure_3D > airship3D_rotation;
  
  shared_ptr< kte::driving_actuator_3D > airship3D_actuator;
  shared_ptr< kte::inertia_3D > airship3D_inertia;
  shared_ptr< kte::kte_map_chain > airship3D_model;
  shared_ptr< kte::mass_matrix_calc > airship3D_mass_calc;
  
  
  try {
    serialization::xml_iarchive in(model_filename);
    in >> airship3D_frame
       >> airship3D_position
       >> airship3D_rotation
       >> airship3D_actuator
       >> airship3D_inertia
       >> airship3D_model
       >> airship3D_mass_calc;
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the airship's model!");
    return 1;
  };
  
  /* initial states */
  try {
    serialization::xml_iarchive in(init_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("initial_motion",(*airship3D_frame));
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 2;
  };
  
  /* inertial data */
  try {
    serialization::xml_iarchive in(inertia_filename);
    double tmp_mass;
    ReaK::mat<double,ReaK::mat_structure::symmetric> tmp_tensor;
    in & RK_SERIAL_LOAD_WITH_ALIAS("mass",tmp_mass)
       & RK_SERIAL_LOAD_WITH_ALIAS("inertia_tensor",tmp_tensor);
    airship3D_inertia->setMass(tmp_mass);
    airship3D_inertia->setInertiaTensor(tmp_tensor);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 2;
  };
  
  /* input disturbance */
  mat<double,mat_structure::diagonal> Qu;
  try {
    serialization::xml_iarchive in(Qu_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("input_disturbance",Qu);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return 2;
  };
  
  /* measurement noise */
  mat<double,mat_structure::diagonal> R;
  try {
    serialization::xml_iarchive in(R_filename);
    in & RK_SERIAL_LOAD_WITH_ALIAS("measurement_noise",R);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return 2;
  };
  
  shared_ptr<ctrl::kte_nl_system> airship3D_system(new ctrl::kte_nl_system("airship3D_system"));
  
  airship3D_system->dofs_3D.push_back(airship3D_frame);
  airship3D_system->inputs.push_back(airship3D_actuator);
  airship3D_system->outputs.push_back(airship3D_position);
  airship3D_system->outputs.push_back(airship3D_rotation);
  airship3D_system->chain = airship3D_model;
  airship3D_system->mass_calc = airship3D_mass_calc;
  
  typedef ctrl::num_int_dtnl_sys< ctrl::kte_nl_system, dormand_prince45_integrator<double> > sys_type;
  sys_type 
    airship3D_dt_sys( airship3D_system, 
                      dormand_prince45_integrator<double>("ode45_integrator",
							  ReaK::vect_n<double>(),
							  0.0,
							  1e-4,
							  weak_ptr< state_rate_function<double> >(),
							  time_step,
							  time_step * 0.000001,
							  1e-3),
		      time_step,
		      "airship3D_dt_sys");
    
  ctrl::airship3D_inv_dt_system mdl_inv_dt("airship3D_invariant_discrete",
					   1.0,
					   mat<double,mat_structure::symmetric>(mat<double,mat_structure::identity>(3)),
					   time_step);
  
  pp::vector_topology< vect_n<double> > mdl_state_space;
  
  recorder::ssv_recorder results(results_filename);
  
  results << "time" << "meas_q0" << "meas_q1" << "meas_q2" << "meas_q3"
		    << "meas_x" << "meas_y" << "meas_z"
		    << "q0" << "q1" << "q2" << "q3" 
		    << "pos_x" << "pos_y" << "pos_z"
                    << recorder::data_recorder::end_name_row;
  
  recorder::ssv_recorder results_dt(results_filename + ".dt.ssv");
  
  results_dt << "time" << "pos_x" << "pos_y" << "pos_y"
                       << "q0" << "q1" << "q2" << "q3" << recorder::data_recorder::end_name_row;  

  sys_type::point_type x(13);
  x[0] = airship3D_frame->Position[0];
  x[1] = airship3D_frame->Position[1];
  x[2] = airship3D_frame->Position[2];
  x[3] = airship3D_frame->Quat[0];
  x[4] = airship3D_frame->Quat[1];
  x[5] = airship3D_frame->Quat[2];
  x[6] = airship3D_frame->Quat[3];
  x[7] = airship3D_frame->Velocity[0];
  x[8] = airship3D_frame->Velocity[1];
  x[9] = airship3D_frame->Velocity[2];
  x[10] = airship3D_frame->AngVelocity[0];
  x[11] = airship3D_frame->AngVelocity[1];
  x[12] = airship3D_frame->AngVelocity[2];
  
  ctrl::airship3D_inv_dt_system::point_type x_dt = x;
  
  std::cout << "Starting simulation..." << std::endl;
  try {
  for(double t = 0.0; t < end_time; t += time_step) {
    
    sys_type::input_type u = sys_type::input_type(6);
    u[0] = var_rnd() * sqrt(Qu(0,0));
    u[1] = var_rnd() * sqrt(Qu(1,1));
    u[2] = var_rnd() * sqrt(Qu(2,2));
    u[3] = var_rnd() * sqrt(Qu(3,3));
    u[4] = var_rnd() * sqrt(Qu(4,4));
    u[5] = var_rnd() * sqrt(Qu(5,5));
    
    x = airship3D_dt_sys.get_next_state(mdl_state_space,x,u,t);
    sys_type::output_type y = airship3D_dt_sys.get_output(mdl_state_space,x,u,t);
    
    x_dt = mdl_inv_dt.get_next_state(mdl_state_space,x_dt,u,t);
    ctrl::airship3D_inv_dt_system::output_type y_dt = mdl_inv_dt.get_output(mdl_state_space,x_dt,u,t);
    //std::cout << y_dt << std::endl;
    
    std::cout << "\r" << std::setw(20) << t; std::cout.flush();
    
    results << t 
            << (y[3] + var_rnd() * sqrt(R(3,3)))
	    << (y[4] + var_rnd() * sqrt(R(4,4)))
	    << (y[5] + var_rnd() * sqrt(R(5,5)))
	    << (y[6] + var_rnd() * sqrt(R(6,6))) 
	    << (y[0] + var_rnd() * sqrt(R(0,0))) 
	    << (y[1] + var_rnd() * sqrt(R(1,1)))
	    << (y[2] + var_rnd() * sqrt(R(2,2)))
	    << y[3] << y[4] << y[5] << y[6] 
            << y[0] << y[1] << y[2] 
            << recorder::data_recorder::end_value_row;
    
    results_dt << t << y_dt[0] << y_dt[1] << y_dt[2]
               << y_dt[3] << y_dt[4] << y_dt[5] << y_dt[6] << recorder::data_recorder::end_value_row;
  };
  
  } catch(impossible_integration& e) {
    std::cout << "Integration was deemed impossible, with message: '" << e.what() << "'" << std::endl;
    return 3;
  } catch(untolerable_integration& e) {
    std::cout << "Integration was deemed untolerable with message: '" << e.what() << "'" << std::endl;
    return 4;
  };
  
  results << recorder::data_recorder::flush;
  results_dt << recorder::data_recorder::flush;
  
  return 0;
};









