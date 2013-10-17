
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

#include "serialization/archiver_factory.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "airship3D_lin_model.hpp"



#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace po = boost::program_options;
namespace fs = boost::filesystem;


#include <iomanip>
#include <fstream>
#include <ctime>



int main(int argc, char** argv) {
  using namespace ReaK;
  
  
  po::options_description generic_options("Generic options");
  generic_options.add_options()
    ("help,h", "produce this help message.")
  ;
  
  po::options_description io_options("I/O options");
  io_options.add_options()
    ("model,m", po::value< std::string >()->default_value("models/airship3D.rkx"), "specify the model filename (default is 'models/airship3D.xml')")
    ("init,i", po::value< std::string >()->default_value("models/airship3D_init.rkx"), "specify the filename for the airship's initial conditions (default is 'models/airship3D_init.xml')")
    ("inertia,I", po::value< std::string >()->default_value("models/airship3D_inertia.rkx"), "specify the filename for the airship's inertial data (default is 'models/airship3D_inertia.xml')")
    ("Q-matrix,Q", po::value< std::string >()->default_value("models/airship3D_Q.rkx"), "specify the filename for the airship's input disturbance covariance matrix (default is 'models/airship3D_Q.xml')")
    ("R-matrix,R", po::value< std::string >()->default_value("models/airship3D_R.rkx"), "specify the filename for the airship's measurement noise covariance matrix (default is 'models/airship3D_R.xml')")
    ("output,o", po::value< std::string >()->default_value("results_airship3D/output_record"), "specify the path and filename (without extension) for the output of the results (default is 'results_airship3D/output_record')")
  ;
  
  po::options_description sim_options("Simulation options");
  sim_options.add_options()
    ("start-time,s", po::value< double >()->default_value(0.0), "start time of the simulation (default is 0.0)")
    ("end-time,e", po::value< double >()->default_value(1.0), "end time of the simulation (default is 1.0)")
    ("time-step,t", po::value< double >()->default_value(0.01), "time-step used for the simulations (default is 0.01)")
    ("sample-hops,r", po::value< std::size_t >()->default_value(1), "time-steps between each output point, i.e., sample-hops x time-step = sample-period (default is 1)")
  ;
  
  po::options_description output_options("Output options (at least one must be set)");
  output_options.add_options()
    ("xml,x", "if set, output resulting trajectory object in a XML file (rkx)")
    ("protobuf,p", "if set, output resulting trajectory object in a protobuf file (pbuf)")
    ("binary,b", "if set, output resulting trajectory object in a binary file (rkb)")
    ("ssv", "if set, output resulting trajectory time-series in a space-separated values file (ssv)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(sim_options).add(output_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help") || (vm.count("xml") + vm.count("protobuf") + vm.count("binary") + vm.count("ssv") < 1)) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  std::string output_path_name = vm["output-path"].as<std::string>();
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  fs::create_directory(output_path_name.c_str());
  
  
  std::string model_filename = vm["model"].as<std::string>();
  if( ! fs::exists( fs::path( model_filename ) ) ) {
    std::cout << "Model file does not exist!" << std::endl;
    return 2;
  };
  
  std::string init_filename = vm["init"].as<std::string>();
  if( ! fs::exists( fs::path( init_filename ) ) ) {
    std::cout << "Initial-conditions file does not exist!" << std::endl;
    return 3;
  };
  
  std::string inertia_filename = vm["inertia"].as<std::string>();
  if( ! fs::exists( fs::path( inertia_filename ) ) ) {
    std::cout << "Inertial-information file does not exist!" << std::endl;
    return 4;
  };
  
  std::string Qu_filename = vm["Q-matrix"].as<std::string>();
  if( ! fs::exists( fs::path( Qu_filename ) ) ) {
    std::cout << "Input disturbance covariance matrix file does not exist!" << std::endl;
    return 5;
  };
  
  std::string R_filename  = vm["R-matrix"].as<std::string>();
  if( ! fs::exists( fs::path( R_filename ) ) ) {
    std::cout << "Measurement noise covariance matrix file does not exist!" << std::endl;
    return 6;
  };
  
  
  
  std::string results_filename(argv[4]);
  
  
  
  double start_time = vm["start-time"].as<double>();
  double end_time   = vm["end-time"].as<double>();
  double time_step  = vm["time-step"].as<double>();
//   std::size_t hops  = vm["sample-hops"].as<std::size_t>();
  
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  
  shared_ptr< frame_3D<double> > airship3D_frame;
  shared_ptr< kte::position_measure_3D > airship3D_position;
  shared_ptr< kte::rotation_measure_3D > airship3D_rotation;
  
  shared_ptr< kte::driving_actuator_3D > airship3D_actuator;
  shared_ptr< kte::inertia_3D > airship3D_inertia;
  shared_ptr< kte::kte_map_chain > airship3D_model;
  shared_ptr< kte::mass_matrix_calc > airship3D_mass_calc;
  
  
  try {
    *(serialization::open_iarchive(model_filename))
      >> airship3D_frame
      >> airship3D_position
      >> airship3D_rotation
      >> airship3D_actuator
      >> airship3D_inertia
      >> airship3D_model
      >> airship3D_mass_calc;
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the airship's model!");
    return 10;
  };
  
  /* initial states */
  try {
    *(serialization::open_iarchive(init_filename))
      & RK_SERIAL_LOAD_WITH_ALIAS("initial_motion",(*airship3D_frame));
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 11;
  };
  
  /* inertial data */
  try {
    double tmp_mass;
    ReaK::mat<double,ReaK::mat_structure::symmetric> tmp_tensor;
    *(serialization::open_iarchive(inertia_filename))
      & RK_SERIAL_LOAD_WITH_ALIAS("mass",tmp_mass)
      & RK_SERIAL_LOAD_WITH_ALIAS("inertia_tensor",tmp_tensor);
    airship3D_inertia->setMass(tmp_mass);
    airship3D_inertia->setInertiaTensor(tmp_tensor);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 12;
  };
  
  /* input disturbance */
  mat<double,mat_structure::diagonal> Qu;
  try {
    *(serialization::open_iarchive(Qu_filename))
      & RK_SERIAL_LOAD_WITH_ALIAS("input_disturbance",Qu);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return 13;
  };
  
  /* measurement noise */
  mat<double,mat_structure::diagonal> R;
  try {
    *(serialization::open_iarchive(R_filename))
      & RK_SERIAL_LOAD_WITH_ALIAS("measurement_noise",R);
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return 14;
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
                      dormand_prince45_integrator<double>(
                        "ode45_integrator",
                        ReaK::vect_n<double>(),
                        start_time, time_step * 1e-2,
                        weak_ptr< state_rate_function<double> >(),
                        time_step, time_step * 0.000001, 1e-4),
                      time_step,
                      "airship3D_dt_sys");
  
  ctrl::airship3D_inv_dt_system mdl_inv_dt(
    "airship3D_invariant_discrete",
    airship3D_inertia->Mass(),
    airship3D_inertia->InertiaTensor(),
    time_step);
  
  pp::vector_topology< vect_n<double> > mdl_state_space;
  
  shared_ptr< std::ostream > ssv_ptr, ssv_dt_ptr;
  if( vm.count("ssv") ) {
    ssv_ptr = shared_ptr< std::ostream >(new std::ofstream((results_filename + ".ssv").c_str()));
    (*ssv_ptr) << "% time meas_q0 meas_q1 meas_q2 meas_q3 meas_x meas_y meas_z q0 q1 q2 q3 pos_x pos_y pos_z" << std::endl;
    ssv_dt_ptr = shared_ptr< std::ostream >(new std::ofstream((results_filename + "_dt.ssv").c_str()));
    (*ssv_dt_ptr) << "% time pos_x pos_y pos_z q0 q1 q2 q3" << std::endl;
  };
  
  
  
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
    
    for(double t = start_time; t < end_time; t += time_step) {
      
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
      
      std::cout << "\r" << std::setw(20) << t; std::cout.flush();
      
      if( vm.count("ssv") ) {
        (*ssv_ptr)  << t << " "
                    << (y[3] + var_rnd() * sqrt(R(3,3))) << " "
                    << (y[4] + var_rnd() * sqrt(R(4,4))) << " "
                    << (y[5] + var_rnd() * sqrt(R(5,5))) << " "
                    << (y[6] + var_rnd() * sqrt(R(6,6))) << " "
                    << (y[0] + var_rnd() * sqrt(R(0,0))) << " "
                    << (y[1] + var_rnd() * sqrt(R(1,1))) << " "
                    << (y[2] + var_rnd() * sqrt(R(2,2))) << " "
                    << y[3] << " " << y[4] << " " << y[5] << " " << y[6] << " "
                    << y[0] << " " << y[1] << " " << y[2] << std::endl;
        
        (*ssv_dt_ptr) << t << " " << y_dt[0] << " " << y_dt[1] << " " << y_dt[2] << " "
                      << y_dt[3] << " " << y_dt[4] << " " << y_dt[5] << " " << y_dt[6] << std::endl;
      };
    };
  
  } catch(impossible_integration& e) {
    std::cout << "Integration was deemed impossible, with message: '" << e.what() << "'" << std::endl;
    return 3;
  } catch(untolerable_integration& e) {
    std::cout << "Integration was deemed untolerable with message: '" << e.what() << "'" << std::endl;
    return 4;
  };
  
  
  return 0;
};









