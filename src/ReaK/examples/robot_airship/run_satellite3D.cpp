
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

#include "serialization/archiver_factory.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

#include "satellite_invar_models.hpp"

#include "topologies/temporal_space.hpp"
#include "interpolation/discrete_point_trajectory.hpp"


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
    ("init,i", po::value< std::string >()->default_value("models/satellite3D_init.rkx"), "specify the filename for the satellite's initial conditions (default is 'models/satellite3D_init.rkx')")
    ("inertia,I", po::value< std::string >()->default_value("models/satellite3D_inertia.rkx"), "specify the filename for the satellite's inertial data (default is 'models/satellite3D_inertia.rkx')")
    ("Q-matrix,Q", po::value< std::string >()->default_value("models/satellite3D_Q.rkx"), "specify the filename for the satellite's input disturbance covariance matrix (default is 'models/satellite3D_Q.rkx')")
    ("R-matrix,R", po::value< std::string >()->default_value("models/satellite3D_R.rkx"), "specify the filename for the satellite's measurement noise covariance matrix (default is 'models/satellite3D_R.rkx')")
    ("IMU-config", po::value< std::string >()->default_value("models/satellite3D_IMU_config.rkx"), "specify the filename for the satellite's IMU configuration data, specifying its placement on the satellite and the inertial / magnetic-field frame it is relative to (default is 'models/satellite3D_IMU_config.rkx')")
    ("output,o", po::value< std::string >()->default_value("sim_results/satellite3D/output_record"), "specify the filename stem (without extension) for the output of the results (default is 'sim_results/satellite3D/output_record')")
    ("generate-files,g", "if set the output will be the generation of all the modeling files (with default values)")
    ("system-output", po::value< std::string >()->default_value("models/satellite3D_system.rkx"), "specify the filename for the output of the satellite system, when 'generate-files' is set (default is 'models/satellite3D_system.rkx')")
  ;
  
  po::options_description sim_options("Simulation options");
  sim_options.add_options()
    ("start-time,s", po::value< double >()->default_value(0.0), "start time of the simulation (default is 0.0)")
    ("end-time,e", po::value< double >()->default_value(1.0), "end time of the simulation (default is 1.0)")
    ("time-step,t", po::value< double >()->default_value(0.01), "time-step used for the simulations (default is 0.01)")
  ;
  
  po::options_description model_options("Modeling options");
  sim_options.add_options()
    ("gyro", "if set, a set of gyros is added to the model (angular velocity measurements). This requires the 'R-matrix' file to contain a 9x9 matrix.")
    ("IMU", "if set, a set of gyros is added to the model (angular velocity, magnetic field, and accelerometer measurements).\
 This requires the 'R-matrix' file to contain a 15x15 matrix. This option also automatically implies the 'midpoint' option.\
 This option will trigger the use of the 'IMU-config' file to obtain the information necessary about the IMU and the Earth's inertial frame.")
    ("midpoint", "if set, a midpoint rule is used to form the system's linearizations.")
  ;
  
  po::options_description output_options("Output options (at least one must be set)");
  output_options.add_options()
    ("xml,x", "if set, output resulting trajectory object in a XML file (rkx)")
    ("protobuf,p", "if set, output resulting trajectory object in a protobuf file (pbuf)")
    ("binary,b", "if set, output resulting trajectory object in a binary file (rkb)")
    ("ssv", "if set, output resulting trajectory time-series in a space-separated values file (ssv)")
  ;
  
  po::options_description cmdline_options;
  cmdline_options.add(generic_options).add(io_options).add(sim_options).add(model_options).add(output_options);
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, cmdline_options), vm);
  po::notify(vm);
  
  if(vm.count("help") || (vm.count("xml") + vm.count("protobuf") + vm.count("binary") + vm.count("ssv") + vm.count("generate-files") < 1)) {
    std::cout << cmdline_options << std::endl;
    return 1;
  };
  
  
  std::string output_path_name = vm["output"].as<std::string>();
  std::string output_stem_name = output_path_name;
  if(output_stem_name[output_stem_name.size()-1] == '/')
    output_stem_name += "output_record";
  else {
    std::size_t p = output_path_name.find_last_of('/');
    if(p == std::string::npos)
      output_path_name = "";
    else
      output_path_name.erase(p);
  };
  while(output_path_name[output_path_name.length()-1] == '/') 
    output_path_name.erase(output_path_name.length()-1, 1);
  
  if(!output_path_name.empty())
    fs::create_directory(output_path_name.c_str());
  
  
  
  double start_time = vm["start-time"].as<double>();
  double end_time   = vm["end-time"].as<double>();
  double time_step  = vm["time-step"].as<double>();
  
  
  boost::variate_generator< boost::minstd_rand, boost::normal_distribution<double> > var_rnd(boost::minstd_rand(static_cast<unsigned int>(time(NULL))), boost::normal_distribution<double>());
  
  
  /* initial states */
  frame_3D<double> initial_motion;
  try {
    
    std::string init_filename = vm["init"].as<std::string>();
    
    if( vm.count("generate-files") ) {
      *(serialization::open_oarchive(init_filename))
        & RK_SERIAL_SAVE_WITH_NAME(initial_motion);
    } else {
      
      if( ! fs::exists( fs::path( init_filename ) ) ) {
        std::cout << "Initial-conditions file does not exist!" << std::endl;
        return 3;
      };
      
      *(serialization::open_iarchive(init_filename))
        & RK_SERIAL_LOAD_WITH_NAME(initial_motion);
    };
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 11;
  };
  
  /* inertial data */
  double mass = 1.0;
  ReaK::mat<double,ReaK::mat_structure::symmetric> inertia_tensor(1.0, 0.0, 0.0, 1.0, 0.0, 1.0);
  try {
    
    std::string inertia_filename = vm["inertia"].as<std::string>();
    
    if( vm.count("generate-files") ) {
      *(serialization::open_oarchive(inertia_filename))
        & RK_SERIAL_SAVE_WITH_NAME(mass)
        & RK_SERIAL_SAVE_WITH_NAME(inertia_tensor);
    } else {
      if( ! fs::exists( fs::path( inertia_filename ) ) ) {
        std::cout << "Inertial-information file does not exist!" << std::endl;
        return 4;
      };
      
      *(serialization::open_iarchive(inertia_filename))
        & RK_SERIAL_LOAD_WITH_NAME(mass)
        & RK_SERIAL_LOAD_WITH_NAME(inertia_tensor);
    };
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the initial conditions!");
    return 12;
  };
  
  /* input disturbance */
  mat<double,mat_structure::diagonal> input_disturbance(6,true);
  try {
    
    std::string Qu_filename = vm["Q-matrix"].as<std::string>();
    
    if( vm.count("generate-files") ) {
      *(serialization::open_oarchive(Qu_filename))
        & RK_SERIAL_SAVE_WITH_NAME(input_disturbance);
    } else {
      if( ! fs::exists( fs::path( Qu_filename ) ) ) {
        std::cout << "Input disturbance covariance matrix file does not exist!" << std::endl;
        return 5;
      };
      
      *(serialization::open_iarchive(Qu_filename))
        & RK_SERIAL_LOAD_WITH_NAME(input_disturbance);
    };
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the input disturbance covariance matrix!");
    return 13;
  };
  
  /* measurement noise */
  std::size_t m_noise_size = 6;
  if( vm.count("gyro") )
    m_noise_size = 9;
  if( vm.count("IMU") )
    m_noise_size = 15;
  mat<double,mat_structure::diagonal> measurement_noise(m_noise_size,true);
  try {
    
    std::string R_filename  = vm["R-matrix"].as<std::string>();
    
    if( vm.count("generate-files") ) {
      *(serialization::open_oarchive(R_filename))
        & RK_SERIAL_SAVE_WITH_NAME(measurement_noise);
    } else {
      if( ! fs::exists( fs::path( R_filename ) ) ) {
        std::cout << "Measurement noise covariance matrix file does not exist!" << std::endl;
        return 6;
      };
      
      *(serialization::open_iarchive(R_filename))
        & RK_SERIAL_LOAD_WITH_NAME(measurement_noise);
    };
  } catch(...) {
    RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
    return 14;
  };
  
  unit_quat<double> IMU_orientation;
  vect<double,3> IMU_location;
  unit_quat<double> earth_orientation;
  vect<double,3> mag_field_direction(1.0,0.0,0.0);
  if( vm.count("IMU") ) {
    try {
      
      std::string IMUconf_filename  = vm["IMU-config"].as<std::string>();
      
      if( vm.count("generate-files") ) {
        *(serialization::open_oarchive(IMUconf_filename))
          & RK_SERIAL_SAVE_WITH_NAME(IMU_orientation)
          & RK_SERIAL_SAVE_WITH_NAME(IMU_location)
          & RK_SERIAL_SAVE_WITH_NAME(earth_orientation)
          & RK_SERIAL_SAVE_WITH_NAME(mag_field_direction);
      } else {
        if( ! fs::exists( fs::path( IMUconf_filename ) ) ) {
          std::cout << "IMU configuration data file does not exist!" << std::endl;
          return 6;
        };
        
        *(serialization::open_iarchive(IMUconf_filename))
          & RK_SERIAL_LOAD_WITH_NAME(IMU_orientation)
          & RK_SERIAL_LOAD_WITH_NAME(IMU_location)
          & RK_SERIAL_LOAD_WITH_NAME(earth_orientation)
          & RK_SERIAL_LOAD_WITH_NAME(mag_field_direction);
      };
    } catch(...) {
      RK_ERROR("An exception occurred during the loading of the measurement noise covariance matrix!");
      return 14;
    };
  };
  
  shared_ptr< ctrl::satellite3D_inv_dt_system > satellite3D_system;
  
  if( vm.count("midpoint") ) {
    
    if( vm.count("IMU") ) {
      satellite3D_system = shared_ptr< ctrl::satellite3D_inv_dt_system >(new ctrl::satellite3D_IMU_imdt_sys(
        "satellite3D_imdt_with_IMU", mass, inertia_tensor, time_step,
        IMU_orientation, IMU_location, earth_orientation, mag_field_direction));
    } else if( vm.count("gyro") ) {
      satellite3D_system = shared_ptr< ctrl::satellite3D_inv_dt_system >(new ctrl::satellite3D_gyro_imdt_sys(
        "satellite3D_imdt_with_gyros", mass, inertia_tensor, time_step));
    } else {
      satellite3D_system = shared_ptr< ctrl::satellite3D_inv_dt_system >(new ctrl::satellite3D_imdt_sys(
        "satellite3D_imdt", mass, inertia_tensor, time_step));
    };
    
  } else {
    
    if( vm.count("gyro") ) {
      satellite3D_system = shared_ptr< ctrl::satellite3D_inv_dt_system >(new ctrl::satellite3D_gyro_inv_dt_system(
        "satellite3D_idt_with_gyros", mass, inertia_tensor, time_step));
    } else {
      satellite3D_system = shared_ptr< ctrl::satellite3D_inv_dt_system >(new ctrl::satellite3D_inv_dt_system(
        "satellite3D_idt", mass, inertia_tensor, time_step));
    };
    
  };
  
  if( vm.count("generate-files") ) {
    
    try {
      
      std::string sys_filename  = vm["system-output"].as<std::string>();
      
      *(serialization::open_oarchive(sys_filename))
        & RK_SERIAL_SAVE_WITH_NAME(satellite3D_system);
      
    } catch(...) {
      RK_ERROR("An exception occurred during the saving the satellite system file!");
      return 14;
    };
    
    return 0;
  };
  
  shared_ptr< std::ostream > ssv_ptr;
  if( vm.count("ssv") ) {
    ssv_ptr = shared_ptr< std::ostream >(new std::ofstream((output_stem_name + ".ssv").c_str()));
    if( vm.count("IMU") ) {
      (*ssv_ptr) << "% time meas_x meas_y meas_z meas_q0 meas_q1 meas_q2 meas_q3" 
                 << " meas_wx meas_wy meas_wz"
                 << " meas_ax meas_ay meas_az"
                 << " meas_mx meas_my meas_mz"
                 << " pos_x pos_y pos_z q0 q1 q2 q3"
                 << " vel_x vel_y vel_z ang_vel_x ang_vel_y ang_vel_z" << std::endl;
    } else if( vm.count("gyro") ) {
      (*ssv_ptr) << "% time meas_x meas_y meas_z meas_q0 meas_q1 meas_q2 meas_q3" 
                 << " meas_wx meas_wy meas_wz"
                 << " pos_x pos_y pos_z q0 q1 q2 q3"
                 << " vel_x vel_y vel_z ang_vel_x ang_vel_y ang_vel_z" << std::endl;
    } else {
      (*ssv_ptr) << "% time meas_x meas_y meas_z meas_q0 meas_q1 meas_q2 meas_q3" 
                 << " pos_x pos_y pos_z q0 q1 q2 q3"  
                 << " vel_x vel_y vel_z ang_vel_x ang_vel_y ang_vel_z" << std::endl;
    };
  };
  
  
  typedef ctrl::satellite3D_inv_dt_system sys_type;
  typedef sys_type::state_space_type sat_state_space_type;
  typedef pp::temporal_space<sat_state_space_type, pp::time_poisson_topology, pp::time_distance_only> sat_temp_space_type;
  typedef pp::discrete_point_trajectory< sat_temp_space_type > sat_traj_type;
  typedef pp::topology_traits< sat_temp_space_type >::point_type temp_point_type;
  
  sat_temp_space_type sat_space(
    "satellite3D_temporal_space",
    sat_state_space_type(),
    pp::time_poisson_topology("satellite3D_time_space", time_step, (end_time - start_time) * 0.5));
  
  shared_ptr< sat_traj_type > traj_ptr;
  if( vm.count("xml") + vm.count("protobuf") + vm.count("binary") > 0 )
    traj_ptr = shared_ptr< sat_traj_type >(new sat_traj_type( shared_ptr< sat_temp_space_type >(&sat_space, null_deleter()) ));
  
  sys_type::point_type x;
  set_frame_3D(x, initial_motion);
  
  const mat<double,mat_structure::diagonal>& R  = measurement_noise;
  const mat<double,mat_structure::diagonal>& Qu = input_disturbance;
  
  double Rq0 = (R(3,3) + R(4,4) + R(5,5)) / 12.0;  // 12 = 3 * 4,  3 for average, 4 for axis-angle to quaternion.
  
  std::cout << "Starting simulation..." << std::endl;
  
  for(double t = start_time; t < end_time; t += time_step) {
    
    sys_type::input_type u = sys_type::input_type(6);
    u[0] = var_rnd() * sqrt(Qu(0,0));
    u[1] = var_rnd() * sqrt(Qu(1,1));
    u[2] = var_rnd() * sqrt(Qu(2,2));
    u[3] = var_rnd() * sqrt(Qu(3,3));
    u[4] = var_rnd() * sqrt(Qu(4,4));
    u[5] = var_rnd() * sqrt(Qu(5,5));
    
    x = satellite3D_system->get_next_state(sat_space.get_space_topology(), x, u, t);
    sys_type::output_type y = satellite3D_system->get_output(sat_space.get_space_topology(), x, u, t);
    
    std::cout << "\r" << std::setw(20) << t; std::cout.flush();
    
    if( vm.count("ssv") ) {
      
      (*ssv_ptr)  << t << " "
                  << (y[0] + var_rnd() * sqrt(R(0,0))) << " "
                  << (y[1] + var_rnd() * sqrt(R(1,1))) << " "
                  << (y[2] + var_rnd() * sqrt(R(2,2))) << " "
                  << (y[3] + var_rnd() * sqrt(Rq0)) << " "
                  << (y[4] + var_rnd() * sqrt(0.25 * R(3,3))) << " "
                  << (y[5] + var_rnd() * sqrt(0.25 * R(4,4))) << " "
                  << (y[6] + var_rnd() * sqrt(0.25 * R(5,5)));
      
      if( vm.count("IMU") || vm.count("gyro") ) {
        (*ssv_ptr)  << " " 
                    << (y[7] + var_rnd() * sqrt(R(6,6))) << " "
                    << (y[8] + var_rnd() * sqrt(R(7,7))) << " "
                    << (y[9] + var_rnd() * sqrt(R(8,8)));
      };
      
      if( vm.count("IMU") ) {
        (*ssv_ptr)  << " " 
                    << (y[10] + var_rnd() * sqrt(R(9,9))) << " "
                    << (y[11] + var_rnd() * sqrt(R(10,10))) << " "
                    << (y[12] + var_rnd() * sqrt(R(11,11))) << " "
                    << (y[13] + var_rnd() * sqrt(R(12,12))) << " "
                    << (y[14] + var_rnd() * sqrt(R(13,13))) << " "
                    << (y[15] + var_rnd() * sqrt(R(14,14)));
      };
      
      
      frame_3D<double> current_motion = get_frame_3D(x);
      
      (*ssv_ptr)  << " " 
                  << current_motion.Position[0] << " " << current_motion.Position[1] << " " << current_motion.Position[2] << " " 
                  << current_motion.Quat[0] << " " << current_motion.Quat[1] << " " << current_motion.Quat[2] << " " << current_motion.Quat[3] << " "
                  << current_motion.Velocity[0] << " " << current_motion.Velocity[1] << " " << current_motion.Velocity[2] << " " 
                  << current_motion.AngVelocity[0] << " " << current_motion.AngVelocity[1] << " " << current_motion.AngVelocity[2]
                  << std::endl;
      
    };
    
    if( vm.count("xml") + vm.count("protobuf") + vm.count("binary") > 0 ) {
      traj_ptr->push_back(temp_point_type(t, x));
    };
    
  };
  
  if( vm.count("xml") ) {
    
    *(serialization::open_oarchive(output_stem_name + ".rkx"))
      & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    
  };
  
  if( vm.count("protobuf") ) {
    
    *(serialization::open_oarchive(output_stem_name + ".pbuf"))
      & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    
  };
  
  if( vm.count("binary") ) {
    
    *(serialization::open_oarchive(output_stem_name + ".rkb"))
      & RK_SERIAL_SAVE_WITH_ALIAS("se3_trajectory", traj_ptr);
    
  };
  
  
  
  return 0;
};









