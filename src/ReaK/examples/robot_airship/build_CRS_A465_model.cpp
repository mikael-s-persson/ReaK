/**
 * \file build_CRS_A465_model.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2010
 */


#include "CRS_A465_geom_model.hpp"
#include "CRS_A465_2D_analog.hpp"


#include "shapes/plane.hpp"
#include "shapes/box.hpp"
#include "shapes/capped_cylinder.hpp"

#include "proximity/proxy_query_model.hpp"

#include "serialization/xml_archiver.hpp"



int main(int argc, char ** argv) {
  
  ReaK::robot_airship::CRS_A465_geom_builder builder;
  builder.create_geom_from_preset();
  
  builder.save_kte_to_file("models/CRS_A465_raw_components.xml");
  builder.save_kte_and_geom("models/CRS_A465_with_geom.xml");
  builder.save_limits_to_file("models/CRS_A465_limits.xml");
  
  
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > robot_proxy = builder.get_proximity_model();
  ReaK::shared_ptr< ReaK::kte::kte_map_chain > kin_chain = builder.get_kinematics_kte_chain();
  
  ReaK::robot_airship::CRS_A465_2D_model_builder builder2D;
  
  builder2D.create_from_preset();
  
  builder2D.save_kte_to_file("models/CRS_A465_2D_raw_components.xml");
  builder2D.save_limits_to_file("models/CRS_A465_2D_limits.xml");
  
  
  
  ReaK::shared_ptr< ReaK::geom::colored_model_3D > MD148_basic_lab = ReaK::shared_ptr< ReaK::geom::colored_model_3D >(new ReaK::geom::colored_model_3D("MD148_basic_lab_render"));
  ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D > MD148_lab_proxy = ReaK::shared_ptr< ReaK::geom::proxy_query_model_3D >(new ReaK::geom::proxy_query_model_3D("MD148_basic_lab_proxy"));
  {
    using namespace ReaK;
    
    shared_ptr< geom::plane > lab_floor = shared_ptr< geom::plane >(new geom::plane("MD148_floor", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(-0.8,-1.0,0.0), quaternion<double>()),
      vect<double,2>(4.0,6.0)));
    
    shared_ptr< geom::plane > lab_n_wall = shared_ptr< geom::plane >(new geom::plane("MD148_north_wall", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(1.2,-1.0,1.5), axis_angle<double>(M_PI * 0.5, vect<double,3>(0.0,-1.0,0.0)).getQuaternion()),
      vect<double,2>(3.0,6.0)));
    
    shared_ptr< geom::plane > lab_w_wall = shared_ptr< geom::plane >(new geom::plane("MD148_west_wall", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(-0.8,2.0,1.5), axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0)).getQuaternion()),
      vect<double,2>(4.0,3.0)));
    
    shared_ptr< geom::box > lab_robot_track = shared_ptr< geom::box >(new geom::box("MD148_robot_track", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(0.0,-1.71,0.15), quaternion<double>()),
      vect<double,3>(0.4,3.42,0.3)));
    
    shared_ptr< geom::capped_cylinder > lab_robot_track_left = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("MD148_robot_track_left", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(0.1,-1.71,0.15), axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0)).getQuaternion()),
      3.42, 0.18));
    
    shared_ptr< geom::capped_cylinder > lab_robot_track_right = shared_ptr< geom::capped_cylinder >(new geom::capped_cylinder("MD148_robot_track_right", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(-0.1,-1.71,0.15), axis_angle<double>(M_PI * 0.5, vect<double,3>(1.0,0.0,0.0)).getQuaternion()),
      3.42, 0.18));
    
    (*MD148_basic_lab)
     .addElement(geom::color(0,0,0),ReaK::shared_ptr< ReaK::geom::coord_arrows_3D >(new ReaK::geom::coord_arrows_3D("global_frame_arrows",shared_ptr< pose_3D<double> >(),pose_3D<double>(),1.0)))
     .addElement(geom::color(0.3,0.3,0.3),lab_floor)
     .addElement(geom::color(0.8,0.8,0.8),lab_n_wall)
     .addElement(geom::color(0.8,0.8,0.8),lab_w_wall)
     .addElement(geom::color(0.35,0.35,0.35),lab_robot_track);
    
    (*MD148_lab_proxy)
     .addShape(lab_floor)
     .addShape(lab_n_wall)
     .addShape(lab_w_wall)
     .addShape(lab_robot_track_left)
     .addShape(lab_robot_track_right);
    
    serialization::xml_oarchive out("models/MD148_lab_model.xml");
    
    out << MD148_basic_lab << MD148_lab_proxy;
    
  };
  
  return 0;
};







