/**
 * \file build_X8_obstacle_courses.cpp
 *
 * This application constructs the KTE-based kinematics models for the CRS A465 manipulator.
 * The inertial information is phony and is only there for completeness of the model but do 
 * not reflect the actual inertial information of the CRS A465.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date February 2013
 */



#include "shapes/plane.hpp"
#include "shapes/box.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "shapes/capped_cylinder.hpp"
#include "shapes/colored_model.hpp"

#include "proximity/proxy_query_model.hpp"

#include "serialization/xml_archiver.hpp"



int main(int argc, char ** argv) {
  
  using namespace ReaK;
  using namespace geom;
  
  shared_ptr< colored_model_3D >     one_building = shared_ptr< colored_model_3D >(new colored_model_3D("one_building_render"));
  shared_ptr< proxy_query_model_3D > one_building_proxy = shared_ptr< proxy_query_model_3D >(new proxy_query_model_3D("one_building_proxy"));
  {
    
    shared_ptr< plane > floor = shared_ptr< plane >(new plane("floor", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(2.5,2.5,0.0), quaternion<double>()),
      vect<double,2>(5.0,5.0)));
    
    shared_ptr< box > building = shared_ptr< box >(new box("building", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(2.5,2.5,2.5), quaternion<double>()),
      vect<double,3>(1.0,1.0,5.0)));
    
    (*one_building)
     .addElement(color(0,0,0),shared_ptr< coord_arrows_3D >(new coord_arrows_3D("global_frame_arrows",shared_ptr< pose_3D<double> >(),pose_3D<double>(),0.3)))
     .addElement(color(0.3,0.3,0.3),floor)
     .addElement(color(0.6,0.6,0.0),building);
    
    (*one_building_proxy)
     .addShape(floor)
     .addShape(building);
    
    serialization::xml_oarchive out("models/one_building.xml");
    
    out << one_building << one_building_proxy;
    
  };
  
  
  shared_ptr< colored_model_3D >     window_crossing = shared_ptr< colored_model_3D >(new colored_model_3D("window_crossing_render"));
  shared_ptr< proxy_query_model_3D > window_crossing_proxy = shared_ptr< proxy_query_model_3D >(new proxy_query_model_3D("window_crossing_proxy"));
  {
    
    shared_ptr< plane > floor = shared_ptr< plane >(new plane("floor", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(2.5,2.5,0.0), quaternion<double>()),
      vect<double,2>(5.0,5.0)));
    
    shared_ptr< box > wall1 = shared_ptr< box >(new box("wall1", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(1.5,0.5,2.5), quaternion<double>()),
      vect<double,3>(0.2,1.0,5.0)));
    
    shared_ptr< box > wall2 = shared_ptr< box >(new box("wall2", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(1.5,1.75,1.25), quaternion<double>()),
      vect<double,3>(0.2,1.5,2.5)));
    
    shared_ptr< box > wall3 = shared_ptr< box >(new box("wall3", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(1.5,1.75,4.5), quaternion<double>()),
      vect<double,3>(0.2,1.5,1.0)));
    
    shared_ptr< box > wall4 = shared_ptr< box >(new box("wall4", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(1.5,3.75,2.5), quaternion<double>()),
      vect<double,3>(0.2,2.5,5.0)));
    
    shared_ptr< box > wall5 = shared_ptr< box >(new box("wall5", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(3.5,1.25,2.5), quaternion<double>()),
      vect<double,3>(0.2,2.5,5.0)));
    
    shared_ptr< box > wall6 = shared_ptr< box >(new box("wall6", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(3.5,3.25,0.5), quaternion<double>()),
      vect<double,3>(0.2,1.5,1.0)));
    
    shared_ptr< box > wall7 = shared_ptr< box >(new box("wall7", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(3.5,3.25,3.75), quaternion<double>()),
      vect<double,3>(0.2,1.5,2.5)));
    
    shared_ptr< box > wall8 = shared_ptr< box >(new box("wall8", 
      shared_ptr< pose_3D<double> >(), 
      pose_3D<double>(weak_ptr< pose_3D<double> >(), vect<double,3>(3.5,4.5,2.5), quaternion<double>()),
      vect<double,3>(0.2,1.0,5.0)));
    
    (*window_crossing)
     .addElement(color(0,0,0),shared_ptr< coord_arrows_3D >(new coord_arrows_3D("global_frame_arrows",shared_ptr< pose_3D<double> >(),pose_3D<double>(),0.3)))
     .addElement(color(0.3,0.3,0.3),floor)
     .addElement(color(0.6,0.6,0.0),wall1)
     .addElement(color(0.6,0.6,0.0),wall2)
     .addElement(color(0.6,0.6,0.0),wall3)
     .addElement(color(0.6,0.6,0.0),wall4)
     .addElement(color(0.6,0.6,0.0),wall5)
     .addElement(color(0.6,0.6,0.0),wall6)
     .addElement(color(0.6,0.6,0.0),wall7)
     .addElement(color(0.6,0.6,0.0),wall8);
    
    (*window_crossing_proxy)
     .addShape(floor)
     .addShape(wall1)
     .addShape(wall2)
     .addShape(wall3)
     .addShape(wall4)
     .addShape(wall5)
     .addShape(wall6)
     .addShape(wall7)
     .addShape(wall8);
    
    serialization::xml_oarchive out("models/window_crossing.xml");
    
    out << window_crossing << window_crossing_proxy;
    
  };
  
  return 0;
};







