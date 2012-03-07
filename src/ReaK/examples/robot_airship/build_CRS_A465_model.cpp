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


#include "CRS_A465_models.hpp"


int main() {
  
  ReaK::robot_airship::CRS_A465_model_builder builder;
  
  builder.create_from_preset();
  
  builder.save_to_file("models/CRS_A465_raw_components.xml");
  
  return 0;
};







