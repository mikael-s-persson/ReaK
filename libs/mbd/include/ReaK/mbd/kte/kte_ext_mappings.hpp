/**
 * \file kte_ext_mappings.hpp
 *
 * This library declares extended KTE-related mappings. Most importantly, the maps for storing 
 * the kinetostatic frame states pertaining to end and intermediate frames in a chain of KTE models.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
 */

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

#ifndef REAK_KTE_EXT_MAPPINGS_HPP
#define REAK_KTE_EXT_MAPPINGS_HPP

#include <ReaK/math/lin_alg/vect_alg.hpp>
#include <ReaK/math/kinetostatics/kinetostatics.hpp>

#include <map>

namespace ReaK {

namespace kte {


struct velocity_coef_gen {
  double v;
  velocity_coef_gen() : v(0.0) { };
};

struct acceleration_coef_gen {
  double a;
  acceleration_coef_gen() : a(0.0) { };
};

struct force_coef_gen {
  double f;
  force_coef_gen() : f(0.0) { };
};


struct velocity_coef_2D {
  vect<double,2> v;
  double omega;
  velocity_coef_2D() : v(), omega(0.0) { };
};

struct acceleration_coef_2D {
  vect<double,2> a;
  double alpha;
  acceleration_coef_2D() : a(), alpha(0.0) { };
};

struct force_coef_2D {
  vect<double,2> f;
  double tau;
  force_coef_2D() : f(), tau(0.0) { };
};


struct velocity_coef_3D {
  vect<double,3> v;
  vect<double,3> omega;
  velocity_coef_3D() : v(), omega() { };
};

struct acceleration_coef_3D {
  vect<double,3> a;
  vect<double,3> alpha;
  acceleration_coef_3D() : a(), alpha() { };
};

struct force_coef_3D {
  vect<double,3> f;
  vect<double,3> tau;
  force_coef_3D() : f(), tau() { };
};


typedef std::map< double*, velocity_coef_gen >     velocity_coef_map_gen;
typedef std::map< double*, acceleration_coef_gen > acceleration_coef_map_gen;
typedef std::map< double*, force_coef_gen >        force_coef_map_gen;
typedef std::map< double*, velocity_coef_2D >      velocity_coef_map_2D;
typedef std::map< double*, acceleration_coef_2D >  acceleration_coef_map_2D;
typedef std::map< double*, force_coef_2D >         force_coef_map_2D;
typedef std::map< double*, velocity_coef_3D >      velocity_coef_map_3D;
typedef std::map< double*, acceleration_coef_3D >  acceleration_coef_map_3D;
typedef std::map< double*, force_coef_3D >         force_coef_map_3D;


/** This typedef defines a map of generalized coordinate state storage associated to a generalized coordinate pointer in the KTE chain. */ 
typedef std::map< shared_ptr< gen_coord<double> >, shared_ptr< gen_coord<double> > > gen_coord_map;
/** This typedef defines a map of frame 2D state storage associated to a frame 2D pointer in the KTE chain. */ 
typedef std::map< shared_ptr< frame_2D<double> >, shared_ptr< frame_2D<double> > > frame_2D_map;
/** This typedef defines a map of frame 3D state storage associated to a frame 3D pointer in the KTE chain. */ 
typedef std::map< shared_ptr< frame_3D<double> >, shared_ptr< frame_3D<double> > > frame_3D_map;

/**
 * This struct is used as a storage repository to take a "flash" of all kinematics and dynamics states 
 * of the kinetostatic frames (2D, 3D and generalized coordinates) at an instant.  
 */
struct frame_storage {
  gen_coord_map gen_coord_mapping; ///< Stores the generalized coordinate states.
  frame_2D_map frame_2D_mapping; ///< Stores the frame 2D states.
  frame_3D_map frame_3D_mapping; ///< Stores the frame 3D states.
  /** Default constructor, non-POD. */
  frame_storage() : gen_coord_mapping(), frame_2D_mapping(), frame_3D_mapping() { };
};

};

};

#endif





