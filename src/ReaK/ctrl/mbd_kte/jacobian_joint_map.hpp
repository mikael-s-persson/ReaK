/**
 * \file jacobian_joint_map.hpp
 * 
 * This library declares types for jacobian mappings of generalized coordinates. These jacobian frames are
 * associated to the motion in a generalized coordinate by the joints that have these generalized coordinates
 * as input.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef JACOBIAN_JOINT_MAP_HPP
#define JACOBIAN_JOINT_MAP_HPP

#include "math/kinetostatics.hpp"
#include "math/motion_jacobians.hpp"
#include <boost/shared_ptr.hpp>
#include <map>

namespace ReaK {

namespace kte {


/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< boost::shared_ptr< gen_coord<double> >, 
                  boost::shared_ptr< jacobian_gen_gen<double> > > jacobian_joint_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< boost::shared_ptr< gen_coord<double> >, 
                  boost::shared_ptr< jacobian_gen_2D<double> > > jacobian_joint_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< boost::shared_ptr< gen_coord<double> >, 
                  boost::shared_ptr< jacobian_gen_3D<double> > > jacobian_joint_map_3D;


/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< boost::shared_ptr< frame_2D<double> >, 
                  boost::shared_ptr< jacobian_2D_gen<double> > > jacobian_joint2D_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< boost::shared_ptr< frame_2D<double> >, 
                  boost::shared_ptr< jacobian_2D_2D<double> > > jacobian_joint2D_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< boost::shared_ptr< frame_2D<double> >, 
                  boost::shared_ptr< jacobian_2D_3D<double> > > jacobian_joint2D_map_3D;

		  
/** This typedef declares a mapping to associate generalized coordinates to their Jacobian generalized coordinate. */
typedef std::map< boost::shared_ptr< frame_3D<double> >, 
                  boost::shared_ptr< jacobian_3D_gen<double> > > jacobian_joint3D_map_gen;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 2D frame. */
typedef std::map< boost::shared_ptr< frame_3D<double> >, 
                  boost::shared_ptr< jacobian_3D_2D<double> > > jacobian_joint3D_map_2D;

/** This typedef declares a mapping to associate generalized coordinates to their Jacobian 3D frame. */
typedef std::map< boost::shared_ptr< frame_3D<double> >, 
                  boost::shared_ptr< jacobian_3D_3D<double> > > jacobian_joint3D_map_3D;



};

};

#endif //JACOBIAN_JOINT_MAP_HPP













