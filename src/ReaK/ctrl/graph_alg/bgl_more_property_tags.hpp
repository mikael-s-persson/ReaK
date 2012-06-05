/**
 * \file bgl_more_property_tags.hpp
 *
 * This library adds a number of additional property-tags that can be used in a path-planning code.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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


#ifndef REAK_BGL_MORE_PROPERTY_TAGS_HPP
#define REAK_BGL_MORE_PROPERTY_TAGS_HPP

#include <boost/graph/properties.hpp>

namespace boost {

  enum vertex_heuristic_t { vertex_heuristic };
  enum vertex_rhs_t { vertex_rhs };
  enum vertex_key_t { vertex_key };
  enum vertex_position_t { vertex_position };
  
  BOOST_INSTALL_PROPERTY(vertex, heuristic);
  BOOST_INSTALL_PROPERTY(vertex, rhs);
  BOOST_INSTALL_PROPERTY(vertex, key);
  BOOST_INSTALL_PROPERTY(vertex, position);

};


#endif










