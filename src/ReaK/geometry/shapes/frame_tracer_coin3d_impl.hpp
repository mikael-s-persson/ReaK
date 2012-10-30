/**
 * \file frame_tracer_coin3d.hpp
 * 
 * This library defines a sampling-based motion/path planning reporter for tracing out the path of a frame 
 * linked to a DK kinematic model. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_FRAME_TRACER_COIN3D_HPP
#define REAK_FRAME_TRACER_COIN3D_HPP

#include "base/defs.hpp"

#include "lin_alg/vect_alg.hpp"


class SoSeparator;   // forward-declare


/** Main namespace for ReaK */
namespace ReaK {


/** Main namespace for ReaK.Path-Planning */
namespace geom {

struct tracer_coin3d_impl_pimpl; // forward-declare

/**
 * This class can be used
 */
struct tracer_coin3d_impl {
  
  tracer_coin3d_impl_pimpl* path_impl;
  
  SoSeparator* get_separator() const;
  
  explicit tracer_coin3d_impl(bool aIsSolution = false);
  ~tracer_coin3d_impl();
  
  void begin_edge(const vect<double,3>& start_point) const;
  
  void add_point(const vect<double,3>& new_point) const;
  
  void end_edge() const;
  
};


};

};


#endif

