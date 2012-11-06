

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

#include "frame_tracer_coin3d_impl.hpp"


#include <Inventor/nodes/SoSeparator.h>
#include <Inventor/nodes/SoBaseColor.h>
#include <Inventor/nodes/SoRotation.h>
#include <Inventor/nodes/SoTranslation.h>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/nodes/SoCoordinate3.h>
#include <Inventor/nodes/SoLineSet.h>


namespace ReaK {

namespace geom {

  

struct tracer_coin3d_impl_pimpl {
  SoSeparator* path_sep;
  SoBaseColor* color;
  SoCoordinate3* coords;
  SoLineSet* ln_set;
  std::vector< vect<double,3> > current_points;
}; 

SoSeparator* tracer_coin3d_impl::get_separator() const {
  return path_impl->path_sep;
};

tracer_coin3d_impl::tracer_coin3d_impl(bool aIsSolution) {
  path_impl = new tracer_coin3d_impl_pimpl;
  path_impl->path_sep = new SoSeparator;
  path_impl->path_sep->ref();
  
  path_impl->color = new SoBaseColor;
  if(aIsSolution)
    path_impl->color->rgb = SbColor(1.0, 0.0, 0.0);
  else
    path_impl->color->rgb = SbColor(1.0, 0.6, 0.0);
  path_impl->path_sep->addChild(path_impl->color);
  
  path_impl->coords = new SoCoordinate3;
  path_impl->path_sep->addChild(path_impl->coords);
    
  path_impl->ln_set = new SoLineSet;
  path_impl->path_sep->addChild(path_impl->ln_set);
  
};

tracer_coin3d_impl::tracer_coin3d_impl(const tracer_coin3d_impl& rhs) {
  path_impl = new tracer_coin3d_impl_pimpl;
  path_impl->path_sep = new SoSeparator;
  path_impl->path_sep->ref();
  
  path_impl->color = new SoBaseColor;
  path_impl->color->rgb = rhs.path_impl->color->rgb;
  path_impl->path_sep->addChild(path_impl->color);
  
  path_impl->coords = new SoCoordinate3;
  path_impl->path_sep->addChild(path_impl->coords);
    
  path_impl->ln_set = new SoLineSet;
  path_impl->path_sep->addChild(path_impl->ln_set);
};

tracer_coin3d_impl& tracer_coin3d_impl::operator=(const tracer_coin3d_impl& rhs) {
  path_impl->path_sep->unref();
  
  path_impl->path_sep = new SoSeparator;
  path_impl->path_sep->ref();
  
  path_impl->color = new SoBaseColor;
  path_impl->color->rgb = rhs.path_impl->color->rgb;
  path_impl->path_sep->addChild(path_impl->color);
  
  path_impl->coords = new SoCoordinate3;
  path_impl->path_sep->addChild(path_impl->coords);
    
  path_impl->ln_set = new SoLineSet;
  path_impl->path_sep->addChild(path_impl->ln_set);
  
  path_impl->current_points.clear();
  
  return *this;
};


tracer_coin3d_impl::~tracer_coin3d_impl() {
  path_impl->path_sep->unref();
  delete path_impl;
};

void tracer_coin3d_impl::begin_edge(const vect<double,3>& start_point) const {
  if(!(path_impl->current_points.empty()))
    end_edge();
  if(std::isnan(start_point[0]) || std::isnan(start_point[1]) || std::isnan(start_point[2]))
    return;
  path_impl->current_points.push_back(start_point);
};

void tracer_coin3d_impl::add_point(const vect<double,3>& new_point) const {
  if(std::isnan(new_point[0]) || std::isnan(new_point[1]) || std::isnan(new_point[2]))
    return;
  path_impl->current_points.push_back(new_point);
};

void tracer_coin3d_impl::end_edge() const {
  std::size_t offset = path_impl->coords->point.getNum();
  
  for(std::size_t i = 0; i < path_impl->current_points.size(); ++i) {
    const vect<double,3>& v = path_impl->current_points[i];
    path_impl->coords->point.set1Value(offset + i, v[0], v[1], v[2]);
  };
  
  std::size_t ln_offset = path_impl->ln_set->numVertices.getNum();
  path_impl->ln_set->numVertices.set1Value(ln_offset, path_impl->current_points.size());
  
  path_impl->current_points.clear();
};




};


};





