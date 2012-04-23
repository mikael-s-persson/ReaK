
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

#include "composite_shape_2D.hpp"

namespace ReaK {

namespace geom {



void composite_shape_2D::render() const {
  for(std::vector< shared_ptr< shape_2D > >::const_iterator it = mShapes.begin(); it != mShapes.end(); ++it)
    (*it)->render();
};
    
composite_shape_2D::composite_shape_2D(const std::string& aName) : shape_2D(aName), mShapes() { };

void RK_CALL composite_shape_2D::save(ReaK::serialization::oarchive& A, unsigned int) const {
  shape_2D::save(A,shape_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mShapes);
};

void RK_CALL composite_shape_2D::load(ReaK::serialization::iarchive& A, unsigned int) {
  shape_2D::load(A,shape_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mShapes);
};



};


};





