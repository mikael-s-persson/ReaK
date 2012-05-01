
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

#include "capped_rectangle.hpp"

namespace ReaK {

namespace geom {


capped_rectangle::capped_rectangle(const std::string& aName,
		                   const shared_ptr< pose_2D<double> >& aAnchor,
		                   const pose_2D<double>& aPose,
		                   const vect<double,2>& aDimensions) :
		                   shape_2D(aName,aAnchor,aPose),
		                   mDimensions(aDimensions) { };
    
    
void RK_CALL capped_rectangle::save(ReaK::serialization::oarchive& A, unsigned int) const {
  shape_2D::save(A,shape_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mDimensions);
};

void RK_CALL capped_rectangle::load(ReaK::serialization::iarchive& A, unsigned int) {
  shape_2D::load(A,shape_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mDimensions);
};




};


};





