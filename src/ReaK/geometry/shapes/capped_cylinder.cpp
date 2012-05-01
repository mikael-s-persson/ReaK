
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

#include "capped_cylinder.hpp"

namespace ReaK {

namespace geom {



capped_cylinder::capped_cylinder(const std::string& aName,
                                 const shared_ptr< pose_3D<double> >& aAnchor,
		                 const pose_3D<double>& aPose,
		                 double aLength,
		                 double aRadius) :
		                 shape_3D(aName,aAnchor,aPose),
		                 mLength(aLength),
		                 mRadius(aRadius) { };
    
void RK_CALL capped_cylinder::save(ReaK::serialization::oarchive& A, unsigned int) const {
  shape_3D::save(A,shape_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mLength)
    & RK_SERIAL_SAVE_WITH_NAME(mRadius);
};

void RK_CALL capped_cylinder::load(ReaK::serialization::iarchive& A, unsigned int) {
  shape_3D::load(A,shape_3D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mLength)
    & RK_SERIAL_LOAD_WITH_NAME(mRadius);
};



};


};





