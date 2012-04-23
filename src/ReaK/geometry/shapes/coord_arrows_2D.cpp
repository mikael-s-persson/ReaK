
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

#include "coord_arrows_2D.hpp"

namespace ReaK {

namespace geom {



coord_arrows_2D::coord_arrows_2D(const std::string& aName,
				 const shared_ptr< pose_2D<double> >& aAnchor,
				 const pose_2D<double>& aPose,
				 double aArrowLength) :
				 geometry_2D(aName,aAnchor,aPose),
				 mArrowLength(aArrowLength) { };
    
void RK_CALL coord_arrows_2D::save(ReaK::serialization::oarchive& A, unsigned int) const {
  geometry_2D::save(A,geometry_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(mArrowLength);
};

void RK_CALL coord_arrows_2D::load(ReaK::serialization::iarchive& A, unsigned int) {
  geometry_2D::load(A,geometry_2D::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(mArrowLength);
};



};


};





