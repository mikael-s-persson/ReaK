
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

#include <ReaK/geometry/shapes/grid_3D.hpp>

namespace ReaK {

namespace geom {


grid_3D::grid_3D( const std::string& aName, const shared_ptr< pose_3D< double > >& aAnchor,
                  const pose_3D< double >& aPose, const vect< double, 3 >& aDimensions,
                  const vect< std::size_t, 3 >& aSquareCounts )
    : geometry_3D( aName, aAnchor, aPose ), mDimensions( aDimensions ), mSquareCounts( aSquareCounts ){};


void RK_CALL grid_3D::save( ReaK::serialization::oarchive& A, unsigned int ) const {
  geometry_3D::save( A, geometry_3D::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_SAVE_WITH_NAME( mDimensions ) & RK_SERIAL_SAVE_WITH_NAME( mSquareCounts );
};

void RK_CALL grid_3D::load( ReaK::serialization::iarchive& A, unsigned int ) {
  geometry_3D::load( A, geometry_3D::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_LOAD_WITH_NAME( mDimensions ) & RK_SERIAL_LOAD_WITH_NAME( mSquareCounts );
};
};
};
