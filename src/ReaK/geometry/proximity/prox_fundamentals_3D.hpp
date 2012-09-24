/**
 * \file prox_fundamentals_3D.hpp
 *
 * This library declares fundamental proximity query methods for 3D.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2012
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

#ifndef REAK_PROX_FUNDAMENTALS_3D_HPP
#define REAK_PROX_FUNDAMENTALS_3D_HPP

#include "proximity_finder_3D.hpp"

#include "shapes/box.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D findProximityBoxToPoint(const shared_ptr< box >& aBox, const vect<double,3>& aPoint);


proximity_record_3D findProximityBoxToLine(const shared_ptr< box >& aBox, const vect<double,3>& aCenter, const vect<double,3>& aTangent, double aHalfLength);



};

};

#endif










