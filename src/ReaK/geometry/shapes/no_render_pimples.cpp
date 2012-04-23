
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

#include "line_seg_2D.hpp"
#include "line_seg_3D.hpp"
#include "coord_arrows_2D.hpp"
#include "coord_arrows_3D.hpp"
#include "grid_2D.hpp"
#include "grid_3D.hpp"
#include "circle.hpp"
#include "ellipse.hpp"
#include "rectangle.hpp"

#include "plane.hpp"
#include "sphere.hpp"
#include "ellipsoid.hpp"
#include "cylinder.hpp"
#include "box.hpp"

namespace ReaK {

namespace geom {

void line_seg_2D::render() const { };

void line_seg_3D::render() const { };

void coord_arrows_2D::render() const { };

void coord_arrows_3D::render() const { };

void grid_2D::render() const { };

void grid_3D::render() const { };

void circle::render() const { };

void ellipse::render() const { };

void rectangle::render() const { };


void plane::render() const { };

void sphere::render() const { };

void ellipsoid::render() const { };

void cylinder::render() const { };

void box::render() const { };

};


};





