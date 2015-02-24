
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

#include <ReaK/geometry/proximity/prox_box_box.hpp>

#include <ReaK/geometry/proximity/prox_fundamentals_3D.hpp>


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D compute_proximity(const box& aBox1, const shape_3D_precompute_pack& aPack1,
                                      const box& aBox2, const shape_3D_precompute_pack& aPack2) {
  return findProximityByGJKEPA(
    box_support_func(aBox1, aPack1.global_pose), 
    box_support_func(aBox2, aPack2.global_pose));
};


proximity_record_3D prox_box_box::computeProximity(const shape_3D_precompute_pack& aPack1, 
                                                   const shape_3D_precompute_pack& aPack2) {
  if((!mBox1) || (!mBox2))
    return proximity_record_3D();
  
  if( aPack1.parent == mBox1 )
    return compute_proximity(*mBox1, aPack1, *mBox2, aPack2);
  else // note, respect the order of packs.
    return compute_proximity(*mBox2, aPack1, *mBox1, aPack2);
};


prox_box_box::prox_box_box(const box* aBox1, const box* aBox2) :
                           proximity_finder_3D(),
                           mBox1(aBox1),
                           mBox2(aBox2) { };


};

};

