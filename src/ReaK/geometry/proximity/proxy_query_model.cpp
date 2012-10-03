
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

#include "proxy_query_model.hpp"

#include "shapes/line_seg_2D.hpp"
#include "shapes/grid_2D.hpp"
#include "shapes/coord_arrows_2D.hpp"
#include "shapes/circle.hpp"
#include "shapes/rectangle.hpp"
#include "shapes/capped_rectangle.hpp"
#include "shapes/composite_shape_2D.hpp"

#include "shapes/line_seg_3D.hpp"
#include "shapes/grid_3D.hpp"
#include "shapes/coord_arrows_3D.hpp"
#include "shapes/plane.hpp"
#include "shapes/sphere.hpp"
#include "shapes/box.hpp"
#include "shapes/cylinder.hpp"
#include "shapes/capped_cylinder.hpp"
#include "shapes/composite_shape_3D.hpp"


#include "prox_circle_circle.hpp"
#include "prox_circle_crect.hpp"
#include "prox_circle_rectangle.hpp"
#include "prox_crect_crect.hpp"
#include "prox_crect_rectangle.hpp"
#include "prox_rectangle_rectangle.hpp"

#include "prox_plane_plane.hpp"
#include "prox_plane_sphere.hpp"
#include "prox_plane_ccylinder.hpp"
#include "prox_plane_cylinder.hpp"
#include "prox_plane_box.hpp"
#include "prox_sphere_sphere.hpp"
#include "prox_sphere_ccylinder.hpp"
#include "prox_sphere_cylinder.hpp"
#include "prox_sphere_box.hpp"
#include "prox_ccylinder_ccylinder.hpp"
#include "prox_ccylinder_cylinder.hpp"     // NOTE: not working.
#include "prox_ccylinder_box.hpp"
#include "prox_cylinder_cylinder.hpp"      // NOTE: not working.
#include "prox_cylinder_box.hpp"           // NOTE: not working.
#include "prox_box_box.hpp"                // NOTE: not working.


namespace ReaK {

namespace geom {



void proxy_query_pair_2D::createProxFinderList() {
  mProxFinders.clear();
  if(!mModel1 || !mModel2)
    return;
  
  // for all shapes in mModel1
  for(std::size_t i = 0; i < mModel1->mShapeList.size(); ++i) {
    if(!mModel1->mShapeList[i])
      continue;
    // for all shapes in mModel2
    for(std::size_t j = 0; j < mModel2->mShapeList.size(); ++j) {
      if(!mModel2->mShapeList[j])
        continue;
      // if one of the model is a circle?
      if((mModel1->mShapeList[i]->getObjectType() == circle::getStaticObjectType()) ||
         (mModel2->mShapeList[j]->getObjectType() == circle::getStaticObjectType())) {
        shared_ptr<circle> ci_geom;
        shared_ptr<shape_2D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == circle::getStaticObjectType()) {
          ci_geom = rtti::rk_static_ptr_cast<circle>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          ci_geom = rtti::rk_static_ptr_cast<circle>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a circle..
        if(other_geom->getObjectType() == circle::getStaticObjectType()) {
          shared_ptr<circle> ci2_geom = rtti::rk_static_ptr_cast<circle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_circle_circle >(new prox_circle_circle(ci_geom, ci2_geom)));
        }
        // if the other is a capped-rectangle..
        else if(other_geom->getObjectType() == capped_rectangle::getStaticObjectType()) {
          shared_ptr<capped_rectangle> cr_geom = rtti::rk_static_ptr_cast<capped_rectangle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_circle_crect >(new prox_circle_crect(ci_geom, cr_geom)));
        }
        // if the other is a rectangle..
        else if(other_geom->getObjectType() == rectangle::getStaticObjectType()) {
          shared_ptr<rectangle> re_geom = rtti::rk_static_ptr_cast<rectangle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_circle_rectangle >(new prox_circle_rectangle(ci_geom, re_geom)));
        };
      }
      // if one of the model is a capped-rectangle?
      else if((mModel1->mShapeList[i]->getObjectType() == capped_rectangle::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == capped_rectangle::getStaticObjectType())) {
        shared_ptr<capped_rectangle> cr_geom;
        shared_ptr<shape_2D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == capped_rectangle::getStaticObjectType()) {
          cr_geom = rtti::rk_static_ptr_cast<capped_rectangle>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          cr_geom = rtti::rk_static_ptr_cast<capped_rectangle>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a capped-rectangle..
        if(other_geom->getObjectType() == capped_rectangle::getStaticObjectType()) {
          shared_ptr<capped_rectangle> cr2_geom = rtti::rk_static_ptr_cast<capped_rectangle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_crect_crect >(new prox_crect_crect(cr_geom, cr2_geom)));
        }
        // if the other is a rectangle..
        else if(other_geom->getObjectType() == rectangle::getStaticObjectType()) {
          shared_ptr<rectangle> re_geom = rtti::rk_static_ptr_cast<rectangle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_crect_rectangle >(new prox_crect_rectangle(cr_geom, re_geom)));
        };
      }
      // if one of the model is a rectangle?
      else if((mModel1->mShapeList[i]->getObjectType() == rectangle::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == rectangle::getStaticObjectType())) {
        shared_ptr<rectangle> re_geom;
        shared_ptr<shape_2D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == rectangle::getStaticObjectType()) {
          re_geom = rtti::rk_static_ptr_cast<rectangle>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          re_geom = rtti::rk_static_ptr_cast<rectangle>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a rectangle..
        if(other_geom->getObjectType() == rectangle::getStaticObjectType()) {
          shared_ptr<rectangle> re2_geom = rtti::rk_static_ptr_cast<rectangle>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_rectangle_rectangle >(new prox_rectangle_rectangle(re_geom, re2_geom)));
        };
      };
    };
  };
  
  return;
};

shared_ptr< proximity_finder_2D > proxy_query_pair_2D::findMinimumDistance() const {
  if(mProxFinders.empty())
    return shared_ptr< proximity_finder_2D >();
  
  std::size_t min_i = 0;
  mProxFinders[0]->computeProximity();
  double min_dist = mProxFinders[0]->getLastResult().mDistance;
  
  for(std::size_t i = 1; i < mProxFinders.size(); ++i) {
    mProxFinders[i]->computeProximity();
    if(min_dist > mProxFinders[i]->getLastResult().mDistance) {
      min_i = i;
      min_dist = mProxFinders[i]->getLastResult().mDistance;
    };
  };
  
  return mProxFinders[min_i];
};
    
bool proxy_query_pair_2D::gatherCollisionPoints(std::vector< proximity_record_2D >& aOutput) const {
  bool collision_found = false;
  
  // TODO: Add a circle-circle collision test on each shape-pair before calling "computeProximity".
  
  for(std::size_t i = 0; i < mProxFinders.size(); ++i) {
    mProxFinders[i]->computeProximity();
    if(mProxFinders[i]->getLastResult().mDistance < 0.0) {
      aOutput.push_back(mProxFinders[i]->getLastResult());
      collision_found = true;
    };
  };
  
  return collision_found;
};






void proxy_query_pair_3D::createProxFinderList() {
  mProxFinders.clear();
  if(!mModel1 || !mModel2)
    return;
  
  // for all shapes in mModel1
  for(std::size_t i = 0; i < mModel1->mShapeList.size(); ++i) {
    if(!mModel1->mShapeList[i])
      continue;
    // for all shapes in mModel2
    for(std::size_t j = 0; j < mModel2->mShapeList.size(); ++j) {
      if(!mModel2->mShapeList[j])
        continue;
      
      // if one of the model is a plane?
      if((mModel1->mShapeList[i]->getObjectType() == plane::getStaticObjectType()) ||
         (mModel2->mShapeList[j]->getObjectType() == plane::getStaticObjectType())) {
        shared_ptr<plane> pl_geom;
        shared_ptr<shape_3D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == plane::getStaticObjectType()) {
          pl_geom = rtti::rk_static_ptr_cast<plane>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          pl_geom = rtti::rk_static_ptr_cast<plane>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a plane..
        if(other_geom->getObjectType() == plane::getStaticObjectType()) {
          shared_ptr<plane> pl2_geom = rtti::rk_static_ptr_cast<plane>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_plane_plane >(new prox_plane_plane(pl_geom, pl2_geom)));
        }
        // if the other is a sphere..
        else if(other_geom->getObjectType() == sphere::getStaticObjectType()) {
          shared_ptr<sphere> sp_geom = rtti::rk_static_ptr_cast<sphere>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_plane_sphere >(new prox_plane_sphere(pl_geom, sp_geom)));
        }
        // if the other is a ccylinder..
        else if(other_geom->getObjectType() == capped_cylinder::getStaticObjectType()) {
          shared_ptr<capped_cylinder> cc_geom = rtti::rk_static_ptr_cast<capped_cylinder>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_plane_ccylinder >(new prox_plane_ccylinder(pl_geom, cc_geom)));
        }
        // if the other is a cylinder..
        else if(other_geom->getObjectType() == cylinder::getStaticObjectType()) {
          shared_ptr<cylinder> cy_geom = rtti::rk_static_ptr_cast<cylinder>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_plane_cylinder >(new prox_plane_cylinder(pl_geom, cy_geom)));
        }
        // if the other is a box..
        else if(other_geom->getObjectType() == box::getStaticObjectType()) {
          shared_ptr<box> bx_geom = rtti::rk_static_ptr_cast<box>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_plane_box >(new prox_plane_box(pl_geom, bx_geom)));
        };
      }
      // if one of the model is a sphere?
      else if((mModel1->mShapeList[i]->getObjectType() == sphere::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == sphere::getStaticObjectType())) {
        shared_ptr<sphere> sp_geom;
        shared_ptr<shape_3D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == sphere::getStaticObjectType()) {
          sp_geom = rtti::rk_static_ptr_cast<sphere>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          sp_geom = rtti::rk_static_ptr_cast<sphere>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a sphere..
        if(other_geom->getObjectType() == sphere::getStaticObjectType()) {
          shared_ptr<sphere> sp2_geom = rtti::rk_static_ptr_cast<sphere>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_sphere_sphere >(new prox_sphere_sphere(sp_geom, sp2_geom)));
        }
        // if the other is a ccylinder..
        else if(other_geom->getObjectType() == capped_cylinder::getStaticObjectType()) {
          shared_ptr<capped_cylinder> cc_geom = rtti::rk_static_ptr_cast<capped_cylinder>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_sphere_ccylinder >(new prox_sphere_ccylinder(sp_geom, cc_geom)));
        }
        // if the other is a cylinder..
        else if(other_geom->getObjectType() == cylinder::getStaticObjectType()) {
          shared_ptr<cylinder> cy_geom = rtti::rk_static_ptr_cast<cylinder>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_sphere_cylinder >(new prox_sphere_cylinder(sp_geom, cy_geom)));
        }
        // if the other is a box..
        else if(other_geom->getObjectType() == box::getStaticObjectType()) {
          shared_ptr<box> bx_geom = rtti::rk_static_ptr_cast<box>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_sphere_box >(new prox_sphere_box(sp_geom, bx_geom)));
        };
      }
      // if one of the model is a ccylinder?
      else if((mModel1->mShapeList[i]->getObjectType() == capped_cylinder::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == capped_cylinder::getStaticObjectType())) {
        shared_ptr<capped_cylinder> cc_geom;
        shared_ptr<shape_3D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == capped_cylinder::getStaticObjectType()) {
          cc_geom = rtti::rk_static_ptr_cast<capped_cylinder>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          cc_geom = rtti::rk_static_ptr_cast<capped_cylinder>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a ccylinder..
        if(other_geom->getObjectType() == capped_cylinder::getStaticObjectType()) {
          shared_ptr<capped_cylinder> cc2_geom = rtti::rk_static_ptr_cast<capped_cylinder>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_ccylinder_ccylinder >(new prox_ccylinder_ccylinder(cc_geom, cc2_geom)));
        }
        // if the other is a cylinder..
        else if(other_geom->getObjectType() == cylinder::getStaticObjectType()) {
          shared_ptr<cylinder> cy_geom = rtti::rk_static_ptr_cast<cylinder>(other_geom);
          //mProxFinders.push_back(shared_ptr< prox_ccylinder_cylinder >(new prox_ccylinder_cylinder(cc_geom, cy_geom)));
        }
        // if the other is a box..
        else if(other_geom->getObjectType() == box::getStaticObjectType()) {
          shared_ptr<box> bx_geom = rtti::rk_static_ptr_cast<box>(other_geom);
          mProxFinders.push_back(shared_ptr< prox_ccylinder_box >(new prox_ccylinder_box(cc_geom, bx_geom)));
        };
      }
      // if one of the model is a cylinder?
      else if((mModel1->mShapeList[i]->getObjectType() == cylinder::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == cylinder::getStaticObjectType())) {
        shared_ptr<cylinder> cy_geom;
        shared_ptr<shape_3D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == cylinder::getStaticObjectType()) {
          cy_geom = rtti::rk_static_ptr_cast<cylinder>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          cy_geom = rtti::rk_static_ptr_cast<cylinder>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a cylinder..
        if(other_geom->getObjectType() == cylinder::getStaticObjectType()) {
          shared_ptr<cylinder> cy2_geom = rtti::rk_static_ptr_cast<cylinder>(other_geom);
          //mProxFinders.push_back(shared_ptr< prox_cylinder_cylinder >(new prox_cylinder_cylinder(cy_geom, cy2_geom)));
        }
        // if the other is a box..
        else if(other_geom->getObjectType() == box::getStaticObjectType()) {
          shared_ptr<box> bx_geom = rtti::rk_static_ptr_cast<box>(other_geom);
          //mProxFinders.push_back(shared_ptr< prox_cylinder_box >(new prox_cylinder_box(cy_geom, bx_geom)));
        };
      }
      // if one of the model is a box?
      else if((mModel1->mShapeList[i]->getObjectType() == box::getStaticObjectType()) ||
              (mModel2->mShapeList[j]->getObjectType() == box::getStaticObjectType())) {
        shared_ptr<box> bx_geom;
        shared_ptr<shape_3D> other_geom;
        if(mModel1->mShapeList[i]->getObjectType() == box::getStaticObjectType()) {
          bx_geom = rtti::rk_static_ptr_cast<box>(mModel1->mShapeList[i]);
          other_geom = mModel2->mShapeList[j];
        } else {
          bx_geom = rtti::rk_static_ptr_cast<box>(mModel2->mShapeList[j]);
          other_geom = mModel1->mShapeList[i];
        };
        // if the other is a box..
        if(other_geom->getObjectType() == box::getStaticObjectType()) {
          shared_ptr<box> bx2_geom = rtti::rk_static_ptr_cast<box>(other_geom);
          //mProxFinders.push_back(shared_ptr< prox_box_box >(new prox_box_box(bx_geom, bx2_geom)));
        };
      };
    };
  };
  
  return;
};

shared_ptr< proximity_finder_3D > proxy_query_pair_3D::findMinimumDistance() const {
  if(mProxFinders.empty())
    return shared_ptr< proximity_finder_3D >();
  
  std::size_t min_i = 0;
  mProxFinders[0]->computeProximity();
  double min_dist = mProxFinders[0]->getLastResult().mDistance;
  
  for(std::size_t i = 1; i < mProxFinders.size(); ++i) {
    mProxFinders[i]->computeProximity();
    if(min_dist > mProxFinders[i]->getLastResult().mDistance) {
      min_i = i;
      min_dist = mProxFinders[i]->getLastResult().mDistance;
    };
  };
  
  return mProxFinders[min_i];
};
    
bool proxy_query_pair_3D::gatherCollisionPoints(std::vector< proximity_record_3D >& aOutput) const {
  bool collision_found = false;
  
  // TODO: Add a sphere-sphere collision test on each shape-pair before calling "computeProximity".
  
  for(std::size_t i = 0; i < mProxFinders.size(); ++i) {
    mProxFinders[i]->computeProximity();
    if(mProxFinders[i]->getLastResult().mDistance < 0.0) {
      aOutput.push_back(mProxFinders[i]->getLastResult());
      collision_found = true;
    };
  };
  
  return collision_found;
};



};


};





