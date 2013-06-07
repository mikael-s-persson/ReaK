
/*
 *    Copyright 2013 Sven Mikael Persson
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


#include "kte_chain_geometry.hpp"


#include "mbd_kte/kte_chain_visitation.hpp"

#include "mbd_kte/rigid_link.hpp"
#include "mbd_kte/inertia.hpp"

#include <string>
#include <map>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


namespace detail {
  
  typedef boost::tuple< kte::rigid_link_2D, kte::inertia_2D > kte_geom_visit_2D_types;
  typedef boost::tuple< kte::rigid_link_3D, kte::inertia_3D > kte_geom_visit_3D_types;
  
  
  struct kte_geom_2D_visitor {
    const kte_chain_geometry_2D* parent;
    
    kte_geom_2D_visitor(const kte_chain_geometry_2D* aParent) : parent(aParent) { };
    
    void operator()(const kte::rigid_link_2D& aObj) {
      std::string obj_name = aObj.getName();
      std::map< std::string, colored_geometry_2D >::const_iterator it = parent->mGeomList.find(obj_name);
      if((it != parent->mGeomList.end()) && (it->second.mGeom))
        it->second.mGeom->setAnchor(aObj.BaseFrame());
    };
    
    void operator()(const kte::inertia_2D& aObj) {
      std::string obj_name = aObj.getName();
      std::map< std::string, colored_geometry_2D >::const_iterator it = parent->mGeomList.find(obj_name);
      if((it != parent->mGeomList.end()) && (it->second.mGeom))
        it->second.mGeom->setAnchor(aObj.CenterOfMass()->mFrame);
    };
    
  };
  
  
  struct kte_geom_3D_visitor {
    const kte_chain_geometry_3D* parent;
    
    kte_geom_3D_visitor(const kte_chain_geometry_3D* aParent) : parent(aParent) { };
    
    void operator()(const kte::rigid_link_3D& aObj) {
      std::string obj_name = aObj.getName();
      std::map< std::string, colored_geometry_3D >::const_iterator it = parent->mGeomList.find(obj_name);
      if((it != parent->mGeomList.end()) && (it->second.mGeom))
        it->second.mGeom->setAnchor(aObj.BaseFrame());
    };
    
    void operator()(const kte::inertia_3D& aObj) {
      std::string obj_name = aObj.getName();
      std::map< std::string, colored_geometry_3D >::const_iterator it = parent->mGeomList.find(obj_name);
      if((it != parent->mGeomList.end()) && (it->second.mGeom))
        it->second.mGeom->setAnchor(aObj.CenterOfMass()->mFrame);
    };
    
  };
  
  
};


shared_ptr< colored_model_2D > kte_chain_geometry_2D::attachGeomToKTEChain(const kte::kte_map_chain& aKTEChain) const {
  // first, visit the KTE chain to associated each geometry to its KTE-related anchor.
  detail::kte_geom_2D_visitor vis = detail::kte_geom_2D_visitor(this);
  kte::visit_kte_chain< detail::kte_geom_visit_2D_types >(vis, aKTEChain);
  
  // then, construct the colored_model_2D from the resulting (attached) geometries:
  shared_ptr< colored_model_2D > result = shared_ptr< colored_model_2D >(new colored_model_2D(aKTEChain.getName() + "_geom"), scoped_deleter());
  for(std::map< std::string, colored_geometry_2D >::const_iterator it = mGeomList.begin(); it != mGeomList.end(); ++it) {
    if(it->second.mGeom) {
      result->addAnchor(it->second.mGeom->getAnchor());
      result->addElement(it->second.mColor, it->second.mGeom);
    };
  };
  return result;
};


shared_ptr< proxy_query_model_2D > kte_chain_geometry_2D::attachProxyToKTEChain(const kte::kte_map_chain& aKTEChain) const {
  // first, visit the KTE chain to associated each geometry to its KTE-related anchor.
  detail::kte_geom_2D_visitor vis = detail::kte_geom_2D_visitor(this);
  kte::visit_kte_chain< detail::kte_geom_visit_2D_types >(vis, aKTEChain);
  
  // then, construct the proxy_query_model_2D from the resulting (attached) geometries:
  shared_ptr< proxy_query_model_2D > result = shared_ptr< proxy_query_model_2D >(new proxy_query_model_2D(aKTEChain.getName() + "_geom"), scoped_deleter());
  for(std::map< std::string, colored_geometry_2D >::const_iterator it = mGeomList.begin(); it != mGeomList.end(); ++it) {
    if(it->second.mGeom) {
      shared_ptr< shape_2D > p_shape = rtti::rk_dynamic_ptr_cast< shape_2D >(it->second.mGeom);
      if(p_shape)
        result->addShape(p_shape);
    };
  };
  return result;
};




shared_ptr< colored_model_3D > kte_chain_geometry_3D::attachGeomToKTEChain(const kte::kte_map_chain& aKTEChain) const {
  detail::kte_geom_3D_visitor vis = detail::kte_geom_3D_visitor(this);
  kte::visit_kte_chain< detail::kte_geom_visit_3D_types >(vis, aKTEChain);
  
  // then, construct the colored_model_3D from the resulting (attached) geometries:
  shared_ptr< colored_model_3D > result = shared_ptr< colored_model_3D >(new colored_model_3D(aKTEChain.getName() + "_geom"), scoped_deleter());
  for(std::map< std::string, colored_geometry_3D >::const_iterator it = mGeomList.begin(); it != mGeomList.end(); ++it) {
    if(it->second.mGeom) {
      result->addAnchor(it->second.mGeom->getAnchor());
      result->addElement(it->second.mColor, it->second.mGeom);
    };
  };
  return result;
};


shared_ptr< proxy_query_model_3D > kte_chain_geometry_3D::attachProxyToKTEChain(const kte::kte_map_chain& aKTEChain) const {
  detail::kte_geom_3D_visitor vis = detail::kte_geom_3D_visitor(this);
  kte::visit_kte_chain< detail::kte_geom_visit_3D_types >(vis, aKTEChain);
  
  // then, construct the proxy_query_model_3D from the resulting (attached) geometries:
  shared_ptr< proxy_query_model_3D > result = shared_ptr< proxy_query_model_3D >(new proxy_query_model_3D(aKTEChain.getName() + "_geom"), scoped_deleter());
  for(std::map< std::string, colored_geometry_3D >::const_iterator it = mGeomList.begin(); it != mGeomList.end(); ++it) {
    if(it->second.mGeom) {
      shared_ptr< shape_3D > p_shape = rtti::rk_dynamic_ptr_cast< shape_3D >(it->second.mGeom);
      if(p_shape)
        result->addShape(p_shape);
    };
  };
  return result;
};



};

};









