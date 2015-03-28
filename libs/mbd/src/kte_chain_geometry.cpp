
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


#include <ReaK/mbd/kte/kte_chain_geometry.hpp>

#include <ReaK/mbd/kte/kte_map_chain.hpp>
#include <ReaK/mbd/kte/kte_chain_visitation.hpp>
#include <ReaK/mbd/kte/rigid_link.hpp>
#include <ReaK/mbd/kte/inertia.hpp>
#include <ReaK/mbd/kte/revolute_joint.hpp>
#include <ReaK/mbd/kte/prismatic_joint.hpp>

#include <string>
#include <map>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace kte {


namespace detail {

typedef boost::tuple< rigid_link_2D, prismatic_joint_2D, revolute_joint_2D, inertia_2D > kte_geom_visit_2D_types;
typedef boost::tuple< rigid_link_3D, prismatic_joint_3D, revolute_joint_3D, inertia_3D > kte_geom_visit_3D_types;


struct kte_geom_2D_visitor {
  const kte_chain_geometry_2D* parent;

  kte_geom_2D_visitor( const kte_chain_geometry_2D* aParent ) : parent( aParent ){};

  template < typename KTEType >
  void connect_to_base_or_end( const KTEType& aObj ) {
    std::string obj_name = aObj.getName();

    std::map< std::string, std::vector< geom::colored_geometry_2D > >::const_iterator it
      = parent->mGeomList.find( obj_name + "_base" );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom )
          it->second[i].mGeom->setAnchor( aObj.BaseFrame() );
      };
    };

    std::map< std::string, std::vector< shared_ptr< geom::shape_2D > > >::const_iterator it_prox
      = parent->mProxyShapeList.find( obj_name + "_base" );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] )
          it_prox->second[i]->setAnchor( aObj.BaseFrame() );
      };
    };

    it = parent->mGeomList.find( obj_name + "_end" );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom )
          it->second[i].mGeom->setAnchor( aObj.EndFrame() );
      };
    };

    it_prox = parent->mProxyShapeList.find( obj_name + "_end" );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] )
          it_prox->second[i]->setAnchor( aObj.EndFrame() );
      };
    };
  };

  void operator()( const rigid_link_2D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const prismatic_joint_2D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const revolute_joint_2D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const inertia_2D& aObj ) {
    std::string obj_name = aObj.getName();

    std::map< std::string, std::vector< geom::colored_geometry_2D > >::const_iterator it
      = parent->mGeomList.find( obj_name );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom )
          it->second[i].mGeom->setAnchor( aObj.CenterOfMass()->mFrame );
      };
    };

    std::map< std::string, std::vector< shared_ptr< geom::shape_2D > > >::const_iterator it_prox
      = parent->mProxyShapeList.find( obj_name );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] )
          it_prox->second[i]->setAnchor( aObj.CenterOfMass()->mFrame );
      };
    };
  };
};


struct kte_geom_3D_visitor {
  const kte_chain_geometry_3D* parent;

  kte_geom_3D_visitor( const kte_chain_geometry_3D* aParent ) : parent( aParent ){};

  template < typename KTEType >
  void connect_to_base_or_end( const KTEType& aObj ) {
    std::string obj_name = aObj.getName();

    std::map< std::string, std::vector< geom::colored_geometry_3D > >::const_iterator it
      = parent->mGeomList.find( obj_name + "_base" );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom ) {
          it->second[i].mGeom->setAnchor( aObj.BaseFrame() );
        };
      };
    };

    std::map< std::string, std::vector< shared_ptr< geom::shape_3D > > >::const_iterator it_prox
      = parent->mProxyShapeList.find( obj_name + "_base" );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] ) {
          it_prox->second[i]->setAnchor( aObj.BaseFrame() );
        };
      };
    };

    it = parent->mGeomList.find( obj_name + "_end" );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom ) {
          it->second[i].mGeom->setAnchor( aObj.EndFrame() );
        };
      };
    };

    it_prox = parent->mProxyShapeList.find( obj_name + "_end" );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] ) {
          it_prox->second[i]->setAnchor( aObj.EndFrame() );
        };
      };
    };
  };

  void operator()( const rigid_link_3D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const prismatic_joint_3D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const revolute_joint_3D& aObj ) { connect_to_base_or_end( aObj ); };

  void operator()( const inertia_3D& aObj ) {
    std::string obj_name = aObj.getName();

    std::map< std::string, std::vector< geom::colored_geometry_3D > >::const_iterator it
      = parent->mGeomList.find( obj_name );
    if( ( it != parent->mGeomList.end() ) && ( it->second.size() ) ) {
      for( std::size_t i = 0; i < it->second.size(); ++i ) {
        if( it->second[i].mGeom )
          it->second[i].mGeom->setAnchor( aObj.CenterOfMass()->mFrame );
      };
    };

    std::map< std::string, std::vector< shared_ptr< geom::shape_3D > > >::const_iterator it_prox
      = parent->mProxyShapeList.find( obj_name );
    if( ( it_prox != parent->mProxyShapeList.end() ) && ( it_prox->second.size() ) ) {
      for( std::size_t i = 0; i < it_prox->second.size(); ++i ) {
        if( it_prox->second[i] )
          it_prox->second[i]->setAnchor( aObj.CenterOfMass()->mFrame );
      };
    };
  };
};
};


std::pair< shared_ptr< geom::colored_model_2D >, shared_ptr< geom::proxy_query_model_2D > >
  kte_chain_geometry_2D::attachToKTEChain( const kte_map_chain& aKTEChain ) const {
  // first, visit the KTE chain to associated each geometry to its KTE-related anchor.
  detail::kte_geom_2D_visitor vis = detail::kte_geom_2D_visitor( this );
  visit_kte_chain< detail::kte_geom_visit_2D_types >( vis, aKTEChain );

  // then, construct the colored_model_2D from the resulting (attached) geometries:
  shared_ptr< geom::colored_model_2D > result_geom = shared_ptr< geom::colored_model_2D >(
    new geom::colored_model_2D( aKTEChain.getName() + "_geom" ), scoped_deleter() );
  for( std::map< std::string, std::vector< geom::colored_geometry_2D > >::const_iterator it = mGeomList.begin();
       it != mGeomList.end(); ++it ) {
    for( std::size_t i = 0; i < it->second.size(); ++i ) {
      if( it->second[i].mGeom ) {
        result_geom->addAnchor( it->second[i].mGeom->getAnchor() );
        result_geom->addElement( it->second[i].mColor, it->second[i].mGeom );
      };
    };
  };

  // then, construct the proxy_query_model_2D from the resulting (attached) geometries:
  shared_ptr< geom::proxy_query_model_2D > result_prox = shared_ptr< geom::proxy_query_model_2D >(
    new geom::proxy_query_model_2D( aKTEChain.getName() + "_proxy" ), scoped_deleter() );
  for( std::map< std::string, std::vector< shared_ptr< geom::shape_2D > > >::const_iterator it
       = mProxyShapeList.begin();
       it != mProxyShapeList.end(); ++it ) {
    for( std::size_t i = 0; i < it->second.size(); ++i ) {
      if( it->second[i] ) {
        result_prox->addShape( it->second[i] );
      };
    };
  };

  return std::pair< shared_ptr< geom::colored_model_2D >, shared_ptr< geom::proxy_query_model_2D > >( result_geom,
                                                                                                      result_prox );
};


std::pair< shared_ptr< geom::colored_model_3D >, shared_ptr< geom::proxy_query_model_3D > >
  kte_chain_geometry_3D::attachToKTEChain( const kte_map_chain& aKTEChain ) const {
  detail::kte_geom_3D_visitor vis = detail::kte_geom_3D_visitor( this );
  visit_kte_chain< detail::kte_geom_visit_3D_types >( vis, aKTEChain );

  // then, construct the colored_model_3D from the resulting (attached) geometries:
  shared_ptr< geom::colored_model_3D > result_geom = shared_ptr< geom::colored_model_3D >(
    new geom::colored_model_3D( aKTEChain.getName() + "_geom" ), scoped_deleter() );
  for( std::map< std::string, std::vector< geom::colored_geometry_3D > >::const_iterator it = mGeomList.begin();
       it != mGeomList.end(); ++it ) {
    for( std::size_t i = 0; i < it->second.size(); ++i ) {
      if( it->second[i].mGeom ) {
        result_geom->addAnchor( it->second[i].mGeom->getAnchor() );
        result_geom->addElement( it->second[i].mColor, it->second[i].mGeom );
      };
    };
  };

  // then, construct the proxy_query_model_3D from the resulting (attached) geometries:
  shared_ptr< geom::proxy_query_model_3D > result_prox = shared_ptr< geom::proxy_query_model_3D >(
    new geom::proxy_query_model_3D( aKTEChain.getName() + "_proxy" ), scoped_deleter() );
  for( std::map< std::string, std::vector< shared_ptr< geom::shape_3D > > >::const_iterator it
       = mProxyShapeList.begin();
       it != mProxyShapeList.end(); ++it ) {
    for( std::size_t i = 0; i < it->second.size(); ++i ) {
      if( it->second[i] ) {
        result_prox->addShape( it->second[i] );
      };
    };
  };

  return std::pair< shared_ptr< geom::colored_model_3D >, shared_ptr< geom::proxy_query_model_3D > >( result_geom,
                                                                                                      result_prox );
};
};
};
