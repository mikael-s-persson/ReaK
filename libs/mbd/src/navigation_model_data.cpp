
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include <ReaK/core/serialization/archiver_factory.hpp>
#include <ReaK/math/optimization/optim_exceptions.hpp>

#include <ReaK/mbd/models/navigation_model_data.hpp>

#include <ReaK/geometry/shapes/colored_model.hpp>
#include <ReaK/geometry/proximity/proxy_query_model.hpp>


namespace ReaK {

namespace kte {


navigation_scenario::navigation_scenario() : named_object() { this->setName( "navigation_scenario" ); };


void navigation_scenario::load_robot( const std::string& fileName ) {

  ( *serialization::open_iarchive( fileName ) ) >> robot_base_frame >> robot_kin_model >> robot_jt_limits
    >> robot_geom_model >> robot_proxy;

  robot_env_proxies.clear();
  create_robot_env_proxies();
};


void navigation_scenario::save_robot( const std::string& fileName ) const {

  ( *serialization::open_oarchive( fileName ) ) << robot_base_frame << robot_kin_model << robot_jt_limits
                                                << robot_geom_model << robot_proxy;
};


void navigation_scenario::load_environment( const std::string& fileName ) {

  shared_ptr< geom::proxy_query_model_3D > env_proxy;
  shared_ptr< geom::colored_model_3D > env_geom_model;

  ( *serialization::open_iarchive( fileName ) ) >> env_geom_model >> env_proxy;

  if( env_geom_model )
    env_geom_models.push_back( env_geom_model );
  if( env_proxy )
    env_proxy_models.push_back( env_proxy );

  create_robot_env_proxies();
};

void navigation_scenario::save_environment( std::size_t id, const std::string& fileName ) const {
  if( id >= env_geom_models.size() )
    return;

  ( *serialization::open_oarchive( fileName ) ) << env_geom_models[id] << env_proxy_models[id];
};


void navigation_scenario::clear_environment() {
  env_geom_models.clear();
  robot_env_proxies.clear();
};


void navigation_scenario::create_robot_env_proxies() {

  if( robot_proxy ) {

    for( std::size_t i = robot_env_proxies.size(); i < env_proxy_models.size(); ++i ) {
      robot_env_proxies.push_back( shared_ptr< geom::proxy_query_pair_3D >( new geom::proxy_query_pair_3D(
        "robot_env_proxy:" + env_proxy_models[i]->getName(), robot_proxy, env_proxy_models[i] ) ) );
    };
  };
};


void RK_CALL navigation_scenario::save( serialization::oarchive& A, unsigned int ) const {
  named_object::save( A, named_object::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_SAVE_WITH_NAME( robot_base_frame ) & RK_SERIAL_SAVE_WITH_NAME( robot_kin_model )
    & RK_SERIAL_SAVE_WITH_NAME( robot_jt_limits ) & RK_SERIAL_SAVE_WITH_NAME( robot_proxy )
    & RK_SERIAL_SAVE_WITH_NAME( robot_geom_model ) & RK_SERIAL_SAVE_WITH_NAME( target_frame )
    & RK_SERIAL_SAVE_WITH_NAME( env_geom_models ) & RK_SERIAL_SAVE_WITH_NAME( env_proxy_models )
    & RK_SERIAL_SAVE_WITH_NAME( robot_env_proxies );
};

void RK_CALL navigation_scenario::load( serialization::iarchive& A, unsigned int ) {
  named_object::load( A, named_object::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_LOAD_WITH_NAME( robot_base_frame ) & RK_SERIAL_LOAD_WITH_NAME( robot_kin_model )
    & RK_SERIAL_LOAD_WITH_NAME( robot_jt_limits ) & RK_SERIAL_LOAD_WITH_NAME( robot_proxy )
    & RK_SERIAL_LOAD_WITH_NAME( robot_geom_model ) & RK_SERIAL_LOAD_WITH_NAME( target_frame )
    & RK_SERIAL_LOAD_WITH_NAME( env_geom_models ) & RK_SERIAL_LOAD_WITH_NAME( env_proxy_models )
    & RK_SERIAL_LOAD_WITH_NAME( robot_env_proxies );
};
};
};
