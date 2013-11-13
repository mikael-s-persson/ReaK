
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

#include "chaser_target_model_data.hpp"

#include "shapes/colored_model.hpp"
#include "proximity/proxy_query_model.hpp"

#include "serialization/archiver_factory.hpp"

#include "optimization/optim_exceptions.hpp"

namespace ReaK {
  
namespace kte {

  
chaser_target_data::chaser_target_data() : named_object() {
  this->setName("chaser_target_model_data");
};


void chaser_target_data::load_chaser(const std::string& fileName) {
  
  (*serialization::open_iarchive(fileName))
    >> chaser_base_frame 
    >> chaser_kin_model 
    >> chaser_jt_limits 
    >> chaser_geom_model 
    >> chaser_proxy;
  
  create_chaser_target_proxy();
  chaser_env_proxies.clear();
  create_chaser_env_proxies();
};


void chaser_target_data::save_chaser(const std::string& fileName) const {
  
  (*serialization::open_oarchive(fileName))
    << chaser_base_frame 
    << chaser_kin_model 
    << chaser_jt_limits 
    << chaser_geom_model 
    << chaser_proxy;
  
};



void chaser_target_data::load_target(const std::string& fileName) {
  
  shared_ptr< frame_3D<double> > target_base;
  
  (*serialization::open_iarchive(fileName))
    >> target_base
    >> target_kin_model
    >> target_frame
    >> target_geom_model
    >> target_proxy;
  
  create_chaser_target_proxy();
  target_env_proxies.clear();
  create_target_env_proxies();
};

void chaser_target_data::save_target(const std::string& fileName) const {
  
  (*serialization::open_oarchive(fileName))
    << shared_ptr< frame_3D<double> >()
    << target_kin_model
    << target_frame
    << target_geom_model
    << target_proxy;
  
};



void chaser_target_data::load_environment(const std::string& fileName) {
  
  shared_ptr< geom::proxy_query_model_3D > env_proxy;
  shared_ptr< geom::colored_model_3D > env_geom_model;
  
  (*serialization::open_iarchive(fileName)) >> env_geom_model >> env_proxy;
  
  if( env_geom_model )
    env_geom_models.push_back(env_geom_model);
  if( env_proxy )
    env_proxy_models.push_back(env_proxy);
  
  create_chaser_env_proxies();
  create_target_env_proxies();
  
};

void chaser_target_data::save_environment(std::size_t id, const std::string& fileName) const {
  if(id >= env_geom_models.size())
    return;
  
  (*serialization::open_oarchive(fileName)) << env_geom_models[id] << env_proxy_models[id];
  
};


void chaser_target_data::clear_environment() {
  env_geom_models.clear();
  chaser_env_proxies.clear();
  target_env_proxies.clear();
};


void chaser_target_data::create_chaser_target_proxy() {
  
  if( chaser_proxy && target_proxy ) {
    chaser_target_proxy = shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D("chaser_target_proxy", chaser_proxy, target_proxy));
  } else if( chaser_target_proxy ) {
    chaser_target_proxy.reset();
  };
  
};


void chaser_target_data::create_chaser_env_proxies() {
  
  if( chaser_proxy ) {
    
    for(std::size_t i = chaser_env_proxies.size(); i < env_proxy_models.size(); ++i) {
      chaser_env_proxies.push_back(
        shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D(
          "chaser_env_proxy:" + env_proxy_models[i]->getName(), 
          chaser_proxy, 
          env_proxy_models[i])));
    };
    
  };
  
};


void chaser_target_data::create_target_env_proxies() {
  
  if( target_proxy ) {
    
    for(std::size_t i = target_env_proxies.size(); i < env_proxy_models.size(); ++i) {
      target_env_proxies.push_back(
        shared_ptr< geom::proxy_query_pair_3D >(new geom::proxy_query_pair_3D(
          "target_env_proxy:" + env_proxy_models[i]->getName(), 
          target_proxy, 
          env_proxy_models[i])));
    };
    
  };
  
};



void RK_CALL chaser_target_data::save(serialization::oarchive& A, unsigned int) const {
  named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_SAVE_WITH_NAME(chaser_base_frame)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_kin_model)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_jt_limits)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_proxy)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_geom_model)
    & RK_SERIAL_SAVE_WITH_NAME(target_kin_model)
    & RK_SERIAL_SAVE_WITH_NAME(target_frame)
    & RK_SERIAL_SAVE_WITH_NAME(target_proxy)
    & RK_SERIAL_SAVE_WITH_NAME(target_geom_model)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_target_proxy)
    & RK_SERIAL_SAVE_WITH_NAME(env_geom_models)
    & RK_SERIAL_SAVE_WITH_NAME(env_proxy_models)
    & RK_SERIAL_SAVE_WITH_NAME(chaser_env_proxies)
    & RK_SERIAL_SAVE_WITH_NAME(target_env_proxies);
};

void RK_CALL chaser_target_data::load(serialization::iarchive& A, unsigned int) {
  named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
  A & RK_SERIAL_LOAD_WITH_NAME(chaser_base_frame)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_kin_model)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_jt_limits)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_proxy)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_geom_model)
    & RK_SERIAL_LOAD_WITH_NAME(target_kin_model)
    & RK_SERIAL_LOAD_WITH_NAME(target_frame)
    & RK_SERIAL_LOAD_WITH_NAME(target_proxy)
    & RK_SERIAL_LOAD_WITH_NAME(target_geom_model)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_target_proxy)
    & RK_SERIAL_LOAD_WITH_NAME(env_geom_models)
    & RK_SERIAL_LOAD_WITH_NAME(env_proxy_models)
    & RK_SERIAL_LOAD_WITH_NAME(chaser_env_proxies)
    & RK_SERIAL_LOAD_WITH_NAME(target_env_proxies);
};



};

};

