/**
 * \file bgl_more_property_maps.hpp
 * 
 * This library provides additional property-maps for Boost Graph Library's property maps.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2012
 */

/*
 *    Copyright 2011 Sven Mikael Persson
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

#ifndef REAK_BGL_MORE_PROPERTY_MAPS_HPP
#define REAK_BGL_MORE_PROPERTY_MAPS_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

    
namespace boost {


template <typename T, typename PropertyType>
class data_member_property_map {
  public:
    typedef T PropertyType::* member_ptr_type;
    typedef data_member_property_map<T,PropertyType> self;
  private:
    member_ptr_type mem_ptr;
  public:
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
    typedef PropertyType key_type;
    typedef lvalue_property_map_tag category;
    
    data_member_property_map(member_ptr_type aMemPtr) : mem_ptr(aMemPtr) { };
    
    const_reference operator[](const key_type& p) const {
      return p.*mem_ptr;
    };
    
    reference operator[](key_type& p) const {
      return p.*mem_ptr;
    };
    
    friend
    const_reference get(const self& m, const key_type& p) {
      return p.*(m.mem_ptr);
    };
    
    friend
    reference get(const self& m, key_type& p) {
      return p.*(m.mem_ptr);
    };
    
    friend
    void put(const self& m, key_type& p, const_reference value) {
      p.*(m.mem_ptr) = value;
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend
    void put(const self& m, key_type& p, value_type&& value) {
      p.*(m.mem_ptr) = std::move(value);
    };
#endif
  
  
};


template <typename OutputMap, typename InputMap>
class composite_property_map {
  public:
    typedef composite_property_map<OutputMap,InputMap> self;
    
  private:
    OutputMap prop_out;
    InputMap prop_in;
  public:
    typedef typename property_traits< OutputMap >::value_type value_type;
    typedef typename property_traits< InputMap >::key_type key_type;
    typedef typename property_traits< OutputMap >::category category;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    
    composite_property_map(OutputMap aPropOut, InnerMap aPropIn) : prop_out(aPropOut), prop_in(aPropIn) { };
    
    reference operator[](const key_type& k) {
      return prop_out[ get(prop_in, k) ];
    };
    
    reference operator[](key_type& k) {
      return prop_out[ get(prop_in, k) ];
    };
    
    friend
    const_reference get(self& m, const key_type& p) {
      return get(prop_out, get(prop_in, p));
    };
    
    friend
    void put(self& m, const key_type& p, const_reference value) {
      put(prop_out, get(prop_in, p), value);
    };
    
    friend
    void put(self& m, key_type& p, const_reference value) {
      put(prop_out, get(prop_in, p), value);
    };
    
#ifdef RK_ENABLE_CXX0X_FEATURES
    friend
    void put(self& m, const key_type& p, value_type&& value) {
      put(prop_out, get(prop_in, p), std::move(value));
    };
    
    friend
    void put(self& m, key_type& p, value_type&& value) {
      put(prop_out, get(prop_in, p), std::move(value));
    };
#endif
};




};


#endif


















