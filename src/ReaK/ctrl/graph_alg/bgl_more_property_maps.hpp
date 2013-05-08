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
#include <boost/property_map/property_map.hpp>

#include "bgl_more_property_tags.hpp"


namespace boost {



template <typename Reference, typename LvaluePropertyMap>
struct subobject_put_get_helper { };

template <typename PropertyMap, typename Reference, typename K>
inline
const typename property_traits<PropertyMap>::value_type& get(const subobject_put_get_helper<Reference, PropertyMap>& pa, const K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
};

template <typename PropertyMap, typename Reference, typename K>
inline
Reference get(const subobject_put_get_helper<Reference, PropertyMap>& pa, K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
};

#ifndef RK_ENABLE_CXX0X_FEATURES

template <typename PropertyMap, typename Reference, typename K, typename V>
inline
void put(const subobject_put_get_helper<Reference, PropertyMap>& pa, K& k, const V& v) {
  static_cast<const PropertyMap&>(pa)[k] = v;
};

#else

template <typename PropertyMap, typename Reference, typename K, typename V>
inline
void put(const subobject_put_get_helper<Reference, PropertyMap>& pa, K& k, V&& v) {
  static_cast<const PropertyMap&>(pa)[k] = std::forward<V>(v);
};

#endif


template <typename Graph, typename PropertyMapTag>
struct whole_bundle_property_map :
    public put_get_helper<
      typename mpl::if_< is_same< PropertyMapTag, vertex_bundle_t>,
        typename Graph::vertex_bundled,
	typename Graph::edge_bundled >::type&,
      whole_bundle_property_map<Graph,PropertyMapTag> > {
  private:
    Graph* pg;
  public:
    typedef is_same< PropertyMapTag, vertex_bundle_t> is_vertex_bundle;
    typedef is_const< Graph > is_const_graph;
    typedef typename mpl::if_< is_vertex_bundle,
      typename Graph::vertex_bundled,
      typename Graph::edge_bundled >::type value_type;
    typedef typename mpl::if_< is_const_graph,
      const value_type&,
      value_type& >::type reference;
    typedef typename mpl::if_< is_vertex_bundle,
      typename graph_traits<Graph>::vertex_descriptor,
      typename graph_traits<Graph>::edge_descriptor >::type key_type;
    typedef typename mpl::if_< is_const_graph,
      readable_property_map_tag,
      lvalue_property_map_tag >::type category;

    whole_bundle_property_map(Graph* aPG) : pg(aPG) { };
    reference operator[](key_type k) const { return (*pg)[k]; };

};


template <typename T, typename Graph, typename PropertyMapTag>
struct propgraph_property_map :
    public put_get_helper< T&, propgraph_property_map<T, Graph, PropertyMapTag> > {
  private:
    Graph* pg;
    PropertyMapTag tag;
  public:
    typedef is_same< typename property_kind<PropertyMapTag>::type, vertex_property_tag> is_vertex_prop;
    typedef is_const< Graph > is_const_graph;
    typedef T value_type;
    typedef T& reference;
    typedef typename mpl::if_< is_vertex_prop,
      typename graph_traits<Graph>::vertex_descriptor,
      typename graph_traits<Graph>::edge_descriptor >::type key_type;
    typedef typename mpl::if_< is_const< T >,
      readable_property_map_tag,
      lvalue_property_map_tag >::type category;

    propgraph_property_map(Graph* aPG, PropertyMapTag aTag) : pg(aPG), tag(aTag) { };
    reference operator[](key_type k) const { return get(tag,*pg,k); };

};


template <typename T>
struct self_property_map :
    public subobject_put_get_helper<T&, self_property_map<T> > {
  typedef self_property_map self;

  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T key_type;
  typedef typename mpl::if_< is_const<T>,
    readable_property_map_tag,
    lvalue_property_map_tag >::type category;

  self_property_map() { };
  reference operator[](reference p) const { return p; };
  const_reference operator[](const_reference p) const { return p; };
};


template <typename T, typename PropertyType>
class data_member_property_map :
    public subobject_put_get_helper<T&, data_member_property_map<T, PropertyType> > {
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
    typedef typename mpl::if_< is_const<T>,
      readable_property_map_tag,
      lvalue_property_map_tag >::type category;

    data_member_property_map(member_ptr_type aMemPtr) : mem_ptr(aMemPtr) { };
    reference operator[](key_type& p) const { return p.*mem_ptr; };
    const_reference operator[](const key_type& p) const { return p.*mem_ptr; };
};

template <typename T, typename PropertyType>
class data_member_property_map< const T, const PropertyType> :
    public subobject_put_get_helper<const T&, data_member_property_map<const T, const PropertyType> > {
  public:
    typedef T value_type;
    typedef const T& reference;
    typedef const T& const_reference;
    typedef const PropertyType key_type;

    typedef T key_type::* member_ptr_type;
    typedef data_member_property_map<const T, const PropertyType> self;
  private:
    member_ptr_type mem_ptr;
  public:
    typedef typename mpl::if_< is_const<T>,
      readable_property_map_tag,
      lvalue_property_map_tag >::type category;

    data_member_property_map(member_ptr_type aMemPtr) : mem_ptr(aMemPtr) { };
    reference operator[](key_type& p) const { return p.*mem_ptr; };
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
    typedef typename property_traits< OutputMap >::reference reference;
    typedef const reference const_reference;

    composite_property_map(OutputMap aPropOut, InputMap aPropIn) : prop_out(aPropOut), prop_in(aPropIn) { };

    reference operator[](const key_type& k) const {
      return prop_out[ prop_in[k] ];
    };

    reference operator[](key_type& k) const {
      return prop_out[ prop_in[k] ];
    };

    friend
    value_type get(const self& m, const key_type& p) {
      return m.prop_out[m.prop_in[p]];
    };

#ifndef RK_ENABLE_CXX0X_FEATURES
    template <typename V>
    friend
    void put(const self& m, const key_type& p, const V& value) {
      put(m.prop_out, m.prop_in[p], value);
    };

    template <typename V>
    friend
    void put(const self& m, key_type& p, const V& value) {
      put(m.prop_out, m.prop_in[p], value);
    };

#else

    template <typename V>
    friend
    void put(const self& m, const key_type& p, V&& value) {
      put(m.prop_out, m.prop_in[p], std::forward<V>(value));
    };

    template <typename V>
    friend
    void put(const self& m, key_type& p, V&& value) {
      put(m.prop_out, m.prop_in[p], std::forward<V>(value));
    };
#endif
};



template <typename T, typename Graph, typename PropertyMapTag>
class bundle_member_property_map :
    public put_get_helper< T&, bundle_member_property_map<T, Graph, PropertyMapTag> > {
  public:
    typedef bundle_member_property_map<T, Graph, PropertyMapTag> self;
    typedef is_same< PropertyMapTag, vertex_bundle_t> is_vertex_bundle;
    typedef is_const< Graph > is_const_graph;
    typedef typename mpl::if_< is_vertex_bundle,
      typename Graph::vertex_bundled,
      typename Graph::edge_bundled >::type bundle_type;
    typedef T bundle_type::* member_ptr_type;
  private:
    Graph* pg;
    member_ptr_type mem_ptr;
  public:
    typedef T value_type;
    typedef T& reference;
    typedef typename mpl::if_< is_vertex_bundle,
      typename graph_traits<Graph>::vertex_descriptor,
      typename graph_traits<Graph>::edge_descriptor >::type key_type;
    typedef typename mpl::if_< is_const_graph,
      readable_property_map_tag,
      lvalue_property_map_tag >::type category;

    bundle_member_property_map(Graph* aPG, member_ptr_type aMemPtr) : pg(aPG), mem_ptr(aMemPtr) { };
    reference operator[](key_type p) const { return (*pg)[p].*mem_ptr; };
};





template <typename BundleMemberMap, typename Graph>
composite_property_map< BundleMemberMap, whole_bundle_property_map< Graph, vertex_bundle_t > > 
  bundle_prop_to_vertex_prop(BundleMemberMap bundle_prop, Graph& g) {
  typedef composite_property_map< BundleMemberMap, whole_bundle_property_map< Graph, vertex_bundle_t > > ResultType;
  return ResultType(bundle_prop, whole_bundle_property_map< Graph, vertex_bundle_t >(&g));
};      


};


#endif


















