/**
 * \file adj_list_tree_overlay.hpp
 *
 * This library provides a class
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
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

#ifndef REAK_ADJ_LIST_TREE_OVERLAY_HPP
#define REAK_ADJ_LIST_TREE_OVERLAY_HPP

#include <ReaK/core/base/defs.hpp>

#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>

#include <vector>
#include <queue>
#include <algorithm>
#include <functional>


// BGL-Extra includes:
#include <boost/graph/more_property_tags.hpp>
#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/graph/bfl_d_ary_tree.hpp>
//#include <boost/graph/tree_adaptor.hpp>
#include <boost/graph/more_property_maps.hpp>

// Pending inclusion in BGL-Extra:
#include "bgl_raw_property_graph.hpp"

#include "tree_organizer_concept.hpp"



namespace ReaK {

namespace graph {


template <typename AdjListOnTreeType>
class alt_tree_view;  // forward-declaration.


template <typename AdjListOnTreeType>
class alt_graph_view;  // forward-declaration.


namespace detail { namespace {


template <typename TreeVertexProperty, typename AdjVertexProperty>
struct tree_vertex_properties_impl {
  typedef tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> self;
  std::size_t partner_node;
  TreeVertexProperty tree_data;
  AdjVertexProperty user_data;

  friend bool operator <(const self& lhs, const self& rhs) {
     return lhs.partner_node < rhs.partner_node;
  };
  friend bool operator <=(const self& lhs, const self& rhs) {
    return lhs.partner_node <= rhs.partner_node;
  };
  friend bool operator >(const self& lhs, const self& rhs) {
    return lhs.partner_node > rhs.partner_node;
  };
  friend bool operator >=(const self& lhs, const self& rhs) {
    return lhs.partner_node >= rhs.partner_node;
  };
  friend bool operator ==(const self& lhs, const self& rhs) {
    return lhs.partner_node == rhs.partner_node;
  };
  friend bool operator !=(const self& lhs, const self& rhs) {
    return lhs.partner_node != rhs.partner_node;
  };
  /*
  friend
  boost::data_member_property_map<const std::size_t, const self>
    get(boost::vertex_key_t, const self&) {
    return boost::data_member_property_map<const std::size_t, const self>(&self::partner_node);
  };

  friend
  boost::data_member_property_map<AdjVertexProperty, self>
    get(boost::vertex_second_bundle_t, const self&) {
    return boost::data_member_property_map<AdjVertexProperty, self>(&self::user_data);
  };*/
};


template <typename TreeVertexProperty, typename AdjVertexProperty>
boost::data_member_property_map<const std::size_t, const tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> >
 get(boost::vertex_key_t, const tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>&) {
  return boost::data_member_property_map<const std::size_t, const tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> >(&tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>::partner_node);
};

template <typename TreeVertexProperty, typename AdjVertexProperty>
boost::data_member_property_map<AdjVertexProperty, tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> >
 get(boost::vertex_second_bundle_t, const tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>&) {
  return boost::data_member_property_map<AdjVertexProperty, tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> >(&tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>::user_data);
};


/*
 * This class template implements an adjacency-list (in BGL-style) that over-shadows (i.e., a multi-graph)
 * a tree to store the underlying vertices according to the layout that the tree uses (e.g., B-tree (BFL),
 * or COB-tree). An additional benefit of this class template is that the B-tree can act without being an
 * indirect to a graph, and the two graphs are always kept consistent, which would be difficult to do externally.
 */
template <typename AdjListType,
          typename VertexProperty,
          typename TreeVertexProperty,
          typename TreeEdgeProperty,
          typename TreeStorageTag>
class adj_list_on_tree {
  public:
    typedef adj_list_on_tree<AdjListType, VertexProperty, TreeVertexProperty, TreeEdgeProperty, TreeStorageTag> self;
    
    typedef tree_vertex_properties_impl<TreeVertexProperty, VertexProperty> tree_vertex_properties;
    typedef TreeEdgeProperty tree_edge_properties;
    
    typedef typename boost::tree_storage< tree_vertex_properties, tree_edge_properties, TreeStorageTag>::type tree_type;
    typedef TreeVertexProperty tree_vertex_bundled;
    typedef TreeEdgeProperty tree_edge_bundled;
    
    typedef AdjListType adj_list_type;
    typedef std::size_t adj_vertex_descriptor;
    
    typedef VertexProperty adj_vertex_bundled;
    typedef typename adj_list_type::edge_bundled adj_edge_bundled;
    
  private:
    adj_list_type m_adj_list;
    tree_type m_tree;
    
  public:
    friend class alt_tree_view< self >;
    friend class alt_graph_view< self >;

};

}; }; // detail






template <typename OutEdgeListS = boost::vecBC,
          typename DirectedS = boost::directedS,
          typename VertexProperty = boost::no_property,
          typename AdjEdgeProperty = boost::no_property,
          typename AdjEdgeListS = boost::vecBC,
          typename TreeStorageTag = boost::bfl_d_ary_tree_storage<2> >
struct adj_list_on_tree_tag {
  
  template <typename TreeVertexProperty, typename TreeEdgeProperty>
  struct alt {
    typedef detail::tree_vertex_properties_impl<TreeVertexProperty, VertexProperty> tree_vertex_properties;
    typedef TreeEdgeProperty tree_edge_properties;
    typedef typename boost::tree_storage< tree_vertex_properties, tree_edge_properties, TreeStorageTag>::type tree_type;
    
    typedef typename boost::graph_traits<tree_type>::vertex_descriptor tree_vertex_descriptor;
    
    struct adj_vertex_properties {
      tree_vertex_descriptor tree_vertex;
    };
    typedef AdjEdgeProperty adj_edge_properties;
    
    typedef boost::adjacency_list_BC< OutEdgeListS, boost::poolBC, DirectedS, adj_vertex_properties, adj_edge_properties > adj_list_type;
    
    typedef detail::adj_list_on_tree<adj_list_type, VertexProperty, TreeVertexProperty, TreeEdgeProperty, TreeStorageTag> type;
  };
  
};


}; // namespace graph

}; // namespace ReaK


namespace boost {

template <typename TreeVertexProperty,
          typename TreeEdgeProperty,
          typename OutEdgeListS,
          typename DirectedS,
          typename VertexProperty,
          typename AdjEdgeProperty,
          typename AdjEdgeListS,
          typename TreeStorageTag>
struct tree_storage< TreeVertexProperty, TreeEdgeProperty,
                     ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> > {
  typedef ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> alt_tag_type;
  typedef ReaK::graph::alt_tree_view< typename alt_tag_type::template alt<TreeVertexProperty,TreeEdgeProperty>::type > type;
};


template <typename OutEdgeListS,
          typename DirectedS,
          typename VertexProperty,
          typename AdjEdgeProperty,
          typename AdjEdgeListS,
          typename TreeStorageTag>
struct tree_storage_traits< ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS, TreeStorageTag> > :
  tree_storage_traits< TreeStorageTag > { };

};





namespace ReaK {

namespace graph {



template <typename AdjListOnTreeType>
class alt_tree_view {
  private:
    shared_ptr< AdjListOnTreeType > m_alt;

  public:
    typedef alt_tree_view< AdjListOnTreeType > self;

    typedef typename AdjListOnTreeType::tree_type tree_type;

    // Graph traits:
    typedef typename boost::graph_traits< tree_type >::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits< tree_type >::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits< tree_type >::directed_category directed_category;
    typedef typename boost::graph_traits< tree_type >::edge_parallel_category edge_parallel_category;
    typedef typename boost::graph_traits< tree_type >::traversal_category traversal_category;

    static vertex_descriptor null_vertex() { return boost::graph_traits< tree_type >::null_vertex(); };

    // IncidenceGraph traits:
    typedef typename boost::graph_traits< tree_type >::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits< tree_type >::degree_size_type degree_size_type;

    // BidirectionalGraph traits:
    typedef typename boost::graph_traits< tree_type >::in_edge_iterator in_edge_iterator;

    // VertexListGraph traits:
    typedef typename boost::graph_traits< tree_type >::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits< tree_type >::vertices_size_type vertices_size_type;

    // EdgeListGraph traits:
    typedef typename boost::graph_traits< tree_type >::edge_iterator edge_iterator;
    typedef typename boost::graph_traits< tree_type >::edges_size_type edges_size_type;

    // AdjacencyGraph traits:
    typedef typename boost::graph_traits< tree_type >::adjacency_iterator adjacency_iterator;
    typedef typename boost::tree_traits< tree_type >::child_vertex_iterator child_vertex_iterator;

    // PropertyGraph traits:
    typedef typename tree_type::edge_property_type edge_property_type;
    typedef typename tree_type::vertex_property_type vertex_property_type;

    typedef typename AdjListOnTreeType::tree_vertex_bundled vertex_bundled;
    typedef typename AdjListOnTreeType::tree_edge_bundled edge_bundled;

    friend class alt_graph_view<AdjListOnTreeType>;

    alt_tree_view() : m_alt(new AdjListOnTreeType()) { };

    explicit alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g);


    const tree_type& get_tree() const { return m_alt->m_tree; };
    tree_type& get_tree() { return m_alt->m_tree; };


    // Bundled Property-map functions (used by the boost::property_map< self, T Bundle::* > classes).

    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return m_alt->m_tree[v_i].tree_data;
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return m_alt->m_tree[v_i].tree_data;
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_alt->m_tree[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_alt->m_tree[e_i];
    };


    friend
    vertex_property_type& get_raw_vertex_property(self& g, const vertex_descriptor& v_i) {
      return g.get_tree()[v_i];
    };
    friend
    const vertex_property_type& get_raw_vertex_property(const self& g, const vertex_descriptor& v_i) {
      return g.get_tree()[v_i];
    };
    friend
    edge_property_type& get_raw_edge_property(self& g, const edge_descriptor& e_i) {
      return g.get_tree()[e_i];
    };
    friend
    const edge_property_type& get_raw_edge_property(const self& g, const edge_descriptor& e_i) {
      return g.get_tree()[e_i];
    };


    friend const vertex_bundled& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };

    friend void put( self& g, const vertex_descriptor& v_i, const vertex_bundled& value) {
      g[v_i] = value;
    };

    friend const edge_bundled& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };

    friend void put( self& g, const edge_descriptor& e_i, const edge_bundled& value) {
      g[e_i] = value;
    };

    template <typename T, typename Bundle>
    friend
    typename boost::property_map< self, T Bundle::* >::type
    get( T Bundle::* p, self& g) {
      return typename boost::property_map< self, T Bundle::* >::type(&g, p);
    };

    template <typename T, typename Bundle>
    friend
    typename boost::property_map< self, T Bundle::* >::const_type
    get( T Bundle::* p, const self& g) {
      return typename boost::property_map< self, T Bundle::* >::const_type(&g, p);
    };




    // MutablePropertyTreeConcept

    vertex_descriptor create_root_impl(const vertex_property_type& vp) {
      vertex_descriptor v_root = create_root(vp, m_alt->m_tree);
      m_alt->m_adj_list[ m_alt->m_tree[v_root].partner_node ].tree_vertex = v_root;
      return v_root;
    };

    std::pair< vertex_descriptor, edge_descriptor>
      add_child_vertex_impl(vertex_descriptor v,
                       const vertex_property_type& vp,
                       const edge_property_type& ep = edge_property_type()) {
      std::pair< vertex_descriptor, edge_descriptor> result = add_child_vertex(v, vp, ep, m_alt->m_tree);
      m_alt->m_adj_list[ m_alt->m_tree[result.first].partner_node ].tree_vertex = result.first;
      return result;
    };

    template <typename OutputIter>
    OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      return remove_branch(v, it_out, m_alt->m_tree);
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor create_root_impl(vertex_property_type&& vp) {
      vertex_descriptor v_root = create_root(std::move(vp), m_alt->m_tree);
      m_alt->m_adj_list[ m_alt->m_tree[v_root].partner_node ].tree_vertex = v_root;
      return v_root;
    };

    std::pair< vertex_descriptor, edge_descriptor>
      add_child_vertex_impl(vertex_descriptor v, vertex_property_type&& vp) {
      std::pair< vertex_descriptor, edge_descriptor> result = add_child_vertex(v, std::move(vp), m_alt->m_tree);
      m_alt->m_adj_list[ m_alt->m_tree[result.first].partner_node ].tree_vertex = result.first;
      return result;
    };

    std::pair< vertex_descriptor, edge_descriptor>
      add_child_vertex_impl(vertex_descriptor v, vertex_property_type&& vp, edge_property_type&& ep) {
      std::pair< vertex_descriptor, edge_descriptor> result = add_child_vertex(v, std::move(vp), std::move(ep), m_alt->m_tree);
      m_alt->m_adj_list[ m_alt->m_tree[result.first].partner_node ].tree_vertex = result.first;
      return result;
    };
#endif


};

};  // end namespace graph

};  // end namespace ReaK


namespace boost {



template <typename AdjListOnTreeType>
struct raw_vertex_property_map< ReaK::graph::alt_tree_view<AdjListOnTreeType> > {
  typedef boost::raw_vertex_propgraph_map<
    typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_property_type,
    ReaK::graph::alt_tree_view<AdjListOnTreeType> > type;
  typedef boost::raw_vertex_propgraph_map<
    const typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_property_type,
    const ReaK::graph::alt_tree_view<AdjListOnTreeType> > const_type;
};

template <typename AdjListOnTreeType>
struct raw_edge_property_map< ReaK::graph::alt_tree_view<AdjListOnTreeType> > {
  typedef boost::raw_edge_propgraph_map<
    typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::edge_property_type,
    ReaK::graph::alt_tree_view<AdjListOnTreeType> > type;
  typedef boost::raw_edge_propgraph_map<
    const typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::edge_property_type,
    const ReaK::graph::alt_tree_view<AdjListOnTreeType> > const_type;
};

template <typename AdjListOnTreeType>
struct raw_vertex_to_bundle_map< ReaK::graph::alt_tree_view<AdjListOnTreeType> > {
  typedef data_member_property_map<
    typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_bundled,
    typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_property_type > type;
  typedef data_member_property_map<
    const typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_bundled,
    const typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_property_type > const_type;
};



template <typename AdjListOnTreeType, typename T, typename Bundle>
struct property_map< ReaK::graph::alt_tree_view<AdjListOnTreeType>, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, ReaK::graph::alt_tree_view<AdjListOnTreeType>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const ReaK::graph::alt_tree_view<AdjListOnTreeType>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};



};  // end namespace boost



namespace ReaK {

namespace graph {


template <typename AdjListOnTreeType>
typename boost::raw_vertex_property_map< alt_tree_view<AdjListOnTreeType> >::type
  get_raw_vertex_property_map(alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_vertex_property_map< alt_tree_view<AdjListOnTreeType> >::type result_type;
  return result_type(&g);
};

template <typename AdjListOnTreeType>
typename boost::raw_vertex_property_map< alt_tree_view<AdjListOnTreeType> >::const_type
  get_raw_vertex_property_map(const alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_vertex_property_map< alt_tree_view<AdjListOnTreeType> >::const_type result_type;
  return result_type(&g);
};

template <typename AdjListOnTreeType>
typename boost::raw_edge_property_map< alt_tree_view<AdjListOnTreeType> >::type
  get_raw_edge_property_map(alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_edge_property_map< alt_tree_view<AdjListOnTreeType> >::type result_type;
  return result_type(&g);
};

template <typename AdjListOnTreeType>
typename boost::raw_edge_property_map< alt_tree_view<AdjListOnTreeType> >::const_type
  get_raw_edge_property_map(const alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_edge_property_map< alt_tree_view<AdjListOnTreeType> >::const_type result_type;
  return result_type(&g);
};


template <typename AdjListOnTreeType>
typename boost::raw_vertex_to_bundle_map< alt_tree_view<AdjListOnTreeType> >::type
  get_raw_vertex_to_bundle_map(alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_vertex_to_bundle_map< alt_tree_view<AdjListOnTreeType> >::type result_type;
  typedef typename alt_tree_view<AdjListOnTreeType>::vertex_property_type VProp;
  return result_type(&VProp::tree_data);
};

template <typename AdjListOnTreeType>
typename boost::raw_vertex_to_bundle_map< alt_tree_view<AdjListOnTreeType> >::const_type
  get_raw_vertex_to_bundle_map(const alt_tree_view<AdjListOnTreeType>& g) {
  typedef typename boost::raw_vertex_to_bundle_map< alt_tree_view<AdjListOnTreeType> >::const_type result_type;
  typedef const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type VProp;
  return result_type(&VProp::tree_data);
};












template <typename AdjListOnTreeType>
class alt_graph_view {
  public:
    typedef alt_graph_view< AdjListOnTreeType > self;

    typedef typename AdjListOnTreeType::adj_list_type graph_type;
    typedef typename AdjListOnTreeType::tree_type tree_type;

    // Graph traits:
    typedef typename boost::graph_traits< graph_type >::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits< graph_type >::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits< graph_type >::directed_category directed_category;
    typedef typename boost::graph_traits< graph_type >::edge_parallel_category edge_parallel_category;
    typedef typename boost::graph_traits< graph_type >::traversal_category traversal_category;

    static vertex_descriptor null_vertex() { return boost::graph_traits< graph_type >::null_vertex(); };

    // IncidenceGraph traits:
    typedef typename boost::graph_traits< graph_type >::out_edge_iterator out_edge_iterator;
    typedef typename boost::graph_traits< graph_type >::degree_size_type degree_size_type;

    // BidirectionalGraph traits:
    typedef typename boost::graph_traits< graph_type >::in_edge_iterator in_edge_iterator;

    // VertexListGraph traits:
    typedef typename boost::graph_traits< graph_type >::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits< graph_type >::vertices_size_type vertices_size_type;

    // EdgeListGraph traits:
    typedef typename boost::graph_traits< graph_type >::edge_iterator edge_iterator;
    typedef typename boost::graph_traits< graph_type >::edges_size_type edges_size_type;

    // AdjacencyGraph traits:
    typedef typename boost::graph_traits< graph_type >::adjacency_iterator adjacency_iterator;

    // PropertyGraph traits:
    typedef typename AdjListOnTreeType::adj_edge_bundled edge_property_type;
    typedef typename AdjListOnTreeType::adj_vertex_bundled vertex_property_type;

    typedef typename AdjListOnTreeType::adj_vertex_bundled vertex_bundled;
    typedef typename AdjListOnTreeType::adj_edge_bundled edge_bundled;

    typedef void graph_bundled;

  private:

    typedef typename boost::graph_traits< tree_type >::vertex_descriptor tree_vertex_desc;
    typedef typename tree_type::vertex_property_type tree_vertex_prop;

    struct mutation_visitor_base {

      mutation_visitor_base() { };
      virtual ~mutation_visitor_base() { };

      virtual void remove_vertex(tree_vertex_desc, alt_tree_view< AdjListOnTreeType >&) const = 0;
      virtual void add_vertex(const tree_vertex_prop&, alt_tree_view< AdjListOnTreeType >&) const = 0;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      virtual void add_vertex(tree_vertex_prop&&, alt_tree_view< AdjListOnTreeType >&) const = 0;
#endif
    };

    template <typename GraphMutationVisitor>
    struct mutation_visitor : mutation_visitor_base {
      GraphMutationVisitor m_vis;

      BOOST_CONCEPT_ASSERT((TreeOrganizerVisitorConcept<GraphMutationVisitor, alt_tree_view< AdjListOnTreeType > >));

      mutation_visitor(GraphMutationVisitor aVis) : m_vis(aVis) { };
      virtual ~mutation_visitor() { };

      virtual void remove_vertex(tree_vertex_desc tv, alt_tree_view< AdjListOnTreeType >& t) const {
        m_vis.remove_vertex(tv, t);
      };
      virtual void add_vertex(const tree_vertex_prop& tvp, alt_tree_view< AdjListOnTreeType >& t) const {
        m_vis.add_vertex(tvp, t);
      };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      virtual void add_vertex(tree_vertex_prop&& tvp, alt_tree_view< AdjListOnTreeType >& t) const {
        m_vis.add_vertex(std::move(tvp), t);
      };
#endif
    };

    shared_ptr< AdjListOnTreeType > m_alt;
    shared_ptr< mutation_visitor_base > m_vis;

  public:

    friend class alt_tree_view<AdjListOnTreeType>;

    template <typename GraphMutationVisitor>
    explicit alt_graph_view(GraphMutationVisitor aVis) :
                            m_alt(new AdjListOnTreeType()),
                            m_vis(new mutation_visitor<GraphMutationVisitor>(aVis)) { };

    template <typename GraphMutationVisitor>
    alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t, GraphMutationVisitor aVis);




    // Bundled Property-map functions (used by the boost::property_map< self, T Bundle::* > classes).

    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return m_alt->m_tree[ m_alt->m_adj_list[v_i].tree_vertex ].user_data;
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return m_alt->m_tree[ m_alt->m_adj_list[v_i].tree_vertex ].user_data;
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_alt->m_adj_list[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_alt->m_adj_list[e_i];
    };

    friend const vertex_bundled& get( const self& g, const vertex_descriptor& v_i) {
      return g[v_i];
    };

    friend void put( self& g, const vertex_descriptor& v_i, const vertex_bundled& value) {
      g[v_i] = value;
    };

    friend const edge_bundled& get( const self& g, const edge_descriptor& e_i) {
      return g[e_i];
    };

    friend void put( self& g, const edge_descriptor& e_i, const edge_bundled& value) {
      g[e_i] = value;
    };

    template <typename T, typename Bundle>
    friend
    typename boost::property_map< self, T Bundle::* >::type
    get( T Bundle::* p, self& g) {
      return typename boost::property_map< self, T Bundle::* >::type(&g, p);
    };

    template <typename T, typename Bundle>
    friend
    typename boost::property_map< self, T Bundle::* >::const_type
    get( T Bundle::* p, const self& g) {
      return typename boost::property_map< self, T Bundle::* >::const_type(&g, p);
    };




    // IncidenceGraph concept

    std::pair<out_edge_iterator, out_edge_iterator> out_edges_impl(vertex_descriptor v) const {
      return out_edges(v, m_alt->m_adj_list);
    };
    vertex_descriptor source_impl(edge_descriptor e) const {
      return source(e, m_alt->m_adj_list);
    };
    vertex_descriptor target_impl(edge_descriptor e) const {
      return target(e, m_alt->m_adj_list);
    };
    degree_size_type out_degree_impl(vertex_descriptor v) const {
      return out_degree(v, m_alt->m_adj_list);
    };

    // BidirectionalGraph concept

    std::pair<in_edge_iterator, in_edge_iterator> in_edges_impl(vertex_descriptor v) const {
      return in_edges(v, m_alt->m_adj_list);
    };
    degree_size_type in_degree_impl(vertex_descriptor v) const {
      return in_degree(v, m_alt->m_adj_list);
    };
    degree_size_type degree_impl(edge_descriptor e) const {
      return degree(e, m_alt->m_adj_list);
    };

    // AdjacencyGraph concept

    std::pair<adjacency_iterator, adjacency_iterator> adjacent_vertices_impl(vertex_descriptor u) const {
      return adjacent_vertices(u, m_alt->m_adj_list);
    };

    // VertexListGraph concept

    std::pair< vertex_iterator, vertex_iterator > vertices_impl() const {
      return vertices(m_alt->m_adj_list);
    };
    vertices_size_type num_vertices_impl() const {
      return num_vertices(m_alt->m_adj_list);
    };

    // EdgeListGraph concept

    std::pair< edge_iterator, edge_iterator > edges_impl() const {
      return edges(m_alt->m_adj_list);
    };
    edges_size_type num_edges_impl() const {
      return num_edges(m_alt->m_adj_list);
    };

    // AdjacencyMatrix concept

    std::pair<edge_descriptor,bool> edge_impl(vertex_descriptor u, vertex_descriptor v) const {
      return edge(u,v,m_alt->m_adj_list);
    };

    // MutableGraph concept

    void clear_vertex_impl(vertex_descriptor v) {
      clear_vertex(v, m_alt->m_adj_list);
    };

    void remove_vertex_impl(vertex_descriptor v) {
      alt_tree_view< AdjListOnTreeType > tree_view(*this);
      m_vis->remove_vertex(m_alt->m_adj_list[v].tree_vertex, tree_view);
      clear_vertex(v, m_alt->m_adj_list);
      remove_vertex(v, m_alt->m_adj_list);
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v) {
      return add_edge(u, v, m_alt->m_adj_list);
    };

    void remove_edge_impl(vertex_descriptor u, vertex_descriptor v) {
      remove_edge(u, v, m_alt->m_adj_list);
    };

    void remove_edge_impl(edge_descriptor e) {
      remove_edge(e, m_alt->m_adj_list);
    };

    void remove_edge_impl(edge_iterator e_iter) {
      remove_edge(e_iter, m_alt->m_adj_list);
    };

    // MutablePropertyGraph concept

    vertex_descriptor add_vertex_impl(const vertex_property_type& vp) {
      vertex_descriptor v = add_vertex(m_alt->m_adj_list);
      typename tree_type::vertex_property_type tree_vp;
      tree_vp.partner_node = v;
      tree_vp.user_data = vp;
      alt_tree_view< AdjListOnTreeType > tree_view(*this);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      m_vis->add_vertex(std::move(tree_vp), tree_view);
#else
      m_vis->add_vertex(tree_vp, tree_view);
#endif
      return v;
    };

    void remove_vertex_impl(vertex_descriptor v, vertex_property_type& vp) {
      vp = m_alt->m_tree[ m_alt->m_adj_list[v].tree_vertex ].user_data;
      remove_vertex_impl(v);
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v, const edge_property_type& ep) {
      return add_edge(u, v, ep, m_alt->m_adj_list);
    };

    void remove_edge_impl(vertex_descriptor u, vertex_descriptor v, edge_property_type& ep) {
      remove_edge(u, v, ep, m_alt->m_adj_list);
    };

    void remove_edge_impl(edge_descriptor e, edge_property_type& ep) {
      remove_edge(e, ep, m_alt->m_adj_list);
    };

    void remove_edge_impl(edge_iterator e_iter, edge_property_type& ep) {
      remove_edge(e_iter, ep, m_alt->m_adj_list);
    };

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_vertex_impl(vertex_property_type&& vp) {
      vertex_descriptor v = add_vertex(m_alt->m_adj_list);
      typename tree_type::vertex_property_type tree_vp;
      tree_vp.partner_node = v;
      tree_vp.user_data = std::move(vp);
      alt_tree_view< AdjListOnTreeType > tree_view(*this);
      m_vis->add_vertex(std::move(tree_vp), tree_view);
      return v;
    };

    std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u, vertex_descriptor v, edge_property_type&& ep) {
      return add_edge(u, v, std::move(ep), m_alt->m_adj_list);
    };
#endif

    // Other useful / unclassified functions:

    void update_vertex(vertex_descriptor v) {
      alt_tree_view< AdjListOnTreeType > tree_view(*this);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      typename tree_type::vertex_property_type tree_vp = std::move(m_alt->m_tree[ m_alt->m_adj_list[v].tree_vertex ]);
      m_vis->remove_vertex(m_alt->m_adj_list[v].tree_vertex, tree_view);
      m_vis->add_vertex(std::move(tree_vp), tree_view);
#else
      typename tree_type::vertex_property_type tree_vp = m_alt->m_tree[ m_alt->m_adj_list[v].tree_vertex ];
      m_vis->remove_vertex(m_alt->m_adj_list[v].tree_vertex, tree_view);
      m_vis->add_vertex(tree_vp, tree_view);
#endif
    };

    
#if 0
    // this could go into NonCompactGraphConcept
    void repack_graph() {
      if(m_alt->m_available_pnodes.empty())
        return;
      // here, the idea is to essentially empty out the graveyard,
      if( boost::is_same< vertices_size_type, vertex_descriptor >::type::value ) {
        // if vertex descriptors are into a random-access container, then a simple remove loop won't cut it.
        // the idea here is to achieve a behavior similar to the std::remove function, but by taking
        // elements from the back of the vertex list and swapping them with the ones in the graveyard.
        using std::swap;
        vertex_descriptor v_end = reinterpret_cast<vertex_descriptor>(num_vertices(m_alt->m_adj_list));
        while(!m_alt->m_available_pnodes.empty()) {
          vertex_descriptor v = m_alt->m_available_pnodes.top();
          if(v > v_end) {
            ++v_end;
            break;
          };
          swap(m_alt->m_adj_list[v], m_alt->m_adj_list[v_end]);
          out_edge_iterator ei, ei_end;
          for(boost::tie(ei,ei_end) = out_edges(v_end, m_alt->m_adj_list); ei != ei_end; ++ei) {
            std::pair<edge_descriptor, bool> e = add_edge(v, target(*ei, m_alt->m_adj_list), m_alt->m_adj_list);
            swap(m_alt->m_adj_list[e], m_alt->m_adj_list[*ei]);
          };
          if( boost::is_same< directed_category, boost::directed_tag >::type::value ) {
            in_edge_iterator in_ei, in_ei_end;
            for(boost::tie(in_ei, in_ei_end) = in_edges(v_end, m_alt->m_adj_list); in_ei != in_ei_end; ++in_ei) {
              std::pair<edge_descriptor, bool> e = add_edge(source(*in_ei, m_alt->m_adj_list), v, m_alt->m_adj_list);
              swap(m_alt->m_adj_list[e], m_alt->m_adj_list[*in_ei]);
            };
          };
          clear_vertex(v_end, m_alt->m_adj_list);
          m_alt->m_available_pnodes.pop();
        };
        while(!m_alt->m_available_pnodes.empty())
          m_alt->m_available_pnodes.pop();
        // at this point v_end is the first invalid vertex remaining.
        vertex_descriptor v_rm = reinterpret_cast<vertex_descriptor>(num_vertices(m_alt->m_adj_list));
        while(v_rm > v_end)
          remove_vertex(--v_rm, m_alt->m_adj_list);

      } else {
        // if vertex descriptors are not random-access, they should not get invalidated by removal of other vertices.
        while(!m_alt->m_available_pnodes.empty()) {
          clear_vertex(m_alt->m_available_pnodes.top(), m_alt->m_adj_list);
          remove_vertex(m_alt->m_available_pnodes.top(), m_alt->m_adj_list);
          m_alt->m_available_pnodes.pop();
        };
      };
      // re-updating all partner-vertex descriptors in the tree.
      {
      vertex_iterator vi, vi_end;
      for(boost::tie(vi,vi_end) = vertices(m_alt->m_adj_list); vi != vi_end; ++vi)
        m_alt->m_tree[ m_alt->m_adj_list[ *vi ].tree_vertex ].partner_node = *vi;
      };
    };
#endif 

};

template <typename AdjListOnTreeType>
alt_tree_view<AdjListOnTreeType>::alt_tree_view(const alt_graph_view< AdjListOnTreeType >& g) : m_alt(g.m_alt) { };

template <typename AdjListOnTreeType>
template <typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType>::alt_graph_view(const alt_tree_view< AdjListOnTreeType >& t,
                                                  GraphMutationVisitor aVis ) :
                                                  m_alt(t.m_alt),
                                                  m_vis(new typename alt_graph_view<AdjListOnTreeType>::template mutation_visitor<GraphMutationVisitor>(aVis)) { };



template <typename AdjListOnTreeType, typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType> make_graph_view(const alt_tree_view< AdjListOnTreeType >& t, GraphMutationVisitor aVis) {
  return alt_graph_view<AdjListOnTreeType>(t,aVis);
};








/******************************************************************************************
 * *************************************************************************************
 *                             Tree View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::out_edge_iterator,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::out_edge_iterator >
  out_edges(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
            const alt_tree_view<AdjListOnTreeType>& g) {
  return out_edges(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  source(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_tree_view<AdjListOnTreeType>& g) {
  return source(e, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  target(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_tree_view<AdjListOnTreeType>& g) {
  return target(e, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type
  out_degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
             const alt_tree_view<AdjListOnTreeType>& g) {
  return out_degree(v, g.get_tree());
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::in_edge_iterator,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::in_edge_iterator >
  in_edges(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
           const alt_tree_view<AdjListOnTreeType>& g) {
  return in_edges(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type
  in_degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
            const alt_tree_view<AdjListOnTreeType>& g) {
  return in_degree(v, g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::degree_size_type
  degree(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_tree_view<AdjListOnTreeType>& g) {
  return degree(e, g.get_tree());
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::adjacency_iterator,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::adjacency_iterator >
  adjacent_vertices(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor u,
                    const alt_tree_view<AdjListOnTreeType>& g) {
  return adjacent_vertices(u, g.get_tree());
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_iterator >
  vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return vertices(g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertices_size_type
  num_vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_vertices(g.get_tree());
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_iterator,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_iterator >
  edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return edges(g.get_tree());
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edges_size_type
  num_edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_edges(g.get_tree());
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor, bool>
  edge(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor u,
       typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
       const alt_tree_view<AdjListOnTreeType>& g) {
  return edge(u,v,g.get_tree());
};

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  get_root_vertex(const alt_tree_view<AdjListOnTreeType>& g) {
  return get_root_vertex(g.get_tree());
};

template <typename AdjListOnTreeType>
std::pair< typename boost::tree_traits< alt_tree_view<AdjListOnTreeType> >::child_vertex_iterator,
           typename boost::tree_traits< alt_tree_view<AdjListOnTreeType> >::child_vertex_iterator >
  child_vertices(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                 const alt_tree_view<AdjListOnTreeType>& g) {
  return child_vertices(v,g.get_tree());
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  create_root(const alt_tree_view<AdjListOnTreeType>& g) {
  return g.create_root_can_only_be_called_with_a_vertex_property();
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor>
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   const alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex_can_only_be_called_with_a_vertex_property(v);
};


template <typename AdjListOnTreeType>
void remove_branch(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   alt_tree_view<AdjListOnTreeType>& g) {
  g.remove_branch_can_only_be_called_with_a_vertex_property_range_for_output(v);
};


/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  create_root(const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp,
              alt_tree_view<AdjListOnTreeType>& g) {
  return g.create_root_impl(vp);
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor>
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp,
                   alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex_impl(v, vp);
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor>
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   const typename alt_tree_view<AdjListOnTreeType>::vertex_property_type& vp,
                   const typename alt_tree_view<AdjListOnTreeType>::edge_property_type& ep,
                   alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex_impl(v, vp, ep);
};


template <typename AdjListOnTreeType, typename OutputIter>
OutputIter remove_branch(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                         OutputIter it_out,
                         alt_tree_view<AdjListOnTreeType>& g) {
  return g.remove_branch_impl(v,it_out);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor
  create_root(typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp,
              alt_tree_view<AdjListOnTreeType>& g) {
  return g.create_root_impl(std::move(vp));
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor>
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp,
                   alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex_impl(v, std::move(vp));
};

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor,
           typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::edge_descriptor>
  add_child_vertex(typename boost::graph_traits< alt_tree_view<AdjListOnTreeType> >::vertex_descriptor v,
                   typename alt_tree_view<AdjListOnTreeType>::vertex_property_type&& vp,
                   typename alt_tree_view<AdjListOnTreeType>::edge_property_type&& ep,
                   alt_tree_view<AdjListOnTreeType>& g) {
  return g.add_child_vertex_impl(v, std::move(vp), std::move(ep));
};
#endif






/******************************************************************************************
 * *************************************************************************************
 *                             Graph View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::out_edge_iterator,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::out_edge_iterator >
  out_edges(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
            const alt_graph_view<AdjListOnTreeType>& g) {
  return g.out_edges_impl(v);
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor
  source(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_graph_view<AdjListOnTreeType>& g) {
  return g.source_impl(e);
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor
  target(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_graph_view<AdjListOnTreeType>& g) {
  return g.target_impl(e);
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::degree_size_type
  out_degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
             const alt_graph_view<AdjListOnTreeType>& g) {
  return g.out_degree_impl(v);
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::in_edge_iterator,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::in_edge_iterator >
  in_edges(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
           const alt_graph_view<AdjListOnTreeType>& g) {
  return g.in_edges_impl(v);
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::degree_size_type
  in_degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
            const alt_graph_view<AdjListOnTreeType>& g) {
  return g.in_degree_impl(v);
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::degree_size_type
  degree(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor e,
         const alt_graph_view<AdjListOnTreeType>& g) {
  return g.degree_impl(e);
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::adjacency_iterator,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::adjacency_iterator >
  adjacent_vertices(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
                    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.adjacent_vertices_impl(u);
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_iterator,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_iterator >
  vertices(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.vertices_impl();
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertices_size_type
  num_vertices(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.num_vertices_impl();
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair< typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_iterator,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_iterator >
  edges(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.edges_impl();
};

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edges_size_type
  num_edges(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.num_edges_impl();
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor, bool>
  edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
       typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
       const alt_graph_view<AdjListOnTreeType>& g) {
  return g.edge_impl(u,v);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor
  add_vertex(alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_vertex_function_can_only_be_called_with_a_vertex_property();
//  return g.add_vertex();
};

template <typename AdjListOnTreeType>
void clear_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
                  alt_graph_view<AdjListOnTreeType>& g) {
  g.clear_vertex_impl(v);
};

template <typename AdjListOnTreeType>
void remove_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
                   alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_vertex_impl(v);
};

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor, bool>
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
           alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_edge_impl(u,v);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
                 typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(u,v);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor e,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(e);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_iterator e_iter,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(e_iter);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor
  add_vertex(const typename alt_graph_view<AdjListOnTreeType>::vertex_property_type& vp,
             alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_vertex_impl(vp);
};

template <typename AdjListOnTreeType>
void remove_vertex(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
                   typename alt_graph_view<AdjListOnTreeType>::vertex_property_type& vp,
                   alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_vertex_impl(v,vp);
};

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor, bool>
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
           const typename alt_graph_view<AdjListOnTreeType>::edge_property_type& ep,
           alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_edge_impl(u,v,ep);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
                 typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
                 typename alt_graph_view<AdjListOnTreeType>::edge_property_type& ep,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(u,v,ep);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor e,
                 typename alt_graph_view<AdjListOnTreeType>::edge_property_type& ep,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(e,ep);
};

template <typename AdjListOnTreeType>
void remove_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_iterator e_iter,
                 typename alt_graph_view<AdjListOnTreeType>::edge_property_type& ep,
                 alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(e_iter,ep);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename AdjListOnTreeType>
typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor
  add_vertex(typename alt_graph_view<AdjListOnTreeType>::vertex_property_type&& vp,
             alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_vertex_impl(std::move(vp));
};

template <typename AdjListOnTreeType>
std::pair<typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::edge_descriptor, bool>
  add_edge(typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor u,
           typename boost::graph_traits< alt_graph_view<AdjListOnTreeType> >::vertex_descriptor v,
           typename alt_graph_view<AdjListOnTreeType>::edge_property_type&& ep,
           alt_graph_view<AdjListOnTreeType>& g) {
  return g.add_edge_impl(u,v,std::move(ep));
};
#endif


};

};



namespace boost {

template <typename TreeVertexProperty, typename AdjVertexProperty>
struct property_map< ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>,
                     vertex_key_t > {
  typedef data_member_property_map<std::size_t, ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> > type;
  typedef data_member_property_map<const std::size_t, const ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> > const_type;
};

template <typename TreeVertexProperty, typename AdjVertexProperty>
struct property_map< ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>,
                     vertex_second_bundle_t > {
  typedef data_member_property_map<AdjVertexProperty, ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> > type;
  typedef data_member_property_map<const AdjVertexProperty, const ReaK::graph::detail::tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty> > const_type;
};




template <typename AdjListOnTreeType, typename T, typename Bundle>
struct property_map< ReaK::graph::alt_graph_view<AdjListOnTreeType>, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename ReaK::graph::alt_graph_view<AdjListOnTreeType>::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, ReaK::graph::alt_graph_view<AdjListOnTreeType>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const ReaK::graph::alt_graph_view<AdjListOnTreeType>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


};


#endif


















