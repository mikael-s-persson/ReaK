
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

#include <iostream>

#include "d_ary_bf_tree.hpp"
#include "d_ary_cob_tree.hpp"
#include <boost/graph/adjacency_list.hpp>
#include "bgl_tree_adaptor.hpp"
#include "linked_tree.hpp"
#include "pooled_adjacency_list.hpp"


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE bgl_trees
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


typedef boost::mpl::list< 
  ReaK::graph::d_ary_bf_tree<int, 4, int>, 
  ReaK::graph::d_ary_cob_tree<int, 4, int>, 
  ReaK::graph::linked_tree<boost::vecS, boost::vecS, int, int>,
  ReaK::graph::linked_tree<boost::listS, boost::vecS, int, int>, 
  ReaK::graph::linked_tree<boost::vecS, boost::listS, int, int>, 
  ReaK::graph::linked_tree<boost::listS, boost::listS, int, int>,
  ReaK::graph::tree_storage<int, int>::type,
  boost::pooled_adjacency_list<boost::bidirectionalS, int, int > > intint_treetest_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( intint_tree_test, TreeType, intint_treetest_types )
{
  typedef typename boost::graph_traits<TreeType>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<TreeType>::vertex_iterator VertexIter;
  typedef typename boost::graph_traits<TreeType>::edge_descriptor Edge;
  typedef typename boost::graph_traits<TreeType>::out_edge_iterator OutEdgeIter;
  typedef typename boost::graph_traits<TreeType>::in_edge_iterator InEdgeIter;
  typedef typename ReaK::graph::tree_traits<TreeType>::child_vertex_iterator ChildVertIter;
  
  TreeType g;
  
  Vertex v_root;
  BOOST_CHECK_NO_THROW( v_root = create_root(g) );
  BOOST_CHECK_NO_THROW( remove_branch(v_root,g) );
  int vp_r = 1;
  BOOST_CHECK_NO_THROW( v_root = create_root(vp_r, g) );
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  
  std::vector<int> props;
  BOOST_CHECK_NO_THROW( remove_branch(v_root, back_inserter(props), g) );
  BOOST_CHECK_EQUAL( props.size(), 1 );
  BOOST_CHECK_EQUAL( props[0], 1 );
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  BOOST_CHECK_NO_THROW( v_root = create_root(std::move(vp_r), g) );
#else
  BOOST_CHECK_NO_THROW( v_root = create_root(vp_r, g) );
#endif
  props.clear();
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  
  int vp_rc[] = {2,3,4,5};
  int ep_rc[] = {1002,1003,1004,1005};
  Vertex v_rc[4];
  Edge e_rc[4];
  for(int i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( boost::tie(v_rc[i],e_rc[i]) = add_child_vertex(v_root,vp_rc[i],ep_rc[i],g) );
  };
  BOOST_CHECK_EQUAL( num_vertices(g), 5 );
  
  int vp_rc1c[] = {6,7,8,9};
  int ep_rc1c[] = {2006,2007,2008,2009};
  Vertex v_rc1c[4];
  Edge e_rc1c[4];
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( boost::tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g) );
  };
  BOOST_CHECK_EQUAL( num_vertices(g), 9 );
  
  
  BOOST_CHECK_NO_THROW( v_root = get_root_vertex(g) );
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1002);
    BOOST_CHECK_EQUAL( e_list[1], 1003);
    BOOST_CHECK_EQUAL( e_list[2], 1004);
    BOOST_CHECK_EQUAL( e_list[3], 1005);
    
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      if(is_vertex_valid(*cvi,g))
        vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 2);
    BOOST_CHECK_EQUAL( vp_list[1], 3);
    BOOST_CHECK_EQUAL( vp_list[2], 4);
    BOOST_CHECK_EQUAL( vp_list[3], 5);
    
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if((is_vertex_valid(*cvi,g)) && (g[*cvi] == 2)) {
        
        BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          if(is_edge_valid(*ei,g)) {
            BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
            e_list2.push_back(g[*ei]);
          };
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 2006);
        BOOST_CHECK_EQUAL( e_list2[1], 2007);
        BOOST_CHECK_EQUAL( e_list2[2], 2008);
        BOOST_CHECK_EQUAL( e_list2[3], 2009);
        
      };
    };
  };
  
  int vp_rc2c[] = {10,11,12,13};
  int ep_rc2c[] = {3010,3011,3012,3013};
  Vertex v_rc2c[4];
  Edge e_rc2c[4];
  for(std::size_t i = 0; i < 4; ++i) {
#ifdef RK_ENABLE_CXX0X_FEATURES
    BOOST_CHECK_NO_THROW( boost::tie(v_rc2c[i],e_rc2c[i]) = add_child_vertex(v_rc[1],std::move(vp_rc2c[i]),std::move(ep_rc2c[i]),g) );
#else
    BOOST_CHECK_NO_THROW( boost::tie(v_rc2c[i],e_rc2c[i]) = add_child_vertex(v_rc[1],vp_rc2c[i],ep_rc2c[i],g) );
#endif
  };
  
  BOOST_CHECK_EQUAL( num_vertices(g), 13 );
  
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_rc[1],g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 3010);
    BOOST_CHECK_EQUAL( e_list[1], 3011);
    BOOST_CHECK_EQUAL( e_list[2], 3012);
    BOOST_CHECK_EQUAL( e_list[3], 3013);
    
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_rc[1],g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      if(is_vertex_valid(*cvi,g))
        vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 10);
    BOOST_CHECK_EQUAL( vp_list[1], 11);
    BOOST_CHECK_EQUAL( vp_list[2], 12);
    BOOST_CHECK_EQUAL( vp_list[3], 13);
  };
  
  remove_branch(v_rc[0],g);
  BOOST_CHECK_EQUAL( num_vertices(g), 8 );
  
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1003);
    BOOST_CHECK_EQUAL( e_list[1], 1004);
    BOOST_CHECK_EQUAL( e_list[2], 1005);
    
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      if(is_vertex_valid(*cvi,g))
        vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 3);
    BOOST_CHECK_EQUAL( vp_list[1], 4);
    BOOST_CHECK_EQUAL( vp_list[2], 5);
    
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if((is_vertex_valid(*cvi,g)) && (g[*cvi] == 3)) {
        
        BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          if(is_edge_valid(*ei,g)) {
            BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
            e_list2.push_back(g[*ei]);
          };
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 3010);
        BOOST_CHECK_EQUAL( e_list2[1], 3011);
        BOOST_CHECK_EQUAL( e_list2[2], 3012);
        BOOST_CHECK_EQUAL( e_list2[3], 3013);
        
      };
    };
  };
  
  remove_branch(v_rc[1], back_inserter(props), g);
  BOOST_CHECK_EQUAL( props.size(), 5 );
  BOOST_CHECK_EQUAL( props[0], 3);  // the first vertex should be the root of the branch.
  std::sort(props.begin(), props.end());
  BOOST_CHECK_EQUAL( props[1], 10);
  BOOST_CHECK_EQUAL( props[2], 11);
  BOOST_CHECK_EQUAL( props[3], 12);
  BOOST_CHECK_EQUAL( props[4], 13);
  
  BOOST_CHECK_EQUAL( num_vertices(g), 3 );
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1004);
    BOOST_CHECK_EQUAL( e_list[1], 1005);
    
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      if(is_vertex_valid(*cvi,g))
        vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 4);
    BOOST_CHECK_EQUAL( vp_list[1], 5);
    
  };
  
  
  BOOST_CHECK_NO_THROW( boost::tie(v_rc[0],e_rc[0]) = add_child_vertex(v_root,vp_rc[0],ep_rc[0],g) );
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( boost::tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g) );
  };
  
  BOOST_CHECK_EQUAL( num_vertices(g), 8 );
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1002);
    BOOST_CHECK_EQUAL( e_list[1], 1004);
    BOOST_CHECK_EQUAL( e_list[2], 1005);
    
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      if(is_vertex_valid(*cvi,g))
        vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 2);
    BOOST_CHECK_EQUAL( vp_list[1], 4);
    BOOST_CHECK_EQUAL( vp_list[2], 5);
    
    BOOST_CHECK_NO_THROW( boost::tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if((is_vertex_valid(*cvi,g)) && (g[*cvi] == 2)) {
        
        BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          if(is_edge_valid(*ei,g)) {
            BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
            e_list2.push_back(g[*ei]);
          };
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 2006);
        BOOST_CHECK_EQUAL( e_list2[1], 2007);
        BOOST_CHECK_EQUAL( e_list2[2], 2008);
        BOOST_CHECK_EQUAL( e_list2[3], 2009);
        
      };
    };
  };
  
  
  std::vector<int> all_vertices;
  VertexIter vi, vi_end;
  for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
    if(is_vertex_valid(*vi,g))
      all_vertices.push_back(g[*vi]);
  std::sort(all_vertices.begin(), all_vertices.end());
  BOOST_CHECK_EQUAL( all_vertices[0], 1 );
  BOOST_CHECK_EQUAL( all_vertices[1], 2 );
  BOOST_CHECK_EQUAL( all_vertices[2], 4 );
  BOOST_CHECK_EQUAL( all_vertices[3], 5 );
  BOOST_CHECK_EQUAL( all_vertices[4], 6 );
  BOOST_CHECK_EQUAL( all_vertices[5], 7 );
  BOOST_CHECK_EQUAL( all_vertices[6], 8 );
  BOOST_CHECK_EQUAL( all_vertices[7], 9 );
  
  BOOST_CHECK_EQUAL( vp_rc1c[1], 7 );
  InEdgeIter iei, iei_end;
  unsigned int j = 0;
  for(boost::tie(iei,iei_end) = in_edges(v_rc1c[1],g); (iei != iei_end) && (j < 4); ++iei, ++j) {
    if(is_edge_valid(*iei,g)) {
      BOOST_CHECK_EQUAL( g[*iei], 2007 );
      BOOST_CHECK_EQUAL( g[source(*iei,g)], 2 );
      BOOST_CHECK_EQUAL( g[target(*iei,g)], 7 );
    };
  };
  
  BOOST_CHECK_EQUAL( vp_rc[0], 2 );
  BOOST_CHECK_EQUAL( vp_rc1c[2], 8 );
  Edge e_rc_rc1c;
  BOOST_CHECK_NO_THROW( e_rc_rc1c = edge(v_rc[0],v_rc1c[2],g).first );
  BOOST_CHECK_EQUAL( g[e_rc_rc1c], 2008 );
  BOOST_CHECK_EQUAL( g[source(e_rc_rc1c,g)], 2 );
  BOOST_CHECK_EQUAL( g[target(e_rc_rc1c,g)], 8 );
  
};




