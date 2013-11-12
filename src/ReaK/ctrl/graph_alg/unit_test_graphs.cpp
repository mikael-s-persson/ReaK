
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

#define BOOST_TEST_MODULE bgl_graphs
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


typedef boost::mpl::list< 
  boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS, int, int>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::bidirectionalS, int, int>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS, int, int>,
  boost::adjacency_list< boost::listS, boost::listS, boost::bidirectionalS, int, int>,
  boost::pooled_adjacency_list<boost::bidirectionalS, int, int > > intint_graphtest_types;

BOOST_AUTO_TEST_CASE_TEMPLATE( intint_bgl_mutable_graph_test, Graph, intint_graphtest_types )
{
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIter;
  typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
  typedef typename boost::graph_traits<Graph>::edge_iterator EdgeIter;
  typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
  typedef typename boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;
  
  Graph g;
  
  Vertex v_root = Vertex();
  BOOST_CHECK_NO_THROW( v_root = add_vertex(g) );
  BOOST_CHECK_EQUAL( num_vertices(g), 1);
  BOOST_CHECK_NO_THROW( remove_vertex(v_root,g) );
  BOOST_CHECK_EQUAL( num_vertices(g), 0);
  
  BOOST_CHECK_NO_THROW( v_root = add_vertex(1, g) );
  g[v_root] = 1;
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  
  int vp_rc[] = {2,3,4,5};
  int ep_rc[] = {1002,1003,1004,1005};
  Vertex v_rc[4];
  Edge e_rc[4];
  for(int i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( v_rc[i] = add_vertex(g) );
    g[ v_rc[i] ] = vp_rc[i];
    BOOST_CHECK_EQUAL( g[ v_rc[i] ], vp_rc[i] );
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW( boost::tie(e_rc[i],edge_added_success) = add_edge(v_root, v_rc[i], g) );
    BOOST_CHECK( edge_added_success );
    g[ e_rc[i] ] = ep_rc[i];
    BOOST_CHECK_EQUAL( g[ e_rc[i] ], ep_rc[i] );
  };
  BOOST_CHECK_EQUAL( num_vertices(g), 5 );
  
  int vp_rc1c[] = {6,7,8,9};
  int ep_rc1c[] = {2006,2007,2008,2009};
  Vertex v_rc1c[4];
  Edge e_rc1c[4];
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( v_rc1c[i] = add_vertex(vp_rc1c[i], g) );
    BOOST_CHECK_EQUAL( g[ v_rc1c[i] ], vp_rc1c[i] );
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW( boost::tie(e_rc1c[i],edge_added_success) = add_edge(v_rc[0], v_rc1c[i], ep_rc1c[i], g) );
    BOOST_CHECK( edge_added_success );
    BOOST_CHECK_EQUAL( g[ e_rc1c[i] ], ep_rc1c[i] );
  };
  BOOST_CHECK_EQUAL( num_vertices(g), 9 );
  
  
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
    
    
    InEdgeIter iei, iei_end;
    BOOST_CHECK_NO_THROW( boost::tie(iei, iei_end) = in_edges(v_rc[0], g) );
    BOOST_CHECK( iei != iei_end );
    BOOST_CHECK_EQUAL( g[*iei], 1002);
    ++iei;
    BOOST_CHECK( iei == iei_end );
    
    
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_rc[0],g) );
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
  
  int vp_rc2c[] = {10,11,12,13};
  int ep_rc2c[] = {3010,3011,3012,3013};
  Vertex v_rc2c[4];
  Edge e_rc2c[4];
  for(std::size_t i = 0; i < 4; ++i) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    BOOST_CHECK_NO_THROW( v_rc2c[i] = add_vertex(std::move(vp_rc2c[i]), g) );
#else
    BOOST_CHECK_NO_THROW( v_rc2c[i] = add_vertex(vp_rc2c[i], g) );
#endif
    bool edge_added_success = false;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    BOOST_CHECK_NO_THROW( boost::tie(e_rc2c[i],edge_added_success) = add_edge(v_rc[1], v_rc2c[i], ep_rc2c[i], g) );
#else
    BOOST_CHECK_NO_THROW( boost::tie(e_rc2c[i],edge_added_success) = add_edge(v_rc[1], v_rc2c[i], ep_rc2c[i], g) );
#endif
    BOOST_CHECK( edge_added_success );
  };
  BOOST_CHECK_EQUAL( num_vertices(g), 13 );
  
  {
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei,ei_end) = out_edges(v_rc[1],g) );
    std::vector<int> e_list;
    std::vector<int> vp_list;
    for(; ei != ei_end; ++ei) {
      if(is_edge_valid(*ei,g)) {
        BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
        e_list.push_back(g[*ei]);
        vp_list.push_back(g[target(*ei,g)]);
      };
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 3010);
    BOOST_CHECK_EQUAL( e_list[1], 3011);
    BOOST_CHECK_EQUAL( e_list[2], 3012);
    BOOST_CHECK_EQUAL( e_list[3], 3013);
    
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 10);
    BOOST_CHECK_EQUAL( vp_list[1], 11);
    BOOST_CHECK_EQUAL( vp_list[2], 12);
    BOOST_CHECK_EQUAL( vp_list[3], 13);
  };
  
  BOOST_CHECK_NO_THROW( clear_vertex(v_rc[0],g) );
  
  BOOST_CHECK_EQUAL( out_degree(v_rc[0], g), 0 );
  BOOST_CHECK_EQUAL( in_degree(v_rc[0], g), 0 );
  BOOST_CHECK_EQUAL( out_degree(v_root, g), 3 );
  BOOST_CHECK_EQUAL( in_degree(v_rc1c[0], g), 0 );
  BOOST_CHECK_EQUAL( in_degree(v_rc1c[1], g), 0 );
  BOOST_CHECK_EQUAL( in_degree(v_rc1c[2], g), 0 );
  BOOST_CHECK_EQUAL( in_degree(v_rc1c[3], g), 0 );
  
  BOOST_CHECK_EQUAL( num_vertices(g), 13 );
  {
    VertexIter vi, vi_end;
    BOOST_CHECK_NO_THROW( boost::tie(vi, vi_end) = vertices(g) );
    std::vector<int> vp_list;
    for(; vi != vi_end; ++vi)
      if( is_vertex_valid(*vi, g) )
        vp_list.push_back( g[*vi] );
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 1 );
    BOOST_CHECK_EQUAL( vp_list[1], 2 );
    BOOST_CHECK_EQUAL( vp_list[2], 3 );
    BOOST_CHECK_EQUAL( vp_list[3], 4 );
    BOOST_CHECK_EQUAL( vp_list[4], 5 );
    BOOST_CHECK_EQUAL( vp_list[5], 6 );
    BOOST_CHECK_EQUAL( vp_list[6], 7 );
    BOOST_CHECK_EQUAL( vp_list[7], 8 );
    BOOST_CHECK_EQUAL( vp_list[8], 9 );
    BOOST_CHECK_EQUAL( vp_list[9], 10 );
    BOOST_CHECK_EQUAL( vp_list[10], 11 );
    BOOST_CHECK_EQUAL( vp_list[11], 12 );
    BOOST_CHECK_EQUAL( vp_list[12], 13 );
  };
  
  BOOST_CHECK_EQUAL( num_edges(g), 7 );
  {
    EdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei, ei_end) = edges(g) );
    std::vector<int> ep_list;
    for(; ei != ei_end; ++ei)
      if( is_edge_valid(*ei, g) )
        ep_list.push_back( g[*ei] );
    std::sort(ep_list.begin(), ep_list.end());
    BOOST_CHECK_EQUAL( ep_list[0], 1003 );
    BOOST_CHECK_EQUAL( ep_list[1], 1004 );
    BOOST_CHECK_EQUAL( ep_list[2], 1005 );
    BOOST_CHECK_EQUAL( ep_list[3], 3010 );
    BOOST_CHECK_EQUAL( ep_list[4], 3011 );
    BOOST_CHECK_EQUAL( ep_list[5], 3012 );
    BOOST_CHECK_EQUAL( ep_list[6], 3013 );
  };
  
  
  
  BOOST_CHECK_NO_THROW( remove_edge(v_rc[1], v_rc2c[2], g) );
  
  BOOST_CHECK_EQUAL( num_vertices(g), 13 );
  {
    VertexIter vi, vi_end;
    BOOST_CHECK_NO_THROW( boost::tie(vi, vi_end) = vertices(g) );
    std::vector<int> vp_list;
    for(; vi != vi_end; ++vi) {
      if( is_vertex_valid(*vi, g) ) {
        vp_list.push_back( g[*vi] );
      };
    };
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 1 );
    BOOST_CHECK_EQUAL( vp_list[1], 2 );
    BOOST_CHECK_EQUAL( vp_list[2], 3 );
    BOOST_CHECK_EQUAL( vp_list[3], 4 );
    BOOST_CHECK_EQUAL( vp_list[4], 5 );
    BOOST_CHECK_EQUAL( vp_list[5], 6 );
    BOOST_CHECK_EQUAL( vp_list[6], 7 );
    BOOST_CHECK_EQUAL( vp_list[7], 8 );
    BOOST_CHECK_EQUAL( vp_list[8], 9 );
    BOOST_CHECK_EQUAL( vp_list[9], 10 );
    BOOST_CHECK_EQUAL( vp_list[10], 11 );
    BOOST_CHECK_EQUAL( vp_list[11], 12 );
    BOOST_CHECK_EQUAL( vp_list[12], 13 );
  };
  
  BOOST_CHECK_EQUAL( num_edges(g), 6 );
  {
    EdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei, ei_end) = edges(g) );
    std::vector<int> ep_list;
    for(; ei != ei_end; ++ei) {
      if( is_edge_valid(*ei, g) ) {
        ep_list.push_back( g[*ei] );
      };
    };
    std::sort(ep_list.begin(), ep_list.end());
    BOOST_CHECK_EQUAL( ep_list[0], 1003 );
    BOOST_CHECK_EQUAL( ep_list[1], 1004 );
    BOOST_CHECK_EQUAL( ep_list[2], 1005 );
    BOOST_CHECK_EQUAL( ep_list[3], 3010 );
    BOOST_CHECK_EQUAL( ep_list[4], 3011 );
    BOOST_CHECK_EQUAL( ep_list[5], 3013 );
  };
  
  
  
  BOOST_CHECK_NO_THROW( remove_edge(e_rc2c[3], g) );
  
  BOOST_CHECK_EQUAL( num_vertices(g), 13 );
  {
    VertexIter vi, vi_end;
    BOOST_CHECK_NO_THROW( boost::tie(vi, vi_end) = vertices(g) );
    std::vector<int> vp_list;
    for(; vi != vi_end; ++vi) {
      if( is_vertex_valid(*vi, g) ) {
        vp_list.push_back( g[*vi] );
      };
    };
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 1 );
    BOOST_CHECK_EQUAL( vp_list[1], 2 );
    BOOST_CHECK_EQUAL( vp_list[2], 3 );
    BOOST_CHECK_EQUAL( vp_list[3], 4 );
    BOOST_CHECK_EQUAL( vp_list[4], 5 );
    BOOST_CHECK_EQUAL( vp_list[5], 6 );
    BOOST_CHECK_EQUAL( vp_list[6], 7 );
    BOOST_CHECK_EQUAL( vp_list[7], 8 );
    BOOST_CHECK_EQUAL( vp_list[8], 9 );
    BOOST_CHECK_EQUAL( vp_list[9], 10 );
    BOOST_CHECK_EQUAL( vp_list[10], 11 );
    BOOST_CHECK_EQUAL( vp_list[11], 12 );
    BOOST_CHECK_EQUAL( vp_list[12], 13 );
  };
  
  BOOST_CHECK_EQUAL( num_edges(g), 5 );
  {
    EdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei, ei_end) = edges(g) );
    std::vector<int> ep_list;
    for(; ei != ei_end; ++ei) {
      if( is_edge_valid(*ei, g) ) {
        ep_list.push_back( g[*ei] );
      };
    };
    std::sort(ep_list.begin(), ep_list.end());
    BOOST_CHECK_EQUAL( ep_list[0], 1003 );
    BOOST_CHECK_EQUAL( ep_list[1], 1004 );
    BOOST_CHECK_EQUAL( ep_list[2], 1005 );
    BOOST_CHECK_EQUAL( ep_list[3], 3010 );
    BOOST_CHECK_EQUAL( ep_list[4], 3011 );
  };
  
  
  
  BOOST_CHECK_NO_THROW( clear_vertex(v_rc2c[0], g) );
  BOOST_CHECK_NO_THROW( remove_vertex(v_rc2c[0], g) );
  
  BOOST_CHECK_EQUAL( num_vertices(g), 12 );
  {
    VertexIter vi, vi_end;
    BOOST_CHECK_NO_THROW( boost::tie(vi, vi_end) = vertices(g) );
    std::vector<int> vp_list;
    for(; vi != vi_end; ++vi) {
      if( is_vertex_valid(*vi, g) ) {
        vp_list.push_back( g[*vi] );
      };
    };
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 1 );
    BOOST_CHECK_EQUAL( vp_list[1], 2 );
    BOOST_CHECK_EQUAL( vp_list[2], 3 );
    BOOST_CHECK_EQUAL( vp_list[3], 4 );
    BOOST_CHECK_EQUAL( vp_list[4], 5 );
    BOOST_CHECK_EQUAL( vp_list[5], 6 );
    BOOST_CHECK_EQUAL( vp_list[6], 7 );
    BOOST_CHECK_EQUAL( vp_list[7], 8 );
    BOOST_CHECK_EQUAL( vp_list[8], 9 );
    BOOST_CHECK_EQUAL( vp_list[9], 11 );
    BOOST_CHECK_EQUAL( vp_list[10], 12 );
    BOOST_CHECK_EQUAL( vp_list[11], 13 );
  };
  
  BOOST_CHECK_EQUAL( num_edges(g), 4 );
  {
    EdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( boost::tie(ei, ei_end) = edges(g) );
    std::vector<int> ep_list;
    for(; ei != ei_end; ++ei) {
      if( is_edge_valid(*ei, g) ) {
        ep_list.push_back( g[*ei] );
      };
    };
    std::sort(ep_list.begin(), ep_list.end());
    BOOST_CHECK_EQUAL( ep_list[0], 1003 );
    BOOST_CHECK_EQUAL( ep_list[1], 1004 );
    BOOST_CHECK_EQUAL( ep_list[2], 1005 );
    BOOST_CHECK_EQUAL( ep_list[3], 3011 );
  };
  
  
};




