
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


template <typename TreeType>
void print_tree(TreeType& g, typename boost::graph_traits<TreeType>::vertex_descriptor v) {
  typedef typename boost::graph_traits<TreeType>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<TreeType>::edge_descriptor Edge;
  typedef typename boost::graph_traits<TreeType>::out_edge_iterator OutEdgeIter;
  typedef typename boost::graph_traits<TreeType>::in_edge_iterator InEdgeIter;
  typedef typename ReaK::graph::tree_traits<TreeType>::child_vertex_iterator ChildVertIter;
  
  std::cout << "Vertex " << g[v] << std::endl;
  
  OutEdgeIter ei, ei_end;
  for(boost::tie(ei,ei_end) = out_edges(v,g); ei != ei_end; ++ei) {
    if(is_edge_valid(*ei,g))
      std::cout << "\tEdge " << g[*ei] << " = Vertex " 
                             << g[source(*ei,g)] << " --> Vertex " 
                             << g[target(*ei,g)] << std::endl;
  };
  
  ChildVertIter cvi, cvi_end;
  std::vector< Vertex > v_list;
  for(boost::tie(cvi,cvi_end) = child_vertices(v,g); cvi != cvi_end; ++cvi) {
    if(is_vertex_valid(*cvi,g))
      v_list.push_back(*cvi);
  };
  for(std::size_t i = 0; i < v_list.size(); ++i)
    print_tree(g, v_list[i]);
  
};


template <typename TreeType>
void test_tree(TreeType& g) {
  typedef typename boost::graph_traits<TreeType>::vertex_descriptor Vertex;
  typedef typename boost::graph_traits<TreeType>::vertex_iterator VertexIter;
  typedef typename boost::graph_traits<TreeType>::edge_descriptor Edge;
  typedef typename boost::graph_traits<TreeType>::out_edge_iterator OutEdgeIter;
  typedef typename boost::graph_traits<TreeType>::in_edge_iterator InEdgeIter;
  
  Vertex v_root = create_root(g);
  remove_branch(v_root,g);
  int vp_r = 1;
  v_root = create_root(vp_r, g);
  std::cout << "Root vp = " << g[v_root] << std::endl;
  
  std::vector<int> props;
  remove_branch(v_root, back_inserter(props), g);
  for(std::size_t i = 0; i < props.size(); ++i)
    std::cout << "vp " << i << " = " << props[i] << std::endl;
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  v_root = create_root(std::move(vp_r), g);
#else
  v_root = create_root(vp_r, g);
#endif
  props.clear();
  std::cout << "Root vp = " << g[v_root] << std::endl;
  
  int vp_rc[] = {2,3,4,5};
  int ep_rc[] = {1002,1003,1004,1005};
  Vertex v_rc[4];
  Edge e_rc[4];
  for(int i = 0; i < 4; ++i)
    boost::tie(v_rc[i],e_rc[i]) = add_child_vertex(v_root,vp_rc[i],ep_rc[i],g);
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  int vp_rc1c[] = {6,7,8,9};
  int ep_rc1c[] = {2006,2007,2008,2009};
  Vertex v_rc1c[4];
  Edge e_rc1c[4];
  for(std::size_t i = 0; i < 4; ++i)
    boost::tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g);
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  print_tree(g,get_root_vertex(g));
  
  int vp_rc2c[] = {10,11,12,13};
  int ep_rc2c[] = {3010,3011,3012,3013};
  Vertex v_rc2c[4];
  Edge e_rc2c[4];
  for(std::size_t i = 0; i < 4; ++i)
    boost::tie(v_rc2c[i],e_rc2c[i]) = add_child_vertex(v_rc[1],std::move(vp_rc2c[i]),std::move(ep_rc2c[i]),g);
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  print_tree(g,get_root_vertex(g));
  
  remove_branch(v_rc[0],g);
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  print_tree(g,get_root_vertex(g));
  
  remove_branch(v_rc[1], back_inserter(props), g);
  for(std::size_t i = 0; i < props.size(); ++i)
    std::cout << "vp " << i << " = " << props[i] << std::endl;
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  print_tree(g,get_root_vertex(g));
  
  boost::tie(v_rc[0],e_rc[0]) = add_child_vertex(v_root,vp_rc[0],ep_rc[0],g);
  for(std::size_t i = 0; i < 4; ++i)
    boost::tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g);
  std::cout << "Size = " << num_vertices(g) << std::endl;
  
  print_tree(g,get_root_vertex(g));
  
  std::cout << "All vertices are: " << std::endl;
  VertexIter vi, vi_end;
  for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
    if(is_vertex_valid(*vi,g))
      std::cout << g[*vi] << std::endl;
  std::cout << " .. that was all the verties." << std::endl;
  
  std::cout << "Vertex " << vp_rc1c[1] << " has the following in-edges:" << std::endl;
  InEdgeIter iei, iei_end;
  unsigned int j = 0;
  for(boost::tie(iei,iei_end) = in_edges(v_rc1c[1],g); (iei != iei_end) && (j < 4); ++iei, ++j) {
    if(is_edge_valid(*iei,g))
      std::cout << "\tEdge " << g[*iei] << " = Vertex " 
                             << g[source(*iei,g)] << " --> Vertex " 
                             << g[target(*iei,g)] << std::endl;
  };
  
  std::cout << "Edge between Vertex " << vp_rc[0] << " and Vertex " << vp_rc1c[2] << " is found as:" << std::endl;
  Edge e_rc_rc1c = edge(v_rc[0],v_rc1c[2],g).first;
  std::cout << "\tEdge " << g[e_rc_rc1c] << " = Vertex "
                         << g[source(e_rc_rc1c,g)] << " --> Vertex "
			 << g[target(e_rc_rc1c,g)] << std::endl;
  
};



int main() {
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- BF-tree  -----------------------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::d_ary_bf_tree<int, 4, int> bf_tree;
  test_tree(bf_tree);
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- COB-tree  ----------------------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::d_ary_cob_tree<int, 4, int> cob_tree;
  test_tree(cob_tree);
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- Linked-tree (vec,vec) ----------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::linked_tree<boost::vecS, boost::vecS, int, int> lk1_tree;
  test_tree(lk1_tree);
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- Linked-tree (list,vec) ---------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::linked_tree<boost::listS, boost::vecS, int, int> lk2_tree;
  test_tree(lk2_tree);
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- Linked-tree (vec,list) ---------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::linked_tree<boost::vecS, boost::listS, int, int> lk3_tree;
  test_tree(lk3_tree);
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- Linked-tree (list,list) --------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::linked_tree<boost::listS, boost::listS, int, int> lk4_tree;
  test_tree(lk4_tree);
  
  
  std::cout << "----------------------------------------------------------------------" << std::endl
            << "------------------------------- Adj-list (list,list) -----------------" << std::endl
            << "----------------------------------------------------------------------" << std::endl;
  ReaK::graph::tree_storage<int, int>::type adj_list_tree;
  test_tree(adj_list_tree);
  
  
  return 0;
};




