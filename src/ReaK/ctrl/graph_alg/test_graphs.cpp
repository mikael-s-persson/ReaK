
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

#include <iostream>
#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/properties.hpp>

#include "path_planning/dvp_layout_adjacency_list.hpp"
#include "topologies/hyperbox_topology.hpp"
#include "lin_alg/vect_alg.hpp"

#include "graph_alg/pooled_adjacency_list.hpp"


typedef ReaK::pp::hyperbox_topology< ReaK::vect<double,2> > TopologyType;
  
typedef TopologyType::point_type PointType;

struct WorldGridVertexProperties {
  PointType pos;
  
  WorldGridVertexProperties(const PointType& aPos = PointType()) : pos(aPos) { };
};

struct WorldGridEdgeProperties {
  double dist;
  
  WorldGridEdgeProperties(double aDist = double(0.0)) : dist(aDist) { };
};

template <typename Graph>
bool test_propgraph_functions(Graph& g) {
  
  typedef typename boost::graph_traits< Graph >::vertex_descriptor Vertex;
  typedef typename boost::graph_traits< Graph >::edge_descriptor Edge;
  
  typedef typename Graph::vertex_property_type VertexProp;
  typedef typename Graph::edge_property_type EdgeProp;
  
  typedef typename boost::property_map< Graph, boost::vertex_raw_prop_to_bundle_t>::type VProp2Bundle;
  typedef typename boost::property_map< Graph, boost::edge_raw_prop_to_bundle_t>::type EProp2Bundle;
  
  VProp2Bundle vp_2_bundle = get(boost::vertex_raw_prop_to_bundle, g);
  EProp2Bundle ep_2_bundle = get(boost::edge_raw_prop_to_bundle, g); RK_UNUSED(ep_2_bundle);
  
  std::vector<WorldGridVertexProperties> pts;
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.0,0.0)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.0,0.2)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.2,0.2)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.6,0.2)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.6,0.4)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.6,0.6)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.6,0.8)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.6,1.0)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.2,1.0)));
  pts.push_back(WorldGridVertexProperties(ReaK::vect<double,2>(0.0,1.0)));
  
  std::vector<Vertex> pts_v;
  for(std::size_t i = 0; i < pts.size(); ++i) {
    VertexProp vp;
    put(vp_2_bundle, vp, pts[i]);
    pts_v.push_back(add_vertex(vp, g));
  };
  
  for(std::size_t i = 0; i < pts_v.size(); ++i) {
    if( g[pts_v[i]].pos != pts[i].pos ) {
      RK_ERROR("The position property of vertex " << i << " was not preserved in the addition!");
      return false;
    };
  };
  
  for(std::size_t i = 0; i < pts_v.size() / 2; ++i) {
    remove_vertex(pts_v[i],g);
  };
  
  for(std::size_t i = pts_v.size() / 2; i < pts_v.size(); ++i) {
    if( get(vp_2_bundle, get(boost::vertex_raw_property, g, pts_v[i])).pos
        != get(boost::vertex_raw_prop_to_bundle, g, pts[i]).pos ) {
      RK_ERROR("The position property of vertex " << i << " was not preserved after the removal!");
      return false;
    };
  };
  
  for(std::size_t i = 0; i < pts_v.size() / 2; ++i) {
    VertexProp vp;
    vp_2_bundle[vp] = pts[i];
    pts_v[i] = add_vertex(std::move(vp),g);
  };
  
  for(std::size_t i = 0; i < pts_v.size(); ++i) {
    if( g[pts_v[i]].pos != pts[i].pos ) {
      RK_ERROR("The position property of vertex " << i << " was not preserved in the second addition!");
      return false;
    };
  };
  
  
  std::vector<WorldGridEdgeProperties> segs;
  segs.push_back(WorldGridEdgeProperties(0.2));  // 0 to 1
  segs.push_back(WorldGridEdgeProperties(0.2));  // 1 to 2
  segs.push_back(WorldGridEdgeProperties(0.4));  // 2 to 3
  segs.push_back(WorldGridEdgeProperties(0.2));  // 3 to 4
  segs.push_back(WorldGridEdgeProperties(0.2));  // 4 to 5
  segs.push_back(WorldGridEdgeProperties(0.2));  // 5 to 6
  segs.push_back(WorldGridEdgeProperties(0.2));  // 6 to 7
  segs.push_back(WorldGridEdgeProperties(0.4));  // 7 to 8
  segs.push_back(WorldGridEdgeProperties(0.2));  // 8 to 9
  segs.push_back(WorldGridEdgeProperties(1.0));  // 9 to 0
  
  std::vector<Edge> segs_e;
  for(std::size_t i = 0; i < segs.size(); ++i) {
    EdgeProp ep;
    put(ep_2_bundle, ep, segs[i]);
    segs_e.push_back(add_edge(pts_v[i], pts_v[(i+1) % pts_v.size()], ep, g).first);
  };
  
  for(std::size_t i = 0; i < segs_e.size(); ++i) {
    if( g[segs_e[i]].dist != segs[i].dist ) {
      RK_ERROR("The dist property of edge " << i << " was not preserved in the edge-addition!");
      return false;
    };
  };
  
  for(std::size_t i = 0; i < segs_e.size(); ++i) {
    if( get(ep_2_bundle, get(boost::edge_raw_property, g, segs_e[i])).dist
        != get(boost::edge_raw_prop_to_bundle, g, segs[i]).dist ) {
      RK_ERROR("The dist property of edge " << i << " was not preserved in the edge-addition!");
      return false;
    };
  };
  
  
  return true;
};


int main() {
  
  ReaK::shared_ptr<TopologyType> m_space = ReaK::shared_ptr<TopologyType>(new TopologyType("",ReaK::vect<double,2>(0.0,0.0),ReaK::vect<double,2>(1.0,1.0)));
  
  typedef boost::data_member_property_map< PointType, WorldGridVertexProperties > PositionMap;
  
  typedef ReaK::pp::dvp_adjacency_list< 
    WorldGridVertexProperties,
    WorldGridEdgeProperties,
    TopologyType, PositionMap,
    2, ReaK::pp::random_vp_chooser,
    ReaK::graph::d_ary_bf_tree_storage<2>,
    boost::vecS, boost::bidirectionalS, boost::listS > WorldPartition2BF;
  
  typedef WorldPartition2BF::adj_list_type WorldGrid2BF;
  
  WorldPartition2BF dvp2(m_space, PositionMap(&WorldGridVertexProperties::pos));
  WorldGrid2BF g2 = dvp2.get_adjacency_list();
  
  if(!test_propgraph_functions(g2))
    return 1;
  
  typedef ReaK::pp::dvp_adjacency_list< 
    WorldGridVertexProperties,
    WorldGridEdgeProperties,
    TopologyType, PositionMap,
    2, ReaK::pp::random_vp_chooser,
    ReaK::graph::d_ary_bf_tree_storage<2>,
    boost::listS, boost::bidirectionalS, boost::listS > WorldPartition2BF_listS;
  
  typedef WorldPartition2BF_listS::adj_list_type WorldGrid2BF_listS;
  
  WorldPartition2BF_listS dvp2_ls(m_space, PositionMap(&WorldGridVertexProperties::pos));
  WorldGrid2BF_listS g2_ls = dvp2_ls.get_adjacency_list();
  
  if(!test_propgraph_functions(g2_ls))
    return 1;
  
  
  typedef boost::pooled_adjacency_list<
    boost::bidirectionalS,
    WorldGridVertexProperties,
    WorldGridEdgeProperties> WorldGridPooled;
  
  WorldGridPooled g_pooled;
  
  if(!test_propgraph_functions(g_pooled))
    return 1;
  
  
  return 0;
};













