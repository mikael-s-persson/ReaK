
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

#ifndef REACHABILITY_SORT_HPP
#define REACHABILITY_SORT_HPP

#include "reachability_space_concept.hpp"

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <map>
#include <set>
#include <algorithm>


namespace ReaK {

namespace pp {


template <typename Graph, typename PositionMap, typename ReachabilityTopology>
class reachability_sorted_set {
  public:
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<Graph>::vertex_iterator VertexIter;
    typedef typename reachability_topology_traits<ReachabilityTopology>::point_type Point;
    
    struct vertex_tuple {
      Vertex u;
      double backward_reach;
      double forward_reach;
      vertex_tuple(Vertex aU, double aBackwardReach, double aForwardReach) : 
                   u(aU), backward_reach(aBackwardReach), forward_reach(aForwardReach) { };
    };
    
    struct backward { };
    struct forward { };
    
    struct vertex_access {
      Vertex& operator()(vertex_tuple& elem) const throw() { return elem.u; };
      const Vertex& operator()(const vertex_tuple& elem) const throw() { return elem.u; };
    };
      
  private:
    Graph& m_g;
    PositionMap m_position;
    const ReachabilityTopology& m_space;
    
    typedef boost::multi_index_container<
      vertex_tuple,
      boost::multi_index::indexed_by<
        boost::multi_index::ordered_non_unique< boost::multi_index::tag<backward>,
          boost::multi_index::member<vertex_tuple,double,&vertex_tuple::backward_reach> >,
	boost::multi_index::ordered_non_unique< boost::multi_index::tag<forward>,
	  boost::multi_index::member<vertex_tuple,double,&vertex_tuple::forward_reach> > > > VertexMultiMap;
   
    typedef VertexMultiMap::index<backward>::type BackwardReachIndex;
    typedef VertexMultiMap::index<forward>::type ForwardReachIndex;
 
    VertexMultiMap m_map;
    
    
  public:
        
    typedef BackwardReachIndex::iterator back_iterator;
    typedef BackwardReachIndex::const_iterator const_back_iterator;
    typedef BackwardReachIndex::reverse_iterator reverse_back_iterator;
    typedef BackwardReachIndex::const_reverse_iterator const_reverse_back_iterator;
        
    typedef ForwardReachIndex::iterator forth_iterator;
    typedef ForwardReachIndex::const_iterator const_forth_iterator;
    typedef ForwardReachIndex::reverse_iterator reverse_forth_iterator;
    typedef ForwardReachIndex::const_reverse_iterator const_reverse_forth_iterator;

    typedef std::pair<const_reverse_back_iterator, const_reverse_back_iterator> backward_range;
    typedef std::pair<const_forth_iterator, const_forth_iterator> forward_range;
    
    
    
    reachability_sorted_set(Graph& g, PositionMap position, const ReachabilityTopology& space) :
                            m_g(g), m_position(position), m_space(space) { 
      boost::function_requires< ReachabilitySpaceConcept<ReachabilityTopology> >();
      VertexIter ui, ui_end;
      for(boost::tie(ui, ui_end) = boost::vertices(m_g); ui != ui_end; ++ui) {
	Point p = boost::get(m_position,*ui);
	m_map.insert(vertex_tuple(*ui, m_space.backward_reach(p), m_space.forward_reach(p)));
      };
      
    };
    
    void swap(reachability_sorted_set<Graph,PositionMap,ReachabilityTopology>& rhs) throw() {
      std::swap(m_map,rhs.m_map); //swap only the map since the other members should be the same.
    };
    
    reachability_sorted_set(const reachability_sorted_set<Graph,PositionMap,ReachabilityTopology>& rhs) :
                            m_g(rhs.m_g), m_position(rhs.m_position), m_space(rhs.m_space), m_map(rhs.m_map) { };
    
    reachability_sorted_set<Graph,PositionMap,ReachabilityTopology>& operator=(const reachability_sorted_set<Graph,PositionMap,ReachabilityTopology>& rhs) {
      reachability_sorted_set<Graph,PositionMap,ReachabilityTopology> tmp(rhs);
      swap(tmp);
      return *this;
    };
    
    void can_reach(const Point& p, std::vector< std::pair<double,Vertex> >& pred_list, 
		   std::size_t max_number = std::numeric_limits<std::size_t>::max(), 
		   double max_radius = std::numeric_limits<double>::infinity()) const {
      double back_p = m_space.backward_reach(p);
      double forth_p = m_space.forward_reach(p);
      double t_p = back_p + forth_p;
      const BackwardReachIndex& back_index = m_map.get<backward>();
      const_back_iterator itb = back_index.begin();
      const_back_iterator itb_end = back_index.upper_bound(back_p); 
      
      for(;itb != itb_end; ++itb) {
	double d_t = t_p - itb->backward_reach - itb->forward_reach;
	if((itb->forward_reach <= forth_p) && (d_t >= 0.0) && (d_t < max_radius)) {
	  double dist = m_space.distance(boost::get(m_position,itb->u),p);
	  if((dist < std::numeric_limits<double>::infinity()) && (dist < max_radius)) {
	    std::pair<double,Vertex> tmp(dist, itb->u);
	    pred_list.insert(std::lower_bound(pred_list.begin(), pred_list.end(), tmp,
					      boost::bind(std::less<double>(),
							  boost::bind(&std::pair<double,Vertex>::first, _1),
							  boost::bind(&std::pair<double,Vertex>::first, _2))), tmp);
	    if(pred_list.size() >= max_number) {
	      pred_list.pop_back();
	      max_radius = pred_list.back().first;
	    };
	  };
	};
      };
      
    };
    
    void reachable_from(const Point& p, std::vector< std::pair<double,Vertex> >& succ_list, 
		        std::size_t max_number = std::numeric_limits<std::size_t>::max(), 
		        double max_radius = std::numeric_limits<double>::infinity()) const {
      double back_p = m_space.backward_reach(p);
      double forth_p = m_space.forward_reach(p);
      double t_p = back_p + forth_p;
      const BackwardReachIndex& back_index = m_map.get<backward>();
      const_back_iterator itb = back_index.lower_bound(back_p); 
      const_back_iterator itb_end = back_index.end();
      
      for(;itb != itb_end; ++itb) {
	double d_t = itb->backward_reach + itb->forward_reach - t_p;
	if((itb->forward_reach >= forth_p) && (d_t >= 0.0) && (d_t < max_radius)) {
	  double dist = m_space.distance(p,boost::get(m_position,itb->u));
	  if((dist < std::numeric_limits<double>::infinity()) && (dist < max_radius)) {
	    std::pair<double,Vertex> tmp(dist, itb->u);
	    succ_list.insert(std::lower_bound(succ_list.begin(), succ_list.end(), tmp,
					      boost::bind(std::less<double>(),
							  boost::bind(&std::pair<double,Vertex>::first, _1),
							  boost::bind(&std::pair<double,Vertex>::first, _2))), tmp);
	    if(succ_list.size() >= max_number) {
	      succ_list.pop_back();
	      max_radius = succ_list.back().first;
	    };
	  };
	};
      };
      
    };
    
    void reachable(const Point& p, std::vector< std::pair<double,Vertex> >& pred_list, 
		                   std::vector< std::pair<double,Vertex> >& succ_list, 
				   std::size_t max_number = std::numeric_limits<std::size_t>::max(), 
		                   double max_radius = std::numeric_limits<double>::infinity()) const {
      double back_p = m_space.backward_reach(p);
      double forth_p = m_space.forward_reach(p);
      double t_p = back_p + forth_p;
      const BackwardReachIndex& back_index = m_map.get<backward>();
      const_back_iterator itb = back_index.begin();
      const_back_iterator itb_end = back_index.upper_bound(back_p); 
      
      double back_max_radius = max_radius;
      for(;itb != itb_end; ++itb) {
	double d_t = t_p - itb->backward_reach - itb->forward_reach;
	if((itb->forward_reach <= forth_p) && (d_t >= 0.0) && (d_t < back_max_radius)) {
	  double dist = m_space.distance(boost::get(m_position,itb->u),p);
	  if((dist < std::numeric_limits<double>::infinity()) && (dist < back_max_radius)) {
	    std::pair<double,Vertex> tmp(dist, itb->u);
	    pred_list.insert(std::lower_bound(pred_list.begin(), pred_list.end(), tmp,
					      boost::bind(std::less<double>(),
							  boost::bind(&std::pair<double,Vertex>::first, _1),
							  boost::bind(&std::pair<double,Vertex>::first, _2))), tmp);
	    if(pred_list.size() >= max_number) {
	      pred_list.pop_back();
	      back_max_radius = pred_list.back().first;
	    };
	  };
	};
      };
      
      itb = itb_end;
      for(;itb != back_index.begin();--itb) {
	if(itb->backward_reach < back_p) {
	  ++itb;
	  break;
	};
      };
      
      itb_end = back_index.end();
      
      double forth_max_radius = max_radius;
      for(;itb != itb_end; ++itb) {
	double d_t = itb->backward_reach + itb->forward_reach - t_p;
	if((itb->forward_reach >= forth_p) && (d_t >= 0.0) && (d_t < forth_max_radius)) {
	  double dist = m_space.distance(p,boost::get(m_position,itb->u));
	  if((dist < std::numeric_limits<double>::infinity()) && (dist < forth_max_radius)) {
	    std::pair<double,Vertex> tmp(dist, itb->u);
	    succ_list.insert(std::lower_bound(succ_list.begin(), succ_list.end(), tmp,
					      boost::bind(std::less<double>(),
							  boost::bind(&std::pair<double,Vertex>::first, _1),
							  boost::bind(&std::pair<double,Vertex>::first, _2))), tmp);
	    if(succ_list.size() >= max_number) {
	      succ_list.pop_back();
	      forth_max_radius = succ_list.back().first;
	    };
	  };
	};
      };
      
    };
    
    
    
    
    
    
    /*********************************** STL SET-like interface *************************************/
    
    typedef Vertex key_type;
    typedef Vertex value_type;
    typedef Vertex& reference;
    typedef const Vertex& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef Vertex* pointer;
    typedef const Vertex* const_pointer;
    
    struct key_compare {
      PositionMap m_position;
      const ReachabilityTopology* m_space;
      key_compare(PositionMap p, const ReachabilityTopology* s) : m_position(p), m_space(s) { };
      bool operator()(Vertex u, Vertex v) const {
	Point p_u = boost::get(m_position,u);
	Point p_v = boost::get(m_position,v);
	return (m_space->backward_reach(p_u) < m_space->backward_reach(p_v));
      };
    };
    typedef key_compare value_compare;
    
    typedef boost::transform_iterator<vertex_access, BackwardReachIndex::iterator> iterator;
    typedef boost::transform_iterator<vertex_access, BackwardReachIndex::const_iterator> const_iterator;
    typedef boost::transform_iterator<vertex_access, BackwardReachIndex::reverse_iterator> reverse_iterator;
    typedef boost::transform_iterator<vertex_access, BackwardReachIndex::const_reverse_iterator> const_reverse_iterator;

    
    iterator begin() { return iterator(m_map.begin(),vertex_access()); };
    const_iterator begin() const { return const_iterator(m_map.begin(),vertex_access()); };
    
    iterator end() { return iterator(m_map.end(),vertex_access()); };
    const_iterator end() const { return const_iterator(m_map.end(),vertex_access()); };
    
    reverse_iterator rbegin() { return reverse_iterator(m_map.rbegin(),vertex_access()); };
    const_reverse_iterator rbegin() const { return const_reverse_iterator(m_map.rbegin(),vertex_access()); };
    
    reverse_iterator rend() { return reverse_iterator(m_map.rend(),vertex_access()); };
    const_reverse_iterator rend() const { return const_reverse_iterator(m_map.rend(),vertex_access()); };
    
    bool empty() const { return m_map.empty(); };
    std::size_t size() const { return m_map.size(); };
    std::size_t max_size() const { return m_map.max_size(); };
    
    std::pair<iterator,bool> insert(Vertex u) {
      Point p = boost::get(m_position,u);
      m_map.insert(vertex_tuple(u, m_space.backward_reach(p), m_space.forward_reach(p)));
    };
    
    iterator insert(iterator pos, Vertex u) {
      Point p = boost::get(m_position,u);
      return iterator(m_map.insert(pos.base(), vertex_tuple(u, m_space.backward_reach(p), m_space.forward_reach(p))),vertex_access());
    };
    
    iterator insert(const_iterator pos, Vertex u) {
      Point p = boost::get(m_position,u);
      return iterator(m_map.insert(pos.base(), vertex_tuple(u, m_space.backward_reach(p), m_space.forward_reach(p))),vertex_access());
    };
    
    template <typename ForwardIter>
    void insert(ForwardIter first, ForwardIter last) {
      for(;first != last; ++first) {
	Point p = boost::get(m_position,*first);
        m_map.insert(vertex_tuple(*first, m_space.backward_reach(p), m_space.forward_reach(p)));
      };
    };

    void erase(iterator pos) {
      m_map.erase(pos.base());
    };
    
    void erase(const_iterator pos) {
      m_map.erase(pos.base());
    };
    
    std::size_t erase(Vertex u) {
      Point p = boost::get(m_position, u);
      BackwardReachIndex::iterator it = m_map.lower_bound(m_space.backward_reach(p));
      std::size_t result = 0;
      while((it != m_map.end()) && (it->u == u)) {
	m_map.erase(it++);
	++result;
      };
      return result;
    };
    
    template <typename ForwardIter>
    void erase(ForwardIter first, ForwardIter last) {
      for(;first != last;++first) erase(*first);
    };
    
    void clear() { m_map.clear(); };
    
    key_compare key_comp() const { return key_compare(m_position,&m_space); };
    value_compare value_comp() const { return key_compare(m_position,&m_space); };
    
    const_iterator find(Vertex u) const {
      Point p = boost::get(m_position, u);
      BackwardReachIndex::const_iterator it = m_map.lower_bound(m_space.backward_reach(p));
      if((it != m_map.end()) && (it->u == u))
	return const_iterator(it,vertex_access());
      return const_iterator(m_map.end(), vertex_access());
    };
    
    std::size_t count(Vertex u) const {
      Point p = boost::get(m_position, u);
      BackwardReachIndex::const_iterator it = m_map.lower_bound(m_space.backward_reach(p));
      std::size_t result = 0;
      while((it != m_map.end()) && (it->u == u)) {
	++result; ++it;
      };
      return result;
    };
    
    const_iterator lower_bound ( Vertex u ) const {
      Point p = boost::get(m_position, u);
      return const_iterator(m_map.lower_bound(m_space.backward_reach(p)),vertex_access());
    };
    
    const_iterator upper_bound ( Vertex u ) const {
      Point p = boost::get(m_position, u);
      return const_iterator(m_map.upper_bound(m_space.backward_reach(p)),vertex_access());
    };
    
    std::pair<const_iterator,const_iterator> equal_range(Vertex u) const {
      Point p = boost::get(m_position, u);
      BackwardReachIndex::const_iterator it = m_map.lower_bound(m_space.backward_reach(p));
      std::pair<const_iterator,const_iterator> result;
      result.first = const_iterator(it,vertex_access());
      while((it != m_map.end()) && (it->u == u)) {
	++it;
      };
      result.second = const_iterator(it,vertex_access());
      return result;
    };
    
    /*********************************** END OF: STL set interface ***********************************/

};



};

};

#endif













