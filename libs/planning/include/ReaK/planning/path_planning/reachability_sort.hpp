/**
 * \file reachability_sort.hpp
 *
 * This library implements a set of vertices whose position lies on a reachability-space
 * and sorts them according to the reachability metrics.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2011
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

#ifndef REAK_REACHABILITY_SORT_HPP
#define REAK_REACHABILITY_SORT_HPP

#include "ReaK/topologies/spaces/reachability_space_concept.hpp"

#include "boost/graph/graph_concepts.hpp"
#include "boost/graph/graph_traits.hpp"

#include "boost/iterator/transform_iterator.hpp"
#include "boost/multi_index/member.hpp"
#include "boost/multi_index/ordered_index.hpp"
#include "boost/multi_index/tag.hpp"
#include "boost/multi_index_container.hpp"

#include <algorithm>
#include <map>
#include <set>

namespace ReaK::pp {

/**
 * This class implements a set of vertices whose position lies on a reachability-space
 * and sorts them according to the reachability metrics. This class allows for very
 * efficient queries of nearest-reachable-neighbors, much more efficient than
 * space partitioning methods (DVP-tree or Kd-tree for example). This class also
 * implements the same interface as the STL set class.
 * \tparam Graph The graph on which the vertices are taken from.
 * \tparam PositionMap The property-map type that can map the vertex descriptors (which should be the value-type of the
 * iterators) to a point (position).
 * \tparam ReachabilityTopology The topology type on which the points can reside, should model the
 * ReachabilitySpaceConcept.
 */
template <typename Graph, typename PositionMap, typename ReachabilityTopology>
class reachability_sorted_set {
 public:
  BOOST_CONCEPT_ASSERT((ReachabilitySpaceConcept<ReachabilityTopology>));

  using self =
      reachability_sorted_set<Graph, PositionMap, ReachabilityTopology>;

  using Vertex = graph::graph_vertex_t<Graph>;
  using Point =
      typename reachability_topology_traits<ReachabilityTopology>::point_type;

  /**
   * This simple POD type stores a vertex descriptor and its backward and forward reach.
   */
  struct vertex_tuple {
    Vertex u;               ///< Holds the vertex descriptor.
    double backward_reach;  ///< Holds the backward reach of the vertex.
    double forward_reach;   ///< Holds the forward reach of the vertex.
                            /**
                            * Parametrized constructor.
                            * \param aU The vertex descriptor.
                            * \param aBackwardReach The backward reach of the vertex.
                            * \param aForwardReach The forward reach of the vertex.
                            */
    vertex_tuple(Vertex aU, double aBackwardReach, double aForwardReach)
        : u(aU), backward_reach(aBackwardReach), forward_reach(aForwardReach) {}
  };

  struct backward {};
  struct forward {};

  struct vertex_access {
    Vertex& operator()(vertex_tuple& elem) const noexcept { return elem.u; }
    const Vertex& operator()(const vertex_tuple& elem) const throw() {
      return elem.u;
    }
  };

 private:
  Graph& m_g;
  PositionMap m_position;
  const ReachabilityTopology& m_space;

  using VertexMultiMap = boost::multi_index_container<
      vertex_tuple,
      boost::multi_index::indexed_by<
          boost::multi_index::ordered_non_unique<
              boost::multi_index::tag<backward>,
              boost::multi_index::member<vertex_tuple, double,
                                         &vertex_tuple::backward_reach>>,
          boost::multi_index::ordered_non_unique<
              boost::multi_index::tag<forward>,
              boost::multi_index::member<vertex_tuple, double,
                                         &vertex_tuple::forward_reach>>>>;

  using BackwardReachIndex = VertexMultiMap::index<backward>::type;
  using ForwardReachIndex = VertexMultiMap::index<forward>::type;

  VertexMultiMap m_map;

 public:
  /** The iterator type to iterate based on the backward reach ordering. */
  using back_iterator = BackwardReachIndex::iterator;
  /** The const-iterator type to iterate based on the backward reach ordering. */
  using const_back_iterator = BackwardReachIndex::const_iterator;
  /** The reverse-iterator type to iterate based on the backward reach ordering. */
  using reverse_back_iterator = BackwardReachIndex::reverse_iterator;
  /** The const-reverse-iterator type to iterate based on the backward reach ordering. */
  using const_reverse_back_iterator =
      BackwardReachIndex::const_reverse_iterator;

  /** The iterator type to iterate based on the forward reach ordering. */
  using forth_iterator = ForwardReachIndex::iterator;
  /** The const-iterator type to iterate based on the forward reach ordering. */
  using const_forth_iterator = ForwardReachIndex::const_iterator;
  /** The reverse-iterator type to iterate based on the forward reach ordering. */
  using reverse_forth_iterator = ForwardReachIndex::reverse_iterator;
  /** The const-reverse-iterator type to iterate based on the forward reach ordering. */
  using const_reverse_forth_iterator =
      ForwardReachIndex::const_reverse_iterator;

  using backward_range =
      std::pair<const_reverse_back_iterator, const_reverse_back_iterator>;
  using forward_range = std::pair<const_forth_iterator, const_forth_iterator>;

  /**
   * Parametrized constructor.
   * \param g The graph from which to take the vertices.
   * \param position The property-map that associates position values to each vertex of the graph.
   * \param space The topology on which the position values reside.
   */
  reachability_sorted_set(Graph& g, PositionMap position,
                          const ReachabilityTopology& space)
      : m_g(g), m_position(position), m_space(space) {
    for (auto [ui, ui_end] = vertices(m_g); ui != ui_end; ++ui) {
      Point p = get(m_position, *ui);
      m_map.insert(vertex_tuple(*ui, m_space.backward_reach(p),
                                m_space.forward_reach(p)));
    }
  }

  /**
   * Standard swap function. Swaps only the sorted set, not the graph, position-map or topology.
   */
  friend void swap(
      reachability_sorted_set<Graph, PositionMap, ReachabilityTopology>&
          rhs) noexcept {
    using std::swap;
    // swap only the map since the other members should be the same.
    std::swap(m_map, rhs.m_map);
  }

  /**
   * Standard copy-constructor.
   */
  reachability_sorted_set(const self& rhs)
      : m_g(rhs.m_g),
        m_position(rhs.m_position),
        m_space(rhs.m_space),
        m_map(rhs.m_map) {}

  /**
   * Standard assignment operator.
   */
  self& operator=(self rhs) {
    swap(*this, rhs);
    return *this;
  }

  /**
   * Obtains the list of all the vertices that can reach the given point (all potential
   * predecessors). Essentially, this function performs a nearest-neighbor queries for
   * all the points that can reach the given point.
   * \param p The point that should be reachable from the predecessors.
   * \param pred_list Stores, as output, the list of predecessors that are the nearest-neighbors (sorted by increasing
   * distance).
   * \param max_number The maximum number of predecessors to find.
   * \param max_radius The maximum reachability radius for the predecessors.
   */
  void can_reach(
      const Point& p, std::vector<std::pair<double, Vertex>>& pred_list,
      std::size_t max_number = std::numeric_limits<std::size_t>::max(),
      double max_radius = std::numeric_limits<double>::infinity()) const {
    double back_p = m_space.backward_reach(p);
    double forth_p = m_space.forward_reach(p);
    double t_p = back_p + forth_p;
    const BackwardReachIndex& back_index = m_map.get<backward>();
    auto itb = back_index.begin();
    auto itb_end = back_index.upper_bound(back_p);

    for (; itb != itb_end; ++itb) {
      double d_t = t_p - itb->backward_reach - itb->forward_reach;
      if ((itb->forward_reach <= forth_p) && (d_t >= 0.0) &&
          (d_t < max_radius)) {
        double dist = m_space.distance(get(m_position, itb->u), p);
        if ((dist < std::numeric_limits<double>::infinity()) &&
            (dist < max_radius)) {
          std::pair tmp(dist, itb->u);
          pred_list.insert(
              std::lower_bound(pred_list.begin(), pred_list.end(), tmp,
                               [](const auto& lhs, const auto& rhs) {
                                 return lhs.first < rhs.first;
                               }),
              tmp);
          if (pred_list.size() >= max_number) {
            pred_list.pop_back();
            max_radius = pred_list.back().first;
          }
        }
      }
    }
  }

  /**
   * Obtains the list of all the vertices that can be reached from the given point (all potential
   * successors). Essentially, this function performs a nearest-neighbor queries for
   * all the points that can be reached from the given point.
   * \param p The point that should be able to reach the successors.
   * \param succ_list Stores, as output, the list of successors that are the nearest-neighbors (sorted by increasing
   * distance).
   * \param max_number The maximum number of successors to find.
   * \param max_radius The maximum reachability radius for the successors.
   */
  void reachable_from(
      const Point& p, std::vector<std::pair<double, Vertex>>& succ_list,
      std::size_t max_number = std::numeric_limits<std::size_t>::max(),
      double max_radius = std::numeric_limits<double>::infinity()) const {
    double back_p = m_space.backward_reach(p);
    double forth_p = m_space.forward_reach(p);
    double t_p = back_p + forth_p;
    const BackwardReachIndex& back_index = m_map.get<backward>();
    auto itb = back_index.lower_bound(back_p);
    auto itb_end = back_index.end();

    for (; itb != itb_end; ++itb) {
      double d_t = itb->backward_reach + itb->forward_reach - t_p;
      if ((itb->forward_reach >= forth_p) && (d_t >= 0.0) &&
          (d_t < max_radius)) {
        double dist = m_space.distance(p, get(m_position, itb->u));
        if ((dist < std::numeric_limits<double>::infinity()) &&
            (dist < max_radius)) {
          std::pair tmp(dist, itb->u);
          succ_list.insert(
              std::lower_bound(succ_list.begin(), succ_list.end(), tmp,
                               [](const auto& lhs, const auto& rhs) {
                                 return lhs.first < rhs.first;
                               }),
              tmp);
          if (succ_list.size() >= max_number) {
            succ_list.pop_back();
            max_radius = succ_list.back().first;
          }
        }
      }
    }
  }

  /**
   * Obtains the list of all the vertices that can reach the given point (all potential predecessors) and
   * be reached from the given point (all potential successors). Essentially, this function performs
   * nearest-neighbor queries in both directions (combines the work of both functions reachable_from
   * and can_reach).
   * \param p The point in question.
   * \param pred_list Stores, as output, the list of predecessors that are the nearest-neighbors (sorted by increasing
   * distance).
   * \param succ_list Stores, as output, the list of successors that are the nearest-neighbors (sorted by increasing
   * distance).
   * \param max_number The maximum number of successors/predecessors to find.
   * \param max_radius The maximum reachability radius for the successors/predecessors.
   */
  void reachable(
      const Point& p, std::vector<std::pair<double, Vertex>>& pred_list,
      std::vector<std::pair<double, Vertex>>& succ_list,
      std::size_t max_number = std::numeric_limits<std::size_t>::max(),
      double max_radius = std::numeric_limits<double>::infinity()) const {
    double back_p = m_space.backward_reach(p);
    double forth_p = m_space.forward_reach(p);
    double t_p = back_p + forth_p;
    const BackwardReachIndex& back_index = m_map.get<backward>();
    const_back_iterator itb = back_index.begin();
    const_back_iterator itb_end = back_index.upper_bound(back_p);

    double back_max_radius = max_radius;
    for (; itb != itb_end; ++itb) {
      double d_t = t_p - itb->backward_reach - itb->forward_reach;
      if ((itb->forward_reach <= forth_p) && (d_t >= 0.0) &&
          (d_t < back_max_radius)) {
        double dist = m_space.distance(get(m_position, itb->u), p);
        if ((dist < std::numeric_limits<double>::infinity()) &&
            (dist < back_max_radius)) {
          std::pair tmp(dist, itb->u);
          pred_list.insert(
              std::lower_bound(pred_list.begin(), pred_list.end(), tmp,
                               [](const auto& lhs, const auto& rhs) {
                                 return lhs.first < rhs.first;
                               }),
              tmp);
          if (pred_list.size() >= max_number) {
            pred_list.pop_back();
            back_max_radius = pred_list.back().first;
          }
        }
      }
    }

    itb = itb_end;
    for (; itb != back_index.begin(); --itb) {
      if (itb->backward_reach < back_p) {
        ++itb;
        break;
      }
    }

    itb_end = back_index.end();

    double forth_max_radius = max_radius;
    for (; itb != itb_end; ++itb) {
      double d_t = itb->backward_reach + itb->forward_reach - t_p;
      if ((itb->forward_reach >= forth_p) && (d_t >= 0.0) &&
          (d_t < forth_max_radius)) {
        double dist = m_space.distance(p, get(m_position, itb->u));
        if ((dist < std::numeric_limits<double>::infinity()) &&
            (dist < forth_max_radius)) {
          std::pair tmp(dist, itb->u);
          succ_list.insert(
              std::lower_bound(succ_list.begin(), succ_list.end(), tmp,
                               [](const auto& lhs, const auto& rhs) {
                                 return lhs.first < rhs.first;
                               }),
              tmp);
          if (succ_list.size() >= max_number) {
            succ_list.pop_back();
            forth_max_radius = succ_list.back().first;
          }
        }
      }
    }
  }

  /*********************************** STL SET-like interface *************************************/

  /** The key-type for the set elements. */
  using key_type = Vertex;
  /** The value-type for the set elements. */
  using value_type = Vertex;
  /** The reference-type for the set elements. */
  using reference = Vertex&;
  /** The const-reference-type for the set elements. */
  using const_reference = const Vertex&;
  /** The size-type for the set. */
  using size_type = std::size_t;
  /** The difference-type for the set. */
  using difference_type = std::ptrdiff_t;
  /** The pointer-type for the set elements. */
  using pointer = Vertex*;
  /** The const-pointer-type for the set elements. */
  using const_pointer = const Vertex*;

  struct key_compare {
    PositionMap m_position;
    const ReachabilityTopology* m_space;
    key_compare(PositionMap p, const ReachabilityTopology* s)
        : m_position(p), m_space(s) {}
    bool operator()(Vertex u, Vertex v) const {
      Point p_u = get(m_position, u);
      Point p_v = get(m_position, v);
      return (m_space->backward_reach(p_u) < m_space->backward_reach(p_v));
    }
  };
  /** The comparison functor for the set elements. */
  using value_compare = key_compare;

  /** The iterator for the set elements. */
  using iterator =
      boost::transform_iterator<vertex_access, BackwardReachIndex::iterator>;
  /** The const-iterator for the set elements. */
  using const_iterator =
      boost::transform_iterator<vertex_access,
                                BackwardReachIndex::const_iterator>;
  /** The reverse-iterator for the set elements. */
  using reverse_iterator =
      boost::transform_iterator<vertex_access,
                                BackwardReachIndex::reverse_iterator>;
  /** The const-reverse-iterator for the set elements. */
  using const_reverse_iterator =
      boost::transform_iterator<vertex_access,
                                BackwardReachIndex::const_reverse_iterator>;

  /**
   * Returns the iterator at the begining of the set.
   */
  auto begin() { return iterator(m_map.begin(), vertex_access()); }
  /**
   * Returns the const-iterator at the begining of the set.
   */
  auto begin() const { return const_iterator(m_map.begin(), vertex_access()); }

  /**
   * Returns the iterator at the one-past-last element of the set.
   */
  auto end() { return iterator(m_map.end(), vertex_access()); }
  /**
   * Returns the const-iterator at the one-past-last element of the set.
   */
  auto end() const { return const_iterator(m_map.end(), vertex_access()); }

  /**
   * Returns the reverse-iterator at the begining of the set (actually the end).
   */
  auto rbegin() { return reverse_iterator(m_map.rbegin(), vertex_access()); }
  /**
   * Returns the const-reverse-iterator at the begining of the set (actually the end).
   */
  auto rbegin() const {
    return const_reverse_iterator(m_map.rbegin(), vertex_access());
  }

  /**
   * Returns the reverse-iterator at the one-past-last element of the set (actually the begin).
   */
  auto rend() { return reverse_iterator(m_map.rend(), vertex_access()); }
  /**
   * Returns the const-reverse-iterator at the one-past-last element of the set (actually the begin).
   */
  auto rend() const {
    return const_reverse_iterator(m_map.rend(), vertex_access());
  }

  /**
   * Checks if the set is empty.
   * \return True if the set is empty.
   */
  bool empty() const { return m_map.empty(); }
  /**
   * Returns the size of the set.
   */
  std::size_t size() const { return m_map.size(); }
  /**
   * Returns the max-size of the set.
   */
  std::size_t max_size() const { return m_map.max_size(); }

  /**
   * Inserts a value into the set.
   * \param u The value to be added to the set.
   * \return A pair that contains the iterator to the added value and a bool telling whether the insertion was
   * successful.
   */
  std::pair<iterator, bool> insert(Vertex u) {
    Point p = get(m_position, u);
    auto [it, insert_succeeded] = m_map.insert(
        vertex_tuple(u, m_space.backward_reach(p), m_space.forward_reach(p)));
    return {iterator(it, vertex_access()), insert_succeeded};
  }

  /**
   * Inserts a value into the set at the given position.
   * \param pos The position in the set where to add the value.
   * \param u The value to be added to the set.
   * \return The iterator to the added value.
   */
  iterator insert(iterator pos, Vertex u) {
    Point p = get(m_position, u);
    return iterator(
        m_map.insert(pos.base(), vertex_tuple(u, m_space.backward_reach(p),
                                              m_space.forward_reach(p))),
        vertex_access());
  }

  /**
   * Inserts a value into the set at the given position.
   * \param pos The position in the set where to add the value.
   * \param u The value to be added to the set.
   * \return The iterator to the added value.
   */
  iterator insert(const_iterator pos, Vertex u) {
    Point p = get(m_position, u);
    return iterator(
        m_map.insert(pos.base(), vertex_tuple(u, m_space.backward_reach(p),
                                              m_space.forward_reach(p))),
        vertex_access());
  }

  /**
   * Inserts a range of values into the set.
   * \param first The first element to add to the set.
   * \param last The one-past-last element to add to the set.
   */
  template <typename ForwardIter>
  void insert(ForwardIter first, ForwardIter last) {
    for (; first != last; ++first) {
      Point p = get(m_position, *first);
      m_map.insert(vertex_tuple(*first, m_space.backward_reach(p),
                                m_space.forward_reach(p)));
    }
  }

  /**
   * Removes a value at a given position from the set.
   * \param pos The position of the element to be removed.
   */
  void erase(iterator pos) { m_map.erase(pos.base()); }

  /**
   * Removes a value at a given position from the set.
   * \param pos The position of the element to be removed.
   */
  void erase(const_iterator pos) { m_map.erase(pos.base()); }

  /**
   * Removes a value from the set.
   * \param u The element to be removed.
   * \return The number of elements removed from the set.
   */
  std::size_t erase(Vertex u) {
    Point p = get(m_position, u);
    auto it = m_map.lower_bound(m_space.backward_reach(p));
    std::size_t result = 0;
    while ((it != m_map.end()) && (it->u == u)) {
      m_map.erase(it++);
      ++result;
    }
    return result;
  }

  /**
   * Removes a range of values from the set.
   * \tparam ForwardIter A forward-iterator type to access the range of values to remove.
   * \param first The first element to remove to the set.
   * \param last The one-past-last element to remove to the set.
   */
  template <typename ForwardIter>
  void erase(ForwardIter first, ForwardIter last) {
    for (; first != last; ++first) {
      erase(*first);
    }
  }

  /**
   * Clears the set.
   */
  void clear() { m_map.clear(); }

  /**
   * Returns the key-comparison functor.
   */
  key_compare key_comp() const { return key_compare(m_position, &m_space); }
  /**
   * Returns the value-comparison functor.
   */
  value_compare value_comp() const { return key_compare(m_position, &m_space); }

  /**
   * Finds a given value in the set.
   * \param u The value to be found in the set.
   * \return The iterator to the given value in the set.
   */
  iterator find(Vertex u) {
    Point p = get(m_position, u);
    auto it = m_map.lower_bound(m_space.backward_reach(p));
    if ((it != m_map.end()) && (it->u == u)) {
      return iterator(it, vertex_access());
    }
    return iterator(m_map.end(), vertex_access());
  }

  /**
   * Finds a given value in the set.
   * \param u The value to be found in the set.
   * \return The iterator to the given value in the set.
   */
  const_iterator find(Vertex u) const {
    Point p = get(m_position, u);
    auto it = m_map.lower_bound(m_space.backward_reach(p));
    if ((it != m_map.end()) && (it->u == u)) {
      return const_iterator(it, vertex_access());
    }
    return const_iterator(m_map.end(), vertex_access());
  }

  /**
   * Count the number of occurrences of a given value in the set.
   * \param u The value to be found in the set.
   * \return The number of occurrences of the given value in the set.
   */
  std::size_t count(Vertex u) const {
    Point p = get(m_position, u);
    auto it = m_map.lower_bound(m_space.backward_reach(p));
    std::size_t result = 0;
    while ((it != m_map.end()) && (it->u == u)) {
      ++result;
      ++it;
    }
    return result;
  }

  /**
   * Finds the first occurrence of a value in the set.
   * \param u The value to be found in the set.
   * \return The iterator to the first occurrence of a value in the set.
   */
  const_iterator lower_bound(Vertex u) const {
    Point p = get(m_position, u);
    return const_iterator(m_map.lower_bound(m_space.backward_reach(p)),
                          vertex_access());
  }

  /**
   * Finds the one-past-last occurrence of a value in the set.
   * \param u The value to be found in the set.
   * \return The iterator to the one-past-last occurrence of a value in the set.
   */
  const_iterator upper_bound(Vertex u) const {
    Point p = get(m_position, u);
    return const_iterator(m_map.upper_bound(m_space.backward_reach(p)),
                          vertex_access());
  }

  /**
   * Finds the range of (first,one-past-last) occurrence of a value in the set
   * \param u The value to be found in the set.
   * \return The range to the (first,one-past-last) occurrence of a value in the set.
   */
  std::pair<const_iterator, const_iterator> equal_range(Vertex u) const {
    Point p = get(m_position, u);
    auto it = m_map.lower_bound(m_space.backward_reach(p));
    auto first = const_iterator(it, vertex_access());
    while ((it != m_map.end()) && (it->u == u)) {
      ++it;
    }
    auto last = const_iterator(it, vertex_access());
    return {first, last};
  }

  /*********************************** END OF: STL set interface ***********************************/
};

}  // namespace ReaK::pp

#endif
