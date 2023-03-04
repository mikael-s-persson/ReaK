/**
 * \file any_knn_synchro.hpp
 *
 * This library defines a type-erasure base-class for K-nearest-neighbor synchronization objects.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
 */

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

#ifndef REAK_ANY_KNN_SYNCHRO_HPP
#define REAK_ANY_KNN_SYNCHRO_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object.hpp>
#include <ReaK/planning/graph_alg/any_graph.hpp>

namespace ReaK::pp {

/**
 * This class can be used as the base for a dynamically polymorphic KNN synchronizer.
 */
class any_knn_synchro {
 public:
  using self = any_knn_synchro;

 protected:
  virtual void added_vertex_impl(graph::any_graph::vertex_descriptor /*unused*/,
                                 graph::any_graph& /*unused*/) const {}
  virtual void removed_vertex_impl(
      graph::any_graph::vertex_descriptor /*unused*/,
      graph::any_graph& /*unused*/) const {}

 public:
  virtual ~any_knn_synchro() = default;

  /**
   * Called to notify the synchronizer that a vertex was just added to the graph.
   * \tparam Vertex The vertex-descriptor type for the graph.
   * \tparam Graph The graph structure type.
   * \param v The vertex that was just added to the graph.
   * \param g The graph to which a vertex was just added.
   */
  template <typename Vertex, typename Graph>
  void added_vertex(Vertex v, Graph& g) const {
    graph::type_erased_graph<Graph> teg(&g);
    this->added_vertex_impl(graph::any_graph::vertex_descriptor(std::any(v)),
                            teg);
  }

  /**
   * Called to notify the synchronizer that a vertex is about to be removed from the graph.
   * \tparam Vertex The vertex-descriptor type for the graph.
   * \tparam Graph The graph structure type.
   * \param v The vertex that is about to be removed from the graph.
   * \param g The graph from which a vertex is about to be removed.
   */
  template <typename Vertex, typename Graph>
  void removed_vertex(Vertex v, Graph& g) const {
    graph::type_erased_graph<Graph> teg(&g);
    this->removed_vertex_impl(graph::any_graph::vertex_descriptor(std::any(v)),
                              teg);
  }
};

/**
 * This class can be used to wrap a generic KNN synchronizer within a dynamically
 * polymorphic KNN synchronizer. This operates on type-erasure via the ReaK::graph::any_graph class.
 * \tparam Graph The graph type.
 * \tparam KNNSynchro The KNN synchronizer object type to be encapsulated by this type-erasure class.
 */
template <typename Graph, typename KNNSynchro>
class type_erased_knn_synchro : public any_knn_synchro {
 public:
  using base_type = any_knn_synchro;
  using self = type_erased_knn_synchro<Graph, KNNSynchro>;

 protected:
  using Vertex = graph::graph_vertex_t<Graph>;

  KNNSynchro synchro;

  void added_vertex_impl(graph::any_graph::vertex_descriptor tev,
                         graph::any_graph& teg) const override {
    synchro.added_vertex(
        std::any_cast<Vertex>(tev),
        static_cast<graph::type_erased_graph<Graph>&>(teg).base());
  }

  void removed_vertex_impl(graph::any_graph::vertex_descriptor tev,
                           graph::any_graph& teg) const override {
    synchro.removed_vertex(
        std::any_cast<Vertex>(tev),
        static_cast<graph::type_erased_graph<Graph>&>(teg).base());
  }

 public:
  explicit type_erased_knn_synchro(KNNSynchro aSynchro = KNNSynchro())
      : any_knn_synchro(), synchro(aSynchro) {}

  ~type_erased_knn_synchro() override = default;
};

}  // namespace ReaK::pp

#endif
