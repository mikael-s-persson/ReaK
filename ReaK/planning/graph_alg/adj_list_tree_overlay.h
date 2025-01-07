/**
 * \file adj_list_tree_overlay.h
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

#ifndef REAK_PLANNING_GRAPH_ALG_ADJ_LIST_TREE_OVERLAY_H_
#define REAK_PLANNING_GRAPH_ALG_ADJ_LIST_TREE_OVERLAY_H_

#include "ReaK/core/base/defs.h"

#include "bagl/adjacency_list.h"
#include "bagl/bfl_d_ary_tree.h"
#include "bagl/graph_concepts.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_traits.h"

#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <vector>

#include "ReaK/planning/graph_alg/tree_organizer_concept.h"

namespace ReaK::graph {

template <typename AdjListOnTreeType>
class alt_tree_view;  // forward-declaration.

template <typename AdjListOnTreeType>
class alt_graph_view;  // forward-declaration.

// Special property tag used to get full adj-list-tree vertex property.
// Only used in partitioning logic to move indexed vertices around.
struct vertex_alt_full_prop_t {
  using kind = bagl::vertex_property_tag;
  static constexpr std::string_view name = "vertex_alt_full_prop";
};
constexpr vertex_alt_full_prop_t vertex_alt_full_prop = {};

namespace alt_detail {

template <typename TreeVertexProperty, typename AdjVertexProperty>
struct tree_vertex_properties_impl {
  using self =
      tree_vertex_properties_impl<TreeVertexProperty, AdjVertexProperty>;
  std::size_t partner_node = std::numeric_limits<std::size_t>::max();
  TreeVertexProperty tree_data;
  AdjVertexProperty user_data;

  tree_vertex_properties_impl(const self&) = default;
  tree_vertex_properties_impl(self&&) noexcept = default;
  self& operator=(const self&) = default;
  self& operator=(self&&) noexcept = default;

  tree_vertex_properties_impl() = default;
  ~tree_vertex_properties_impl() = default;

  tree_vertex_properties_impl(const TreeVertexProperty& rhs) : tree_data(rhs) {}
  tree_vertex_properties_impl(TreeVertexProperty&& rhs) noexcept
      : tree_data(std::move(rhs)) {}

  tree_vertex_properties_impl(const AdjVertexProperty& rhs) : user_data(rhs) {}
  tree_vertex_properties_impl(AdjVertexProperty&& rhs) noexcept
      : user_data(std::move(rhs)) {}

  friend bool operator<(const self& lhs, const self& rhs) {
    return lhs.partner_node < rhs.partner_node;
  }
  friend bool operator<=(const self& lhs, const self& rhs) {
    return lhs.partner_node <= rhs.partner_node;
  }
  friend bool operator>(const self& lhs, const self& rhs) {
    return lhs.partner_node > rhs.partner_node;
  }
  friend bool operator>=(const self& lhs, const self& rhs) {
    return lhs.partner_node >= rhs.partner_node;
  }
  friend bool operator==(const self& lhs, const self& rhs) {
    return lhs.partner_node == rhs.partner_node;
  }
  friend bool operator!=(const self& lhs, const self& rhs) {
    return lhs.partner_node != rhs.partner_node;
  }
};

/*
 * This class template implements an adjacency-list that over-shadows (i.e., a multi-graph)
 * a tree to store the underlying vertices according to the layout that the tree uses (e.g., B-tree (BFL),
 * or COB-tree). An additional benefit of this class template is that the B-tree can act without being an
 * indirect to a graph, and the two graphs are always kept consistent, which would be difficult to do externally.
 */
template <typename AdjListType, typename VertexProperty,
          typename TreeVertexProperty, typename TreeEdgeProperty,
          typename TreeStorageTag>
class adj_list_on_tree {
 public:
  using self = adj_list_on_tree<AdjListType, VertexProperty, TreeVertexProperty,
                                TreeEdgeProperty, TreeStorageTag>;

  using tree_vertex_properties =
      tree_vertex_properties_impl<TreeVertexProperty, VertexProperty>;
  using tree_edge_properties = TreeEdgeProperty;

  using tree_type =
      typename bagl::tree_storage<tree_vertex_properties, tree_edge_properties,
                                  TreeStorageTag>::type;

  using tree_vertex_property_type = TreeVertexProperty;
  using tree_edge_property_type = TreeEdgeProperty;

  using adj_list_type = AdjListType;
  using adj_vertex_descriptor = std::size_t;

  using adj_vertex_property_type = VertexProperty;
  using adj_edge_property_type = bagl::edge_property_type<adj_list_type>;

 private:
  adj_list_type adj_list_;
  tree_type tree_;

 public:
  friend class alt_tree_view<self>;
  friend class alt_graph_view<self>;
};

}  // namespace alt_detail

template <typename OutEdgeListS = bagl::vec_s,
          typename DirectedS = bagl::directed_s,
          typename VertexProperty = bagl::no_property,
          typename AdjEdgeProperty = bagl::no_property,
          typename AdjEdgeListS = bagl::vec_s,
          typename TreeStorageTag = bagl::bfl_d_ary_tree_storage<2>>
struct adj_list_on_tree_tag {

  template <typename TreeVertexProperty, typename TreeEdgeProperty>
  struct alt {
    using tree_vertex_properties =
        alt_detail::tree_vertex_properties_impl<TreeVertexProperty,
                                                VertexProperty>;
    using tree_edge_properties = TreeEdgeProperty;
    using tree_type =
        typename bagl::tree_storage<tree_vertex_properties,
                                    tree_edge_properties, TreeStorageTag>::type;

    using tree_vertex_descriptor = bagl::graph_vertex_descriptor_t<tree_type>;

    struct adj_vertex_properties {
      tree_vertex_descriptor tree_vertex;
    };
    using adj_edge_properties = AdjEdgeProperty;

    using adj_list_type =
        bagl::adjacency_list<OutEdgeListS, bagl::pool_s, DirectedS,
                             adj_vertex_properties, adj_edge_properties>;

    using type = alt_detail::adj_list_on_tree<adj_list_type, VertexProperty,
                                              TreeVertexProperty,
                                              TreeEdgeProperty, TreeStorageTag>;
  };
};

}  // namespace ReaK::graph

namespace bagl {

template <typename TreeVertexProperty, typename TreeEdgeProperty,
          typename OutEdgeListS, typename DirectedS, typename VertexProperty,
          typename AdjEdgeProperty, typename AdjEdgeListS,
          typename TreeStorageTag>
struct tree_storage<TreeVertexProperty, TreeEdgeProperty,
                    ReaK::graph::adj_list_on_tree_tag<
                        OutEdgeListS, DirectedS, VertexProperty,
                        AdjEdgeProperty, AdjEdgeListS, TreeStorageTag>> {
  using alt_tag_type =
      ReaK::graph::adj_list_on_tree_tag<OutEdgeListS, DirectedS, VertexProperty,
                                        AdjEdgeProperty, AdjEdgeListS,
                                        TreeStorageTag>;
  using type = ReaK::graph::alt_tree_view<typename alt_tag_type::template alt<
      TreeVertexProperty, TreeEdgeProperty>::type>;
};

template <typename OutEdgeListS, typename DirectedS, typename VertexProperty,
          typename AdjEdgeProperty, typename AdjEdgeListS,
          typename TreeStorageTag>
struct tree_storage_traits<ReaK::graph::adj_list_on_tree_tag<
    OutEdgeListS, DirectedS, VertexProperty, AdjEdgeProperty, AdjEdgeListS,
    TreeStorageTag>> : tree_storage_traits<TreeStorageTag> {};

}  // namespace bagl

namespace ReaK::graph {

struct alt_tree_view_tag {};

template <typename AdjListOnTreeType>
class alt_tree_view {
 private:
  std::shared_ptr<AdjListOnTreeType> alt_;

 public:
  using self = alt_tree_view<AdjListOnTreeType>;

  using tree_type = typename AdjListOnTreeType::tree_type;

  // Graph traits:
  using vertex_descriptor = bagl::graph_vertex_descriptor_t<tree_type>;
  using edge_descriptor = bagl::graph_edge_descriptor_t<tree_type>;
  using directed_category = bagl::graph_directed_category_t<tree_type>;
  using edge_parallel_category =
      bagl::graph_edge_parallel_category_t<tree_type>;
  using traversal_category = bagl::graph_traversal_category_t<tree_type>;

  static vertex_descriptor null_vertex() {
    return bagl::graph_traits<tree_type>::null_vertex();
  };

  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;

  // PropertyGraph traits:
  using graph_property_type = bagl::no_property;
  using edge_property_type =
      typename AdjListOnTreeType::tree_edge_property_type;
  using vertex_property_type =
      typename AdjListOnTreeType::tree_vertex_property_type;

  using graph_bundled = bagl::no_property;
  using vertex_bundled =
      bagl::lookup_one_property_t<vertex_property_type, bagl::vertex_bundle_t>;
  using edge_bundled =
      bagl::lookup_one_property_t<edge_property_type, bagl::edge_bundle_t>;

  using graph_tag = alt_tree_view_tag;

  friend class alt_graph_view<AdjListOnTreeType>;

  alt_tree_view() : alt_(std::make_shared<AdjListOnTreeType>()) {}

  explicit alt_tree_view(const alt_graph_view<AdjListOnTreeType>& g);

  const tree_type& get_tree() const { return alt_->tree_; }
  tree_type& get_tree() { return alt_->tree_; }

  static auto alt_key_map() {
    using RawVProp = typename AdjListOnTreeType::tree_vertex_properties;
    return bagl::make_data_member_property_map(&RawVProp::partner_node);
  }
  static auto alt_adj_data_map() {
    using RawVProp = typename AdjListOnTreeType::tree_vertex_properties;
    return bagl::make_data_member_property_map(&RawVProp::user_data);
  }
  static auto alt_tree_data_map() {
    using RawVProp = typename AdjListOnTreeType::tree_vertex_properties;
    return bagl::make_data_member_property_map(&RawVProp::tree_data);
  }

  // Bundled Property-map functions (used by the bagl::property_map< self, T Bundle::* > classes).

  vertex_bundled& operator[](const vertex_descriptor& v) {
    return get_property_value(alt_->tree_[v].tree_data, bagl::vertex_bundle);
  }
  const vertex_bundled& operator[](const vertex_descriptor& v) const {
    return get_property_value(alt_->tree_[v].tree_data, bagl::vertex_bundle);
  }
  edge_bundled& operator[](const edge_descriptor& e) {
    return get_property_value(alt_->tree_[e], bagl::edge_bundle);
  }
  const edge_bundled& operator[](const edge_descriptor& e) const {
    return get_property_value(alt_->tree_[e], bagl::edge_bundle);
  }

  auto operator[](bagl::graph_bundle_t /*unused*/) const {
    return bagl::no_property{};
  }

  // MutablePropertyTreeConcept

  template <typename VertexOIter, typename EdgeOIter>
  std::pair<VertexOIter, EdgeOIter> remove_branch_impl(vertex_descriptor v,
                                                       VertexOIter vit_out,
                                                       EdgeOIter eit_out) {
    return remove_branch(v, alt_->tree_, vit_out, eit_out);
  }

  template <typename... VPArgs>
  vertex_descriptor create_root_impl(VPArgs&&... vp_args) {
    auto v_root = create_root(alt_->tree_, std::forward<VPArgs>(vp_args)...);
    alt_->adj_list_[alt_->tree_[v_root].partner_node].tree_vertex = v_root;
    return v_root;
  }

  template <typename VProp>
  std::tuple<vertex_descriptor, edge_descriptor, bool> add_child_impl(
      vertex_descriptor v, VProp&& vp) {
    auto result = add_child(v, alt_->tree_, std::forward<VProp>(vp));
    alt_->adj_list_[alt_->tree_[std::get<0>(result)].partner_node].tree_vertex =
        std::get<0>(result);
    return result;
  }

  template <typename VProp, typename EProp>
  std::tuple<vertex_descriptor, edge_descriptor, bool> add_child_impl(
      vertex_descriptor v, VProp&& vp, EProp&& ep) {
    auto result = add_child(v, alt_->tree_, std::move(vp), std::move(ep));
    alt_->adj_list_[alt_->tree_[std::get<0>(result)].partner_node].tree_vertex =
        std::get<0>(result);
    return result;
  }
};

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
const auto& get(
    const alt_tree_view<AdjListOnTreeType>& g,
    typename alt_tree_view<AdjListOnTreeType>::vertex_descriptor v) {
  return g[v];
}
template <typename AdjListOnTreeType>
auto& get(alt_tree_view<AdjListOnTreeType>& g,
          typename alt_tree_view<AdjListOnTreeType>::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <typename AdjListOnTreeType>
const auto& get(
    const alt_tree_view<AdjListOnTreeType>& g,
    const typename alt_tree_view<AdjListOnTreeType>::edge_descriptor& e) {
  return g[e];
}
template <typename AdjListOnTreeType>
auto& get(alt_tree_view<AdjListOnTreeType>& g,
          const typename alt_tree_view<AdjListOnTreeType>::edge_descriptor& e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <typename AdjListOnTreeType>
const auto& get(const alt_tree_view<AdjListOnTreeType>& g,
                bagl::graph_bundle_t /*unused*/) {
  return g[bagl::graph_bundle];
}
template <typename AdjListOnTreeType>
auto& get(alt_tree_view<AdjListOnTreeType>& g,
          bagl::graph_bundle_t /*unused*/) {
  return g[bagl::graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <typename AdjListOnTreeType, typename VProp>
void put(alt_tree_view<AdjListOnTreeType>& g,
         typename alt_tree_view<AdjListOnTreeType>::vertex_descriptor v,
         VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <typename AdjListOnTreeType, typename EProp>
void put(alt_tree_view<AdjListOnTreeType>& g,
         const typename alt_tree_view<AdjListOnTreeType>::edge_descriptor& e,
         EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
auto& get_property(
    alt_tree_view<AdjListOnTreeType>& g,
    typename alt_tree_view<AdjListOnTreeType>::vertex_descriptor v) {
  return g.get_tree()[v].tree_data;
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
const auto& get_property(
    const alt_tree_view<AdjListOnTreeType>& g,
    typename alt_tree_view<AdjListOnTreeType>::vertex_descriptor v) {
  return g.get_tree()[v].tree_data;
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <typename AdjListOnTreeType>
auto& get_property(
    alt_tree_view<AdjListOnTreeType>& g,
    const typename alt_tree_view<AdjListOnTreeType>::edge_descriptor& e) {
  return get_property(g.get_tree(), e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <typename AdjListOnTreeType>
const auto& get_property(
    const alt_tree_view<AdjListOnTreeType>& g,
    const typename alt_tree_view<AdjListOnTreeType>::edge_descriptor& e) {
  return get_property(g.get_tree(), e);
}

template <typename AdjListOnTreeType>
auto get_property(const alt_tree_view<AdjListOnTreeType>& g,
                  bagl::graph_all_t /*unused*/) {
  return bagl::no_property{};
}

template <typename AdjListOnTreeType, typename Tag>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
get_property(alt_tree_view<AdjListOnTreeType>& /*unused*/,
             Tag /*unused*/) = delete;

template <typename AdjListOnTreeType, typename Tag>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
get_property(const alt_tree_view<AdjListOnTreeType>& /*unused*/,
             Tag /*unused*/) = delete;

template <typename AdjListOnTreeType, typename Tag, typename T>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
set_property(alt_tree_view<AdjListOnTreeType>& /*unused*/, Tag /*unused*/,
             T&& /*unused*/) = delete;

}  // namespace ReaK::graph

namespace bagl {

template <typename AdjListOnTreeType>
struct property_map<ReaK::graph::alt_tree_view<AdjListOnTreeType>,
                    ReaK::graph::vertex_alt_full_prop_t>
    : property_map<
          typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::tree_type,
          vertex_all_t> {};

}  // namespace bagl

namespace ReaK::graph {

template <typename AdjListOnTreeType>
auto get(vertex_alt_full_prop_t /*unused*/,
         alt_tree_view<AdjListOnTreeType>& g) {
  return get(bagl::vertex_all, g.get_tree());
}

template <typename AdjListOnTreeType>
auto get(vertex_alt_full_prop_t /*unused*/,
         const alt_tree_view<AdjListOnTreeType>& g) {
  return get(bagl::vertex_all, g.get_tree());
}

template <typename AdjListOnTreeType, typename Key>
decltype(auto) get(vertex_alt_full_prop_t /*unused*/,
                   const alt_tree_view<AdjListOnTreeType>& g, const Key& k) {
  return get_property(g.get_tree(), k);
}

template <typename AdjListOnTreeType, typename Key>
decltype(auto) get(vertex_alt_full_prop_t /*unused*/,
                   alt_tree_view<AdjListOnTreeType>& g, const Key& k) {
  return get_property(g.get_tree(), k);
}

template <typename AdjListOnTreeType, typename Key, typename Value>
void put(vertex_alt_full_prop_t /*unused*/, alt_tree_view<AdjListOnTreeType>& g,
         const Key& k, Value&& val) {
  get_property(g.get_tree(), k) = std::forward<Value>(val);
}

struct alt_tree_view_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename bagl::property_value<Property, Tag>::type;

    using type = bagl::tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type =
        bagl::tagged_in_property_property_map<value_type, const Graph, Tag>;
  };
};

}  // namespace ReaK::graph

namespace bagl {

/* specializations used by bagl/properties.h */
template <>
struct vertex_property_selector<ReaK::graph::alt_tree_view_tag> {
  using type = ReaK::graph::alt_tree_view_property_selector;
};

template <>
struct edge_property_selector<ReaK::graph::alt_tree_view_tag> {
  using type = ReaK::graph::alt_tree_view_property_selector;
};

template <>
struct graph_property_selector<ReaK::graph::alt_tree_view_tag> {
  using type = ReaK::graph::alt_tree_view_property_selector;
};

}  // namespace bagl

namespace ReaK::graph {

template <typename AdjListOnTreeType, typename Property>
auto get(Property p, alt_tree_view<AdjListOnTreeType>& g) {
  using Map = bagl::property_map_t<alt_tree_view<AdjListOnTreeType>, Property>;
  return Map(&g, p);
}

template <typename AdjListOnTreeType, typename Property>
auto get(Property p, const alt_tree_view<AdjListOnTreeType>& g) {
  using Map =
      bagl::property_map_const_t<alt_tree_view<AdjListOnTreeType>, Property>;
  return Map(&g, p);
}

template <typename AdjListOnTreeType, typename Property, typename Key>
decltype(auto) get(Property p, const alt_tree_view<AdjListOnTreeType>& g,
                   const Key& k) {
  return bagl::get_property_value(get_property(g, k), p);
}

template <typename AdjListOnTreeType, typename Property, typename Key>
decltype(auto) get(Property p, alt_tree_view<AdjListOnTreeType>& g,
                   const Key& k) {
  return bagl::get_property_value(get_property(g, k), p);
}

template <typename AdjListOnTreeType, typename Property, typename Key,
          typename Value>
void put(Property p, alt_tree_view<AdjListOnTreeType>& g, const Key& k,
         Value&& val) {
  bagl::get_property_value(get_property(g, k), p) = std::forward<Value>(val);
}

}  // namespace ReaK::graph

namespace bagl {

template <typename AdjListOnTreeType, typename T, typename Bundle>
struct property_map<ReaK::graph::alt_tree_view<AdjListOnTreeType>,
                    T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v = std::is_convertible_v<
      typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::vertex_bundled*,
      non_const_Bundle*>;
  static constexpr bool is_edge_bundle_v = std::is_convertible_v<
      typename ReaK::graph::alt_tree_view<AdjListOnTreeType>::edge_bundled*,
      non_const_Bundle*>;
  using tag_type = std::conditional_t<
      is_vertex_bundle_v, vertex_bundle_t,
      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<
      non_const_T, ReaK::graph::alt_tree_view<AdjListOnTreeType>, tag_type>;
  using const_type = bundle_member_property_map<
      const non_const_T, const ReaK::graph::alt_tree_view<AdjListOnTreeType>,
      tag_type>;
};

}  // namespace bagl

namespace ReaK::graph {

template <typename AdjListOnTreeType, typename T, typename Bundle>
auto get(T Bundle::*p, alt_tree_view<AdjListOnTreeType>& g) {
  return typename bagl::property_map<alt_tree_view<AdjListOnTreeType>,
                                     T Bundle::*>::type{&g, p};
}

template <typename AdjListOnTreeType, typename T, typename Bundle>
auto get(T Bundle::*p, const alt_tree_view<AdjListOnTreeType>& g) {
  return typename bagl::property_map<alt_tree_view<AdjListOnTreeType>,
                                     T Bundle::*>::const_type{&g, p};
}

template <typename AdjListOnTreeType, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p,
                               const alt_tree_view<AdjListOnTreeType>& g,
                               const Key& k) {
  return (g[k]).*p;
}

template <typename AdjListOnTreeType, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, alt_tree_view<AdjListOnTreeType>& g, const Key& k,
         T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

struct alt_graph_view_tag {};

template <typename AdjListOnTreeType>
class alt_graph_view {
 public:
  using self = alt_graph_view<AdjListOnTreeType>;

  using graph_type = typename AdjListOnTreeType::adj_list_type;
  using tree_type = typename AdjListOnTreeType::tree_type;

  // Graph traits:
  using vertex_descriptor = bagl::graph_vertex_descriptor_t<graph_type>;
  using edge_descriptor = bagl::graph_edge_descriptor_t<graph_type>;
  using directed_category = bagl::graph_directed_category_t<graph_type>;
  using edge_parallel_category =
      bagl::graph_edge_parallel_category_t<graph_type>;
  using traversal_category = bagl::graph_traversal_category_t<graph_type>;

  static vertex_descriptor null_vertex() {
    return bagl::graph_traits<graph_type>::null_vertex();
  }

  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;

  // PropertyGraph traits:
  using graph_property_type = bagl::no_property;
  using vertex_property_type =
      typename AdjListOnTreeType::adj_vertex_property_type;
  using edge_property_type = typename AdjListOnTreeType::adj_edge_property_type;

  using graph_bundled = bagl::no_property;
  using vertex_bundled =
      bagl::lookup_one_property_t<vertex_property_type, bagl::vertex_bundle_t>;
  using edge_bundled =
      bagl::lookup_one_property_t<edge_property_type, bagl::edge_bundle_t>;

  using graph_tag = alt_graph_view_tag;

 private:
  using tree_vertex_desc = bagl::graph_vertex_descriptor_t<tree_type>;
  using tree_vertex_prop = typename tree_type::vertex_property_type;

  struct mutation_visitor_base {

    mutation_visitor_base() = default;
    virtual ~mutation_visitor_base() = default;

    virtual void remove_vertex(tree_vertex_desc,
                               alt_tree_view<AdjListOnTreeType>&) const = 0;
    virtual void add_vertex(alt_tree_view<AdjListOnTreeType>&,
                            const tree_vertex_prop&) const = 0;
    virtual void add_vertex(alt_tree_view<AdjListOnTreeType>&,
                            tree_vertex_prop&&) const = 0;
  };

  template <typename GraphMutationVisitor>
  struct mutation_visitor : mutation_visitor_base {
    GraphMutationVisitor vis;

    static_assert(TreeOrganizerVisitor<GraphMutationVisitor,
                                       alt_tree_view<AdjListOnTreeType>>);

    explicit mutation_visitor(GraphMutationVisitor a_vis) : vis(a_vis) {}
    ~mutation_visitor() override = default;

    void remove_vertex(tree_vertex_desc tv,
                       alt_tree_view<AdjListOnTreeType>& t) const override {
      vis.remove_vertex(tv, t);
    }
    void add_vertex(alt_tree_view<AdjListOnTreeType>& t,
                    const tree_vertex_prop& tvp) const override {
      vis.add_vertex(t, tvp);
    }
    void add_vertex(alt_tree_view<AdjListOnTreeType>& t,
                    tree_vertex_prop&& tvp) const override {
      vis.add_vertex(t, std::move(tvp));
    }
  };

  std::shared_ptr<AdjListOnTreeType> alt_;
  std::shared_ptr<mutation_visitor_base> vis_;

 public:
  friend class alt_tree_view<AdjListOnTreeType>;

  template <typename GraphMutationVisitor>
  explicit alt_graph_view(GraphMutationVisitor vis)
      : alt_(std::make_shared<AdjListOnTreeType>()),
        vis_(std::make_shared<mutation_visitor<GraphMutationVisitor>>(vis)) {}

  template <typename GraphMutationVisitor>
  alt_graph_view(const alt_tree_view<AdjListOnTreeType>& t,
                 GraphMutationVisitor vis);

  // Bundled Property-map functions (used by the bagl::property_map< self, T Bundle::* > classes).

  vertex_bundled& operator[](const vertex_descriptor& v) {
    return bagl::get_property_value(
        alt_->tree_[alt_->adj_list_[v].tree_vertex].user_data,
        bagl::vertex_bundle);
  }
  const vertex_bundled& operator[](const vertex_descriptor& v) const {
    return bagl::get_property_value(
        alt_->tree_[alt_->adj_list_[v].tree_vertex].user_data,
        bagl::vertex_bundle);
  }
  edge_bundled& operator[](const edge_descriptor& e) {
    return bagl::get_property_value(alt_->adj_list_[e], bagl::edge_bundle);
  }
  const edge_bundled& operator[](const edge_descriptor& e) const {
    return bagl::get_property_value(alt_->adj_list_[e], bagl::edge_bundle);
  }

  auto operator[](bagl::graph_bundle_t /*unused*/) const {
    return bagl::no_property{};
  }

  // Returns a reference to the vertex-property associated to the given vertex descriptor.
  vertex_property_type& get_property_impl(vertex_descriptor v) {
    return alt_->tree_[alt_->adj_list_[v].tree_vertex].user_data;
  }

  // Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  const vertex_property_type& get_property_impl(vertex_descriptor v) const {
    return alt_->tree_[alt_->adj_list_[v].tree_vertex].user_data;
  }

  // Returns a reference to the edge-property associated to the given edge descriptor.
  edge_property_type& get_property_impl(const edge_descriptor& e) {
    return get_property(alt_->adj_list_, e);
  }

  // Returns a const-reference to the edge-property associated to the given edge descriptor.
  const edge_property_type& get_property_impl(const edge_descriptor& e) const {
    return get_property(alt_->adj_list_, e);
  }

  // IncidenceGraph concept

  auto out_edges_impl(vertex_descriptor v) const {
    return out_edges(v, alt_->adj_list_);
  }
  vertex_descriptor source_impl(edge_descriptor e) const {
    return source(e, alt_->adj_list_);
  }
  vertex_descriptor target_impl(edge_descriptor e) const {
    return target(e, alt_->adj_list_);
  }
  degree_size_type out_degree_impl(vertex_descriptor v) const {
    return out_degree(v, alt_->adj_list_);
  }

  // BidirectionalGraph concept

  auto in_edges_impl(vertex_descriptor v) const {
    return in_edges(v, alt_->adj_list_);
  }
  degree_size_type in_degree_impl(vertex_descriptor v) const {
    return in_degree(v, alt_->adj_list_);
  }
  degree_size_type degree_impl(edge_descriptor e) const {
    return degree(e, alt_->adj_list_);
  }

  // AdjacencyGraph concept

  auto adjacent_vertices_impl(vertex_descriptor u) const {
    return adjacent_vertices(u, alt_->adj_list_);
  }

  // VertexListGraph concept

  auto vertices_impl() const { return vertices(alt_->adj_list_); }
  vertices_size_type num_vertices_impl() const {
    return num_vertices(alt_->adj_list_);
  }

  // EdgeListGraph concept

  auto edges_impl() const { return edges(alt_->adj_list_); }
  edges_size_type num_edges_impl() const { return num_edges(alt_->adj_list_); }

  // AdjacencyMatrix concept

  std::pair<edge_descriptor, bool> edge_impl(vertex_descriptor u,
                                             vertex_descriptor v) const {
    return edge(u, v, alt_->adj_list_);
  }

  // MutableGraph concept

  void clear_vertex_impl(vertex_descriptor v) {
    clear_vertex(v, alt_->adj_list_);
  }

  void remove_vertex_impl(vertex_descriptor v) {
    alt_tree_view<AdjListOnTreeType> tree_view(*this);
    vis_->remove_vertex(alt_->adj_list_[v].tree_vertex, tree_view);
    clear_vertex(v, alt_->adj_list_);
    remove_vertex(v, alt_->adj_list_);
  }

  void remove_edge_impl(vertex_descriptor u, vertex_descriptor v) {
    remove_edge(u, v, alt_->adj_list_);
  }

  void remove_edge_impl(edge_descriptor e) { remove_edge(e, alt_->adj_list_); }

  // MutablePropertyGraph concept

  template <typename... VPArgs>
  vertex_descriptor add_vertex_impl(VPArgs&&... vp_args) {
    vertex_descriptor v = add_vertex(alt_->adj_list_);
    typename tree_type::vertex_property_type tree_vp;
    tree_vp.partner_node = v;
    tree_vp.user_data = vertex_property_type{std::forward<VPArgs>(vp_args)...};
    alt_tree_view<AdjListOnTreeType> tree_view(*this);
    vis_->add_vertex(tree_view, std::move(tree_vp));
    return v;
  }

  void remove_vertex_impl(vertex_descriptor v, vertex_property_type& vp) {
    vp = alt_->tree_[alt_->adj_list_[v].tree_vertex].user_data;
    remove_vertex_impl(v);
  }

  template <typename... EPArgs>
  std::pair<edge_descriptor, bool> add_edge_impl(vertex_descriptor u,
                                                 vertex_descriptor v,
                                                 EPArgs&&... ep_args) {
    return add_edge(u, v, alt_->adj_list_, std::forward<EPArgs>(ep_args)...);
  }

  template <typename EProp>
  void remove_edge_impl(vertex_descriptor u, vertex_descriptor v, EProp* ep) {
    remove_edge(u, v, alt_->adj_list_, ep);
  }

  template <typename EProp>
  void remove_edge_impl(edge_descriptor e, EProp* ep) {
    remove_edge(e, alt_->adj_list_, ep);
  }

  // Other useful / unclassified functions:

  void update_vertex(vertex_descriptor v) {
    alt_tree_view<AdjListOnTreeType> tree_view(*this);
    typename tree_type::vertex_property_type tree_vp =
        std::move(alt_->tree_[alt_->adj_list_[v].tree_vertex]);
    vis_->remove_vertex(alt_->adj_list_[v].tree_vertex, tree_view);
    vis_->add_vertex(tree_view, std::move(tree_vp));
  }
};

template <typename AdjListOnTreeType>
alt_tree_view<AdjListOnTreeType>::alt_tree_view(
    const alt_graph_view<AdjListOnTreeType>& g)
    : alt_(g.alt_) {}

template <typename AdjListOnTreeType>
template <typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType>::alt_graph_view(
    const alt_tree_view<AdjListOnTreeType>& t, GraphMutationVisitor vis)
    : alt_(t.alt_),
      vis_(
          std::make_shared<typename alt_graph_view<AdjListOnTreeType>::
                               template mutation_visitor<GraphMutationVisitor>>(
              vis)) {}

template <typename AdjListOnTreeType, typename GraphMutationVisitor>
alt_graph_view<AdjListOnTreeType> make_graph_view(
    const alt_tree_view<AdjListOnTreeType>& t, GraphMutationVisitor vis) {
  return alt_graph_view<AdjListOnTreeType>(t, vis);
}

/******************************************************************************************
 * *************************************************************************************
 *                             Tree View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto out_edges(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return out_edges(v, g.get_tree());
}

template <typename AdjListOnTreeType>
auto source(bagl::graph_edge_descriptor_t<alt_tree_view<AdjListOnTreeType>> e,
            const alt_tree_view<AdjListOnTreeType>& g) {
  return source(e, g.get_tree());
}

template <typename AdjListOnTreeType>
auto target(bagl::graph_edge_descriptor_t<alt_tree_view<AdjListOnTreeType>> e,
            const alt_tree_view<AdjListOnTreeType>& g) {
  return target(e, g.get_tree());
}

template <typename AdjListOnTreeType>
auto out_degree(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return out_degree(v, g.get_tree());
}

/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto in_edges(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return in_edges(v, g.get_tree());
}

template <typename AdjListOnTreeType>
auto in_degree(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return in_degree(v, g.get_tree());
}

template <typename AdjListOnTreeType>
auto degree(bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
            const alt_tree_view<AdjListOnTreeType>& g) {
  return degree(v, g.get_tree());
}

/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto adjacent_vertices(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> u,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return adjacent_vertices(u, g.get_tree());
}

/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return vertices(g.get_tree());
}

template <typename AdjListOnTreeType>
auto num_vertices(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_vertices(g.get_tree());
}

/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return edges(g.get_tree());
}

template <typename AdjListOnTreeType>
auto num_edges(const alt_tree_view<AdjListOnTreeType>& g) {
  return num_edges(g.get_tree());
}

/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto edge(bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> u,
          bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
          const alt_tree_view<AdjListOnTreeType>& g) {
  return edge(u, v, g.get_tree());
}

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
auto get_root_vertex(const alt_tree_view<AdjListOnTreeType>& g) {
  return get_root_vertex(g.get_tree());
}

template <typename AdjListOnTreeType>
auto child_vertices(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) {
  return child_vertices(v, g.get_tree());
}

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType>
auto add_child(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    const alt_tree_view<AdjListOnTreeType>& g) = delete;

template <typename AdjListOnTreeType>
void remove_branch(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    alt_tree_view<AdjListOnTreeType>& g) = delete;

/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/

template <typename AdjListOnTreeType, typename... VPArgs>
auto create_root(alt_tree_view<AdjListOnTreeType>& g, VPArgs&&... vp_args) {
  return g.create_root_impl(std::forward<VPArgs>(vp_args)...);
}

template <typename AdjListOnTreeType, typename VProp>
auto add_child(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    alt_tree_view<AdjListOnTreeType>& g, VProp&& vp) {
  return g.add_child_impl(v, std::forward<VProp>(vp));
}

template <typename AdjListOnTreeType, typename VProp, typename EProp>
auto add_child(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    alt_tree_view<AdjListOnTreeType>& g, VProp&& vp, EProp&& ep) {
  return g.add_child_impl(v, vp, ep);
}

template <typename AdjListOnTreeType, typename OutputIter>
OutputIter remove_branch(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    alt_tree_view<AdjListOnTreeType>& g, OutputIter it_out) {
  return g
      .remove_branch_impl(v, it_out,
                          bagl::container_detail::ignore_output_iter())
      .first;
}

template <typename AdjListOnTreeType, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(
    bagl::graph_vertex_descriptor_t<alt_tree_view<AdjListOnTreeType>> v,
    alt_tree_view<AdjListOnTreeType>& g, VertexOIter vit_out,
    EdgeOIter eit_out) {
  return g.remove_branch_impl(v, vit_out, eit_out);
}

/******************************************************************************************
 * *************************************************************************************
 *                             Graph View functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto out_edges(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.out_edges_impl(v);
}

template <typename AdjListOnTreeType>
auto source(bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>> e,
            const alt_graph_view<AdjListOnTreeType>& g) {
  return g.source_impl(e);
}

template <typename AdjListOnTreeType>
auto target(bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>> e,
            const alt_graph_view<AdjListOnTreeType>& g) {
  return g.target_impl(e);
}

template <typename AdjListOnTreeType>
auto out_degree(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.out_degree_impl(v);
}

/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto in_edges(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.in_edges_impl(v);
}

template <typename AdjListOnTreeType>
auto in_degree(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.in_degree_impl(v);
}

template <typename AdjListOnTreeType>
auto degree(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.degree_impl(v);
}

/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto adjacent_vertices(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> u,
    const alt_graph_view<AdjListOnTreeType>& g) {
  return g.adjacent_vertices_impl(u);
}

/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto vertices(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.vertices_impl();
}

template <typename AdjListOnTreeType>
auto num_vertices(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.num_vertices_impl();
}

/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto edges(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.edges_impl();
}

template <typename AdjListOnTreeType>
auto num_edges(const alt_graph_view<AdjListOnTreeType>& g) {
  return g.num_edges_impl();
}

/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
auto edge(bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> u,
          bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
          const alt_graph_view<AdjListOnTreeType>& g) {
  return g.edge_impl(u, v);
}

/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType>
void clear_vertex(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g) {
  g.clear_vertex_impl(v);
}

template <typename AdjListOnTreeType>
void remove_vertex(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_vertex_impl(v);
}

template <typename AdjListOnTreeType>
void remove_edge(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> u,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(u, v);
}

template <typename AdjListOnTreeType>
void remove_edge(
    bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>> e,
    alt_graph_view<AdjListOnTreeType>& g) {
  g.remove_edge_impl(e);
}

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <typename AdjListOnTreeType, typename... VPArgs>
auto add_vertex(alt_graph_view<AdjListOnTreeType>& g, VPArgs&&... vp_args) {
  return g.add_vertex_impl(std::forward<VPArgs>(vp_args)...);
}

template <typename AdjListOnTreeType, typename VProp>
void remove_vertex(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g, VProp* vp) {
  g.remove_vertex_impl(v, vp);
}

template <typename AdjListOnTreeType, typename... EPArgs>
auto add_edge(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> u,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g, EPArgs&&... ep_args) {
  return g.add_edge_impl(u, v, std::forward<EPArgs>(ep_args)...);
}

template <typename AdjListOnTreeType, typename EProp>
void remove_edge(
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> u,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
    alt_graph_view<AdjListOnTreeType>& g, EProp* ep) {
  g.remove_edge_impl(u, v, ep);
}

template <typename AdjListOnTreeType, typename EProp>
void remove_edge(
    bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>> e,
    alt_graph_view<AdjListOnTreeType>& g, EProp* ep) {
  g.remove_edge_impl(e, ep);
}

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
const auto& get(
    const alt_graph_view<AdjListOnTreeType>& g,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v) {
  return g[v];
}
template <typename AdjListOnTreeType>
auto& get(
    alt_graph_view<AdjListOnTreeType>& g,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <typename AdjListOnTreeType>
const auto& get(
    const alt_graph_view<AdjListOnTreeType>& g,
    const bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>>& e) {
  return g[e];
}
template <typename AdjListOnTreeType>
auto& get(
    alt_graph_view<AdjListOnTreeType>& g,
    const bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>>& e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <typename AdjListOnTreeType>
const auto& get(const alt_graph_view<AdjListOnTreeType>& g,
                bagl::graph_bundle_t /*unused*/) {
  return g[bagl::graph_bundle];
}
template <typename AdjListOnTreeType>
auto& get(alt_graph_view<AdjListOnTreeType>& g,
          bagl::graph_bundle_t /*unused*/) {
  return g[bagl::graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <typename AdjListOnTreeType, typename VProp>
void put(alt_graph_view<AdjListOnTreeType>& g,
         bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v,
         VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <typename AdjListOnTreeType, typename EProp>
void put(
    alt_graph_view<AdjListOnTreeType>& g,
    const bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>>& e,
    EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
auto& get_property(
    alt_graph_view<AdjListOnTreeType>& g,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v) {
  return g.get_property_impl(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <typename AdjListOnTreeType>
const auto& get_property(
    const alt_graph_view<AdjListOnTreeType>& g,
    bagl::graph_vertex_descriptor_t<alt_graph_view<AdjListOnTreeType>> v) {
  return g.get_property_impl(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <typename AdjListOnTreeType>
auto& get_property(
    alt_graph_view<AdjListOnTreeType>& g,
    const bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>>& e) {
  return g.get_property_impl(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <typename AdjListOnTreeType>
const auto& get_property(
    const alt_graph_view<AdjListOnTreeType>& g,
    const bagl::graph_edge_descriptor_t<alt_graph_view<AdjListOnTreeType>>& e) {
  return g.get_property_impl(e);
}

template <typename AdjListOnTreeType>
auto get_property(const alt_graph_view<AdjListOnTreeType>& g,
                  bagl::graph_all_t /*unused*/) {
  return bagl::no_property{};
}

template <typename AdjListOnTreeType, typename Tag>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
get_property(alt_graph_view<AdjListOnTreeType>& /*unused*/,
             Tag /*unused*/) = delete;

template <typename AdjListOnTreeType, typename Tag>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
get_property(const alt_graph_view<AdjListOnTreeType>& /*unused*/,
             Tag /*unused*/) = delete;

template <typename AdjListOnTreeType, typename Tag, typename T>
std::enable_if_t<
    std::is_same_v<bagl::property_kind_t<Tag>, bagl::graph_property_tag>>
set_property(alt_graph_view<AdjListOnTreeType>& /*unused*/, Tag /*unused*/,
             T&& /*unused*/) = delete;

struct alt_graph_view_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename bagl::property_value<Property, Tag>::type;

    using type = bagl::tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type =
        bagl::tagged_in_property_property_map<value_type, const Graph, Tag>;
  };
};

}  // namespace ReaK::graph

namespace bagl {

template <typename AdjListOnTreeType, typename T, typename Bundle>
struct property_map<ReaK::graph::alt_graph_view<AdjListOnTreeType>,
                    T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v = std::is_convertible_v<
      typename ReaK::graph::alt_graph_view<AdjListOnTreeType>::vertex_bundled*,
      non_const_Bundle*>;
  static constexpr bool is_edge_bundle_v = std::is_convertible_v<
      typename ReaK::graph::alt_graph_view<AdjListOnTreeType>::edge_bundled*,
      non_const_Bundle*>;
  using tag_type = std::conditional_t<
      is_vertex_bundle_v, vertex_bundle_t,
      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<
      non_const_T, ReaK::graph::alt_graph_view<AdjListOnTreeType>, tag_type>;
  using const_type = bundle_member_property_map<
      const non_const_T, const ReaK::graph::alt_graph_view<AdjListOnTreeType>,
      tag_type>;
};

/* specializations used by bagl/properties.h */
template <>
struct vertex_property_selector<ReaK::graph::alt_graph_view_tag> {
  using type = ReaK::graph::alt_graph_view_property_selector;
};

template <>
struct edge_property_selector<ReaK::graph::alt_graph_view_tag> {
  using type = ReaK::graph::alt_graph_view_property_selector;
};

template <>
struct graph_property_selector<ReaK::graph::alt_graph_view_tag> {
  using type = ReaK::graph::alt_graph_view_property_selector;
};

}  // namespace bagl

namespace ReaK::graph {

template <typename AdjListOnTreeType, typename Property>
auto get(Property p, alt_graph_view<AdjListOnTreeType>& g) {
  using Map = bagl::property_map_t<alt_graph_view<AdjListOnTreeType>, Property>;
  return Map(&g, p);
}

template <typename AdjListOnTreeType, typename Property>
auto get(Property p, const alt_graph_view<AdjListOnTreeType>& g) {
  using Map =
      bagl::property_map_const_t<alt_graph_view<AdjListOnTreeType>, Property>;
  return Map(&g, p);
}

template <typename AdjListOnTreeType, typename Property, typename Key>
decltype(auto) get(Property p, const alt_graph_view<AdjListOnTreeType>& g,
                   const Key& k) {
  return bagl::get_property_value(get_property(g, k), p);
}

template <typename AdjListOnTreeType, typename Property, typename Key>
decltype(auto) get(Property p, alt_graph_view<AdjListOnTreeType>& g,
                   const Key& k) {
  return bagl::get_property_value(get_property(g, k), p);
}

template <typename AdjListOnTreeType, typename Property, typename Key,
          typename Value>
void put(Property p, alt_graph_view<AdjListOnTreeType>& g, const Key& k,
         Value&& val) {
  bagl::get_property_value(get_property(g, k), p) = std::forward<Value>(val);
}

template <typename AdjListOnTreeType, typename T, typename Bundle>
auto get(T Bundle::*p, alt_graph_view<AdjListOnTreeType>& g) {
  return typename bagl::property_map<alt_graph_view<AdjListOnTreeType>,
                                     T Bundle::*>::type{&g, p};
}

template <typename AdjListOnTreeType, typename T, typename Bundle>
auto get(T Bundle::*p, const alt_graph_view<AdjListOnTreeType>& g) {
  return typename bagl::property_map<alt_graph_view<AdjListOnTreeType>,
                                     T Bundle::*>::const_type{&g, p};
}

template <typename AdjListOnTreeType, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p,
                               const alt_graph_view<AdjListOnTreeType>& g,
                               const Key& k) {
  return (g[k]).*p;
}

template <typename AdjListOnTreeType, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, alt_graph_view<AdjListOnTreeType>& g, const Key& k,
         T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_ADJ_LIST_TREE_OVERLAY_H_
