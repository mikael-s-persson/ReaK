/**
 * \file avl_tree_detail.h
 *
 * This library implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN)
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 *
 * \todo Currently, there is no support for set/map (i.e., there is no enforcement of unique-keys).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date April 2013
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

#ifndef REAK_PLANNING_GRAPH_ALG_AVL_TREE_DETAIL_H_
#define REAK_PLANNING_GRAPH_ALG_AVL_TREE_DETAIL_H_

#include "ReaK/core/base/defs.h"

#include "bagl/bst_inorder_iterator.h"
#include "bagl/graph_concepts.h"
#include "bagl/property_map.h"
#include "bagl/tree_adaptor.h"

#include <iterator>
#include <queue>
#include <vector>

namespace ReaK::graph {

namespace avl_detail {

enum avl_container_style {
  avl_set_style,
  avl_multiset_style,
  avl_map_style,
  avl_multimap_style
};

template <typename TreeType, typename Compare,
          avl_container_style ContainerStyle>
struct avl_tree_helper {
  using value_type = bagl::vertex_bundle_type<TreeType>;
  using key_type = value_type;
  using mapped_type = value_type;

  using value_compare = Compare;
  using key_compare = Compare;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  static const bool allow_duplicates = false;

  static const key_type& value_to_key(const value_type& rhs) { return rhs; }

  static const mapped_type& value_to_mapped(const value_type& rhs) {
    return rhs;
  }
  static mapped_type& value_to_mapped(value_type& rhs) { return rhs; }

  template <typename U, typename V>
  static value_type keymap_to_value(U&& k, V&& /*unused*/) {
    return value_type(std::forward<U>(k));
  }
  static value_type keymap_to_value(const key_type& k,
                                    const mapped_type& /*unused*/) {
    return value_type(k);
  }
};

template <typename TreeType, typename Compare>
struct avl_tree_helper<TreeType, Compare, avl_multiset_style> {
  using value_type = bagl::vertex_bundle_type<TreeType>;
  using key_type = value_type;
  using mapped_type = value_type;

  using value_compare = Compare;
  using key_compare = Compare;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  static const bool allow_duplicates = true;

  static const key_type& value_to_key(const value_type& rhs) { return rhs; }

  static const mapped_type& value_to_mapped(const value_type& rhs) {
    return rhs;
  }
  static mapped_type& value_to_mapped(value_type& rhs) { return rhs; }

  template <typename U, typename V>
  static value_type keymap_to_value(U&& k, V&& /*unused*/) {
    return value_type(std::forward<U>(k));
  }
  static value_type keymap_to_value(const key_type& k,
                                    const mapped_type& /*unused*/) {
    return value_type(k);
  }
};

template <typename TreeType, typename Compare>
struct avl_tree_helper<TreeType, Compare, avl_map_style> {
  using value_type = bagl::vertex_bundle_type<TreeType>;
  using key_type = typename value_type::first_type;
  using mapped_type = typename value_type::second_type;

  struct value_compare {
    Compare compare_;

    explicit value_compare(Compare aComp) : compare_(aComp) {}

    bool operator()(const value_type& lhs, const value_type& rhs) const {
      return compare_(lhs.first, rhs.first);
    }
    bool operator()(const value_type& lhs, const key_type& rhs) const {
      return compare_(lhs.first, rhs);
    }
    bool operator()(const key_type& lhs, const value_type& rhs) const {
      return compare_(lhs, rhs.first);
    }
    bool operator()(const key_type& lhs, const key_type& rhs) const {
      return compare_(lhs, rhs);
    }
  };
  using key_compare = value_compare;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  static const bool allow_duplicates = false;

  static const key_type& value_to_key(const value_type& rhs) {
    return rhs.first;
  }

  static const mapped_type& value_to_mapped(const value_type& rhs) {
    return rhs.second;
  }
  static mapped_type& value_to_mapped(value_type& rhs) { return rhs.second; }

  template <typename U, typename V>
  static value_type keymap_to_value(U&& k, V&& m) {
    return value_type(std::forward<U>(k), std::forward<V>(m));
  }
  static value_type keymap_to_value(const key_type& k, const mapped_type& m) {
    return value_type(k, m);
  }
};

template <typename TreeType, typename Compare>
struct avl_tree_helper<TreeType, Compare, avl_multimap_style> {
  using value_type = bagl::vertex_bundle_type<TreeType>;
  using key_type = typename value_type::first_type;
  using mapped_type = typename value_type::second_type;

  struct value_compare {
    Compare compare_;

    explicit value_compare(Compare aComp) : compare_(aComp){};

    bool operator()(const value_type& lhs, const value_type& rhs) const {
      return compare_(lhs.first, rhs.first);
    }
    bool operator()(const value_type& lhs, const key_type& rhs) const {
      return compare_(lhs.first, rhs);
    }
    bool operator()(const key_type& lhs, const value_type& rhs) const {
      return compare_(lhs, rhs.first);
    }
    bool operator()(const key_type& lhs, const key_type& rhs) const {
      return compare_(lhs, rhs);
    }
  };
  using key_compare = value_compare;

  using reference = value_type&;
  using const_reference = const value_type&;
  using pointer = value_type*;
  using const_pointer = const value_type*;

  static const bool allow_duplicates = true;

  static const key_type& value_to_key(const value_type& rhs) {
    return rhs.first;
  }

  static const mapped_type& value_to_mapped(const value_type& rhs) {
    return rhs.second;
  }
  static mapped_type& value_to_mapped(value_type& rhs) { return rhs.second; }

  template <typename U, typename V>
  static value_type keymap_to_value(U&& k, V&& m) {
    return value_type(std::forward<U>(k), std::forward<V>(m));
  }
  static value_type keymap_to_value(const key_type& k, const mapped_type& m) {
    return value_type(k, m);
  }
};

}  // namespace avl_detail

/**
 * This class implements the details of a AVL Tree that
 * allows for O(logN) time look-ups in a strict weak ordering, with amortized O(logN)
 * insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * \tparam TreeType The tree type to be used to store the entries of this AVL tree.
 * \tparam Compare The comparison functor type that can compare two elements.
 */
template <typename TreeType, typename Compare,
          avl_detail::avl_container_style ContainerStyle>
class avl_tree_impl {
 public:
  using self = avl_tree_impl<TreeType, Compare, ContainerStyle>;

  using helper_type =
      avl_detail::avl_tree_helper<TreeType, Compare, ContainerStyle>;

  using value_type = typename helper_type::value_type;
  using key_type = typename helper_type::key_type;
  using mapped_type = typename helper_type::mapped_type;
  using value_compare = typename helper_type::value_compare;
  using key_compare = typename helper_type::key_compare;

  using reference = typename helper_type::reference;
  using const_reference = typename helper_type::const_reference;
  using pointer = typename helper_type::pointer;
  using const_pointer = typename helper_type::const_pointer;

  using iterator = bagl::bst_inorder_iterator<const TreeType, const value_type>;
  using const_iterator = iterator;

  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using allocator_type = std::allocator<value_type>;

 private:
  using tree_indexer = TreeType;

  using vertex_type = bagl::graph_vertex_descriptor_t<tree_indexer>;
  using edge_type = bagl::graph_edge_descriptor_t<tree_indexer>;

  using vertex_property = bagl::vertex_property_type<tree_indexer>;
  using edge_property = bagl::edge_property_type<tree_indexer>;

  tree_indexer tree_;      ///< Tree storage.
  vertex_type root_;       ///< Root node of the tree.
  value_compare compare_;  ///< The comparison functor.

  struct construction_task {
    using PropIter = typename std::vector<vertex_property>::iterator;

    vertex_type node;
    PropIter first;
    PropIter last;
  };

  /* NOTE Invalidates vertices */
  /* NOTE This is a non-recursive version of the construct-node algorithm */
  /* Does not require persistent vertices */
  /* This is the main tree construction function. It takes the vertices in the iterator range and organizes them
   * as a sub-tree below the aParentNode node. */
  void construct_node(vertex_type node,
                      typename construction_task::PropIter first,
                      typename construction_task::PropIter last) {
    using std::swap;

    std::queue<construction_task> tasks;
    tasks.emplace(node, first, last);

    // breadth-first construction of the tree from node and everything below it.
    while (!tasks.empty()) {
      construction_task cur_task = tasks.front();
      tasks.pop();

      // partition the elements around the median element:
      auto cur_task_median =
          cur_task.first + std::distance(cur_task.first, cur_task.last) / 2;
      std::nth_element(cur_task.first, cur_task_median, cur_task.last,
                       compare_);

      tree_[cur_task.node] = std::move(*cur_task_median);

      // first take care of the [first, median) range.
      if (cur_task.first == cur_task_median) {
        continue;  // this is a leaf node.
      }
      auto [new_v, new_e, new_added] = add_child(cur_task.node, tree_);
      tasks.emplace(new_v, cur_task.first, cur_task_median);

      // then, take care of the (median, last) range.
      ++cur_task_median;  // exclude the median.
      if (cur_task_median == cur_task.last) {
        continue;  // there is no right-child.
      }
      std::tie(new_v, new_e, new_added) = add_child(cur_task.node, tree_);
      tasks.emplace(new_v, cur_task_median, cur_task.last);
    }
  }

  bool has_left_child(vertex_type u) const { return out_degree(u, tree_) >= 1; }
  bool has_right_child(vertex_type u) const {
    return out_degree(u, tree_) >= 2;
  }

  vertex_type get_left_child(vertex_type u) const {
    return *(children(u, tree_).begin());
  }
  vertex_type get_right_child(vertex_type u) const {
    return *(std::next(children(u, tree_).begin()));
  }

  vertex_type find_lower_bound(const key_type& k, vertex_type cur_node) const {
    vertex_type cur_lb = tree_indexer::null_vertex();
    while (true) {
      if (!compare_(tree_[cur_node], k)) {
        cur_lb = cur_node;
        // descend on the left side.
        if (!has_left_child(cur_node)) {
          break;
        }
        cur_node = get_left_child(cur_node);
      } else {
        // descend on the right side.
        if (!has_right_child(cur_node)) {
          break;
        }
        cur_node = get_right_child(cur_node);
      }
    }
    return cur_lb;
  }

  vertex_type find_upper_bound(const key_type& k, vertex_type cur_node) const {
    vertex_type cur_ub = tree_indexer::null_vertex();
    while (true) {
      if (compare_(k, tree_[cur_node])) {
        cur_ub = cur_node;
        // descend on the left side.
        if (!has_left_child(cur_node)) {
          break;
        }
        cur_node = get_left_child(cur_node);
      } else {
        // descend on the right side.
        if (!has_right_child(cur_node)) {
          break;
        }
        cur_node = get_right_child(cur_node);
      }
    }
    return cur_ub;
  }

  std::pair<std::size_t, std::size_t> get_minmax_depth(vertex_type node) const {
    std::pair<std::size_t, std::size_t> result(
        std::numeric_limits<std::size_t>::max(), 0);

    // assume that a breadth-first traversal is more efficient:
    std::queue<std::pair<vertex_type, std::size_t>> tasks;
    tasks.emplace(node, 0);
    while (!tasks.empty()) {
      auto [cur_node, cur_depth] = tasks.front();
      tasks.pop();
      if (!has_left_child(cur_node)) {
        if (cur_depth < result.first) {
          result.first = cur_depth;
        }
        if (cur_depth > result.second) {
          result.second = cur_depth;
        }
        continue;
      };
      ++cur_depth;
      tasks.emplace(get_left_child(cur_node), cur_depth);
      if (has_right_child(cur_node)) {
        tasks.emplace(get_right_child(cur_node), cur_depth);
      }
    }

    return result;
  }

  // collect all the keys below and including the given root-node. Appends all the keys to the list.
  void collect_nodes(vertex_type root_node,
                     std::vector<vertex_property>& vp_list) {
    // assume that a breadth-first traversal is more efficient:
    std::queue<vertex_type> tasks;
    tasks.push(root_node);
    while (!tasks.empty()) {
      vertex_type cur_task = tasks.front();
      tasks.pop();
      vp_list.emplace_back(std::move(tree_[cur_task]));
      if (!has_left_child(cur_task)) {
        continue;
      }
      tasks.push(get_left_child(cur_task));
      if (has_right_child(cur_task)) {
        tasks.push(get_right_child(cur_task));
      }
    }
    if (has_right_child(root_node)) {
      remove_branch(get_right_child(root_node), tree_);
    }
    if (has_left_child(root_node)) {
      remove_branch(get_left_child(root_node), tree_);
    }
  }

  // re-balance the sub-tree rooted at the given node. Uses a "collect and re-construct" approach.
  void rebalance_subtree(vertex_type root_node) {
    // collect and remove.
    std::vector<vertex_property> collected_nodes;
    collect_nodes(root_node, collected_nodes);
    // re-construct:
    construct_node(root_node, collected_nodes.begin(), collected_nodes.end());
  }

  // re-balance the sub-tree rooted at the given node. Uses a "collect and re-construct" approach.
  void rebalance_and_insert_subtree(vertex_type root_node,
                                    vertex_property&& added) {
    // collect and remove.
    std::vector<vertex_property> collected_nodes;
    collected_nodes.emplace_back(std::move(added));
    collect_nodes(root_node, collected_nodes);
    // re-construct:
    construct_node(root_node, collected_nodes.begin(), collected_nodes.end());
  }

  std::pair<vertex_type, bool> find_imbalance_impl(vertex_type u,
                                                   int depth_change) const {
    // if we were to change the depth below u by the given depth-change, then would it
    // create an imbalance somewhere? if yes, what vertex stems the sub-tree that must be re-balanced?
    std::pair<std::size_t, std::size_t> c_depth = get_minmax_depth(u);
    // assume c_depth.first == c_depth.second, i.e., the sub-tree is full and balanced (otherwise, what's the point of
    // calling this function):
    c_depth.first += depth_change;
    c_depth.second += depth_change;
    vertex_type orig_u = u;

    auto ei_rg = in_edges(u, tree_);
    auto ei = ei_rg.begin();
    while (ei != ei_rg.end()) {  // until you hit the root node.
      vertex_type p = source(*ei, tree_);
      auto vi_rg = children(p, tree_);
      std::pair<std::size_t, std::size_t> o_depth;
      if (*vi_rg.begin() != u) {
        o_depth = get_minmax_depth(*vi_rg.begin());
      } else if (out_degree(p, tree_) > 1) {
        o_depth = get_minmax_depth(*std::next(vi_rg.begin()));
      } else {
        o_depth = std::pair<std::size_t, std::size_t>(0, 0);
      }

      // check if the other branch already has compatible depth-bounds:
      if ((c_depth.first >= o_depth.first) &&
          (c_depth.second <= o_depth.second)) {
        return {orig_u, false};
      }

      // if not, then check if the other branch already has incompatible depth-bounds:
      if ((o_depth.first < c_depth.second - 1) ||
          (o_depth.second > c_depth.first + 1)) {
        // then, the sub-tree spanned from "p" will be unbalanced:
        return {p, true};
      }

      // otherwise, keep on looking:
      if (o_depth.first < c_depth.first) {
        c_depth.first = o_depth.first;
      }
      if (o_depth.second > c_depth.second) {
        c_depth.second = o_depth.second;
      }
      c_depth.first += 1;  // move up one level.
      c_depth.second += 1;
      u = p;
      ei_rg = in_edges(u, tree_);
      ei = ei_rg.begin();
    }
    // could not find an imbalance anywhere.
    return {orig_u, false};
  }

  std::pair<iterator, bool> insert_after_terminal_impl(vertex_type before_u,
                                                       vertex_property&& vp) {
    using std::swap;
    if (out_degree(before_u, tree_)) {
      // then the before node is terminal but not a leaf, safe to insert:
      auto [new_v, new_e, new_added] =
          add_child(before_u, tree_, std::move(vp));
      return {iterator(&tree_, new_v), true};
    }
    // otherwise, a new child will be needed below the 'before_u' node.
    // find an unbalanced sub-tree.
    auto [imbal_u, imbalanced] = find_imbalance_impl(before_u, 1);
    if (imbalanced) {
      // must re-balance and insert at sub-tree from imbal_u:
      key_type tmp_k = helper_type::value_to_key(vp);
      rebalance_and_insert_subtree(imbal_u, std::move(vp));
      return {iterator(&tree_, find_lower_bound(tmp_k, imbal_u)), true};
    }
    auto [new_v, new_e, new_added] = add_child(before_u, tree_, std::move(vp));
    swap(tree_[new_v], tree_[before_u]);
    return {iterator(&tree_, before_u), true};
  }

  vertex_type insert_in_range_impl(iterator tvi, iterator tvi_end,
                                   vertex_property&& vp) {
    // lets try to find a sweet spot (balance-wise) to place the new node:
    vertex_type tmp_u = tvi.base();
    for (; tvi != tvi_end; ++tvi) {
      if (out_degree(tvi.base(), tree_) == 1) {
        // then, this is definitely a good spot:
        auto [new_v, new_e, new_added] =
            add_child(tvi.base(), tree_, std::move(vp));
        return new_v;
      }
    }
    // if we reach this point, then, try to go up to some re-balanceable node.
    auto dummy_status = bagl::bst_traversal_status::on_right_branch;
    bagl::bst_detail::bst_move_up_to_next(tree_, tmp_u, dummy_status);
    if (tmp_u != root_) {
      auto [imbal_u, imbalanced] = find_imbalance_impl(tmp_u, 1);
      if (imbalanced) {
        // must re-balance and insert at sub-tree from imbal_u:
        key_type tmp_k = helper_type::value_to_key(vp);
        rebalance_and_insert_subtree(imbal_u, std::move(vp));
        return find_lower_bound(tmp_k, imbal_u);
      }
      rebalance_and_insert_subtree(tmp_u, std::move(vp));
      return tmp_u;
    }
    return tree_indexer::null_vertex();
  }

  /*
   * The vertex before_u is the "lower-bound" (either the first equal element or the first greater element if value is
   * not
   * already in the tree).
   * The vertex after_u is the "upper-bound" (the first greater element).
   * The value is the property of the vertex to be inserted.
   */
  std::pair<iterator, bool> insert_impl(vertex_type before_u,
                                        vertex_type after_u,
                                        vertex_property vp) {  // by-value.
    using std::swap;

    if ((after_u == tree_indexer::null_vertex()) &&
        (before_u != tree_indexer::null_vertex())) {
      // this means that all the vertices remaining after before_u are all equal to the given value.
      vertex_type tmp_u = insert_in_range_impl(
          iterator(&tree_, before_u), iterator::end(&tree_), std::move(vp));
      if (tmp_u == tree_indexer::null_vertex()) {
        // this will cause it to be inserted at the end (see below).
        before_u = tmp_u;
      } else {
        return {iterator(&tree_, tmp_u), true};
      }
    }

    if ((before_u == tree_indexer::null_vertex()) &&
        (after_u == tree_indexer::null_vertex())) {
      // this means that the vertex must be appended to the end. (value is greater than all elements in the tree)
      before_u = bagl::bst_detail::bst_go_down_right(tree_, root_);
      return insert_after_terminal_impl(before_u, std::move(vp));
    }

    // then, the special case where the equal range is empty.
    if (before_u == after_u) {
      if (out_degree(after_u, tree_) == 0) {
        // this means we could add a left-child to this node, see if it causes an imbalance:
        auto [imbal_u, imbalanced] = find_imbalance_impl(after_u, 1);
        if (imbalanced) {
          // must re-balance and insert at sub-tree from imbal_u:
          key_type tmp_k = helper_type::value_to_key(vp);
          rebalance_and_insert_subtree(imbal_u, std::move(vp));
          return {iterator(&tree_, find_lower_bound(tmp_k, imbal_u)), true};
        }
        auto [new_v, new_e, new_added] =
            add_child(after_u, tree_, std::move(vp));
        return {iterator(&tree_, new_v), true};
      }

      if (out_degree(after_u, tree_) == 1) {
        // this means that the left-child must be a leaf (otherwise it would be imbalanced), so, we can rotate the tree.
        auto [new_v, new_e, new_added] =
            add_child(after_u, tree_, std::move(vp));
        swap(tree_[after_u], tree_[new_v]);
        return {iterator(&tree_, *(children(after_u, tree_).begin())), true};
      }

      // else, after_u is full, lets try after the previous node (which must be a leaf or non-full node):
      before_u = bagl::bst_detail::bst_go_down_right(
          tree_, *(children(after_u, tree_).begin()));
      return insert_after_terminal_impl(before_u, std::move(vp));
    }

    // then, the general case: the equal range is somewhere in the middle of the set.
    vertex_type tmp_u = insert_in_range_impl(
        iterator(&tree_, before_u), iterator(&tree_, after_u), std::move(vp));
    if (tmp_u == tree_indexer::null_vertex()) {
      // this case is impossible.
      return {iterator::end(&tree_), false};
    }
    return {iterator(&tree_, tmp_u), true};
    // at this point, the vertex must have been inserted one way or another (worst-case: the whole tree got
    // re-constructed).
  }

  /*
   * Erase all vertices in the range [before, after)  (excl. after).
   * The vertex before_it is the "lower-bound" (either the first equal element or the first greater element if value is
   * not
   * already in the tree).
   * The vertex after_it is the "upper-bound" (the first greater element).
   */
  iterator erase_impl(iterator before_it, iterator after_it) {
    iterator before_after_it = after_it;
    --before_after_it;
    if ((before_it == after_it) ||
        (root_ == bagl::graph_traits<tree_indexer>::null_vertex())) {
      return before_it;
    }

    key_type tmp_after_key = {};
    bool after_at_end = false;
    if (after_it != iterator::end(&tree_)) {
      tmp_after_key = helper_type::value_to_key(tree_[after_it.base()]);
    } else {
      after_at_end = true;
    }

    // find the biggest possible sub-tree that contains all elements from the range.
    vertex_type mid_root = root_;
    while (true) {
      // if root is less than lower-bound than go right.
      if (compare_(tree_[mid_root], tree_[before_it.base()])) {
        mid_root = get_right_child(mid_root);
        continue;
      }
      // if root is greater than upper-bound than go left.
      if (compare_(tree_[before_after_it.base()], tree_[mid_root])) {
        mid_root = get_left_child(mid_root);
        continue;
      }
      // if the root is the highest root that is neither less than 'before' nor greater than 'after'
      break;
    }

    // collect and remove.
    std::vector<vertex_property> collected_nodes;
    iterator it_near(&tree_,
                     bagl::bst_detail::bst_go_down_left(tree_, mid_root));
    iterator it_far(&tree_,
                    bagl::bst_detail::bst_go_down_right(tree_, mid_root));

    while (true) {
      if (it_near == before_it) {
        if (it_far.base() == before_after_it.base()) {
          break;
        }
        it_near = after_it;
        continue;
      }
      collected_nodes.push_back(std::move(tree_[it_near.base()]));
      if (it_near == it_far) {
        break;
      }
      ++it_near;
    }

    if (has_right_child(mid_root)) {
      remove_branch(get_right_child(mid_root), tree_);
    }
    remove_branch(get_left_child(mid_root), tree_);

    if ((root_ == mid_root) && (collected_nodes.empty())) {
      remove_branch(root_, tree_);
      root_ = bagl::graph_traits<tree_indexer>::null_vertex();
      return iterator::end(&tree_);
    }

    if (collected_nodes.empty()) {
      vertex_type tmp_parent =
          source(*in_edges(mid_root, tree_).begin(), tree_);
      remove_branch(mid_root, tree_);
      mid_root = tmp_parent;
    } else {
      // re-construct:
      construct_node(mid_root, collected_nodes.begin(), collected_nodes.end());
    }

    while (true) {
      std::pair<vertex_type, bool> imbal_u = find_imbalance_impl(mid_root, 0);
      if (!imbal_u.second) {
        break;
      }
      rebalance_subtree(imbal_u.first);
      mid_root = imbal_u.first;
    }

    if (after_at_end) {
      return iterator::end(&tree_);
    }
    return iterator(&tree_, find_lower_bound(tmp_after_key, mid_root));
  }

 public:
  /**
   * Creates a AVL-tree with no elements.
   */
  explicit avl_tree_impl(const allocator_type& /*unused*/ = allocator_type())
      : tree_(),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        compare_(Compare()) {}

  /**
   * Creates a AVL-tree with no elements.
   */
  explicit avl_tree_impl(const Compare& comp,
                         const allocator_type& /*unused*/ = allocator_type())
      : tree_(),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        compare_(comp) {}

  /**
   * Builds a AVL-tree from an iterator range.
   * \param first An input iterator (start of range).
   * \param last An input iterator (end of range).
   * \param comp A comparison functor.
   */
  template <typename InputIterator>
  avl_tree_impl(InputIterator first, InputIterator last,
                const Compare& comp = Compare(),
                const allocator_type& /*unused*/ = allocator_type())
      : tree_(),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        compare_(comp) {
    root_ = create_root(tree_);
    std::vector<value_type> v(first, last);
    construct_node(root_, v.begin(), v.end());
  }

  /**
   * Builds a AVL-tree from an initializer_list.
   * \param vlist An std::initializer_list.
   * \param comp A comparison functor.
   */
  avl_tree_impl(std::initializer_list<value_type> vlist,
                const Compare& comp = Compare(),
                const allocator_type& /*unused*/ = allocator_type())
      : tree_(),
        root_(bagl::graph_traits<tree_indexer>::null_vertex()),
        compare_(comp) {
    root_ = create_root(tree_);
    std::vector<value_type> v(vlist);
    construct_node(root_, v.begin(), v.end());
  }

  /**
   * Checks if the AVL-tree is empty.
   * \return True if the AVL-tree is empty.
   */
  bool empty() const noexcept { return (num_vertices(tree_) == 0); }
  /**
   * Returns the size of the AVL-tree (the number of vertices it contains).
   * \return The size of the AVL-tree (the number of vertices it contains).
   */
  size_type size() const noexcept { return num_vertices(tree_); }
  /**
   * Returns the maximum size of the AVL-tree.
   * \return The maximum size of the AVL-tree.
   */
  size_type max_size() const noexcept {
    return std::numeric_limits<size_type>::max();
  }

  /**
   * Standard swap function.
   */
  void swap(self& rhs) {
    using std::swap;
    swap(this->tree_, rhs.tree_);
    swap(this->root_, rhs.root_);
    swap(this->compare_, rhs.compare_);
  }

  /**
   * Standard swap function.
   */
  friend void swap(self& lhs, self& rhs) { lhs.swap(rhs); }

  /**
   * Returns the depth of the tree.
   * \return The depth of the tree.
   */
  size_type depth() const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return 0;
    }
    return get_minmax_depth(root_).second;
  }

  /** Returns the comparison object with which the tree was constructed.  */
  key_compare key_comp() const { return compare_; }
  /** Returns the comparison object with which the tree was constructed.  */
  value_compare value_comp() const { return compare_; }

  /** Returns the allocator object with which the tree was constructed.  */
  allocator_type get_allocator() const { return allocator_type(); }

  /**
   * Returns a read-only (constant) iterator that points to the first element in the set.
   * Iteration is done in ascending order according to the keys.
   */
  const_iterator begin() const { return const_iterator::begin(&tree_); }

  /**
   * Returns a read-only (constant) iterator that points one past the last element in the set.
   * Iteration is done in ascending order according to the keys.
   */
  const_iterator end() const { return const_iterator::end(&tree_); }

  /**
   * Returns a read-only (constant) iterator that points to the first element in the set.
   * Iteration is done in ascending order according to the keys.
   */
  const_iterator cbegin() const { return begin(); }
  /**
   * Returns a read-only (constant) iterator that points one past the last element in the set.
   * Iteration is done in ascending order according to the keys.
   */
  const_iterator cend() const { return end(); }

  /**
   * Returns a read-only (constant) reverse iterator that points to the last element in the set.
   * Iteration is done in descending order according to the keys.
   */
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(const_iterator::end(&tree_));
  }

  /**
   * Returns a read-only (constant) reverse iterator that points to the last element in the set.
   * Iteration is done in descending order according to the keys.
   */
  const_reverse_iterator rend() const {
    return const_reverse_iterator(const_iterator::begin(&tree_));
  }

  /**
   * Returns a read-only (constant) reverse iterator that points to the last element in the set.
   * Iteration is done in descending order according to the keys.
   */
  const_reverse_iterator crbegin() const { return rbegin(); }
  /**
   * Returns a read-only (constant) reverse iterator that points to the last element in the set.
   * Iteration is done in descending order according to the keys.
   */
  const_reverse_iterator crend() const { return rend(); }

  /**
   * Tries to locate a key in a set.
   * \param k Key to be located.
   * \return Iterator pointing to sought-after element, or end() if not found.
   */
  const_iterator find(const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return const_iterator::end(&tree_);
    }
    vertex_type u = find_lower_bound(k, root_);
    if (!compare_(k, tree_[u])) {
      return const_iterator(&tree_, u);
    }
    return const_iterator::end(&tree_);
  }

  /**
   * Finds the beginning of a subsequence matching given key.
   * \param k Key to be located.
   * \return Iterator pointing to first element equal to or greater than key, or end().
   */
  const_iterator lower_bound(const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return const_iterator::end(&tree_);
    }
    vertex_type u = find_lower_bound(k, root_);
    if (u == tree_indexer::null_vertex()) {
      return end();
    }
    return const_iterator(&tree_, u);
  }

  /**
   * Finds the end of a subsequence matching given key.
   * \param k Key to be located.
   * \return Iterator pointing to the first element greater than key, or end().
   */
  const_iterator upper_bound(const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return const_iterator::end(&tree_);
    }
    vertex_type u = find_upper_bound(k, root_);
    if (u == tree_indexer::null_vertex()) {
      return end();
    }
    return const_iterator(&tree_, u);
  }

  /**
   * Finds a subsequence matching given key.
   * \param k Key to be located.
   * \return Pair of iterators that possibly points to the subsequence matching given key.
   */
  std::pair<const_iterator, const_iterator> equal_range(
      const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return {const_iterator::end(&tree_), const_iterator::end(&tree_)};
    }
    return {lower_bound(k), upper_bound(k)};
  }

  /**
   * Finds the number of keys.
   * \param k Key to be located.
   * \return Number of elements with specified key.
   */
  size_type count(const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      return 0;
    }
    auto [first, last] = equal_range(k);
    return std::distance(first, last);
  }

  /**
   * Attempts to insert an element into the set.
   * \param vp Element to be inserted.
   * \return A pair, of which the first element is an iterator that points to the possibly inserted element,
   *         and the second is a bool that is true if the element was actually inserted.
   */
  std::pair<iterator, bool> insert(const value_type& vp) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      root_ = create_root(tree_, vp);
      return {begin(), true};
    }
    key_type k = helper_type::value_to_key(vp);
    return insert_impl(find_lower_bound(k, root_), find_upper_bound(k, root_),
                       vp);
  }

  /**
   * Attempts to insert an element into the set.
   * \param vp Element to be inserted.
   * \return A pair, of which the first element is an iterator that points to the possibly inserted element,
   *         and the second is a bool that is true if the element was actually inserted.
   */
  template <typename P>
  std::pair<iterator, bool> insert(P&& vp) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      root_ = create_root(tree_, std::forward<P>(vp));
      return {begin(), true};
    }
    key_type k = helper_type::value_to_key(vp);
    return insert_impl(find_lower_bound(k, root_), find_upper_bound(k, root_),
                       value_type(std::forward<P>(vp)));
  }

  /**
   * Attempts to insert an element into the set.
   * \param pos An iterator that serves as a hint as to where the element should be inserted.
   * \param vp Element to be inserted.
   * \return An iterator that points to the element with key of x (may or may not be the element passed in).
   */
  iterator insert(const_iterator pos, const value_type& vp) {
    if (pos == end()) {
      return insert(vp).first;
    }
    if (compare_(*pos, vp)) {  // if position is less than value.
      const_iterator p2 = pos;
      ++p2;
      if (!compare_(*p2, vp)) {  // if p2 if not less than value.
        // find first iterator that is greater than value:
        const_iterator p3 = p2;
        while ((p3 != end()) && (!compare_(vp, *p3))) {
          ++p3;
        }
        return insert_impl(p2.base(), p3.base(), vp).first;
      }
      // else, it means that the hint is wrong:
      return insert(vp).first;
    }
    // else, position is not less than value.
    if (pos == begin()) {
      return insert_impl(pos.base(), pos.base(), vp).first;
    }
    const_iterator p4 = pos;
    --p4;
    if (!compare_(*p4, vp)) {
      // hint is wrong:
      return insert(vp).first;
    }
    // else, p4 is less than value. Find first iterator that is greater than value:
    const_iterator p5 = pos;
    while ((p5 != end()) && (!compare_(vp, *p5))) {
      ++p5;
    }
    return insert_impl(pos.base(), p5.base(), vp).first;
  }

  /**
   * Attempts to insert an element into the set.
   * \param pos An iterator that serves as a hint as to where the element should be inserted.
   * \param vp Element to be inserted.
   * \return An iterator that points to the element with key of x (may or may not be the element passed in).
   */
  template <typename P>
  iterator insert(const_iterator pos, P&& vp) {
    if (pos == end()) {
      return insert(std::forward<P>(vp)).first;
    }
    if (compare_(*pos, vp)) {  // if position is less than value.
      const_iterator p2 = pos;
      ++p2;
      if (!compare_(*p2, vp)) {  // if p2 if not less than value.
        // find first iterator that is greater than value:
        const_iterator p3 = p2;
        while ((p3 != end()) && (!compare_(vp, *p3))) {
          ++p3;
        }
        return insert_impl(p2.base(), p3.base(),
                           value_type(std::forward<P>(vp)))
            .first;
      }
      // else, it means that the hint is wrong:
      return insert(std::forward<P>(vp)).first;
    }
    // else, position is not less than value.
    if (pos == begin()) {
      return insert_impl(pos.base(), pos.base(),
                         value_type(std::forward<P>(vp)))
          .first;
    }
    const_iterator p4 = pos;
    --p4;
    if (!compare_(*p4, vp)) {
      // hint is wrong:
      return insert(std::forward<P>(vp)).first;
    }
    // else, p4 is less than value. Find first iterator that is greater than value:
    const_iterator p5 = pos;
    while ((p5 != end()) && (!compare_(vp, *p5))) {
      ++p5;
    }
    return insert_impl(pos.base(), p5.base(), value_type(std::forward<P>(vp)))
        .first;
  }

  /**
   * A template function that attempts to insert a range of elements.
   * \tparam InputIterator An input-iterator type that can be dereferenced to a value-type rvalue (or rvalue-ref
   * (C++11)).
   * \param first Iterator pointing to the start of the range to be inserted.
   * \param last Iterator pointing to the end of the range.
   */
  template <typename InputIterator>
  void insert(InputIterator first, InputIterator last) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      root_ = create_root(tree_);
      std::vector<value_type> v(first, last);
      construct_node(root_, v.begin(), v.end());
      return;
    }
    for (; first != last; ++first) {
      insert(*first);
    }
  }

  /**
   * Attempts to insert a list of elements into the set.
   * \param aList A std::initializer_list of elements to be inserted.
   */
  void insert(std::initializer_list<value_type> vlist) {
    insert(vlist.begin(), vlist.end());
  }

  /**
   * Attempts to emplace-create an element into the set.
   * \param args The constructor arguments required to create the element.
   * \return A pair, of which the first element is an iterator that points to the possibly inserted element,
   *         and the second is a bool that is true if the element was actually inserted.
   * \note This function does not really create the element emplace since it must first form the element to find
   *       its appropriate position in the set.
   */
  template <typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    value_type tmp_val(std::forward<Args>(args)...);
    return insert(std::move(tmp_val));
  }

  /**
   * Attempts to emplace-create an element into the set.
   * \param pos An iterator that serves as a hint as to where the element should be inserted.
   * \param args The constructor arguments required to create the element.
   * \return A pair, of which the first element is an iterator that points to the possibly inserted element,
   *         and the second is a bool that is true if the element was actually inserted.
   * \note This function does not really create the element emplace since it must first form the element to find
   *       its appropriate position in the set.
   */
  template <typename... Args>
  iterator emplace_hint(const_iterator pos, Args&&... args) {
    value_type tmp_val(std::forward<Args>(args)...);
    return insert(pos, std::move(tmp_val));
  }

  /**
   * Erases a [first,last) range of elements from a set.
   * \param first Iterator pointing to the start of the range to be erased.
   * \param last Iterator pointing to the end of the range to be erased.
   * \return New iterator to last.
   */
  iterator erase(const_iterator first, const_iterator last) {
    if (first == last) {
      return last;
    }
    return erase_impl(first, last);
  }

  /**
   * Erases an element from a set.
   * \param pos An iterator pointing to the element to be erased.
   * \return An iterator pointing to the element immediately following position
   *         prior to the element being erased. If no such element exists, end() is returned.
   */
  iterator erase(const_iterator pos) {
    const_iterator p2 = pos;
    ++p2;
    return erase(pos, p2);
  }

  /**
   * Erases elements according to the provided key.
   * \param k Key of element to be erased.
   * \return The number of elements erased.
   */
  size_type erase(const key_type& k) {
    std::pair<const_iterator, const_iterator> er = equal_range(k);
    size_type result = std::distance(er.first, er.second);
    if (result == 0) {
      return 0;
    }
    erase(er.first, er.second);
    return result;
  }

  /**
 * Erases all elements in a set.
 */
  void clear() noexcept {
    if (num_vertices(tree_) != 0) {
      remove_branch(root_, tree_);
      root_ = bagl::graph_traits<tree_indexer>::null_vertex();
    }
  }

  /**
   * Subscript ( [] ) access to map data. Allows for easy lookup with the subscript ( [] ) operator.
   * Returns data associated with the key specified in subscript. If the key does not exist, a pair
   * with that key is created using default values, which is then returned.
   * \param k The key for which data should be retrieved.
   * \return A reference to the data of the (key,data) pair.
   * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
   */
  mapped_type& operator[](const key_type& k) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      iterator it =
          insert(helper_type::keymap_to_value(k, mapped_type())).first;
      return helper_type::value_to_mapped(tree_[it.base()]);
    }
    vertex_type u = find_lower_bound(k, root_);
    if (!compare_(k, tree_[u])) {
      return helper_type::value_to_mapped(tree_[u]);
    }  // insert the key:
    iterator it = insert(const_iterator(&tree_, u),
                         helper_type::keymap_to_value(k, mapped_type()));
    return helper_type::value_to_mapped(tree_[it.base()]);
  }

  /**
   * Subscript ( [] ) access to map data. Allows for easy lookup with the subscript ( [] ) operator.
   * Returns data associated with the key specified in subscript. If the key does not exist, a pair
   * with that key is created using default values, which is then returned.
   * \param k The key for which data should be retrieved.
   * \return A reference to the data of the (key,data) pair.
   * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
   */
  mapped_type& operator[](key_type&& k) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      iterator it =
          insert(helper_type::keymap_to_value(std::move(k), mapped_type()))
              .first;
      return helper_type::value_to_mapped(tree_[it.base()]);
    }
    vertex_type u = find_lower_bound(k, root_);
    if (!compare_(k, tree_[u])) {
      return helper_type::value_to_mapped(tree_[u]);
    }  // insert the key:
    iterator it =
        insert(const_iterator(&tree_, u),
               helper_type::keymap_to_value(std::move(k), mapped_type()));
    return helper_type::value_to_mapped(tree_[it.base()]);
  }

  /**
   * Access to map data.
   * \param k The key for which data should be retrieved.
   * \return A reference to the data whose key is equivalent to k, if such a data is present in the map.
   * \throw std::out_of_range If no such data is present.
   * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
   */
  mapped_type& at(const key_type& k) {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      throw std::out_of_range("The key does not match an element of the map!");
    }
    vertex_type u = find_lower_bound(k, root_);
    if (!compare_(k, tree_[u])) {
      return helper_type::value_to_mapped(tree_[u]);
    }
    throw std::out_of_range("The key does not match an element of the map!");
  }

  /**
   * Access to map const data.
   * \param k The key for which data should be retrieved.
   * \return A const reference to the data whose key is equivalent to k, if such a data is present in the map.
   * \throw std::out_of_range If no such data is present.
   * \note If used on a multimap or multiset, this function will use the first match (lower-bound).
   */
  const mapped_type& at(const key_type& k) const {
    if (root_ == bagl::graph_traits<tree_indexer>::null_vertex()) {
      throw std::out_of_range("The key does not match an element of the map!");
    }
    vertex_type u = find_lower_bound(k, root_);
    if (!compare_(k, tree_[u])) {
      return helper_type::value_to_mapped(tree_[u]);
    }
    throw std::out_of_range("The key does not match an element of the map!");
  }

  /* for private use only (useful to keep a map / set synchronized to another dynamic structure. */
  struct mutation_visitor {
    self* parent_;

    explicit mutation_visitor(self* p) : parent_(p) {}

    void remove_vertex(vertex_type v, tree_indexer& /*unused*/) const {
      parent_->erase(v);
    }

    void add_vertex(tree_indexer& /*unused*/, const vertex_property& vp) const {
      parent_->insert(vp);
    }

    void add_vertex(tree_indexer& /*unused*/, vertex_property&& vp) const {
      parent_->insert(std::move(vp));
    }
  };
};

}  // namespace ReaK::graph

#endif  // REAK_PLANNING_GRAPH_ALG_AVL_TREE_DETAIL_H_
