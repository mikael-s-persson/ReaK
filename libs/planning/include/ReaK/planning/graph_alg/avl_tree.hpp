/**
 * \file avl_tree.hpp
 *
 * This library implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
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

#ifndef REAK_AVL_TREE_HPP
#define REAK_AVL_TREE_HPP

#include <ReaK/core/base/defs.hpp>

#include "avl_tree_detail.hpp"

// BGL-Extra includes:
#include <boost/graph/bfl_d_ary_tree.hpp>
#include <boost/graph/vebl_d_ary_tree.hpp>

namespace ReaK::graph {

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::set class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam T Type of the elements. Each element in a set container is also uniquely identified by this
 *           value (each value is itself also the element's key). Aliased as member types avlbfl_set::key_type and
 * avlbfl_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool.
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T>>
struct avlbfl_set : public avl_tree_impl<boost::bfl_d_ary_tree<2, T>, Compare,
                                         detail::avl_set_style> {
  using base_impl_type = avl_tree_impl<boost::bfl_d_ary_tree<2, T>, Compare,
                                       detail::avl_set_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlbfl_set(const allocator_type&) : base_impl_type() {}

  avlbfl_set() : avlbfl_set(allocator_type()) {}

  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlbfl_set(const Compare& comp,
                      const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree set from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbfl_set(InputIterator aFirst, InputIterator aLast,
             const Compare& comp = Compare(),
             const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree set from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlbfl_set(std::initializer_list<value_type> aList,
             const Compare& comp = Compare(),
             const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::multiset class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam T Type of the elements. Each element in a set container is also uniquely identified by this
 *           value (each value is itself also the element's key). Aliased as member types avlbfl_multiset::key_type and
 * avlbfl_multiset::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool.
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T>>
struct avlbfl_multiset
    : public avl_tree_impl<boost::bfl_d_ary_tree<2, T>, Compare,
                           detail::avl_multiset_style> {
  using base_impl_type = avl_tree_impl<boost::bfl_d_ary_tree<2, T>, Compare,
                                       detail::avl_multiset_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlbfl_multiset(const allocator_type&) : base_impl_type() {}

  avlbfl_multiset() : avlbfl_multiset(allocator_type()) {}

  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlbfl_multiset(const Compare& comp,
                           const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree multiset from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbfl_multiset(InputIterator aFirst, InputIterator aLast,
                  const Compare& comp = Compare(),
                  const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree multiset from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlbfl_multiset(std::initializer_list<value_type> aList,
                  const Compare& comp = Compare(),
                  const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::map class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam Key Type of the keys. Each element in a map is uniquely identified by its key value. Aliased as
 *             member type avlbfl_map::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as
 *           member type avlbfl_map::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values,
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key>>
struct avlbfl_map
    : public avl_tree_impl<boost::bfl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                           detail::avl_map_style> {
  using base_impl_type =
      avl_tree_impl<boost::bfl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                    detail::avl_map_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlbfl_map(const allocator_type&) : base_impl_type() {}

  avlbfl_map() : avlbfl_map(allocator_type()) {}

  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlbfl_map(const Compare& comp,
                      const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree map from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbfl_map(InputIterator aFirst, InputIterator aLast,
             const Compare& comp = Compare(),
             const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree map from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlbfl_map(std::initializer_list<value_type> aList,
             const Compare& comp = Compare(),
             const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::multimap class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam Key Type of the keys. Each element in a map is uniquely identified by its key value. Aliased as
 *             member type avlbfl_multimap::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as
 *           member type avlbfl_multimap::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values,
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key>>
struct avlbfl_multimap
    : public avl_tree_impl<boost::bfl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                           detail::avl_multimap_style> {
  using base_impl_type =
      avl_tree_impl<boost::bfl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                    detail::avl_multimap_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlbfl_multimap(const allocator_type&) : base_impl_type() {}

  avlbfl_multimap() : avlbfl_multimap(allocator_type()) {}

  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlbfl_multimap(const Compare& comp,
                           const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree multimap from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbfl_multimap(InputIterator aFirst, InputIterator aLast,
                  const Compare& comp = Compare(),
                  const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree multimap from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlbfl_multimap(std::initializer_list<value_type> aList,
                  const Compare& comp = Compare(),
                  const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::set class template.
 * The AVL tree is stored in a cache-oblivious B-tree layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam T Type of the elements. Each element in a set container is also uniquely identified by this
 *           value (each value is itself also the element's key). Aliased as member types avlvebl_set::key_type and
 * avlvebl_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool.
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T>>
struct avlvebl_set : public avl_tree_impl<boost::vebl_d_ary_tree<2, T>, Compare,
                                          detail::avl_set_style> {
  using base_impl_type = avl_tree_impl<boost::vebl_d_ary_tree<2, T>, Compare,
                                       detail::avl_set_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlvebl_set(const allocator_type&) : base_impl_type() {}

  avlvebl_set() : avlvebl_set(allocator_type()) {}

  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlvebl_set(const Compare& comp,
                       const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree set from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlvebl_set(InputIterator aFirst, InputIterator aLast,
              const Compare& comp = Compare(),
              const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree set from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlvebl_set(std::initializer_list<value_type> aList,
              const Compare& comp = Compare(),
              const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::multiset class template.
 * The AVL tree is stored in a cache-oblivious B-tree layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam T Type of the elements. Each element in a set container is also uniquely identified by this
 *           value (each value is itself also the element's key). Aliased as member types avlvebl_set::key_type and
 * avlvebl_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool.
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T>>
struct avlvebl_multiset
    : public avl_tree_impl<boost::vebl_d_ary_tree<2, T>, Compare,
                           detail::avl_multiset_style> {
  using base_impl_type = avl_tree_impl<boost::vebl_d_ary_tree<2, T>, Compare,
                                       detail::avl_multiset_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlvebl_multiset(const allocator_type&) : base_impl_type() {}

  avlvebl_multiset() : avlvebl_multiset(allocator_type()) {}

  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlvebl_multiset(const Compare& comp,
                            const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree multiset from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlvebl_multiset(InputIterator aFirst, InputIterator aLast,
                   const Compare& comp = Compare(),
                   const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree multiset from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlvebl_multiset(std::initializer_list<value_type> aList,
                   const Compare& comp = Compare(),
                   const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::map class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam Key Type of the keys. Each element in a map is uniquely identified by its key value. Aliased as
 *             member type avlvebl_map::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as
 *           member type avlvebl_map::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values,
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key>>
struct avlvebl_map
    : public avl_tree_impl<boost::vebl_d_ary_tree<2, std::pair<Key, T>>,
                           Compare, detail::avl_map_style> {
  using base_impl_type =
      avl_tree_impl<boost::vebl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                    detail::avl_map_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlvebl_map(const allocator_type&) : base_impl_type() {}

  avlvebl_map() : avlvebl_map(allocator_type()) {}

  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlvebl_map(const Compare& comp,
                       const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree map from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlvebl_map(InputIterator aFirst, InputIterator aLast,
              const Compare& comp = Compare(),
              const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree map from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlvebl_map(std::initializer_list<value_type> aList,
              const Compare& comp = Compare(),
              const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering,
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::multimap class template.
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam Key Type of the keys. Each element in a map is uniquely identified by its key value. Aliased as
 *             member type avlvebl_multimap::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as
 *           member type avlvebl_multimap::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values,
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key>>
struct avlvebl_multimap
    : public avl_tree_impl<boost::vebl_d_ary_tree<2, std::pair<Key, T>>,
                           Compare, detail::avl_multimap_style> {
  using base_impl_type =
      avl_tree_impl<boost::vebl_d_ary_tree<2, std::pair<Key, T>>, Compare,
                    detail::avl_multimap_style>;
  using value_type = typename base_impl_type::value_type;
  using allocator_type = typename base_impl_type::allocator_type;

  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlvebl_multimap(const allocator_type&) : base_impl_type() {}

  avlvebl_multimap() : avlvebl_multimap(allocator_type()) {}

  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlvebl_multimap(const Compare& comp,
                            const allocator_type& = allocator_type())
      : base_impl_type(comp) {}

  /**
    * Builds a AVL-tree multimap from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlvebl_multimap(InputIterator aFirst, InputIterator aLast,
                   const Compare& comp = Compare(),
                   const allocator_type& = allocator_type())
      : base_impl_type(aFirst, aLast, comp) {}

  /**
    * Builds a AVL-tree multimap from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */
  avlvebl_multimap(std::initializer_list<value_type> aList,
                   const Compare& comp = Compare(),
                   const allocator_type& = allocator_type())
      : base_impl_type(aList, comp) {}
};

}  // namespace ReaK::graph

#endif
