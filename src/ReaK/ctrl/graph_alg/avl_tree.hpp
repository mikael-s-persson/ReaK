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

#include "base/defs.hpp"

#include "graph_alg/avl_tree_detail.hpp"

#include "graph_alg/d_ary_bf_tree.hpp"
#include "graph_alg/d_ary_cob_tree.hpp"


namespace ReaK {

namespace graph {
  

/**
 * This class implements an AVL Tree that allows for O(logN) time look-ups in a strict weak ordering, 
 * with amortized O(logN) insertion-deletion. A AVL-tree is a self-balancing binary search tree which 
 * only requires a strict weak ordering (comparison function, analogous to less-than).
 * This class provides essentially the same interface as the standard std::set class template. 
 * The AVL tree is stored in a breadth-first layout in contiguous memory, and thus, the 
 * performance in enhanced in terms for cache-locality, at the amortized expense of copies during 
 * re-balancing and re-allocation. Note also, that insertion/deletions generally invalidate existing iterators.
 * \tparam T Type of the elements. Each element in a set container is also uniquely identified by this 
 *           value (each value is itself also the element's key). Aliased as member types avlbf_set::key_type and avlbf_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool. 
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall 
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T> >
struct avlbf_set : public avl_tree_impl< d_ary_bf_tree<T, 2>, Compare, detail::avl_set_style > {
  typedef avl_tree_impl< d_ary_bf_tree<T, 2>, Compare, detail::avl_set_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlbf_set(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlbf_set(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree set from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbf_set(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree set from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlbf_set(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aList, comp) { };
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
 *           value (each value is itself also the element's key). Aliased as member types avlbf_multiset::key_type and avlbf_multiset::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool. 
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall 
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T> >
struct avlbf_multiset : public avl_tree_impl< d_ary_bf_tree<T, 2>, Compare, detail::avl_multiset_style > {
  typedef avl_tree_impl< d_ary_bf_tree<T, 2>, Compare, detail::avl_multiset_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlbf_multiset(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlbf_multiset(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree multiset from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbf_multiset(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree multiset from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlbf_multiset(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aList, comp) { };
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
 *             member type avlbf_map::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as 
 *           member type avlbf_map::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The 
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values, 
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key> >
struct avlbf_map : public avl_tree_impl< d_ary_bf_tree< std::pair<Key,T>, 2>, Compare, detail::avl_map_style > {
  typedef avl_tree_impl< d_ary_bf_tree< std::pair<Key,T>, 2>, Compare, detail::avl_map_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlbf_map(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlbf_map(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree map from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbf_map(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree map from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlbf_map(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aList, comp) { };
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
 *             member type avlbf_multimap::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as 
 *           member type avlbf_multimap::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The 
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values, 
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key> >
struct avlbf_multimap : public avl_tree_impl< d_ary_bf_tree< std::pair<Key,T>, 2>, Compare, detail::avl_multimap_style > {
  typedef avl_tree_impl< d_ary_bf_tree< std::pair<Key,T>, 2>, Compare, detail::avl_multimap_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlbf_multimap(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlbf_multimap(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree multimap from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlbf_multimap(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                 base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree multimap from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlbf_multimap(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                 base_impl_type(aList, comp) { };
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
 *           value (each value is itself also the element's key). Aliased as member types avlcob_set::key_type and avlcob_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool. 
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall 
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T> >
struct avlcob_set : public avl_tree_impl< d_ary_cob_tree<T, 2>, Compare, detail::avl_set_style > {
  typedef avl_tree_impl< d_ary_cob_tree<T, 2>, Compare, detail::avl_set_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlcob_set(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree set with no elements.
    */
  explicit avlcob_set(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree set from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlcob_set(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree set from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlcob_set(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
            base_impl_type(aList, comp) { };
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
 *           value (each value is itself also the element's key). Aliased as member types avlcob_set::key_type and avlcob_set::value_type.
 * \tparam Compare A binary predicate that takes two arguments of the same type as the elements and returns a bool. 
 *                 The expression comp(a,b), where comp is an object of this type and a and b are key values, shall 
 *                 return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename T, typename Compare = std::less<T> >
struct avlcob_multiset : public avl_tree_impl< d_ary_cob_tree<T, 2>, Compare, detail::avl_multiset_style > {
  typedef avl_tree_impl< d_ary_cob_tree<T, 2>, Compare, detail::avl_multiset_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlcob_multiset(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree multiset with no elements.
    */
  explicit avlcob_multiset(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree multiset from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlcob_multiset(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  base_impl_type(aFirst, aLast, comp) { };
           
  /**
    * Builds a AVL-tree multiset from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlcob_multiset(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  base_impl_type(aList, comp) { };
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
 *             member type avlcob_map::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as 
 *           member type avlcob_map::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The 
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values, 
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key> >
struct avlcob_map : public avl_tree_impl< d_ary_cob_tree< std::pair<Key,T>, 2>, Compare, detail::avl_map_style > {
  typedef avl_tree_impl< d_ary_cob_tree< std::pair<Key,T>, 2>, Compare, detail::avl_map_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlcob_map(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree map with no elements.
    */
  explicit avlcob_map(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree map from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlcob_map(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
             base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree map from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlcob_map(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
             base_impl_type(aList, comp) { };
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
 *             member type avlcob_multimap::key_type.
 * \tparam T Type of the mapped value. Each element in a map stores some data as its mapped value. Aliased as 
 *           member type avlcob_multimap::mapped_type.
 * \tparam Compare A binary predicate that takes two element keys as arguments and returns a bool. The 
 *                 expression comp(a,b), where comp is an object of this type and a and b are key values, 
 *                 shall return true if a is considered to go before b in the strict weak ordering the function defines.
 */
template <typename Key, typename T, typename Compare = std::less<Key> >
struct avlcob_multimap : public avl_tree_impl< d_ary_cob_tree< std::pair<Key,T>, 2>, Compare, detail::avl_multimap_style > {
  typedef avl_tree_impl< d_ary_cob_tree< std::pair<Key,T>, 2>, Compare, detail::avl_multimap_style > base_impl_type;
  typedef typename base_impl_type::value_type value_type;
  typedef typename base_impl_type::allocator_type allocator_type;
  
  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlcob_multimap(const allocator_type& = allocator_type()) : base_impl_type() { };
  
  /**
    * Creates a AVL-tree multimap with no elements.
    */
  explicit avlcob_multimap(const Compare& comp, const allocator_type& = allocator_type()) : base_impl_type(comp) { };
  
  /**
    * Builds a AVL-tree multimap from an iterator range.
    * \param aFirst An input iterator (start of range).
    * \param aLast An input iterator (end of range).
    * \param comp A comparison functor.
    */
  template <typename InputIterator>
  avlcob_multimap(InputIterator aFirst, InputIterator aLast, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  base_impl_type(aFirst, aLast, comp) { };
        
  /**
    * Builds a AVL-tree multimap from an initializer_list.
    * \param aList An std::initializer_list.
    * \param comp A comparison functor.
    */  
  avlcob_multimap(std::initializer_list< value_type > aList, const Compare& comp = Compare(), const allocator_type& = allocator_type()) : 
                  base_impl_type(aList, comp) { };
};




};

};


#endif


















