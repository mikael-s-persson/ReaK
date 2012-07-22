/**
 * \file intro_sort.hpp
 * 
 * This library provides a generic intro-sort function.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2012
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

#ifndef REAK_INTRO_SORT_HPP
#define REAK_INTRO_SORT_HPP

#include "base/defs.hpp"

#include <algorithm>
#include <iterator>
#include <functional>

#include "quick_sort.hpp"

#include "insertion_sort.hpp"
#include "heap_sort.hpp"

#ifndef RK_ENABLE_CXX0X_FEATURES
#include <boost/bind.hpp>
#endif

namespace ReaK {
  
/** This is the namespace for all ReaK sorting algorithms implementations. */
namespace sorting {


namespace detail {

template <typename RandomAccessIter, typename Compare, typename PivotChooser>
void intro_sort_impl(RandomAccessIter first, RandomAccessIter last, Compare comp, PivotChooser choose_pivot) {
  int depth_count = 2 * std::log2(last - first);
  while(first != last) {
    std::size_t dist = last - first;
    if(dist < 2)
      return;
    if(dist < 50)
      return insertion_sort(first,last,comp);
    choose_pivot(first, last, comp);
    RandomAccessIter before_last = last - 1;
#ifdef RK_ENABLE_CXX0X_FEATURES
    RandomAccessIter pivot = std::partition(first, before_last, [&](decltype(*first) x) -> bool { return comp(x, *before_last); });
#else
    RandomAccessIter pivot = std::partition(first, before_last, boost::bind(comp, _1, *before_last));
#endif
    std::iter_swap(pivot, before_last);
    if( pivot - first < last - pivot ) {
      intro_sort_impl(first, pivot, comp, choose_pivot);
      first = pivot + 1;
    } else {
      intro_sort_impl(pivot + 1, last, comp, choose_pivot);
      last = pivot;
    };
    if(--depth_count < 0)
      return heap_sort(first,last,comp);
  };
};

}; // detail





/**
 * This function performs an intro sort on a given range of elements, and with the 
 * given comparison functor and the given pivot chooser (e.g., random_pivot, 
 * median_of_3_pivots, median_of_3_random_pivots, or first_pivot).
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename Iter, typename Compare, typename PivotChooser>
inline
void intro_sort(Iter first, Iter last, Compare comp, PivotChooser choose_pivot) {
  detail::intro_sort_impl(first,last,comp,choose_pivot);
};

/**
 * This function performs an intro sort on a given range of elements, and with the 
 * given comparison functor.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename Iter, typename Compare>
inline
void intro_sort(Iter first, Iter last, Compare comp) {
  detail::intro_sort_impl(first,last,comp,median_of_3_pivots());
};

/**
 * This function performs an intro sort on a given range of elements, and by using the
 * less-than operator on the elements of the range.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 */
template <typename RandomAccessIter>
inline
void intro_sort(RandomAccessIter first, RandomAccessIter last) {
  intro_sort(first, last, std::less< typename std::iterator_traits<RandomAccessIter>::value_type >());
};


};

};

#endif



