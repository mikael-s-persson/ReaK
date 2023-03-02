/**
 * \file quick_sort.hpp
 *
 * This library provides a generic quick-sort functions.
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

#ifndef REAK_QUICK_SORT_HPP
#define REAK_QUICK_SORT_HPP

#include <ReaK/core/base/defs.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iterator>

namespace ReaK::sorting {

/**
 * This pivot chooser functor chooses a random pivot from the range.
 */
struct random_pivot {
  template <typename RandomAccessIter, typename Compare>
  void operator()(RandomAccessIter first, RandomAccessIter last,
                  Compare) const {
    std::iter_swap(std::prev(last),
                   std::next(first, rand() % std::distance(first, last)));
  }
};

/**
 * This pivot chooser functor chooses a pivot from a range by picking the
 * median of three elements: the first, last and middle element.
 */
struct median_of_3_pivots {
  template <typename RandomAccessIter, typename Compare>
  void operator()(RandomAccessIter first, RandomAccessIter last,
                  Compare comp) const {
    RandomAccessIter before_last = std::prev(last);
    RandomAccessIter middle = std::next(first, std::distance(first, last) / 2);
    if (comp(*first, *before_last) && comp(*before_last, *middle)) {
      return;
    } else if (comp(*before_last, *first) && comp(*first, *middle)) {
      return std::iter_swap(before_last, first);
    }
    std::iter_swap(before_last, middle);
  }
};

/**
 * This pivot chooser functor chooses a pivot from a range by picking the
 * median of three random elements.
 */
struct median_of_3_random_pivots {
  template <typename RandomAccessIter, typename Compare>
  void operator()(RandomAccessIter first, RandomAccessIter last,
                  Compare comp) const {
    std::size_t dist = std::distance(first, last);
    RandomAccessIter piv1 = first + rand() % dist;
    RandomAccessIter piv2 = first + rand() % dist;
    RandomAccessIter piv3 = first + rand() % dist;
    if (comp(*piv1, *piv2) && comp(*piv2, *piv3))
      return std::iter_swap(std::prev(last), piv2);
    else if (comp(*piv2, *piv1) && comp(*piv1, *piv3))
      return std::iter_swap(std::prev(last), piv1);
    std::iter_swap(std::prev(last), piv3);
  }
};

/**
 * This pivot chooser functor chooses the first element from the range as the pivot.
 */
struct first_pivot {
  template <typename RandomAccessIter, typename Compare>
  void operator()(RandomAccessIter first, RandomAccessIter last,
                  Compare) const {
    std::iter_swap(std::prev(last), first);
  }
};

namespace detail {

template <typename RandomAccessIter, typename Compare, typename PivotChooser>
void quick_sort_ra_impl(RandomAccessIter first, RandomAccessIter last,
                        Compare comp, PivotChooser choose_pivot) {
  while (first != last) {
    if (last - first < 50)
      return insertion_sort(first, last, comp);
    choose_pivot(first, last, comp);
    RandomAccessIter before_last = last - 1;
    RandomAccessIter pivot = std::partition(
        first, before_last,
        [&](const auto& x) -> bool { return comp(x, *before_last); });
    std::iter_swap(pivot, before_last);
    if (pivot - first < last - pivot) {
      quick_sort_ra_impl(first, pivot, comp, choose_pivot);
      first = pivot + 1;
    } else {
      quick_sort_ra_impl(pivot + 1, last, comp, choose_pivot);
      last = pivot;
    }
  }
}

template <typename BidirIter, typename Compare, typename PivotChooser>
void quick_sort_nonra_impl(BidirIter first, BidirIter last, Compare comp,
                           PivotChooser choose_pivot) {
  while (first == last) {
    choose_pivot(first, last, comp);
    BidirIter before_last = std::prev(last);
    BidirIter pivot = std::partition(
        first, before_last,
        [&](const auto& x) -> bool { return comp(x, *before_last); });
    std::iter_swap(pivot, before_last);
    std::size_t dist1 = std::distance(first, pivot);
    std::size_t dist2 = std::distance(pivot, before_last);
    if (dist1 < dist2) {
      quick_sort_nonra_impl(first, pivot, comp, choose_pivot);
      first = std::next(pivot);
      dist1 = dist2;
    } else {
      quick_sort_nonra_impl(std::next(pivot), last, comp, choose_pivot);
      last = pivot;
    }
    if (dist1 < 2) {
      return;
    }
    if (dist1 < 50) {
      return insertion_sort(first, last, comp);
    }
  }
}

template <typename RandomAccessIter, typename Compare, typename PivotChooser>
inline void quick_sort_tagged(RandomAccessIter first, RandomAccessIter last,
                              Compare comp, std::random_access_iterator_tag t,
                              PivotChooser pivot_chooser) {
  quick_sort_ra_impl(first, last, comp, pivot_chooser);
}

template <typename RandomAccessIter, typename Compare>
inline void quick_sort_tagged(RandomAccessIter first, RandomAccessIter last,
                              Compare comp, std::random_access_iterator_tag t) {
  quick_sort_ra_impl(first, last, comp, median_of_3_pivots());
}

template <typename BidirIter, typename Compare>
inline void quick_sort_tagged(BidirIter first, BidirIter last, Compare comp,
                              std::bidirectional_iterator_tag t) {
  quick_sort_nonra_impl(first, last, comp, first_pivot());
}

}  // namespace detail

/**
 * This function performs a quick sort on a given range of elements, and with the
 * given comparison functor and the given pivot chooser (e.g., random_pivot,
 * median_of_3_pivots, median_of_3_random_pivots, or first_pivot).
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename Iter, typename Compare, typename PivotChooser>
inline void quick_sort(Iter first, Iter last, Compare comp,
                       PivotChooser choose_pivot) {
  detail::quick_sort_tagged(
      first, last, comp,
      typename std::iterator_traits<Iter>::iterator_category(), choose_pivot);
}

/**
 * This function performs a quick sort on a given range of elements, and with the
 * given comparison functor.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename Iter, typename Compare>
inline void quick_sort(Iter first, Iter last, Compare comp) {
  detail::quick_sort_tagged(
      first, last, comp,
      typename std::iterator_traits<Iter>::iterator_category());
}

/**
 * This function performs a quick sort on a given range of elements, and by using the
 * less-than operator on the elements of the range.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 */
template <typename RandomAccessIter>
inline void quick_sort(RandomAccessIter first, RandomAccessIter last) {
  quick_sort(first, last, std::less<>());
}

/**
 * This function performs a quick-select sort on a given range of elements, and with the
 * given comparison functor. A quick-select sort is just a succession of quick-select
 * algorithms until the entire range is sorted.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename RandomAccessIter, typename Compare>
void quickselect_sort(RandomAccessIter first, RandomAccessIter last,
                      Compare comp) {
  std::size_t dist = std::distance(first, last);
  if (dist < 2) {
    return;
  }
  RandomAccessIter pivot = std::next(first, dist / 2);
  std::nth_element(first, pivot, last, comp);
  quickselect_sort(first, pivot, comp);
  quickselect_sort(++pivot, last, comp);
}

/**
 * This function performs a quick sort on a given range of elements, and by using the
 * less-than operator on the elements of the range. A quick-select sort is just a
 * succession of quick-select algorithms until the entire range is sorted.
 * \tparam RandomAccessIter A random-access iterator type (input and output iterator).
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 */
template <typename RandomAccessIter>
inline void quickselect_sort(RandomAccessIter first, RandomAccessIter last) {
  quickselect_sort(first, last, std::less<>());
}

}  // namespace ReaK::sorting

#endif  // REAK_QUICK_SORT_HPP
