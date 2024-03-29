/**
 * \file merge_sort.h
 *
 * This library provides a generic merge-sort function.
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

#ifndef REAK_CORE_SORTING_MERGE_SORT_H_
#define REAK_CORE_SORTING_MERGE_SORT_H_

#include <algorithm>
#include <functional>
#include <iterator>

namespace ReaK::sorting {

/**
 * This function performs a merge sort on a given range of elements, and with the
 * given comparison functor.
 * \tparam BidirIter A bidirectional iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename BidirIter, typename Compare>
void merge_sort(BidirIter first, BidirIter last, Compare comp) {
  auto dist = std::distance(first, last) / 2;
  if (dist) {
    BidirIter middle = std::next(first, dist);
    merge_sort(first, middle, comp);
    merge_sort(middle, last, comp);
    std::inplace_merge(first, middle, last, comp);
  }
}

/**
 * This function performs a merge sort on a given range of elements, and by using the
 * less-than operator on the elements of the range.
 * \tparam BidirIter A bidirectional iterator type (input and output iterator).
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 */
template <typename BidirIter>
inline void merge_sort(BidirIter first, BidirIter last) {
  merge_sort(first, last, std::less<>());
}

}  // namespace ReaK::sorting

#endif  // REAK_CORE_SORTING_MERGE_SORT_H_
