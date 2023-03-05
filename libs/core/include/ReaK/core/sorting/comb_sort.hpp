/**
 * \file comb_sort.hpp
 *
 * This library provides a generic comb-sort function.
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

#ifndef REAK_COMB_SORT_HPP
#define REAK_COMB_SORT_HPP

#include "ReaK/core/base/defs.hpp"

#include <algorithm>
#include <functional>
#include <iterator>

namespace ReaK::sorting {

/**
 * This function performs a comb sort on a given range of elements, and with the
 * given comparison functor.
 * \tparam ForwardIter A foward iterator type (input and output iterator).
 * \tparam Compare A comparison functor type that can produce a bool value to order two elements.
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 * \param comp The comparison functor to use to determine the order of elements.
 */
template <typename ForwardIter, typename Compare>
void comb_sort(ForwardIter first, ForwardIter last, Compare comp) {
  std::size_t gap = std::distance(first, last);
  bool swapped = true;
  while ((gap > 1) || swapped) {
    if (gap > 1) {
      gap = (gap * 4) / 5;
    }
    swapped = false;
    ForwardIter current_next = first;
    std::advance(current_next, gap);
    for (ForwardIter current = first; current_next != last;
         ++current, ++current_next) {
      if (comp(*current_next, *current) && (swapped = true)) {
        std::iter_swap(current, current_next);
      }
    }
  }
}

/**
 * This function performs a comb sort on a given range of elements, and by using the
 * less-than operator on the elements of the range.
 * \tparam ForwardIter A foward iterator type (input and output iterator).
 * \param first The start of the range to be sorted.
 * \param last One element past the end of the range to be sorted.
 */
template <typename ForwardIter>
inline void comb_sort(ForwardIter first, ForwardIter last) {
  comb_sort(first, last, std::less<>());
}

}  // namespace ReaK::sorting

#endif  // REAK_COMB_SORT_HPP
