
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

#include "ReaK/core/sorting/bubble_sort.hpp"
#include "ReaK/core/sorting/comb_sort.hpp"
#include "ReaK/core/sorting/heap_sort.hpp"
#include "ReaK/core/sorting/insertion_sort.hpp"
#include "ReaK/core/sorting/intro_sort.hpp"
#include "ReaK/core/sorting/merge_sort.hpp"
#include "ReaK/core/sorting/quick_sort.hpp"
#include "ReaK/core/sorting/selection_sort.hpp"
#include "ReaK/core/sorting/shell_sort.hpp"

#include <cmath>
#include <ctime>

#include "gtest/gtest.h"

TEST(SortingTests, AllTests) {
  using namespace ReaK;
  using namespace sorting;

  std::vector<int> orig_values(100);
  for (std::size_t i = 0; i < 100; ++i) {
    orig_values[i] = std::rand() % 10000;
  }

  std::vector<int> ref_values = orig_values;
  std::sort(ref_values.begin(), ref_values.end());

  std::vector<int> test_values = orig_values;
  EXPECT_NO_THROW(selection_sort(test_values.begin(), test_values.end()));
  bool is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "selection_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(insertion_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "insertion_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(bubble_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "bubble_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(merge_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "merge_sort algorithm";

#if 0
  test_values = orig_values;
  EXPECT_NO_THROW(shell_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for(std::size_t i = 0; i < 100; ++i) {
    if(test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "shell_sort algorithm";
#endif

  test_values = orig_values;
  EXPECT_NO_THROW(comb_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "comb_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(heap_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "heap_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(quick_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "quick_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(quick_sort(test_values.begin(), test_values.end(),
                             std::less<>(), random_pivot()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "quick_sort random_pivot algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(quick_sort(test_values.begin(), test_values.end(),
                             std::less<>(), first_pivot()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "quick_sort first_pivot algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(quickselect_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "quickselect_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(intro_sort(test_values.begin(), test_values.end()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "intro_sort algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(intro_sort(test_values.begin(), test_values.end(),
                             std::less<>(), random_pivot()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "intro_sort random_pivot algorithm";

  test_values = orig_values;
  EXPECT_NO_THROW(intro_sort(test_values.begin(), test_values.end(),
                             std::less<>(), first_pivot()));
  is_correct = true;
  for (std::size_t i = 0; i < 100; ++i) {
    if (test_values[i] != ref_values[i]) {
      is_correct = false;
      break;
    }
  }
  EXPECT_TRUE(is_correct) << "intro_sort first_pivot algorithm";
}
