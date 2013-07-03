
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

#include "selection_sort.hpp"
#include "insertion_sort.hpp"
#include "bubble_sort.hpp"
#include "shell_sort.hpp"
#include "comb_sort.hpp"
#include "merge_sort.hpp"
#include "heap_sort.hpp"
#include "quick_sort.hpp"
#include "intro_sort.hpp"

#include "base/chrono_incl.hpp"

#include <ctime>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>

const std::size_t RANGE = 100000;
std::size_t REPEATS = 10;

struct rand_int_generator {
  std::size_t m_range;
  rand_int_generator(std::size_t aRange) : m_range(aRange) { };
  int operator()() const { return rand() % m_range; };
};

void generate_random_seq(int* first, int* last) {
  std::generate(first, last, rand_int_generator(RANGE));
};

void generate_almost_sorted(int* first, int* last) {
  std::generate(first, last, rand_int_generator(RANGE));
  std::sort(first,last);
  std::size_t dist = (last - first);
  for(std::size_t i = 0; i < dist / 10; ++i)
    first[rand() % dist] = rand() % RANGE;
};

void generate_almost_reversed(int* first, int* last) {
  std::generate(first, last, rand_int_generator(RANGE));
  std::sort(first,last,std::greater<int>());
  std::size_t dist = (last - first);
  for(std::size_t i = 0; i < dist / 10; ++i)
    first[rand() % dist] = rand() % RANGE;
};

void generate_partly_sorted(int* first, int* last) {
  std::generate(first, last, rand_int_generator(RANGE));
  std::sort(first,last);
  std::size_t dist = (last - first);
  for(std::size_t i = dist - dist / 10; i < dist; ++i)
    first[i] = rand() % RANGE;
};

template <typename SortFunction, typename GenFunction>
std::size_t sort_and_get_us(std::size_t SIZE, SortFunction func, GenFunction gen) {
  
  using namespace ReaKaux::chrono;
  
  int* tmp_array = new int[SIZE];
  
  high_resolution_clock::time_point t_start = high_resolution_clock::now();
  high_resolution_clock::time_point t_end = t_start;
  auto t_duration = t_end - t_start;
  
  for(std::size_t i = 0; i < REPEATS; ++i) {
    gen(tmp_array, tmp_array + SIZE);
    t_start = high_resolution_clock::now();
    func(tmp_array, tmp_array + SIZE);
    t_end = high_resolution_clock::now();
    
    t_duration = t_duration + (t_end - t_start);
  };
  
  delete[] tmp_array;
  
  return (duration_cast<microseconds>(t_duration).count() / REPEATS);
};

int main(int argc, char** argv) {
  
  using namespace ReaK::sorting;
  
  std::size_t max_size = 1000000;
  if(argc > 1)
    std::stringstream(argv[1]) >> max_size;
  
  if(argc > 2)
    std::stringstream(argv[2]) >> REPEATS;
  
  std::srand(static_cast<unsigned int>(std::time(0)));
  
  typedef void (*GenFuncType)(int*,int*);
  
  GenFuncType generators[] = {&generate_random_seq,
                              &generate_almost_sorted,
                              &generate_almost_reversed,
                              &generate_partly_sorted};
  
  for(std::size_t i = 0; i < 4; ++i) {
    for(std::size_t j = 10; j < max_size; j *= 2) {
      std::cout << std::setw(10) << j;
      
      if(j < 10000) {
        std::cout << std::setw(10) << sort_and_get_us(j, selection_sort<int*>, generators[i]);
        
        std::cout << std::setw(10) << sort_and_get_us(j, insertion_sort<int*>, generators[i]);
        
        std::cout << std::setw(10) << sort_and_get_us(j, bubble_sort<int*>, generators[i]);
      } else {
        std::cout << std::setw(10) << 0;
        std::cout << std::setw(10) << 0;
        std::cout << std::setw(10) << 0;
      };
      
      std::cout << std::setw(10) << sort_and_get_us(j, std::sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, merge_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, shell_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, comb_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, heap_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, quick_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, std::bind(quick_sort<int*,std::less<int>,random_pivot>,std::placeholders::_1,std::placeholders::_2,std::less<int>(),random_pivot()), generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, std::bind(quick_sort<int*,std::less<int>,first_pivot>,std::placeholders::_1,std::placeholders::_2,std::less<int>(),first_pivot()), generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, quickselect_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, intro_sort<int*>, generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, std::bind(intro_sort<int*,std::less<int>,random_pivot>,std::placeholders::_1,std::placeholders::_2,std::less<int>(),random_pivot()), generators[i]);
      
      std::cout << std::setw(10) << sort_and_get_us(j, std::bind(intro_sort<int*,std::less<int>,first_pivot>,std::placeholders::_1,std::placeholders::_2,std::less<int>(),first_pivot()), generators[i]);
      
      std::cout << std::endl;
    };
  };
  std::cout << std::flush;
  
  return 0;
};




