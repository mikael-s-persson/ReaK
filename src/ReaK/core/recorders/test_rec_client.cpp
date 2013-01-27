
/*
 *    Copyright 2011 Sven Mikael Persson
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

#include "tcp_recorder.hpp"

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace ReaK;
using namespace recorder;

int main(int argc, char** argv) {
  
  std::string ep_name = "127.0.0.1:17017";
  if(argc == 2)
    ep_name = argv[1];
  
  {
  tcp_extractor input_rec(ep_name);
  std::cout << "Connection established!" << std::endl;
  
  unsigned int col_count = input_rec.getColCount();
  std::cout << "Obtained the following names: ";
  for(unsigned int i = 0; i < col_count; ++i) {
    std::string tmp;
    input_rec >> tmp;
    std::cout << tmp << " ";
  };
  std::cout << std::endl;
  
  while(true) {
    try {
      for(unsigned int i = 0; i < col_count; ++i) {
        double tmp;
        input_rec >> tmp;
        std::cout << std::setw(10) << std::setprecision(4) << tmp;
      };
      input_rec >> data_extractor::end_value_row;
      std::cout << std::endl;
    } catch(...) {
      break;
    };
  };
  };
  
  return 0;
};












