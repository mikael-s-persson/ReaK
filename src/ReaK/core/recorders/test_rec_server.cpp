
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

using namespace ReaK;
using namespace recorder;

int main(int argc, char** argv) {
  
  std::string portnum = "17017";
  if(argc == 2)
    portnum = argv[1];
  
  {
  tcp_recorder output_rec(portnum);
  std::cout << "Connection established." << std::endl;
  
  output_rec << "x" << "sin(x)" << "cos(x)" << data_recorder::end_name_row;
  std::cout << "Names sent." << std::endl;
  
  for(double x = 0.0; x < 10.0; x += 0.5) {
    output_rec << x << std::sin(x) << std::cos(x) << data_recorder::end_value_row;
  };
  output_rec << data_recorder::flush;
  std::cout << "Done!" << std::endl;
  char c;
//   std::cin >> c;
  };
  
  return 0;
};












