
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

#include <ReaK/core/base/defs.hpp>
#include <ReaK/math/lin_alg/mat_alg.hpp>

#include <iostream>
#include <fstream>
#include <cstdio>


/*
 * This program is just a little test program that is used temporarily to test little
 * bits of code here and there related to the lin-alg library. In other words, if I
 * don't know if a particular expression is going to compile and compute correctly,
 * I can just test it in this "sand-box".
 */

int main() {
  using namespace ReaK;

  mat< double, mat_structure::rectangular > m( 2, 2 );
  m( 0, 0 ) = 1.0;
  m( 0, 1 ) = 0.0;
  m( 1, 0 ) = 0.0;
  m( 1, 1 ) = 1.0;

  std::cout << m << std::endl;

  return 0;
};
