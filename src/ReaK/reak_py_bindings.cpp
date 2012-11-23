/**
 *\file reak_py_bindings.cpp
 *
 * This source file defines the list of exports for the python bindings on the ReaK platform.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date June 2012
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


#include <boost/python.hpp>

namespace PyReaK {

void export_base();
void export_vect_alg();
void export_kinetostatics();
void export_mbd_kte();
void export_kte_models();

};

using namespace boost::python;


BOOST_PYTHON_MODULE(libreak_py) {
  
  using namespace PyReaK;
  
  export_base();
  export_vect_alg();
  export_kinetostatics();
  export_mbd_kte();
  export_kte_models();
  
};

















