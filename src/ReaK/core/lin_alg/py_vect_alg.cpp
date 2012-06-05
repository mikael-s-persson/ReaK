/**
 *\file py_vect_alg.cpp
 *
 * This source file defines export functions for the python bindings on vect_alg classes 
 * of the ReaK platform.
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



#include "base/defs.hpp"

#include "vect_alg.hpp"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "base/py_fixes.hpp"

#include <sstream>

namespace PyReaK {
  
 
template <typename Vector>
double vect_getitem(const Vector& v, std::size_t i) {
  return v[i]; 
};

template <typename Vector>
void vect_setitem(Vector& v, std::size_t i, double d) {
  v[i] = d;
};

template <typename Vector>
std::string vect_to_string(const Vector& v) {
  std::stringstream ss;
  ss << v;
  return ss.str();
};


void export_vect_alg() {

  using namespace boost::python;
  
  class_< ReaK::vect<double,2>,
          bases< ReaK::serialization::serializable >
        >("Vector2D")
    .def(init<double,double>())
    .def(self + self)
    .def(self - self)
    .def(-self)
    .def(self * self)
    .def(self += self)
    .def(self -= self)
    .def(self * double())
    .def(double() * self)
    .def(self / double())
    .def(self *= double())
    .def(self /= double())
    .def("__str__",vect_to_string< ReaK::vect<double,2> >)
    .def("__len__",&ReaK::vect<double,2>::size)
    .def("__getitem__",vect_getitem< ReaK::vect<double,2> >)
    .def("__setitem__",vect_setitem< ReaK::vect<double,2> >);
  
  class_< ReaK::vect<double,3>,
          bases< ReaK::serialization::serializable >
        >("Vector3D")
    .def(init<double,double,double>())
    .def(self + self)
    .def(self - self)
    .def(-self)
    .def(self * self)
    .def(self += self)
    .def(self -= self)
    .def(self * double())
    .def(double() * self)
    .def(self / double())
    .def(self *= double())
    .def(self /= double())
    .def("__str__",vect_to_string< ReaK::vect<double,3> >)
    .def("__len__",&ReaK::vect<double,3>::size)
    .def("__getitem__",vect_getitem< ReaK::vect<double,3> >)
    .def("__setitem__",vect_setitem< ReaK::vect<double,3> >);
  
  class_< ReaK::vect<double,4>,
          bases< ReaK::serialization::serializable >
        >("Vector4D")
    .def(init<double,double,double,double>())
    .def(self + self)
    .def(self - self)
    .def(-self)
    .def(self * self)
    .def(self += self)
    .def(self -= self)
    .def(self * double())
    .def(double() * self)
    .def(self / double())
    .def(self *= double())
    .def(self /= double())
    .def("__str__",vect_to_string< ReaK::vect<double,4> >)
    .def("__len__",&ReaK::vect<double,4>::size)
    .def("__getitem__",vect_getitem< ReaK::vect<double,4> >)
    .def("__setitem__",vect_setitem< ReaK::vect<double,4> >);
  
  class_< ReaK::vect_n<double, std::allocator<double> >,
          bases< ReaK::serialization::serializable >
        >("VectorND")
    .def(init<double,double,double>())
    .def(init<double,double,double,double>())
    .def(init<double,double,double,double,double>())
    .def(init<double,double,double,double,double,double>())
    .def(init<double,double,double,double,double,double,double>())
    .def(init<double,double,double,double,double,double,double,double>())
    .def(init<double,double,double,double,double,double,double,double,double>())
    .def(init<double,double,double,double,double,double,double,double,double,double>())
    .def(init<double,double,double,double,double,double,double,double,double,double,
	      double>())
    .def(init<double,double,double,double,double,double,double,double,double,double,
	      double,double>())
    .def(init<double,double,double,double,double,double,double,double,double,double,
	      double,double,double>())
    .def(init<double,double,double,double,double,double,double,double,double,double,
	      double,double,double,double>())
    .def(self + self)
    .def(self - self)
    .def(-self)
    .def(self * self)
    .def(self += self)
    .def(self -= self)
    .def(self * double())
    .def(double() * self)
    .def(self / double())
    .def(self *= double())
    .def(self /= double())
    .def("__str__",vect_to_string< ReaK::vect_n<double, std::allocator<double> > >)
    .def("__len__",&ReaK::vect_n<double, std::allocator<double> >::size)
    .def("__getitem__",vect_getitem< ReaK::vect_n<double, std::allocator<double> > >)
    .def("__setitem__",vect_setitem< ReaK::vect_n<double, std::allocator<double> > >)
    .def("resize", &ReaK::vect_n<double, std::allocator<double> >::resize);
  
  
};


};
















