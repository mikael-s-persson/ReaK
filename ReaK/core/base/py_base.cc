/**
 *\file py_base.cpp
 *
 * This source file defines export functions for the python bindings on base classes
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


#include "ReaK/core/base/named_object.h"
#include "ReaK/core/base/shared_object.h"

#include "boost/python.hpp"

using namespace boost::python;

namespace PyReaK {

void export_base() {

  class_<ReaK::typed_object, std::shared_ptr<ReaK::typed_object>,
         boost::noncopyable>("TypedObj", no_init);

  class_<ReaK::serializable, bases<ReaK::typed_object>,
         std::shared_ptr<ReaK::serializable>, boost::noncopyable>(
      "SerializableObj", no_init);

  class_<ReaK::shared_object, bases<ReaK::serializable>,
         std::shared_ptr<ReaK::shared_object>, boost::noncopyable>("SharedObj",
                                                                   no_init);

  class_<ReaK::named_interface, std::shared_ptr<ReaK::named_interface>,
         boost::noncopyable>("NamedInterface", no_init)
      .def("get_name", pure_virtual(&ReaK::named_interface::get_name),
           return_value_policy<copy_const_reference>())
      .def("set_name", pure_virtual(&ReaK::named_interface::set_name));

  class_<ReaK::named_object, bases<ReaK::shared_object, ReaK::named_interface>,
         std::shared_ptr<ReaK::named_object>>("NamedObj");
  def("create_named_obj", std::make_shared<ReaK::named_object>);

  implicitly_convertible<std::shared_ptr<ReaK::named_object>,
                         std::shared_ptr<ReaK::shared_object>>();
  implicitly_convertible<std::shared_ptr<ReaK::named_object>,
                         std::shared_ptr<ReaK::named_interface>>();
  implicitly_convertible<std::shared_ptr<ReaK::shared_object>,
                         std::shared_ptr<ReaK::serializable>>();
  implicitly_convertible<std::shared_ptr<ReaK::serializable>,
                         std::shared_ptr<ReaK::typed_object>>();
}

}  // namespace PyReaK

BOOST_PYTHON_MODULE(libreak_py_base) {
  using namespace PyReaK;
  export_base();
}
