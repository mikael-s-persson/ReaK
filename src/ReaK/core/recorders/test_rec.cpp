
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

#include "ssv_recorder.hpp"
#include "tsv_recorder.hpp"
#include "serialization/xml_archiver.hpp"

using namespace ReaK;
using namespace recorder;
using namespace serialization;

int main() {
  
  
  std::cout << RK_HERE << " reached!" << std::endl;
  {
  ssv_recorder output_rec("test_data.sscdat");
  std::cout << RK_HERE << " reached!" << std::endl;
  output_rec << "x" << "2*x" << "x^2" << data_recorder::end_name_row;
  std::cout << RK_HERE << " reached!" << std::endl;
  for(double x=0;x < 10; x += 0.5) {
    output_rec << x << 2*x << x*x << data_recorder::end_value_row;
  };
  std::cout << RK_HERE << " reached!" << std::endl;
  output_rec << data_recorder::flush;
  
  xml_oarchive output_arc("ssv_rec_archive.xml");
  output_arc << output_rec;
  };
  
  std::cout << RK_HERE << " reached!" << std::endl;
  {
#ifndef RK_ENABLE_CXX0X_FEATURES
  boost::shared_ptr<serializable> output_rec;
#else
  std::shared_ptr<serializable> output_rec;
#endif
  
  {
  std::cout << RK_HERE << " reached!" << std::endl;
  xml_iarchive input_arc("ssv_rec_archive.xml");
  std::cout << RK_HERE << " reached!" << std::endl;
  input_arc >> output_rec;
  };
#ifndef RK_ENABLE_CXX0X_FEATURES
  boost::shared_ptr<ssv_recorder> output_rec_ssv = rtti::rk_dynamic_ptr_cast<ssv_recorder>(output_rec);
#else
  std::shared_ptr<ssv_recorder> output_rec_ssv = rtti::rk_dynamic_ptr_cast<ssv_recorder>(output_rec);
#endif
  
  std::cout << RK_HERE << " reached!" << std::endl;
  for(double x=10;x < 20; x += 0.5) {
    (*output_rec_ssv) << x << 2*x << x*x << data_recorder::end_value_row;
  };
  std::cout << RK_HERE << " reached!" << std::endl;
  (*output_rec_ssv) << data_recorder::flush;
  std::cout << RK_HERE << " reached!" << std::endl;
  };
  std::cout << RK_HERE << " reached!" << std::endl;
  return 0;
};












