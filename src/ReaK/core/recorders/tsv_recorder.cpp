
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

#include "tsv_recorder.hpp"

namespace ReaK {

namespace recorder {

void tsv_recorder::writeRow() {
  if((output_file.is_open()) && (rowCount > 0) && (colCount > 0)) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    output_file << std::endl;
    output_file << values_rm.front();
    values_rm.pop();
    for(unsigned int i=1;i<colCount;++i) {
      output_file << "\t" << values_rm.front();
      values_rm.pop();
    };
    --rowCount;
  };
};

void tsv_recorder::writeNames() {
  output_file << "%";
  std::vector<std::string>::iterator it = names.begin();
  for(;it != names.end(); ++it)
    output_file << "\t" << (*it);
  output_file.flush();
};


bool tsv_extractor::readRow() {
  if((input_file.is_open()) && (colCount > 0)) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    std::string temp;
    std::getline(input_file,temp,'\n');
    std::stringstream ss(temp);
    for(unsigned int i=0;i<colCount;++i) {
      double tmp = 0;
      ss >> tmp;
      if(!ss)
	return false;
      values_rm.push(tmp);
    };
  };
  return true;
};


};

};





