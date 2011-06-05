
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

#include "bin_recorder.hpp"

namespace ReaK {

namespace recorder {


void bin_recorder::writeRow() {
  if((output_file.is_open()) && (rowCount > 0) && (colCount > 0)) {
    for(unsigned int i=0;i<colCount;++i) {
      double tmp(values_rm.front());
      output_file.write(reinterpret_cast<char*>(&tmp),sizeof(double));
      values_rm.pop();
    };
    --rowCount;
  };
};

void bin_recorder::writeNames() {
  output_file.write(reinterpret_cast<char*>(&colCount),sizeof(unsigned int));
  std::vector<std::string>::iterator it = names.begin();
  for(;it != names.end(); ++it)
    output_file.write(it->c_str(),it->size() + 1);
};

void bin_recorder::setFileName(const std::string& aFileName) {
  if(colCount != 0) {
    *this << close;

    boost::unique_lock< boost::mutex > lock_here(access_mutex);
    if(output_file.is_open())
      output_file.close();
    output_file.open(aFileName.c_str(),std::ios_base::out | std::ios_base::binary);
    fileName = aFileName;
    if(output_file.is_open()) {
      colCount = names.size();
      writeNames();
    };
  } else {
    boost::unique_lock< boost::mutex > lock_here(access_mutex);
    if(output_file.is_open())
      output_file.close();
    output_file.open(aFileName.c_str(),std::ios_base::out | std::ios_base::binary);
    fileName = aFileName;
  };
};



bool bin_extractor::readRow() {
  if((input_file.is_open()) && (colCount > 0)) {
    for(unsigned int i=0;i<colCount;++i) {
      double tmp = 0;
      input_file.read(reinterpret_cast<char*>(&tmp),sizeof(double));
      if(!input_file)
	return false;
      values_rm.push(tmp);
    };
  };
  return true;
};

bool bin_extractor::readNames() {
  input_file.read(reinterpret_cast<char*>(&colCount),sizeof(unsigned int));
  char temp[128];
  for(unsigned int i = 0; i < colCount; ++i) {
    char* temp_ptr = temp;
    while((temp_ptr < temp + 128) && (input_file.read(temp_ptr,1)) && (*temp_ptr != '\0'))
      ++temp_ptr;
    if((temp_ptr >= temp + 128) || (!input_file))
      return false;
    names.push_back(std::string(temp));
  };
  return true;
};

bool bin_extractor::loadFile(const std::string& aFileName) {
  if(colCount != 0) {
    *this >> close;

    boost::unique_lock< boost::mutex > lock_here(access_mutex);
    if(input_file.is_open())
      input_file.close();
    input_file.open(aFileName.c_str(),std::ios_base::in | std::ios_base::binary);
    fileName = aFileName;
  } else {
    boost::unique_lock< boost::mutex > lock_here(access_mutex);
    if(input_file.is_open())
      input_file.close();
    input_file.open(aFileName.c_str(),std::ios_base::in | std::ios_base::binary);
    fileName = aFileName;
  };
  return true;
};



};


};



