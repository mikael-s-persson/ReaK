
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

namespace ReaK {

namespace recorder {


void ssv_recorder::writeRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((output_file.is_open()) && (rowCount > 0) && (colCount > 0)) {
    output_file << std::endl;
    output_file << values_rm.front();
    values_rm.pop();
    for(unsigned int i=1;i<colCount;++i) {
      output_file << " " << values_rm.front();
      values_rm.pop();
    };
    --rowCount;
  };
};

void ssv_recorder::writeNames() {
  output_file << "%";
  std::vector<std::string>::iterator it = names.begin();
  for(;it != names.end(); ++it)
    output_file << " " << (*it);
  output_file.flush();
};

void ssv_recorder::setFileName(const std::string& aFileName) {
  if(colCount != 0) {
    *this << close;

    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    if(output_file.is_open())
      output_file.close();
    output_file.open(aFileName.c_str());
    fileName = aFileName;
    if(output_file.is_open()) {
      colCount = names.size();
      writeNames();
    };
  } else {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    if(output_file.is_open())
      output_file.close();
    output_file.open(aFileName.c_str());
    fileName = aFileName;
  };
};




bool ssv_extractor::readRow() {
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

bool ssv_extractor::readNames() {
  std::string temp;
  std::getline(input_file,temp,'\n');
  std::stringstream ss(temp);
  std::string temp_name;
  ss >> temp_name; //ignore the first %
  while(ss >> temp_name) {
    names.push_back(temp_name);
  };
  colCount = names.size();
  return true;
};

bool ssv_extractor::loadFile(const std::string& aFileName) {
  if(colCount != 0) {
    *this >> close;

    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    if(input_file.is_open())
      input_file.close();
    input_file.open(aFileName.c_str());
    fileName = aFileName;
  } else {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    if(input_file.is_open())
      input_file.close();
    input_file.open(aFileName.c_str());
    fileName = aFileName;
  };
  return true;
};



};


};





