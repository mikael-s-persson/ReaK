
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

#include <ReaK/core/recorders/ssv_recorder.hpp>

namespace ReaK {

namespace recorder {


void ssv_recorder::writeRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((out_stream) && (*out_stream) && (rowCount > 0) && (colCount > 0)) {
    (*out_stream) << std::endl;
    (*out_stream) << values_rm.front();
    values_rm.pop();
    for(unsigned int i = 1; i < colCount; ++i) {
      (*out_stream) << " " << values_rm.front();
      values_rm.pop();
    };
    --rowCount;
  };
};

void ssv_recorder::writeNames() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((!out_stream) || (!(*out_stream)))
    return;
  (*out_stream) << "%";
  std::vector<std::string>::iterator it = names.begin();
  for(;it != names.end(); ++it)
    (*out_stream) << " " << (*it);
  out_stream->flush();
};

void ssv_recorder::setStreamImpl(const shared_ptr<std::ostream>& aStreamPtr) {
  if(colCount != 0) {
    *this << close;
    if((aStreamPtr) && (*aStreamPtr)) {
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      out_stream = aStreamPtr;
      out_stream->setf(std::ios::scientific, std::ios::floatfield);
      out_stream->precision(11);
      colCount = names.size();
      lock_here.unlock();
      writeNames();
    };
  } else {
    if((aStreamPtr) && (*aStreamPtr)) {
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      out_stream = aStreamPtr;
      out_stream->setf(std::ios::scientific, std::ios::floatfield);
      out_stream->precision(11);
    };
  };
};



bool ssv_extractor::readRow() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((in_stream) && (*in_stream) && (colCount > 0)) {
    std::string temp;
    std::getline(*in_stream, temp, '\n');
    if(!(*in_stream))
      return false;
    std::stringstream ss(temp);
    for(unsigned int i = 0; i < colCount; ++i) {
      double tmp = 0;
      ss >> tmp;
      if(!ss)
        return false;
      values_rm.push(tmp);
    };
  };
  if((in_stream) && !(*in_stream))
    return false;
  return true;
};

bool ssv_extractor::readNames() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((!in_stream) || (!(*in_stream)))
    return false;
  std::string temp;
  std::getline(*in_stream,temp,'\n');
  std::stringstream ss(temp);
  std::string temp_name;
  ss >> temp_name; //ignore the first %
  while(ss >> temp_name) {
    names.push_back(temp_name);
  };
  colCount = names.size();
  return true;
};

void ssv_extractor::setStreamImpl(const shared_ptr<std::istream>& aStreamPtr) {
  if(colCount != 0)
    *this >> close;
  if((aStreamPtr) && (*aStreamPtr)) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    in_stream = aStreamPtr;
    lock_here.unlock();
    readNames();
  };
};



};


};





