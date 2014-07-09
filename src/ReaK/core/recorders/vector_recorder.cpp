
/*
 *    Copyright 2014 Sven Mikael Persson
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

#include <ReaK/core/recorders/vector_recorder.hpp>

namespace ReaK {

namespace recorder {


vector_recorder::vector_recorder() : data_recorder(), vec_data(NULL) { };

vector_recorder::vector_recorder(std::vector< std::vector<double> >* aVecData) : data_recorder(), vec_data(aVecData) {
  setFileName("");
};

vector_recorder::~vector_recorder() { };

void vector_recorder::writeRow() {
  if(!vec_data)
    return;
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if((rowCount > 0) && (colCount > 0)) {
    std::vector<double> v_tmp(colCount, 0.0);
    for(std::size_t i = 0; i < colCount; ++i) {
      v_tmp[i] = values_rm.front();
      values_rm.pop();
    };
    vec_data->push_back(v_tmp);
    --rowCount;
  };
};

void vector_recorder::writeNames() { };

void vector_recorder::setFileName(const std::string& aFileName) { };

void vector_recorder::setVecData(std::vector< std::vector<double> >* aVecData) {
  vec_data = aVecData;
};





vector_extractor::vector_extractor() : data_extractor(), vec_data(NULL), cur_vec_index(0) { };

vector_extractor::vector_extractor(const std::vector< std::vector<double> >* aVecData) : data_extractor(), vec_data(aVecData), cur_vec_index(0) {
  
};

vector_extractor::~vector_extractor() {};


void vector_extractor::addName(const std::string& s) { 
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  names.push_back(s); 
  ++colCount;
};


bool vector_extractor::readRow() {
  if(!vec_data || (cur_vec_index >= vec_data->size()))
    return false;
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  if(colCount > 0) {
    const std::vector<double>& v_tmp = (*vec_data)[cur_vec_index];
    for(std::size_t i = 0; (i < colCount) && (i < v_tmp.size()); ++i)
      values_rm.push(v_tmp[i]);
    ++cur_vec_index;
  };
  return true;
};

bool vector_extractor::readNames() {
  return true;
};

void vector_extractor::setFileName(const std::string& aFileName) { };

void vector_extractor::setVecData(const std::vector< std::vector<double> >* aVecData) {
  vec_data = aVecData;
  cur_vec_index = 0;
};


};


};



