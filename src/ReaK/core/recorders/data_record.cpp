
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

#include "data_record.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>

namespace ReaK {

namespace recorder {



void data_recorder::record_process::operator()() {
  std::size_t currentIter = 0;
  std::size_t iterStep = 1;
  boost::posix_time::ptime last_time = boost::posix_time::microsec_clock::local_time();
  while(parent.colCount != 0) {
    {
      if((parent.rowCount > parent.maxBufferSize) || (currentIter % iterStep == 0))
	parent.writeRow();
    };
    if(currentIter == 1000) {
      boost::posix_time::ptime current_time = boost::posix_time::microsec_clock::local_time();
      boost::posix_time::time_duration dt = current_time - last_time;
      std::size_t numSamples = static_cast<std::size_t>(1 + (parent.flushSampleRate * dt.total_microseconds()) / 1000000);
      iterStep = 1000 / numSamples;
      last_time = current_time;
      currentIter = 0;
    };
    ++currentIter;
    ReaKaux::this_thread::yield();
  };
};


data_recorder::~data_recorder() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  colCount = 0;
  if(writing_thread) {
    lock_here.unlock();
    writing_thread->join();
    lock_here.lock();
    writing_thread = ReaK::shared_ptr<ReaKaux::thread>();
  };
};


data_recorder& data_recorder::operator <<(double value) {
  if(colCount != 0) {
    if(currentColumn < colCount) {
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      values_rm.push(value);
      ++currentColumn;
    } else
      throw out_of_bounds();
  };
  return *this;
};

data_recorder& data_recorder::operator <<(const std::string& name) {
  if(colCount == 0) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    names.push_back(name);
  };
  return *this;
};

data_recorder& data_recorder::operator <<(flag some_flag) {
  if(some_flag == end_name_row) {
    if(colCount != 0)
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    colCount = names.size();
    writeNames();
    currentColumn = 0;
    rowCount = 0;
    writing_thread = ReaK::shared_ptr<ReaKaux::thread>(new ReaKaux::thread(record_process(*this)));
  } else if(some_flag == end_value_row) {
    if(colCount == 0)
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    for(;currentColumn < colCount;++currentColumn)
      values_rm.push(0.0);
    currentColumn = 0;
    ++rowCount;
  } else if(some_flag == flush) {
    //flush all data right away... normally would be done at the closure or pause...
    //not while doing other things because the function will not return until this is done.
    while(rowCount > 0)
      writeRow();
  } else if(some_flag == close) {
    //flush and stop thread.
    while(rowCount > 0)
      writeRow();
    colCount = 0;
    if(writing_thread) {
      writing_thread->join();
      writing_thread = ReaK::shared_ptr<ReaKaux::thread>();
    };
  };
  return *this;
};













void data_extractor::extract_process::operator()() {
  std::size_t currentIter = 0;
  std::size_t iterStep = 1;
  boost::posix_time::ptime last_time = boost::posix_time::microsec_clock::local_time();
  while(parent.colCount != 0) {
    {
      if((parent.values_rm.size() < parent.minBufferSize * parent.colCount) || 
	 (currentIter % iterStep == 0))
	if(!parent.readRow())
	  break;
    };
    if(currentIter == 1000) {
      boost::posix_time::ptime current_time = boost::posix_time::microsec_clock::local_time();
      boost::posix_time::time_duration dt = current_time - last_time;
      std::size_t numSamples = static_cast<std::size_t>(1 + (parent.flushSampleRate * dt.total_microseconds()) / 1000000);
      iterStep = 1000 / numSamples;
      last_time = current_time;
      currentIter = 0;
    };
    ++currentIter;
    ReaKaux::this_thread::yield();
  };
};


data_extractor::~data_extractor() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
  colCount = 0;
  if((reading_thread) && (reading_thread->joinable())) {
    lock_here.unlock();
    reading_thread->join();
    lock_here.lock();
    reading_thread.reset();
  };
};


data_extractor& data_extractor::operator >>(double& value) {
  if(colCount != 0) {
    if(currentColumn < colCount) {
      if(values_rm.empty()) {
	while((values_rm.size() < minBufferSize) && (readRow())) ;
	if(values_rm.empty())
	  throw end_of_record();
      };
      ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
      value = values_rm.front();
      values_rm.pop();
      ++currentColumn;
    } else
      throw out_of_bounds();
  } else
    throw end_of_record();
  return *this;
};

data_extractor& data_extractor::operator >>(std::string& name) {
  if(currentNameCol < colCount) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    name = names[currentNameCol++];
  } else
    throw out_of_bounds();
  return *this;
};

data_extractor& data_extractor::operator >>(flag some_flag) {
  if(some_flag == end_value_row) {
    if(colCount == 0)
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    for(;currentColumn < colCount;++currentColumn)
      values_rm.pop();
    currentColumn = 0;
  } else if(some_flag == advance) {
    //flush all data right away... normally would be done at the closure or pause...
    //not while doing other things because the function will not return until this is done.
    while((values_rm.size() < minBufferSize) && (readRow())) ;
  } else if(some_flag == close) {
    //flush and stop thread.
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here(access_mutex);
    while(!values_rm.empty())
      values_rm.pop();
    colCount = 0;
    if(reading_thread) {
      lock_here.unlock();
      if(reading_thread->joinable())
        reading_thread->join();
      lock_here.lock();
      reading_thread.reset();
    };
  };
  return *this;
};

void data_extractor::setFileName(const std::string& aFileName) {
  colCount = 0;
  names.resize(0);
  while(!values_rm.empty())
    values_rm.pop();
  if((loadFile(aFileName)) && (readNames())) {
    fileName = aFileName;
    currentColumn = 0;
    reading_thread = ReaK::shared_ptr<ReaKaux::thread>(new ReaKaux::thread(extract_process(*this)));
  };
};




};


};








