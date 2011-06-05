
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

//#include "base/rk_named_object.hpp"

namespace ReaK {

namespace recorder {



void data_recorder::record_process::operator()() {
  unsigned int currentIter = 0;
  unsigned int iterStep = 1;
  boost::posix_time::ptime last_time = boost::posix_time::microsec_clock::local_time();
  while(parent.colCount != 0) {
    {
      boost::unique_lock< boost::mutex > lock_here(parent.access_mutex);
      if((parent.rowCount > parent.maxBufferSize) || (currentIter % iterStep == 0))
	parent.writeRow();
    };
    if(currentIter == 1000) {
      boost::posix_time::ptime current_time = boost::posix_time::microsec_clock::local_time();
      boost::posix_time::time_duration dt = current_time - last_time;
      unsigned int numSamples = 1 + (parent.flushSampleRate * dt.total_microseconds()) / 1000000;
      iterStep = 1000 / numSamples;
      last_time = current_time;
      currentIter = 0;
    };
    ++currentIter;
    boost::this_thread::yield();
  };
};


data_recorder::~data_recorder() {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  colCount = 0;
  if((writing_thread) && (writing_thread->joinable())) {
    lock_here.unlock();
    if(writing_thread->get_id() != boost::this_thread::get_id())
      writing_thread->join();
    lock_here.lock();
    writing_thread.reset();
  };
};


data_recorder& data_recorder::operator <<(double value) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(colCount != 0) {
    if(currentColumn < colCount) {
      values_rm.push(value);
      ++currentColumn;
    } else
      throw out_of_bounds();
  };
  return *this;
};

data_recorder& data_recorder::operator <<(const std::string& name) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(colCount == 0) {
    names.push_back(name);
  };
  return *this;
};

data_recorder& data_recorder::operator <<(flag some_flag) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(some_flag == end_name_row) {
    if(colCount != 0)
      throw improper_flag();
    colCount = names.size();
    writeNames();
    currentColumn = 0;
    rowCount = 0;
    writing_thread = boost::shared_ptr<boost::thread>(new boost::thread(record_process(*this)));
  } else if(some_flag == end_value_row) {
    if(colCount == 0)
      throw improper_flag();
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
      lock_here.unlock();
      if(writing_thread->joinable())
        writing_thread->join();
      lock_here.lock();
      writing_thread.reset();
    };
  };
  return *this;
};













void data_extractor::extract_process::operator()() {
  unsigned int currentIter = 0;
  unsigned int iterStep = 1;
  boost::posix_time::ptime last_time = boost::posix_time::microsec_clock::local_time();
  while(parent.colCount != 0) {
    {
      boost::unique_lock< boost::mutex > lock_here(parent.access_mutex);
      if((parent.values_rm.size() < parent.minBufferSize * parent.colCount) || 
	 (currentIter % iterStep == 0))
	if(!parent.readRow())
	  break;
    };
    if(currentIter == 1000) {
      boost::posix_time::ptime current_time = boost::posix_time::microsec_clock::local_time();
      boost::posix_time::time_duration dt = current_time - last_time;
      unsigned int numSamples = 1 + (parent.flushSampleRate * dt.total_microseconds()) / 1000000;
      iterStep = 1000 / numSamples;
      last_time = current_time;
      currentIter = 0;
    };
    ++currentIter;
    boost::this_thread::yield();
  };
};


data_extractor::~data_extractor() {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  colCount = 0;
  if((reading_thread) && (reading_thread->joinable())) {
    lock_here.unlock();
    if(reading_thread->get_id() != boost::this_thread::get_id())
      reading_thread->join();
    lock_here.lock();
    reading_thread.reset();
  };
};


data_extractor& data_extractor::operator >>(double& value) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(colCount != 0) {
    if(currentColumn < colCount) {
      value = values_rm.front();
      values_rm.pop();
      ++currentColumn;
    } else
      throw out_of_bounds();
  };
  return *this;
};

data_extractor& data_extractor::operator >>(std::string& name) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(colCount == 0) {
    names.push_back(name);
  };
  return *this;
};

data_extractor& data_extractor::operator >>(flag some_flag) {
  boost::unique_lock< boost::mutex > lock_here(access_mutex);
  if(some_flag == end_value_row) {
    if(colCount == 0)
      throw improper_flag();
    for(;currentColumn < colCount;++currentColumn)
      values_rm.pop();
    currentColumn = 0;
  } else if(some_flag == advance) {
    //flush all data right away... normally would be done at the closure or pause...
    //not while doing other things because the function will not return until this is done.
    while((values_rm.size() < minBufferSize) && (readRow())) ;
  } else if(some_flag == close) {
    //flush and stop thread.
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
    reading_thread = boost::shared_ptr<boost::thread>(new boost::thread(extract_process(*this)));
  };
};




};


};








