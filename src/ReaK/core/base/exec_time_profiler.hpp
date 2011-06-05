
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

#ifndef RK_EXEC_TIME_PROFILER_HPP
#define RK_EXEC_TIME_PROFILER_HPP

#include "defs.hpp"

#include <vector>
#include <fstream>
#include <ios>

#include <boost/date_time.hpp>

namespace ReaK {


#if RK_VERBOSITY>7
class exec_time_profiler {
  private:
    std::vector< boost::posix_time::ptime > mTimeTable;
    std::vector< int > mIDTable;
    
    std::ofstream mOutputStream;
    
    uint mCurrentIndex;
    uint mTotalCount;
    
  public:
    
    exec_time_profiler(const std::string& aFileName) : mTimeTable(),
                                                       mIDTable(),
                                                       mOutputStream(aFileName.c_str(),std::ios_base::out | std::ios_base::trunc),
                                                       mCurrentIndex(0),
                                                       mTotalCount(0) { };
    virtual ~exec_time_profiler() { mOutputStream.flush(); };
    
    void markTime(int Index) {
      if(mTotalCount > 0) {
        if(mCurrentIndex < mTotalCount)
          mTimeTable[mCurrentIndex] = boost::posix_time::microsec_clock::local_time();
      } else {
        mTimeTable.push_back(boost::posix_time::microsec_clock::local_time());
        mIDTable.push_back(Index);
      };
      mCurrentIndex++;
    };
    
    void flush() {
      if(mTotalCount == 0) {
        mTotalCount = mTimeTable.size();
        for(uint i=0;i<mTotalCount-1;++i)
          mOutputStream << mIDTable[i] << " to " << mIDTable[i+1] << "\t";
      };
      mCurrentIndex = 0;
      mOutputStream << std::endl;
      for(uint i=0;i<mTotalCount-1;++i) {
        boost::posix_time::time_duration delta_t = mTimeTable[i+1] - mTimeTable[i];
        mOutputStream << delta_t.total_microseconds() << "\t";
      };
      mOutputStream.flush();
    };

};
#else
class exec_time_profiler {
  public:
    
    exec_time_profiler(const std::string& aFileName) { };
    virtual ~exec_time_profiler() { };
    
    void markTime(int Index) { };
    
    void flush() { };

};
#endif


};


#define RK_EXEC_TIME_PROFILER_CREATE ReaK::exec_time_profiler mMyExecTimeProfiler;

#define RK_EXEC_TIME_PROFILER_INIT(CONTAINER_NAME) mMyExecTimeProfiler(std::string("profiles/") + CONTAINER_NAME + "_exec_times.csv")

#define RK_EXEC_TIME_PROFILER_MARK mMyExecTimeProfiler.markTime(__LINE__);

#define RK_EXEC_TIME_PROFILER_FLUSH mMyExecTimeProfiler.flush();


#endif //RK_EXEC_TIME_PROFILER_HPP














