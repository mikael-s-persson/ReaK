
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

#include <ReaK/core/recorders/data_record.hpp>

#include <ReaK/core/base/chrono_incl.hpp>

#include <fstream>

namespace ReaK {

namespace recorder {


namespace ch = ReaKaux::chrono;
typedef ch::high_resolution_clock hrc;

void data_recorder::closeRecordProcess() {
  while( rowCount > 0 )
    writeRow();
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  colCount = 0;
  if( writing_thread ) {
    lock_here.unlock();
    if( writing_thread->joinable() )
      writing_thread->join();
    lock_here.lock();
    writing_thread.reset();
  };
};

void data_recorder::record_process::operator()() {
  hrc::time_point last_time = hrc::now();
  while( parent.colCount != 0 ) {
    do {
      parent.writeRow();
    } while( parent.rowCount > parent.maxBufferSize );

    if( parent.flushSampleRate == 0 ) {
      ReaKaux::this_thread::yield();
    } else {
      double time_period = 1.0 / parent.flushSampleRate;
      hrc::time_point time_to_reach = last_time + ch::duration_cast< hrc::duration >(
                                                    ch::duration< double, ReaKaux::ratio< 1, 1 > >( time_period ) );
      ReaKaux::this_thread::sleep_until( time_to_reach );
    };
    last_time = hrc::now();
  };
};


data_recorder::~data_recorder() { closeRecordProcess(); };


data_recorder& data_recorder::operator<<( double value ) {
  if( colCount == 0 )
    throw end_of_record();
  if( currentColumn >= colCount )
    throw out_of_bounds();

  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  values_rm.push( value );
  ++currentColumn;
  return *this;
};

data_recorder& data_recorder::operator<<( const named_value_row& values ) {
  if( colCount == 0 )
    return *this;
  for( std::vector< double >::const_iterator it = values.values.begin(); it != values.values.end(); ++it )
    ( *this ) << ( *it );
  return ( *this ) << end_value_row;
};

data_recorder& data_recorder::operator<<( const std::string& name ) {
  if( colCount == 0 ) {
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
    names.push_back( name );
  };
  return *this;
};

data_recorder& data_recorder::operator<<( flag some_flag ) {
  if( some_flag == end_name_row ) {
    if( colCount != 0 )
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
    colCount = names.size();
    for( std::size_t i = 0; i < colCount; ++i )
      named_indices[names[i]] = i;
    lock_here.unlock();
    writeNames();
    lock_here.lock();
    currentColumn = 0;
    rowCount = 0;
    writing_thread = ReaK::shared_ptr< ReaKaux::thread >( new ReaKaux::thread( record_process( *this ) ) );
  } else if( some_flag == end_value_row ) {
    if( colCount == 0 )
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
    for( ; currentColumn < colCount; ++currentColumn )
      values_rm.push( 0.0 );
    currentColumn = 0;
    ++rowCount;
  } else if( some_flag == flush ) {
    // flush all data right away... normally would be done at the closure or pause...
    // not while doing other things because the function will not return until this is done.
    while( rowCount > 0 )
      writeRow();
  } else if( some_flag == close ) {
    // flush and stop thread.
    closeRecordProcess();
  };
  return *this;
};

void data_recorder::setFileName( const std::string& aFilename ) {
  shared_ptr< std::ofstream > file_out( new std::ofstream( aFilename.c_str() ) );
  if( file_out->is_open() )
    setStreamImpl( file_out );
};

void RK_CALL data_recorder::save( serialization::oarchive& A, unsigned int ) const {
  shared_object::save( A, shared_object::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_SAVE_WITH_NAME( static_cast< unsigned int >( colCount ) ) & RK_SERIAL_SAVE_WITH_NAME( flushSampleRate )
    & RK_SERIAL_SAVE_WITH_NAME( maxBufferSize ) & RK_SERIAL_SAVE_WITH_NAME( names );
};

void RK_CALL data_recorder::load( serialization::iarchive& A, unsigned int ) {
  closeRecordProcess();
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  shared_object::load( A, shared_object::getStaticObjectType()->TypeVersion() );
  unsigned int aColCount;
  A& RK_SERIAL_LOAD_WITH_ALIAS( "colCount", aColCount ) & RK_SERIAL_LOAD_WITH_NAME( flushSampleRate )
    & RK_SERIAL_LOAD_WITH_NAME( maxBufferSize ) & RK_SERIAL_LOAD_WITH_NAME( names );
  colCount = aColCount;
  for( std::size_t i = 0; i < colCount; ++i )
    named_indices[names[i]] = i;
  rowCount = 0;
  currentColumn = 0;
  values_rm = std::queue< double >();
  lock_here.unlock();
  writing_thread = ReaK::shared_ptr< ReaKaux::thread >( new ReaKaux::thread( record_process( *this ) ) );
};


void data_extractor::closeExtractProcess() {
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  while( !values_rm.empty() )
    values_rm.pop();
  colCount = 0;
  if( reading_thread ) {
    lock_here.unlock();
    if( reading_thread->joinable() )
      reading_thread->join();
    lock_here.lock();
    reading_thread.reset();
  };
};

void data_extractor::extract_process::operator()() {
  hrc::time_point last_time = hrc::now();
  while( parent.colCount != 0 ) {
    do {
      if( !parent.readRow() ) {
        return;
      };
    } while( parent.values_rm.size() < parent.minBufferSize * parent.colCount );
    if( parent.flushSampleRate == 0 ) {
      ReaKaux::this_thread::yield();
    } else {
      double time_period = 1.0 / parent.flushSampleRate;
      hrc::time_point time_to_reach = last_time + ch::duration_cast< hrc::duration >(
                                                    ch::duration< double, ReaKaux::ratio< 1, 1 > >( time_period ) );
      ReaKaux::this_thread::sleep_until( time_to_reach );
    };
    last_time = hrc::now();
  };
};


data_extractor::~data_extractor() { closeExtractProcess(); };


data_extractor& data_extractor::operator>>( double& value ) {
  if( colCount == 0 )
    throw end_of_record();
  if( currentColumn >= colCount )
    throw out_of_bounds();
  while( values_rm.size() < colCount - currentColumn ) {
    // This second size-test may look redundant,
    //  but values could have been extracted between the last test and this one
    //  (and failures have occurred because of this problem!)
    if( !readRow() && ( values_rm.size() < colCount - currentColumn ) )
      throw end_of_record();
  };
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  value = values_rm.front();
  values_rm.pop();
  ++currentColumn;
  return *this;
};

data_extractor& data_extractor::operator>>( std::string& name ) {
  if( currentNameCol >= colCount )
    throw out_of_bounds();

  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  name = names[currentNameCol++];
  return *this;
};

data_extractor& data_extractor::operator>>( named_value_row& values ) {
  if( colCount == 0 )
    return *this;
  for( std::vector< double >::iterator it = values.values.begin(); it != values.values.end(); ++it )
    ( *this ) >> ( *it );
  return ( *this ) >> end_value_row;
};

data_extractor& data_extractor::operator>>( flag some_flag ) {
  if( some_flag == end_value_row ) {
    if( colCount == 0 )
      throw improper_flag();
    ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
    for( ; currentColumn < colCount; ++currentColumn )
      values_rm.pop();
    currentColumn = 0;
  } else if( some_flag == advance ) {
    // flush all data right away... normally would be done at the closure or pause...
    // not while doing other things because the function will not return until this is done.
    std::size_t last_count = values_rm.size();
    do {
      last_count = values_rm.size();
      if( !readRow() )
        break;
    } while( ( values_rm.size() > last_count ) && ( values_rm.size() < minBufferSize * colCount ) );
  } else if( some_flag == close ) {
    // flush and stop thread.
    closeExtractProcess();
  };
  return *this;
};

void data_extractor::setStreamWrappedCall( const shared_ptr< std::istream >& aStreamPtr ) {
  setStreamImpl( aStreamPtr );
  for( std::size_t i = 0; i < colCount; ++i )
    named_indices[names[i]] = i;
  currentColumn = 0;
  reading_thread = ReaK::shared_ptr< ReaKaux::thread >( new ReaKaux::thread( extract_process( *this ) ) );
};

void data_extractor::setFileName( const std::string& aFilename ) {
  shared_ptr< std::ifstream > file_in( new std::ifstream( aFilename.c_str() ) );
  if( file_in->is_open() )
    setStreamWrappedCall( file_in );
};


void RK_CALL data_extractor::save( serialization::oarchive& A, unsigned int ) const {
  shared_object::save( A, shared_object::getStaticObjectType()->TypeVersion() );
  A& RK_SERIAL_SAVE_WITH_NAME( static_cast< unsigned int >( colCount ) ) & RK_SERIAL_SAVE_WITH_NAME( flushSampleRate )
    & RK_SERIAL_SAVE_WITH_NAME( minBufferSize ) & RK_SERIAL_SAVE_WITH_NAME( names );
};

void RK_CALL data_extractor::load( serialization::iarchive& A, unsigned int ) {
  closeExtractProcess();
  ReaKaux::unique_lock< ReaKaux::mutex > lock_here( access_mutex );
  shared_object::load( A, shared_object::getStaticObjectType()->TypeVersion() );
  unsigned int aColCount;
  A& RK_SERIAL_LOAD_WITH_ALIAS( "colCount", aColCount ) & RK_SERIAL_LOAD_WITH_NAME( flushSampleRate )
    & RK_SERIAL_LOAD_WITH_NAME( minBufferSize ) & RK_SERIAL_LOAD_WITH_NAME( names );
  colCount = aColCount;
  for( std::size_t i = 0; i < names.size(); ++i )
    named_indices[names[i]] = i;
  currentColumn = 0;
  currentNameCol = 0;
  values_rm = std::queue< double >();
  lock_here.unlock();
  reading_thread = ReaK::shared_ptr< ReaKaux::thread >( new ReaKaux::thread( extract_process( *this ) ) );
};
};
};
