/**
 * \file protobuf_archiver.hpp
 *
 * This library declares the class for a protobuf archive to which an object hierarchy
 * can be serialized to and from.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_PROTOBUF_ARCHIVER_HPP
#define REAK_PROTOBUF_ARCHIVER_HPP

#include "archiver.hpp"

#include <iostream>
#include <stack>
#include <string>
#include <vector>
#include <map>


namespace ReaK {

namespace serialization {

/**
 * Protobuf input archive.
 */
class protobuf_iarchive : public iarchive {
private:
  shared_ptr< std::istream > file_stream;

protected:
  virtual iarchive& RK_CALL load_serializable_ptr( serializable_shared_pointer& Item );

  virtual iarchive& RK_CALL load_serializable_ptr( const std::pair< std::string, serializable_shared_pointer& >& Item );

  virtual iarchive& RK_CALL load_serializable( serializable& Item );

  virtual iarchive& RK_CALL load_serializable( const std::pair< std::string, serializable& >& Item );

  virtual iarchive& RK_CALL load_char( char& i );

  virtual iarchive& RK_CALL load_char( const std::pair< std::string, char& >& i );

  virtual iarchive& RK_CALL load_unsigned_char( unsigned char& u );

  virtual iarchive& RK_CALL load_unsigned_char( const std::pair< std::string, unsigned char& >& u );

  virtual iarchive& RK_CALL load_int( int& i );

  virtual iarchive& RK_CALL load_int( const std::pair< std::string, int& >& i );

  virtual void load_varint( unsigned int& u );

  virtual iarchive& RK_CALL load_unsigned_int( unsigned int& u );

  virtual iarchive& RK_CALL load_unsigned_int( const std::pair< std::string, unsigned int& >& u );

  virtual iarchive& RK_CALL load_float( float& f );

  virtual iarchive& RK_CALL load_float( const std::pair< std::string, float& >& f );

  virtual iarchive& RK_CALL load_double( double& d );

  virtual iarchive& RK_CALL load_double( const std::pair< std::string, double& >& d );

  virtual iarchive& RK_CALL load_bool( bool& b );

  virtual iarchive& RK_CALL load_bool( const std::pair< std::string, bool& >& b );

  virtual iarchive& RK_CALL load_string( std::string& s );

  virtual iarchive& RK_CALL load_string( const std::pair< std::string, std::string& >& s );

public:
  protobuf_iarchive( const std::string& FileName );
  protobuf_iarchive( std::istream& aStream );
  virtual ~protobuf_iarchive();
};

/**
 * Protobuf output archive.
 */
class protobuf_oarchive : public oarchive {
private:
  shared_ptr< std::ostream > file_stream;
  std::stack< unsigned int > field_IDs;
  std::stack< unsigned int > repeat_state;

protected:
  virtual oarchive& RK_CALL
    saveToNewArchive_impl( const serializable_shared_pointer& Item, const std::string& FileName );

  virtual oarchive& RK_CALL
    saveToNewArchiveNamed_impl( const std::pair< std::string, const serializable_shared_pointer& >& Item,
                                const std::string& FileName );

  virtual oarchive& RK_CALL save_serializable_ptr( const serializable_shared_pointer& Item );

  virtual oarchive& RK_CALL
    save_serializable_ptr( const std::pair< std::string, const serializable_shared_pointer& >& Item );

  virtual oarchive& RK_CALL save_serializable( const serializable& Item );

  virtual oarchive& RK_CALL save_serializable( const std::pair< std::string, const serializable& >& Item );

  virtual oarchive& RK_CALL save_char( char i );

  virtual oarchive& RK_CALL save_char( const std::pair< std::string, char >& i );

  virtual oarchive& RK_CALL save_unsigned_char( unsigned char u );

  virtual oarchive& RK_CALL save_unsigned_char( const std::pair< std::string, unsigned char >& u );

  virtual oarchive& RK_CALL save_int( int i );

  virtual oarchive& RK_CALL save_int( const std::pair< std::string, int >& i );

  void save_varint( unsigned int u );

  virtual oarchive& RK_CALL save_unsigned_int( unsigned int u );

  virtual oarchive& RK_CALL save_unsigned_int( const std::pair< std::string, unsigned int >& u );

  virtual oarchive& RK_CALL save_float( float f );

  virtual oarchive& RK_CALL save_float( const std::pair< std::string, float >& f );

  virtual oarchive& RK_CALL save_double( double d );

  virtual oarchive& RK_CALL save_double( const std::pair< std::string, double >& d );

  virtual oarchive& RK_CALL save_bool( bool b );

  virtual oarchive& RK_CALL save_bool( const std::pair< std::string, bool >& b );

  virtual oarchive& RK_CALL save_string( const std::string& s );

  virtual oarchive& RK_CALL save_string( const std::pair< std::string, const std::string& >& s );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName, const std::string& s );

  virtual void RK_CALL finish_repeated_field();

  virtual void RK_CALL start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2 );

  virtual void RK_CALL
    start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2, const std::string& s );

  virtual void RK_CALL finish_repeated_pair();


public:
  protobuf_oarchive( const std::string& FileName );
  protobuf_oarchive( std::ostream& aStream );
  virtual ~protobuf_oarchive();
};


/**
 * Protobuf scheme constructor.
 */
class protobuf_schemer : public oarchive {
private:
  shared_ptr< std::ostream > file_stream;
  std::stack< unsigned int > field_IDs;
  std::stack< unsigned int > repeat_state;

  std::vector< std::string > schemes;
  std::map< std::string, std::size_t > scheme_map;

protected:
  unsigned int get_chunk_hdr();

  virtual oarchive& RK_CALL
    saveToNewArchive_impl( const serializable_shared_pointer& Item, const std::string& FileName );

  virtual oarchive& RK_CALL
    saveToNewArchiveNamed_impl( const std::pair< std::string, const serializable_shared_pointer& >& Item,
                                const std::string& FileName );

  virtual oarchive& RK_CALL save_serializable_ptr( const serializable_shared_pointer& Item );

  virtual oarchive& RK_CALL
    save_serializable_ptr( const std::pair< std::string, const serializable_shared_pointer& >& Item );

  virtual oarchive& RK_CALL save_serializable( const serializable& Item );

  virtual oarchive& RK_CALL save_serializable( const std::pair< std::string, const serializable& >& Item );

  virtual oarchive& RK_CALL save_char( char i );

  virtual oarchive& RK_CALL save_char( const std::pair< std::string, char >& i );

  virtual oarchive& RK_CALL save_unsigned_char( unsigned char u );

  virtual oarchive& RK_CALL save_unsigned_char( const std::pair< std::string, unsigned char >& u );

  virtual oarchive& RK_CALL save_int( int i );

  virtual oarchive& RK_CALL save_int( const std::pair< std::string, int >& i );

  virtual oarchive& RK_CALL save_unsigned_int( unsigned int u );

  virtual oarchive& RK_CALL save_unsigned_int( const std::pair< std::string, unsigned int >& u );

  virtual oarchive& RK_CALL save_float( float f );

  virtual oarchive& RK_CALL save_float( const std::pair< std::string, float >& f );

  virtual oarchive& RK_CALL save_double( double d );

  virtual oarchive& RK_CALL save_double( const std::pair< std::string, double >& d );

  virtual oarchive& RK_CALL save_bool( bool b );

  virtual oarchive& RK_CALL save_bool( const std::pair< std::string, bool >& b );

  virtual oarchive& RK_CALL save_string( const std::string& s );

  virtual oarchive& RK_CALL save_string( const std::pair< std::string, const std::string& >& s );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName, const std::string& s );

  virtual void RK_CALL finish_repeated_field();

  virtual void RK_CALL start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2 );

  virtual void RK_CALL
    start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2, const std::string& s );

  virtual void RK_CALL finish_repeated_pair();

public:
  void print_schemes( std::ostream& aStream );

  protobuf_schemer();
  virtual ~protobuf_schemer();
};


}; // serialization

}; // ReaK

#endif
