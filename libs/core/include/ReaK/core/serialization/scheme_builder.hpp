/**
 * \file scheme_builder.hpp
 *
 * This library declares the class for a creating type schemes that represent the fields contained in the
 * serialization of a given type.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date November 2012
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

#ifndef REAK_SCHEME_BUILDER_HPP
#define REAK_SCHEME_BUILDER_HPP

#include "archiver.hpp"

#include <stack>
#include <string>
#include <utility>
#include <map>

namespace ReaK {

namespace serialization {

class serializable_obj_scheme;
class type_scheme;

std::map< std::string, shared_ptr< type_scheme > >& get_global_schemes();

/**
 * Protobuf scheme constructor.
 */
class scheme_builder : public oarchive {
private:
  std::stack< shared_ptr< serializable_obj_scheme > > field_stack;
  std::stack< std::pair< std::string, std::string > > value_name_stack;

protected:
  template < typename T >
  void save_primitive( const std::string& aName );

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

  virtual oarchive& RK_CALL save_int( std::ptrdiff_t i );

  virtual oarchive& RK_CALL save_int( const std::pair< std::string, std::ptrdiff_t >& i );

  virtual oarchive& RK_CALL save_unsigned_int( std::size_t u );

  virtual oarchive& RK_CALL save_unsigned_int( const std::pair< std::string, std::size_t >& u );

  virtual oarchive& RK_CALL save_float( float f );

  virtual oarchive& RK_CALL save_float( const std::pair< std::string, float >& f );

  virtual oarchive& RK_CALL save_double( double d );

  virtual oarchive& RK_CALL save_double( const std::pair< std::string, double >& d );

  virtual oarchive& RK_CALL save_bool( bool b );

  virtual oarchive& RK_CALL save_bool( const std::pair< std::string, bool >& b );

  virtual oarchive& RK_CALL save_string( const std::string& s );

  virtual oarchive& RK_CALL save_string( const std::pair< std::string, const std::string& >& s );

  virtual void RK_CALL signal_polymorphic_field( const std::string& aBaseTypeName, const unsigned int* aTypeID,
                                                 const std::string& aFieldName );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName );

  virtual void RK_CALL start_repeated_field( const std::string& aTypeName, const std::string& s );

  virtual void RK_CALL finish_repeated_field();

  virtual void RK_CALL start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2 );

  virtual void RK_CALL
    start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2, const std::string& s );

  virtual void RK_CALL finish_repeated_pair();

public:
  scheme_builder();
  virtual ~scheme_builder();
};


}; // serialization

}; // ReaK

#endif
