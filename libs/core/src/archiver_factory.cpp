
/*
 *    Copyright 2013 Sven Mikael Persson
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

#include <ReaK/core/serialization/archiver_factory.hpp>

#include <ReaK/core/serialization/bin_archiver.hpp>
#include <ReaK/core/serialization/xml_archiver.hpp>
#include <ReaK/core/serialization/protobuf_archiver.hpp>

#include <ios>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Serialization */
namespace serialization {


shared_ptr< iarchive > open_iarchive( const std::string& aFileName ) {

  std::string fileExt = aFileName.substr( aFileName.find_last_of( '.' ) + 1 );

  if( ( fileExt == "rkx" ) || ( fileExt == "xml" ) )
    return shared_ptr< iarchive >( new xml_iarchive( aFileName ) );
  else if( ( fileExt == "rkb" ) || ( fileExt == "bin" ) )
    return shared_ptr< iarchive >( new bin_iarchive( aFileName ) );
  else if( fileExt == "pbuf" )
    return shared_ptr< iarchive >( new protobuf_iarchive( aFileName ) );
  else
    throw std::ios_base::failure( "Sorry, this file-type is not supported!" );
};

shared_ptr< oarchive > open_oarchive( const std::string& aFileName ) {

  std::string fileExt = aFileName.substr( aFileName.find_last_of( '.' ) + 1 );

  if( ( fileExt == "rkx" ) || ( fileExt == "xml" ) )
    return shared_ptr< oarchive >( new xml_oarchive( aFileName ) );
  else if( ( fileExt == "rkb" ) || ( fileExt == "bin" ) )
    return shared_ptr< oarchive >( new bin_oarchive( aFileName ) );
  else if( fileExt == "pbuf" )
    return shared_ptr< oarchive >( new protobuf_oarchive( aFileName ) );
  else
    throw std::ios_base::failure( "Sorry, this file-type is not supported!" );
};
};
};
