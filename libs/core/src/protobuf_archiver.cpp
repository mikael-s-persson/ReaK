
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

#include <ReaK/core/serialization/protobuf_archiver.hpp>


#include <ReaK/core/base/shared_object.hpp>
#include <ReaK/core/rtti/rtti.hpp>

#include <ReaK/core/serialization/archiving_exceptions.hpp>

#include <string>
#include <fstream>
#include <sstream>


#include <cstdint>


namespace ReaK {

namespace serialization {


namespace {

union float_to_ulong {
  float f;
  std::uint32_t ui32;
};

union double_to_ulong {
  double d;
  std::uint64_t ui64;
  std::uint32_t ui32[2];
};

union llong_to_ulong {
  std::int64_t i64;
  std::uint64_t ui64;
  std::uint32_t ui32[2];
};

union ulong_to_uword {
  std::uint32_t ui32;
  std::uint16_t ui16[2];
  std::uint8_t ui8[4];
};


void le2h_1ui32( std::uint32_t& value ) {
#if RK_BYTE_ORDER == RK_ORDER_BIG_ENDIAN
  ulong_to_uword tmp;
  tmp.ui32 = value;
  std::uint8_t tmp_b = tmp.ui8[0];
  tmp.ui8[0] = tmp.ui8[3];
  tmp.ui8[3] = tmp_b;
  tmp_b = tmp.ui8[1];
  tmp.ui8[1] = tmp.ui8[2];
  tmp.ui8[2] = tmp_b;
  value = tmp.ui32;
#elif RK_BYTE_ORDER == RK_ORDER_PDP_ENDIAN
  ulong_to_uword tmp;
  tmp.ui32 = value;
  std::uint16_t tmp2 = tmp.ui16[0];
  tmp.ui16[0] = tmp.ui16[1];
  tmp.ui16[1] = tmp2;
  value = tmp.ui32;
#endif
};

template < typename UnionT >
void le2h_2ui32( UnionT& value ) {
#if RK_BYTE_ORDER == RK_ORDER_BIG_ENDIAN
  le2h_1ui32( value.ui32[0] );
  le2h_1ui32( value.ui32[1] );
  std::uint32_t tmp = value.ui32[0];
  value.ui32[0] = value.ui32[1];
  value.ui32[1] = tmp;
#endif
};
};


static const char bad_stream_msg[]
  = "Protobuf archive could not be loaded due to an unexpected failure or bad state of the input stream!";
static const char unexpected_eof_msg[] = "Protobuf archive is corrupt! Unexpectingly reached end-of-file!";
static const char corrupt_header_msg[] = "Protobuf archive has a corrupt header!";
static const char unknown_version_msg[] = "Protobuf archive is of an unknown file version!";

static const char bad_out_stream_msg[]
  = "Protobuf archive could not be saved due to an unexpected failure or bad state of the output stream!";


protobuf_iarchive::protobuf_iarchive( const std::string& FileName ) {

  file_stream = shared_ptr< std::istream >( new std::ifstream( FileName.c_str(), std::ios::binary | std::ios::in ) );

  std::string header;
  protobuf_iarchive::load_string( header );
  unsigned int version;
  protobuf_iarchive::load_unsigned_int( version );

  if( !( header == "reak_serialization::protobuf_archive" ) )
    throw std::ios_base::failure( corrupt_header_msg );
  if( version != 2 )
    throw std::ios_base::failure( unknown_version_msg );
};


protobuf_iarchive::protobuf_iarchive( std::istream& aStream ) {

  file_stream = shared_ptr< std::istream >( &aStream, null_deleter() );

  std::string header;
  protobuf_iarchive::load_string( header );
  unsigned int version;
  protobuf_iarchive::load_unsigned_int( version );

  if( !( header == "reak_serialization::protobuf_archive" ) )
    throw std::ios_base::failure( corrupt_header_msg );
  if( version != 2 )
    throw std::ios_base::failure( unknown_version_msg );
};

protobuf_iarchive::~protobuf_iarchive(){};

iarchive& RK_CALL protobuf_iarchive::load_serializable_ptr( serializable_shared_pointer& Item ) {
  archive_object_header hdr;
  Item = serializable_shared_pointer();

  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 2 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading serializable object pointer should "
          "have wire-type 2. Got chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec
       << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };

  protobuf_iarchive::load_varint( hdr.size );
  std::streampos start_pos = file_stream->tellg();
  if( start_pos < 0 )
    throw std::ios_base::failure( bad_stream_msg );

  std::vector< unsigned int > typeIDvect;
  unsigned int i;
  do {
    *this >> i;
    typeIDvect.push_back( i );
  } while( i != 0 );

  *this >> hdr.type_version >> hdr.object_ID >> hdr.is_external;

  if( hdr.object_ID == 0 ) {
    // Item already null.
    std::streampos end_pos = file_stream->tellg();
    if( hdr.size + start_pos != end_pos )
      if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
        throw std::ios_base::failure( bad_stream_msg );
    return *this;
  };
  if( ( hdr.object_ID < mObjRegistry.size() ) && ( mObjRegistry[hdr.object_ID] ) ) {
    Item = mObjRegistry[hdr.object_ID];
    std::streampos end_pos = file_stream->tellg();
    if( hdr.size + start_pos != end_pos )
      if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
        throw std::ios_base::failure( bad_stream_msg );
    return *this;
  };

  if( hdr.is_external ) {
    std::string ext_filename;
    *this >> ext_filename;
    std::streampos end_pos = file_stream->tellg();
    if( hdr.size + start_pos != end_pos )
      if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
        throw std::ios_base::failure( bad_stream_msg );

    protobuf_iarchive a( ext_filename ); // if this throws, let it propagate up (no point catching and throwing).
    a >> Item;

    return *this;
  };

  // Find the class in question in the repository.
  rtti::so_type* p = rtti::so_type_repo::getInstance().findType( &( typeIDvect[0] ) );
  if( ( !p ) || ( p->TypeVersion() < hdr.type_version ) ) {
    std::streampos end_pos = file_stream->tellg();
    if( hdr.size + start_pos != end_pos )
      if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
        throw std::ios_base::failure( bad_stream_msg );
    throw unsupported_type( unsupported_type::not_found_in_repo, &( typeIDvect[0] ) );
  };
  ReaK::shared_ptr< shared_object > po( p->CreateObject() );
  if( !po ) {
    std::streampos end_pos = file_stream->tellg();
    if( hdr.size + start_pos != end_pos )
      if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
        throw std::ios_base::failure( bad_stream_msg );
    throw unsupported_type( unsupported_type::could_not_create, &( typeIDvect[0] ) );
  };

  Item = po;
  if( hdr.object_ID < mObjRegistry.size() )
    mObjRegistry[hdr.object_ID] = Item;
  else if( hdr.object_ID == mObjRegistry.size() )
    mObjRegistry.push_back( Item ); // in theory, only this condition should occur
  else if( hdr.object_ID > mObjRegistry.size() ) {
    mObjRegistry.resize( hdr.object_ID + 1 );
    mObjRegistry[hdr.object_ID] = Item;
  };

  Item->load( *this, hdr.type_version );

  std::streampos end_pos = file_stream->tellg();
  if( hdr.size + start_pos != end_pos )
    if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
      throw std::ios_base::failure( bad_stream_msg );

  return *this;
};

iarchive& RK_CALL
  protobuf_iarchive::load_serializable_ptr( const std::pair< std::string, serializable_shared_pointer& >& Item ) {
  return protobuf_iarchive::load_serializable_ptr( Item.second );
};

iarchive& RK_CALL protobuf_iarchive::load_serializable( serializable& Item ) {
  archive_object_header hdr;

  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 2 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading serializable object should have "
          "wire-type 2. Got chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec
       << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };

  protobuf_iarchive::load_varint( hdr.size );
  std::streampos start_pos = file_stream->tellg();

  std::vector< unsigned int > typeIDvect;
  unsigned int i;
  do {
    *this >> i;
    typeIDvect.push_back( i );
  } while( i != 0 );

  *this >> hdr.type_version;

  Item.load( *this, hdr.type_version );
  std::streampos end_pos = file_stream->tellg();

  if( hdr.size + start_pos != end_pos )
    if( !file_stream->seekg( start_pos + std::streampos( hdr.size ) ) )
      throw std::ios_base::failure( bad_stream_msg );

  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_serializable( const std::pair< std::string, serializable& >& Item ) {
  return protobuf_iarchive::load_serializable( Item.second );
};

iarchive& RK_CALL protobuf_iarchive::load_char( char& i ) {
  int il = 0;
  protobuf_iarchive::load_int( il );
  i = il;
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_char( const std::pair< std::string, char& >& i ) {
  return protobuf_iarchive::load_char( i.second );
};

iarchive& RK_CALL protobuf_iarchive::load_unsigned_char( unsigned char& u ) {
  unsigned int ul = 0;
  protobuf_iarchive::load_unsigned_int( ul );
  u = ul;
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_unsigned_char( const std::pair< std::string, unsigned char& >& u ) {
  return protobuf_iarchive::load_unsigned_char( u.second );
};

iarchive& RK_CALL protobuf_iarchive::load_int( int& i ) {
  unsigned int u = 0;
  protobuf_iarchive::load_unsigned_int( u );
  i = ( u >> 1 ) ^ ( -static_cast< int >( u & 1 ) );
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_int( const std::pair< std::string, int& >& i ) {
  return protobuf_iarchive::load_int( i.second );
};

void protobuf_iarchive::load_varint( unsigned int& u ) {
  u = 0;
  std::uint8_t tmp = 0;
  if( !file_stream->read( reinterpret_cast< char* >( &tmp ), 1 ) )
    throw std::ios_base::failure( unexpected_eof_msg );
  std::uint8_t shifts = 0;
  u = tmp & 0x7F;
  while( tmp & 0x80 ) {
    if( !file_stream->read( reinterpret_cast< char* >( &tmp ), 1 ) )
      throw std::ios_base::failure( unexpected_eof_msg );
    shifts += 7;
    u |= ( tmp & 0x7F ) << shifts;
  };
};

iarchive& RK_CALL protobuf_iarchive::load_unsigned_int( unsigned int& u ) {
  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 0 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading varint should have wire-type 0. Got "
          "chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };

  protobuf_iarchive::load_varint( u );
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_unsigned_int( const std::pair< std::string, unsigned int& >& u ) {
  return protobuf_iarchive::load_unsigned_int( u.second );
};

iarchive& RK_CALL protobuf_iarchive::load_float( float& f ) {
  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 5 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading float should have wire-type 5. Got "
          "chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };
  float_to_ulong tmp;
  if( !file_stream->read( reinterpret_cast< char* >( &tmp ), sizeof( float ) ) )
    throw std::ios_base::failure( unexpected_eof_msg );
  le2h_1ui32( tmp.ui32 );
  f = tmp.f;
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_float( const std::pair< std::string, float& >& f ) {
  return protobuf_iarchive::load_float( f.second );
};

iarchive& RK_CALL protobuf_iarchive::load_double( double& d ) {
  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 1 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading double should have wire-type 1. Got "
          "chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };
  double_to_ulong tmp;
  if( !file_stream->read( reinterpret_cast< char* >( &tmp ), sizeof( double ) ) )
    throw std::ios_base::failure( unexpected_eof_msg );
  le2h_2ui32( tmp );
  d = tmp.d;
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_double( const std::pair< std::string, double& >& d ) {
  return protobuf_iarchive::load_double( d.second );
};

iarchive& RK_CALL protobuf_iarchive::load_bool( bool& b ) {
  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 0 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading bool should have wire-type 0. Got "
          "chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };
  char tmp = 0;
  if( !file_stream->read( &tmp, 1 ) )
    throw std::ios_base::failure( unexpected_eof_msg );
  b = ( tmp ? true : false );
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_bool( const std::pair< std::string, bool& >& b ) {
  return protobuf_iarchive::load_bool( b.second );
};

iarchive& RK_CALL protobuf_iarchive::load_string( std::string& s ) {
  unsigned int chunk_hdr;
  protobuf_iarchive::load_varint( chunk_hdr );
  if( ( chunk_hdr & 0x07 ) != 2 ) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg( std::ios_base::beg );
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! Loading string should have wire-type 2. Got "
          "chunk-ID: " << std::hex << chunk_hdr << " at offset " << std::dec << ( current_pos - start_pos ) << ".";
    throw std::ios_base::failure( ss.str() );
  };
  unsigned int u;
  protobuf_iarchive::load_varint( u );
  s.resize( u );
  if( !file_stream->read( &s[0], u ) )
    throw std::ios_base::failure( unexpected_eof_msg );
  return *this;
};

iarchive& RK_CALL protobuf_iarchive::load_string( const std::pair< std::string, std::string& >& s ) {
  return protobuf_iarchive::load_string( s.second );
};


protobuf_oarchive::protobuf_oarchive( const std::string& FileName ) {
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  file_stream = shared_ptr< std::ostream >( new std::ofstream( FileName.c_str(), std::ios::binary | std::ios::out ) );

  std::string header = "reak_serialization::protobuf_archive";
  protobuf_oarchive::save_string( header );
  unsigned int version = 2;
  protobuf_oarchive::save_unsigned_int( version );
};

protobuf_oarchive::protobuf_oarchive( std::ostream& aStream ) {
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  file_stream = shared_ptr< std::ostream >( &aStream, null_deleter() );

  std::string header = "reak_serialization::protobuf_archive";
  protobuf_oarchive::save_string( header );
  unsigned int version = 2;
  protobuf_oarchive::save_unsigned_int( version );
};

protobuf_oarchive::~protobuf_oarchive(){};

oarchive& RK_CALL
  protobuf_oarchive::saveToNewArchive_impl( const serializable_shared_pointer& Item, const std::string& FileName ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 2; // wire-type 2: length-delimited.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );

  shared_ptr< std::ostream > tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  shared_ptr< std::stringstream > str_stream( new std::stringstream( std::ios::binary | std::ios::out ) );
  file_stream = str_stream;
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  archive_object_header hdr;
  bool already_saved( false );
  const unsigned int* type_ID = NULL;

  if( Item ) {
    std::map< serializable_shared_pointer, unsigned int >::const_iterator it = mObjRegMap.find( Item );

    if( it != mObjRegMap.end() ) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back( Item );
      mObjRegMap[Item] = hdr.object_ID;
    };

    rtti::so_type* obj_type = Item->getObjectType();
    type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = true;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    already_saved = true;
  };

  protobuf_oarchive::start_repeated_field( "unsigned int" );
  while( ( type_ID ) && ( *type_ID ) ) {
    protobuf_oarchive::save_unsigned_int( *type_ID );
    ++type_ID;
  };
  protobuf_oarchive::save_unsigned_int( 0 );
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int( hdr.type_version );
  protobuf_oarchive::save_unsigned_int( hdr.object_ID );
  protobuf_oarchive::save_bool( hdr.is_external );

  if( !already_saved ) {
    protobuf_oarchive::save_string( FileName );

    protobuf_oarchive a( FileName );
    a << Item;
  };

  file_stream = tmp_str_ptr;
  field_IDs.pop();
  repeat_state.pop();

  str_stream->seekp( 0, std::ios::end );
  hdr.size = static_cast< unsigned int >( str_stream->tellp() );
  str_stream->clear();

  protobuf_oarchive::save_varint( hdr.size );
  if( !( *file_stream << str_stream->str() ) )
    throw std::ios_base::failure( bad_out_stream_msg );

  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.

  return *this;
};

oarchive& RK_CALL protobuf_oarchive::saveToNewArchiveNamed_impl(
  const std::pair< std::string, const serializable_shared_pointer& >& Item, const std::string& FileName ) {
  return protobuf_oarchive::saveToNewArchive_impl( Item.second, FileName );
};


oarchive& RK_CALL protobuf_oarchive::save_serializable_ptr( const serializable_shared_pointer& Item ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 2; // wire-type 2: length-delimited.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );

  shared_ptr< std::ostream > tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  shared_ptr< std::stringstream > str_stream( new std::stringstream( std::ios::binary | std::ios::out ) );
  file_stream = str_stream;
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  archive_object_header hdr;
  bool already_saved( false );
  const unsigned int* type_ID = NULL;

  if( Item ) {
    std::map< serializable_shared_pointer, unsigned int >::const_iterator it = mObjRegMap.find( Item );

    if( it != mObjRegMap.end() ) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back( Item );
      mObjRegMap[Item] = hdr.object_ID;
    };

    rtti::so_type* obj_type = Item->getObjectType();
    type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = false;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    already_saved = true;
  };

  protobuf_oarchive::start_repeated_field( "unsigned int" );
  while( ( type_ID ) && ( *type_ID ) ) {
    protobuf_oarchive::save_unsigned_int( *type_ID );
    ++type_ID;
  };
  protobuf_oarchive::save_unsigned_int( 0 );
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int( hdr.type_version );
  protobuf_oarchive::save_unsigned_int( hdr.object_ID );
  protobuf_oarchive::save_bool( hdr.is_external );

  if( !already_saved ) {
    Item->save( *this, hdr.type_version );
  };

  file_stream = tmp_str_ptr;
  field_IDs.pop();
  repeat_state.pop();

  str_stream->seekp( 0, std::ios::end );
  hdr.size = static_cast< unsigned int >( str_stream->tellp() );
  str_stream->clear();

  protobuf_oarchive::save_varint( hdr.size );
  if( !( *file_stream << str_stream->str() ) )
    throw std::ios_base::failure( bad_out_stream_msg );

  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.

  return *this;
};


oarchive& RK_CALL
  protobuf_oarchive::save_serializable_ptr( const std::pair< std::string, const serializable_shared_pointer& >& Item ) {
  return protobuf_oarchive::save_serializable_ptr( Item.second );
};

oarchive& RK_CALL protobuf_oarchive::save_serializable( const serializable& Item ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 2; // wire-type 2: length-delimited.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );

  shared_ptr< std::ostream > tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  shared_ptr< std::stringstream > str_stream( new std::stringstream( std::ios::binary | std::ios::out ) );
  file_stream = str_stream;
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  archive_object_header hdr;
  const unsigned int* type_ID = Item.getObjectType()->TypeID_begin();
  hdr.type_version = Item.getObjectType()->TypeVersion();
  hdr.object_ID = 0;
  hdr.is_external = false;
  hdr.size = 0;

  protobuf_oarchive::start_repeated_field( "unsigned int" );
  while( *type_ID ) {
    protobuf_oarchive::save_unsigned_int( *type_ID );
    ++type_ID;
  };
  protobuf_oarchive::save_unsigned_int( 0 );
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int( hdr.type_version );

  Item.save( *this, hdr.type_version );

  file_stream = tmp_str_ptr;
  field_IDs.pop();
  repeat_state.pop();

  str_stream->seekp( 0, std::ios::end );
  hdr.size = static_cast< unsigned int >( str_stream->tellp() );
  str_stream->clear();

  protobuf_oarchive::save_varint( hdr.size );
  if( !( *file_stream << str_stream->str() ) )
    throw std::ios_base::failure( bad_out_stream_msg );

  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.

  return *this;
};

oarchive& RK_CALL protobuf_oarchive::save_serializable( const std::pair< std::string, const serializable& >& Item ) {
  return protobuf_oarchive::save_serializable( Item.second );
};

oarchive& RK_CALL protobuf_oarchive::save_char( char i ) { return protobuf_oarchive::save_int( i ); };

oarchive& RK_CALL protobuf_oarchive::save_char( const std::pair< std::string, char >& i ) {
  return protobuf_oarchive::save_char( i.second );
};

oarchive& RK_CALL protobuf_oarchive::save_unsigned_char( unsigned char u ) {
  return protobuf_oarchive::save_unsigned_int( u );
};

oarchive& RK_CALL protobuf_oarchive::save_unsigned_char( const std::pair< std::string, unsigned char >& u ) {
  return protobuf_oarchive::save_unsigned_char( u.second );
};

oarchive& RK_CALL protobuf_oarchive::save_int( int i ) {
  protobuf_oarchive::save_unsigned_int( ( i << 1 ) ^ ( i >> ( sizeof( int ) * 8 - 1 ) ) );
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_int( const std::pair< std::string, int >& i ) {
  return protobuf_oarchive::save_int( i.second );
};

void protobuf_oarchive::save_varint( unsigned int u ) {
  std::uint8_t buf[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // 80-bits, supports at most a 64-bit varint.
  std::uint8_t* pbuf = buf;
  *pbuf = ( u & 0x7F );
  u >>= 7;
  while( u ) {
    *pbuf |= 0x80; // set first msb because there is more to come.
    pbuf++;
    *pbuf = ( u & 0x7F );
    u >>= 7;
  };
  if( !file_stream->write( reinterpret_cast< char* >( buf ), pbuf - buf + 1 ) )
    throw std::ios_base::failure( bad_out_stream_msg );
};

oarchive& RK_CALL protobuf_oarchive::save_unsigned_int( unsigned int u ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ); // wire-type 0: varint.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );
  protobuf_oarchive::save_varint( u );
  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_unsigned_int( const std::pair< std::string, unsigned int >& u ) {
  return protobuf_oarchive::save_unsigned_int( u.second );
};


oarchive& RK_CALL protobuf_oarchive::save_float( float f ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 5; // wire-type 5: 32-bit.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );
  float_to_ulong tmp = {f};
  le2h_1ui32( tmp.ui32 );
  if( !file_stream->write( reinterpret_cast< char* >( &tmp ), sizeof( float_to_ulong ) ) )
    throw std::ios_base::failure( bad_out_stream_msg );
  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_float( const std::pair< std::string, float >& f ) {
  return protobuf_oarchive::save_float( f.second );
};


oarchive& RK_CALL protobuf_oarchive::save_double( double d ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 1; // wire-type 1: 64-bit.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );
  double_to_ulong tmp = {d};
  le2h_2ui32( tmp );
  if( !file_stream->write( reinterpret_cast< char* >( &tmp ), sizeof( double_to_ulong ) ) )
    throw std::ios_base::failure( bad_out_stream_msg );
  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_double( const std::pair< std::string, double >& d ) {
  return protobuf_oarchive::save_double( d.second );
};


oarchive& RK_CALL protobuf_oarchive::save_bool( bool b ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ); // wire-type 0: varint.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );
  char tmp = 0;
  if( b )
    tmp = 1;
  if( !file_stream->put( tmp ) )
    throw std::ios_base::failure( bad_out_stream_msg );
  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_bool( const std::pair< std::string, bool >& b ) {
  return protobuf_oarchive::save_bool( b.second );
};


oarchive& RK_CALL protobuf_oarchive::save_string( const std::string& s ) {
  unsigned int chunk_hdr = ( field_IDs.top() << 3 ) | 2; // wire-type 2: length-delimited.
  if( repeat_state.top() & 0x02 ) {
    if( repeat_state.top() & 0x04 )
      chunk_hdr += 8;
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint( chunk_hdr );
  unsigned int u = s.length();
  protobuf_oarchive::save_varint( u );
  if( !file_stream->write( s.data(), s.length() ) )
    throw std::ios_base::failure( bad_out_stream_msg );
  if( !repeat_state.top() )
    ++( field_IDs.top() ); // increment the field ID.
  return *this;
};


oarchive& RK_CALL protobuf_oarchive::save_string( const std::pair< std::string, const std::string& >& s ) {
  return protobuf_oarchive::save_string( s.second );
};


void RK_CALL protobuf_oarchive::start_repeated_field( const std::string& ) { repeat_state.push( 1 ); };

void RK_CALL protobuf_oarchive::start_repeated_field( const std::string&, const std::string& ) {
  repeat_state.push( 1 );
};

void RK_CALL protobuf_oarchive::finish_repeated_field() {
  repeat_state.pop();
  ++( field_IDs.top() );
};

void RK_CALL protobuf_oarchive::start_repeated_pair( const std::string&, const std::string& ) {
  repeat_state.push( 3 );
};

void RK_CALL protobuf_oarchive::start_repeated_pair( const std::string&, const std::string&, const std::string& ) {
  repeat_state.push( 3 );
};

void RK_CALL protobuf_oarchive::finish_repeated_pair() {
  repeat_state.pop();
  field_IDs.top() += 2;
};


protobuf_schemer::protobuf_schemer() {
  field_IDs.push( 0 );
  repeat_state.push( 0 );

  file_stream = shared_ptr< std::ostream >( new std::stringstream() );
};

protobuf_schemer::~protobuf_schemer(){};


void protobuf_schemer::print_schemes( std::ostream& aStream ) {
  for( std::size_t i = 0; i < schemes.size(); ++i )
    aStream << schemes[i] << std::endl << std::endl;
};


unsigned int protobuf_schemer::get_chunk_hdr() {
  unsigned int chunk_hdr = 0;
  if( repeat_state.top() & 0x08 )
    return ~chunk_hdr;
  else
    chunk_hdr = field_IDs.top();
  return chunk_hdr;
};


oarchive& RK_CALL
  protobuf_schemer::saveToNewArchive_impl( const serializable_shared_pointer& Item, const std::string& FileName ) {
  return protobuf_schemer::save_serializable_ptr(
    std::pair< std::string, const serializable_shared_pointer& >( "item", Item ) );
};

oarchive& RK_CALL protobuf_schemer::saveToNewArchiveNamed_impl(
  const std::pair< std::string, const serializable_shared_pointer& >& Item, const std::string& FileName ) {
  return protobuf_schemer::save_serializable_ptr( Item );
};


oarchive& RK_CALL protobuf_schemer::save_serializable_ptr( const serializable_shared_pointer& Item ) {
  return protobuf_schemer::save_serializable_ptr(
    std::pair< std::string, const serializable_shared_pointer& >( "item_ptr", Item ) );
};


oarchive& RK_CALL
  protobuf_schemer::save_serializable_ptr( const std::pair< std::string, const serializable_shared_pointer& >& Item ) {
  if( !Item.second )
    return *this;

  unsigned int chunk_hdr = get_chunk_hdr();

#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  constexpr auto tname = rtti::get_type_id< serializable_shared_pointer >::type_name;
  std::string aObjTypeName = tname.to_string();
#else
  std::string aObjTypeName = rtti::get_type_id< serializable_shared_pointer >::type_name();
#endif
  aObjTypeName += "<" + Item.second->getObjectType()->TypeName() + ">";
  std::map< serializable_shared_pointer, unsigned int >::const_iterator it = mObjRegMap.find( Item.second );

  if( it == mObjRegMap.end() ) {
    shared_ptr< std::ostream > tmp_str_ptr = file_stream;

    shared_ptr< std::stringstream > str_stream( new std::stringstream( std::ios::binary | std::ios::out ) );
    file_stream = str_stream;
    field_IDs.push( 4 );
    repeat_state.push( 0 );

    *file_stream << "message " << aObjTypeName << " {" << std::endl << "  repeated uint32 type_ID = 0;" << std::endl
                 << "  required uint32 version = 1;" << std::endl << "  required luid32 object_ID = 2;" << std::endl
                 << "  required bool is_external = 3;" << std::endl;

    mObjRegistry.push_back( Item.second );
    mObjRegMap[Item.second] = mObjRegistry.size() - 1;

    Item.second->save( *this, Item.second->getObjectType()->TypeVersion() );
    *file_stream << "}" << std::endl;

    file_stream = tmp_str_ptr;
    field_IDs.pop();
    repeat_state.pop();

    std::map< std::string, std::size_t >::iterator itm = scheme_map.find( aObjTypeName );
    if( itm == scheme_map.end() ) {
      schemes.push_back( str_stream->str() );
      scheme_map[aObjTypeName] = schemes.size() - 1;
    } else {
      std::string s_tmp = str_stream->str();
      if( schemes[itm->second].length() < s_tmp.length() )
        schemes[itm->second] = s_tmp;
    };
  };

  if( ~chunk_hdr == 0 )
    return *this;

  *file_stream << "  required " << aObjTypeName << " " << Item.first << " = " << chunk_hdr << ";" << std::endl;

  field_IDs.top() += 1;

  return *this;
};

oarchive& RK_CALL protobuf_schemer::save_serializable( const serializable& Item ) {
  return protobuf_schemer::save_serializable( std::pair< std::string, const serializable& >( "item", Item ) );
};

oarchive& RK_CALL protobuf_schemer::save_serializable( const std::pair< std::string, const serializable& >& Item ) {
  unsigned int chunk_hdr = get_chunk_hdr();

  shared_ptr< std::ostream > tmp_str_ptr = file_stream;

  shared_ptr< std::stringstream > str_stream( new std::stringstream() );
  file_stream = str_stream;
  field_IDs.push( 2 );
  repeat_state.push( 0 );

  *file_stream << "message " << Item.second.getObjectType()->TypeName() << " {" << std::endl
               << "  repeated uint32 type_ID = 0;" << std::endl << "  required uint32 version = 1;" << std::endl;

  Item.second.save( *this, Item.second.getObjectType()->TypeVersion() );

  *file_stream << "}" << std::endl;

  file_stream = tmp_str_ptr;
  field_IDs.pop();
  repeat_state.pop();

  std::map< std::string, std::size_t >::iterator itm = scheme_map.find( Item.second.getObjectType()->TypeName() );
  if( itm == scheme_map.end() ) {
    schemes.push_back( str_stream->str() );
    scheme_map[Item.second.getObjectType()->TypeName()] = schemes.size() - 1;
  } else {
    std::string s_tmp = str_stream->str();
    if( schemes[itm->second].length() < s_tmp.length() )
      schemes[itm->second] = s_tmp;
  };

  if( ~chunk_hdr == 0 )
    return *this;

  *file_stream << "  required " << Item.second.getObjectType()->TypeName() << " " << Item.first << " = " << chunk_hdr
               << ";" << std::endl;

  field_IDs.top() += 1;

  return *this;
};

oarchive& RK_CALL protobuf_schemer::save_char( char i ) {
  return protobuf_schemer::save_char( std::pair< std::string, char >( "i", i ) );
};

oarchive& RK_CALL protobuf_schemer::save_char( const std::pair< std::string, char >& i ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required int32 " << i.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};

oarchive& RK_CALL protobuf_schemer::save_unsigned_char( unsigned char u ) {
  return protobuf_schemer::save_unsigned_char( std::pair< std::string, unsigned char >( "u", u ) );
};

oarchive& RK_CALL protobuf_schemer::save_unsigned_char( const std::pair< std::string, unsigned char >& u ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required uint32 " << u.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};

oarchive& RK_CALL protobuf_schemer::save_int( int i ) {
  return protobuf_schemer::save_int( std::pair< std::string, int >( "i", i ) );
};


oarchive& RK_CALL protobuf_schemer::save_int( const std::pair< std::string, int >& i ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required int32 " << i.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};

oarchive& RK_CALL protobuf_schemer::save_unsigned_int( unsigned int u ) {
  return protobuf_schemer::save_unsigned_int( std::pair< std::string, unsigned int >( "u", u ) );
};


oarchive& RK_CALL protobuf_schemer::save_unsigned_int( const std::pair< std::string, unsigned int >& u ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required uint32 " << u.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};


oarchive& RK_CALL protobuf_schemer::save_float( float f ) {
  return protobuf_schemer::save_float( std::pair< std::string, float >( "f", f ) );
};


oarchive& RK_CALL protobuf_schemer::save_float( const std::pair< std::string, float >& f ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required float " << f.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};


oarchive& RK_CALL protobuf_schemer::save_double( double d ) {
  return protobuf_schemer::save_double( std::pair< std::string, double >( "d", d ) );
};


oarchive& RK_CALL protobuf_schemer::save_double( const std::pair< std::string, double >& d ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required double " << d.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};


oarchive& RK_CALL protobuf_schemer::save_bool( bool b ) {
  return protobuf_schemer::save_bool( std::pair< std::string, bool >( "b", b ) );
};


oarchive& RK_CALL protobuf_schemer::save_bool( const std::pair< std::string, bool >& b ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required bool " << b.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};


oarchive& RK_CALL protobuf_schemer::save_string( const std::string& s ) {
  return protobuf_schemer::save_string( std::pair< std::string, const std::string& >( "str", s ) );
};


oarchive& RK_CALL protobuf_schemer::save_string( const std::pair< std::string, const std::string& >& s ) {
  unsigned int chunk_hdr = get_chunk_hdr();
  if( ~chunk_hdr == 0 )
    return *this;
  *file_stream << "  required string " << s.first << " = " << chunk_hdr << ";" << std::endl;
  field_IDs.top() += 1;
  return *this;
};


void RK_CALL protobuf_schemer::start_repeated_field( const std::string& aTypeName ) {
  repeat_state.push( 9 );
  *file_stream << "  repeated " << aTypeName << " value = " << field_IDs.top() << ";" << std::endl;
};

void RK_CALL protobuf_schemer::start_repeated_field( const std::string& aTypeName, const std::string& aName ) {
  repeat_state.push( 9 );
  *file_stream << "  repeated " << aTypeName << " " << aName << " = " << field_IDs.top() << ";" << std::endl;
};

void RK_CALL protobuf_schemer::finish_repeated_field() {
  repeat_state.pop();
  field_IDs.top() += 1;
};

void RK_CALL protobuf_schemer::start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2 ) {
  repeat_state.push( 11 );
  *file_stream << "  repeated " << aTypeName1 << " map_key = " << field_IDs.top() << ";" << std::endl << "  repeated "
               << aTypeName2 << " map_value = " << ( field_IDs.top() + 1 ) << ";" << std::endl;
};

void RK_CALL protobuf_schemer::start_repeated_pair( const std::string& aTypeName1, const std::string& aTypeName2,
                                                    const std::string& aName ) {
  repeat_state.push( 11 );
  *file_stream << "  repeated " << aTypeName1 << " " << aName << "_key = " << field_IDs.top() << ";" << std::endl
               << "  repeated " << aTypeName2 << " " << aName << "_value = " << ( field_IDs.top() + 1 ) << ";"
               << std::endl;
};

void RK_CALL protobuf_schemer::finish_repeated_pair() {
  repeat_state.pop();
  field_IDs.top() += 2;
};


}; // serialization


}; // ReaK
