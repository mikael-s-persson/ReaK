
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

#include <ReaK/core/serialization/objtree_archiver.hpp>

#include <ReaK/core/base/named_object.hpp>
#include <ReaK/core/rtti/rtti.hpp>

#include <ReaK/core/serialization/scheme_builder.hpp>
#include <ReaK/core/serialization/type_schemes.hpp>

#include <ReaK/core/serialization/archiving_exceptions.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <map>

#include <algorithm>
#include <cctype>

namespace ReaK {

namespace serialization {


std::string::iterator xml_field_editor::mark_field( std::string::iterator it_prev, std::string::iterator it_end,
                                                    const std::string& fld_name,
                                                    const shared_ptr< type_scheme >& scheme ) {
  if( !scheme )
    return it_prev;

  std::string& xml_src = ( *( p_parent->get_object_graph() ) )[node].xml_src;

  if( scheme->is_single_field() ) { // must be a primitive scheme.

    std::string test_seq = "<" + fld_name;
    std::string::iterator it = std::search( it_prev, it_end, test_seq.begin(), test_seq.end() );
    if( it == it_end )
      return it_prev;
    src_markers.push_back( it - xml_src.begin() );
    field_schemes.push_back( scheme );
    field_names.push_back( fld_name );
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 = std::search( it, it_end, test_seq.begin(), test_seq.end() );
    if( it_end2 != it_end )
      it_end2 += test_seq.length();
    it_prev = it_end2;

  } else if( scheme->getObjectType() == serializable_obj_scheme::getStaticObjectType() ) {

    std::string test_seq = "<" + fld_name;
    std::string::iterator it = std::search( it_prev, it_end, test_seq.begin(), test_seq.end() );
    if( it == it_end )
      return it_prev;
    it = std::find( it, it_end, '>' );
    while( ( it != it_end ) && std::isspace( *it ) )
      ++it;
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 = std::search( it, it_end, test_seq.begin(), test_seq.end() );

    for( std::size_t i = 0; i < scheme->get_field_count(); ++i ) {
      std::pair< std::string, shared_ptr< type_scheme > > fld = scheme->get_field( i );
      it = mark_field( it, it_end2, fld.first, fld.second );
    };

    if( it_end2 != it_end )
      it_end2 += test_seq.length();
    it_prev = it_end2;

  } else if( scheme->getObjectType() == serializable_ptr_scheme::getStaticObjectType() ) {

    std::string test_seq = "<" + fld_name;
    std::string::iterator it = std::search( it_prev, it_end, test_seq.begin(), test_seq.end() );
    if( it == it_end )
      return it_prev;
    src_markers.push_back( it - xml_src.begin() );
    field_schemes.push_back( scheme );
    field_names.push_back( fld_name );
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 = std::search( it, it_end, test_seq.begin(), test_seq.end() );
    if( it_end2 != it_end )
      it_end2 += test_seq.length();
    it_prev = it_end2;

  } else if( scheme->getObjectType() == vector_type_scheme::getStaticObjectType() ) {

    std::string test_seq = "<" + fld_name + "_count";
    std::string::iterator it = std::search( it_prev, it_end, test_seq.begin(), test_seq.end() );
    if( it == it_end )
      return it_prev;
    src_markers.push_back( it - xml_src.begin() );
    field_schemes.push_back( shared_ptr< type_scheme >( new primitive_scheme< unsigned int >() ) );
    field_names.push_back( fld_name + "_count" );
    test_seq = "</" + fld_name + "_count>";
    std::string::iterator it_end2 = std::search( it, it_end, test_seq.begin(), test_seq.end() );
    if( it_end2 != it_end )
      it_end2 += test_seq.length();
    it_prev = it_end2;

    std::size_t i = 0;
    std::stringstream ss;
    ss << fld_name << "_q[" << i << "]";
    it = mark_field( it_prev, it_end, ss.str(), scheme->get_field( 0 ).second );
    while( it != it_prev ) {
      it_prev = it;
      ++i;
      ss.str( "" );
      ss << fld_name << "_q[" << i << "]";
      it = mark_field( it_prev, it_end, ss.str(), scheme->get_field( 0 ).second );
    };

  } else if( scheme->getObjectType() == map_type_scheme::getStaticObjectType() ) {

    std::string test_seq = "<" + fld_name + "_count";
    std::string::iterator it = std::search( it_prev, it_end, test_seq.begin(), test_seq.end() );
    if( it == it_end )
      return it_prev;
    src_markers.push_back( it - xml_src.begin() );
    field_schemes.push_back( shared_ptr< type_scheme >( new primitive_scheme< unsigned int >() ) );
    field_names.push_back( fld_name + "_count" );
    test_seq = "</" + fld_name + "_count>";
    std::string::iterator it_end2 = std::search( it, it_end, test_seq.begin(), test_seq.end() );
    if( it_end2 != it_end )
      it_end2 += test_seq.length();
    it_prev = it_end2;

    std::size_t i = 0;
    std::stringstream ss;
    ss << fld_name << "_key[" << i << "]";
    it = mark_field( it_prev, it_end, ss.str(), scheme->get_field( 0 ).second );
    while( it != it_prev ) {
      it_prev = it;
      ss.str( "" );
      ss << fld_name << "_value[" << i << "]";
      it = mark_field( it_prev, it_end, ss.str(), scheme->get_field( 1 ).second );
      it_prev = it;
      ++i;
      ss.str( "" );
      ss << fld_name << "_key[" << i << "]";
      it = mark_field( it_prev, it_end, ss.str(), scheme->get_field( 0 ).second );
    };
  };
  return it_prev;
};

std::size_t xml_field_editor::get_field_index( const std::string& aName ) const {
  for( std::size_t i = 0; i < src_markers.size(); ++i ) {
    if( aName == field_names[i] )
      return i;
  };
  return src_markers.size();
};

shared_ptr< type_scheme > xml_field_editor::get_type_scheme() const {
  shared_ptr< serializable > cur_ptr = ( *( p_parent->get_object_graph() ) )[node].p_obj;
  if( !cur_ptr )
    return shared_ptr< type_scheme >();
  std::string t_name = cur_ptr->getObjectType()->TypeName();
  if( t_name == "" )
    return shared_ptr< type_scheme >();
  std::map< std::string, shared_ptr< type_scheme > >::iterator itm = get_global_schemes().find( t_name );
  if( ( itm == get_global_schemes().end() ) || ( !( itm->second ) ) )
    return shared_ptr< type_scheme >();
  return itm->second;
};

const std::string& xml_field_editor::get_complete_src() const {
  return ( *( p_parent->get_object_graph() ) )[node].xml_src;
};

void xml_field_editor::set_complete_src( const std::string& aXMLSrc ) {
  ( *( p_parent->get_object_graph() ) )[node].xml_src = aXMLSrc;
  p_parent->ot_input_arc.load_current_from_node( node );
  shared_ptr< serializable > cur_ptr = ( *( p_parent->get_object_graph() ) )[node].p_obj;
  src_markers.clear();
  field_schemes.clear();
  field_names.clear();
  if( !cur_ptr )
    return;
  cur_ptr->load( p_parent->ot_input_arc, cur_ptr->getObjectType()->TypeVersion() );

  std::string t_name = cur_ptr->getObjectType()->TypeName();
  if( t_name == "" )
    return;
  std::map< std::string, shared_ptr< type_scheme > >::iterator itm = get_global_schemes().find( t_name );
  if( ( itm == get_global_schemes().end() ) || ( !( itm->second ) ) )
    return;
  std::string::iterator it_prev = ( *( p_parent->get_object_graph() ) )[node].xml_src.begin();
  for( std::size_t i = 0; i < itm->second->get_field_count(); ++i ) {
    std::pair< std::string, shared_ptr< type_scheme > > fld = itm->second->get_field( i );
    it_prev = mark_field( it_prev, ( *( p_parent->get_object_graph() ) )[node].xml_src.end(), fld.first, fld.second );
  };
};

std::string xml_field_editor::get_object_name( object_node_desc aNode ) const {
  return p_parent->get_object_name( aNode );
};

std::string xml_field_editor::get_object_name() const { return p_parent->get_object_name( node ); };

xml_field_editor::xml_field_editor( objtree_editor* aParent, object_node_desc aNode )
    : p_parent( aParent ), node( aNode ), src_markers(), field_schemes(), field_names() {

  shared_ptr< serializable > cur_ptr = ( *( p_parent->get_object_graph() ) )[node].p_obj;
  if( !cur_ptr )
    return;
  std::string t_name = cur_ptr->getObjectType()->TypeName();
  if( t_name == "" )
    return;
  std::map< std::string, shared_ptr< type_scheme > >::iterator itm = get_global_schemes().find( t_name );
  if( ( itm == get_global_schemes().end() ) || ( !( itm->second ) ) )
    return;
  std::string::iterator it_prev = ( *( p_parent->get_object_graph() ) )[node].xml_src.begin();
  for( std::size_t i = 0; i < itm->second->get_field_count(); ++i ) {
    std::pair< std::string, shared_ptr< type_scheme > > fld = itm->second->get_field( i );
    it_prev = mark_field( it_prev, ( *( p_parent->get_object_graph() ) )[node].xml_src.end(), fld.first, fld.second );
  };
};

std::size_t xml_field_editor::get_total_field_count() const { return src_markers.size(); };

std::pair< std::string, shared_ptr< type_scheme > > xml_field_editor::get_field( std::size_t aIndex ) const {
  return std::pair< std::string, shared_ptr< type_scheme > >( field_names[aIndex], field_schemes[aIndex] );
};

std::string xml_field_editor::get_field_src( std::size_t aIndex ) const {
  if( aIndex >= src_markers.size() )
    return "";
  else if( aIndex + 1 == src_markers.size() )
    return std::string( ( *( p_parent->get_object_graph() ) )[node].xml_src.begin() + src_markers[aIndex],
                        ( *( p_parent->get_object_graph() ) )[node].xml_src.end() );
  else
    return std::string( ( *( p_parent->get_object_graph() ) )[node].xml_src.begin() + src_markers[aIndex],
                        ( *( p_parent->get_object_graph() ) )[node].xml_src.begin() + src_markers[aIndex + 1] );
};

std::string xml_field_editor::get_field_src( const std::string& aName ) const {
  return get_field_src( get_field_index( aName ) );
};

std::string xml_field_editor::get_field_value( std::size_t aIndex ) const {
  std::string& xml_src = ( *( p_parent->get_object_graph() ) )[node].xml_src;
  std::string::iterator it = xml_src.begin() + src_markers[aIndex];
  std::string test_str;
  if( field_schemes[aIndex]->getObjectType() == serializable_ptr_scheme::getStaticObjectType() )
    test_str = "object_ID=\"";
  else
    test_str = ">\"";

  it = std::search( it, xml_src.end(), test_str.begin(), test_str.end() );
  if( it == xml_src.end() )
    return "";
  it += test_str.length();
  std::string::iterator it_end = std::find( it, xml_src.end(), '\"' );
  if( field_schemes[aIndex]->getObjectType() == serializable_ptr_scheme::getStaticObjectType() ) {
    object_node_desc fld_node = 0;
    std::stringstream( std::string( it, it_end ) ) >> fld_node;
    return get_object_name( fld_node );
  } else
    return std::string( it, it_end );
};

std::string xml_field_editor::get_field_value( const std::string& aName ) const {
  return get_field_value( get_field_index( aName ) );
};

void xml_field_editor::set_field_value( std::size_t aIndex, const std::string& aValue ) {
  std::string& xml_src = ( *( p_parent->get_object_graph() ) )[node].xml_src;
  std::string::iterator it = xml_src.begin() + src_markers[aIndex];
  std::string test_str;
  if( field_schemes[aIndex]->getObjectType() == serializable_ptr_scheme::getStaticObjectType() )
    test_str = "object_ID=\"";
  else
    test_str = ">\"";
  it = std::search( it, xml_src.end(), test_str.begin(), test_str.end() );
  if( it == xml_src.end() )
    return;
  it += test_str.length();
  std::string::iterator it_end = std::find( it, xml_src.end(), '\"' );
  std::size_t orig_len = it_end - it;
  std::string new_xml_src( xml_src.begin(), it );
  std::ptrdiff_t len_diff = aValue.length() - orig_len;

  if( field_schemes[aIndex]->getObjectType() == serializable_ptr_scheme::getStaticObjectType() ) {
    object_node_desc orig_node, new_node;
    std::stringstream( std::string( it, it_end ) ) >> orig_node;
    {
      std::string test_str3 = "(ID:";
      std::string::const_iterator it3 = std::search( aValue.begin(), aValue.end(), test_str3.begin(), test_str3.end() );
      if( it3 == aValue.end() )
        new_node = 0;
      else {
        it3 += test_str3.length();
        std::string::const_iterator it3_end = std::find( it3, aValue.end(), ')' );
        std::stringstream( std::string( it3, it3_end ) ) >> new_node;
      };
    };
    // first, check if the original node appeared anywhere else in the same xml-source.
    std::size_t orig_count = 0;
    {
      std::stringstream ss2;
      ss2 << "object_ID=\"" << orig_node << "\"";
      std::string test_str2 = ss2.str();
      std::string::iterator it2 = std::search( xml_src.begin(), xml_src.end(), test_str2.begin(), test_str2.end() );
      while( it2 != xml_src.end() ) {
        ++orig_count;
        ++it2;
        it2 = std::search( it2, xml_src.end(), test_str2.begin(), test_str2.end() );
      };
    };
    // then, check if the new node appears anywhere in the xml-source already.
    std::size_t new_count = 0;
    {
      std::stringstream ss2;
      ss2 << "object_ID=\"" << new_node << "\"";
      std::string test_str2 = ss2.str();
      std::string::iterator it2 = std::search( xml_src.begin(), xml_src.end(), test_str2.begin(), test_str2.end() );
      while( it2 != xml_src.end() ) {
        ++orig_count;
        ++it2;
        it2 = std::search( it2, xml_src.end(), test_str2.begin(), test_str2.end() );
      };
    };
    if( new_count == 0 ) {
      if( orig_count == 1 ) // the nodes must be swapped.
        p_parent->replace_child( node, new_node, orig_node );
      else // the new-node must be linked to the parent without removing the existing link.
        p_parent->create_child( node, new_node );
    } else if( orig_count == 1 )
      p_parent->sever_child( node, orig_node );

    std::stringstream ss4;
    ss4 << new_node;
    std::string new_node_str = ss4.str();
    len_diff = new_node_str.length() - orig_len;
    new_xml_src.append( new_node_str );
  } else {
    new_xml_src.append( aValue );
  };
  new_xml_src.append( it_end, xml_src.end() );
  for( std::size_t j = aIndex; j < src_markers.size(); ++j )
    src_markers[j] += len_diff;
  xml_src = std::move( new_xml_src );
  p_parent->ot_input_arc.load_current_from_node( node );
  shared_ptr< serializable > cur_ptr = ( *( p_parent->get_object_graph() ) )[node].p_obj;
  if( cur_ptr )
    cur_ptr->load( p_parent->ot_input_arc, cur_ptr->getObjectType()->TypeVersion() );
};

void xml_field_editor::set_field_value( const std::string& aName, const std::string& aValue ) {
  set_field_value( get_field_index( aName ), aValue );
};


void xml_field_editor::set_field_newptr( std::size_t aIndex, const shared_ptr< serializable >& aNewPtr ) {
  if( field_schemes[aIndex]->getObjectType() != serializable_ptr_scheme::getStaticObjectType() )
    return;
  std::string& xml_src = ( *( p_parent->get_object_graph() ) )[node].xml_src;

  std::string::iterator it = xml_src.begin() + src_markers[aIndex];
  std::string test_str = "object_ID=\"";
  it = std::search( it, xml_src.end(), test_str.begin(), test_str.end() );
  if( it == xml_src.end() )
    return;
  it += test_str.length();
  std::string::iterator it_end = std::find( it, xml_src.end(), '\"' );

  object_node_desc orig_node = 0;
  std::stringstream( std::string( it, it_end ) ) >> orig_node;
  // first, check if the original node appeared anywhere else in the same xml-source.
  std::size_t orig_count = 0;
  {
    std::stringstream ss2;
    ss2 << "object_ID=\"" << orig_node << "\"";
    std::string test_str2 = ss2.str();
    std::string::iterator it2 = std::search( xml_src.begin(), xml_src.end(), test_str2.begin(), test_str2.end() );
    while( it2 != xml_src.end() ) {
      ++orig_count;
      ++it2;
      it2 = std::search( it2, xml_src.end(), test_str2.begin(), test_str2.end() );
    };
  };

  object_node_desc new_node = 0;
  if( orig_count == 1 ) // the original node must be replaced.
    new_node = p_parent->add_new_object( aNewPtr, node, orig_node );
  else // the new-node must be linked to the parent without removing the existing link.
    new_node = p_parent->add_new_object( aNewPtr, node );

  std::size_t orig_len = it_end - it;
  std::stringstream ss3;
  ss3 << new_node;
  std::string aValue = ss3.str();
  std::ptrdiff_t len_diff = aValue.length() - orig_len;
  std::string new_xml_src( xml_src.begin(), it );
  new_xml_src.append( aValue );
  new_xml_src.append( it_end, xml_src.end() );
  for( std::size_t j = aIndex; j < src_markers.size(); ++j )
    src_markers[j] += len_diff;
  xml_src = std::move( new_xml_src );
  p_parent->ot_input_arc.load_current_from_node( node );
  shared_ptr< serializable > cur_ptr = ( *( p_parent->get_object_graph() ) )[node].p_obj;
  if( cur_ptr )
    cur_ptr->load( p_parent->ot_input_arc, cur_ptr->getObjectType()->TypeVersion() );
};

void xml_field_editor::set_field_newptr( const std::string& aName, const shared_ptr< serializable >& aNewPtr ) {
  set_field_newptr( get_field_index( aName ), aNewPtr );
};


char objtree_iarchive::getNextChar() {
  char c = '\0';
  current_ss->get( c );
  while( ( c == ' ' ) || ( c == '\t' ) || ( c == '\n' ) || ( c == '\r' ) )
    current_ss->get( c );
  return c;
};

std::string objtree_iarchive::readToken() {
  std::string result;
  char c = getNextChar();
  if( c != '<' )
    return result;
  c = getNextChar();
  if( ( c == '!' ) || ( c == '?' ) ) {
    char line_str[512];
    current_ss->getline( line_str, 512 );
    return readToken();
  };
  while( c != '>' ) {
    result += c;
    current_ss->get( c );
  };
  return result;
};


void objtree_iarchive::skipToEndToken( const std::string& name ) {
  std::string token = readToken();
  trimStr( token );
  while( token != "/" + name ) {
    token = readToken();
    trimStr( token );
  };
};

void objtree_iarchive::trimStr( std::string& s ) {
  unsigned int i = 0;
  for( ; ( ( i < s.size() ) && ( ( s[i] == ' ' ) || ( s[i] == '\t' ) || ( s[i] == '\n' ) || ( s[i] == '\r' ) ) ); ++i )
    ;
  std::string result;
  for( ; ( ( i < s.size() ) && ( s[i] != ' ' ) && ( s[i] != '\t' ) && ( s[i] != '\n' ) && ( s[i] != '\r' ) ); ++i )
    result += s[i];
  s = result;
};

bool objtree_iarchive::readNamedValue( const std::string& value_name, std::string& value_str ) {
  std::string token = readToken();
  trimStr( token );
  if( ( value_name.empty() ) || ( token != value_name ) )
    return false;

  char c = '\0';
  current_ss->get( c );
  while( c != '\"' )
    current_ss->get( c );

  value_str.clear();
  current_ss->get( c );
  while( c != '\"' ) {
    value_str += c;
    current_ss->get( c );
  };

  token = readToken();
  unsigned int i = 0;
  for( ; ( ( i < token.size() ) && ( token[i] != '/' ) ); ++i )
    ;
  std::string tmp;
  ++i;
  for( ; ( ( i < token.size() ) && ( token[i] != value_name[0] ) ); ++i )
    ;
  for( ; ( ( i < token.size() ) && ( tmp.size() < value_name.size() ) ); ++i )
    tmp += token[i];
  if( tmp != value_name )
    return false;
  return true;
};

archive_object_header objtree_iarchive::readHeader( const std::string& obj_name,
                                                    std::vector< unsigned int >& outTypeID ) {
  archive_object_header result;
  outTypeID.clear();

  std::string token = readToken();
  if( token.empty() )
    return result;

  std::string name;
  unsigned int i = 0;
  for( ; ( ( i < token.size() ) && ( token[i] == ' ' ) ); ++i )
    ;
  for( ; ( ( i < token.size() ) && ( token[i] != ' ' ) ); ++i )
    name += token[i];

  if( ( name != obj_name ) || ( i == token.size() ) )
    return result;

  std::map< std::string, std::string > values;
  while( i < token.size() ) {
    for( ; ( ( i < token.size() )
             && ( ( token[i] == ' ' ) || ( token[i] == '\t' ) || ( token[i] == '\n' ) || ( token[i] == '\r' ) ) );
         ++i )
      ;
    std::string value_key;
    for( ; ( ( i < token.size() ) && ( token[i] != ' ' ) && ( token[i] != '=' ) ); ++i )
      value_key += token[i];
    std::string value_str;
    for( ; ( ( i < token.size() ) && ( token[i] != '\"' ) ); ++i )
      ;
    ++i;
    for( ; ( ( i < token.size() ) && ( token[i] != '\"' ) ); ++i )
      value_str += token[i];
    ++i;
    values[value_key] = value_str;
  };

  std::string IDstr = values["type_ID"];
  if( !IDstr.empty() ) {
    for( i = 0; i < IDstr.size(); ++i ) {
      std::string numstr;
      for( ; ( ( i < IDstr.size() ) && ( IDstr[i] != '.' ) ); ++i )
        numstr += IDstr[i];
      outTypeID.push_back( strtoul( numstr.c_str(), nullptr, 0 ) );
    };
  };

  if( values["version"].empty() )
    result.type_version = 0;
  else
    result.type_version = strtoul( values["version"].c_str(), nullptr, 0 );

  if( values["object_ID"].empty() )
    result.object_ID = 0;
  else
    result.object_ID = strtoul( values["object_ID"].c_str(), nullptr, 0 );

  result.is_external = false;

  return result;
};


void objtree_iarchive::load_current_from_node( object_node_desc aNode ) {
  current_ss = shared_ptr< std::stringstream >( new std::stringstream( ( *obj_graph )[aNode].xml_src ) );
};

objtree_iarchive::objtree_iarchive( const shared_ptr< object_graph >& aObjGraph, object_node_desc aRoot )
    : obj_graph( aObjGraph ), obj_graph_root( aRoot ) {

  if( num_vertices( *obj_graph ) == 0 )
    obj_graph_root
      = add_vertex( *obj_graph ); // add a root node. This case doesn't make much sense, it means the graph is empty.

  current_ss = shared_ptr< std::stringstream >( new std::stringstream( ( *obj_graph )[obj_graph_root].xml_src ) );
};

objtree_iarchive::~objtree_iarchive(){};


iarchive& RK_CALL objtree_iarchive::load_serializable_ptr( serializable_shared_pointer& Item ) {
  return objtree_iarchive::load_serializable_ptr(
    std::pair< std::string, serializable_shared_pointer& >( "Item", Item ) );
};

iarchive& RK_CALL
  objtree_iarchive::load_serializable_ptr( const std::pair< std::string, serializable_shared_pointer& >& Item ) {
  Item.second = serializable_shared_pointer();
  typedef boost::graph_traits< object_graph >::vertex_descriptor Vertex;

  std::vector< unsigned int > typeID;
  archive_object_header hdr = readHeader( Item.first, typeID );
  if( ( typeID.empty() ) || ( hdr.type_version == 0 ) || ( hdr.object_ID == 0 ) ) {
    skipToEndToken( Item.first );
    return *this;
  };

  if( ( hdr.object_ID < num_vertices( *obj_graph ) )
      && ( ( *obj_graph )[static_cast< Vertex >( hdr.object_ID )].p_obj ) ) {
    Item.second = ( *obj_graph )[static_cast< Vertex >( hdr.object_ID )].p_obj;
    skipToEndToken( Item.first );

    // re-read the xml source
    shared_ptr< std::stringstream > tmp_ss = current_ss;
    current_ss = shared_ptr< std::stringstream >(
      new std::stringstream( ( *obj_graph )[static_cast< Vertex >( hdr.object_ID )].xml_src ) );

    Item.second->load( *this, hdr.type_version );

    current_ss = tmp_ss;
  };

  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_serializable( serializable& Item ) {
  return objtree_iarchive::load_serializable( std::pair< std::string, serializable& >( "Item", Item ) );
};

iarchive& RK_CALL objtree_iarchive::load_serializable( const std::pair< std::string, serializable& >& Item ) {
  archive_object_header hdr;

  std::vector< unsigned int > typeID;
  hdr = readHeader( Item.first, typeID );
  if( ( hdr.type_ID == nullptr ) || ( hdr.type_version == 0 ) ) {
    skipToEndToken( Item.first );
    return *this;
  };

  Item.second.load( *this, hdr.type_version );

  skipToEndToken( Item.first );
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_char( char& i ) {
  return objtree_iarchive::load_char( std::pair< std::string, char& >( "char", i ) );
};

iarchive& RK_CALL objtree_iarchive::load_char( const std::pair< std::string, char& >& i ) {
  std::string value_str;
  if( readNamedValue( i.first, value_str ) ) {
    if( value_str.empty() )
      i.second = 0;
    else {
      int temp;
      std::stringstream( value_str ) >> temp;
      i.second = char( temp );
    };
  } else
    i.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_char( unsigned char& u ) {
  return objtree_iarchive::load_unsigned_char( std::pair< std::string, unsigned char& >( "unsigned_char", u ) );
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_char( const std::pair< std::string, unsigned char& >& u ) {
  std::string value_str;
  if( readNamedValue( u.first, value_str ) ) {
    if( value_str.empty() )
      u.second = 0;
    else {
      unsigned int temp;
      std::stringstream( value_str ) >> temp;
      u.second = char( temp );
    };
  } else
    u.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_int( std::ptrdiff_t& i ) {
  return objtree_iarchive::load_int( std::pair< std::string, std::ptrdiff_t& >( "int", i ) );
};

iarchive& RK_CALL objtree_iarchive::load_int( const std::pair< std::string, std::ptrdiff_t& >& i ) {
  std::string value_str;
  if( readNamedValue( i.first, value_str ) ) {
    if( value_str.empty() )
      i.second = 0;
    else
      std::stringstream( value_str ) >> i.second;
  } else
    i.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_int( std::size_t& u ) {
  return objtree_iarchive::load_unsigned_int( std::pair< std::string, std::size_t& >( "unsigned_int", u ) );
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_int( const std::pair< std::string, std::size_t& >& u ) {
  std::string value_str;
  if( readNamedValue( u.first, value_str ) ) {
    if( value_str.empty() )
      u.second = 0;
    else
      std::stringstream( value_str ) >> u.second;
  } else
    u.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_float( float& f ) {
  return objtree_iarchive::load_float( std::pair< std::string, float& >( "real", f ) );
};

iarchive& RK_CALL objtree_iarchive::load_float( const std::pair< std::string, float& >& f ) {
  std::string value_str;
  if( readNamedValue( f.first, value_str ) ) {
    if( value_str.empty() )
      f.second = 0;
    else
      std::stringstream( value_str ) >> f.second;
  } else
    f.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_double( double& d ) {
  return objtree_iarchive::load_double( std::pair< std::string, double& >( "real", d ) );
};

iarchive& RK_CALL objtree_iarchive::load_double( const std::pair< std::string, double& >& d ) {
  std::string value_str;
  if( readNamedValue( d.first, value_str ) ) {
    if( value_str.empty() )
      d.second = 0;
    else
      std::stringstream( value_str ) >> d.second;
  } else
    d.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_bool( bool& b ) {
  return objtree_iarchive::load_bool( std::pair< std::string, bool& >( "bool", b ) );
};

iarchive& RK_CALL objtree_iarchive::load_bool( const std::pair< std::string, bool& >& b ) {
  std::string value_str;
  if( readNamedValue( b.first, value_str ) ) {
    if( value_str.empty() )
      b.second = false;
    else if( value_str == "true" )
      b.second = true;
    else
      b.second = false;
  } else
    b.second = false;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_string( std::string& s ) {
  return objtree_iarchive::load_string( std::pair< std::string, std::string& >( "string", s ) );
};

iarchive& RK_CALL objtree_iarchive::load_string( const std::pair< std::string, std::string& >& s ) {
  readNamedValue( s.first, s.second );
  return *this;
};


void objtree_oarchive::register_new_object( object_node_desc aNode ) {
  mObjRegMap[( *obj_graph )[aNode].p_obj] = static_cast< unsigned int >( aNode );
};

void objtree_oarchive::unregister_object( object_node_desc aNode ) { mObjRegMap.erase( ( *obj_graph )[aNode].p_obj ); };

void objtree_oarchive::save_current_stream() { ( *obj_graph )[current_node].xml_src = current_ss->str(); };

void objtree_oarchive::load_current_from_node( object_node_desc aNode ) {
  current_ss = shared_ptr< std::stringstream >( new std::stringstream( ( *obj_graph )[aNode].xml_src ) );
  current_node = aNode;
};

void objtree_oarchive::fresh_current_node( object_node_desc aNode ) {
  current_ss = shared_ptr< std::stringstream >( new std::stringstream() );
  current_node = aNode;
};

objtree_oarchive::objtree_oarchive( const shared_ptr< object_graph >& aObjGraph, object_node_desc aRoot )
    : obj_graph( aObjGraph ), obj_graph_root( aRoot ) {
  current_ss = shared_ptr< std::stringstream >( new std::stringstream() );

  boost::graph_traits< object_graph >::vertex_iterator vi, vi_end;
  for( boost::tie( vi, vi_end ) = vertices( *obj_graph ); vi != vi_end; ++vi )
    if( ( *obj_graph )[*vi].p_obj )
      mObjRegMap[( *obj_graph )[*vi].p_obj] = static_cast< unsigned int >( *vi );
  current_node = obj_graph_root;
};

objtree_oarchive::~objtree_oarchive() { save_current_stream(); };


oarchive& RK_CALL
  objtree_oarchive::saveToNewArchive_impl( const serializable_shared_pointer& Item, const std::string& ) {
  return objtree_oarchive::save_serializable_ptr(
    std::pair< std::string, const serializable_shared_pointer& >( "Item", Item ) );
};

oarchive& RK_CALL objtree_oarchive::saveToNewArchiveNamed_impl(
  const std::pair< std::string, const serializable_shared_pointer& >& Item, const std::string& ) {
  return objtree_oarchive::save_serializable_ptr( Item );
};


oarchive& RK_CALL objtree_oarchive::save_serializable_ptr( const serializable_shared_pointer& Item ) {
  return objtree_oarchive::save_serializable_ptr(
    std::pair< std::string, const serializable_shared_pointer& >( "Item", Item ) );
};


oarchive& RK_CALL
  objtree_oarchive::save_serializable_ptr( const std::pair< std::string, const serializable_shared_pointer& >& Item ) {


  if( Item.second ) {
    std::map< serializable_shared_pointer, std::size_t >::const_iterator it = mObjRegMap.find( Item.second );
    std::size_t object_ID = 0;
    if( it != mObjRegMap.end() ) {
      object_ID = it->second;
    } else {
      object_ID = static_cast< std::size_t >( add_vertex( *obj_graph ) );
      ( *obj_graph )[object_ID].p_obj = Item.second;
      mObjRegMap[Item.second] = object_ID;
    };

    if( !edge( current_node, object_ID, *obj_graph ).second )
      add_edge( current_node, object_ID, *obj_graph );

    rtti::so_type* obj_type = Item.second->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    unsigned int type_version = obj_type->TypeVersion();

    ( *current_ss ) << "<" << Item.first << " type_ID=\"";
    while( *type_ID ) {
      ( *current_ss ) << *type_ID << ".";
      ++type_ID;
    };
    ( *current_ss ) << "0\" version=\"" << type_version << "\" object_ID=\"" << object_ID << "\"></" << Item.first
                    << ">" << std::endl;

    shared_ptr< std::stringstream > tmp = current_ss;
    current_ss = shared_ptr< std::stringstream >( new std::stringstream() );

    boost::graph_traits< object_graph >::vertex_descriptor tmp_v = current_node;
    current_node = object_ID;

    Item.second->save( *this, type_version );

    // grab the resulting string.
    ( *obj_graph )[object_ID].xml_src = current_ss->str();

    current_ss = tmp;
    current_node = tmp_v;
  } else {
    ( *current_ss ) << "<" << Item.first << " type_ID=\"0\" version=\"0\" object_ID=\"0\"></" << Item.first << ">"
                    << std::endl;
  };

  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_serializable( const serializable& Item ) {
  return *this & std::pair< std::string, const serializable& >( "Item", Item );
};

oarchive& RK_CALL objtree_oarchive::save_serializable( const std::pair< std::string, const serializable& >& Item ) {
  archive_object_header hdr;
  const unsigned int* type_ID = Item.second.getObjectType()->TypeID_begin();
  hdr.type_version = Item.second.getObjectType()->TypeVersion();
  hdr.object_ID = 0;
  hdr.size = 0;
  hdr.is_external = false;

  ( *current_ss ) << "<" << Item.first << " type_ID=\"";
  while( *type_ID ) {
    ( *current_ss ) << *type_ID << ".";
    ++type_ID;
  };
  ( *current_ss ) << "0\" version=\"" << hdr.type_version << "\">" << std::endl;

  Item.second.save( *this, hdr.type_version );

  ( *current_ss ) << "</" << Item.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_char( char i ) {
  return objtree_oarchive::save_char( std::pair< std::string, char >( "char", i ) );
};

oarchive& RK_CALL objtree_oarchive::save_char( const std::pair< std::string, char >& i ) {
  ( *current_ss ) << "<" << i.first << ">\"" << static_cast< int >( i.second ) << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_char( unsigned char u ) {
  return objtree_oarchive::save_unsigned_char( std::pair< std::string, unsigned char >( "unsigned_char", u ) );
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_char( const std::pair< std::string, unsigned char >& u ) {
  ( *current_ss ) << "<" << u.first << ">\"" << static_cast< unsigned int >( u.second ) << "\"</" << u.first << ">"
                  << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_int( std::ptrdiff_t i ) {
  return objtree_oarchive::save_int( std::pair< std::string, std::ptrdiff_t >( "int", i ) );
};


oarchive& RK_CALL objtree_oarchive::save_int( const std::pair< std::string, std::ptrdiff_t >& i ) {
  ( *current_ss ) << "<" << i.first << ">\"" << i.second << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_int( std::size_t u ) {
  return objtree_oarchive::save_unsigned_int( std::pair< std::string, std::size_t >( "unsigned_int", u ) );
};


oarchive& RK_CALL objtree_oarchive::save_unsigned_int( const std::pair< std::string, std::size_t >& u ) {
  ( *current_ss ) << "<" << u.first << ">\"" << u.second << "\"</" << u.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_float( float f ) {
  return objtree_oarchive::save_float( std::pair< std::string, float >( "real", f ) );
};


oarchive& RK_CALL objtree_oarchive::save_float( const std::pair< std::string, float >& f ) {
  ( *current_ss ) << "<" << f.first << ">\"" << f.second << "\"</" << f.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_double( double d ) {
  return objtree_oarchive::save_double( std::pair< std::string, double >( "real", d ) );
};


oarchive& RK_CALL objtree_oarchive::save_double( const std::pair< std::string, double >& d ) {
  ( *current_ss ) << "<" << d.first << ">\"" << d.second << "\"</" << d.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_bool( bool b ) {
  return objtree_oarchive::save_bool( std::pair< std::string, bool >( "bool", b ) );
};


oarchive& RK_CALL objtree_oarchive::save_bool( const std::pair< std::string, bool >& b ) {
  ( *current_ss ) << "<" << b.first << ">\"" << ( b.second ? "true" : "false" ) << "\"</" << b.first << ">"
                  << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_string( const std::string& s ) {
  return objtree_oarchive::save_string( std::pair< std::string, const std::string& >( "string", s ) );
};


oarchive& RK_CALL objtree_oarchive::save_string( const std::pair< std::string, const std::string& >& s ) {
  ( *current_ss ) << "<" << s.first << ">\"" << s.second << "\"</" << s.first << ">" << std::endl;
  return *this;
};


objtree_editor::objtree_editor()
    : obj_graph( new object_graph() ), obj_graph_root( add_vertex( *obj_graph ) ),
      ot_output_arc( obj_graph, obj_graph_root ), ot_input_arc( obj_graph, obj_graph_root ), obj_graph_graveyard(){};

objtree_editor::objtree_editor( const shared_ptr< object_graph >& aObjGraph, object_node_desc aRoot )
    : obj_graph( aObjGraph ), obj_graph_root( aRoot ), ot_output_arc( obj_graph, obj_graph_root ),
      ot_input_arc( obj_graph, obj_graph_root ), obj_graph_graveyard(){};


object_node_desc objtree_editor::add_new_object( const shared_ptr< serializable >& aNewObj, object_node_desc aParent,
                                                 object_node_desc aOldChild ) {
  if( !aNewObj )
    return obj_graph_root;
  object_node_desc result;
  if( obj_graph_graveyard.empty() )
    result = add_vertex( object_graph_node( aNewObj, "" ), *obj_graph );
  else {
    result = obj_graph_graveyard.top();
    obj_graph_graveyard.pop();
    ( *obj_graph )[result].p_obj = aNewObj;
    ( *obj_graph )[result].xml_src = "";
  };
  if( aOldChild ) {
    remove_edge( aParent, aOldChild, *obj_graph );
    remove_object( aOldChild ); // if this was the last in-edge on aOldChild.
  };
  add_edge( aParent, result, *obj_graph );
  ot_output_arc.register_new_object( result );
  ot_output_arc.fresh_current_node( result );
  aNewObj->save( ot_output_arc, aNewObj->getObjectType()->TypeVersion() );
  ot_output_arc.save_current_stream();
  return result;
};

void objtree_editor::remove_object( object_node_desc aNode ) {
  if( in_degree( aNode, *obj_graph ) )
    return;
  std::vector< object_node_desc > v;
  v.reserve( out_degree( aNode, *obj_graph ) );
  boost::graph_traits< object_graph >::out_edge_iterator eo, eo_end;
  for( boost::tie( eo, eo_end ) = out_edges( aNode, *obj_graph ); eo != eo_end; ++eo )
    v.push_back( target( *eo, *obj_graph ) );
  clear_vertex( aNode, *obj_graph );
  for( std::vector< object_node_desc >::iterator it = v.begin(); it != v.end(); ++it )
    remove_object( *it ); // this will only really have an effect if the node has no other in-edge.
  ot_output_arc.unregister_object( aNode );
  ( *obj_graph )[aNode].p_obj = shared_ptr< serializable >();
  ( *obj_graph )[aNode].xml_src = "";
  obj_graph_graveyard.push( aNode );
};

void objtree_editor::replace_child( object_node_desc aParent, object_node_desc aNewChild, object_node_desc aOldChild ) {
  if( ( aNewChild ) && ( !edge( aParent, aNewChild, *obj_graph ).second ) )
    add_edge( aParent, aNewChild, *obj_graph );
  if( aOldChild ) {
    remove_edge( aParent, aOldChild, *obj_graph );
    remove_object( aOldChild ); // this will only really have an effect if the node has no other in-edge.
  };
};

std::string objtree_editor::get_object_name( object_node_desc aNode ) const {
  return get_objtree_name( *obj_graph, aNode );
};

std::vector< std::string > objtree_editor::get_objects_derived_from( rtti::so_type* aType ) const {
  boost::graph_traits< serialization::object_graph >::vertex_iterator vi, vi_end;
  boost::tie( vi, vi_end ) = vertices( *obj_graph );
  std::vector< std::string > result;
  for( ; vi != vi_end; ++vi ) {
    shared_ptr< serializable > p_obj = ( *obj_graph )[*vi].p_obj;
    if( ( p_obj ) && ( p_obj->castTo( aType ) ) )
      result.push_back( get_object_name( *vi ) );
  };
  return result;
};


std::string get_objtree_name( const object_graph& obj_graph, object_node_desc node_id ) {
  shared_ptr< named_object > item_ptr = rtti::rk_dynamic_ptr_cast< named_object >( obj_graph[node_id].p_obj );
  std::stringstream ss;
  if( item_ptr )
    ss << item_ptr->getName() << " (ID:" << node_id << ")";
  else
    ss << "Object (ID:" << node_id << ")";
  return ss.str();
};

object_node_desc get_objtree_node_id( const object_graph& obj_graph, const std::string& obj_name ) {
  std::string node_id_str = obj_name.substr( obj_name.find( "(ID:" ) + 4 );
  node_id_str = node_id_str.substr( 0, node_id_str.find( ")" ) );
  std::stringstream ss( node_id_str );
  std::size_t result = 0;
  ss >> result;
  return result;
};


}; // serialization

}; // ReaK
