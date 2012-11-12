
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

#include "objtree_archiver.hpp"


#include "base/shared_object.hpp"

#include <string>
#include <fstream>
#include <iomanip>

#include <algorithm>


namespace ReaK {

namespace serialization {





char objtree_iarchive::getNextChar() {
  char c;
  current_ss->get(c);
  while((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r'))
    current_ss->get(c);
  return c;
};

std::string objtree_iarchive::readToken() {
  std::string result;
  char c = getNextChar();
  if(c != '<')
    return result;
  c = getNextChar();
  if((c == '!') || (c == '?')) {
    char line_str[512];
    current_ss->getline(line_str,512);
    return readToken();
  };
  while(c != '>') {
    result += c;
    current_ss->get(c);
  };
  return result;
};


void objtree_iarchive::skipToEndToken(const std::string& name) {
  std::string token = readToken();
  trimStr(token);
  while(token != "/" + name) {
    token = readToken();
    trimStr(token);
  };
};

void objtree_iarchive::trimStr(std::string& s) {
  unsigned int i=0;
  for(;((i < s.size()) && ((s[i] == ' ') || (s[i] == '\t') || (s[i] == '\n') || (s[i] == '\r')));++i) ;
  std::string result;
  for(;((i < s.size()) && (s[i] != ' ') && (s[i] != '\t') && (s[i] != '\n') && (s[i] != '\r'));++i)
    result += s[i];
  s = result;
};

bool objtree_iarchive::readNamedValue(const std::string& value_name,std::string& value_str) {
  std::string token = readToken();
  trimStr(token);
  if((value_name.empty()) || (token != value_name))
    return false;
  
  char c;
  current_ss->get(c);
  while(c != '\"')
    current_ss->get(c);
  
  value_str.clear();
  current_ss->get(c);
  while(c != '\"') {
    value_str += c;
    current_ss->get(c);
  };
  
  token = readToken();
  unsigned int i=0;
  for(;((i<token.size()) && (token[i] != '/'));++i) ;
  std::string tmp;
  ++i;
  for(;((i<token.size()) && (token[i] != value_name[0]));++i) ;
  for(;((i<token.size()) && (tmp.size() < value_name.size()));++i)
    tmp += token[i];
  if(tmp != value_name)
    return false;
  return true;
};

archive_object_header objtree_iarchive::readHeader(const std::string& obj_name) {
  archive_object_header result;
  
  std::string token = readToken();
  if(token.empty())
    return result;
  
  std::string name;
  unsigned int i=0;
  for(;((i < token.size()) && (token[i] == ' '));++i) ;
  for(;((i < token.size()) && (token[i] != ' '));++i)
    name += token[i];
  
  if((name != obj_name) || (i == token.size()))
    return result;
  
  std::map<std::string,std::string> values;
  while(i < token.size()) {
    for(;((i < token.size()) && ((token[i] == ' ') || (token[i] == '\t') || (token[i] == '\n') || (token[i] == '\r')));++i) ;
    std::string value_key;
    for(;((i < token.size()) && (token[i] != ' ') && (token[i] != '='));++i)
      value_key += token[i];
    std::string value_str;
    for(;((i < token.size()) && (token[i] != '\"'));++i) ;
    ++i;
    for(;((i < token.size()) && (token[i] != '\"'));++i)
      value_str += token[i];
    ++i;
    values[value_key] = value_str;
  };
  
  std::string IDstr = values["type_ID"];
  if(IDstr.empty())
    result.type_ID = NULL;
  else {
    std::vector<unsigned int> nums;
    for(i=0;i<IDstr.size();++i) {
      std::string numstr;
      for(;((i < IDstr.size()) && (IDstr[i] != '.'));++i) 
        numstr += IDstr[i];
      nums.push_back(strtoul(numstr.c_str(),NULL,0));
    };
    result.type_ID = new unsigned int[nums.size()];
    std::copy(nums.begin(),nums.end(),result.type_ID);
  };
  
  if(values["version"].empty())
    result.type_version = 0;
  else
    result.type_version = strtoul(values["version"].c_str(),NULL,0);
  
  if(values["object_ID"].empty())
    result.object_ID = 0;
  else
    result.object_ID = strtoul(values["object_ID"].c_str(),NULL,0);
  
  result.is_external = false;
  
  return result;
};

objtree_iarchive::objtree_iarchive(const shared_ptr< object_graph >& aObjGraph) : obj_graph(aObjGraph) {
  
  current_ss = shared_ptr< std::stringstream >(new std::stringstream());
  
  if(num_vertices(*obj_graph) == 0)
    add_vertex(*obj_graph);  // add a root node. This case doesn't make much sense, it means the graph is empty.
  
};

objtree_iarchive::~objtree_iarchive() { };


iarchive& RK_CALL objtree_iarchive::load_serializable_ptr(serializable_shared_pointer& Item) {
  return objtree_iarchive::load_serializable_ptr(std::pair<std::string, serializable_shared_pointer& >("Item",Item));
};

iarchive& RK_CALL objtree_iarchive::load_serializable_ptr(const std::pair<std::string, serializable_shared_pointer& >& Item) {
  Item.second = serializable_shared_pointer();
  typedef boost::graph_traits< object_graph >::vertex_descriptor Vertex;
  
  archive_object_header hdr = readHeader(Item.first);
  if((hdr.type_ID == NULL) || 
     (hdr.type_version == 0) || 
     (hdr.object_ID == 0)) {
    skipToEndToken(Item.first);
    delete[] hdr.type_ID;
    return *this;
  };
  
  if((hdr.object_ID < num_vertices(*obj_graph)) && 
     ((*obj_graph)[static_cast<Vertex>(hdr.object_ID)].p_obj)) {
    Item.second = (*obj_graph)[static_cast<Vertex>(hdr.object_ID)].p_obj;
    skipToEndToken(Item.first);
    delete[] hdr.type_ID;
    
    // re-read the xml source
    shared_ptr< std::stringstream > tmp_ss = current_ss;
    current_ss = shared_ptr< std::stringstream >(new std::stringstream((*obj_graph)[static_cast<Vertex>(hdr.object_ID)].xml_src));
    
    Item.second->load(*this, hdr.type_version);
    
    current_ss = tmp_ss;
  };
  
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_serializable(serializable& Item) {
  return objtree_iarchive::load_serializable(std::pair<std::string, serializable&>("Item",Item));
};

iarchive& RK_CALL objtree_iarchive::load_serializable(const std::pair<std::string, serializable& >& Item) {
  archive_object_header hdr;

  hdr = readHeader(Item.first);
  if((hdr.type_ID == NULL) || (hdr.type_version == 0)) {
    skipToEndToken(Item.first);
    return *this;
  };

  Item.second.load(*this,hdr.type_version);

  skipToEndToken(Item.first);
  delete[] hdr.type_ID;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_char(char& i) {
  return objtree_iarchive::load_char(std::pair<std::string, char& >("char",i));
};

iarchive& RK_CALL objtree_iarchive::load_char(const std::pair<std::string, char& >& i) {
  std::string value_str;
  if(readNamedValue(i.first,value_str)) {
    if(value_str.empty())
      i.second = 0;
    else {
      int temp;
      std::stringstream(value_str) >> temp;
      i.second = char(temp);
    };
  } else
    i.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_char(unsigned char& u) {
  return objtree_iarchive::load_unsigned_char(std::pair<std::string, unsigned char& >("unsigned char",u));
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_char(const std::pair<std::string, unsigned char& >& u) {
  std::string value_str;
  if(readNamedValue(u.first,value_str)) {
    if(value_str.empty())
      u.second = 0;
    else {
      unsigned int temp;
      std::stringstream(value_str) >> temp;
      u.second = char(temp);
    };
  } else
    u.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_int(int& i) {
  return objtree_iarchive::load_int(std::pair<std::string, int& >("int",i));
};

iarchive& RK_CALL objtree_iarchive::load_int(const std::pair<std::string, int& >& i) {
  std::string value_str;
  if(readNamedValue(i.first,value_str)) {
    if(value_str.empty())
      i.second = 0;
    else
      std::stringstream(value_str) >> i.second;
  } else
    i.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_int(unsigned int& u) {
  return objtree_iarchive::load_unsigned_int(std::pair<std::string, unsigned int& >("unsigned int",u));
};

iarchive& RK_CALL objtree_iarchive::load_unsigned_int(const std::pair<std::string, unsigned int& >& u) {
  std::string value_str;
  if(readNamedValue(u.first,value_str)) {
    if(value_str.empty())
      u.second = 0;
    else
      std::stringstream(value_str) >> u.second;
  } else
    u.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_float(float& f) {
  return objtree_iarchive::load_float(std::pair<std::string, float& >("real",f));
};

iarchive& RK_CALL objtree_iarchive::load_float(const std::pair<std::string, float& >& f) {
  std::string value_str;
  if(readNamedValue(f.first,value_str)) {
    if(value_str.empty())
      f.second = 0;
    else
      std::stringstream(value_str) >> f.second;
  } else
    f.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_double(double& d) {
  return objtree_iarchive::load_double(std::pair<std::string, double& >("real",d));
};

iarchive& RK_CALL objtree_iarchive::load_double(const std::pair<std::string, double& >& d) {
  std::string value_str;
  if(readNamedValue(d.first,value_str)) {
    if(value_str.empty())
      d.second = 0;
    else
      std::stringstream(value_str) >> d.second;
  } else
    d.second = 0;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_bool(bool& b) {
  return objtree_iarchive::load_bool(std::pair<std::string, bool& >("bool",b));
};

iarchive& RK_CALL objtree_iarchive::load_bool(const std::pair<std::string, bool& >& b) {
  std::string value_str;
  if(readNamedValue(b.first,value_str)) {
    if(value_str.empty())
      b.second = false;
    else if(value_str == "true")
      b.second = true;
    else
      b.second = false;
  } else
    b.second = false;
  return *this;
};

iarchive& RK_CALL objtree_iarchive::load_string(std::string& s) {
  return objtree_iarchive::load_string(std::pair<std::string, std::string& >("string",s));
};

iarchive& RK_CALL objtree_iarchive::load_string(const std::pair<std::string, std::string& >& s) {
  readNamedValue(s.first,s.second);
  return *this;
};









objtree_oarchive::objtree_oarchive(const shared_ptr< object_graph >& aObjGraph) : obj_graph(aObjGraph) {
  current_ss = shared_ptr< std::stringstream >(new std::stringstream());
  
  boost::graph_traits< object_graph >::vertex_iterator vi, vi_end;
  for(boost::tie(vi, vi_end) = vertices(*obj_graph); vi != vi_end; ++vi)
    mObjRegMap[(*obj_graph)[*vi].p_obj] = static_cast<unsigned int>(*vi);
  current_node = add_vertex(*obj_graph);
};

objtree_oarchive::~objtree_oarchive() { };


oarchive& RK_CALL objtree_oarchive::saveToNewArchive_impl(const serializable_shared_pointer& Item, const std::string&) {
  return objtree_oarchive::save_serializable_ptr(std::pair<std::string, const serializable_shared_pointer& >("Item",Item));
};

oarchive& RK_CALL objtree_oarchive::saveToNewArchiveNamed_impl(const std::pair<std::string, const serializable_shared_pointer& >& Item, const std::string&) {
  return objtree_oarchive::save_serializable_ptr(Item);
};


oarchive& RK_CALL objtree_oarchive::save_serializable_ptr(const serializable_shared_pointer& Item) {
  return objtree_oarchive::save_serializable_ptr(std::pair<std::string, const serializable_shared_pointer& >("Item",Item));
};


oarchive& RK_CALL objtree_oarchive::save_serializable_ptr(const std::pair<std::string, const serializable_shared_pointer& >& Item) {
  
  
  if(Item.second) {
    std::map< serializable_shared_pointer, unsigned int>::const_iterator it = mObjRegMap.find(Item.second);
    unsigned int object_ID = 0;
    if(it != mObjRegMap.end()) {
      object_ID = it->second;
    } else {
      object_ID = static_cast<unsigned int>(add_vertex(*obj_graph));
      (*obj_graph)[object_ID].p_obj = Item.second;
      mObjRegMap[Item.second] = object_ID;
    };
    
    add_edge(current_node, object_ID, *obj_graph);
    
    rtti::so_type::shared_pointer obj_type = Item.second->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    unsigned int type_version = obj_type->TypeVersion();
    
    (*current_ss) << "<" << Item.first << " type_ID=\"";
    while(*type_ID) {
      (*current_ss) << *type_ID << "."; 
      ++type_ID;
    };
    (*current_ss) << "0\" version=\"" << type_version << "\" object_ID=\"" << object_ID << "\"></" << Item.first << ">" << std::endl;
    
    shared_ptr< std::stringstream > tmp = current_ss;
    current_ss = shared_ptr< std::stringstream >(new std::stringstream());
    
    boost::graph_traits< object_graph >::vertex_descriptor tmp_v = current_node;
    current_node = object_ID;
    
    Item.second->save(*this,type_version);
    
    // grab the resulting string.
    (*obj_graph)[object_ID].xml_src = current_ss->str();
    
    current_ss = tmp;
    current_node = tmp_v;
  } else {
    (*current_ss) << "<" << Item.first << " type_ID=\"0\" version=\"0\" object_ID=\"0\"></" << Item.first << ">" << std::endl;
  };
  
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_serializable(const serializable& Item) {
  return *this & std::pair<std::string, const serializable& >("Item",Item);
};

oarchive& RK_CALL objtree_oarchive::save_serializable(const std::pair<std::string, const serializable& >& Item) {
  archive_object_header hdr;
  const unsigned int* type_ID = Item.second.getObjectType()->TypeID_begin();
  hdr.type_version = Item.second.getObjectType()->TypeVersion();
  hdr.object_ID = 0;
  hdr.size = 0;
  hdr.is_external = false;

  (*current_ss) << "<" << Item.first << " type_ID=\""; 
  while(*type_ID) {
    (*current_ss) << *type_ID << ".";
    ++type_ID;
  };
  (*current_ss) << "0\" version=\"" << hdr.type_version << "\">" << std::endl;

  Item.second.save(*this,hdr.type_version);
  
  (*current_ss) << "</" << Item.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_char(char i) {
  return objtree_oarchive::save_char(std::pair<std::string, char >("char",i));
};

oarchive& RK_CALL objtree_oarchive::save_char(const std::pair<std::string, char >& i) {
  (*current_ss) << "<" << i.first << ">\"" << static_cast<int>(i.second) << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_char(unsigned char u) {
  return objtree_oarchive::save_unsigned_char(std::pair<std::string, unsigned char >("unsigned char",u));
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_char(const std::pair<std::string, unsigned char >& u) {
  (*current_ss) << "<" << u.first << ">\"" << static_cast<unsigned int>(u.second) << "\"</" << u.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_int(int i) {
  return objtree_oarchive::save_int(std::pair<std::string, int >("int",i));
};


oarchive& RK_CALL objtree_oarchive::save_int(const std::pair<std::string, int >& i) {
  (*current_ss) << "<" << i.first << ">\"" << i.second << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL objtree_oarchive::save_unsigned_int(unsigned int u) {
  return objtree_oarchive::save_unsigned_int(std::pair<std::string, unsigned int >("unsigned int",u));
};


oarchive& RK_CALL objtree_oarchive::save_unsigned_int(const std::pair<std::string, unsigned int >& u) {
  (*current_ss) << "<" << u.first << ">\"" << u.second << "\"</" << u.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_float(float f) {
  return objtree_oarchive::save_float(std::pair<std::string, float >("real",f));
};


oarchive& RK_CALL objtree_oarchive::save_float(const std::pair<std::string, float >& f) {
  (*current_ss) << "<" << f.first << ">\"" << f.second << "\"</" << f.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_double(double d) {
  return objtree_oarchive::save_double(std::pair<std::string, double >("real",d));
};


oarchive& RK_CALL objtree_oarchive::save_double(const std::pair<std::string, double >& d) {
  (*current_ss) << "<" << d.first << ">\"" << d.second << "\"</" << d.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_bool(bool b) {
  return objtree_oarchive::save_bool(std::pair<std::string, bool >("bool",b));
};


oarchive& RK_CALL objtree_oarchive::save_bool(const std::pair<std::string, bool >& b) {
  (*current_ss) << "<" << b.first << ">\"" << (b.second ? "true" : "false") << "\"</" << b.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL objtree_oarchive::save_string(const std::string& s) {
  return objtree_oarchive::save_string(std::pair<std::string, const std::string& >("string",s));
};


oarchive& RK_CALL objtree_oarchive::save_string(const std::pair<std::string, const std::string& >& s) {
  (*current_ss) << "<" << s.first << ">\"" << s.second << "\"</" << s.first << ">" << std::endl;
  return *this;
};






}; //serialization


}; //ReaK


