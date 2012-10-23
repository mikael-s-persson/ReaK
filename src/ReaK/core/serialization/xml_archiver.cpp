
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

#include "xml_archiver.hpp"

#include <iostream>
#include <iterator>

#include "base/shared_object.hpp"

#include <cstdlib>

namespace ReaK {

namespace serialization {


char xml_iarchive::getNextChar() {
  char c;
  file_stream.get(c);
  while((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r'))
    file_stream.get(c);
  return c;
};

std::string xml_iarchive::readToken() {
  std::string result;
  char c = getNextChar();
  if(c != '<')
    return result;
  c = getNextChar();
  if((c == '!') || (c == '?')) {
    char line_str[512];
    file_stream.getline(line_str,512);
    return readToken();
  };
  while(c != '>') {
    result += c;
    file_stream.get(c);
  };
  return result;
};


void xml_iarchive::skipToEndToken(const std::string& name) {
  std::string token = readToken();
  trimStr(token);
  while(token != "/" + name) {
    token = readToken();
    trimStr(token);
  };
};

void xml_iarchive::trimStr(std::string& s) {
  unsigned int i=0;
  for(;((i < s.size()) && ((s[i] == ' ') || (s[i] == '\t') || (s[i] == '\n') || (s[i] == '\r')));++i) ;
  std::string result;
  for(;((i < s.size()) && (s[i] != ' ') && (s[i] != '\t') && (s[i] != '\n') && (s[i] != '\r'));++i)
    result += s[i];
  s = result;
};

bool xml_iarchive::readNamedValue(const std::string& value_name,std::string& value_str) {
  std::string token = readToken();
  trimStr(token);
  if((value_name.empty()) || (token != value_name))
    return false;

  char c;
  file_stream.get(c);
  while(c != '\"')
    file_stream.get(c);

  value_str.clear();
  file_stream.get(c);
  while(c != '\"') {
    value_str += c;
    file_stream.get(c);
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

archive_object_header xml_iarchive::readHeader(const std::string& obj_name) {
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

  if(values["is_external"].empty())
    result.is_external = false;
  else if(values["is_external"] == "true")
    result.is_external = true;
  else
    result.is_external = false;

  return result;
};

xml_iarchive::xml_iarchive(const std::string& FileName) {
  file_stream.open(FileName.c_str());

  archive_object_header global_hdr = readHeader("reak_serialization");
  if(global_hdr.type_version != 2)
    throw std::ios_base::failure("ReaK XML Archive is of an unknown version!");

};

xml_iarchive::~xml_iarchive() {
  file_stream.close();
};


iarchive& RK_CALL xml_iarchive::load_serializable_ptr(serializable_shared_pointer& Item) {
  return xml_iarchive::load_serializable_ptr(std::pair<std::string, serializable_shared_pointer& >("Item",Item));
};

iarchive& RK_CALL xml_iarchive::load_serializable_ptr(const std::pair<std::string, serializable_shared_pointer& >& Item) {
  archive_object_header hdr;
  Item.second = serializable_shared_pointer();
  
  hdr = readHeader(Item.first);
  if((hdr.type_ID == NULL) || (hdr.type_version == 0) || (hdr.object_ID == 0)) {
    skipToEndToken(Item.first);
    delete[] hdr.type_ID;
    return *this;
  };

  if((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    Item.second = mObjRegistry[hdr.object_ID];
    skipToEndToken(Item.first);
    delete[] hdr.type_ID;
    return *this;
  };

  if(hdr.is_external) {
    std::string ext_filename;
    xml_iarchive::load_string(std::pair<std::string, std::string& >("filename", ext_filename));

    skipToEndToken(Item.first);

    try {
      xml_iarchive a(ext_filename);
      a & Item;
    } catch(std::ios_base::failure e) {
      RK_ERROR("Could not load object: " << Item.first << " - failed with message: " << e.what());
      Item.second = serializable_shared_pointer();
    };

    delete[] hdr.type_ID;
    return *this;
  };

  //Find the class in question in the repository.
  rtti::so_type::weak_pointer p( rtti::so_type_repo::getInstance().findType(hdr.type_ID) );
  delete[] hdr.type_ID;
  if((p.expired()) || (p.lock()->TypeVersion() < hdr.type_version)) {
    skipToEndToken(Item.first);
    return *this;
  };
  ReaK::shared_ptr<shared_object> po(p.lock()->CreateObject());
  if(!po) {
    skipToEndToken(Item.first);
    return *this;
  };

  Item.second = po;
  if(hdr.object_ID < mObjRegistry.size())
    mObjRegistry[hdr.object_ID] = Item.second;
  else if(hdr.object_ID == mObjRegistry.size())
    mObjRegistry.push_back(Item.second);                //in theory, only this condition should occur
  else if(hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = Item.second;
  };

  Item.second->load(*this,hdr.type_version);

  skipToEndToken(Item.first);
  return *this;
};

iarchive& RK_CALL xml_iarchive::load_serializable(serializable& Item) {
  return xml_iarchive::load_serializable(std::pair<std::string, serializable&>("Item",Item));
};

iarchive& RK_CALL xml_iarchive::load_serializable(const std::pair<std::string, serializable& >& Item) {
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

iarchive& RK_CALL xml_iarchive::load_char(char& i) {
  return xml_iarchive::load_char(std::pair<std::string, char& >("char",i));
};

iarchive& RK_CALL xml_iarchive::load_char(const std::pair<std::string, char& >& i) {
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

iarchive& RK_CALL xml_iarchive::load_unsigned_char(unsigned char& u) {
  return xml_iarchive::load_unsigned_char(std::pair<std::string, unsigned char& >("unsigned char",u));
};

iarchive& RK_CALL xml_iarchive::load_unsigned_char(const std::pair<std::string, unsigned char& >& u) {
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

iarchive& RK_CALL xml_iarchive::load_int(int& i) {
  return xml_iarchive::load_int(std::pair<std::string, int& >("int",i));
};

iarchive& RK_CALL xml_iarchive::load_int(const std::pair<std::string, int& >& i) {
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

iarchive& RK_CALL xml_iarchive::load_unsigned_int(unsigned int& u) {
  return xml_iarchive::load_unsigned_int(std::pair<std::string, unsigned int& >("unsigned int",u));
};

iarchive& RK_CALL xml_iarchive::load_unsigned_int(const std::pair<std::string, unsigned int& >& u) {
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

iarchive& RK_CALL xml_iarchive::load_float(float& f) {
  return xml_iarchive::load_float(std::pair<std::string, float& >("real",f));
};

iarchive& RK_CALL xml_iarchive::load_float(const std::pair<std::string, float& >& f) {
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

iarchive& RK_CALL xml_iarchive::load_double(double& d) {
  return xml_iarchive::load_double(std::pair<std::string, double& >("real",d));
};

iarchive& RK_CALL xml_iarchive::load_double(const std::pair<std::string, double& >& d) {
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

iarchive& RK_CALL xml_iarchive::load_bool(bool& b) {
  return xml_iarchive::load_bool(std::pair<std::string, bool& >("bool",b));
};

iarchive& RK_CALL xml_iarchive::load_bool(const std::pair<std::string, bool& >& b) {
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

iarchive& RK_CALL xml_iarchive::load_string(std::string& s) {
  return xml_iarchive::load_string(std::pair<std::string, std::string& >("string",s));
};

iarchive& RK_CALL xml_iarchive::load_string(const std::pair<std::string, std::string& >& s) {
  readNamedValue(s.first,s.second);
  return *this;
};











xml_oarchive::xml_oarchive(const std::string& FileName) {
  file_stream.open(FileName.c_str());

  file_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>" << std::endl;
  file_stream << "<!DOCTYPE reak_serialization>" << std::endl;
  file_stream << "<reak_serialization version=\"2\">" << std::endl;
  tabulation = 0;
};

xml_oarchive::~xml_oarchive() {

  file_stream << "</reak_serialization>" << std::endl;

  file_stream.close();
};

void xml_oarchive::addTabulations() {
  for(unsigned int i=0;i<tabulation;++i) {
    file_stream << "    ";
  };
};


oarchive& RK_CALL xml_oarchive::saveToNewArchive_impl(const serializable_shared_pointer& Item, const std::string& FileName) {
  return xml_oarchive::saveToNewArchiveNamed_impl(std::pair<std::string, const serializable_shared_pointer& >("Item",Item),FileName);
};

oarchive& RK_CALL xml_oarchive::saveToNewArchiveNamed_impl(const std::pair<std::string, const serializable_shared_pointer& >& Item, const std::string& FileName) {
  archive_object_header hdr;
  bool already_saved(false);

  if(Item.second) {
    std::map< serializable_shared_pointer, unsigned int>::const_iterator it = mObjRegMap.find(Item.second);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item.second);
      mObjRegMap[Item.second] = hdr.object_ID;
    };

    rtti::so_type::shared_pointer obj_type = Item.second->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = true;
    hdr.size = 0;

    addTabulations();
    file_stream << "<" << Item.first << " type_ID=\"";
    while(*type_ID) {
      file_stream << *type_ID << ".";
      ++type_ID;
    };
    file_stream << "0\" version=\"" << hdr.type_version
	        << "\" object_ID=\"" << hdr.object_ID
	        << "\" is_external=\"true\">" << std::endl;
  } else {
    already_saved = true;

    addTabulations();
    file_stream << "<" << Item.first << " type_ID=\"0\" version=\"0\" object_ID=\"0\" is_external=\"false\">" << std::endl;
  };

  if(!already_saved) {
    tabulation++;
    addTabulations();
    file_stream << "<filename>\"" << FileName << "\"</filename>" << std::endl;
    tabulation--;

    xml_oarchive a(FileName);
    a & Item;
  };

  addTabulations();
  file_stream << "</" << Item.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_serializable_ptr(const serializable_shared_pointer& Item) {
  return *this & std::pair<std::string, const serializable_shared_pointer& >("Item",Item);
};


oarchive& RK_CALL xml_oarchive::save_serializable_ptr(const std::pair<std::string, const serializable_shared_pointer& >& Item) {

  archive_object_header hdr;
  bool already_saved(false);

  if(Item.second) {
    std::map< serializable_shared_pointer, unsigned int>::const_iterator it = mObjRegMap.find(Item.second);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item.second);
      mObjRegMap[Item.second] = hdr.object_ID;
    };

    rtti::so_type::shared_pointer obj_type = Item.second->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = false;
    hdr.size = 0;

    addTabulations();
    file_stream << "<" << Item.first << " type_ID=\"";
    while(*type_ID) {
      file_stream << *type_ID << "."; 
      ++type_ID;
    };
    file_stream << "0\" version=\"" << hdr.type_version
	        << "\" object_ID=\"" << hdr.object_ID
	        << "\" is_external=\"false\">" << std::endl;
  } else {
    already_saved = true;

    addTabulations();
    file_stream << "<" << Item.first << " type_ID=\"0\" version=\"0\" object_ID=\"0\" is_external=\"false\">" << std::endl;
  };

  if(!already_saved) {

    tabulation++;
    Item.second->save(*this,hdr.type_version);
    tabulation--;

  };

  addTabulations();
  file_stream << "</" << Item.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_serializable(const serializable& Item) {
  return *this & std::pair<std::string, const serializable& >("Item",Item);
};

oarchive& RK_CALL xml_oarchive::save_serializable(const std::pair<std::string, const serializable& >& Item) {
  archive_object_header hdr;
  const unsigned int* type_ID = Item.second.getObjectType()->TypeID_begin();
  hdr.type_version = Item.second.getObjectType()->TypeVersion();
  hdr.object_ID = 0;
  hdr.size = 0;
  hdr.is_external = false;

  addTabulations();
  file_stream << "<" << Item.first << " type_ID=\""; 
  while(*type_ID) {
    file_stream << *type_ID << ".";
    ++type_ID;
  };
  file_stream << "0\" version=\"" << hdr.type_version << "\">" << std::endl;

  tabulation++;
  Item.second.save(*this,hdr.type_version);
  tabulation--;

  addTabulations();
  file_stream << "</" << Item.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_char(char i) {
  return xml_oarchive::save_char(std::pair<std::string, char >("char",i));
};

oarchive& RK_CALL xml_oarchive::save_char(const std::pair<std::string, char >& i) {
  addTabulations();
  file_stream << "<" << i.first << ">\"" << static_cast<int>(i.second) << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL xml_oarchive::save_unsigned_char(unsigned char u) {
  return xml_oarchive::save_unsigned_char(std::pair<std::string, unsigned char >("unsigned char",u));
};

oarchive& RK_CALL xml_oarchive::save_unsigned_char(const std::pair<std::string, unsigned char >& u) {
  addTabulations();
  file_stream << "<" << u.first << ">\"" << static_cast<unsigned int>(u.second) << "\"</" << u.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL xml_oarchive::save_int(int i) {
  return xml_oarchive::save_int(std::pair<std::string, int >("int",i));
};


oarchive& RK_CALL xml_oarchive::save_int(const std::pair<std::string, int >& i) {
  addTabulations();
  file_stream << "<" << i.first << ">\"" << i.second << "\"</" << i.first << ">" << std::endl;
  return *this;
};

oarchive& RK_CALL xml_oarchive::save_unsigned_int(unsigned int u) {
  return xml_oarchive::save_unsigned_int(std::pair<std::string, unsigned int >("unsigned int",u));
};


oarchive& RK_CALL xml_oarchive::save_unsigned_int(const std::pair<std::string, unsigned int >& u) {
  addTabulations();
  file_stream << "<" << u.first << ">\"" << u.second << "\"</" << u.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_float(float f) {
  return xml_oarchive::save_float(std::pair<std::string, float >("real",f));
};


oarchive& RK_CALL xml_oarchive::save_float(const std::pair<std::string, float >& f) {
  addTabulations();
  file_stream << "<" << f.first << ">\"" << f.second << "\"</" << f.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_double(double d) {
  return xml_oarchive::save_double(std::pair<std::string, double >("real",d));
};


oarchive& RK_CALL xml_oarchive::save_double(const std::pair<std::string, double >& d) {
  addTabulations();
  file_stream << "<" << d.first << ">\"" << d.second << "\"</" << d.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_bool(bool b) {
  return xml_oarchive::save_bool(std::pair<std::string, bool >("bool",b));
};


oarchive& RK_CALL xml_oarchive::save_bool(const std::pair<std::string, bool >& b) {
  addTabulations();
  file_stream << "<" << b.first << ">\"" << (b.second ? "true" : "false") << "\"</" << b.first << ">" << std::endl;
  return *this;
};


oarchive& RK_CALL xml_oarchive::save_string(const std::string& s) {
  return xml_oarchive::save_string(std::pair<std::string, const std::string& >("string",s));
};


oarchive& RK_CALL xml_oarchive::save_string(const std::pair<std::string, const std::string& >& s) {
  addTabulations();
  file_stream << "<" << s.first << ">\"" << s.second << "\"</" << s.first << ">" << std::endl;
  return *this;
};



}; //serialization


}; //ReaK

