
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

#include "bin_archiver.hpp"


#include "base/shared_object.hpp"

#include <string>
#include <fstream>
#include <iomanip>

#include <algorithm>

namespace ReaK {

namespace serialization {



bin_iarchive::bin_iarchive(const std::string& FileName) {

  file_stream.open(FileName.c_str(),std::ios::binary | std::ios::in);

  std::string header;
  *this >> header;
  unsigned int version;
  *this >> version;

  if(!(header == "reak_serialization::bin_archive"))
    throw std::ios_base::failure("Binary Archive has a corrupt header!");
  if(version != 2)
    throw std::ios_base::failure("Binary Archive is of an unknown file version!");

};

bin_iarchive::~bin_iarchive() {

  file_stream.close();
};

iarchive& RK_CALL bin_iarchive::load_serializable_ptr(boost::shared_ptr<serializable>& Item) {
  archive_object_header hdr;

  std::vector<unsigned int> typeIDvect;
  unsigned int i;
  do {
    *this >> i;
    typeIDvect.push_back(i);
  } while(i != 0);
  
  *this >> hdr.type_version >> hdr.object_ID >> hdr.is_external >> hdr.size;

  if((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    Item = mObjRegistry[hdr.object_ID];
    file_stream.ignore(hdr.size);
    return *this;
  };

  if(hdr.is_external) {
    std::string ext_filename;
    std::streampos start_pos = file_stream.tellg();
    *this >> ext_filename;
    std::streampos end_pos = file_stream.tellg();
    if (hdr.size + start_pos != end_pos)
      file_stream.seekg(start_pos + std::streampos(hdr.size));

    bin_iarchive a(ext_filename);
    a >> Item;

    return *this;
  };

  hdr.type_ID = new unsigned int[typeIDvect.size()];
  std::copy(typeIDvect.begin(),typeIDvect.end(),hdr.type_ID);
  //Find the class in question in the repository.
  boost::weak_ptr<rtti::so_type> p( rtti::so_type_repo::getInstance().findType(hdr.type_ID) );
  delete[] hdr.type_ID;
  if((p.expired()) || (p.lock()->TypeVersion() < hdr.type_version)) {
    file_stream.ignore(hdr.size);
    Item = boost::shared_ptr<serializable>();
    return *this;
  };
  boost::shared_ptr<shared_object> po(p.lock()->CreateObject());
  if(!po) {
    file_stream.ignore(hdr.size);
    Item = boost::shared_ptr<serializable>();
    return *this;
  };

  Item = po;
  if(hdr.object_ID < mObjRegistry.size())
    mObjRegistry[hdr.object_ID] = Item;
  else if(hdr.object_ID == mObjRegistry.size())
    mObjRegistry.push_back(Item);                //in theory, only this condition should occur
  else if(hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = Item;
  };

  std::streampos start_pos = file_stream.tellg();
  Item->load(*this,hdr.type_version);
  std::streampos end_pos = file_stream.tellg();

  if (hdr.size + start_pos != end_pos)
    file_stream.seekg(start_pos + std::streampos(hdr.size));

  return *this;
};

iarchive& RK_CALL bin_iarchive::load_serializable_ptr(const std::pair<std::string, boost::shared_ptr<serializable>& >& Item) {
  return bin_iarchive::load_serializable_ptr(Item.second);
};

iarchive& RK_CALL bin_iarchive::load_serializable(serializable& Item) {
  archive_object_header hdr;
  
  std::vector<unsigned int> typeIDvect;
  unsigned int i;
  do {
    *this >> i;
    typeIDvect.push_back(i);
  } while(i != 0);
  
  *this >> hdr.type_version >> hdr.size;

  std::streampos start_pos = file_stream.tellg();
  Item.load(*this,hdr.type_version);
  std::streampos end_pos = file_stream.tellg();

  if (hdr.size + start_pos != end_pos)
    file_stream.seekg(start_pos + std::streampos(hdr.size));

  return *this;
};

iarchive& RK_CALL bin_iarchive::load_serializable(const std::pair<std::string, serializable& >& Item) {
  return bin_iarchive::load_serializable(Item.second);
};

iarchive& RK_CALL bin_iarchive::load_int(int& i) {
  file_stream.read(reinterpret_cast<char*>(&i),sizeof(int));
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_int(const std::pair<std::string, int& >& i) {
  return bin_iarchive::load_int(i.second);
};

iarchive& RK_CALL bin_iarchive::load_unsigned_int(unsigned int& u) {
  file_stream.read(reinterpret_cast<char*>(&u),sizeof(unsigned int));
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_unsigned_int(const std::pair<std::string, unsigned int& >& u) {
  return bin_iarchive::load_unsigned_int(u.second);
};

iarchive& RK_CALL bin_iarchive::load_float(float& f) {
  file_stream.read(reinterpret_cast<char*>(&f),sizeof(float));
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_float(const std::pair<std::string, float& >& f) {
  return bin_iarchive::load_float(f.second);
};

iarchive& RK_CALL bin_iarchive::load_double(double& d) {
  file_stream.read(reinterpret_cast<char*>(&d),sizeof(double));
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_double(const std::pair<std::string, double& >& d) {
  return bin_iarchive::load_double(d.second);
};

iarchive& RK_CALL bin_iarchive::load_bool(bool& b) {
  file_stream.read(reinterpret_cast<char*>(&b),sizeof(bool));
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_bool(const std::pair<std::string, bool& >& b) {
  return bin_iarchive::load_bool(b.second);
};

iarchive& RK_CALL bin_iarchive::load_string(std::string& s) {
  std::getline(file_stream,s,'\0');
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_string(const std::pair<std::string, std::string& >& s) {
  return bin_iarchive::load_string(s.second);
};











bin_oarchive::bin_oarchive(const std::string& FileName) {
  file_stream.open(FileName.c_str(),std::ios::binary | std::ios::out);

  char header[] = "reak_serialization::bin_archive";
  file_stream.write(header,std::strlen(header)+1);
  unsigned int version = 2;
  file_stream.write(reinterpret_cast<char*>(&version),sizeof(unsigned int));


};

bin_oarchive::~bin_oarchive() {

  file_stream.close();
};

oarchive& RK_CALL bin_oarchive::saveToNewArchive_impl(const boost::shared_ptr<serializable>& Item, const std::string& FileName) {
  archive_object_header hdr;
  bool already_saved(false);
  const unsigned int* type_ID = NULL;

  if(Item) {
    std::map< boost::shared_ptr<serializable>, unsigned int>::const_iterator it = mObjRegMap.find(Item);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item);
      mObjRegMap[Item] = hdr.object_ID;
    };

    boost::shared_ptr<rtti::so_type> obj_type = Item->getObjectType();
    type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = true;
    hdr.size = 0;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  };

  while((type_ID) && (*type_ID)) {
    bin_oarchive::save_unsigned_int(*type_ID);
    ++type_ID;
  };
  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if(already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream.tellp();
    bin_oarchive::save_unsigned_int(hdr.size);
    bin_oarchive::save_string(FileName);
    std::streampos end_pos = file_stream.tellp();
    typedef unsigned int tmp_uint;
    hdr.size = tmp_uint(end_pos - size_pos - sizeof(tmp_uint));
    file_stream.seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream.seekp(end_pos);

    bin_oarchive a(FileName);
    a << Item;
  };

  return *this;
};

oarchive& RK_CALL bin_oarchive::saveToNewArchiveNamed_impl(const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item, const std::string& FileName) {
  return bin_oarchive::saveToNewArchive_impl(Item.second,FileName);
};


oarchive& RK_CALL bin_oarchive::save_serializable_ptr(const boost::shared_ptr<serializable>& Item) {
  archive_object_header hdr;
  bool already_saved(false);
  const unsigned int* type_ID = NULL;

  if(Item) {
    std::map< boost::shared_ptr<serializable>, unsigned int>::const_iterator it = mObjRegMap.find(Item);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item);
      mObjRegMap[Item] = hdr.object_ID;
    };

    boost::shared_ptr<rtti::so_type> obj_type = Item->getObjectType();
    type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = false;
    hdr.size = 0;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  };
  
  while((type_ID) && (*type_ID)) {
    bin_oarchive::save_unsigned_int(*type_ID);
    ++type_ID;
  };
  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if(already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream.tellp();
    bin_oarchive::save_unsigned_int(hdr.size);

    Item->save(*this,hdr.type_version);

    std::streampos end_pos = file_stream.tellp();
    typedef unsigned int tmp_uint;
    hdr.size = tmp_uint(end_pos - size_pos - sizeof(tmp_uint));
    file_stream.seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream.seekp(end_pos);
  };

  return *this;
};


oarchive& RK_CALL bin_oarchive::save_serializable_ptr(const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item) {
  return bin_oarchive::save_serializable_ptr(Item.second);
};

oarchive& RK_CALL bin_oarchive::save_serializable(const serializable& Item) {
  archive_object_header hdr;

  const unsigned int* type_ID = Item.getObjectType()->TypeID_begin();
  hdr.type_version = Item.getObjectType()->TypeVersion();
  hdr.object_ID = 0;
  hdr.is_external = false;
  hdr.size = 0;
  
  while(*type_ID) {
    bin_oarchive::save_unsigned_int(*type_ID);
    ++type_ID;
  };
  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);

  std::streampos size_pos = file_stream.tellp();
  bin_oarchive::save_unsigned_int(hdr.size);

  Item.save(*this,hdr.type_version);

  std::streampos end_pos = file_stream.tellp();
  typedef unsigned int tmp_uint;
  hdr.size = tmp_uint(end_pos - size_pos - sizeof(tmp_uint));
  file_stream.seekp(size_pos);
  bin_oarchive::save_unsigned_int(hdr.size);
  file_stream.seekp(end_pos);

  return *this;
};

oarchive& RK_CALL bin_oarchive::save_serializable(const std::pair<std::string, const serializable& >& Item) {
  return bin_oarchive::save_serializable(Item.second);
};

oarchive& RK_CALL bin_oarchive::save_int(int i) {
  file_stream.write(reinterpret_cast<char*>(&i),sizeof(int));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_int(const std::pair<std::string, int >& i) {
  return bin_oarchive::save_int(i.second);
};

oarchive& RK_CALL bin_oarchive::save_unsigned_int(unsigned int u) {
  file_stream.write(reinterpret_cast<char*>(&u),sizeof(unsigned int));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_unsigned_int(const std::pair<std::string, unsigned int >& u) {
  return bin_oarchive::save_unsigned_int(u.second);
};


oarchive& RK_CALL bin_oarchive::save_float(float f) {
  file_stream.write(reinterpret_cast<char*>(&f),sizeof(float));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_float(const std::pair<std::string, float >& f) {
  return bin_oarchive::save_float(f.second);
};


oarchive& RK_CALL bin_oarchive::save_double(double d) {
  file_stream.write(reinterpret_cast<char*>(&d),sizeof(double));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_double(const std::pair<std::string, double >& d) {
  return bin_oarchive::save_double(d.second);
};


oarchive& RK_CALL bin_oarchive::save_bool(bool b) {
  file_stream.write(reinterpret_cast<char*>(&b),sizeof(bool));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_bool(const std::pair<std::string, bool >& b) {
  return bin_oarchive::save_bool(b.second);
};


oarchive& RK_CALL bin_oarchive::save_string(const std::string& s) {
  file_stream.write(s.c_str(),std::strlen(s.c_str())+1);
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_string(const std::pair<std::string, const std::string& >& s) {
  return bin_oarchive::save_string(s.second);
};



}; //serialization


}; //ReaK


