
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

#include <ReaK/core/serialization/bin_archiver.hpp>


#include <ReaK/core/base/shared_object.hpp>
#include <ReaK/core/rtti/so_type.hpp>
#include <ReaK/core/rtti/so_type_repo.hpp>

#include <ReaK/core/serialization/archiving_exceptions.hpp>

#include <ReaK/core/base/endian_conversions.hpp>


#include <string>
#include <map>
#include <vector>

#include <fstream>

namespace ReaK {

namespace serialization {


bin_iarchive::bin_iarchive(const std::string& FileName) {
  
  file_stream = shared_ptr< std::istream >(new std::ifstream(FileName.c_str(), std::ios::binary | std::ios::in));
  
  std::string header;
  *this >> header;
  unsigned int version;
  *this >> version;

  if(!(header == "reak_serialization::bin_archive"))
    throw std::ios_base::failure("Binary Archive has a corrupt header!");
  if(version != 2)
    throw std::ios_base::failure("Binary Archive is of an unknown file version!");

};

bin_iarchive::bin_iarchive(std::istream& aStream) {
  
  file_stream = shared_ptr< std::istream >(&aStream, null_deleter());
  
  std::string header;
  *this >> header;
  unsigned int version;
  *this >> version;

  if(!(header == "reak_serialization::bin_archive"))
    throw std::ios_base::failure("Binary Archive has a corrupt header!");
  if(version != 2)
    throw std::ios_base::failure("Binary Archive is of an unknown file version!");

};


bin_iarchive::~bin_iarchive() {};



iarchive& RK_CALL bin_iarchive::load_serializable_ptr(serializable_shared_pointer& Item) {
  archive_object_header hdr;
  Item = serializable_shared_pointer();

  std::vector<unsigned int> typeIDvect;
  unsigned int i;
  do {
    bin_iarchive::load_unsigned_int(i);
    typeIDvect.push_back(i);
  } while(i != 0);
  
  bin_iarchive::load_unsigned_int(hdr.type_version);
  bin_iarchive::load_unsigned_int(hdr.object_ID);
  bin_iarchive::load_bool(hdr.is_external);
  bin_iarchive::load_unsigned_int(hdr.size);
  
  if(hdr.object_ID == 0) {
    return *this;
  };
  if((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    Item = mObjRegistry[hdr.object_ID];
    file_stream->ignore(hdr.size);
    return *this;
  };

  if(hdr.is_external) {
    std::string ext_filename;
    std::streampos start_pos = file_stream->tellg();
    *this >> ext_filename;
    std::streampos end_pos = file_stream->tellg();
    if (hdr.size + start_pos != end_pos)
      file_stream->seekg(start_pos + std::streampos(hdr.size));

    bin_iarchive a(ext_filename);  // if this throws, let it propagate up (no point catching and throwing).
    a >> Item;

    return *this;
  };
  
  //Find the class in question in the repository.
  rtti::so_type::weak_pointer p( rtti::so_type_repo::getInstance().findType(&(typeIDvect[0])) );
  if((p.expired()) || (p.lock()->TypeVersion() < hdr.type_version)) {
    file_stream->ignore(hdr.size);
    throw unsupported_type(unsupported_type::not_found_in_repo, &(typeIDvect[0]));
  };
  ReaK::shared_ptr<shared_object> po(p.lock()->CreateObject());
  if(!po) {
    file_stream->ignore(hdr.size);
    throw unsupported_type(unsupported_type::could_not_create, &(typeIDvect[0]));
  };

  Item = po;
  if(hdr.object_ID < mObjRegistry.size())  // somehow this object-ID was previously skipped over.
    mObjRegistry[hdr.object_ID] = Item;
  else if(hdr.object_ID == mObjRegistry.size())
    mObjRegistry.push_back(Item);                //in theory, only this condition should occur
  else if(hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = Item;
  };

  std::streampos start_pos = file_stream->tellg();
  Item->load(*this,hdr.type_version);
  std::streampos end_pos = file_stream->tellg();

  if (hdr.size + start_pos != end_pos)
    file_stream->seekg(start_pos + std::streampos(hdr.size));

  return *this;
};

iarchive& RK_CALL bin_iarchive::load_serializable_ptr(const std::pair<std::string, serializable_shared_pointer& >& Item) {
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

  std::streampos start_pos = file_stream->tellg();
  Item.load(*this,hdr.type_version);
  std::streampos end_pos = file_stream->tellg();

  if (hdr.size + start_pos != end_pos)
    file_stream->seekg(start_pos + std::streampos(hdr.size));

  return *this;
};

iarchive& RK_CALL bin_iarchive::load_serializable(const std::pair<std::string, serializable& >& Item) {
  return bin_iarchive::load_serializable(Item.second);
};

iarchive& RK_CALL bin_iarchive::load_char(char& i) {
  file_stream->read(reinterpret_cast<char*>(&i),1);
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_char(const std::pair<std::string, char& >& i) {
  return bin_iarchive::load_char(i.second);
};

iarchive& RK_CALL bin_iarchive::load_unsigned_char(unsigned char& u) {
  file_stream->read(reinterpret_cast<char*>(&u),1);
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_unsigned_char(const std::pair<std::string, unsigned char& >& u) {
  return bin_iarchive::load_unsigned_char(u.second);
};

iarchive& RK_CALL bin_iarchive::load_int(int& i) {
  llong_to_ulong tmp; 
  file_stream->read(reinterpret_cast<char*>(&tmp),sizeof(llong_to_ulong));
  ntoh_2ui32(tmp);
  i = static_cast<int>(tmp.i64);
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_int(const std::pair<std::string, int& >& i) {
  return bin_iarchive::load_int(i.second);
};

iarchive& RK_CALL bin_iarchive::load_unsigned_int(unsigned int& u) {
  llong_to_ulong tmp; 
  file_stream->read(reinterpret_cast<char*>(&tmp),sizeof(llong_to_ulong));
  ntoh_2ui32(tmp);
  u = static_cast<unsigned int>(tmp.ui64);
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_unsigned_int(const std::pair<std::string, unsigned int& >& u) {
  return bin_iarchive::load_unsigned_int(u.second);
};

iarchive& RK_CALL bin_iarchive::load_float(float& f) {
  float_to_ulong tmp; 
  file_stream->read(reinterpret_cast<char*>(&tmp),sizeof(float_to_ulong));
  ntoh_1ui32(tmp);
  f = tmp.f;
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_float(const std::pair<std::string, float& >& f) {
  return bin_iarchive::load_float(f.second);
};

iarchive& RK_CALL bin_iarchive::load_double(double& d) {
  double_to_ulong tmp; 
  file_stream->read(reinterpret_cast<char*>(&tmp),sizeof(double_to_ulong));
  ntoh_2ui32(tmp);
  d = tmp.d;
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_double(const std::pair<std::string, double& >& d) {
  return bin_iarchive::load_double(d.second);
};

iarchive& RK_CALL bin_iarchive::load_bool(bool& b) {
  char tmp = 0;
  file_stream->read(&tmp,1);
  b = (tmp ? true : false);
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_bool(const std::pair<std::string, bool& >& b) {
  return bin_iarchive::load_bool(b.second);
};

iarchive& RK_CALL bin_iarchive::load_string(std::string& s) {
  std::getline(*file_stream,s,'\0');
  return *this;
};

iarchive& RK_CALL bin_iarchive::load_string(const std::pair<std::string, std::string& >& s) {
  return bin_iarchive::load_string(s.second);
};











bin_oarchive::bin_oarchive(const std::string& FileName) {
  
  file_stream = shared_ptr< std::ostream >(new std::ofstream(FileName.c_str(), std::ios::binary | std::ios::out));
  
  *this << std::string("reak_serialization::bin_archive");
  unsigned int version = 2;
  *this << version;
  
};

bin_oarchive::bin_oarchive(std::ostream& aStream) {
  
  file_stream = shared_ptr< std::ostream >(&aStream, null_deleter());
  
  *this << std::string("reak_serialization::bin_archive");
  unsigned int version = 2;
  *this << version;
  
};

bin_oarchive::~bin_oarchive() { };

oarchive& RK_CALL bin_oarchive::saveToNewArchive_impl(const serializable_shared_pointer& Item, const std::string& FileName) {
  archive_object_header hdr;
  bool already_saved(false);

  if(Item) {
    std::map< serializable_shared_pointer, unsigned int>::const_iterator it = mObjRegMap.find(Item);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item);
      mObjRegMap[Item] = hdr.object_ID;
    };

    rtti::so_type::shared_pointer obj_type = Item->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = true;
    hdr.size = 0;
    while((type_ID) && (*type_ID)) {
      bin_oarchive::save_unsigned_int(*type_ID);
      ++type_ID;
    };
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  };

  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if(already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream->tellp();
    bin_oarchive::save_unsigned_int(hdr.size);
    
    std::streampos start_pos = file_stream->tellp();
    bin_oarchive::save_string(FileName);
    std::streampos end_pos = file_stream->tellp();
    typedef unsigned int tmp_uint;
    hdr.size = tmp_uint(end_pos - start_pos);
    file_stream->seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream->seekp(end_pos);
    
    bin_oarchive a(FileName);
    a << Item;
  };

  return *this;
};

oarchive& RK_CALL bin_oarchive::saveToNewArchiveNamed_impl(const std::pair<std::string, const serializable_shared_pointer& >& Item, const std::string& FileName) {
  return bin_oarchive::saveToNewArchive_impl(Item.second,FileName);
};


oarchive& RK_CALL bin_oarchive::save_serializable_ptr(const serializable_shared_pointer& Item) {
  archive_object_header hdr;
  bool already_saved(false);

  if(Item) {
    std::map< serializable_shared_pointer, unsigned int>::const_iterator it = mObjRegMap.find(Item);

    if(it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(Item);
      mObjRegMap[Item] = hdr.object_ID;
    };

    rtti::so_type::shared_pointer obj_type = Item->getObjectType();
    const unsigned int* type_ID = obj_type->TypeID_begin();
    hdr.type_version = obj_type->TypeVersion();
    hdr.is_external = false;
    hdr.size = 0;
    while((type_ID) && (*type_ID)) {
      bin_oarchive::save_unsigned_int(*type_ID);
      ++type_ID;
    };
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  };
  
  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if(already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream->tellp();
    bin_oarchive::save_unsigned_int(hdr.size);

    std::streampos start_pos = file_stream->tellp();
    Item->save(*this,hdr.type_version);
    std::streampos end_pos = file_stream->tellp();

    typedef unsigned int tmp_uint;
    hdr.size = tmp_uint(end_pos - start_pos);
    file_stream->seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream->seekp(end_pos);
  };

  return *this;
};


oarchive& RK_CALL bin_oarchive::save_serializable_ptr(const std::pair<std::string, const serializable_shared_pointer& >& Item) {
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
  
  std::streampos size_pos = file_stream->tellp();
  bin_oarchive::save_unsigned_int(hdr.size);
  
  std::streampos start_pos = file_stream->tellp();
  Item.save(*this,hdr.type_version);
  std::streampos end_pos = file_stream->tellp();

  typedef unsigned int tmp_uint;
  hdr.size = tmp_uint(end_pos - start_pos);
  file_stream->seekp(size_pos);
  bin_oarchive::save_unsigned_int(hdr.size);
  file_stream->seekp(end_pos);

  return *this;
};

oarchive& RK_CALL bin_oarchive::save_serializable(const std::pair<std::string, const serializable& >& Item) {
  return bin_oarchive::save_serializable(Item.second);
};

oarchive& RK_CALL bin_oarchive::save_char(char i) {
  file_stream->write(reinterpret_cast<char*>(&i),1);
  return *this;
};

oarchive& RK_CALL bin_oarchive::save_char(const std::pair<std::string, char >& i) {
  return bin_oarchive::save_char(i.second);
};

oarchive& RK_CALL bin_oarchive::save_unsigned_char(unsigned char u) {
  file_stream->write(reinterpret_cast<char*>(&u),1);
  return *this;
};

oarchive& RK_CALL bin_oarchive::save_unsigned_char(const std::pair<std::string, unsigned char >& u) {
  return bin_oarchive::save_unsigned_char(u.second);
};

oarchive& RK_CALL bin_oarchive::save_int(int i) {
  llong_to_ulong tmp; tmp.i64 = i;
  hton_2ui32(tmp);
  file_stream->write(reinterpret_cast<char*>(&tmp),sizeof(llong_to_ulong));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_int(const std::pair<std::string, int >& i) {
  return bin_oarchive::save_int(i.second);
};

oarchive& RK_CALL bin_oarchive::save_unsigned_int(unsigned int u) {
  llong_to_ulong tmp; tmp.ui64 = u;
  hton_2ui32(tmp);
  file_stream->write(reinterpret_cast<char*>(&tmp),sizeof(llong_to_ulong));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_unsigned_int(const std::pair<std::string, unsigned int >& u) {
  return bin_oarchive::save_unsigned_int(u.second);
};


oarchive& RK_CALL bin_oarchive::save_float(float f) {
  float_to_ulong tmp = { f };
  hton_1ui32(tmp);
  file_stream->write(reinterpret_cast<char*>(&tmp),sizeof(float_to_ulong));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_float(const std::pair<std::string, float >& f) {
  return bin_oarchive::save_float(f.second);
};


oarchive& RK_CALL bin_oarchive::save_double(double d) {
  double_to_ulong tmp = { d };
  hton_2ui32(tmp);
  file_stream->write(reinterpret_cast<char*>(&tmp),sizeof(double_to_ulong));
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_double(const std::pair<std::string, double >& d) {
  return bin_oarchive::save_double(d.second);
};


oarchive& RK_CALL bin_oarchive::save_bool(bool b) {
  char tmp = 0;
  if(b) tmp = 1;
  file_stream->write(&tmp,1);
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_bool(const std::pair<std::string, bool >& b) {
  return bin_oarchive::save_bool(b.second);
};


oarchive& RK_CALL bin_oarchive::save_string(const std::string& s) {
  file_stream->write(s.c_str(), s.length() + 1 );
  return *this;
};


oarchive& RK_CALL bin_oarchive::save_string(const std::pair<std::string, const std::string& >& s) {
  return bin_oarchive::save_string(s.second);
};



}; //serialization


}; //ReaK


