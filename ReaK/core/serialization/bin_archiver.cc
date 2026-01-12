
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

#include "ReaK/core/serialization/bin_archiver.h"

#include "ReaK/core/base/shared_object.h"
#include "ReaK/core/rtti/rtti.h"

#include "ReaK/core/serialization/archiving_exceptions.h"

#include "ReaK/core/base/endian_conversions.h"

#include <map>
#include <string>
#include <vector>

#include <fstream>

namespace ReaK::serialization {

bin_iarchive::bin_iarchive(const std::string& file_name) {

  file_stream = std::make_shared<std::ifstream>(
      file_name.c_str(), std::ios::binary | std::ios::in);

  std::string header;
  *this >> header;
  unsigned int version = 0;
  *this >> version;

  if (!(header == "reak_serialization::bin_archive")) {
    throw std::ios_base::failure("Binary Archive has a corrupt header!");
  }
  if (version != 2) {
    throw std::ios_base::failure(
        "Binary Archive is of an unknown file version!");
  }
}

bin_iarchive::bin_iarchive(std::istream& stream) {

  file_stream = std::shared_ptr<std::istream>(&stream, null_deleter());

  std::string header;
  *this >> header;
  unsigned int version = 0;
  *this >> version;

  if (!(header == "reak_serialization::bin_archive")) {
    throw std::ios_base::failure("Binary Archive has a corrupt header!");
  }
  if (version != 2) {
    throw std::ios_base::failure(
        "Binary Archive is of an unknown file version!");
  }
}

bin_iarchive::~bin_iarchive() = default;

iarchive& bin_iarchive::load_serializable_ptr(
    serializable_shared_pointer& item) {
  archive_object_header hdr;
  item = serializable_shared_pointer();

  std::vector<unsigned int> type_id_vect;
  unsigned int i = 1;
  while (i != 0) {
    std::uint64_t tmp_i = 0;
    bin_iarchive::load_unsigned_int(tmp_i);
    i = static_cast<unsigned int>(tmp_i);
    type_id_vect.push_back(i);
  }

  std::uint64_t tmp_tv = 0;
  bin_iarchive::load_unsigned_int(tmp_tv);
  hdr.type_version = static_cast<unsigned int>(tmp_tv);
  bin_iarchive::load_unsigned_int(hdr.object_ID);
  bin_iarchive::load_bool(hdr.is_external);
  bin_iarchive::load_unsigned_int(hdr.size);

  if (hdr.object_ID == 0) {
    return *this;
  }
  if ((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    item = mObjRegistry[hdr.object_ID];
    file_stream->ignore(hdr.size);
    return *this;
  }

  if (hdr.is_external) {
    std::string ext_filename;
    std::streampos start_pos = file_stream->tellg();
    *this >> ext_filename;
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      file_stream->seekg(start_pos + std::streampos(hdr.size));
    }

    bin_iarchive a(
        ext_filename);  // if this throws, let it propagate up (no point catching and throwing).
    a >> item;

    return *this;
  }

  // Find the class in question in the repository.
  rtti::so_type* p =
      rtti::so_type_repo::get_instance().find_type(type_id_vect.data());
  if ((p == nullptr) || (p->version() < hdr.type_version)) {
    file_stream->ignore(hdr.size);
    throw unsupported_type(unsupported_type::not_found_in_repo,
                           type_id_vect.data());
  }
  std::shared_ptr<shared_object> po(p->create_object());
  if (!po) {
    file_stream->ignore(hdr.size);
    throw unsupported_type(unsupported_type::could_not_create,
                           type_id_vect.data());
  }

  item = po;
  if (hdr.object_ID <
      mObjRegistry
          .size()) {  // somehow this object-ID was previously skipped over.
    mObjRegistry[hdr.object_ID] = item;
  } else if (hdr.object_ID == mObjRegistry.size()) {
    mObjRegistry.push_back(
        item);  // in theory, only this condition should occur
  } else if (hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = item;
  }

  std::streampos start_pos = file_stream->tellg();
  item->load(*this, hdr.type_version);
  std::streampos end_pos = file_stream->tellg();

  if (std::streampos(hdr.size) + start_pos != end_pos) {
    file_stream->seekg(start_pos + std::streampos(hdr.size));
  }

  return *this;
}

iarchive& bin_iarchive::load_serializable_ptr(
    const std::pair<std::string, serializable_shared_pointer&>& item) {
  return bin_iarchive::load_serializable_ptr(item.second);
}

iarchive& bin_iarchive::load_serializable(serializable& item) {
  archive_object_header hdr;

  std::vector<unsigned int> type_id_vect;
  unsigned int i = 1;
  while (i != 0) {
    *this >> i;
    type_id_vect.push_back(i);
  }

  *this >> hdr.type_version >> hdr.size;

  std::streampos start_pos = file_stream->tellg();
  item.load(*this, hdr.type_version);
  std::streampos end_pos = file_stream->tellg();

  if (std::streampos(hdr.size) + start_pos != end_pos) {
    file_stream->seekg(start_pos + std::streampos(hdr.size));
  }

  return *this;
}

iarchive& bin_iarchive::load_serializable(
    const std::pair<std::string, serializable&>& item) {
  return bin_iarchive::load_serializable(item.second);
}

iarchive& bin_iarchive::load_char(char& i) {
  file_stream->read(reinterpret_cast<char*>(&i), 1);
  return *this;
}

iarchive& bin_iarchive::load_char(const std::pair<std::string, char&>& i) {
  return bin_iarchive::load_char(i.second);
}

iarchive& bin_iarchive::load_unsigned_char(unsigned char& u) {
  file_stream->read(reinterpret_cast<char*>(&u), 1);
  return *this;
}

iarchive& bin_iarchive::load_unsigned_char(
    const std::pair<std::string, unsigned char&>& u) {
  return bin_iarchive::load_unsigned_char(u.second);
}

iarchive& bin_iarchive::load_int(std::int64_t& i) {
  file_stream->read(reinterpret_cast<char*>(&i), sizeof(std::int64_t));
  ntoh_any(i);
  return *this;
}

iarchive& bin_iarchive::load_int(
    const std::pair<std::string, std::int64_t&>& i) {
  return bin_iarchive::load_int(i.second);
}

iarchive& bin_iarchive::load_unsigned_int(std::uint64_t& u) {
  file_stream->read(reinterpret_cast<char*>(&u), sizeof(std::uint64_t));
  ntoh_any(u);
  return *this;
}

iarchive& bin_iarchive::load_unsigned_int(
    const std::pair<std::string, std::uint64_t&>& u) {
  return bin_iarchive::load_unsigned_int(u.second);
}

iarchive& bin_iarchive::load_float(float& f) {
  file_stream->read(reinterpret_cast<char*>(&f), sizeof(float));
  ntoh_any(f);
  return *this;
}

iarchive& bin_iarchive::load_float(const std::pair<std::string, float&>& f) {
  return bin_iarchive::load_float(f.second);
}

iarchive& bin_iarchive::load_double(double& d) {
  file_stream->read(reinterpret_cast<char*>(&d), sizeof(double));
  ntoh_any(d);
  return *this;
}

iarchive& bin_iarchive::load_double(const std::pair<std::string, double&>& d) {
  return bin_iarchive::load_double(d.second);
}

iarchive& bin_iarchive::load_bool(bool& b) {
  char tmp = 0;
  file_stream->read(&tmp, 1);
  b = (tmp != 0);
  return *this;
}

iarchive& bin_iarchive::load_bool(const std::pair<std::string, bool&>& b) {
  return bin_iarchive::load_bool(b.second);
}

iarchive& bin_iarchive::load_string(std::string& s) {
  std::getline(*file_stream, s, '\0');
  return *this;
}

iarchive& bin_iarchive::load_string(
    const std::pair<std::string, std::string&>& s) {
  return bin_iarchive::load_string(s.second);
}

bin_oarchive::bin_oarchive(const std::string& file_name) {
  file_stream = std::make_shared<std::ofstream>(
      file_name.c_str(), std::ios::binary | std::ios::out);

  *this << std::string("reak_serialization::bin_archive");
  unsigned int version = 2;
  *this << version;
}

bin_oarchive::bin_oarchive(std::ostream& stream) {
  file_stream = std::shared_ptr<std::ostream>(&stream, null_deleter());

  *this << std::string("reak_serialization::bin_archive");
  unsigned int version = 2;
  *this << version;
}

bin_oarchive::~bin_oarchive() = default;

oarchive& bin_oarchive::save_to_new_archive_impl(
    const serializable_shared_pointer& item, const std::string& file_name) {
  archive_object_header hdr;
  bool already_saved(false);

  if (item) {
    auto it = mObjRegMap.find(item);

    if (it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(item);
      mObjRegMap[item] = hdr.object_ID;
    }

    rtti::so_type* obj_type = item->get_object_type();
    const unsigned int* type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = true;
    hdr.size = 0;
    while (((type_id) != nullptr) && ((*type_id) != 0U)) {
      bin_oarchive::save_unsigned_int(*type_id);
      ++type_id;
    }
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  }

  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if (already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream->tellp();
    bin_oarchive::save_unsigned_int(hdr.size);

    std::streampos start_pos = file_stream->tellp();
    bin_oarchive::save_string(file_name);
    std::streampos end_pos = file_stream->tellp();
    hdr.size = end_pos - start_pos;
    file_stream->seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream->seekp(end_pos);

    bin_oarchive a(file_name);
    a << item;
  }

  return *this;
}

oarchive& bin_oarchive::save_to_new_archive_named_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& item,
    const std::string& file_name) {
  return bin_oarchive::save_to_new_archive_impl(item.second, file_name);
}

oarchive& bin_oarchive::save_serializable_ptr(
    const serializable_shared_pointer& item) {
  archive_object_header hdr;
  bool already_saved(false);

  if (item) {
    auto it = mObjRegMap.find(item);

    if (it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(item);
      mObjRegMap[item] = hdr.object_ID;
    }

    rtti::so_type* obj_type = item->get_object_type();
    const unsigned int* type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = false;
    hdr.size = 0;
    while (((type_id) != nullptr) && ((*type_id) != 0U)) {
      bin_oarchive::save_unsigned_int(*type_id);
      ++type_id;
    }
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    hdr.size = 0;
    already_saved = true;
  }

  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);
  bin_oarchive::save_unsigned_int(hdr.object_ID);
  bin_oarchive::save_bool(hdr.is_external);

  if (already_saved) {
    bin_oarchive::save_unsigned_int(hdr.size);
  } else {
    std::streampos size_pos = file_stream->tellp();
    bin_oarchive::save_unsigned_int(hdr.size);

    std::streampos start_pos = file_stream->tellp();
    item->save(*this, hdr.type_version);
    std::streampos end_pos = file_stream->tellp();

    hdr.size = end_pos - start_pos;
    file_stream->seekp(size_pos);
    bin_oarchive::save_unsigned_int(hdr.size);
    file_stream->seekp(end_pos);
  }

  return *this;
}

oarchive& bin_oarchive::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& item) {
  return bin_oarchive::save_serializable_ptr(item.second);
}

oarchive& bin_oarchive::save_serializable(const serializable& item) {
  archive_object_header hdr;

  const unsigned int* type_id = item.get_object_type()->id_begin();
  hdr.type_version = item.get_object_type()->version();
  hdr.object_ID = 0;
  hdr.is_external = false;
  hdr.size = 0;

  while (*type_id != 0U) {
    bin_oarchive::save_unsigned_int(*type_id);
    ++type_id;
  }
  bin_oarchive::save_unsigned_int(0);

  bin_oarchive::save_unsigned_int(hdr.type_version);

  std::streampos size_pos = file_stream->tellp();
  bin_oarchive::save_unsigned_int(hdr.size);

  std::streampos start_pos = file_stream->tellp();
  item.save(*this, hdr.type_version);
  std::streampos end_pos = file_stream->tellp();

  hdr.size = end_pos - start_pos;
  file_stream->seekp(size_pos);
  bin_oarchive::save_unsigned_int(hdr.size);
  file_stream->seekp(end_pos);

  return *this;
}

oarchive& bin_oarchive::save_serializable(
    const std::pair<std::string, const serializable&>& item) {
  return bin_oarchive::save_serializable(item.second);
}

oarchive& bin_oarchive::save_char(char i) {
  file_stream->write(reinterpret_cast<char*>(&i), 1);
  return *this;
}

oarchive& bin_oarchive::save_char(const std::pair<std::string, char>& i) {
  return bin_oarchive::save_char(i.second);
}

oarchive& bin_oarchive::save_unsigned_char(unsigned char u) {
  file_stream->write(reinterpret_cast<char*>(&u), 1);
  return *this;
}

oarchive& bin_oarchive::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  return bin_oarchive::save_unsigned_char(u.second);
}

oarchive& bin_oarchive::save_int(std::int64_t i) {
  hton_any(i);
  file_stream->write(reinterpret_cast<char*>(&i), sizeof(std::int64_t));
  return *this;
}

oarchive& bin_oarchive::save_int(
    const std::pair<std::string, std::int64_t>& i) {
  return bin_oarchive::save_int(i.second);
}

oarchive& bin_oarchive::save_unsigned_int(std::uint64_t u) {
  hton_any(u);
  file_stream->write(reinterpret_cast<char*>(&u), sizeof(std::uint64_t));
  return *this;
}

oarchive& bin_oarchive::save_unsigned_int(
    const std::pair<std::string, std::uint64_t>& u) {
  return bin_oarchive::save_unsigned_int(u.second);
}

oarchive& bin_oarchive::save_float(float f) {
  hton_any(f);
  file_stream->write(reinterpret_cast<char*>(&f), sizeof(float));
  return *this;
}

oarchive& bin_oarchive::save_float(const std::pair<std::string, float>& f) {
  return bin_oarchive::save_float(f.second);
}

oarchive& bin_oarchive::save_double(double d) {
  hton_any(d);
  file_stream->write(reinterpret_cast<char*>(&d), sizeof(double));
  return *this;
}

oarchive& bin_oarchive::save_double(const std::pair<std::string, double>& d) {
  return bin_oarchive::save_double(d.second);
}

oarchive& bin_oarchive::save_bool(bool b) {
  char tmp = 0;
  if (b) {
    tmp = 1;
  }
  file_stream->write(&tmp, 1);
  return *this;
}

oarchive& bin_oarchive::save_bool(const std::pair<std::string, bool>& b) {
  return bin_oarchive::save_bool(b.second);
}

oarchive& bin_oarchive::save_string(const std::string& s) {
  file_stream->write(s.c_str(), s.length() + 1);
  return *this;
}

oarchive& bin_oarchive::save_string(
    const std::pair<std::string, const std::string&>& s) {
  return bin_oarchive::save_string(s.second);
}

}  // namespace ReaK::serialization
