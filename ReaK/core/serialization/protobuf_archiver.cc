
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

#include "ReaK/core/serialization/protobuf_archiver.h"

#include "ReaK/core/base/endian_conversions.h"
#include "ReaK/core/base/shared_object.h"
#include "ReaK/core/rtti/rtti.h"

#include "ReaK/core/serialization/archiving_exceptions.h"

#include <bit>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <string>

namespace ReaK::serialization {

namespace {

const char* const bad_stream_msg =
    "Protobuf archive could not be loaded due to an unexpected failure or bad "
    "state of the input stream!";
const char* const unexpected_eof_msg =
    "Protobuf archive is corrupt! Unexpectingly reached end-of-file!";
const char* const corrupt_header_msg = "Protobuf archive has a corrupt header!";
const char* const unknown_version_msg =
    "Protobuf archive is of an unknown file version!";
const char* const bad_out_stream_msg =
    "Protobuf archive could not be saved due to an unexpected failure or bad "
    "state of the output stream!";

}  // namespace

protobuf_iarchive::protobuf_iarchive(const std::string& file_name) {

  file_stream = std::make_shared<std::ifstream>(
      file_name.c_str(), std::ios::binary | std::ios::in);

  std::string header;
  protobuf_iarchive::load_string(header);
  std::uint64_t version = 0;
  protobuf_iarchive::load_unsigned_int(version);

  if (!(header == "reak_serialization::protobuf_archive")) {
    throw std::ios_base::failure(corrupt_header_msg);
  }
  if (version != 2) {
    throw std::ios_base::failure(unknown_version_msg);
  }
};

protobuf_iarchive::protobuf_iarchive(std::istream& stream) {

  file_stream = std::shared_ptr<std::istream>(&stream, null_deleter());

  std::string header;
  protobuf_iarchive::load_string(header);
  std::uint64_t version = 0;
  protobuf_iarchive::load_unsigned_int(version);

  if (!(header == "reak_serialization::protobuf_archive")) {
    throw std::ios_base::failure(corrupt_header_msg);
  }
  if (version != 2) {
    throw std::ios_base::failure(unknown_version_msg);
  }
}

protobuf_iarchive::~protobuf_iarchive() = default;
;

iarchive& protobuf_iarchive::load_serializable_ptr(
    serializable_shared_pointer& item) {
  archive_object_header hdr;
  item = serializable_shared_pointer();

  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 2) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading serializable object pointer should "
          "have wire-type 2. Got chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }

  protobuf_iarchive::load_varint(hdr.size);
  std::streampos start_pos = file_stream->tellg();
  if (start_pos < 0) {
    throw std::ios_base::failure(bad_stream_msg);
  }

  std::vector<unsigned int> type_id_vect;
  unsigned int i = 0;
  do {
    *this >> i;
    type_id_vect.push_back(i);
  } while (i != 0);

  *this >> hdr.type_version >> hdr.object_ID >> hdr.is_external;

  if (hdr.object_ID == 0) {
    // item already null.
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
        throw std::ios_base::failure(bad_stream_msg);
      }
    }
    return *this;
  }
  if ((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    item = mObjRegistry[hdr.object_ID];
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
        throw std::ios_base::failure(bad_stream_msg);
      }
    }
    return *this;
  }

  if (hdr.is_external) {
    std::string ext_filename;
    *this >> ext_filename;
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
        throw std::ios_base::failure(bad_stream_msg);
      }
    }

    protobuf_iarchive a(
        ext_filename);  // if this throws, let it propagate up (no point catching and throwing).
    a >> item;

    return *this;
  }

  // Find the class in question in the repository.
  rtti::so_type* p =
      rtti::so_type_repo::get_instance().find_type(type_id_vect.data());
  if ((p == nullptr) || (p->version() < hdr.type_version)) {
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
        throw std::ios_base::failure(bad_stream_msg);
      }
    }
    throw unsupported_type(unsupported_type::not_found_in_repo,
                           type_id_vect.data());
  }
  std::shared_ptr<shared_object> po(p->create_object());
  if (!po) {
    std::streampos end_pos = file_stream->tellg();
    if (std::streampos(hdr.size) + start_pos != end_pos) {
      if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
        throw std::ios_base::failure(bad_stream_msg);
      }
    }
    throw unsupported_type(unsupported_type::could_not_create,
                           type_id_vect.data());
  }

  item = po;
  if (hdr.object_ID < mObjRegistry.size()) {
    mObjRegistry[hdr.object_ID] = item;
  } else if (hdr.object_ID == mObjRegistry.size()) {
    mObjRegistry.push_back(
        item);  // in theory, only this condition should occur
  } else if (hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = item;
  }

  item->load(*this, hdr.type_version);

  std::streampos end_pos = file_stream->tellg();
  if (std::streampos(hdr.size) + start_pos != end_pos) {
    if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
      throw std::ios_base::failure(bad_stream_msg);
    }
  }

  return *this;
}

iarchive& protobuf_iarchive::load_serializable_ptr(
    const std::pair<std::string, serializable_shared_pointer&>& item) {
  return protobuf_iarchive::load_serializable_ptr(item.second);
}

iarchive& protobuf_iarchive::load_serializable(serializable& item) {
  archive_object_header hdr;

  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 2) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading serializable object should have "
          "wire-type 2. Got chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }

  protobuf_iarchive::load_varint(hdr.size);
  std::streampos start_pos = file_stream->tellg();

  std::vector<unsigned int> type_id_vect;
  unsigned int i = 1;
  while (i != 0) {
    *this >> i;
    type_id_vect.push_back(i);
  }

  *this >> hdr.type_version;

  item.load(*this, hdr.type_version);
  std::streampos end_pos = file_stream->tellg();

  if (std::streampos(hdr.size) + start_pos != end_pos) {
    if (!file_stream->seekg(start_pos + std::streampos(hdr.size))) {
      throw std::ios_base::failure(bad_stream_msg);
    }
  }

  return *this;
}

iarchive& protobuf_iarchive::load_serializable(
    const std::pair<std::string, serializable&>& item) {
  return protobuf_iarchive::load_serializable(item.second);
}

iarchive& protobuf_iarchive::load_char(char& i) {
  std::int64_t il = 0;
  protobuf_iarchive::load_int(il);
  i = static_cast<char>(il);
  return *this;
}

iarchive& protobuf_iarchive::load_char(const std::pair<std::string, char&>& i) {
  return protobuf_iarchive::load_char(i.second);
}

iarchive& protobuf_iarchive::load_unsigned_char(unsigned char& u) {
  std::uint64_t ul = 0;
  protobuf_iarchive::load_unsigned_int(ul);
  u = static_cast<unsigned char>(ul);
  return *this;
}

iarchive& protobuf_iarchive::load_unsigned_char(
    const std::pair<std::string, unsigned char&>& u) {
  return protobuf_iarchive::load_unsigned_char(u.second);
}

iarchive& protobuf_iarchive::load_int(std::int64_t& i) {
  std::uint64_t u = 0;
  protobuf_iarchive::load_unsigned_int(u);
  i = (u >> 1) ^ (-static_cast<std::int64_t>(u & 1));
  return *this;
}

iarchive& protobuf_iarchive::load_int(
    const std::pair<std::string, std::int64_t&>& i) {
  return protobuf_iarchive::load_int(i.second);
}

void protobuf_iarchive::load_varint(std::uint64_t& u) {
  u = 0;
  std::uint8_t tmp = 0;
  if (!file_stream->read(reinterpret_cast<char*>(&tmp), 1)) {
    throw std::ios_base::failure(unexpected_eof_msg);
  }
  std::uint8_t shifts = 0;
  u = tmp & 0x7F;
  while ((tmp & 0x80) != 0) {
    if (!file_stream->read(reinterpret_cast<char*>(&tmp), 1)) {
      throw std::ios_base::failure(unexpected_eof_msg);
    }
    shifts += 7;
    u |= (tmp & 0x7F) << shifts;
  }
}

iarchive& protobuf_iarchive::load_unsigned_int(std::uint64_t& u) {
  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 0) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading varint should have wire-type 0. Got "
          "chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }

  protobuf_iarchive::load_varint(u);
  return *this;
}

iarchive& protobuf_iarchive::load_unsigned_int(
    const std::pair<std::string, std::uint64_t&>& u) {
  return protobuf_iarchive::load_unsigned_int(u.second);
}

iarchive& protobuf_iarchive::load_float(float& f) {
  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 5) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading float should have wire-type 5. Got "
          "chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }
  if (!file_stream->read(reinterpret_cast<char*>(&f), sizeof(float))) {
    throw std::ios_base::failure(unexpected_eof_msg);
  }
  from_endian<std::endian::little>(f);
  return *this;
}

iarchive& protobuf_iarchive::load_float(
    const std::pair<std::string, float&>& f) {
  return protobuf_iarchive::load_float(f.second);
}

iarchive& protobuf_iarchive::load_double(double& d) {
  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 1) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading double should have wire-type 1. Got "
          "chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }
  if (!file_stream->read(reinterpret_cast<char*>(&d), sizeof(double))) {
    throw std::ios_base::failure(unexpected_eof_msg);
  }
  from_endian<std::endian::little>(d);
  return *this;
}

iarchive& protobuf_iarchive::load_double(
    const std::pair<std::string, double&>& d) {
  return protobuf_iarchive::load_double(d.second);
}

iarchive& protobuf_iarchive::load_bool(bool& b) {
  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 0) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading bool should have wire-type 0. Got "
          "chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }
  char tmp = 0;
  if (!file_stream->read(&tmp, 1)) {
    throw std::ios_base::failure(unexpected_eof_msg);
  }
  b = (tmp != 0);
  return *this;
}

iarchive& protobuf_iarchive::load_bool(const std::pair<std::string, bool&>& b) {
  return protobuf_iarchive::load_bool(b.second);
}

iarchive& protobuf_iarchive::load_string(std::string& s) {
  std::uint64_t chunk_hdr = 0;
  protobuf_iarchive::load_varint(chunk_hdr);
  if ((chunk_hdr & 0x07) != 2) {
    std::streampos current_pos = file_stream->tellg();
    file_stream->seekg(std::ios_base::beg);
    std::streampos start_pos = file_stream->tellg();
    std::stringstream ss;
    ss << "Protobuf archive is inconsistent with requested read operation! "
          "Loading string should have wire-type 2. Got "
          "chunk-ID: "
       << std::hex << chunk_hdr << " at offset " << std::dec
       << (current_pos - start_pos) << ".";
    throw std::ios_base::failure(ss.str());
  }
  std::uint64_t u = 0;
  protobuf_iarchive::load_varint(u);
  s.resize(u);
  if (!file_stream->read(s.data(), u)) {
    throw std::ios_base::failure(unexpected_eof_msg);
  }
  return *this;
}

iarchive& protobuf_iarchive::load_string(
    const std::pair<std::string, std::string&>& s) {
  return protobuf_iarchive::load_string(s.second);
}

protobuf_oarchive::protobuf_oarchive(const std::string& file_name) {
  field_ids.push(0);
  repeat_state.push(0);

  file_stream = std::make_shared<std::ofstream>(
      file_name.c_str(), std::ios::binary | std::ios::out);

  std::string header = "reak_serialization::protobuf_archive";
  protobuf_oarchive::save_string(header);
  std::uint64_t version = 2;
  protobuf_oarchive::save_unsigned_int(version);
}

protobuf_oarchive::protobuf_oarchive(std::ostream& stream) {
  field_ids.push(0);
  repeat_state.push(0);

  file_stream = std::shared_ptr<std::ostream>(&stream, null_deleter());

  std::string header = "reak_serialization::protobuf_archive";
  protobuf_oarchive::save_string(header);
  std::uint64_t version = 2;
  protobuf_oarchive::save_unsigned_int(version);
}

protobuf_oarchive::~protobuf_oarchive() = default;

oarchive& protobuf_oarchive::save_to_new_archive_impl(
    const serializable_shared_pointer& item, const std::string& file_name) {
  std::uint64_t chunk_hdr =
      (field_ids.top() << 3) | 2;  // wire-type 2: length-delimited.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint(chunk_hdr);

  std::shared_ptr<std::ostream> tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  auto str_stream =
      std::make_shared<std::stringstream>(std::ios::binary | std::ios::out);
  file_stream = str_stream;
  field_ids.push(0);
  repeat_state.push(0);

  archive_object_header hdr;
  bool already_saved(false);
  const unsigned int* type_id = nullptr;

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
    type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = true;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    already_saved = true;
  }

  protobuf_oarchive::start_repeated_field("unsigned int");
  while (((type_id) != nullptr) && ((*type_id) != 0U)) {
    protobuf_oarchive::save_unsigned_int(*type_id);
    ++type_id;
  }
  protobuf_oarchive::save_unsigned_int(0);
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int(hdr.type_version);
  protobuf_oarchive::save_unsigned_int(hdr.object_ID);
  protobuf_oarchive::save_bool(hdr.is_external);

  if (!already_saved) {
    protobuf_oarchive::save_string(file_name);

    protobuf_oarchive a(file_name);
    a << item;
  }

  file_stream = tmp_str_ptr;
  field_ids.pop();
  repeat_state.pop();

  str_stream->seekp(0, std::ios::end);
  hdr.size = static_cast<std::uint64_t>(str_stream->tellp());
  str_stream->clear();

  protobuf_oarchive::save_varint(hdr.size);
  if (!(*file_stream << str_stream->str())) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }

  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }

  return *this;
}

oarchive& protobuf_oarchive::save_to_new_archive_named_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& item,
    const std::string& file_name) {
  return protobuf_oarchive::save_to_new_archive_impl(item.second, file_name);
}

oarchive& protobuf_oarchive::save_serializable_ptr(
    const serializable_shared_pointer& item) {
  std::uint64_t chunk_hdr =
      (field_ids.top() << 3) | 2;  // wire-type 2: length-delimited.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint(chunk_hdr);

  std::shared_ptr<std::ostream> tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  auto str_stream =
      std::make_shared<std::stringstream>(std::ios::binary | std::ios::out);
  file_stream = str_stream;
  field_ids.push(0);
  repeat_state.push(0);

  archive_object_header hdr;
  bool already_saved(false);
  const unsigned int* type_id = nullptr;

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
    type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = false;
  } else {
    hdr.type_version = 0;
    hdr.object_ID = 0;
    hdr.is_external = false;
    already_saved = true;
  }

  protobuf_oarchive::start_repeated_field("unsigned int");
  while (((type_id) != nullptr) && ((*type_id) != 0U)) {
    protobuf_oarchive::save_unsigned_int(*type_id);
    ++type_id;
  }
  protobuf_oarchive::save_unsigned_int(0);
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int(hdr.type_version);
  protobuf_oarchive::save_unsigned_int(hdr.object_ID);
  protobuf_oarchive::save_bool(hdr.is_external);

  if (!already_saved) {
    item->save(*this, hdr.type_version);
  }

  file_stream = tmp_str_ptr;
  field_ids.pop();
  repeat_state.pop();

  str_stream->seekp(0, std::ios::end);
  hdr.size = static_cast<std::uint64_t>(str_stream->tellp());
  str_stream->clear();

  protobuf_oarchive::save_varint(hdr.size);
  if (!(*file_stream << str_stream->str())) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }

  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }

  return *this;
}

oarchive& protobuf_oarchive::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& item) {
  return protobuf_oarchive::save_serializable_ptr(item.second);
}

oarchive& protobuf_oarchive::save_serializable(const serializable& item) {
  std::uint64_t chunk_hdr =
      (field_ids.top() << 3) | 2;  // wire-type 2: length-delimited.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint(chunk_hdr);

  std::shared_ptr<std::ostream> tmp_str_ptr = file_stream;

  // TODO: Fix the exception-safety of this buffer-swapping code (ensure a swap back on exception):
  auto str_stream =
      std::make_shared<std::stringstream>(std::ios::binary | std::ios::out);
  file_stream = str_stream;
  field_ids.push(0);
  repeat_state.push(0);

  archive_object_header hdr;
  const unsigned int* type_id = item.get_object_type()->id_begin();
  hdr.type_version = item.get_object_type()->version();
  hdr.object_ID = 0;
  hdr.is_external = false;
  hdr.size = 0;

  protobuf_oarchive::start_repeated_field("unsigned int");
  while (*type_id != 0U) {
    protobuf_oarchive::save_unsigned_int(*type_id);
    ++type_id;
  }
  protobuf_oarchive::save_unsigned_int(0);
  protobuf_oarchive::finish_repeated_field();

  protobuf_oarchive::save_unsigned_int(hdr.type_version);

  item.save(*this, hdr.type_version);

  file_stream = tmp_str_ptr;
  field_ids.pop();
  repeat_state.pop();

  str_stream->seekp(0, std::ios::end);
  hdr.size = static_cast<std::uint64_t>(str_stream->tellp());
  str_stream->clear();

  protobuf_oarchive::save_varint(hdr.size);
  if (!(*file_stream << str_stream->str())) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }

  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }

  return *this;
}

oarchive& protobuf_oarchive::save_serializable(
    const std::pair<std::string, const serializable&>& item) {
  return protobuf_oarchive::save_serializable(item.second);
}

oarchive& protobuf_oarchive::save_char(char i) {
  return protobuf_oarchive::save_int(i);
}

oarchive& protobuf_oarchive::save_char(const std::pair<std::string, char>& i) {
  return protobuf_oarchive::save_char(i.second);
}

oarchive& protobuf_oarchive::save_unsigned_char(unsigned char u) {
  return protobuf_oarchive::save_unsigned_int(u);
}

oarchive& protobuf_oarchive::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  return protobuf_oarchive::save_unsigned_char(u.second);
}

oarchive& protobuf_oarchive::save_int(std::int64_t i) {
  protobuf_oarchive::save_unsigned_int((i << 1) ^
                                       (i >> (sizeof(std::int64_t) * 8 - 1)));
  return *this;
}

oarchive& protobuf_oarchive::save_int(
    const std::pair<std::string, std::int64_t>& i) {
  return protobuf_oarchive::save_int(i.second);
}

void protobuf_oarchive::save_varint(std::uint64_t u) {
  std::array<std::uint8_t, 10> buf = {
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0};  // 80-bits, supports at most a 64-bit varint.
  std::uint8_t* pbuf = buf.data();
  *pbuf = (u & 0x7F);
  u >>= 7;
  while (u != 0U) {
    *pbuf |= 0x80;  // set first msb because there is more to come.
    pbuf++;
    *pbuf = (u & 0x7F);
    u >>= 7;
  }
  if (!file_stream->write(reinterpret_cast<char*>(buf.data()),
                          pbuf - buf.data() + 1)) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }
}

oarchive& protobuf_oarchive::save_unsigned_int(std::uint64_t u) {
  std::uint64_t chunk_hdr = (field_ids.top() << 3);  // wire-type 0: varint.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  };
  protobuf_oarchive::save_varint(chunk_hdr);
  protobuf_oarchive::save_varint(u);
  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }
  return *this;
}

oarchive& protobuf_oarchive::save_unsigned_int(
    const std::pair<std::string, std::uint64_t>& u) {
  return protobuf_oarchive::save_unsigned_int(u.second);
}

oarchive& protobuf_oarchive::save_float(float f) {
  std::uint64_t chunk_hdr = (field_ids.top() << 3) | 5;  // wire-type 5: 32-bit.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  }
  protobuf_oarchive::save_varint(chunk_hdr);
  to_endian<std::endian::little>(f);
  if (!file_stream->write(reinterpret_cast<char*>(&f), sizeof(float))) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }
  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }
  return *this;
}

oarchive& protobuf_oarchive::save_float(
    const std::pair<std::string, float>& f) {
  return protobuf_oarchive::save_float(f.second);
}

oarchive& protobuf_oarchive::save_double(double d) {
  std::uint64_t chunk_hdr = (field_ids.top() << 3) | 1;  // wire-type 1: 64-bit.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  }
  protobuf_oarchive::save_varint(chunk_hdr);
  to_endian<std::endian::little>(d);
  if (!file_stream->write(reinterpret_cast<char*>(&d), sizeof(double))) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }
  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }
  return *this;
}

oarchive& protobuf_oarchive::save_double(
    const std::pair<std::string, double>& d) {
  return protobuf_oarchive::save_double(d.second);
}

oarchive& protobuf_oarchive::save_bool(bool b) {
  std::uint64_t chunk_hdr = (field_ids.top() << 3);  // wire-type 0: varint.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  }
  protobuf_oarchive::save_varint(chunk_hdr);
  char tmp = 0;
  if (b) {
    tmp = 1;
  }
  if (!file_stream->put(tmp)) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }
  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }
  return *this;
}

oarchive& protobuf_oarchive::save_bool(const std::pair<std::string, bool>& b) {
  return protobuf_oarchive::save_bool(b.second);
}

oarchive& protobuf_oarchive::save_string(const std::string& s) {
  std::uint64_t chunk_hdr =
      (field_ids.top() << 3) | 2;  // wire-type 2: length-delimited.
  if ((repeat_state.top() & 0x02) != 0U) {
    if ((repeat_state.top() & 0x04) != 0U) {
      chunk_hdr += 8;
    }
    repeat_state.top() ^= 0x04;
  }
  protobuf_oarchive::save_varint(chunk_hdr);
  std::uint64_t u = s.length();
  protobuf_oarchive::save_varint(u);
  if (!file_stream->write(s.data(), s.length())) {
    throw std::ios_base::failure(bad_out_stream_msg);
  }
  if (repeat_state.top() == 0U) {
    ++(field_ids.top());  // increment the field ID.
  }
  return *this;
};

oarchive& protobuf_oarchive::save_string(
    const std::pair<std::string, const std::string&>& s) {
  return protobuf_oarchive::save_string(s.second);
}

void protobuf_oarchive::start_repeated_field(const std::string& /*type_name*/) {
  repeat_state.push(1);
}

void protobuf_oarchive::start_repeated_field(const std::string& /*type_name*/,
                                             const std::string& /*name*/) {
  repeat_state.push(1);
}

void protobuf_oarchive::finish_repeated_field() {
  repeat_state.pop();
  ++(field_ids.top());
}

void protobuf_oarchive::start_repeated_pair(
    const std::string& /*type_name_1*/, const std::string& /*type_name_2*/) {
  repeat_state.push(3);
}

void protobuf_oarchive::start_repeated_pair(const std::string& /*type_name_1*/,
                                            const std::string& /*type_name_2*/,
                                            const std::string& /*name*/) {
  repeat_state.push(3);
}

void protobuf_oarchive::finish_repeated_pair() {
  repeat_state.pop();
  field_ids.top() += 2;
}

protobuf_schemer::protobuf_schemer() {
  field_ids.push(0);
  repeat_state.push(0);

  file_stream = std::make_shared<std::stringstream>();
}

protobuf_schemer::~protobuf_schemer() = default;
;

void protobuf_schemer::print_schemes(std::ostream& stream) {
  for (auto& scheme : schemes) {
    stream << scheme << "\n\n";
  }
}

std::uint64_t protobuf_schemer::get_chunk_hdr() {
  std::uint64_t chunk_hdr = 0;
  if ((repeat_state.top() & 0x08) != 0U) {
    return ~chunk_hdr;
  }
  chunk_hdr = field_ids.top();

  return chunk_hdr;
}

oarchive& protobuf_schemer::save_to_new_archive_impl(
    const serializable_shared_pointer& item, const std::string& /*file_name*/) {
  return protobuf_schemer::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("item", item));
}

oarchive& protobuf_schemer::save_to_new_archive_named_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& item,
    const std::string& /*file_name*/) {
  return protobuf_schemer::save_serializable_ptr(item);
}

oarchive& protobuf_schemer::save_serializable_ptr(
    const serializable_shared_pointer& item) {
  return protobuf_schemer::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("item_ptr",
                                                                 item));
}

oarchive& protobuf_schemer::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& item) {
  if (!item.second) {
    return *this;
  }

  std::uint64_t chunk_hdr = get_chunk_hdr();

  constexpr auto tname =
      rtti::get_type_id<serializable_shared_pointer>::type_name;
  std::string obj_type_name = std::string(tname);
  obj_type_name += "<" + item.second->get_object_type()->name() + ">";
  auto it = mObjRegMap.find(item.second);

  if (it == mObjRegMap.end()) {
    std::shared_ptr<std::ostream> tmp_str_ptr = file_stream;

    auto str_stream =
        std::make_shared<std::stringstream>(std::ios::binary | std::ios::out);
    file_stream = str_stream;
    field_ids.push(4);
    repeat_state.push(0);

    *file_stream << "message " << obj_type_name << " {\n"
                 << "  repeated uint32 type_ID = 0;\n"
                 << "  required uint32 version = 1;\n"
                 << "  required luid32 object_ID = 2;\n"
                 << "  required bool is_external = 3;\n";

    mObjRegistry.push_back(item.second);
    mObjRegMap[item.second] = mObjRegistry.size() - 1;

    item.second->save(*this, item.second->get_object_type()->version());
    *file_stream << "}\n";

    file_stream = tmp_str_ptr;
    field_ids.pop();
    repeat_state.pop();

    auto itm = scheme_map.find(obj_type_name);
    if (itm == scheme_map.end()) {
      schemes.push_back(str_stream->str());
      scheme_map[obj_type_name] = schemes.size() - 1;
    } else {
      std::string s_tmp = str_stream->str();
      if (schemes[itm->second].length() < s_tmp.length()) {
        schemes[itm->second] = s_tmp;
      }
    }
  }

  if (~chunk_hdr == 0) {
    return *this;
  }

  *file_stream << "  required " << obj_type_name << " " << item.first << " = "
               << chunk_hdr << ";\n";

  field_ids.top() += 1;

  return *this;
}

oarchive& protobuf_schemer::save_serializable(const serializable& item) {
  return protobuf_schemer::save_serializable(
      std::pair<std::string, const serializable&>("item", item));
}

oarchive& protobuf_schemer::save_serializable(
    const std::pair<std::string, const serializable&>& item) {
  std::uint64_t chunk_hdr = get_chunk_hdr();

  std::shared_ptr<std::ostream> tmp_str_ptr = file_stream;

  auto str_stream = std::make_shared<std::stringstream>();
  file_stream = str_stream;
  field_ids.push(2);
  repeat_state.push(0);

  *file_stream << "message " << item.second.get_object_type()->name() << " {\n"
               << "  repeated uint32 type_ID = 0;\n"
               << "  required uint32 version = 1;\n";

  item.second.save(*this, item.second.get_object_type()->version());

  *file_stream << "}\n";

  file_stream = tmp_str_ptr;
  field_ids.pop();
  repeat_state.pop();

  auto itm = scheme_map.find(item.second.get_object_type()->name());
  if (itm == scheme_map.end()) {
    schemes.push_back(str_stream->str());
    scheme_map[item.second.get_object_type()->name()] = schemes.size() - 1;
  } else {
    std::string s_tmp = str_stream->str();
    if (schemes[itm->second].length() < s_tmp.length()) {
      schemes[itm->second] = s_tmp;
    }
  }

  if (~chunk_hdr == 0) {
    return *this;
  }

  *file_stream << "  required " << item.second.get_object_type()->name() << " "
               << item.first << " = " << chunk_hdr << ";\n";

  field_ids.top() += 1;

  return *this;
}

oarchive& protobuf_schemer::save_char(char i) {
  return protobuf_schemer::save_char(std::pair<std::string, char>("i", i));
}

oarchive& protobuf_schemer::save_char(const std::pair<std::string, char>& i) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required int32 " << i.first << " = " << chunk_hdr << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_unsigned_char(unsigned char u) {
  return protobuf_schemer::save_unsigned_char(
      std::pair<std::string, unsigned char>("u", u));
}

oarchive& protobuf_schemer::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required uint32 " << u.first << " = " << chunk_hdr
               << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_int(std::int64_t i) {
  return protobuf_schemer::save_int(
      std::pair<std::string, std::int64_t>("i", i));
}

oarchive& protobuf_schemer::save_int(
    const std::pair<std::string, std::int64_t>& i) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required int32 " << i.first << " = " << chunk_hdr << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_unsigned_int(std::uint64_t u) {
  return protobuf_schemer::save_unsigned_int(
      std::pair<std::string, std::uint64_t>("u", u));
}

oarchive& protobuf_schemer::save_unsigned_int(
    const std::pair<std::string, std::uint64_t>& u) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required uint32 " << u.first << " = " << chunk_hdr
               << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_float(float f) {
  return protobuf_schemer::save_float(std::pair<std::string, float>("f", f));
}

oarchive& protobuf_schemer::save_float(const std::pair<std::string, float>& f) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required float " << f.first << " = " << chunk_hdr << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_double(double d) {
  return protobuf_schemer::save_double(std::pair<std::string, double>("d", d));
}

oarchive& protobuf_schemer::save_double(
    const std::pair<std::string, double>& d) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required double " << d.first << " = " << chunk_hdr
               << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_bool(bool b) {
  return protobuf_schemer::save_bool(std::pair<std::string, bool>("b", b));
}

oarchive& protobuf_schemer::save_bool(const std::pair<std::string, bool>& b) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required bool " << b.first << " = " << chunk_hdr << ";\n";
  field_ids.top() += 1;
  return *this;
}

oarchive& protobuf_schemer::save_string(const std::string& s) {
  return protobuf_schemer::save_string(
      std::pair<std::string, const std::string&>("str", s));
}

oarchive& protobuf_schemer::save_string(
    const std::pair<std::string, const std::string&>& s) {
  std::uint64_t chunk_hdr = get_chunk_hdr();
  if (~chunk_hdr == 0) {
    return *this;
  }
  *file_stream << "  required string " << s.first << " = " << chunk_hdr
               << ";\n";
  field_ids.top() += 1;
  return *this;
}

void protobuf_schemer::start_repeated_field(const std::string& type_name) {
  repeat_state.push(9);
  *file_stream << "  repeated " << type_name << " value = " << field_ids.top()
               << ";\n";
}

void protobuf_schemer::start_repeated_field(const std::string& type_name,
                                            const std::string& name) {
  repeat_state.push(9);
  *file_stream << "  repeated " << type_name << " " << name << " = "
               << field_ids.top() << ";\n";
}

void protobuf_schemer::finish_repeated_field() {
  repeat_state.pop();
  field_ids.top() += 1;
}

void protobuf_schemer::start_repeated_pair(const std::string& type_name_1,
                                           const std::string& type_name_2) {
  repeat_state.push(11);
  *file_stream << "  repeated " << type_name_1
               << " map_key = " << field_ids.top() << ";\n"
               << "  repeated " << type_name_2
               << " map_value = " << (field_ids.top() + 1) << ";\n";
}

void protobuf_schemer::start_repeated_pair(const std::string& type_name_1,
                                           const std::string& type_name_2,
                                           const std::string& name) {
  repeat_state.push(11);
  *file_stream << "  repeated " << type_name_1 << " " << name
               << "_key = " << field_ids.top() << ";\n"
               << "  repeated " << type_name_2 << " " << name
               << "_value = " << (field_ids.top() + 1) << ";\n";
}

void protobuf_schemer::finish_repeated_pair() {
  repeat_state.pop();
  field_ids.top() += 2;
}

}  // namespace ReaK::serialization
