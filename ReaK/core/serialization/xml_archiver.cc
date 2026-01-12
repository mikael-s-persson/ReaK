
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

#include "ReaK/core/serialization/xml_archiver.h"

#include "ReaK/core/base/shared_object.h"
#include "ReaK/core/rtti/rtti.h"

#include "ReaK/core/serialization/archiving_exceptions.h"

#include <fstream>
#include <map>
#include <vector>

#include <cstdlib>

namespace ReaK::serialization {

char xml_iarchive::get_next_char() {
  char c = 0;
  file_stream->get(c);
  while ((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r')) {
    file_stream->get(c);
  }
  return c;
}

std::string xml_iarchive::read_token() {
  std::string result;
  char c = get_next_char();
  if (c != '<') {
    return result;
  }
  c = get_next_char();
  if ((c == '!') || (c == '?')) {
    std::array<char, 512> line_str{};
    file_stream->getline(line_str.data(), line_str.size());
    return read_token();
  }
  while (c != '>') {
    result += c;
    file_stream->get(c);
  }
  return result;
}

void xml_iarchive::skip_to_end_token(const std::string& name) {
  std::string token = read_token();
  trim_str(token);
  while (token != "/" + name) {
    token = read_token();
    trim_str(token);
  }
}

void xml_iarchive::trim_str(std::string& s) {
  unsigned int i = 0;
  for (; ((i < s.size()) && ((s[i] == ' ') || (s[i] == '\t') ||
                             (s[i] == '\n') || (s[i] == '\r')));
       ++i) {}
  std::string result;
  for (; ((i < s.size()) && (s[i] != ' ') && (s[i] != '\t') && (s[i] != '\n') &&
          (s[i] != '\r'));
       ++i) {
    result += s[i];
  }
  s = result;
}

bool xml_iarchive::read_named_value(const std::string& value_name,
                                    std::string& value_str) {
  std::string token = read_token();
  trim_str(token);
  if ((value_name.empty()) || (token != value_name)) {
    return false;
  }

  char c = 0;
  file_stream->get(c);
  while (c != '\"') {
    file_stream->get(c);
  }

  value_str.clear();
  file_stream->get(c);
  while (c != '\"') {
    value_str += c;
    file_stream->get(c);
  }

  token = read_token();
  unsigned int i = 0;
  for (; ((i < token.size()) && (token[i] != '/')); ++i) {}
  std::string tmp;
  ++i;
  for (; ((i < token.size()) && (token[i] != value_name[0])); ++i) {}
  for (; ((i < token.size()) && (tmp.size() < value_name.size())); ++i) {
    tmp += token[i];
  }
  return tmp == value_name;
}

archive_object_header xml_iarchive::read_header(
    const std::string& obj_name, std::vector<std::uint32_t>& out_type_id) {
  archive_object_header result;
  out_type_id.clear();

  std::string token = read_token();
  if (token.empty()) {
    return result;
  }

  std::string name;
  unsigned int i = 0;
  for (; ((i < token.size()) && (token[i] == ' ')); ++i) {}
  for (; ((i < token.size()) && (token[i] != ' ')); ++i) {
    name += token[i];
  }

  if ((name != obj_name) || (i == token.size())) {
    return result;
  }

  std::map<std::string, std::string> values;
  while (i < token.size()) {
    for (; ((i < token.size()) && ((token[i] == ' ') || (token[i] == '\t') ||
                                   (token[i] == '\n') || (token[i] == '\r')));
         ++i) {}
    std::string value_key;
    for (; ((i < token.size()) && (token[i] != ' ') && (token[i] != '='));
         ++i) {
      value_key += token[i];
    }
    std::string value_str;
    for (; ((i < token.size()) && (token[i] != '\"')); ++i) {}
    ++i;
    for (; ((i < token.size()) && (token[i] != '\"')); ++i) {
      value_str += token[i];
    }
    ++i;
    values[value_key] = value_str;
  }

  std::string id_str = values["type_ID"];
  if (!id_str.empty()) {
    for (i = 0; i < id_str.size(); ++i) {
      std::string numstr;
      for (; ((i < id_str.size()) && (id_str[i] != '.')); ++i) {
        numstr += id_str[i];
      }
      out_type_id.push_back(strtoul(numstr.c_str(), nullptr, 0));
    }
  }

  if (values["version"].empty()) {
    result.type_version = 0;
  } else {
    result.type_version = strtoul(values["version"].c_str(), nullptr, 0);
  }

  if (values["object_ID"].empty()) {
    result.object_ID = 0;
  } else {
    result.object_ID = strtoul(values["object_ID"].c_str(), nullptr, 0);
  }

  if (values["is_external"].empty()) {
    result.is_external = false;
  } else if (values["is_external"] == "true") {
    result.is_external = true;
  } else {
    result.is_external = false;
  }

  return result;
}

xml_iarchive::xml_iarchive(const std::string& file_name) {

  file_stream =
      std::make_shared<std::ifstream>(file_name.c_str(), std::ios::in);

  std::vector<std::uint32_t> type_id;
  archive_object_header global_hdr = read_header("reak_serialization", type_id);
  if (global_hdr.type_version != 2) {
    throw std::ios_base::failure("ReaK XML Archive is of an unknown version!");
  }
}

xml_iarchive::xml_iarchive(std::istream& stream) {

  file_stream = std::shared_ptr<std::istream>(&stream, null_deleter());

  std::vector<std::uint32_t> type_id;
  archive_object_header global_hdr = read_header("reak_serialization", type_id);
  if (global_hdr.type_version != 2) {
    throw std::ios_base::failure("ReaK XML Archive is of an unknown version!");
  }
}

xml_iarchive::~xml_iarchive() = default;

iarchive& xml_iarchive::load_serializable_ptr(
    serializable_shared_pointer& item) {
  return xml_iarchive::load_serializable_ptr(
      std::pair<std::string, serializable_shared_pointer&>("Item", item));
}

iarchive& xml_iarchive::load_serializable_ptr(
    const std::pair<std::string, serializable_shared_pointer&>& item) {
  archive_object_header hdr;
  item.second = serializable_shared_pointer();

  std::vector<std::uint32_t> type_id;
  hdr = read_header(item.first, type_id);
  if ((type_id.empty()) || (hdr.type_version == 0) || (hdr.object_ID == 0)) {
    skip_to_end_token(item.first);
    return *this;
  }

  if ((hdr.object_ID < mObjRegistry.size()) && (mObjRegistry[hdr.object_ID])) {
    item.second = mObjRegistry[hdr.object_ID];
    skip_to_end_token(item.first);
    return *this;
  }

  if (hdr.is_external) {
    std::string ext_filename;
    xml_iarchive::load_string(
        std::pair<std::string, std::string&>("filename", ext_filename));

    skip_to_end_token(item.first);

    xml_iarchive a(
        ext_filename);  // if this throws, let it propagate up (no point catching and throwing).
    a & item;

    return *this;
  }

  // Find the class in question in the repository.
  rtti::so_type* p =
      rtti::so_type_repo::get_instance().find_type(type_id.data());
  if ((p == nullptr) || (p->version() < hdr.type_version)) {
    skip_to_end_token(item.first);
    throw unsupported_type(unsupported_type::not_found_in_repo, type_id.data());
  }
  std::shared_ptr<shared_object> po(p->create_object());
  if (!po) {
    skip_to_end_token(item.first);
    throw unsupported_type(unsupported_type::could_not_create, type_id.data());
  }

  item.second = po;
  if (hdr.object_ID < mObjRegistry.size()) {
    mObjRegistry[hdr.object_ID] = item.second;
  } else if (hdr.object_ID == mObjRegistry.size()) {
    mObjRegistry.push_back(
        item.second);  // in theory, only this condition should occur
  } else if (hdr.object_ID > mObjRegistry.size()) {
    mObjRegistry.resize(hdr.object_ID + 1);
    mObjRegistry[hdr.object_ID] = item.second;
  }

  item.second->load(*this, hdr.type_version);

  skip_to_end_token(item.first);
  return *this;
}

iarchive& xml_iarchive::load_serializable(serializable& item) {
  return xml_iarchive::load_serializable(
      std::pair<std::string, serializable&>("Item", item));
}

iarchive& xml_iarchive::load_serializable(
    const std::pair<std::string, serializable&>& item) {
  archive_object_header hdr;

  std::vector<std::uint32_t> type_id;
  hdr = read_header(item.first, type_id);
  if ((type_id.empty()) || (hdr.type_version == 0)) {
    skip_to_end_token(item.first);
    return *this;
  }

  item.second.load(*this, hdr.type_version);

  skip_to_end_token(item.first);
  return *this;
}

iarchive& xml_iarchive::load_char(char& i) {
  return xml_iarchive::load_char(std::pair<std::string, char&>("char", i));
}

iarchive& xml_iarchive::load_char(const std::pair<std::string, char&>& i) {
  std::string value_str;
  if (read_named_value(i.first, value_str)) {
    if (value_str.empty()) {
      i.second = 0;
    } else {
      int temp = 0;
      std::stringstream(value_str) >> temp;
      i.second = static_cast<char>(temp);
    }
  } else {
    i.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_unsigned_char(unsigned char& u) {
  return xml_iarchive::load_unsigned_char(
      std::pair<std::string, unsigned char&>("unsigned_char", u));
}

iarchive& xml_iarchive::load_unsigned_char(
    const std::pair<std::string, unsigned char&>& u) {
  std::string value_str;
  if (read_named_value(u.first, value_str)) {
    if (value_str.empty()) {
      u.second = 0;
    } else {
      unsigned int temp = 0;
      std::stringstream(value_str) >> temp;
      u.second = static_cast<char>(temp);
    }
  } else {
    u.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_int(std::int64_t& i) {
  return xml_iarchive::load_int(
      std::pair<std::string, std::int64_t&>("int", i));
}

iarchive& xml_iarchive::load_int(
    const std::pair<std::string, std::int64_t&>& i) {
  std::string value_str;
  if (read_named_value(i.first, value_str)) {
    if (value_str.empty()) {
      i.second = 0;
    } else {
      std::stringstream(value_str) >> i.second;
    }
  } else {
    i.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_unsigned_int(std::uint64_t& u) {
  return xml_iarchive::load_unsigned_int(
      std::pair<std::string, std::uint64_t&>("unsigned_int", u));
}

iarchive& xml_iarchive::load_unsigned_int(
    const std::pair<std::string, std::uint64_t&>& u) {
  std::string value_str;
  if (read_named_value(u.first, value_str)) {
    if (value_str.empty()) {
      u.second = 0;
    } else {
      std::stringstream(value_str) >> u.second;
    }
  } else {
    u.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_float(float& f) {
  return xml_iarchive::load_float(std::pair<std::string, float&>("real", f));
}

iarchive& xml_iarchive::load_float(const std::pair<std::string, float&>& f) {
  std::string value_str;
  if (read_named_value(f.first, value_str)) {
    if (value_str.empty()) {
      f.second = 0;
    } else {
      std::stringstream(value_str) >> f.second;
    }
  } else {
    f.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_double(double& d) {
  return xml_iarchive::load_double(std::pair<std::string, double&>("real", d));
}

iarchive& xml_iarchive::load_double(const std::pair<std::string, double&>& d) {
  std::string value_str;
  if (read_named_value(d.first, value_str)) {
    if (value_str.empty()) {
      d.second = 0;
    } else {
      std::stringstream(value_str) >> d.second;
    }
  } else {
    d.second = 0;
  }
  return *this;
}

iarchive& xml_iarchive::load_bool(bool& b) {
  return xml_iarchive::load_bool(std::pair<std::string, bool&>("bool", b));
}

iarchive& xml_iarchive::load_bool(const std::pair<std::string, bool&>& b) {
  std::string value_str;
  if (read_named_value(b.first, value_str)) {
    if (value_str.empty()) {
      b.second = false;
    } else if (value_str == "true") {
      b.second = true;
    } else {
      b.second = false;
    }
  } else {
    b.second = false;
  }
  return *this;
}

iarchive& xml_iarchive::load_string(std::string& s) {
  return xml_iarchive::load_string(
      std::pair<std::string, std::string&>("string", s));
}

iarchive& xml_iarchive::load_string(
    const std::pair<std::string, std::string&>& s) {
  read_named_value(s.first, s.second);
  return *this;
}

xml_oarchive::xml_oarchive(const std::string& file_name) : tabulation(0) {

  file_stream =
      std::make_shared<std::ofstream>(file_name.c_str(), std::ios::out);

  (*file_stream)
      << R"(<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>)";
  (*file_stream) << "\n<!DOCTYPE reak_serialization>\n";
  (*file_stream) << "<reak_serialization version=\"2\">\n";
}

xml_oarchive::xml_oarchive(std::ostream& stream) : tabulation(0) {

  file_stream = std::shared_ptr<std::ostream>(&stream, null_deleter());

  (*file_stream)
      << R"(<?xml version="1.0" encoding="UTF-8" standalone="yes" ?>)";
  (*file_stream) << "\n<!DOCTYPE reak_serialization>\n";
  (*file_stream) << "<reak_serialization version=\"2\">\n";
}

xml_oarchive::~xml_oarchive() {
  (*file_stream) << "</reak_serialization>\n";
}

void xml_oarchive::add_tabulations() {
  for (unsigned int i = 0; i < tabulation; ++i) {
    (*file_stream) << "    ";
  }
}

oarchive& xml_oarchive::save_to_new_archive_impl(
    const serializable_shared_pointer& item, const std::string& file_name) {
  return xml_oarchive::save_to_new_archive_named_impl(
      std::pair<std::string, const serializable_shared_pointer&>("Item", item),
      file_name);
}

oarchive& xml_oarchive::save_to_new_archive_named_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& item,
    const std::string& file_name) {
  archive_object_header hdr;
  bool already_saved(false);

  if (item.second) {
    auto it = mObjRegMap.find(item.second);

    if (it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(item.second);
      mObjRegMap[item.second] = hdr.object_ID;
    }

    rtti::so_type* obj_type = item.second->get_object_type();
    const std::uint32_t* type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = true;
    hdr.size = 0;

    add_tabulations();
    (*file_stream) << "<" << item.first << " type_ID=\"";
    while (*type_id != 0U) {
      (*file_stream) << *type_id << ".";
      ++type_id;
    }
    (*file_stream) << "0\" version=\"" << hdr.type_version << "\" object_ID=\""
                   << hdr.object_ID << R"(" is_external="true">)"
                   << "\n";
  } else {
    already_saved = true;

    add_tabulations();
    (*file_stream)
        << "<" << item.first
        << R"( type_ID="0" version="0" object_ID="0" is_external="false">)"
        << "\n";
  }

  if (!already_saved) {
    tabulation++;
    add_tabulations();
    (*file_stream) << "<filename>\"" << file_name << "\"</filename>\n";
    tabulation--;

    xml_oarchive a(file_name);
    a & item;
  }

  add_tabulations();
  (*file_stream) << "</" << item.first << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_serializable_ptr(
    const serializable_shared_pointer& item) {
  return *this & std::pair<std::string, const serializable_shared_pointer&>(
                     "Item", item);
}

oarchive& xml_oarchive::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& item) {

  archive_object_header hdr;
  bool already_saved(false);

  if (item.second) {
    auto it = mObjRegMap.find(item.second);

    if (it != mObjRegMap.end()) {
      hdr.object_ID = it->second;
      already_saved = true;
    } else {
      hdr.object_ID = mObjRegistry.size();
      mObjRegistry.push_back(item.second);
      mObjRegMap[item.second] = hdr.object_ID;
    }

    rtti::so_type* obj_type = item.second->get_object_type();
    const std::uint32_t* type_id = obj_type->id_begin();
    hdr.type_version = obj_type->version();
    hdr.is_external = false;
    hdr.size = 0;

    add_tabulations();
    (*file_stream) << "<" << item.first << " type_ID=\"";
    while (*type_id != 0U) {
      (*file_stream) << *type_id << ".";
      ++type_id;
    }
    (*file_stream) << "0\" version=\"" << hdr.type_version << "\" object_ID=\""
                   << hdr.object_ID << R"(" is_external="false">)"
                   << "\n";
  } else {
    already_saved = true;

    add_tabulations();
    (*file_stream)
        << "<" << item.first
        << R"( type_ID="0" version="0" object_ID="0" is_external="false">)"
        << "\n";
  }

  if (!already_saved) {

    tabulation++;
    item.second->save(*this, hdr.type_version);
    tabulation--;
  }

  add_tabulations();
  (*file_stream) << "</" << item.first << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_serializable(const serializable& item) {
  return *this & std::pair<std::string, const serializable&>("Item", item);
}

oarchive& xml_oarchive::save_serializable(
    const std::pair<std::string, const serializable&>& item) {
  archive_object_header hdr;
  const std::uint32_t* type_id = item.second.get_object_type()->id_begin();
  hdr.type_version = item.second.get_object_type()->version();
  hdr.object_ID = 0;
  hdr.size = 0;
  hdr.is_external = false;

  add_tabulations();
  (*file_stream) << "<" << item.first << " type_ID=\"";
  while (*type_id != 0U) {
    (*file_stream) << *type_id << ".";
    ++type_id;
  }
  (*file_stream) << "0\" version=\"" << hdr.type_version << "\">\n";

  tabulation++;
  item.second.save(*this, hdr.type_version);
  tabulation--;

  add_tabulations();
  (*file_stream) << "</" << item.first << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_char(char i) {
  return xml_oarchive::save_char(std::pair<std::string, char>("char", i));
}

oarchive& xml_oarchive::save_char(const std::pair<std::string, char>& i) {
  add_tabulations();
  (*file_stream) << "<" << i.first << ">\"" << static_cast<int>(i.second)
                 << "\"</" << i.first << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_unsigned_char(unsigned char u) {
  return xml_oarchive::save_unsigned_char(
      std::pair<std::string, unsigned char>("unsigned_char", u));
}

oarchive& xml_oarchive::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  add_tabulations();
  (*file_stream) << "<" << u.first << ">\""
                 << static_cast<unsigned int>(u.second) << "\"</" << u.first
                 << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_int(std::int64_t i) {
  return xml_oarchive::save_int(std::pair<std::string, std::int64_t>("int", i));
}

oarchive& xml_oarchive::save_int(
    const std::pair<std::string, std::int64_t>& i) {
  add_tabulations();
  (*file_stream) << "<" << i.first << ">\"" << i.second << "\"</" << i.first
                 << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_unsigned_int(std::uint64_t u) {
  return xml_oarchive::save_unsigned_int(
      std::pair<std::string, std::uint64_t>("unsigned_int", u));
}

oarchive& xml_oarchive::save_unsigned_int(
    const std::pair<std::string, std::uint64_t>& u) {
  add_tabulations();
  (*file_stream) << "<" << u.first << ">\"" << u.second << "\"</" << u.first
                 << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_float(float f) {
  return xml_oarchive::save_float(std::pair<std::string, float>("real", f));
}

oarchive& xml_oarchive::save_float(const std::pair<std::string, float>& f) {
  add_tabulations();
  (*file_stream) << "<" << f.first << ">\"" << f.second << "\"</" << f.first
                 << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_double(double d) {
  return xml_oarchive::save_double(std::pair<std::string, double>("real", d));
}

oarchive& xml_oarchive::save_double(const std::pair<std::string, double>& d) {
  add_tabulations();
  (*file_stream) << "<" << d.first << ">\"" << d.second << "\"</" << d.first
                 << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_bool(bool b) {
  return xml_oarchive::save_bool(std::pair<std::string, bool>("bool", b));
}

oarchive& xml_oarchive::save_bool(const std::pair<std::string, bool>& b) {
  add_tabulations();
  (*file_stream) << "<" << b.first << ">\"" << (b.second ? "true" : "false")
                 << "\"</" << b.first << ">\n";
  return *this;
}

oarchive& xml_oarchive::save_string(const std::string& s) {
  return xml_oarchive::save_string(
      std::pair<std::string, const std::string&>("string", s));
}

oarchive& xml_oarchive::save_string(
    const std::pair<std::string, const std::string&>& s) {
  add_tabulations();
  (*file_stream) << "<" << s.first << ">\"" << s.second << "\"</" << s.first
                 << ">\n";
  return *this;
}

}  // namespace ReaK::serialization
