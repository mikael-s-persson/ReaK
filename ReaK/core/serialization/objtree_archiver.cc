
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

#include "ReaK/core/serialization/objtree_archiver.h"

#include "ReaK/core/base/named_object.h"
#include "ReaK/core/rtti/rtti.h"

#include "ReaK/core/serialization/scheme_builder.h"
#include "ReaK/core/serialization/type_schemes.h"

#include <map>
#include <sstream>
#include <string>

#include <algorithm>
#include <cctype>
#include <utility>

namespace ReaK::serialization {

std::string::iterator xml_field_editor::mark_field(
    std::string::iterator it_prev, std::string::iterator it_end,
    const std::string& fld_name, const std::shared_ptr<type_scheme>& scheme) {
  if (!scheme) {
    return it_prev;
  }

  std::string& xml_src = (*(p_parent->get_object_graph()))[node].xml_src;

  if (scheme->is_single_field()) {  // must be a primitive scheme.
    std::string test_seq = "<" + fld_name;
    std::string::iterator it =
        std::search(it_prev, it_end, test_seq.begin(), test_seq.end());
    if (it == it_end) {
      return it_prev;
    }
    src_markers.push_back(it - xml_src.begin());
    field_schemes.push_back(scheme);
    field_names.push_back(fld_name);
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 =
        std::search(it, it_end, test_seq.begin(), test_seq.end());
    if (it_end2 != it_end) {
      it_end2 += test_seq.length();
    }
    it_prev = it_end2;
  } else if (scheme->get_object_type() ==
             serializable_obj_scheme::get_static_object_type()) {
    std::string test_seq = "<" + fld_name;
    std::string::iterator it =
        std::search(it_prev, it_end, test_seq.begin(), test_seq.end());
    if (it == it_end) {
      return it_prev;
    }
    it = std::find(it, it_end, '>');
    while ((it != it_end) && (std::isspace(*it) != 0)) {
      ++it;
    }
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 =
        std::search(it, it_end, test_seq.begin(), test_seq.end());

    for (std::uint64_t i = 0; i < scheme->get_field_count(); ++i) {
      std::pair<std::string, std::shared_ptr<type_scheme>> fld =
          scheme->get_field(i);
      it = mark_field(it, it_end2, fld.first, fld.second);
    }

    if (it_end2 != it_end) {
      it_end2 += test_seq.length();
    }
    it_prev = it_end2;
  } else if (scheme->get_object_type() ==
             serializable_ptr_scheme::get_static_object_type()) {
    std::string test_seq = "<" + fld_name;
    std::string::iterator it =
        std::search(it_prev, it_end, test_seq.begin(), test_seq.end());
    if (it == it_end) {
      return it_prev;
    }
    src_markers.push_back(it - xml_src.begin());
    field_schemes.push_back(scheme);
    field_names.push_back(fld_name);
    test_seq = "</" + fld_name + ">";
    std::string::iterator it_end2 =
        std::search(it, it_end, test_seq.begin(), test_seq.end());
    if (it_end2 != it_end) {
      it_end2 += test_seq.length();
    }
    it_prev = it_end2;
  } else if (scheme->get_object_type() ==
             vector_type_scheme::get_static_object_type()) {
    std::string test_seq = "<" + fld_name + "_count";
    std::string::iterator it =
        std::search(it_prev, it_end, test_seq.begin(), test_seq.end());
    if (it == it_end) {
      return it_prev;
    }
    src_markers.push_back(it - xml_src.begin());
    field_schemes.push_back(
        std::shared_ptr<type_scheme>(new primitive_scheme<std::uint64_t>()));
    field_names.push_back(fld_name + "_count");
    test_seq = "</" + fld_name + "_count>";
    std::string::iterator it_end2 =
        std::search(it, it_end, test_seq.begin(), test_seq.end());
    if (it_end2 != it_end) {
      it_end2 += test_seq.length();
    }
    it_prev = it_end2;

    std::uint64_t i = 0;
    std::stringstream ss;
    ss << fld_name << "_q[" << i << "]";
    it = mark_field(it_prev, it_end, ss.str(), scheme->get_field(0).second);
    while (it != it_prev) {
      it_prev = it;
      ++i;
      ss.str("");
      ss << fld_name << "_q[" << i << "]";
      it = mark_field(it_prev, it_end, ss.str(), scheme->get_field(0).second);
    }
  } else if (scheme->get_object_type() ==
             map_type_scheme::get_static_object_type()) {
    std::string test_seq = "<" + fld_name + "_count";
    std::string::iterator it =
        std::search(it_prev, it_end, test_seq.begin(), test_seq.end());
    if (it == it_end) {
      return it_prev;
    }
    src_markers.push_back(it - xml_src.begin());
    field_schemes.push_back(
        std::shared_ptr<type_scheme>(new primitive_scheme<std::uint64_t>()));
    field_names.push_back(fld_name + "_count");
    test_seq = "</" + fld_name + "_count>";
    std::string::iterator it_end2 =
        std::search(it, it_end, test_seq.begin(), test_seq.end());
    if (it_end2 != it_end) {
      it_end2 += test_seq.length();
    }
    it_prev = it_end2;

    std::uint64_t i = 0;
    std::stringstream ss;
    ss << fld_name << "_key[" << i << "]";
    it = mark_field(it_prev, it_end, ss.str(), scheme->get_field(0).second);
    while (it != it_prev) {
      it_prev = it;
      ss.str("");
      ss << fld_name << "_value[" << i << "]";
      it = mark_field(it_prev, it_end, ss.str(), scheme->get_field(1).second);
      it_prev = it;
      ++i;
      ss.str("");
      ss << fld_name << "_key[" << i << "]";
      it = mark_field(it_prev, it_end, ss.str(), scheme->get_field(0).second);
    }
  }
  return it_prev;
}

std::uint64_t xml_field_editor::get_field_index(const std::string& name) const {
  for (std::uint64_t i = 0; i < src_markers.size(); ++i) {
    if (name == field_names[i]) {
      return i;
    }
  }
  return src_markers.size();
}

std::shared_ptr<type_scheme> xml_field_editor::get_type_scheme() const {
  std::shared_ptr<serializable> cur_ptr =
      (*(p_parent->get_object_graph()))[node].p_obj;
  if (!cur_ptr) {
    return {};
  }
  std::string t_name = cur_ptr->get_object_type()->name();
  if (t_name.empty()) {
    return {};
  }
  auto itm = get_global_schemes().find(t_name);
  if ((itm == get_global_schemes().end()) || (!(itm->second))) {
    return {};
  }
  return itm->second;
}

const std::string& xml_field_editor::get_complete_src() const {
  return (*(p_parent->get_object_graph()))[node].xml_src;
}

void xml_field_editor::set_complete_src(const std::string& xml_src) {
  (*(p_parent->get_object_graph()))[node].xml_src = xml_src;
  p_parent->ot_input_arc.load_current_from_node(node);
  std::shared_ptr<serializable> cur_ptr =
      (*(p_parent->get_object_graph()))[node].p_obj;
  src_markers.clear();
  field_schemes.clear();
  field_names.clear();
  if (!cur_ptr) {
    return;
  }
  cur_ptr->load(p_parent->ot_input_arc, cur_ptr->get_object_type()->version());

  std::string t_name = cur_ptr->get_object_type()->name();
  if (t_name.empty()) {
    return;
  }
  auto itm = get_global_schemes().find(t_name);
  if ((itm == get_global_schemes().end()) || (!(itm->second))) {
    return;
  }
  std::string::iterator it_prev =
      (*(p_parent->get_object_graph()))[node].xml_src.begin();
  for (std::uint64_t i = 0; i < itm->second->get_field_count(); ++i) {
    std::pair<std::string, std::shared_ptr<type_scheme>> fld =
        itm->second->get_field(i);
    it_prev = mark_field(it_prev,
                         (*(p_parent->get_object_graph()))[node].xml_src.end(),
                         fld.first, fld.second);
  }
}

std::string xml_field_editor::get_object_name(object_node_desc a_node) const {
  return p_parent->get_object_name(a_node);
}

std::string xml_field_editor::get_object_name() const {
  return p_parent->get_object_name(node);
}

xml_field_editor::xml_field_editor(objtree_editor* a_parent,
                                   object_node_desc a_node)
    : p_parent(a_parent), node(a_node) {
  std::shared_ptr<serializable> cur_ptr =
      (*(p_parent->get_object_graph()))[node].p_obj;
  if (!cur_ptr) {
    return;
  }
  std::string t_name = cur_ptr->get_object_type()->name();
  if (t_name.empty()) {
    return;
  }
  auto itm = get_global_schemes().find(t_name);
  if ((itm == get_global_schemes().end()) || (!(itm->second))) {
    return;
  }
  std::string::iterator it_prev =
      (*(p_parent->get_object_graph()))[node].xml_src.begin();
  for (std::uint64_t i = 0; i < itm->second->get_field_count(); ++i) {
    std::pair<std::string, std::shared_ptr<type_scheme>> fld =
        itm->second->get_field(i);
    it_prev = mark_field(it_prev,
                         (*(p_parent->get_object_graph()))[node].xml_src.end(),
                         fld.first, fld.second);
  }
}

std::uint64_t xml_field_editor::get_total_field_count() const {
  return src_markers.size();
}

std::pair<std::string, std::shared_ptr<type_scheme>>
xml_field_editor::get_field(std::uint64_t index) const {
  return {field_names[index], field_schemes[index]};
}

std::string xml_field_editor::get_field_src(std::uint64_t index) const {
  if (index >= src_markers.size()) {
    return "";
  }
  if (index + 1 == src_markers.size()) {
    return {(*(p_parent->get_object_graph()))[node].xml_src.begin() +
                src_markers[index],
            (*(p_parent->get_object_graph()))[node].xml_src.end()};
  }
  return {(*(p_parent->get_object_graph()))[node].xml_src.begin() +
              src_markers[index],
          (*(p_parent->get_object_graph()))[node].xml_src.begin() +
              src_markers[index + 1]};
}

std::string xml_field_editor::get_field_src(const std::string& name) const {
  return get_field_src(get_field_index(name));
}

std::string xml_field_editor::get_field_value(std::uint64_t index) const {
  std::string& xml_src = (*(p_parent->get_object_graph()))[node].xml_src;
  std::string::iterator it = xml_src.begin() + src_markers[index];
  std::string test_str;
  if (field_schemes[index]->get_object_type() ==
      serializable_ptr_scheme::get_static_object_type()) {
    test_str = "object_ID=\"";
  } else {
    test_str = ">\"";
  }

  it = std::search(it, xml_src.end(), test_str.begin(), test_str.end());
  if (it == xml_src.end()) {
    return "";
  }
  it += test_str.length();
  std::string::iterator it_end = std::find(it, xml_src.end(), '\"');
  if (field_schemes[index]->get_object_type() ==
      serializable_ptr_scheme::get_static_object_type()) {
    object_node_desc fld_node = 0;
    std::stringstream(std::string(it, it_end)) >> fld_node;
    return get_object_name(fld_node);
  }
  return {it, it_end};
}

std::string xml_field_editor::get_field_value(const std::string& name) const {
  return get_field_value(get_field_index(name));
}

void xml_field_editor::set_field_value(std::uint64_t index,
                                       const std::string& value) {
  std::string& xml_src = (*(p_parent->get_object_graph()))[node].xml_src;
  std::string::iterator it = xml_src.begin() + src_markers[index];
  std::string test_str;
  if (field_schemes[index]->get_object_type() ==
      serializable_ptr_scheme::get_static_object_type()) {
    test_str = "object_ID=\"";
  } else {
    test_str = ">\"";
  }
  it = std::search(it, xml_src.end(), test_str.begin(), test_str.end());
  if (it == xml_src.end()) {
    return;
  }
  it += test_str.length();
  std::string::iterator it_end = std::find(it, xml_src.end(), '\"');
  std::uint64_t orig_len = it_end - it;
  std::string new_xml_src(xml_src.begin(), it);
  std::int64_t len_diff = value.length() - orig_len;

  if (field_schemes[index]->get_object_type() ==
      serializable_ptr_scheme::get_static_object_type()) {
    object_node_desc orig_node = 0;
    object_node_desc new_node = 0;
    std::stringstream(std::string(it, it_end)) >> orig_node;
    {
      std::string test_str3 = "(ID:";
      std::string::const_iterator it3 = std::search(
          value.begin(), value.end(), test_str3.begin(), test_str3.end());
      if (it3 == value.end()) {
        new_node = 0;
      } else {
        it3 += test_str3.length();
        std::string::const_iterator it3_end = std::find(it3, value.end(), ')');
        std::stringstream(std::string(it3, it3_end)) >> new_node;
      }
    }
    // first, check if the original node appeared anywhere else in the same xml-source.
    std::uint64_t orig_count = 0;
    {
      std::stringstream ss2;
      ss2 << "object_ID=\"" << orig_node << "\"";
      std::string test_str2 = ss2.str();
      std::string::iterator it2 = std::search(
          xml_src.begin(), xml_src.end(), test_str2.begin(), test_str2.end());
      while (it2 != xml_src.end()) {
        ++orig_count;
        ++it2;
        it2 =
            std::search(it2, xml_src.end(), test_str2.begin(), test_str2.end());
      }
    }
    // then, check if the new node appears anywhere in the xml-source already.
    std::uint64_t new_count = 0;
    {
      std::stringstream ss2;
      ss2 << "object_ID=\"" << new_node << "\"";
      std::string test_str2 = ss2.str();
      std::string::iterator it2 = std::search(
          xml_src.begin(), xml_src.end(), test_str2.begin(), test_str2.end());
      while (it2 != xml_src.end()) {
        ++orig_count;
        ++it2;
        it2 =
            std::search(it2, xml_src.end(), test_str2.begin(), test_str2.end());
      }
    }
    if (new_count == 0) {
      if (orig_count == 1) {  // the nodes must be swapped.
        p_parent->replace_child(node, new_node, orig_node);
      } else {  // the new-node must be linked to the parent without removing the existing link.
        p_parent->create_child(node, new_node);
      }
    } else if (orig_count == 1) {
      p_parent->sever_child(node, orig_node);
    }

    std::stringstream ss4;
    ss4 << new_node;
    std::string new_node_str = ss4.str();
    len_diff = new_node_str.length() - orig_len;
    new_xml_src.append(new_node_str);
  } else {
    new_xml_src.append(value);
  }
  new_xml_src.append(it_end, xml_src.end());
  for (std::uint64_t j = index; j < src_markers.size(); ++j) {
    src_markers[j] += len_diff;
  }
  xml_src = std::move(new_xml_src);
  p_parent->ot_input_arc.load_current_from_node(node);
  std::shared_ptr<serializable> cur_ptr =
      (*(p_parent->get_object_graph()))[node].p_obj;
  if (cur_ptr) {
    cur_ptr->load(p_parent->ot_input_arc,
                  cur_ptr->get_object_type()->version());
  }
}

void xml_field_editor::set_field_value(const std::string& name,
                                       const std::string& value) {
  set_field_value(get_field_index(name), value);
}

void xml_field_editor::set_field_newptr(
    std::uint64_t index, const std::shared_ptr<serializable>& new_ptr) {
  if (field_schemes[index]->get_object_type() !=
      serializable_ptr_scheme::get_static_object_type()) {
    return;
  }
  std::string& xml_src = (*(p_parent->get_object_graph()))[node].xml_src;

  std::string::iterator it = xml_src.begin() + src_markers[index];
  std::string test_str = "object_ID=\"";
  it = std::search(it, xml_src.end(), test_str.begin(), test_str.end());
  if (it == xml_src.end()) {
    return;
  }
  it += test_str.length();
  std::string::iterator it_end = std::find(it, xml_src.end(), '\"');

  object_node_desc orig_node = 0;
  std::stringstream(std::string(it, it_end)) >> orig_node;
  // first, check if the original node appeared anywhere else in the same xml-source.
  std::uint64_t orig_count = 0;
  {
    std::stringstream ss2;
    ss2 << "object_ID=\"" << orig_node << "\"";
    std::string test_str2 = ss2.str();
    std::string::iterator it2 = std::search(xml_src.begin(), xml_src.end(),
                                            test_str2.begin(), test_str2.end());
    while (it2 != xml_src.end()) {
      ++orig_count;
      ++it2;
      it2 = std::search(it2, xml_src.end(), test_str2.begin(), test_str2.end());
    }
  }

  object_node_desc new_node = 0;
  if (orig_count == 1) {  // the original node must be replaced.
    new_node = p_parent->add_new_object(new_ptr, node, orig_node);
  } else {  // the new-node must be linked to the parent without removing the existing link.
    new_node = p_parent->add_new_object(new_ptr, node);
  }

  std::uint64_t orig_len = it_end - it;
  std::stringstream ss3;
  ss3 << new_node;
  std::string value = ss3.str();
  std::int64_t len_diff = value.length() - orig_len;
  std::string new_xml_src(xml_src.begin(), it);
  new_xml_src.append(value);
  new_xml_src.append(it_end, xml_src.end());
  for (std::uint64_t j = index; j < src_markers.size(); ++j) {
    src_markers[j] += len_diff;
  }
  xml_src = std::move(new_xml_src);
  p_parent->ot_input_arc.load_current_from_node(node);
  std::shared_ptr<serializable> cur_ptr =
      (*(p_parent->get_object_graph()))[node].p_obj;
  if (cur_ptr) {
    cur_ptr->load(p_parent->ot_input_arc,
                  cur_ptr->get_object_type()->version());
  }
}

void xml_field_editor::set_field_newptr(
    const std::string& name, const std::shared_ptr<serializable>& new_ptr) {
  set_field_newptr(get_field_index(name), new_ptr);
}

char objtree_iarchive::get_next_char() {
  char c = '\0';
  current_ss->get(c);
  while ((c == ' ') || (c == '\t') || (c == '\n') || (c == '\r')) {
    current_ss->get(c);
  }
  return c;
}

std::string objtree_iarchive::read_token() {
  std::string result;
  char c = get_next_char();
  if (c != '<') {
    return result;
  }
  c = get_next_char();
  if ((c == '!') || (c == '?')) {
    std::array<char, 512> line_str{};
    current_ss->getline(line_str.data(), line_str.size());
    return read_token();
  }
  while (c != '>') {
    result += c;
    current_ss->get(c);
  }
  return result;
}

void objtree_iarchive::skip_to_end_token(const std::string& name) {
  std::string token = read_token();
  trim_str(token);
  while (token != "/" + name) {
    token = read_token();
    trim_str(token);
  }
}

void objtree_iarchive::trim_str(std::string& s) {
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

bool objtree_iarchive::read_named_value(const std::string& value_name,
                                        std::string& value_str) {
  std::string token = read_token();
  trim_str(token);
  if ((value_name.empty()) || (token != value_name)) {
    return false;
  }

  char c = '\0';
  current_ss->get(c);
  while (c != '\"') {
    current_ss->get(c);
  }

  value_str.clear();
  current_ss->get(c);
  while (c != '\"') {
    value_str += c;
    current_ss->get(c);
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

archive_object_header objtree_iarchive::read_header(
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

  result.is_external = false;

  return result;
}

void objtree_iarchive::load_current_from_node(object_node_desc a_node) {
  current_ss =
      std::make_shared<std::stringstream>((*obj_graph)[a_node].xml_src);
}

objtree_iarchive::objtree_iarchive(std::shared_ptr<object_graph> a_obj_graph,
                                   object_node_desc a_root)
    : obj_graph(std::move(a_obj_graph)), obj_graph_root(a_root) {

  if (num_vertices(*obj_graph) == 0) {
    obj_graph_root = add_vertex(
        *obj_graph);  // add a root node. This case doesn't make much sense, it means the graph is empty.
  }

  current_ss =
      std::make_shared<std::stringstream>((*obj_graph)[obj_graph_root].xml_src);
}

objtree_iarchive::~objtree_iarchive() = default;

iarchive& objtree_iarchive::load_serializable_ptr(
    serializable_shared_pointer& item) {
  return objtree_iarchive::load_serializable_ptr(
      std::pair<std::string, serializable_shared_pointer&>("Item", item));
}

iarchive& objtree_iarchive::load_serializable_ptr(
    const std::pair<std::string, serializable_shared_pointer&>& item) {
  item.second = serializable_shared_pointer();
  using Vertex = bagl::graph_vertex_descriptor_t<object_graph>;

  std::vector<std::uint32_t> type_id;
  archive_object_header hdr = read_header(item.first, type_id);
  if ((type_id.empty()) || (hdr.type_version == 0) || (hdr.object_ID == 0)) {
    skip_to_end_token(item.first);
    return *this;
  }

  if ((hdr.object_ID < num_vertices(*obj_graph)) &&
      ((*obj_graph)[static_cast<Vertex>(hdr.object_ID)].p_obj)) {
    item.second = (*obj_graph)[static_cast<Vertex>(hdr.object_ID)].p_obj;
    skip_to_end_token(item.first);

    // re-read the xml source
    std::shared_ptr<std::stringstream> tmp_ss = current_ss;
    current_ss = std::make_shared<std::stringstream>(
        (*obj_graph)[static_cast<Vertex>(hdr.object_ID)].xml_src);

    item.second->load(*this, hdr.type_version);

    current_ss = tmp_ss;
  }

  return *this;
}

iarchive& objtree_iarchive::load_serializable(serializable& item) {
  return objtree_iarchive::load_serializable(
      std::pair<std::string, serializable&>("Item", item));
}

iarchive& objtree_iarchive::load_serializable(
    const std::pair<std::string, serializable&>& item) {
  archive_object_header hdr;

  std::vector<std::uint32_t> type_id;
  hdr = read_header(item.first, type_id);
  if ((hdr.type_ID == nullptr) || (hdr.type_version == 0)) {
    skip_to_end_token(item.first);
    return *this;
  }

  item.second.load(*this, hdr.type_version);

  skip_to_end_token(item.first);
  return *this;
}

iarchive& objtree_iarchive::load_char(char& i) {
  return objtree_iarchive::load_char(std::pair<std::string, char&>("char", i));
}

iarchive& objtree_iarchive::load_char(const std::pair<std::string, char&>& i) {
  std::string value_str;
  if (read_named_value(i.first, value_str)) {
    if (value_str.empty()) {
      i.second = 0;
    } else {
      int temp = 0;
      std::stringstream(value_str) >> temp;
      i.second = static_cast<char>(temp);
    };
  } else {
    i.second = 0;
  }
  return *this;
}

iarchive& objtree_iarchive::load_unsigned_char(unsigned char& u) {
  return objtree_iarchive::load_unsigned_char(
      std::pair<std::string, unsigned char&>("unsigned_char", u));
}

iarchive& objtree_iarchive::load_unsigned_char(
    const std::pair<std::string, unsigned char&>& u) {
  std::string value_str;
  if (read_named_value(u.first, value_str)) {
    if (value_str.empty()) {
      u.second = 0;
    } else {
      std::uint64_t temp = 0;
      std::stringstream(value_str) >> temp;
      u.second = static_cast<char>(temp);
    }
  } else {
    u.second = 0;
  }
  return *this;
}

iarchive& objtree_iarchive::load_int(std::int64_t& i) {
  return objtree_iarchive::load_int(
      std::pair<std::string, std::int64_t&>("int", i));
}

iarchive& objtree_iarchive::load_int(
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

iarchive& objtree_iarchive::load_unsigned_int(std::uint64_t& u) {
  return objtree_iarchive::load_unsigned_int(
      std::pair<std::string, std::uint64_t&>("unsigned_int", u));
}

iarchive& objtree_iarchive::load_unsigned_int(
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

iarchive& objtree_iarchive::load_float(float& f) {
  return objtree_iarchive::load_float(
      std::pair<std::string, float&>("real", f));
}

iarchive& objtree_iarchive::load_float(
    const std::pair<std::string, float&>& f) {
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

iarchive& objtree_iarchive::load_double(double& d) {
  return objtree_iarchive::load_double(
      std::pair<std::string, double&>("real", d));
}

iarchive& objtree_iarchive::load_double(
    const std::pair<std::string, double&>& d) {
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

iarchive& objtree_iarchive::load_bool(bool& b) {
  return objtree_iarchive::load_bool(std::pair<std::string, bool&>("bool", b));
}

iarchive& objtree_iarchive::load_bool(const std::pair<std::string, bool&>& b) {
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

iarchive& objtree_iarchive::load_string(std::string& s) {
  return objtree_iarchive::load_string(
      std::pair<std::string, std::string&>("string", s));
}

iarchive& objtree_iarchive::load_string(
    const std::pair<std::string, std::string&>& s) {
  read_named_value(s.first, s.second);
  return *this;
}

void objtree_oarchive::register_new_object(object_node_desc a_node) {
  mObjRegMap[(*obj_graph)[a_node].p_obj] = static_cast<std::uint64_t>(a_node);
}

void objtree_oarchive::unregister_object(object_node_desc a_node) {
  mObjRegMap.erase((*obj_graph)[a_node].p_obj);
}

void objtree_oarchive::save_current_stream() {
  (*obj_graph)[current_node].xml_src = current_ss->str();
}

void objtree_oarchive::load_current_from_node(object_node_desc a_node) {
  current_ss =
      std::make_shared<std::stringstream>((*obj_graph)[a_node].xml_src);
  current_node = a_node;
}

void objtree_oarchive::fresh_current_node(object_node_desc a_node) {
  current_ss = std::make_shared<std::stringstream>();
  current_node = a_node;
}

objtree_oarchive::objtree_oarchive(std::shared_ptr<object_graph> a_obj_graph,
                                   object_node_desc a_root)
    : obj_graph(std::move(a_obj_graph)),
      obj_graph_root(a_root),
      current_node(obj_graph_root) {
  current_ss = std::make_shared<std::stringstream>();

  for (auto v : vertices(*obj_graph)) {
    if ((*obj_graph)[v].p_obj) {
      mObjRegMap[(*obj_graph)[v].p_obj] = static_cast<std::uint64_t>(v);
    }
  }
}

objtree_oarchive::~objtree_oarchive() {
  save_current_stream();
}

oarchive& objtree_oarchive::save_to_new_archive_impl(
    const serializable_shared_pointer& item, const std::string& /*file_name*/) {
  return objtree_oarchive::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("Item", item));
}

oarchive& objtree_oarchive::save_to_new_archive_named_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& item,
    const std::string& /*file_name*/) {
  return objtree_oarchive::save_serializable_ptr(item);
}

oarchive& objtree_oarchive::save_serializable_ptr(
    const serializable_shared_pointer& item) {
  return objtree_oarchive::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("Item", item));
}

oarchive& objtree_oarchive::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& item) {

  if (item.second) {
    auto it = mObjRegMap.find(item.second);
    std::uint64_t object_id = 0;
    if (it != mObjRegMap.end()) {
      object_id = it->second;
    } else {
      object_id = static_cast<std::uint64_t>(add_vertex(*obj_graph));
      (*obj_graph)[object_id].p_obj = item.second;
      mObjRegMap[item.second] = object_id;
    }

    if (!edge(current_node, object_id, *obj_graph).second) {
      add_edge(current_node, object_id, *obj_graph);
    }

    rtti::so_type* obj_type = item.second->get_object_type();
    const std::uint32_t* type_id = obj_type->id_begin();
    std::uint32_t type_version = obj_type->version();

    (*current_ss) << "<" << item.first << " type_ID=\"";
    while (*type_id != 0U) {
      (*current_ss) << *type_id << ".";
      ++type_id;
    }
    (*current_ss) << "0\" version=\"" << type_version << "\" object_ID=\""
                  << object_id << "\"></" << item.first << ">\n";

    std::shared_ptr<std::stringstream> tmp = current_ss;
    current_ss = std::make_shared<std::stringstream>();

    auto tmp_v = current_node;
    current_node = object_id;

    item.second->save(*this, type_version);

    // grab the resulting string.
    (*obj_graph)[object_id].xml_src = current_ss->str();

    current_ss = tmp;
    current_node = tmp_v;
  } else {
    (*current_ss) << "<" << item.first
                  << R"( type_ID="0" version="0" object_ID="0"></)"
                  << item.first << ">\n";
  }

  return *this;
}

oarchive& objtree_oarchive::save_serializable(const serializable& item) {
  return *this & std::pair<std::string, const serializable&>("Item", item);
}

oarchive& objtree_oarchive::save_serializable(
    const std::pair<std::string, const serializable&>& item) {
  archive_object_header hdr;
  const std::uint32_t* type_id = item.second.get_object_type()->id_begin();
  hdr.type_version = item.second.get_object_type()->version();
  hdr.object_ID = 0;
  hdr.size = 0;
  hdr.is_external = false;

  (*current_ss) << "<" << item.first << " type_ID=\"";
  while (*type_id != 0U) {
    (*current_ss) << *type_id << ".";
    ++type_id;
  };
  (*current_ss) << "0\" version=\"" << hdr.type_version << "\">\n";

  item.second.save(*this, hdr.type_version);

  (*current_ss) << "</" << item.first << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_char(char i) {
  return objtree_oarchive::save_char(std::pair<std::string, char>("char", i));
}

oarchive& objtree_oarchive::save_char(const std::pair<std::string, char>& i) {
  (*current_ss) << "<" << i.first << ">\"" << static_cast<int>(i.second)
                << "\"</" << i.first << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_unsigned_char(unsigned char u) {
  return objtree_oarchive::save_unsigned_char(
      std::pair<std::string, unsigned char>("unsigned_char", u));
}

oarchive& objtree_oarchive::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  (*current_ss) << "<" << u.first << ">\""
                << static_cast<std::uint32_t>(u.second) << "\"</" << u.first
                << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_int(std::int64_t i) {
  return objtree_oarchive::save_int(
      std::pair<std::string, std::int64_t>("int", i));
}

oarchive& objtree_oarchive::save_int(
    const std::pair<std::string, std::int64_t>& i) {
  (*current_ss) << "<" << i.first << ">\"" << i.second << "\"</" << i.first
                << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_unsigned_int(std::uint64_t u) {
  return objtree_oarchive::save_unsigned_int(
      std::pair<std::string, std::uint64_t>("unsigned_int", u));
}

oarchive& objtree_oarchive::save_unsigned_int(
    const std::pair<std::string, std::uint64_t>& u) {
  (*current_ss) << "<" << u.first << ">\"" << u.second << "\"</" << u.first
                << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_float(float f) {
  return objtree_oarchive::save_float(std::pair<std::string, float>("real", f));
}

oarchive& objtree_oarchive::save_float(const std::pair<std::string, float>& f) {
  (*current_ss) << "<" << f.first << ">\"" << f.second << "\"</" << f.first
                << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_double(double d) {
  return objtree_oarchive::save_double(
      std::pair<std::string, double>("real", d));
}

oarchive& objtree_oarchive::save_double(
    const std::pair<std::string, double>& d) {
  (*current_ss) << "<" << d.first << ">\"" << d.second << "\"</" << d.first
                << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_bool(bool b) {
  return objtree_oarchive::save_bool(std::pair<std::string, bool>("bool", b));
}

oarchive& objtree_oarchive::save_bool(const std::pair<std::string, bool>& b) {
  (*current_ss) << "<" << b.first << ">\"" << (b.second ? "true" : "false")
                << "\"</" << b.first << ">\n";
  return *this;
}

oarchive& objtree_oarchive::save_string(const std::string& s) {
  return objtree_oarchive::save_string(
      std::pair<std::string, const std::string&>("string", s));
}

oarchive& objtree_oarchive::save_string(
    const std::pair<std::string, const std::string&>& s) {
  (*current_ss) << "<" << s.first << ">\"" << s.second << "\"</" << s.first
                << ">\n";
  return *this;
}

objtree_editor::objtree_editor()
    : obj_graph(new object_graph()),
      obj_graph_root(add_vertex(*obj_graph)),
      ot_output_arc(obj_graph, obj_graph_root),
      ot_input_arc(obj_graph, obj_graph_root) {}

objtree_editor::objtree_editor(std::shared_ptr<object_graph> a_obj_graph,
                               object_node_desc a_root)
    : obj_graph(std::move(a_obj_graph)),
      obj_graph_root(a_root),
      ot_output_arc(obj_graph, obj_graph_root),
      ot_input_arc(obj_graph, obj_graph_root) {}

object_node_desc objtree_editor::add_new_object(
    const std::shared_ptr<serializable>& new_obj, object_node_desc a_parent,
    object_node_desc old_child) {
  if (!new_obj) {
    return obj_graph_root;
  }
  object_node_desc result = 0;
  if (obj_graph_graveyard.empty()) {
    result = add_vertex(*obj_graph, object_graph_node(new_obj, ""));
  } else {
    result = obj_graph_graveyard.top();
    obj_graph_graveyard.pop();
    (*obj_graph)[result].p_obj = new_obj;
    (*obj_graph)[result].xml_src = "";
  }
  if (old_child != 0U) {
    remove_edge(a_parent, old_child, *obj_graph);
    remove_object(old_child);  // if this was the last in-edge on old_child.
  }
  add_edge(a_parent, result, *obj_graph);
  ot_output_arc.register_new_object(result);
  ot_output_arc.fresh_current_node(result);
  new_obj->save(ot_output_arc, new_obj->get_object_type()->version());
  ot_output_arc.save_current_stream();
  return result;
}

void objtree_editor::remove_object(object_node_desc a_node) {
  if (in_degree(a_node, *obj_graph) != 0U) {
    return;
  }
  std::vector<object_node_desc> v;
  v.reserve(out_degree(a_node, *obj_graph));
  for (auto e : out_edges(a_node, *obj_graph)) {
    v.push_back(target(e, *obj_graph));
  }
  clear_vertex(a_node, *obj_graph);
  for (std::uint64_t& it : v) {
    remove_object(
        it);  // this will only really have an effect if the node has no other in-edge.
  }
  ot_output_arc.unregister_object(a_node);
  (*obj_graph)[a_node].p_obj = std::shared_ptr<serializable>();
  (*obj_graph)[a_node].xml_src = "";
  obj_graph_graveyard.push(a_node);
}

void objtree_editor::replace_child(object_node_desc a_parent,
                                   object_node_desc new_child,
                                   object_node_desc old_child) {
  if (((new_child) != 0U) && (!edge(a_parent, new_child, *obj_graph).second)) {
    add_edge(a_parent, new_child, *obj_graph);
  }
  if (old_child != 0U) {
    remove_edge(a_parent, old_child, *obj_graph);
    // this will only really have an effect if the node has no other in-edge.
    remove_object(old_child);
  }
}

std::string objtree_editor::get_object_name(object_node_desc a_node) const {
  return get_objtree_name(*obj_graph, a_node);
}

std::vector<std::string> objtree_editor::get_objects_derived_from(
    rtti::so_type* aType) const {
  std::vector<std::string> result;
  for (auto v : vertices(*obj_graph)) {
    std::shared_ptr<serializable> p_obj = (*obj_graph)[v].p_obj;
    if ((p_obj) && ((p_obj->cast_to(aType)) != nullptr)) {
      result.push_back(get_object_name(v));
    }
  }
  return result;
}

std::string get_objtree_name(const object_graph& obj_graph,
                             object_node_desc node_id) {
  std::shared_ptr<named_object> item_ptr =
      rtti::rk_dynamic_ptr_cast<named_object>(obj_graph[node_id].p_obj);
  std::stringstream ss;
  if (item_ptr) {
    ss << item_ptr->get_name() << " (ID:" << node_id << ")";
  } else {
    ss << "Object (ID:" << node_id << ")";
  }
  return ss.str();
}

object_node_desc get_objtree_node_id(const object_graph& /*obj_graph*/,
                                     const std::string& obj_name) {
  std::string node_id_str = obj_name.substr(obj_name.find("(ID:") + 4);
  node_id_str = node_id_str.substr(0, node_id_str.find(')'));
  std::stringstream ss(node_id_str);
  std::uint64_t result = 0;
  ss >> result;
  return result;
}

}  // namespace ReaK::serialization
