/**
 * \file objtree_archiver.h
 *
 * This library declares the class for a creating type schemes that represent the fields contained in the
 * serialization of a given type.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date November 2012
 */

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

#ifndef REAK_CORE_SERIALIZATION_OBJTREE_ARCHIVER_H_
#define REAK_CORE_SERIALIZATION_OBJTREE_ARCHIVER_H_

#include "ReaK/core/serialization/archiver.h"

#include <queue>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"

#include "bagl/adjacency_list.h"

namespace ReaK {

namespace rtti {
class so_type;
}  // namespace rtti

namespace serialization {

class objtree_editor;
class type_scheme;

struct object_graph_node {
  std::shared_ptr<serializable> p_obj;
  std::string xml_src;

  explicit object_graph_node(
      std::shared_ptr<serializable> pobj = std::shared_ptr<serializable>(),
      std::string xmlsrc = "")
      : p_obj(std::move(pobj)), xml_src(std::move(xmlsrc)){};
};

using object_graph =
    bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::bidirectional_s,
                         object_graph_node>;
using object_node_desc = bagl::graph_vertex_descriptor_t<object_graph>;

class xml_field_editor {
 private:
  objtree_editor* p_parent;
  object_node_desc node;
  std::vector<std::uint64_t> src_markers;
  std::vector<std::shared_ptr<type_scheme>> field_schemes;
  std::vector<std::string> field_names;

  std::string::iterator mark_field(
      std::string::iterator it_prev, std::string::iterator it_end,
      const std::string& fld_name,
      const std::shared_ptr<ReaK::serialization::type_scheme>& scheme);
  [[nodiscard]] std::uint64_t get_field_index(const std::string& name) const;
  [[nodiscard]] std::string get_object_name(object_node_desc a_node) const;

 public:
  friend class objtree_editor;

  [[nodiscard]] std::shared_ptr<type_scheme> get_type_scheme() const;
  [[nodiscard]] const std::string& get_complete_src() const;
  void set_complete_src(const std::string& xml_src);
  [[nodiscard]] std::string get_object_name() const;

  xml_field_editor(objtree_editor* a_parent, object_node_desc a_node);

  [[nodiscard]] std::uint64_t get_total_field_count() const;
  [[nodiscard]] std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::uint64_t index) const;

  [[nodiscard]] std::string get_field_src(std::uint64_t index) const;
  [[nodiscard]] std::string get_field_src(const std::string& name) const;

  [[nodiscard]] std::string get_field_value(std::uint64_t index) const;
  [[nodiscard]] std::string get_field_value(const std::string& name) const;

  void set_field_value(std::uint64_t index, const std::string& value);
  void set_field_value(const std::string& name, const std::string& value);

  void set_field_newptr(std::uint64_t index,
                        const std::shared_ptr<serializable>& new_ptr);
  void set_field_newptr(const std::string& name,
                        const std::shared_ptr<serializable>& new_ptr);
};

/**
 * An input archive used to generate an object tree.
 */
class objtree_iarchive : public iarchive {
 private:
  std::shared_ptr<object_graph> obj_graph;
  std::shared_ptr<std::stringstream> current_ss;
  object_node_desc obj_graph_root;

  char get_next_char();
  std::string read_token();
  void skip_to_end_token(const std::string& name);
  static void trim_str(std::string& s);
  bool read_named_value(const std::string& value_name, std::string& value_str);
  archive_object_header read_header(const std::string& obj_name,
                                    std::vector<std::uint32_t>& out_type_id);

 protected:
  iarchive& load_serializable_ptr(serializable_shared_pointer& item) override;

  iarchive& load_serializable_ptr(
      const std::pair<std::string, serializable_shared_pointer&>& item)
      override;

  iarchive& load_serializable(serializable& item) override;

  iarchive& load_serializable(
      const std::pair<std::string, serializable&>& item) override;

  iarchive& load_char(char& i) override;

  iarchive& load_char(const std::pair<std::string, char&>& i) override;

  iarchive& load_unsigned_char(unsigned char& u) override;

  iarchive& load_unsigned_char(
      const std::pair<std::string, unsigned char&>& u) override;

  iarchive& load_int(std::int64_t& i) override;

  iarchive& load_int(const std::pair<std::string, std::int64_t&>& i) override;

  iarchive& load_unsigned_int(std::uint64_t& u) override;

  iarchive& load_unsigned_int(
      const std::pair<std::string, std::uint64_t&>& u) override;

  iarchive& load_float(float& f) override;

  iarchive& load_float(const std::pair<std::string, float&>& f) override;

  iarchive& load_double(double& d) override;

  iarchive& load_double(const std::pair<std::string, double&>& d) override;

  iarchive& load_bool(bool& b) override;

  iarchive& load_bool(const std::pair<std::string, bool&>& b) override;

  iarchive& load_string(std::string& s) override;

  iarchive& load_string(const std::pair<std::string, std::string&>& s) override;

  void load_current_from_node(object_node_desc a_node);

 public:
  friend class objtree_editor;
  friend class xml_field_editor;

  [[nodiscard]] std::shared_ptr<object_graph> get_object_graph() const {
    return obj_graph;
  }
  [[nodiscard]] object_node_desc get_root_node() const {
    return obj_graph_root;
  }

  explicit objtree_iarchive(std::shared_ptr<object_graph> a_obj_graph,
                            object_node_desc a_root = object_node_desc{});
  ~objtree_iarchive() override;
};

/**
 * XML output archive.
 */
class objtree_oarchive : public oarchive {
 private:
  std::shared_ptr<object_graph> obj_graph;
  object_node_desc obj_graph_root;
  std::shared_ptr<std::stringstream> current_ss;
  object_node_desc current_node;

 protected:
  oarchive& save_to_new_archive_impl(const serializable_shared_pointer& item,
                                     const std::string& file_name) override;

  oarchive& save_to_new_archive_named_impl(
      const std::pair<std::string, const serializable_shared_pointer&>& item,
      const std::string& file_name) override;

  oarchive& save_serializable_ptr(
      const serializable_shared_pointer& item) override;

  oarchive& save_serializable_ptr(
      const std::pair<std::string, const serializable_shared_pointer&>& item)
      override;

  oarchive& save_serializable(const serializable& item) override;

  oarchive& save_serializable(
      const std::pair<std::string, const serializable&>& item) override;

  oarchive& save_char(char i) override;

  oarchive& save_char(const std::pair<std::string, char>& i) override;

  oarchive& save_unsigned_char(unsigned char u) override;

  oarchive& save_unsigned_char(
      const std::pair<std::string, unsigned char>& u) override;

  oarchive& save_int(std::int64_t i) override;

  oarchive& save_int(const std::pair<std::string, std::int64_t>& i) override;

  oarchive& save_unsigned_int(std::uint64_t u) override;

  oarchive& save_unsigned_int(
      const std::pair<std::string, std::uint64_t>& u) override;

  oarchive& save_float(float f) override;

  oarchive& save_float(const std::pair<std::string, float>& f) override;

  oarchive& save_double(double d) override;

  oarchive& save_double(const std::pair<std::string, double>& d) override;

  oarchive& save_bool(bool b) override;

  oarchive& save_bool(const std::pair<std::string, bool>& b) override;

  oarchive& save_string(const std::string& s) override;

  oarchive& save_string(
      const std::pair<std::string, const std::string&>& s) override;

  void register_new_object(object_node_desc a_node);
  void unregister_object(object_node_desc a_node);
  void save_current_stream();
  void load_current_from_node(object_node_desc a_node);
  void fresh_current_node(object_node_desc a_node);

 public:
  friend class objtree_editor;

  std::shared_ptr<object_graph> get_object_graph() const { return obj_graph; }
  object_node_desc get_root_node() const { return obj_graph_root; }

  void set_current_node(object_node_desc a_node) { current_node = a_node; }
  object_node_desc get_current_node() const { return current_node; }

  explicit objtree_oarchive(std::shared_ptr<object_graph> a_obj_graph,
                            object_node_desc a_root = object_node_desc{});
  ~objtree_oarchive() override;
};

/**
 * This class acts as a manager or container for an object graph that is constantly
 * kept in sync with the objects that it handles. In other words, it overlays the
 * implicit object-graph formed by the relationships between the objects currently
 * loaded in software with explicit "object inspection" capabilities such as those
 * necessary for a generic editor for the objects.
 */
class objtree_editor {
 private:
  std::shared_ptr<object_graph> obj_graph;
  object_node_desc obj_graph_root;
  objtree_oarchive ot_output_arc;
  objtree_iarchive ot_input_arc;
  std::priority_queue<object_node_desc> obj_graph_graveyard;

 public:
  friend class xml_field_editor;

  objtree_editor(const objtree_editor&) = delete;             // non-copyable.
  objtree_editor& operator=(const objtree_editor&) = delete;  // non-assignable.

  [[nodiscard]] const std::shared_ptr<object_graph>& get_object_graph() const {
    return obj_graph;
  }
  [[nodiscard]] object_node_desc get_root_node() const {
    return obj_graph_root;
  }

  objtree_editor();
  explicit objtree_editor(std::shared_ptr<object_graph> a_obj_graph,
                          object_node_desc a_root = object_node_desc{});

  /**
   * This function adds a new object node to the object graph. In this version, the new object has a
   * given parent and is assumed to replace the 'old_child' node as child of the given parent. This
   * has the effect of severing that parent-child connection of the old child (unless null) and replacing
   * it with a parent-child connection for the new node.
   * \param new_obj The new object node to add to the object graph.
   * \return The vertex descriptor of the node within the object-graph.
   */
  object_node_desc add_new_object(
      const std::shared_ptr<serializable>& new_obj, object_node_desc a_parent,
      object_node_desc old_child = object_node_desc{});

  /**
   * This function adds a new object node to the object graph. In this version, the new object has no
   * parent, and is thus defaulted as being a child of the root node.
   * \param new_obj The new object node to add to the object graph.
   * \return The vertex descriptor of the node within the object-graph.
   */
  object_node_desc add_new_object(
      const std::shared_ptr<serializable>& new_obj) {
    return add_new_object(new_obj, obj_graph_root);
  }

  /**
   * This function attempts to remove the given object from the object graph. Note that this operation
   * cannot be performed unless the object has been severed from all its parent-child connections.
   * With shared-ownership semantics, this function should, in principle, trigger the actual deletion
   * of the object referred to by the node being deleted.
   * \param a_node The node to be removed from the graph.
   */
  void remove_object(object_node_desc a_node);

  /**
   * This function replaces (or reroutes) the object graph such that a new child node replaces the old child
   * in a parent-child connection (edge of the graph). Both the new or old child could be null nodes, meaning
   * that either a connection is severed or created, respectively.
   * \param a_parent The parent of which a child is swapped for another.
   * \param new_child The new child that will come and replace the current (old) child. If the new child is null, then
   * the parent-child connection is simply severed.
   * \param old_child The current (or old) child to be replaced. If the old child is null, then a new parent-child
   * connection is created.
   */
  void replace_child(object_node_desc a_parent, object_node_desc new_child,
                     object_node_desc old_child);

  /**
   * This function severs the parent-child connection (edge of the graph).
   * \param a_parent The parent from which a child is severed.
   * \param old_child The current (or old) child to be severed.
   */
  void sever_child(object_node_desc a_parent, object_node_desc old_child) {
    replace_child(a_parent, object_node_desc{}, old_child);
  }

  /**
   * This function creates a parent-child connection (edge of the graph).
   * \param a_parent The parent to which a child is added.
   * \param new_child The new child that will become a child of the parent node.
   */
  void create_child(object_node_desc a_parent, object_node_desc new_child) {
    replace_child(a_parent, new_child, object_node_desc{});
  }

  /**
   * This function returns a field-editor linked to a given node in the object-tree.
   * \param a_node The node to which to link the newly created field-editor.
   * \return A field-editor linked to the given node.
   */
  xml_field_editor create_field_editor(object_node_desc a_node) {
    return {this, a_node};
  }

  std::string get_object_name(object_node_desc a_node) const;

  /**
   * This function returns the list of objects in the object graph which are derived from
   * the given type identifier pointer. The vector of strings returned by this function contain
   * object names in the form 'Object_Name (ID:23)' such that they can be used to set object-pointer
   * fields using a xml_field_editor object.
   * \param aType The type identifier of the base-class of which the objects are sought.
   * \return The list of all objects in the object-graph which meet the criteria, with object names compatible with the
   * 'set_field_value' function in xml_field_editor.
   */
  std::vector<std::string> get_objects_derived_from(rtti::so_type* aType) const;
};

std::string get_objtree_name(const object_graph& obj_graph,
                             object_node_desc node_id);

object_node_desc get_objtree_node_id(const object_graph& obj_graph,
                                     const std::string& obj_name);

}  // namespace serialization
}  // namespace ReaK

#endif  // REAK_CORE_SERIALIZATION_OBJTREE_ARCHIVER_H_
