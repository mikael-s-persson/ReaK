
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

#include "ReaK/core/rtti/so_type.h"

#include <algorithm>
#include <iostream>
#include <set>
#include <utility>

namespace ReaK::rtti {

static unsigned int* create_dummy_so_type_id(const unsigned int* id) {
  using SizeType = unsigned int;
  SizeType id_length = 0;
  while ((id != nullptr) && (id[id_length] != 0)) {
    ++id_length;
  }
  ++id_length;
  SizeType* new_id = so_type::create_type_id(id_length);
  for (SizeType i = 0; i < id_length - 1; ++i) {
    new_id[i] = id[i];
  }
  new_id[id_length - 1] = 0;
  return new_id;
}

class so_type_impl : public so_type {
 public:
  static bool compare_ptr(so_type_impl* t1, so_type_impl* t2) {
    if (t1 == nullptr) {
      return true;
    }
    if (t2 == nullptr) {
      return false;
    }

    const unsigned int* pid1 = t1->type_id;
    const unsigned int* pid2 = t2->type_id;
    while (((*pid1) != 0U) && ((*pid2) != 0U)) {
      if (*pid1 < *pid2) {
        return true;
      }
      if (*pid1 > *pid2) {
        return false;
      }
      ++pid1;
      ++pid2;
    }
    return (*pid1 == 0) && (*pid2 > 0);
  }

  using compare_ptr_t = bool (*)(so_type_impl*, so_type_impl*);

  unsigned int type_version;
  unsigned int* type_id;
  std::string type_name;
  construct_ptr construct;

  std::set<so_type_impl*, compare_ptr_t> descendants;
  std::set<so_type_impl*, compare_ptr_t> ancestors;

  using iter = std::set<so_type_impl*, compare_ptr_t>::iterator;

  so_type_impl(unsigned int a_version, unsigned int* a_id,
               std::string_view a_name, construct_ptr a_construct)
      : type_version(a_version),
        type_id(a_id),
        type_name(a_name),
        construct(a_construct),
        descendants(&compare_ptr),
        ancestors(&compare_ptr) {}

  ~so_type_impl() { delete[] type_id; }

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type_impl* find_descendant_impl(const unsigned int* id) {
    if (compare_equal(id, this->type_id)) {
      return this;
    }

    if (descendants.empty()) {
      return nullptr;
    }

    unsigned int* d_type_id = create_dummy_so_type_id(id);
    so_type_impl d(0, d_type_id, "Root", nullptr);
    auto it = descendants.lower_bound(&d);

    if ((it != descendants.end()) && (compare_equal((*it)->type_id, id))) {
      return (*it);
    }

    for (it = descendants.begin(); it != descendants.end(); ++it) {
      so_type_impl* p = (*it)->find_descendant_impl(id);
      if (p != nullptr) {
        return p;
      }
    }

    return nullptr;
  }

  /// This function gets the number of direct descendants of this.
  [[nodiscard]] unsigned int get_descendant_count_impl() const {
    return static_cast<unsigned int>(descendants.size());
  }

  /// This function gets a Type record by index in the direct descendants of this.
  [[nodiscard]] so_type_impl* get_descendant_impl(unsigned int aIndex) const {
    if (aIndex >= descendants.size()) {
      return nullptr;
    }

    auto it = descendants.begin();
    std::advance(it, aIndex);
    return *it;
  }

  /// This function checks if a typeID is parent to this.
  so_type_impl* find_ancestor_impl(const unsigned int* id) {
    if (compare_equal(id, this->type_id)) {
      return this;
    }

    if (ancestors.empty()) {
      return nullptr;
    }

    unsigned int* d_type_id = create_dummy_so_type_id(id);
    so_type_impl d(0, d_type_id, "Root", nullptr);
    auto it = ancestors.lower_bound(&d);

    if ((it != ancestors.end()) && ((*it) != nullptr) &&
        (compare_equal((*it)->type_id, id))) {
      return *it;
    }

    for (it = ancestors.begin(); it != ancestors.end(); ++it) {
      if (*it != nullptr) {
        so_type_impl* p = (*it)->find_ancestor_impl(id);
        if (p != nullptr) {
          return p;
        }
      }
    }

    return nullptr;
  }

  so_type_impl* add_descendant_impl(so_type_impl* tp) {
    auto it = descendants.lower_bound(tp);
    if ((it != descendants.end()) &&
        (compare_equal((*it)->type_id, tp->type_id))) {
      if ((*it)->type_version < tp->type_version) {
        descendants.erase(it);
        descendants.insert(tp);
      } else {
        return *it;
      }
    } else {
      descendants.insert(it, tp);
    }
    return tp;
  }

  so_type_impl* add_ancestor_impl(so_type_impl* tp) {
    if (tp != nullptr) {
      auto it = ancestors.lower_bound(tp);
      if ((it != ancestors.end()) && ((*it) != nullptr) &&
          (compare_equal((*it)->type_id, tp->type_id))) {
        if ((*it)->type_version < tp->type_version) {
          ancestors.erase(it);
          ancestors.insert(tp);
        } else {
          return *it;
        }
      } else {
        ancestors.insert(tp);
      }
      tp->add_descendant_impl(this);
    }
    return tp;
  }

  /// This function inserts this into a global repo.
  void insert_to_repo_impl(so_type_impl* repo) {
    // Update all ancestors if there are any.
    for (auto it = ancestors.begin(); it != ancestors.end();) {
      so_type_impl* p = repo->find_descendant_impl((*it)->type_id);
      if ((p != nullptr) && (p != *it)) {
        ancestors.erase(it++);
        ancestors.insert(p);
        p->add_descendant_impl(
            this);  // Register as descendant of that ancestor.
      } else {
        if (p == nullptr) {
          p = (*it);
          if (p != nullptr) {
            p->insert_to_repo_impl(repo);
          }
        }
        ++it;
      }
    }

    // Add to global repo if there is no ancestors to this.
    if ((*(this->type_id) != 0) && (ancestors.empty())) {
      repo->add_descendant_impl(this);
    }

    // insert all the descendants.
    for (auto* desc : descendants) {
      if (desc != nullptr) {
        desc->insert_to_repo_impl(repo);
      }
    }
  }
};

so_type* so_type::create_type_info(unsigned int a_version, unsigned int* a_id,
                                   std::string_view a_name,
                                   construct_ptr a_construct) {
  return new so_type_impl(a_version, a_id, a_name, a_construct);
}

so_type_ptr create_dummy_so_type(const unsigned int* id) {
  using SizeType = unsigned int;
  SizeType id_length = 0;
  while ((id != nullptr) && (id[id_length] != 0)) {
    ++id_length;
  }
  ++id_length;
  SizeType* new_id = so_type::create_type_id(id_length);
  for (SizeType i = 0; i < id_length - 1; ++i) {
    new_id[i] = id[i];
  }
  new_id[id_length - 1] = 0;
  return so_type_ptr{new so_type_impl(0, new_id, "Root", nullptr)};
}

so_type_ptr::~so_type_ptr() {
  delete static_cast<so_type_impl*>(ptr);
}

unsigned int* so_type::create_type_id(unsigned int id_length) {
  using SizeType = unsigned int;
  // this is just to ensure that new / delete are in same TU
  return new SizeType[id_length];
}

bool so_type::compare_equal(const unsigned int* pid1,
                            const unsigned int* pid2) {
  while (((*pid1) != 0U) && ((*pid2) != 0U)) {
    if (*pid1 != *pid2) {
      return false;
    }
    ++pid1;
    ++pid2;
  }
  return *pid1 == *pid2;
}

/// This function adds a Descendant of this.
so_type* so_type::add_descendant(so_type* tp) {
  return static_cast<so_type_impl*>(this)->add_descendant_impl(
      static_cast<so_type_impl*>(tp));
}

so_type* so_type::add_ancestor(so_type* tp) {
  return static_cast<so_type_impl*>(this)->add_ancestor_impl(
      static_cast<so_type_impl*>(tp));
}

/// This function finds a TypeID in the descendants (recusively) of this.
so_type* so_type::find_descendant(const unsigned int* id) {
  return static_cast<so_type_impl*>(this)->find_descendant_impl(id);
}

/// This function finds a TypeID in the descendants (recusively) of this.
so_type* so_type::find_descendant(so_type* tp) {
  return static_cast<so_type_impl*>(this)->find_descendant_impl(tp->id_begin());
}

/// This function gets the number of direct descendants of this.
unsigned int so_type::get_direct_descendant_count() {
  return static_cast<so_type_impl*>(this)->get_descendant_count_impl();
}

/// This function gets a type record by index in the direct descendants of this.
so_type* so_type::get_direct_descendant(unsigned int id) {
  return static_cast<so_type_impl*>(this)->get_descendant_impl(id);
}

/// This function checks if a typeID is parent to this.
so_type* so_type::find_ancestor(const unsigned int* id) {
  return static_cast<so_type_impl*>(this)->find_ancestor_impl(id);
}

/// This function checks if a typeID is parent to this.
so_type* so_type::find_ancestor(so_type* tp) {
  return static_cast<so_type_impl*>(this)->find_ancestor_impl(tp->id_begin());
}

/// This function inserts this into a global repo.
void so_type::insert_to_repo(so_type* repo) {
  static_cast<so_type_impl*>(this)->insert_to_repo_impl(
      static_cast<so_type_impl*>(repo));
}

const unsigned int* so_type::id_begin() const {
  return static_cast<const so_type_impl*>(this)->type_id;
}

unsigned int so_type::version() const {
  return static_cast<const so_type_impl*>(this)->type_version;
}

const std::string& so_type::name() const {
  return static_cast<const so_type_impl*>(this)->type_name;
}

std::shared_ptr<shared_object> so_type::create_object() const {
  if (static_cast<const so_type_impl*>(this)->construct != nullptr) {
    return static_cast<const so_type_impl*>(this)->construct();
  }
  return {};
}

bool so_type::is_concrete() const {
  return (static_cast<const so_type_impl*>(this)->construct != nullptr);
}

}  // namespace ReaK::rtti
