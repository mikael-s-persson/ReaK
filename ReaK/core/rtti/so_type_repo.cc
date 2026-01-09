
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

#include "ReaK/core/rtti/so_type_repo.h"

#include <iostream>

namespace ReaK::rtti {

so_type_repo& so_type_repo::get_instance() {
  static const unsigned int t = 0;
  static so_type_ptr static_typemap = create_dummy_so_type(&t);
  static so_type_repo rk_shared_obj_type_repo(static_typemap.ptr);
  return rk_shared_obj_type_repo;
}

so_type_repo::so_type_repo(so_type* type_map)
    : type_map_(type_map), next_(this), prev_(this) {}

so_type_repo::~so_type_repo() {
  if (next_ != this) {
    prev_->next_ = next_;
    next_->prev_ =
        prev_;  // take 'this' out of the ring (if not empty, of course).
  }
}

bool so_type_repo::is_in_ring(so_type_repo* repo) {
  so_type_repo* p = this;
  while (p->next_ != this) {
    if (p->next_ == repo) {
      return true;
    }
    p = p->next_;
  }
  return false;
}

/// This function merges a repo with this.
void so_type_repo::merge(so_type_repo* repo) {
  if (repo == nullptr) {
    return;
  }
  if (is_in_ring(repo)) {
    return;
  }
  if (next_ ==
      this) {  // this means that 'this' is a unique instance (not really a ring). 'this' can simply be linked
               // into repo.
    prev_ = repo->prev_;
    repo->prev_ = this;
    next_ = repo;
    prev_->next_ = this;
  } else if (
      repo->next_ ==
      repo) {  // this means that repo is a unique instance (not really a ring). repo can
               // simply be linked into 'this'.
    repo->prev_ = prev_;
    prev_ = repo;
    repo->next_ = this;
    repo->prev_->next_ = repo;
  } else {  // both repositories are distinct rings. Merge-in all the elements of repo's ring into 'this's ring.
    so_type_repo* p = repo;
    while (p->next_ != repo) {
      so_type_repo* pn = p->next_;
      p->next_ = pn->next_;
      pn->next_->prev_ = p;
      pn->next_ = pn;
      pn->prev_ = pn;  // this takes pn out of the ring.
      merge(pn);       // now merge this single instance into 'this's ring.
    }
    // at this point, all that remains is repo as a unique instance.
    repo->prev_ = prev_;
    prev_ = repo;
    repo->next_ = this;
    repo->prev_->next_ = repo;
  }
}

/// This function finds a TypeID in the descendants (recusively) of this.

so_type* so_type_repo::find_type(const unsigned int* id) const {
  so_type* result = type_map_->find_descendant(id);
  const so_type_repo* p = this;
  while ((p->next_ != this) && (result == nullptr)) {
    p = p->next_;
    result = p->type_map_->find_descendant(id);
  }
  return result;
}

so_type* so_type_repo::find_type(so_type* tp) const {
  return find_type(tp->id_begin());
}

/// This function adds a type to the repo.
so_type* so_type_repo::add_type(so_type* tp) {
  so_type* r = find_type(tp);
  if (((r) != nullptr) && (r->version() > tp->version())) {
    return r;
  }

  return type_map_->add_descendant(tp);
}

}  // namespace ReaK::rtti
