/**
 * \file so_type_repo.h
 *
 * This library implements a mergable cross-module singleton that acts as the repository for
 * all registered shared-object type identifiers. This is the core of the RTTI system.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date january 2010
 */

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

#ifndef REAK_CORE_RTTI_SO_TYPE_REPO_H_
#define REAK_CORE_RTTI_SO_TYPE_REPO_H_

#include "ReaK/core/rtti/so_type.h"

namespace ReaK::rtti {

/**
 * This class declares the interface for the repository of shared object types.
 */
class so_type_repo {
 private:
  so_type* type_map_;

  explicit so_type_repo(so_type* type_map);

  so_type_repo* next_;
  so_type_repo* prev_;

 protected:
  bool is_in_ring(so_type_repo* repo);

 public:
  so_type_repo(const so_type_repo&) = delete;
  so_type_repo& operator=(const so_type_repo&) = delete;

  ~so_type_repo();

  /// This function merges a repo with this.
  void merge(so_type_repo* repo);

  /// This function finds a id in the descendants (recusively) of this.
  so_type* find_type(const unsigned int* id) const;

  /// This function finds a type in the descendants (recusively) of this.
  so_type* find_type(so_type* tp) const;

  /// This function adds a type to the repo.
  so_type* add_type(so_type* tp);

  static so_type_repo& get_instance();
};

}  // namespace ReaK::rtti

#endif  // REAK_CORE_RTTI_SO_TYPE_REPO_H_
