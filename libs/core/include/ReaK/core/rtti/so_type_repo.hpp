/**
 * \file so_type_repo.hpp
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

#ifndef REAK_SO_TYPE_REPO_HPP
#define REAK_SO_TYPE_REPO_HPP

#include <ReaK/core/base/defs.hpp>

#include "so_type.hpp"

namespace ReaK::rtti {

/**
 * This class declares the interface for the repository of shared object types.
 */
class so_type_repo {
 private:
  so_type* mTypeMap;

  explicit so_type_repo(so_type* aTypeMap);

  so_type_repo* next;
  so_type_repo* prev;

 protected:
  bool isInRing(so_type_repo* aRepo);

 public:
  so_type_repo(const so_type_repo&) = delete;
  so_type_repo& operator=(const so_type_repo&) = delete;

  ~so_type_repo();

  /// This function merges a repo with this.
  void merge(so_type_repo* aRepo);

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type* findType(const unsigned int* aTypeID) const;

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type* findType(so_type* aTypeID) const;

  /// This function adds a type to the repo.
  so_type* addType(so_type* aTypeID);

  static so_type_repo& getInstance();
};

/// Global function to access the shared-object repository.
so_type_repo& getRKSharedObjTypeRepo();

}  // namespace ReaK::rtti

#endif
