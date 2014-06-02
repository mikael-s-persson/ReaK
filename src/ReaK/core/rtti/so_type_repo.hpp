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

#include <map>
#include <vector>
#include <string>

#include "so_type.hpp"

namespace ReaK {

namespace rtti {


/**
 * This class declares the interface for the repository of shared object types.
 */
class so_type_repo : public shared_object_base {
private:
  so_type_repo(const so_type_repo&);
  so_type_repo& operator=(const so_type_repo&);
  
  so_type* mTypeMap;

  so_type_repo(so_type* aTypeMap);
  
  so_type_repo* next;
  so_type_repo* prev;
  
protected:
  
  virtual bool RK_CALL isInRing(so_type_repo* aRepo);
  
public:
  virtual void RK_CALL destroy() { delete this; };
  
  virtual ~so_type_repo();

  ///This function merges a repo with this.
  virtual void RK_CALL merge( so_type_repo* aRepo );

  ///This function finds a TypeID in the descendants (recusively) of this.
  virtual so_type::weak_pointer RK_CALL findType (const unsigned int* aTypeID ) const;
  
  ///This function finds a TypeID in the descendants (recusively) of this.
  virtual so_type::weak_pointer RK_CALL findType (const so_type::shared_pointer& aTypeID ) const;

  ///This function adds a type to the repo.
  virtual so_type::weak_pointer addType(const so_type::shared_pointer& aTypeID);

  static so_type_repo& getInstance();
};

/// Global function to access the shared-object repository.
so_type_repo& getRKSharedObjTypeRepo();


};

};

#endif







