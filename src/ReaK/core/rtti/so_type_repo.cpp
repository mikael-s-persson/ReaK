
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

#include <ReaK/core/rtti/so_type_repo.hpp>

#include <iostream>


namespace ReaK {

namespace rtti {

  
so_type_repo& so_type_repo::getInstance() {
  static const unsigned int t = 0;
  static so_type* static_TypeMap = new detail::dummy_so_type(&t);
  static so_type_repo RKSharedObjTypeRepo(static_TypeMap);
  return RKSharedObjTypeRepo;
};

so_type_repo& getRKSharedObjTypeRepo() {
  return so_type_repo::getInstance();
};

so_type_repo::so_type_repo(so_type* aTypeMap) : shared_object_base(), mTypeMap(aTypeMap) { next = this; prev = this; };

so_type_repo::~so_type_repo() {
  if(next != this) {
    prev->next = next;
    next->prev = prev; //take 'this' out of the ring (if not empty, of course).
  };
  
  delete mTypeMap;
};

bool RK_CALL so_type_repo::isInRing(so_type_repo* aRepo) {
  so_type_repo* p = this;
  while(p->next != this) {
    if(p->next == aRepo)
      return true;
    p = p->next;
  };
  return false;
};


///This function merges a repo with this.
void RK_CALL so_type_repo::merge( so_type_repo* aRepo ) {
  if(aRepo == NULL)
    return;
  if(isInRing(aRepo))
    return;
  if(next == this) { //this means that 'this' is a unique instance (not really a ring). 'this' can simply be linked into aRepo.
    prev = aRepo->prev;
    aRepo->prev = this;
    next = aRepo;
    prev->next = this;
  } else if(aRepo->next == aRepo) { //this means that aRepo is a unique instance (not really a ring). aRepo can simply be linked into 'this'.
    aRepo->prev = prev;
    prev = aRepo;
    aRepo->next = this;
    aRepo->prev->next = aRepo;
  } else { //both repositories are distinct rings. Merge-in all the elements of aRepo's ring into 'this's ring.
    so_type_repo* p = aRepo;
    while(p->next != aRepo) {
      so_type_repo* pn = p->next;
      p->next = pn->next;
      pn->next->prev = p;
      pn->next = pn;
      pn->prev = pn; //this takes pn out of the ring.
      merge(pn); //now merge this single instance into 'this's ring.
    };
    //at this point, all that remains is aRepo as a unique instance.
    aRepo->prev = prev;
    prev = aRepo;
    aRepo->next = this;
    aRepo->prev->next = aRepo;
  };
};

///This function finds a TypeID in the descendants (recusively) of this.

so_type::weak_pointer RK_CALL so_type_repo::findType(const unsigned int* aTypeID ) const {
  so_type::weak_pointer result = mTypeMap->findDescendant(aTypeID);
  const so_type_repo* p = this;
  while((p->next != this) && (!result.lock())) {
    p = p->next;
    result = p->mTypeMap->findDescendant(aTypeID);
  };
  return result;
};

so_type::weak_pointer RK_CALL so_type_repo::findType(const so_type::shared_pointer& aTypeID ) const {
  return findType(aTypeID->TypeID_begin());
};

///This function adds a type to the repo.
so_type::weak_pointer so_type_repo::addType(const so_type::shared_pointer& aTypeID) {
  so_type::weak_pointer r = findType(aTypeID);
  if((r.lock()) && (r.lock()->TypeVersion() > aTypeID->TypeVersion()))
    return r;

  return mTypeMap->addDescendant( aTypeID );
};

};

};




