
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

#include <ReaK/core/rtti/so_type.hpp>
#include <ReaK/core/rtti/typed_object.hpp>

#include <iostream>


namespace ReaK {

namespace rtti {


bool so_type_impl::compare_shared(const so_type::shared_pointer& t1, const so_type::shared_pointer& t2) {
  if(t1) {
    if(t2) {
      const unsigned int* pid1 = t1->TypeID_begin();
      const unsigned int* pid2 = t2->TypeID_begin();
      while((*pid1) && (*pid2)) {
        if(*pid1 < *pid2)
          return true;
        else if(*pid1 > *pid2)
          return false;
        ++pid1; ++pid2;
      };
      if((*pid1 == 0) && (*pid2 > 0))
        return true;
      else 
        return false;
    } else
      return false;
  } else 
    return true;
};
  
bool so_type_impl::compare_weak(const so_type::weak_pointer& t1, const so_type::weak_pointer& t2) {
  return compare_shared(t1.lock(),t2.lock());
};

// bool so_type_impl::compare_equal(const unsigned int* pid1, const unsigned int* pid2) {
//   while((*pid1) && (*pid2)) {
//     if(*pid1 != *pid2)
//       return false;
//     ++pid1; ++pid2;
//   };
//   if(*pid1 != *pid2)
//     return false;
//   return true;
// };
  
  
so_type_impl::so_type_impl() : so_type(), 
                               mDescendants(&compare_shared), 
                               mAncestors(&compare_weak) { };
  

so_type::shared_pointer RK_CALL so_type_impl::addDescendant(const so_type::shared_pointer& aObj ) {
  std::set< so_type::shared_pointer, so_type_impl::compare_shared_t>::iterator it = mDescendants.lower_bound(aObj);
  if((it != mDescendants.end()) && 
     (compare_equal((*it)->TypeID_begin(),aObj->TypeID_begin()))) {
    if((*it)->TypeVersion() < aObj->TypeVersion()) {
      mDescendants.erase(it);
      mDescendants.insert(aObj);
    } else 
      return *it;
  } else
    mDescendants.insert(it,aObj);
  return aObj;
};

so_type::weak_pointer RK_CALL so_type_impl::addAncestor(so_type::shared_pointer& aThis, const so_type::weak_pointer& aObj) {
  so_type::shared_pointer p = aObj.lock();
  if(p) {
    std::set< so_type::weak_pointer, so_type_impl::compare_weak_t>::iterator it = mAncestors.lower_bound(aObj);
    if((it != mAncestors.end()) &&
       ((*it).lock()) &&
       (compare_equal((*it).lock()->TypeID_begin(),p->TypeID_begin()))) {
      if((*it).lock()->TypeVersion() < p->TypeVersion()) {
        mAncestors.erase(it);
        mAncestors.insert(aObj);
      } else 
        return *it;
    } else
      mAncestors.insert(aObj);
    p->addDescendant(aThis);
  };
  return aObj;
};

  
///This function finds a TypeID in the descendants (recusively) of this.
so_type::weak_pointer RK_CALL so_type_impl::findDescendant_impl(const unsigned int* aTypeID ) const {
  if( compare_equal(aTypeID, this->TypeID_begin()) )
    return so_type::weak_pointer(rk_static_ptr_cast<so_type>(mThis));
  
  if(mDescendants.empty())
    return so_type::weak_pointer();
  
  detail::dummy_so_type d(aTypeID);
  std::set< so_type::shared_pointer, so_type_impl::compare_shared_t>::const_iterator it = mDescendants.lower_bound(so_type::shared_pointer(&d,null_deleter()));
  
  if((it != mDescendants.end()) && (compare_equal((*it)->TypeID_begin(),aTypeID))) 
    return *it;
  
  for(it = mDescendants.begin(); it != mDescendants.end(); ++it) {
    so_type::weak_pointer p = (*it)->findDescendant(aTypeID);
    if(p.lock())
      return p;
  };
  
  return so_type::weak_pointer();
};


///This function gets the number of direct descendants of this.
unsigned int RK_CALL so_type_impl::getDescendantCount_impl() const {
  return mDescendants.size();
};

///This function gets a Type record by index in the direct descendants of this.
so_type::shared_pointer RK_CALL so_type_impl::getDescendant_impl(unsigned int aIndex) const {
  if(aIndex >= mDescendants.size())
    return so_type::shared_pointer();
  
  std::set< so_type::shared_pointer, so_type_impl::compare_shared_t>::const_iterator it = mDescendants.begin();
  std::advance(it, aIndex);
  return *it;
};


///This function checks if a typeID is parent to this.
so_type::weak_pointer RK_CALL so_type_impl::findAncestor_impl(const unsigned int* aTypeID ) const {
  if( compare_equal(aTypeID, this->TypeID_begin()) ) 
    return rk_static_ptr_cast<so_type>(mThis);
    
  if(mAncestors.empty())
    return so_type::weak_pointer();
  
  detail::dummy_so_type d(aTypeID);
  std::set< so_type::weak_pointer, so_type_impl::compare_weak_t>::const_iterator it = mAncestors.lower_bound(so_type::shared_pointer(&d,null_deleter()));
  
  if((it != mAncestors.end()) && (it->lock()) && (compare_equal(it->lock()->TypeID_begin(),aTypeID))) 
    return *it;
  
  for(it = mAncestors.begin(); it != mAncestors.end(); ++it) {
    if(it->lock()) {
      so_type::weak_pointer p = it->lock()->findAncestor(aTypeID);
      if(p.lock())
        return p;
    };
  };
  
  return so_type::weak_pointer();
};
  
///This function inserts this into a global repo.
void RK_CALL so_type_impl::insertToRepo(const so_type::shared_pointer& aThis,so_type::shared_pointer& aRepo) {
  //Update all ancestors if there are any.
  std::set<so_type::weak_pointer, so_type_impl::compare_weak_t >::iterator it = mAncestors.begin();
  for(;it != mAncestors.end();) {
    so_type::shared_pointer p( aRepo->findDescendant(it->lock()) );
    if(p) {
      mAncestors.erase(it++);
      mAncestors.insert(p);
      p->addDescendant(aThis); //Register as descendant of that ancestor.
    } else {
      p = it->lock();
      if(p)
        p->insertToRepo(p,aRepo);
      ++it;
    };
  };

  //Add to global repo if there is no ancestors to this.
  if((*(TypeID_begin()) != 0) && (mAncestors.empty()))
    aRepo->addDescendant(aThis);

  //insert all the descendants.
  std::set<so_type::shared_pointer, so_type_impl::compare_shared_t >::iterator itd = mDescendants.begin();
  for(;itd != mDescendants.end(); ++itd) {
    if(*itd)
      (*itd)->insertToRepo(*itd,aRepo);
  };

};



};

};









