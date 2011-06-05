
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

#include "so_type.hpp"

#include <iostream>


namespace ReaK {

namespace rtti {


bool so_type_impl::compare_shared(const boost::shared_ptr< so_type >& t1, const boost::shared_ptr< so_type >& t2) {
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
  
bool so_type_impl::compare_weak(const boost::weak_ptr< so_type >& t1, const boost::weak_ptr< so_type >& t2) {
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
  

boost::shared_ptr<so_type> RK_CALL so_type_impl::addDescendant(const boost::shared_ptr<so_type>& aObj ) {
  std::set< boost::shared_ptr<so_type>, so_type_impl::compare_shared_t>::iterator it = mDescendants.lower_bound(aObj);
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

boost::weak_ptr<so_type> RK_CALL so_type_impl::addAncestor(boost::shared_ptr<so_type>& aThis, const boost::weak_ptr<so_type>& aObj) {
  boost::shared_ptr<so_type> p = aObj.lock();
  if(p) {
    std::set< boost::weak_ptr<so_type>, so_type_impl::compare_weak_t>::iterator it = mAncestors.lower_bound(aObj);
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
boost::weak_ptr<so_type> RK_CALL so_type_impl::findDescendant_impl(const unsigned int* aTypeID ) const {
  if( compare_equal(aTypeID, this->TypeID_begin()) ) 
    return boost::static_pointer_cast<so_type>(mThis);
  
  if(mDescendants.empty())
    return boost::weak_ptr<so_type>();
  
  detail::dummy_so_type d(aTypeID);
  std::set< boost::shared_ptr<so_type>, so_type_impl::compare_shared_t>::const_iterator it = mDescendants.lower_bound(boost::shared_ptr<so_type>(&d,null_deleter()));
  
  if((it != mDescendants.end()) && (compare_equal((*it)->TypeID_begin(),aTypeID))) 
    return *it;
  
  for(it = mDescendants.begin(); it != mDescendants.end(); ++it) {
    boost::weak_ptr<so_type> p = (*it)->findDescendant(aTypeID);
    if(p.lock())
      return p;
  };
  
  return boost::weak_ptr<so_type>();
};

///This function checks if a typeID is parent to this.
boost::weak_ptr<so_type> RK_CALL so_type_impl::findAncestor_impl(const unsigned int* aTypeID ) const {
  if( compare_equal(aTypeID, this->TypeID_begin()) ) 
    return boost::static_pointer_cast<so_type>(mThis);
  
  if(mAncestors.empty())
    return boost::weak_ptr<so_type>();
  
  detail::dummy_so_type d(aTypeID);
  std::set< boost::weak_ptr<so_type>, so_type_impl::compare_weak_t>::const_iterator it = mAncestors.lower_bound(boost::shared_ptr<so_type>(&d,null_deleter()));
  
  if((it != mAncestors.end()) && (it->lock()) && (compare_equal(it->lock()->TypeID_begin(),aTypeID))) 
    return *it;
  
  for(it = mAncestors.begin(); it != mAncestors.end(); ++it) {
    if(it->lock()) {
      boost::weak_ptr<so_type> p = it->lock()->findAncestor(aTypeID);
      if(p.lock())
	return p;
    };
  };
  
  return boost::weak_ptr<so_type>();
};
  
///This function inserts this into a global repo.
void RK_CALL so_type_impl::insertToRepo(const boost::shared_ptr<so_type>& aThis,boost::shared_ptr<so_type>& aRepo) {
  //Update all ancestors if there are any.
  std::set<boost::weak_ptr<so_type>, so_type_impl::compare_weak_t >::iterator it = mAncestors.begin();
  for(;it != mAncestors.end();) {
    boost::shared_ptr<so_type> p( aRepo->findDescendant(it->lock()) );
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
  std::set<boost::shared_ptr<so_type>, so_type_impl::compare_shared_t >::iterator itd = mDescendants.begin();
  for(;itd != mDescendants.end(); ++itd) {
    if(*itd)
      (*itd)->insertToRepo(*itd,aRepo);
  };

};



};

};









