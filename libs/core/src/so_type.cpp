
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
#include <ReaK/core/base/typed_object.hpp>

#include <iostream>
#include <algorithm>
#include <set>


namespace ReaK {

namespace rtti {


static unsigned int* create_dummy_so_type_id( const unsigned int* aTypeID ) {
  typedef unsigned int SizeType;
  SizeType TypeIDLength = 0;
  while( aTypeID && ( aTypeID[TypeIDLength] != 0 ) )
    ++TypeIDLength;
  ++TypeIDLength;
  SizeType* typeID = so_type::createTypeID( TypeIDLength );
  for( SizeType i = 0; i < TypeIDLength - 1; ++i )
    typeID[i] = aTypeID[i];
  typeID[TypeIDLength - 1] = 0;
  return typeID;
};


class so_type_impl : public so_type {
public:
  static bool compare_ptr( so_type_impl* t1, so_type_impl* t2 ) {
    if( !t1 )
      return true;
    if( !t2 )
      return false;

    const unsigned int* pid1 = t1->mTypeID;
    const unsigned int* pid2 = t2->mTypeID;
    while( ( *pid1 ) && ( *pid2 ) ) {
      if( *pid1 < *pid2 )
        return true;
      else if( *pid1 > *pid2 )
        return false;
      ++pid1;
      ++pid2;
    };
    if( ( *pid1 == 0 ) && ( *pid2 > 0 ) )
      return true;
    else
      return false;
  };

  typedef bool ( *compare_ptr_t )( so_type_impl*, so_type_impl* );

  unsigned int mTypeVersion;
  unsigned int* mTypeID;
  std::string mTypeName;
  construct_ptr mConstruct;

  std::set< so_type_impl*, compare_ptr_t > mDescendants;
  std::set< so_type_impl*, compare_ptr_t > mAncestors;

  typedef std::set< so_type_impl*, compare_ptr_t >::iterator iter;

  so_type_impl( unsigned int aTypeVersion, unsigned int* aTypeID, std::string aTypeName, construct_ptr aConstruct )
      : mTypeVersion( aTypeVersion ), mTypeID( aTypeID ), mTypeName( aTypeName ), mConstruct( aConstruct ),
        mDescendants( &compare_ptr ), mAncestors( &compare_ptr ){};

  ~so_type_impl() { delete[] mTypeID; };

  /// This function finds a TypeID in the descendants (recusively) of this.
  so_type_impl* findDescendant_impl( const unsigned int* aTypeID ) {
    if( compare_equal( aTypeID, this->mTypeID ) )
      return this;

    if( mDescendants.empty() )
      return nullptr;

    unsigned int* d_typeID = create_dummy_so_type_id( aTypeID );
    so_type_impl d( 0, d_typeID, "Root", nullptr );
    iter it = mDescendants.lower_bound( &d );

    if( ( it != mDescendants.end() ) && ( compare_equal( ( *it )->mTypeID, aTypeID ) ) )
      return ( *it );

    for( it = mDescendants.begin(); it != mDescendants.end(); ++it ) {
      so_type_impl* p = ( *it )->findDescendant_impl( aTypeID );
      if( p )
        return p;
    };

    return nullptr;
  };

  /// This function gets the number of direct descendants of this.
  unsigned int getDescendantCount_impl() { return static_cast< unsigned int >( mDescendants.size() ); };

  /// This function gets a Type record by index in the direct descendants of this.
  so_type_impl* getDescendant_impl( unsigned int aIndex ) {
    if( aIndex >= mDescendants.size() )
      return nullptr;

    iter it = mDescendants.begin();
    std::advance( it, aIndex );
    return *it;
  };

  /// This function checks if a typeID is parent to this.
  so_type_impl* findAncestor_impl( const unsigned int* aTypeID ) {
    if( compare_equal( aTypeID, this->mTypeID ) )
      return this;

    if( mAncestors.empty() )
      return nullptr;

    unsigned int* d_typeID = create_dummy_so_type_id( aTypeID );
    so_type_impl d( 0, d_typeID, "Root", nullptr );
    iter it = mAncestors.lower_bound( &d );

    if( ( it != mAncestors.end() ) && ( *it ) && ( compare_equal( ( *it )->mTypeID, aTypeID ) ) )
      return *it;

    for( it = mAncestors.begin(); it != mAncestors.end(); ++it ) {
      if( *it ) {
        so_type_impl* p = ( *it )->findAncestor_impl( aTypeID );
        if( p )
          return p;
      };
    };

    return nullptr;
  };

  so_type_impl* addDescendant_impl( so_type_impl* aObj ) {
    iter it = mDescendants.lower_bound( aObj );
    if( ( it != mDescendants.end() ) && ( compare_equal( ( *it )->mTypeID, aObj->mTypeID ) ) ) {
      if( ( *it )->mTypeVersion < aObj->mTypeVersion ) {
        mDescendants.erase( it );
        mDescendants.insert( aObj );
      } else
        return *it;
    } else
      mDescendants.insert( it, aObj );
    return aObj;
  };

  so_type_impl* addAncestor_impl( so_type_impl* aObj ) {
    if( aObj ) {
      iter it = mAncestors.lower_bound( aObj );
      if( ( it != mAncestors.end() ) && ( *it ) && ( compare_equal( ( *it )->mTypeID, aObj->mTypeID ) ) ) {
        if( ( *it )->mTypeVersion < aObj->mTypeVersion ) {
          mAncestors.erase( it );
          mAncestors.insert( aObj );
        } else
          return *it;
      } else
        mAncestors.insert( aObj );
      aObj->addDescendant_impl( this );
    };
    return aObj;
  };

  /// This function inserts this into a global repo.
  void insertToRepo_impl( so_type_impl* aRepo ) {
    // Update all ancestors if there are any.
    for( iter it = mAncestors.begin(); it != mAncestors.end(); ) {
      so_type_impl* p = aRepo->findDescendant_impl( ( *it )->mTypeID );
      if( p && ( p != *it ) ) {
        mAncestors.erase( it++ );
        mAncestors.insert( p );
        p->addDescendant_impl( this ); // Register as descendant of that ancestor.
      } else {
        if( !p ) {
          p = ( *it );
          if( p )
            p->insertToRepo_impl( aRepo );
        };
        ++it;
      };
    };

    // Add to global repo if there is no ancestors to this.
    if( ( *( this->mTypeID ) != 0 ) && ( mAncestors.empty() ) )
      aRepo->addDescendant_impl( this );

    // insert all the descendants.
    for( iter itd = mDescendants.begin(); itd != mDescendants.end(); ++itd ) {
      if( *itd )
        ( *itd )->insertToRepo_impl( aRepo );
    };
  };
};


so_type* so_type::createTypeInfo( unsigned int aTypeVersion, unsigned int* aTypeID, const std::string& aTypeName,
                                  construct_ptr aConstruct ) {
  return new so_type_impl( aTypeVersion, aTypeID, aTypeName, aConstruct );
};

so_type_ptr create_dummy_so_type( const unsigned int* aTypeID ) {
  typedef unsigned int SizeType;
  SizeType TypeIDLength = 0;
  while( aTypeID && ( aTypeID[TypeIDLength] != 0 ) )
    ++TypeIDLength;
  ++TypeIDLength;
  SizeType* typeID = so_type::createTypeID( TypeIDLength );
  for( SizeType i = 0; i < TypeIDLength - 1; ++i )
    typeID[i] = aTypeID[i];
  typeID[TypeIDLength - 1] = 0;
  return so_type_ptr( new so_type_impl( 0, typeID, "Root", nullptr ) );
};


so_type_ptr::~so_type_ptr() { delete static_cast< so_type_impl* >( ptr ); };

unsigned int* so_type::createTypeID( unsigned int aTypeIDSize ) {
  typedef unsigned int SizeType;
  return new SizeType[aTypeIDSize]; // this is just to ensure that new / delete are in same TU
};


bool so_type::compare_equal( const unsigned int* pid1, const unsigned int* pid2 ) {
  while( ( *pid1 ) && ( *pid2 ) ) {
    if( *pid1 != *pid2 )
      return false;
    ++pid1;
    ++pid2;
  };
  if( *pid1 != *pid2 )
    return false;
  return true;
};

/// This function adds a Descendant of this.
so_type* so_type::addDescendant( so_type* aObj ) {
  return static_cast< so_type_impl* >( this )->addDescendant_impl( static_cast< so_type_impl* >( aObj ) );
};

so_type* so_type::addAncestor( so_type* aObj ) {
  return static_cast< so_type_impl* >( this )->addAncestor_impl( static_cast< so_type_impl* >( aObj ) );
};

/// This function finds a TypeID in the descendants (recusively) of this.
so_type* so_type::findDescendant( const unsigned int* aTypeID ) {
  return static_cast< so_type_impl* >( this )->findDescendant_impl( aTypeID );
};

/// This function finds a TypeID in the descendants (recusively) of this.
so_type* so_type::findDescendant( so_type* aTypeID ) {
  return static_cast< so_type_impl* >( this )->findDescendant_impl( aTypeID->TypeID_begin() );
};

/// This function gets the number of direct descendants of this.
unsigned int so_type::getDirectDescendantCount() {
  return static_cast< so_type_impl* >( this )->getDescendantCount_impl();
};

/// This function gets a type record by index in the direct descendants of this.
so_type* so_type::getDirectDescendant( unsigned int aIndex ) {
  return static_cast< so_type_impl* >( this )->getDescendant_impl( aIndex );
};

/// This function checks if a typeID is parent to this.
so_type* so_type::findAncestor( const unsigned int* aTypeID ) {
  return static_cast< so_type_impl* >( this )->findAncestor_impl( aTypeID );
};

/// This function checks if a typeID is parent to this.
so_type* so_type::findAncestor( so_type* aTypeID ) {
  return static_cast< so_type_impl* >( this )->findAncestor_impl( aTypeID->TypeID_begin() );
};

/// This function inserts this into a global repo.
void so_type::insertToRepo( so_type* aRepo ) {
  static_cast< so_type_impl* >( this )->insertToRepo_impl( static_cast< so_type_impl* >( aRepo ) );
};

const unsigned int* so_type::TypeID_begin() const { return static_cast< const so_type_impl* >( this )->mTypeID; };

unsigned int so_type::TypeVersion() const { return static_cast< const so_type_impl* >( this )->mTypeVersion; };

const std::string& so_type::TypeName() const { return static_cast< const so_type_impl* >( this )->mTypeName; };

ReaK::shared_ptr< shared_object > so_type::CreateObject() const {
  if( static_cast< const so_type_impl* >( this )->mConstruct )
    return static_cast< const so_type_impl* >( this )->mConstruct();
  else
    return ReaK::shared_ptr< shared_object >();
};

bool so_type::isConcrete() const { return ( static_cast< const so_type_impl* >( this )->mConstruct != nullptr ); };
};
};
