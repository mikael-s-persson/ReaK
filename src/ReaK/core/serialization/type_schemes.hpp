/**
 * \file type_schemes.hpp
 *
 * This library declares the class to represent the type scheme of a given type (i.e., what fields it contains when serialized).
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

#ifndef REAK_TYPE_SCHEMES_HPP
#define REAK_TYPE_SCHEMES_HPP

#include "archiver.hpp"

#include "base/serializable.hpp"

#include <fstream>
#include <stack>


namespace ReaK {

namespace serialization {


/**
 * This class is the base class for schemes that represent a given type.
 * \note Type scheme classes are themselves serialization, meaning that a set of schemes can be saved / loaded to 
 * a file (or any stream).
 */
class type_scheme : public shared_object {
  protected:
    std::string m_type_name;
    const unsigned int* m_type_ID_ptr;
    unsigned int m_type_version;
    
  public:
    
    type_scheme(const std::string& aTypeName, 
                const unsigned int* aTypeIDPtr, 
                unsigned int aTypeVersion = 1) :
                shared_object(),
                m_type_name(aTypeName),
                m_type_ID_ptr(aTypeIDPtr),
                m_type_version(aTypeVersion) { };
    
    virtual ~type_scheme() { };
    
    const std::string& get_type_name() const { return m_type_name; };
    const unsigned int* get_type_ID() const { return m_type_ID_ptr; };
    unsigned int get_type_version() const { return m_type_version; };
    
    virtual bool is_single_field() const { return true; };
    
    virtual std::size_t get_field_count() const { return 0; };
    
    virtual std::pair< std::string, shared_ptr< type_scheme > > get_field(std::size_t i) const = 0;
    
    
};

/**
 * This derived class is used to represent a primitive type (double, int, char, etc.).
 */
template <typename T>
class primitive_scheme : public type_scheme {
  public:
    
    primitive_scheme() : type_scheme(rtti::get_type_id<T>::type_name(), NULL) {
      unsigned int* tmp_ptr = new unsigned int[2];
      tmp_ptr[0] = rtti::get_type_id<T>::ID;
      tmp_ptr[1] = 0;
      m_type_ID_ptr = tmp_ptr;
    };
    
    virtual ~primitive_scheme() {
      delete[] m_type_ID_ptr;
    };
    
    virtual std::pair< std::string, shared_ptr< type_scheme > > get_field(std::size_t i) const {
      return std::pair< std::string, shared_ptr< type_scheme > >("", shared_ptr< type_scheme >());
    };
    
    
  
};


/**
 * This derived class is used to represent the scheme for a serializable object (anything derived from ReaK::serialization::serializable).
 */
class serializable_obj_scheme : public type_scheme {
  protected:
    std::vector< std::pair< std::string, shared_ptr< type_scheme > > > m_fields;
    
  public:
    
    serializable_obj_scheme(const std::string& aTypeName, 
                            const unsigned int* aTypeIDPtr,
                            unsigned int aTypeVersion) : 
                            type_scheme(aTypeName, aTypeIDPtr, aTypeVersion),
                            m_fields() { };
    
    virtual bool is_single_field() const { return false; };
    
    virtual std::size_t get_field_count() const { return m_fields.size(); };
    
    virtual std::pair< std::string, shared_ptr< type_scheme > > get_field(std::size_t i) const {
      return m_fields[i];
    };
    
    void add_field(const std::string& aFieldName, const shared_ptr< type_scheme >& aTypeScheme) {
      m_fields.push_back(std::pair< std::string, shared_ptr< type_scheme >(aFieldName, aTypeScheme));
    };
    
};

/**
 * This derived class is used to represent the scheme for a serializable pointer to 
 * an object (anything derived from ReaK::serialization::serializable). This special 
 * scheme is required in order to deal with the object-ID field which allows for cross-linking
 * the different shared_ptr that point to the same object within an object hierarchy.
 */
class serializable_ptr_scheme : public type_scheme {
  protected:
    shared_ptr< type_scheme > m_object_ID_scheme;
    
  public:
    
    serializable_ptr_scheme(const std::string& aTypeName, 
                            const unsigned int* aTypeIDPtr,
                            unsigned int aTypeVersion,
                            const shared_ptr< type_scheme >& aObjectIDScheme) : 
                            type_scheme(aTypeName, aTypeIDPtr, aTypeVersion),
                            m_object_ID_scheme(aObjectIDScheme) { };
    
    virtual bool is_single_field() const { return false; };
    
    virtual std::size_t get_field_count() const { return 1; };
    
    virtual std::pair< std::string, shared_ptr< type_scheme > > get_field(std::size_t i) const {
      if(i == 0)
        return std::pair< std::string, shared_ptr< type_scheme > >("object_ID", m_object_ID_scheme);
      else
        return std::pair< std::string, shared_ptr< type_scheme > >("", shared_ptr< type_scheme >());
    };
    
};
  
  

}; //serialization

}; //ReaK

#endif






