/**
 * \file archiver.hpp
 *
 * This library declares the base class for an archive to which an object hierarchy
 * can be serialized to. The serialization framework is constructed such that any
 * object type or primitive data type can be saved to and loaded from an archive of any
 * type with simple operators (<< and &).
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

#ifndef ARCHIVER_HPP
#define ARCHIVER_HPP



#include "base/shared_object_base.hpp"

#include "rtti/typed_object.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/preprocessor/stringize.hpp>

#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <list>
#include <iterator>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Serialization */
namespace serialization {

class serializable;


/// This function constructs a name-value-pair for saving purposes.
template <class T>
inline std::pair< std::string, typename ReaK::rtti::get_type_id< T >::save_type > make_save_nvp(const std::string& s, const T& v) {
  return std::pair< std::string, typename ReaK::rtti::get_type_id< T >::save_type >(s,v);
};

/// This function constructs a name-value-pair for loading purposes.
template <class T>
inline std::pair< std::string, typename ReaK::rtti::get_type_id< T >::load_type > make_load_nvp(const std::string& s, T& v) {
  return std::pair< std::string, typename ReaK::rtti::get_type_id< T >::load_type >(s,v);
};

/// Define used for easy association of a value and a name that is the actual variable name, for saving.
#define RK_SERIAL_SAVE_WITH_NAME(VARIABLE) ReaK::serialization::make_save_nvp(std::string(BOOST_PP_STRINGIZE(VARIABLE)),VARIABLE)
/// Define used for easy association of a value and a name that is the actual variable name, for loading.
#define RK_SERIAL_LOAD_WITH_NAME(VARIABLE) ReaK::serialization::make_load_nvp(std::string(BOOST_PP_STRINGIZE(VARIABLE)),VARIABLE)

/// Define used for easy association of a value and a name that is an alias name, for saving.
#define RK_SERIAL_SAVE_WITH_ALIAS(ALIAS_STR,VARIABLE) ReaK::serialization::make_save_nvp(ALIAS_STR,VARIABLE)
/// Define used for easy association of a value and a name that is an alias name, for loading.
#define RK_SERIAL_LOAD_WITH_ALIAS(ALIAS_STR,VARIABLE) ReaK::serialization::make_load_nvp(ALIAS_STR,VARIABLE)



/**
 * This structure holds the information that is put as the header for a serializable object.
 */
struct archive_object_header {
  unsigned int* type_ID; ///< A unique number sequence identifying the object's type.
  unsigned int type_version; ///< The version number of the object's type.
  unsigned int object_ID; ///< A unique number identifying the object instance, uniqueness is guaranteed to within the life-span of the archive.
  bool is_external; ///< Flag that identifies whether the object is serialized in an external archive.
  unsigned int size; ///< Size of the object's data, relevant to a binary archive (to know the space to skip if object is unknown, unsupported).
  archive_object_header() :
      type_ID(NULL),
      type_version(0),
      object_ID(0),
      is_external(false),
      size(0) { };
};

/**
 * This class is the basis for all archives. This class is internal and should not be used (use iarchive or oarchive).
 */
class archive : public shared_object_base {
  private:
    archive(const archive&);
    archive& operator=(const archive&);
  protected:
    std::vector< boost::shared_ptr<serializable> >& mObjRegistry; ///< Holds the registry of objects already loaded.
  public:
    virtual void RK_CALL destroy() { delete this; };

    archive() : shared_object_base(), mObjRegistry(*(new std::vector< boost::shared_ptr<serializable> >())) { };
    
    virtual ~archive() { delete &mObjRegistry; };
};

/**
 * This class is the base class (interface) for all input archive from which an object hierarchy can
 * be reconstructed. The derived class is responsible for registering the object IDs of object instances
 * loaded from the archive in order to link repeated occurrances of the object in the archive.
 */
class iarchive : public archive {
  protected:
    /// Loading a serializable object.
    virtual iarchive& RK_CALL load_serializable_ptr(boost::shared_ptr<serializable>& Item) = 0;

    /// Loading a serializable object with a name.
    virtual iarchive& RK_CALL load_serializable_ptr(const std::pair<std::string, boost::shared_ptr<serializable>& >& Item) = 0;

    /// Loading a serializable object by reference.
    virtual iarchive& RK_CALL load_serializable(serializable& Item) = 0;

    /// Loading a serializable object by reference with a name.
    virtual iarchive& RK_CALL load_serializable(const std::pair<std::string, serializable& >& Item) = 0;

    /// Loading an integer value.
    virtual iarchive& RK_CALL load_int(int& i) = 0;

    /// Loading an integer value with a name.
    virtual iarchive& RK_CALL load_int(const std::pair<std::string, int& >& i) = 0;

    /// Loading an unsigned integer value.
    virtual iarchive& RK_CALL load_unsigned_int(unsigned int& u) = 0;

    /// Loading an unsigned integer value with a name.
    virtual iarchive& RK_CALL load_unsigned_int(const std::pair<std::string, unsigned int& >& u) = 0;

    /// Loading a float value.
    virtual iarchive& RK_CALL load_float(float& f) = 0;

    /// Loading a float value with a name.
    virtual iarchive& RK_CALL load_float(const std::pair<std::string, float& >& f) = 0;

    /// Loading a double value.
    virtual iarchive& RK_CALL load_double(double& d) = 0;

    /// Loading a double value with a name.
    virtual iarchive& RK_CALL load_double(const std::pair<std::string, double& >& d) = 0;

    /// Loading a boolean value.
    virtual iarchive& RK_CALL load_bool(bool& b) = 0;

    /// Loading a boolean value with a name.
    virtual iarchive& RK_CALL load_bool(const std::pair<std::string, bool& >& b) = 0;

    /// Loading a string value.
    virtual iarchive& RK_CALL load_string(std::string& s) = 0;

    /// Loading a string value with a name.
    virtual iarchive& RK_CALL load_string(const std::pair<std::string, std::string& >& s) = 0;

  public:
    /// Default constructor.
    iarchive() : archive() { };

    /// Loading a serializable object.
    friend iarchive& RK_CALL operator >>(iarchive& in, boost::shared_ptr<serializable>& Item) {
      return in.load_serializable_ptr(Item);
    };

    /// Loading a serializable object with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, boost::shared_ptr<serializable>& >& Item) {
      return in.load_serializable_ptr(Item);
    };

    /// Loading a serializable object by reference.
    friend iarchive& RK_CALL operator >>(iarchive& in, serializable& Item) {
      return in.load_serializable(Item);
    };

    /// Loading a serializable object by reference with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, serializable& >& Item) {
      return in.load_serializable(Item);
    };

    /// Loading an integer value.
    friend iarchive& RK_CALL operator >>(iarchive& in, int& i) {
      return in.load_int(i);
    };

    /// Loading an integer value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, int& >& i) {
      return in.load_int(i);
    };

    /// Loading an unsigned integer value.
    friend iarchive& RK_CALL operator >>(iarchive& in, unsigned int& u) {
      return in.load_unsigned_int(u);
    };

    /// Loading an unsigned integer value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, unsigned int& >& u) {
      return in.load_unsigned_int(u);
    };
    
    /// Loading an unsigned integer value.
    friend iarchive& RK_CALL operator >>(iarchive& in, long unsigned int& u) {
      unsigned int tmp;
      in.load_unsigned_int(tmp);
      u = tmp;
      return in;
    };

    /// Loading an unsigned integer value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, long unsigned int& >& u) {
      unsigned int tmp;
      std::pair<std::string, unsigned int&> p(u.first,tmp);
      in.load_unsigned_int(p);
      u.second = p.second;
      return in;
    };

    /// Loading a float value.
    friend iarchive& RK_CALL operator >>(iarchive& in, float& f) {
      return in.load_float(f);
    };

    /// Loading a float value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, float& >& f) {
      return in.load_float(f);
    };

    /// Loading a double value.
    friend iarchive& RK_CALL operator >>(iarchive& in, double& d) {
      return in.load_double(d);
    };

    /// Loading a double value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, double& >& d) {
      return in.load_double(d);
    };

    /// Loading a boolean value.
    friend iarchive& RK_CALL operator >>(iarchive& in, bool& b) {
      return in.load_bool(b);
    };

    /// Loading a boolean value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, bool& >& b) {
      return in.load_bool(b);
    };

    /// Loading a string value.
    friend iarchive& RK_CALL operator >>(iarchive& in, std::string& s) {
      return in.load_string(s);
    };

    /// Loading a string value with a name.
    friend iarchive& RK_CALL operator &(iarchive& in, const std::pair<std::string, std::string& >& s) {
      return in.load_string(s);
    };

    /// Loading a serializable object as a templated pointer.
    template <typename T>
    friend iarchive& operator >>(iarchive& in, boost::shared_ptr<T>& Item) {
      boost::shared_ptr<serializable> tmp;
      in >> tmp;
      if(tmp)
        Item = rtti::rk_dynamic_ptr_cast<T>(tmp);
      else
	Item = boost::shared_ptr<T>();
      return in;
    };

    /// Loading a serializable object as a templated pointer with a name.
    template <typename T>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, boost::shared_ptr<T>& >& Item) {
      boost::shared_ptr<serializable> tmp;
      in & std::pair<std::string, boost::shared_ptr<serializable>& >(Item.first, tmp);
      if(tmp)
        Item.second = rtti::rk_dynamic_ptr_cast<T>(tmp);
      else
	Item.second = boost::shared_ptr<T>();
      return in;
    };

    /// Loading a serializable object as a templated weak pointer.
    template <typename T>
    friend iarchive& operator >>(iarchive& in, boost::weak_ptr<T>& Item) {
      boost::shared_ptr<serializable> tmp;
      in >> tmp;
      if(tmp)
        Item = rtti::rk_dynamic_ptr_cast<T>(tmp);
      else
	Item = boost::weak_ptr<T>();
      return in;
    };

    /// Loading a serializable object as a templated weak pointer with a name.
    template <typename T>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, boost::weak_ptr<T>& >& Item) {
      boost::shared_ptr<serializable> tmp;
      in & std::pair<std::string, boost::shared_ptr<serializable>& >(Item.first, tmp);
      if(tmp)
        Item.second = rtti::rk_dynamic_ptr_cast<T>(tmp);
      else
	Item.second = boost::weak_ptr<T>();
      return in;
    };
    
    template <typename T>
    friend iarchive& operator >>(iarchive& in, T& Item) {
      serializable& tmp = Item;
      in >> tmp;
      return in;
    };
    
    template <typename T>
    friend iarchive& operator &(iarchive& in, const std::pair< std::string, T& >& Item) {
      in & std::pair<std::string, serializable&>(Item.first, Item.second);
      return in;
    };

    /// Loading a STL vector of templated entries.
    template <typename T, typename Allocator>
    friend iarchive& operator >>(iarchive& in, std::vector<T,Allocator>& v) {
      unsigned int count;
      in >> count;
      v.resize(count);
      for(unsigned int i=0;i<count;++i)
	in >> v[i];
      return in;
    };

    /// Loading a STL vector of templated entries with a name.
    template <typename T, typename Allocator>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, std::vector<T,Allocator>& >& v) {
      unsigned int count;
      in & RK_SERIAL_LOAD_WITH_ALIAS(v.first + "_count", count);
      v.second.resize(count);
      for(unsigned int i=0;i<count;++i) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i << "]";
	in & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), v.second[i]);
      };
      return in;
    };

    /// Loading a STL list of templated entries.
    template <typename T, typename Allocator>
    friend iarchive& operator >>(iarchive& in, std::list<T,Allocator>& v) {
      unsigned int count;
      in >> count;
      v.resize(count);
      typename std::list<T,Allocator>::iterator it = v.begin();
      for(;it!=v.end();++it)
	in >> (*it);
      return in;
    };

    /// Loading a STL list of templated entries with a name.
    template <typename T, typename Allocator>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, std::list<T,Allocator>& >& v) {
      unsigned int count;
      in & RK_SERIAL_LOAD_WITH_ALIAS(v.first + "_count", count);
      v.second.resize(count);
      typename std::list<T,Allocator>::iterator it = v.second.begin();
      for(unsigned int i=0;it!=v.second.end();++it) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i++ << "]";
	in & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), (*it));
      };
      return in;
    };

    /// Loading a STL map of templated entries.
    template <typename Key, typename T, typename Compare, typename Allocator>
    friend iarchive& operator >>(iarchive& in, std::map<Key,T,Compare,Allocator>& m) {
      unsigned int count;
      in >> count;
      for(unsigned int i=0;i<count;++i) {
	Key value_key;
	in >> value_key;
	in >> m[value_key];
      };
      return in;
    };

    /// Loading a STL map of templated entries with a name.
    template <typename Key, typename T, typename Compare, typename Allocator>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, std::map<Key,T,Compare,Allocator>& >& m) {
      unsigned int count;
      in & RK_SERIAL_LOAD_WITH_ALIAS(m.first + "_count", count);
      for(unsigned int i=0;i<count;++i) {
	std::stringstream key_s_stream;
	key_s_stream << m.first << "_key[" << i << "]";
	Key value_key;
	in & RK_SERIAL_LOAD_WITH_ALIAS(key_s_stream.str(), value_key);
	std::stringstream value_s_stream;
	value_s_stream << m.first << "_value[" << i << "]";
	in & RK_SERIAL_LOAD_WITH_ALIAS(value_s_stream.str(), m.second[value_key]);
      };
      return in;
    };
    
    /// Loading a STL set of templated entries.
    template <typename T, typename Compare, typename Allocator>
    friend iarchive& operator >>(iarchive& in, std::set<T,Compare,Allocator>& v) {
      unsigned int count;
      in >> count;
      v.resize(count);
      typename std::set<T,Compare,Allocator>::iterator it = v.begin();
      for(;it!=v.end();++it)
	in >> (*it);
      return in;
    };

    /// Loading a STL set of templated entries with a name.
    template <typename T, typename Compare, typename Allocator>
    friend iarchive& operator &(iarchive& in, const std::pair<std::string, std::set<T,Compare,Allocator>& >& v) {
      unsigned int count;
      in & RK_SERIAL_LOAD_WITH_ALIAS(v.first + "_count", count);
      v.second.resize(count);
      typename std::set<T,Compare,Allocator>::iterator it = v.second.begin();
      for(unsigned int i=0;it!=v.second.end();++it) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i++ << "]";
	in & RK_SERIAL_LOAD_WITH_ALIAS(s_stream.str(), (*it));
      };
      return in;
    };

};

/**
 * This class is the base class (interface) for all output archive to which an object hierarchy can
 * be serialized. The derived class is responsible for registering the object IDs of object instances
 * saved to the archive in order to link repeated occurrances of the object in the archive.
 */
class oarchive : public archive {
  protected:
    std::map< boost::shared_ptr<serializable>, unsigned int >& mObjRegMap;

    /// Saving a serializable object to an external archive.
    virtual oarchive& RK_CALL saveToNewArchive_impl(const boost::shared_ptr<serializable>& Item, const std::string& FileName) = 0;

    /// Saving a serializable object with a name to an external archive.
    virtual oarchive& RK_CALL saveToNewArchiveNamed_impl(const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item, const std::string& FileName) = 0;

    /// Saving a serializable object.
    virtual oarchive& RK_CALL save_serializable_ptr(const boost::shared_ptr<serializable>& Item) = 0;

    /// Saving a serializable object with a name.
    virtual oarchive& RK_CALL save_serializable_ptr(const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item) = 0;

    /// Saving a serializable object by reference.
    virtual oarchive& RK_CALL save_serializable(const serializable& Item) = 0;

    /// Saving a serializable object by reference with a name.
    virtual oarchive& RK_CALL save_serializable(const std::pair<std::string, const serializable& >& Item) = 0;

    /// Saving an integer value.
    virtual oarchive& RK_CALL save_int(int i) = 0;

    /// Saving an integer value with a name.
    virtual oarchive& RK_CALL save_int(const std::pair<std::string, int >& i) = 0;

    /// Saving an unsigned integer value.
    virtual oarchive& RK_CALL save_unsigned_int(unsigned int u) = 0;

    /// Saving an unsigned integer value with a name.
    virtual oarchive& RK_CALL save_unsigned_int(const std::pair<std::string, unsigned int >& u) = 0;

    /// Saving a float value.
    virtual oarchive& RK_CALL save_float(float f) = 0;

    /// Saving a float value with a name.
    virtual oarchive& RK_CALL save_float(const std::pair<std::string, float >& f) = 0;

    /// Saving a double value.
    virtual oarchive& RK_CALL save_double(double d) = 0;

    /// Saving a double value with a name.
    virtual oarchive& RK_CALL save_double(const std::pair<std::string, double >& d) = 0;

    /// Saving a boolean value.
    virtual oarchive& RK_CALL save_bool(bool b) = 0;

    /// Saving a boolean value with a name.
    virtual oarchive& RK_CALL save_bool(const std::pair<std::string, bool >& b) = 0;

    /// Saving a string value.
    virtual oarchive& RK_CALL save_string(const std::string& s) = 0;

    /// Saving a string value with a name.
    virtual oarchive& RK_CALL save_string(const std::pair<std::string, const std::string& >& s) = 0;
    
  public:
    /// Default constructor.
    oarchive() : archive(), mObjRegMap(*(new std::map< boost::shared_ptr<serializable>, unsigned int >())) { };
    
    virtual ~oarchive() { delete &mObjRegMap; };

    /// Saving a serializable object to an external archive.
    oarchive& RK_CALL saveToNewArchive(boost::shared_ptr<serializable> Item, const std::string& FileName) {
      return saveToNewArchive_impl(Item,FileName);
    };

    /// Saving a serializable object with a name to an external archive.
    oarchive& RK_CALL saveToNewArchiveNamed(const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item, const std::string& FileName) {
      return saveToNewArchiveNamed_impl(Item,FileName);
    };

    /// Saving a serializable object.
    friend oarchive& RK_CALL operator <<(oarchive& out, boost::shared_ptr<serializable> Item) {
      return out.save_serializable_ptr(Item);
    };

    /// Saving a serializable object with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, const boost::shared_ptr<serializable>& >& Item) {
      return out.save_serializable_ptr(Item);
    };

    /// Saving a serializable object by reference.
    friend oarchive& RK_CALL operator <<(oarchive& out, const serializable& Item) {
      return out.save_serializable(Item);
    };

    /// Saving a serializable object by reference with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, const serializable& >& Item) {
      return out.save_serializable(Item);
    };

    /// Saving an integer value.
    friend oarchive& RK_CALL operator <<(oarchive& out, int i) {
      return out.save_int(i);
    };

    /// Saving an integer value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, int >& i) {
      return out.save_int(i);
    };

    /// Saving an unsigned integer value.
    friend oarchive& RK_CALL operator <<(oarchive& out, unsigned int u) {
      return out.save_unsigned_int(u);
    };

    /// Saving an unsigned integer value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, unsigned int >& u) {
      return out.save_unsigned_int(u);
    };
    
    /// Saving an unsigned integer value.
    friend oarchive& RK_CALL operator <<(oarchive& out, const long unsigned int& u) {
      unsigned int tmp = u;
      return out.save_unsigned_int(tmp);
    };

    /// Saving an unsigned integer value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, const long unsigned int& >& u) {
      return out.save_unsigned_int(std::pair<std::string, unsigned int>(u.first,u.second));
    };

    /// Saving a float value.
    friend oarchive& RK_CALL operator <<(oarchive& out, float f) {
      return out.save_float(f);
    };

    /// Saving a float value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, float >& f) {
      return out.save_float(f);
    };

    /// Saving a double value.
    friend oarchive& RK_CALL operator <<(oarchive& out, double d) {
      return out.save_double(d);
    };

    /// Saving a double value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, double >& d) {
      return out.save_double(d);
    };

    /// Saving a boolean value.
    friend oarchive& RK_CALL operator <<(oarchive& out, bool b) {
      return out.save_bool(b);
    };

    /// Saving a boolean value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, bool >& b) {
      return out.save_bool(b);
    };

    /// Saving a string value.
    friend oarchive& RK_CALL operator <<(oarchive& out, const std::string& s) {
      return out.save_string(s);
    };

    /// Saving a string value with a name.
    friend oarchive& RK_CALL operator &(oarchive& out, const std::pair<std::string, const std::string& >& s) {
      return out.save_string(s);
    };

    /// Saving a serializable object as a templated pointer.
    template <typename T>
    friend oarchive& operator <<(oarchive& out, const boost::shared_ptr<T>& Item) {
      boost::shared_ptr<serializable> tmp;
      if(Item)
        tmp = rtti::rk_dynamic_ptr_cast<serializable>(Item);
      return out << tmp;
    };

    /// Saving a serializable object as a templated pointer with a name.
    template <typename T>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const boost::shared_ptr<T>& >& Item) {
      boost::shared_ptr<serializable> tmp;
      if(Item.second)
        tmp = rtti::rk_dynamic_ptr_cast<serializable>(Item.second);
      return out & std::pair<std::string, const boost::shared_ptr<serializable>& >(Item.first, tmp);
    };

    /// Saving a serializable object as a templated weak pointer.
    template <typename T>
    friend oarchive& operator <<(oarchive& out, const boost::weak_ptr<T>& Item) {
      boost::shared_ptr<serializable> tmp;
      if(!Item.expired())
        tmp = rtti::rk_dynamic_ptr_cast<serializable>(Item.lock());
      return out << tmp;
    };

    /// Saving a serializable object as a templated weak pointer with a name.
    template <typename T>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const boost::weak_ptr<T>& >& Item) {
      boost::shared_ptr<serializable> tmp;
      if(!Item.second.expired())
        tmp = rtti::rk_dynamic_ptr_cast<serializable>(Item.second.lock());
      return out & std::pair<std::string, const boost::shared_ptr<serializable>& >(Item.first, tmp);
    };
    
    /// Saving any object for which there are no other matching overloads (assumed to be a serializable object).
    template <typename T>
    friend oarchive& operator <<(oarchive& in, const T& Item) {
      const serializable& tmp = Item;
      in << tmp;
      return in;
    };
    
    /// Saving any object with name for which there are no other matching overloads (assumed to be a serializable object).
    template <typename T>
    friend oarchive& operator &(oarchive& in, const std::pair< std::string, const T& >& Item) {
      in & std::pair<std::string, const serializable&>(Item.first, Item.second);
      return in;
    };

    /// Saving a STL vector of templated entries.
    template <typename T, typename Allocator>
    friend oarchive& operator <<(oarchive& out, const std::vector<T,Allocator>& v) {
      unsigned int count = v.size();
      out << count;
      for(unsigned int i=0;i<count;++i)
	out << v[i];
      return out;
    };

    /// Saving a STL vector of templated entries with a name.
    template <typename T, typename Allocator>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const std::vector<T,Allocator>& >& v) {
      unsigned int count = v.second.size();
      out & RK_SERIAL_SAVE_WITH_ALIAS(v.first + "_count", count);
      for(unsigned int i=0;i<count;++i) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i << "]";
	out & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), v.second[i]);
      };
      return out;
    };

    /// Saving a STL list of templated entries.
    template <typename T, typename Allocator>
    friend oarchive& operator <<(oarchive& out, const std::list<T,Allocator>& v) {
      unsigned int count = v.size();
      out << count;
      typename std::list<T,Allocator>::const_iterator it = v.begin();
      for(;it!=v.end();++it)
	out << (*it);
      return out;
    };

    /// Saving a STL list of templated entries with a name.
    template <typename T, typename Allocator>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const std::list<T,Allocator>& >& v) {
      unsigned int count = v.second.size();
      out & RK_SERIAL_SAVE_WITH_ALIAS(v.first + "_count", count);
      typename std::list<T,Allocator>::const_iterator it = v.second.begin();
      for(unsigned int i=0;it!=v.second.end();++it) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i++ << "]";
	out & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), (*it));
      };
      return out;
    };

    /// Saving a STL map of templated entries.
    template <typename Key,typename T,typename Compare,typename Allocator>
    friend oarchive& operator <<(oarchive& out, const std::map<Key,T,Compare,Allocator>& m) {
      unsigned int count = m.size();
      out << count;
      typename std::map<Key,T,Compare,Allocator>::const_iterator it = m.begin();
      for(;it != m.end();it++)
	out << it->first << it->second;
      return out;
    };

    /// Saving a STL map of templated entries with a name.
    template <typename Key,typename T,typename Compare,typename Allocator>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const std::map<Key,T,Compare,Allocator>& >& m) {
      unsigned int count = m.second.size();
      out & std::pair<std::string, unsigned int >(m.first + "_count", count);
      typename std::map<Key,T,Compare,Allocator>::const_iterator it = m.second.begin();
      for(unsigned int i=0;it != m.second.end();it++,++i) {
	std::stringstream key_s_stream;
	key_s_stream << m.first << "_key[" << i << "]";
	out & RK_SERIAL_SAVE_WITH_ALIAS(key_s_stream.str(),it->first);
	//(*this) & std::pair<std::string, Key >(key_s_stream.str(), it->first);
	std::stringstream value_s_stream;
	value_s_stream << m.first << "_value[" << i << "]";
	out & RK_SERIAL_SAVE_WITH_ALIAS(value_s_stream.str(), it->second);
      };
      return out;
    };

    /// Saving a STL set of templated entries.
    template <typename T, typename Compare, typename Allocator>
    friend oarchive& operator <<(oarchive& out, const std::set<T,Compare,Allocator>& v) {
      unsigned int count = v.size();
      out << count;
      typename std::set<T,Compare,Allocator>::const_iterator it = v.begin();
      for(;it!=v.end();++it)
	out << (*it);
      return out;
    };

    /// Saving a STL set of templated entries with a name.
    template <typename T, typename Compare, typename Allocator>
    friend oarchive& operator &(oarchive& out, const std::pair<std::string, const std::set<T,Compare,Allocator>& >& v) {
      unsigned int count = v.second.size();
      out & RK_SERIAL_SAVE_WITH_ALIAS(v.first + "_count", count);
      typename std::set<T,Compare,Allocator>::const_iterator it = v.second.begin();
      for(unsigned int i=0;it!=v.second.end();++it) {
	std::stringstream s_stream;
	s_stream << v.first << "_q[" << i++ << "]";
	out & RK_SERIAL_SAVE_WITH_ALIAS(s_stream.str(), (*it));
      };
      return out;
    };
};



};

};

#endif





