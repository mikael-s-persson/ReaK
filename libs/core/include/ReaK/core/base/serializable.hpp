/**
 * \file serializable.hpp
 * This library declares the serializable interface which means that all objects descending
 * this interface can be serialized (saved) and reconstructed from serial data (loaded).
 *
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date february 2010
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

#ifndef REAK_SERIALIZABLE_HPP
#define REAK_SERIALIZABLE_HPP

#include "defs.hpp"
#include "typed_object.hpp"

#include <ReaK/core/serialization/archiver.hpp>

namespace ReaK {

/**
 * This class is the interface to be implemented in order to make a class serializable, i.e.
 * saved and loaded from a serial archive. The idea is to save and load all data in the
 * exact same order, with care for version number as well. All primitive types, STL arrays
 * (map and vector) and any other object of the ReaK platform or shared_ptr to them are
 * in principle serializable as well, requiring very little effort by the programmer of a
 * new class to implement the save() and load() methods (simply stacking the data members
 * onto the input/output archive is all that is required really, see classes in the control
 * branch for examples).
 */
class serializable : public typed_object {
public:
  virtual ~serializable(){};

  /**
   * This method saves the content of the object to a serial archive of any type.
   *
   * \pre any valid object state.
   * \post the object state has been saved to the archive such that it can be completely reconstructed at load-time.
   * \param A any type of output archive.
   * \param Version the version of this object that is to be saved (always the latest version).
   */
  virtual void RK_CALL save( serialization::oarchive& A, unsigned int Version ) const {
    RK_UNUSED( A );
    RK_UNUSED( Version );
  };

  /**
   * This method loads the content of the object from a serial archive of any type.
   *
   * \pre any object state.
   * \post the object state has been loaded from the archive such that it is completely reconstructed and valid.
   * \param A any type of input archive.
   * \param Version the version of this object that was saved (it is the user's responsability to maintain backward
   *compatibility as much as desired).
   */
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int Version ) {
    RK_UNUSED( A );
    RK_UNUSED( Version );
  };

  RK_RTTI_MAKE_ABSTRACT_1BASE( serializable, 0x80000000, 1, "serializable", typed_object )
};

/**
 * This class can be used to make any class (WrappedClass) appear like a serializable object without 
 * requiring any inheritance on that class. However, it is the user's responsability to assign the 
 * class' meta-data either through the RK_RTTI_REGISTER_CLASS_ID macro or through a manual definition 
 * of the rtti::get_type_id and rtti::get_type_info traits.
 */
template <typename WrappedClass>
class fake_serializable : public serializable {
private:
  WrappedClass* p_obj;
  
  typedef boost::remove_cv<WrappedClass>::type decayed_wrapped_type;
  
public:
  
  explicit fake_serializable(WrappedClass* aPObj) BOOST_NOEXCEPT : p_obj(aPObj) {};
  
  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int aVersion) const {
    serialize(A, *p_obj, aVersion);
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int aVersion) {
    deserialize(A, *p_obj, aVersion);
  };

  RK_RTTI_REGISTER_CLASS_1BASE( decayed_wrapped_type, 1, serializable );
  
};

namespace serialization {

template < typename T >
struct get_fake_serializable_version {
  BOOST_STATIC_CONSTANT(unsigned int, value = 1);
};

template < typename T >
typename boost::disable_if< boost::is_convertible< T&, serializable& >, 
iarchive& >::type operator>>( iarchive& in, T& t ) {
  fake_serializable< T, get_fake_serializable_version<T>::value > fs(&m);
  return in >> fs;
};

template < typename T >
typename boost::disable_if< boost::is_convertible< T&, serializable& >, 
iarchive& >::type operator&( iarchive& in, const std::pair< std::string, T& >& m ) {
  fake_serializable< T, get_fake_serializable_version<T>::value > fs(&(m.second));
  return in & std::pair< std::string, serializable& >(m.first, fs);
};

template < typename T >
typename boost::disable_if< boost::is_convertible< T&, serializable& >, 
oarchive& >::type operator<<( oarchive& out, const T& t ) {
  fake_serializable< const T, get_fake_serializable_version<T>::value > fs(&t);
  return out << fs;
};

template < typename T >
typename boost::disable_if< boost::is_convertible< T&, serializable& >, 
oarchive& >::type operator&( oarchive& out, const std::pair< std::string, const T& >& t ) {
  fake_serializable< const T, get_fake_serializable_version<T>::value > fs(&(t.second));
  return out & std::pair< std::string, const serializable& >(t.first, fs);
};

};

};

#endif
