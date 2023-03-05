/**
 * \file type_schemes.hpp
 *
 * This library declares the class to represent the type scheme of a given type (i.e., what fields it contains when
 *serialized).
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

#include "ReaK/core/serialization/archiver.hpp"

#include "ReaK/core/base/serializable.hpp"
#include "ReaK/core/base/shared_object.hpp"

#include "ReaK/core/rtti/rtti.hpp"

#include <string>
#include <utility>
#include <vector>

namespace ReaK::serialization {

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
  type_scheme(std::string aTypeName, const unsigned int* aTypeIDPtr,
              unsigned int aTypeVersion = 1)
      : shared_object(),
        m_type_name(std::move(aTypeName)),
        m_type_ID_ptr(aTypeIDPtr),
        m_type_version(aTypeVersion) {}

  ~type_scheme() override = default;

  /**
   * This function returns the name of the type represented by this type-scheme.
   * \return The name of the type represented by this type-scheme.
   */
  const std::string& get_type_name() const { return m_type_name; }

  /**
   * This function returns the ID of the type represented by this type-scheme, as a null-terminated list of unsigned
   * integers.
   * \return The ID of the type represented by this type-scheme, as a null-terminated list of unsigned integers.
   */
  const unsigned int* get_type_ID() const { return m_type_ID_ptr; }

  /**
   * This function returns the version of the type represented by this type-scheme, as an unsigned integers.
   * \return The version of the type represented by this type-scheme, as an unsigned integers.
   */
  unsigned int get_type_version() const { return m_type_version; }

  /**
   * This function determines if the type-scheme is a single-field scheme (or not).
   * \return TRUE if the type-scheme is a single-field scheme (e.g., a primitive type).
   */
  virtual bool is_single_field() const { return true; }

  /**
   * This function returns the number of fields in this type-scheme.
   * \return The number of fields in this type-scheme.
   */
  virtual std::size_t get_field_count() const { return 0; }

  /**
   * This function returns the field for a given index in the type-scheme's list of fields.
   * \param i The index of the field in the type-scheme's list of fields.
   * \return The field for a given index in the type-scheme's list of fields.
   */
  virtual std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const = 0;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    ReaK::shared_object::save(
        A, shared_object::getStaticObjectType()->TypeVersion());
    std::vector<unsigned int> ID_vect;
    const unsigned int* tmp = m_type_ID_ptr;
    while (((tmp) != nullptr) && ((*tmp) != 0U)) {
      ID_vect.push_back(*tmp);
      ++tmp;
    }
    ID_vect.push_back(0);
    A& RK_SERIAL_SAVE_WITH_ALIAS("TypeName", m_type_name) &
        RK_SERIAL_SAVE_WITH_ALIAS("TypeID", ID_vect) &
        RK_SERIAL_SAVE_WITH_ALIAS("TypeVersion", m_type_version);
  }

  void load(serialization::iarchive& A, unsigned int /*Version*/) override {
    ReaK::shared_object::load(
        A, shared_object::getStaticObjectType()->TypeVersion());
    std::vector<unsigned int> ID_vect;
    A& RK_SERIAL_LOAD_WITH_ALIAS("TypeName", m_type_name) &
        RK_SERIAL_LOAD_WITH_ALIAS("TypeID", ID_vect) &
        RK_SERIAL_LOAD_WITH_ALIAS("TypeVersion", m_type_version);
    rtti::so_type* tmp_wptr =
        rtti::getRKSharedObjTypeRepo().findType(ID_vect.data());
    if (tmp_wptr == nullptr) {
      m_type_ID_ptr = nullptr;
    } else {
      m_type_ID_ptr = tmp_wptr->TypeID_begin();
    }
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(type_scheme, 0x81300001, 1, "type_scheme",
                              shared_object)
};

/**
 * This derived class is used to represent a primitive type (double, int, char, etc.).
 */
template <typename T>
class primitive_scheme : public type_scheme {
 public:
  using self = primitive_scheme<T>;

  primitive_scheme() : type_scheme("", nullptr) {
    constexpr auto tname = rtti::get_type_id<T>::type_name;
    m_type_name = std::string(tname);
    auto* tmp_ptr = new unsigned int[2];
    tmp_ptr[0] = rtti::get_type_id<T>::ID;
    tmp_ptr[1] = 0;
    m_type_ID_ptr = tmp_ptr;
  }

  ~primitive_scheme() override { delete[] m_type_ID_ptr; };

  std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const override {
    return {"", std::shared_ptr<type_scheme>()};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    type_scheme::save(A, type_scheme::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    type_scheme::load(A, type_scheme::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0x81300002, 1, "primitive_scheme",
                              type_scheme)
};

/**
 * This derived class is used to represent the scheme for a vector of objects.
 */
class vector_type_scheme : public type_scheme {
 protected:
  std::shared_ptr<type_scheme> m_value_type;

 public:
  explicit vector_type_scheme(
      const std::string& aTypeName = "",
      std::shared_ptr<type_scheme> aValueType = std::shared_ptr<type_scheme>())
      : type_scheme(aTypeName, nullptr), m_value_type(std::move(aValueType)) {}

  bool is_single_field() const override { return false; }

  std::size_t get_field_count() const override { return 1; }

  std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const override {
    if (i == 0) {
      return {"q", m_value_type};
    }
    return {"", std::shared_ptr<type_scheme>()};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    type_scheme::save(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_ALIAS("ValueType", m_value_type);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    type_scheme::load(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_ALIAS("ValueType", m_value_type);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(vector_type_scheme, 0x81300003, 1,
                              "vector_type_scheme", type_scheme)
};

/**
 * This derived class is used to represent the scheme for a map of objects.
 */
class map_type_scheme : public type_scheme {
 protected:
  std::shared_ptr<type_scheme> m_key_type;
  std::shared_ptr<type_scheme> m_value_type;

 public:
  explicit map_type_scheme(
      const std::string& aTypeName = "",
      std::shared_ptr<type_scheme> aKeyType = std::shared_ptr<type_scheme>(),
      std::shared_ptr<type_scheme> aValueType = std::shared_ptr<type_scheme>())
      : type_scheme(aTypeName, nullptr),
        m_key_type(std::move(aKeyType)),
        m_value_type(std::move(aValueType)) {}

  bool is_single_field() const override { return false; }

  std::size_t get_field_count() const override { return 2; }

  std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const override {
    if (i == 0) {
      return {"key", m_key_type};
    }
    if (i == 1) {
      return {"value", m_value_type};
    }
    return {"", std::shared_ptr<type_scheme>()};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    type_scheme::save(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_ALIAS("KeyType", m_key_type) &
        RK_SERIAL_SAVE_WITH_ALIAS("ValueType", m_value_type);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    type_scheme::load(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_ALIAS("KeyType", m_key_type) &
        RK_SERIAL_LOAD_WITH_ALIAS("ValueType", m_value_type);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(map_type_scheme, 0x81300004, 1, "map_type_scheme",
                              type_scheme)
};

/**
 * This derived class is used to represent the scheme for a serializable object (anything derived from
 * ReaK::serialization::serializable).
 */
class serializable_obj_scheme : public type_scheme {
 protected:
  std::vector<std::pair<std::string, std::shared_ptr<type_scheme>>> m_fields;

 public:
  explicit serializable_obj_scheme(const std::string& aTypeName = "",
                                   const unsigned int* aTypeIDPtr = nullptr,
                                   unsigned int aTypeVersion = 1)
      : type_scheme(aTypeName, aTypeIDPtr, aTypeVersion) {}

  bool is_single_field() const override { return false; }

  std::size_t get_field_count() const override { return m_fields.size(); }

  std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const override {
    return m_fields[i];
  }

  /**
   * Adds a field to the type-scheme's list of fields.
   * \param aFieldName The name of the field in the type's scheme.
   * \param aTypeScheme The type-scheme defining the type of the new field.
   */
  void add_field(const std::string& aFieldName,
                 const std::shared_ptr<type_scheme>& aTypeScheme) {
    m_fields.emplace_back(aFieldName, aTypeScheme);
  }

  /**
   * Pops the last field off of the type-scheme's list of fields.
   */
  void pop_last_field() { m_fields.pop_back(); }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    type_scheme::save(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_ALIAS("Fields", m_fields);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    type_scheme::load(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_ALIAS("Fields", m_fields);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(serializable_obj_scheme, 0x81300005, 1,
                              "serializable_obj_scheme", type_scheme)
};

/**
 * This derived class is used to represent the scheme for a serializable pointer to
 * an object (anything derived from ReaK::serializable). This special
 * scheme is required in order to deal with the object-ID field which allows for cross-linking
 * the different std::shared_ptr that point to the same object within an object hierarchy.
 */
class serializable_ptr_scheme : public type_scheme {
 protected:
  std::shared_ptr<type_scheme> m_object_ID_scheme;

 public:
  explicit serializable_ptr_scheme(
      const std::string& aTypeName = "",
      const unsigned int* aTypeIDPtr = nullptr, unsigned int aTypeVersion = 1,
      std::shared_ptr<type_scheme> aObjectIDScheme =
          std::shared_ptr<type_scheme>())
      : type_scheme(aTypeName, aTypeIDPtr, aTypeVersion),
        m_object_ID_scheme(std::move(aObjectIDScheme)) {}

  bool is_single_field() const override { return false; }

  std::size_t get_field_count() const override { return 1; }

  std::pair<std::string, std::shared_ptr<type_scheme>> get_field(
      std::size_t i) const override {
    if (i == 0) {
      return {"object_ID", m_object_ID_scheme};
    }
    return {"", std::shared_ptr<type_scheme>()};
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    type_scheme::save(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_ALIAS("ObjIDField", m_object_ID_scheme);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    type_scheme::load(A, type_scheme::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_ALIAS("ObjIDField", m_object_ID_scheme);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(serializable_ptr_scheme, 0x81300006, 1,
                              "serializable_ptr_scheme", type_scheme)
};

}  // namespace ReaK::serialization

#endif
