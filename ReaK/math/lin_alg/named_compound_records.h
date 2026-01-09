/**
 * \file named_compound_records.h
 *
 * This library declares the class for representing compound records obtained by mapping from a
 * recorded row from data-recorder. The data can travel both ways (input or output).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_MATH_LIN_ALG_NAMED_COMPOUND_RECORDS_H_
#define REAK_MATH_LIN_ALG_NAMED_COMPOUND_RECORDS_H_


#include "ReaK/core/recorders/data_record.h"

#include "ReaK/math/lin_alg/arithmetic_tuple.h"
#include "ReaK/math/lin_alg/vect_alg.h"

#include <exception>
#include <map>
#include <string>
#include <vector>

namespace ReaK::recorder {

/**
 * This class represents a compound record from a data-recorder row. This class essentially
 * wraps an instance of named_value_row and implements variable re-mappings and templated
 * get-set functions (via a proxy) such that one can do, for example:
\code{.cpp}
  named_compound_record dr;
  dr.add("position_vector")("px")("py")("pz")
    .default_value(vect<double,3>(1.0, 2.0, 3.0));
  dr.add("velocity_vector")("vx")("vy")("vz")
    .default_value(vect<double,3>(0.1, 0.2, 0.3));

  vect<double,3> p(4.0, 3.0, 2.0);
  vect<double,3> v(0.4, 0.3, 0.2);

  data_recorder& data_out = ...;
  dr.nvr_data = data_out.get_fresh_named_value_row();
  dr["position_vector"] = p;
  dr["velocity_vector"] = v;
  data_out << dr.nvr_data;

  data_extractor& data_in = ...;
  dr.nvr_data = data_in.get_fresh_named_value_row();
  data_in >> dr.nvr_data;
  p = dr["position_vector"].as< vect<double,3> >();
  v = dr["velocity_vector"].as< vect<double,3> >();
\endcode
 */
class named_compound_record {
 private:
  struct compound_def {
    std::vector<std::string> names;
    vect_n<double> default_values;
    compound_def() : names(), default_values(){};
  };
  std::map<std::string, compound_def> compounds;

 public:
  /// This data member stores the named-value-row associated with this compound-record.
  named_value_row nvr_data;

  /**
   * Default constructor.
   */
  explicit named_compound_record(const named_value_row& aNameRow)
      : compounds(), nvr_data(aNameRow) {}

  /**
   * This class is a proxy for the conversions to and from an internal vector representation
   * of the compound data to some desired destination / source type.
   */
  struct compound_vect {
    vect_n<double> data;

    compound_vect() : data() {}
    // needed to overcome the templated operator / constructor:
    compound_vect(const compound_vect& rhs) : data(rhs.data) {}
    compound_vect& operator=(const compound_vect& rhs) {
      data = rhs.data;
      return *this;
    }

    compound_vect(compound_vect&& rhs) : data(std::move(rhs.data)) {}
    compound_vect& operator=(compound_vect&& rhs) {
      data = std::move(rhs.data);
      return *this;
    }

    /**
     * This function converts the internally-recorded value of the compound to the
     * desired destination type (using the from_vect conversion function).
     * \tparam T Any type that can be created from a vector with the from_vect function template (ADL'd).
     * \return The value of the compound object, converted to the specified destination type.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    template <typename T>
    T as() const {
      using ReaK::from_vect;
      return from_vect<T>(data);
    }

    /**
     * This function write to the internally-recorded value of the compound from an
     * object of a source type that is convertible to a vector of values (with to_vect (ADL'd function template)).
     * \tparam T Any type that is convertible to a vector of doubles via a to_vect function template (ADL'd).
     * \param aValue The value of the compound object, to be converted to a vector using to_vect.
     * \return A reference to this proxy object.
     */
    template <typename T>
    compound_vect& operator=(const T& v) {
      using ReaK::to_vect;
      data = to_vect<double>(v);
    }

    /**
     * This constructor creates the proxy object by initializing its internally-recorded value for the compound from an
     * object of a source type that is convertible to a vector of values (with to_vect (ADL'd function template)).
     * \tparam T Any type that is convertible to a vector of doubles via a to_vect function template (ADL'd).
     * \param aValue The value of the compound object, to be converted to a vector using to_vect.
     * \return A reference to this proxy object.
     */
    template <typename T>
    explicit compound_vect(const T& v) : data() {
      *this = v;
    }
  };

  /**
   * This function write to the internally-recorded value of a compound from a proxy
   * object that contains the vector-stored values converted from a source type.
   * \param aName The name of the compound to be written to.
   * \param aValue The value of the compound object as a proxy object.
   * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
   */
  void set(const std::string& aName, const compound_vect& aValue) {
    std::map<std::string, compound_def>::iterator c_it = compounds.find(aName);
    if (c_it == compounds.end())
      throw out_of_bounds();
    try {
      for (std::size_t i = 0;
           (i < c_it->second.names.size()) && (i < aValue.data.size()); ++i)
        nvr_data[c_it->second.names[i]] = aValue.data[i];
    } catch (...) {
      throw;
    }
  }

  /**
   * This function reads off the internally-recorded value of a compound to a proxy
   * object that can be used to convert to the desired destination type.
   * \param aName The name of the compound to be read off.
   * \return The value of the compound object, as a proxy object to be converted to a destination type.
   * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
   */
  compound_vect get(const std::string& aName) const {
    std::map<std::string, compound_def>::const_iterator c_it =
        compounds.find(aName);
    if (c_it == compounds.end())
      throw out_of_bounds();
    compound_vect result;
    result.data = c_it->second.default_values;
    try {
      for (std::size_t i = 0;
           (i < c_it->second.names.size()) && (i < result.data.size()); ++i)
        result.data[i] = nvr_data[c_it->second.names[i]];
    } catch (...) {
      throw;
    }
  }

  /**
   * This class is a proxy to provide read-write access to a compound object,
   * by using a templated-conversion from the internal vector representation.
   */
  struct compound_proxy {
    named_compound_record* parent;
    std::string name;
    compound_proxy(named_compound_record* aParent, const std::string& aName)
        : parent(aParent), name(aName) {}

    /**
     * This function write to the internally-recorded value of the compound from an
     * object of a source type that is convertible to a vector of values (with to_vect (ADL'd function template)).
     * \tparam T Any type that is convertible to a vector of doubles via a to_vect function template (ADL'd).
     * \param aValue The value of the compound object, to be converted to a vector using to_vect.
     * \return A reference to this proxy object.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    template <typename T>
    const compound_proxy& operator=(const T& aValue) const {
      parent->set(name, compound_vect(aValue));
      return *this;
    }

    /**
     * This function converts the internally-recorded value of the compound to the
     * desired destination type (using the from_vect conversion function).
     * \tparam T Any type that can be created from a vector with the from_vect function template (ADL'd).
     * \return The value of the compound object, converted to the specified destination type.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    template <typename T>
    T as() const {
      return parent->get(name).as<T>();
    }

    /**
     * This function write to the internally-recorded value of the compound from a proxy
     * object that contains the vector-stored values converted from a source type.
     * \param aValue The value of the compound object as a proxy object.
     * \return A reference to this proxy object.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    const compound_proxy& operator=(const compound_vect& aValue) const {
      parent->set(name, aValue);
      return *this;
    }

    /**
     * This function reads off the internally-recorded value of the compound to a proxy
     * object that can be used to convert to the desired destination type.
     * \return The value of the compound object, as a proxy object to be converted to a destination type.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    operator compound_vect() const { return parent->get(name); }
  };

  /**
   * This class is a proxy to provide read-only access to a compound object,
   * by using a templated-conversion from the internal vector representation.
   */
  struct compound_const_proxy {
    const named_compound_record* parent;
    std::string name;
    compound_const_proxy(const named_compound_record* aParent,
                         const std::string& aName)
        : parent(aParent), name(aName) {}

    /**
     * This function converts the internally-recorded value of the compound to the
     * desired destination type (using the from_vect conversion function).
     * \tparam T Any type that can be created from a vector with the from_vect function template (ADL'd).
     * \return The value of the compound object, converted to the specified destination type.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    template <typename T>
    T as() const {
      return parent->get(name).as<T>();
    }

    /**
     * This function reads off the internally-recorded value of the compound to a proxy
     * object that can be used to convert to the desired destination type.
     * \return The value of the compound object, as a proxy object to be converted to a destination type.
     * \throw out_of_bounds If the named-value-record does not have the required fields for this compound type.
     */
    operator compound_vect() const { return parent->get(name); }
  };

  /**
   * This function provides read-write access to a compound of a given name.
   * \param aName The name of the compound to gain read-write access to.
   * \return A proxy object that can be used to read-off or write-to (by conversion from an
   *         internal vector representation) the value of the compound object.
   */
  compound_proxy operator[](const std::string& aName) {
    return compound_proxy(this, aName);
  }

  /**
   * This function provides read-only access to a compound of a given name.
   * \param aName The name of the compound to gain read-only access to.
   * \return A proxy object that can be used to read-off (by conversion from an
   *         internal vector representation) the value of the compound object.
   */
  compound_const_proxy operator[](const std::string& aName) const {
    return compound_const_proxy(this, aName);
  }

  // TODO implement the construction code for creating the compounds.

  /**
   * This is a proxy class that allows the creation of compounds (field names and default-value).
   */
  struct compound_builder {
    compound_def& c_def;
    explicit compound_builder(compound_def& aCDef) : c_def(aCDef) {}

    /**
     * This function adds named-value to the current compound.
     * \param aName The name of the new named-value to be added to the current compound.
     * \return A reference to this proxy.
     */
    compound_builder& operator()(const std::string& aName) {
      c_def.names.push_back(aName);
      if (c_def.default_values.size() < c_def.names.size())
        c_def.default_values.resize(c_def.names.size(), 0.0);
      return *this;
    }

    /**
     * This function sets the default-value for the current compound.
     * \tparam T Any type that is convertible to a vector of doubles via a to_vect function template (ADL'd).
     * \param aData The new default value for the current compound.
     * \return A reference to this proxy.
     */
    template <typename T>
    compound_builder& default_value(const T& aData) {
      using ReaK::to_vect;
      c_def.default_values = to_vect<double>(aData);
      if (c_def.default_values.size() < c_def.names.size())
        c_def.default_values.resize(c_def.names.size(), 0.0);
      return *this;
    }
  };

  /**
   * This function creates a new compound definition (or reuses existing one) and
   * returns a proxy that can be used to add named-values (in order) and a
   * default-value to the compound.
   * \param aName The name of the compound to be created (or reused).
   * \return A proxy object that can be used to add named-values (in order) and
   *         a default-value to the compound.
   */
  compound_builder add(const std::string& aName) {
    return compound_builder(compounds[aName]);
  }
};

}  // namespace ReaK::recorder

#endif  // REAK_MATH_LIN_ALG_NAMED_COMPOUND_RECORDS_H_
