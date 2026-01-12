/**
 * \file scheme_builder.h
 *
 * This library declares the class for a creating type schemes that represent the fields contained in the
 * serialization of a given type.
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

#ifndef REAK_CORE_SERIALIZATION_SCHEME_BUILDER_H_
#define REAK_CORE_SERIALIZATION_SCHEME_BUILDER_H_

#include "ReaK/core/serialization/archiver.h"

#include <map>
#include <stack>
#include <string>
#include <utility>

namespace ReaK::serialization {

class serializable_obj_scheme;
class type_scheme;

std::map<std::string, std::shared_ptr<type_scheme>>& get_global_schemes();

/**
 * Protobuf scheme constructor.
 */
class scheme_builder : public oarchive {
 private:
  std::stack<std::shared_ptr<serializable_obj_scheme>> field_stack;
  std::stack<std::pair<std::string, std::string>> value_name_stack;

 protected:
  template <typename T>
  void save_primitive(const std::string& name);

  oarchive& save_to_new_archive_impl(const serializable_shared_pointer& item,
                                     const std::string& file_name) override;

  oarchive& save_to_new_archive_named_impl(
      const std::pair<std::string, const serializable_shared_pointer&>& item,
      const std::string& file_name) override;

  oarchive& save_serializable_ptr(
      const serializable_shared_pointer& item) override;

  oarchive& save_serializable_ptr(
      const std::pair<std::string, const serializable_shared_pointer&>& item)
      override;

  oarchive& save_serializable(const serializable& item) override;

  oarchive& save_serializable(
      const std::pair<std::string, const serializable&>& item) override;

  oarchive& save_char(char i) override;

  oarchive& save_char(const std::pair<std::string, char>& i) override;

  oarchive& save_unsigned_char(unsigned char u) override;

  oarchive& save_unsigned_char(
      const std::pair<std::string, unsigned char>& u) override;

  oarchive& save_int(std::ptrdiff_t i) override;

  oarchive& save_int(const std::pair<std::string, std::ptrdiff_t>& i) override;

  oarchive& save_unsigned_int(std::size_t u) override;

  oarchive& save_unsigned_int(
      const std::pair<std::string, std::size_t>& u) override;

  oarchive& save_float(float f) override;

  oarchive& save_float(const std::pair<std::string, float>& f) override;

  oarchive& save_double(double d) override;

  oarchive& save_double(const std::pair<std::string, double>& d) override;

  oarchive& save_bool(bool b) override;

  oarchive& save_bool(const std::pair<std::string, bool>& b) override;

  oarchive& save_string(const std::string& s) override;

  oarchive& save_string(
      const std::pair<std::string, const std::string&>& s) override;

  void signal_polymorphic_field(const std::string& base_type_name,
                                const std::uint32_t* type_id,
                                const std::string& field_name) override;

  void start_repeated_field(const std::string& type_name) override;

  void start_repeated_field(const std::string& type_name,
                            const std::string& name) override;

  void finish_repeated_field() override;

  void start_repeated_pair(const std::string& type_name_1,
                           const std::string& type_name_2) override;

  void start_repeated_pair(const std::string& type_name_1,
                           const std::string& type_name_2,
                           const std::string& name) override;

  void finish_repeated_pair() override;

 public:
  scheme_builder();
  ~scheme_builder() override;
};

}  // namespace ReaK::serialization

#endif  // REAK_CORE_SERIALIZATION_SCHEME_BUILDER_H_
