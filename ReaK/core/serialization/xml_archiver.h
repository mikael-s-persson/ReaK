/**
 * \file xml_archiver.h
 *
 * This library declares the class for an XML archive to which an object hierarchy
 * can be serialized to and from.
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

#ifndef REAK_CORE_SERIALIZATION_XML_ARCHIVER_H_
#define REAK_CORE_SERIALIZATION_XML_ARCHIVER_H_

#include "ReaK/core/serialization/archiver.h"

#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ReaK::serialization {

/**
 * XML input archive.
 */
class xml_iarchive : public iarchive {
 private:
  std::shared_ptr<std::istream> file_stream;
  std::string storage;

  char getNextChar();
  std::string readToken();
  void skipToEndToken(const std::string& name);
  static void trimStr(std::string& s);
  bool readNamedValue(const std::string& value_name, std::string& value_str);
  archive_object_header readHeader(const std::string& obj_name,
                                   std::vector<unsigned int>& outTypeID);

 protected:
  iarchive& load_serializable_ptr(serializable_shared_pointer& Item) override;

  iarchive& load_serializable_ptr(
      const std::pair<std::string, serializable_shared_pointer&>& Item)
      override;

  iarchive& load_serializable(serializable& Item) override;

  iarchive& load_serializable(
      const std::pair<std::string, serializable&>& Item) override;

  iarchive& load_char(char& i) override;

  iarchive& load_char(const std::pair<std::string, char&>& i) override;

  iarchive& load_unsigned_char(unsigned char& u) override;

  iarchive& load_unsigned_char(
      const std::pair<std::string, unsigned char&>& u) override;

  iarchive& load_int(std::ptrdiff_t& i) override;

  iarchive& load_int(const std::pair<std::string, std::ptrdiff_t&>& i) override;

  iarchive& load_unsigned_int(std::size_t& u) override;

  iarchive& load_unsigned_int(
      const std::pair<std::string, std::size_t&>& u) override;

  iarchive& load_float(float& f) override;

  iarchive& load_float(const std::pair<std::string, float&>& f) override;

  iarchive& load_double(double& d) override;

  iarchive& load_double(const std::pair<std::string, double&>& d) override;

  iarchive& load_bool(bool& b) override;

  iarchive& load_bool(const std::pair<std::string, bool&>& b) override;

  iarchive& load_string(std::string& s) override;

  iarchive& load_string(const std::pair<std::string, std::string&>& s) override;

 public:
  explicit xml_iarchive(const std::string& FileName);
  explicit xml_iarchive(std::istream& stream);
  ~xml_iarchive() override;
};

/**
 * XML output archive.
 */
class xml_oarchive : public oarchive {
 private:
  std::shared_ptr<std::ostream> file_stream;
  unsigned int tabulation;

  void addTabulations();

 protected:
  oarchive& saveToNewArchive_impl(const serializable_shared_pointer& Item,
                                  const std::string& FileName) override;

  oarchive& saveToNewArchiveNamed_impl(
      const std::pair<std::string, const serializable_shared_pointer&>& Item,
      const std::string& FileName) override;

  oarchive& save_serializable_ptr(
      const serializable_shared_pointer& Item) override;

  oarchive& save_serializable_ptr(
      const std::pair<std::string, const serializable_shared_pointer&>& Item)
      override;

  oarchive& save_serializable(const serializable& Item) override;

  oarchive& save_serializable(
      const std::pair<std::string, const serializable&>& Item) override;

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

 public:
  explicit xml_oarchive(const std::string& FileName);
  explicit xml_oarchive(std::ostream& stream);
  ~xml_oarchive() override;
};

}  // namespace ReaK::serialization

#endif  // REAK_CORE_SERIALIZATION_XML_ARCHIVER_H_
