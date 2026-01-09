/**
 * \file ascii_recorder.h
 *
 * This library declares the class for data recording to an ASCII file (like Matlab's 'save -ascii' function).
 * Here, "data" is meant as columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date July 2014
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

#ifndef REAK_CORE_RECORDERS_ASCII_RECORDER_H_
#define REAK_CORE_RECORDERS_ASCII_RECORDER_H_

#include "ReaK/core/recorders/data_record.h"

#include <iosfwd>
#include <memory>
#include <utility>
#include <string>

namespace ReaK::recorder {

/**
 * This class handles file IO operations for a space-separated-values data record.
 */
class ascii_recorder : public data_recorder {
 protected:
  void write_row() override;
  void write_names() override;
  void set_stream_impl(
      const std::shared_ptr<std::ostream>& stream_ptr) override;

 public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_recorder() : delimiter(" ") {}

  /**
   * Constructor that opens a file with name file_name.
   */
  explicit ascii_recorder(const std::string& file_name,
                          std::string aDelimiter = " ")
      : delimiter(std::move(aDelimiter)) {
    set_file_name(file_name);
  }

  /**
   * Destructor, closes the file.
   */
  ~ascii_recorder() override = default;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_recorder::save(A, data_recorder::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(delimiter);
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_recorder::load(A, data_recorder::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(delimiter);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(ascii_recorder, 0x81100006, 1, "ascii_recorder",
                              data_recorder)
};

/**
 * This class handles file IO operations for a space-separated-values data extractor.
 */
class ascii_extractor : public data_extractor {
 protected:
  bool read_row() override;
  bool read_names() override;
  void set_stream_impl(
      const std::shared_ptr<std::istream>& stream_ptr) override;

 public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_extractor() : delimiter(" ") {}

  /**
   * Constructor that opens a file with name file_name.
   */
  explicit ascii_extractor(const std::string& file_name,
                           std::string aDelimiter = " ")
      : delimiter(std::move(aDelimiter)) {
    set_file_name(file_name);
  }

  /**
   * Destructor, closes the file.
   */
  ~ascii_extractor() override = default;
  ;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_extractor::save(A,
                         data_extractor::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(delimiter);
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_extractor::load(A,
                         data_extractor::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(delimiter);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(ascii_extractor, 0x81200006, 1, "ascii_extractor",
                              data_extractor)
};

}  // namespace ReaK::recorder

#endif  // REAK_CORE_RECORDERS_ASCII_RECORDER_H_
