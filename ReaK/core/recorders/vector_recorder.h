/**
 * \file vector_recorder.h
 *
 * This library declares a class for data recording to a vector of vectors of values.
 * Here, "data" is meant as columns of floating-point (double) records of data, such as
 * simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2014
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

#ifndef REAK_CORE_RECORDERS_VECTOR_RECORDER_H_
#define REAK_CORE_RECORDERS_VECTOR_RECORDER_H_

#include "ReaK/core/recorders/data_record.h"

#include <vector>

namespace ReaK::recorder {

/**
 * This class handles file IO operations for a raw binary udp-ip stream.
 * The UDP/IP stream is raw in the sense that names are not communicated through the
 * stream (i.e., no meta-data), and the communication is connection-less, i.e., the
 * recorder just spits out the rows of values.
 */
class vector_recorder : public data_recorder {
 protected:
  std::vector<std::vector<double>>* vec_data;

  void write_row() override;
  void write_names() override;
  void set_stream_impl(
      const std::shared_ptr<std::ostream>& stream_ptr) override {}

 public:
  /**
   * Default constructor.
   */
  vector_recorder();

  /**
   * Constructor that opens a file with name file_name.
   */
  explicit vector_recorder(std::vector<std::vector<double>>* a_vec_data);

  /**
   * Destructor, closes the file.
   */
  ~vector_recorder() override;

  void set_vec_data(std::vector<std::vector<double>>* a_vec_data);

  std::vector<std::vector<double>>* get_vec_data() const { return vec_data; }

  void set_file_name(const std::string& file_name) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_recorder::save(A, data_recorder::get_static_object_type()->version());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_recorder::load(A, data_recorder::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(vector_recorder, 0x81100008, 1, "vector_recorder",
                              data_recorder)
};

/**
 * This class handles file IO operations for a raw binary udp-ip stream.
 * The UDP/IP stream is raw in the sense that names are not communicated through the
 * stream (i.e., no meta-data), and the communication is connection-less, i.e., the
 * extractor just reads out the rows of values received, the meta-data is provided by
 * the user before doing any read operations.
 */
class vector_extractor : public data_extractor {
 protected:
  const std::vector<std::vector<double>>* vec_data;
  std::size_t cur_vec_index;

  bool read_row() override;
  bool read_names() override;
  void set_stream_impl(
      const std::shared_ptr<std::istream>& stream_ptr) override {}

 public:
  void add_name(const std::string& s);

  /**
   * Default constructor.
   */
  vector_extractor();

  /**
   * Constructor that opens a file with name file_name.
   */
  explicit vector_extractor(const std::vector<std::vector<double>>* a_vec_data);

  /**
   * Destructor, closes the file.
   */
  ~vector_extractor() override;

  void set_vec_data(const std::vector<std::vector<double>>* a_vec_data);

  const std::vector<std::vector<double>>& get_vec_data() const {
    return *vec_data;
  }

  void set_file_name(const std::string& file_name) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_extractor::save(A,
                         data_extractor::get_static_object_type()->version());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_extractor::load(A,
                         data_extractor::get_static_object_type()->version());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(vector_extractor, 0x81200008, 1,
                              "vector_extractor", data_extractor)
};

}  // namespace ReaK::recorder

#endif  // REAK_CORE_RECORDERS_VECTOR_RECORDER_H_
