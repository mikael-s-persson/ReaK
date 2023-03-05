/**
 * \file vector_recorder.hpp
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

#ifndef REAK_VECTOR_RECORDER_HPP
#define REAK_VECTOR_RECORDER_HPP

#include "ReaK/core/recorders/data_record.hpp"

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

  void writeRow() override;
  void writeNames() override;
  void setStreamImpl(const std::shared_ptr<std::ostream>& aStreamPtr) override {
  }

 public:
  /**
   * Default constructor.
   */
  vector_recorder();

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit vector_recorder(std::vector<std::vector<double>>* aVecData);

  /**
   * Destructor, closes the file.
   */
  ~vector_recorder() override;

  void setVecData(std::vector<std::vector<double>>* aVecData);

  std::vector<std::vector<double>>* getVecData() const { return vec_data; }

  void setFileName(const std::string& aFileName) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_recorder::save(A, data_recorder::getStaticObjectType()->TypeVersion());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_recorder::load(A, data_recorder::getStaticObjectType()->TypeVersion());
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

  bool readRow() override;
  bool readNames() override;
  void setStreamImpl(const std::shared_ptr<std::istream>& aStreamPtr) override {
  }

 public:
  void addName(const std::string& s);

  /**
   * Default constructor.
   */
  vector_extractor();

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit vector_extractor(const std::vector<std::vector<double>>* aVecData);

  /**
   * Destructor, closes the file.
   */
  ~vector_extractor() override;

  void setVecData(const std::vector<std::vector<double>>* aVecData);

  const std::vector<std::vector<double>>& getVecData() const {
    return *vec_data;
  }

  void setFileName(const std::string& aFilename) override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_extractor::save(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_extractor::load(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(vector_extractor, 0x81200008, 1,
                              "vector_extractor", data_extractor)
};

}  // namespace ReaK::recorder

#endif  // REAK_VECTOR_RECORDER_HPP
