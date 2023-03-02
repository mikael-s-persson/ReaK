/**
 * \file bin_recorder.hpp
 *
 * This library declares the class for data recording to a binary file. Here, "data" is meant as
 * columns of floating-point (double) records of data, such as simulation results for example.
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

#ifndef BIN_RECORDER_HPP
#define BIN_RECORDER_HPP

#include "data_record.hpp"

#include <fstream>

namespace ReaK::recorder {

/**
 * This class handles file IO operations for a binary data record.
 */
class bin_recorder : public data_recorder {
 protected:
  void writeRow() override;
  void writeNames() override;
  void setStreamImpl(const std::shared_ptr<std::ostream>& aStreamPtr) override;

 public:
  /**
   * Default constructor.
   */
  bin_recorder() = default;

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit bin_recorder(const std::string& aFileName) {
    setFileName(aFileName);
  }

  /**
   * Destructor, closes the file.
   */
  ~bin_recorder() override = default;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_recorder::save(A, data_recorder::getStaticObjectType()->TypeVersion());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_recorder::load(A, data_recorder::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(bin_recorder, 0x81100004, 1, "bin_recorder",
                              data_recorder)
};

/**
 * This class handles file IO operations for a binary data extractor.
 */
class bin_extractor : public data_extractor {
 protected:
  bool readRow() override;
  bool readNames() override;
  void setStreamImpl(const std::shared_ptr<std::istream>& aStreamPtr) override;

 public:
  /**
   * Default constructor.
   */
  bin_extractor() = default;

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit bin_extractor(const std::string& aFileName) {
    setFileName(aFileName);
  }

  /**
   * Destructor, closes the file.
   */
  ~bin_extractor() override = default;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_extractor::save(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_extractor::load(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(bin_extractor, 0x81200004, 1, "bin_extractor",
                              data_extractor)
};

}  // namespace ReaK::recorder

#endif
