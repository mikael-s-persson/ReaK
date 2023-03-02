/**
 * \file ascii_recorder.hpp
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

#ifndef ASCII_RECORDER_HPP
#define ASCII_RECORDER_HPP

#include "data_record.hpp"

#include <fstream>
#include <utility>

namespace ReaK::recorder {

/**
 * This class handles file IO operations for a space-separated-values data record.
 */
class ascii_recorder : public data_recorder {
 protected:
  void writeRow() override;
  void writeNames() override;
  void setStreamImpl(const std::shared_ptr<std::ostream>& aStreamPtr) override;

 public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_recorder() : delimiter(" ") {}

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit ascii_recorder(const std::string& aFileName,
                          std::string aDelimiter = " ")
      : delimiter(std::move(aDelimiter)) {
    setFileName(aFileName);
  }

  /**
   * Destructor, closes the file.
   */
  ~ascii_recorder() override = default;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_recorder::save(A, data_recorder::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(delimiter);
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_recorder::load(A, data_recorder::getStaticObjectType()->TypeVersion());
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
  bool readRow() override;
  bool readNames() override;
  void setStreamImpl(const std::shared_ptr<std::istream>& aStreamPtr) override;

 public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_extractor() : delimiter(" ") {}

  /**
   * Constructor that opens a file with name aFileName.
   */
  explicit ascii_extractor(const std::string& aFileName,
                           std::string aDelimiter = " ")
      : delimiter(std::move(aDelimiter)) {
    setFileName(aFileName);
  }

  /**
   * Destructor, closes the file.
   */
  ~ascii_extractor() override = default;
  ;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    data_extractor::save(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(delimiter);
  }
  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    data_extractor::load(A,
                         data_extractor::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(delimiter);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(ascii_extractor, 0x81200006, 1, "ascii_extractor",
                              data_extractor)
};

}  // namespace ReaK::recorder

#endif
