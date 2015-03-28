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

namespace ReaK {

namespace recorder {


/**
 * This class handles file IO operations for a space-separated-values data record.
 */
class ascii_recorder : public data_recorder {
protected:
  virtual void writeRow();
  virtual void writeNames();
  virtual void setStreamImpl( const shared_ptr< std::ostream >& aStreamPtr );

public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_recorder() : data_recorder(), delimiter( " " ){};

  /**
   * Constructor that opens a file with name aFileName.
   */
  ascii_recorder( const std::string& aFileName, const std::string& aDelimiter = " " )
      : data_recorder(), delimiter( aDelimiter ) {
    setFileName( aFileName );
  };

  /**
   * Destructor, closes the file.
   */
  virtual ~ascii_recorder(){};

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    data_recorder::save( A, data_recorder::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( delimiter );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    data_recorder::load( A, data_recorder::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( delimiter );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( ascii_recorder, 0x81100006, 1, "ascii_recorder", data_recorder )
};


/**
 * This class handles file IO operations for a space-separated-values data extractor.
 */
class ascii_extractor : public data_extractor {
protected:
  virtual bool readRow();
  virtual bool readNames();
  virtual void setStreamImpl( const shared_ptr< std::istream >& aStreamPtr );

public:
  std::string delimiter;

  /**
   * Default constructor.
   */
  ascii_extractor() : data_extractor(), delimiter( " " ){};

  /**
   * Constructor that opens a file with name aFileName.
   */
  ascii_extractor( const std::string& aFileName, const std::string& aDelimiter = " " )
      : data_extractor(), delimiter( aDelimiter ) {
    setFileName( aFileName );
  };

  /**
   * Destructor, closes the file.
   */
  virtual ~ascii_extractor(){};

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    data_extractor::save( A, data_extractor::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_SAVE_WITH_NAME( delimiter );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    data_extractor::load( A, data_extractor::getStaticObjectType()->TypeVersion() );
    A& RK_SERIAL_LOAD_WITH_NAME( delimiter );
  };

  RK_RTTI_MAKE_CONCRETE_1BASE( ascii_extractor, 0x81200006, 1, "ascii_extractor", data_extractor )
};
};
};


#endif
