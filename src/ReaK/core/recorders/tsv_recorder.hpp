/**
 * \file tsv_recorder.hpp
 *
 * This library declares the class for data recording to a tab-separated-values file. Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
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

#ifndef TSV_RECORDER_HPP
#define TSV_RECORDER_HPP

#include "ssv_recorder.hpp"

namespace ReaK {

namespace recorder {
  
  
/**
 * This class handles file IO operations for a tab-separated-values data record.
 */
class tsv_recorder : public ssv_recorder {
  protected:
    virtual void writeRow();
    virtual void writeNames();

  public:

    /**
     * Default constructor.
     */
    tsv_recorder() : ssv_recorder() { };

    /**
     * Constructor that opens a file with name aFileName.
     */
    tsv_recorder(const std::string& aFileName) : ssv_recorder(aFileName) { };

    /**
     * Destructor, closes the file.
     */
    virtual ~tsv_recorder() { };

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ssv_recorder::save(A,ssv_recorder::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ssv_recorder::load(A,ssv_recorder::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(tsv_recorder,0x81100003,1,"tsv_recorder",ssv_recorder)
};





/**
 * This class handles file IO operations for a tab-separated-values data extractor.
 */
class tsv_extractor : public ssv_extractor {
  protected:
    virtual bool readRow();
  public:

    /**
     * Default constructor.
     */
    tsv_extractor() : ssv_extractor() { };

    /**
     * Constructor that opens a file with name aFileName.
     */
    tsv_extractor(const std::string& aFileName) : ssv_extractor(aFileName) { };

    /**
     * Destructor, closes the file.
     */
    virtual ~tsv_extractor() {  };

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ssv_extractor::save(A,ssv_extractor::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ssv_extractor::load(A,ssv_extractor::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(tsv_extractor,0x81200003,1,"tsv_extractor",ssv_extractor)
};




};

};

#endif





