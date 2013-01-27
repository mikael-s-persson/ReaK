/**
 * \file tcp_recorder.hpp
 *
 * This library declares the class for data recording to a binary tcp-ip stream. Here, 
 * "data" is meant as columns of floating-point (double) records of data, such as 
 * simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date January 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_TCP_RECORDER_HPP
#define REAK_TCP_RECORDER_HPP

#include "data_record.hpp"

namespace ReaK {

namespace recorder {
  

class tcp_server_impl;
class tcp_client_impl;


/**
 * This class handles file IO operations for a binary tcp-ip stream.
 */
class tcp_recorder : public data_recorder {
  protected:
    virtual void writeRow();
    virtual void writeNames();

    shared_ptr<tcp_server_impl> pimpl;
  public:
    
    /**
     * Default constructor.
     */
    tcp_recorder();
    
    /**
     * Constructor that opens a file with name aFileName.
     */
    tcp_recorder(const std::string& aFileName);
    
    /**
     * Destructor, closes the file.
     */
    virtual ~tcp_recorder();

    virtual void setFileName(const std::string& aFileName);

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      data_recorder::save(A,data_recorder::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      data_recorder::load(A,data_recorder::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(tcp_recorder,0x81100005,1,"tcp_recorder",data_recorder)
};



/**
 * This class handles file IO operations for a binary data extractor.
 */
class tcp_extractor : public data_extractor {
  protected:
    virtual bool readRow();
    virtual bool readNames();
    virtual bool loadFile(const std::string& aFileName);
    
    shared_ptr<tcp_client_impl> pimpl;
  public:

    /**
     * Default constructor.
     */
    tcp_extractor();

    /**
     * Constructor that opens a file with name aFileName.
     */
    tcp_extractor(const std::string& aFileName);

    /**
     * Destructor, closes the file.
     */
    virtual ~tcp_extractor();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      data_extractor::save(A,data_extractor::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      data_extractor::load(A,data_extractor::getStaticObjectType()->TypeVersion());
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(tcp_extractor,0x81200005,1,"tcp_extractor",data_extractor)
};



};


};


#endif









