/**
 * \file udp_recorder.hpp
 *
 * This library declares the class for data recording to a binary udp-ip stream. Here, 
 * "data" is meant as columns of floating-point (double) records of data, such as 
 * simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date January 2014
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

#ifndef REAK_UDP_RECORDER_HPP
#define REAK_UDP_RECORDER_HPP

#include "data_record.hpp"

namespace ReaK {

namespace recorder {
  

class udp_server_impl;
class udp_client_impl;


/**
 * This class handles file IO operations for a binary udp-ip stream.
 */
class udp_recorder : public data_recorder {
  protected:
    virtual void writeRow();
    virtual void writeNames();
    virtual void setStreamImpl(const shared_ptr<std::ostream>& aStreamPtr) { };
    
    shared_ptr<udp_server_impl> pimpl;
  public:
    
    /**
     * Default constructor.
     */
    udp_recorder();
    
    /**
     * Constructor that opens a file with name aFileName.
     */
    udp_recorder(const std::string& aFileName);
    
    /**
     * Destructor, closes the file.
     */
    virtual ~udp_recorder();
    
    virtual void setFileName(const std::string& aFileName);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      data_recorder::save(A,data_recorder::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      data_recorder::load(A,data_recorder::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(udp_recorder,0x81100006,1,"udp_recorder",data_recorder)
};



/**
 * This class handles file IO operations for a binary udp-ip stream.
 */
class udp_extractor : public data_extractor {
  protected:
    virtual bool readRow();
    virtual bool readNames();
    virtual void setStreamImpl(const shared_ptr<std::istream>& aStreamPtr) { };
    
    shared_ptr<udp_client_impl> pimpl;
  public:
    
    /**
     * Default constructor.
     */
    udp_extractor();
    
    /**
     * Constructor that opens a file with name aFileName.
     */
    udp_extractor(const std::string& aFileName);
    
    /**
     * Destructor, closes the file.
     */
    virtual ~udp_extractor();
    
    virtual void setFileName(const std::string& aFilename);
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      data_extractor::save(A,data_extractor::getStaticObjectType()->TypeVersion());
    };
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      data_extractor::load(A,data_extractor::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(udp_extractor,0x81200006,1,"udp_extractor",data_extractor)
};



};


};


#endif









