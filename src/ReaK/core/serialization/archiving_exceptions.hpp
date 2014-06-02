/**
 * \file archiving_exceptions.hpp
 *
 * This library declares the exception classes to represent errors that could occur during the 
 * archiving (input or output) processes.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date November 2013
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

#ifndef REAK_ARCHIVING_EXCEPTIONS_HPP
#define REAK_ARCHIVING_EXCEPTIONS_HPP

#include <ReaK/core/base/defs.hpp>

#include <exception>
#include <string>
#include <sstream>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK's Serialization */
namespace serialization {


/**
 * This exception class represents the case when an archiver is trying to deal with a 
 * particular serializable type (which it should find in the RTTI-repository) and 
 * encountered a problem that causes it to consider this particular type to be "unsupported".
 * For example, if an input-archiver is trying to load an object, it must find and 
 * use a suitable factory function for that type, and if that fails, it will report 
 * the error condition through this exception.
 */
class unsupported_type : public std::exception {
  private:
    std::string mMessage;
  public:
    enum problem { not_found_in_repo, could_not_create };
    
    /**
     * Parametrized constructor.
     * \param aIssue The type of issue encountered.
     * \param aTypeID The type-ID pointer of the type that is deemed "unsupported".
     */
    unsupported_type(problem aIssue, unsigned int* aTypeID) {
      std::stringstream ss;
      ss << "Unsupported type! ";
      switch(aIssue) {
        case not_found_in_repo:
          ss << "Could not find the required type in the rtti-repository! ";
          break;
        case could_not_create:
          ss << "Could not create an object of the required type! ";
          break;
      };
      ss << "TypeID:" << std::hex;
      while((*aTypeID) != 0) {
        ss << " " << *aTypeID;
        ++aTypeID;
      };
      mMessage = ss.str();
    };
    
    virtual ~unsupported_type() BOOST_NOEXCEPT_OR_NOTHROW { };
    
    virtual const char* what() const BOOST_NOEXCEPT_OR_NOTHROW {
      return mMessage.c_str();
    };
  
};







};

};

#endif





