/**
 * \file exceptions.hpp
 *
 * This library declares exception classes for the ReaK.RPC library. 
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2014
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


#ifndef REAK_RPC_EXCEPTIONS_HPP
#define REAK_RPC_EXCEPTIONS_HPP

#include <stdexcept>

namespace ReaK {

namespace rpc {

namespace detail {
  class remote_function;
};


class publishing_mismatch : public std::exception {
  private:
    std::string msg;
  public:
    
    /**
     * Default constructor.
     * \param aAction The action that was attempted when the mismatch was detected.
     * \param aGiven The remote function that was attempted to be published / unpublished / whatever.
     * \param aExisting The remote function that already existed and conflicted with the given remote function.
     */
    publishing_mismatch(const std::string& aAction, detail::remote_function* aGiven, detail::remote_function* aExisting);
    /**
     * Destructor.
     */
    ~publishing_mismatch() throw() {};
    
    /**
     * Gets the error message.
     * \return c_string of the error message.
     */
    const char* what() const throw() {
      return msg.c_str();
    };
    
};


class communication_error : public std::exception {
  private:
    std::string msg;
  public:
    
    /**
     * Default constructor.
     * \param aMsg The message associated to the communication error.
     */
    communication_error(std::string aMsg = "") : msg(std::move(aMsg)) { };
    /**
     * Destructor.
     */
    ~communication_error() throw() {};
    
    void set_message(std::string aMsg) { msg = std::move(aMsg); };
    
    /**
     * Gets the error message.
     * \return c_string of the error message.
     */
    const char* what() const throw() {
      return msg.c_str();
    };
    
};

class unrecognized_function : public communication_error {
  public:
    
    /**
     * Default constructor.
     * \param aPFunc The remote function that was not recognized.
     */
    unrecognized_function(detail::remote_function* aPFunc);
    /**
     * Destructor.
     */
    ~unrecognized_function() throw() {};
    
};


class marshalling_error : public communication_error {
  public:
    
    /**
     * Default constructor.
     * \param aPFunc The remote function that could not be marshalled correctly.
     * \param aSerialExceptMsg The marshalling / serializing exception message obtained when attempting the operation.
     */
    marshalling_error(detail::remote_function* aPFunc, const std::string& aSerialExceptMsg);
    /**
     * Destructor.
     */
    ~marshalling_error() throw() {};
    
};


};

};


#endif



