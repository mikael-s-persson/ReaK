
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

#ifndef SSS_EXCEPTIONS_HPP
#define SSS_EXCEPTIONS_HPP

#include <exception>
#include <string>

namespace ReaK {

namespace ctrl { 


class system_incoherency : public std::exception {
  public:
    std::string message; ///< Message string that identifies the singular matrix.

    /**
     * Default constructor.
     * \param aMessage the message corresponding to the incoherency.
     */
    system_incoherency(const std::string& aMessage) : message(std::string("State space system is incoherent, with message '") + aMessage + "'") { };
    /**
     * Destructor.
     */
    ~system_incoherency() throw() {};

    /**
     * Gets the error message.
     * \return c_string of the error message.
     */
    const char* what() const throw() {
      return message.c_str();
    };
    
};
  
  
  

};

};

#endif





