/**
 * \file archiver_factory.hpp
 *
 * This library declares abstract factory functions that can be used to create the correct 
 * archive object depending on the file-extension of the file-name provided.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date August 2013
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

#ifndef REAK_ARCHIVER_FACTORY_HPP
#define REAK_ARCHIVER_FACTORY_HPP

#include "archiver.hpp"

namespace ReaK {

namespace serialization {

/**
 * This function creates a file-archive from a given file-name by recognizing its file-extension
 * and creating the correct archiver (if available).
 * \param aFileName The complete file name of the file to open for read operations (iarchive).
 * \return A pointer to an archive object that will be capable of correctly reading from the specified file.
 */
shared_ptr< iarchive > open_iarchive(const std::string& aFileName);

/**
 * This function creates a file-archive from a given file-name by recognizing its file-extension
 * and creating the correct archiver (if available).
 * \param aFileName The complete file name of the file to open for write operations (oarchive).
 * \return A pointer to an archive object that will be capable of correctly writing to the specified file.
 */
shared_ptr< oarchive > open_oarchive(const std::string& aFileName);


};

};

#endif





