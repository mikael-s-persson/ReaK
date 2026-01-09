
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

#include "ReaK/core/serialization/archiver_factory.h"

#include "ReaK/core/serialization/bin_archiver.h"
#include "ReaK/core/serialization/protobuf_archiver.h"
#include "ReaK/core/serialization/xml_archiver.h"

#include <ios>
#include <memory>

/** Main namespace for ReaK */
namespace ReaK::serialization {

std::shared_ptr<iarchive> open_iarchive(const std::string& file_name) {

  std::string file_ext = file_name.substr(file_name.find_last_of('.') + 1);

  if ((file_ext == "rkx") || (file_ext == "xml")) {
    return std::make_shared<xml_iarchive>(file_name);
  }
  if ((file_ext == "rkb") || (file_ext == "bin")) {
    return std::make_shared<bin_iarchive>(file_name);
  }
  if (file_ext == "pbuf") {
    return std::make_shared<protobuf_iarchive>(file_name);
  }
  throw std::ios_base::failure("Sorry, this file-type is not supported!");
};

std::shared_ptr<oarchive> open_oarchive(const std::string& file_name) {

  std::string file_ext = file_name.substr(file_name.find_last_of('.') + 1);

  if ((file_ext == "rkx") || (file_ext == "xml")) {
    return std::make_shared<xml_oarchive>(file_name);
  }
  if ((file_ext == "rkb") || (file_ext == "bin")) {
    return std::make_shared<bin_oarchive>(file_name);
  }
  if (file_ext == "pbuf") {
    return std::make_shared<protobuf_oarchive>(file_name);
  }
  throw std::ios_base::failure("Sorry, this file-type is not supported!");
};

}  // namespace ReaK::serialization
