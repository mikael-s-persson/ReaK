
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

#include "ReaK/core/recorders/bin_recorder.h"

namespace ReaK::recorder {

void bin_recorder::writeRow() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((out_stream) && (*out_stream) && (rowCount > 0) && (colCount > 0)) {
    for (unsigned int i = 0; i < colCount; ++i) {
      double tmp(values_rm.front());
      out_stream->write(reinterpret_cast<char*>(&tmp), sizeof(double));
      values_rm.pop();
    }
    --rowCount;
  }
}

void bin_recorder::writeNames() {
  if ((!out_stream) || (!(*out_stream))) {
    return;
  }
  unsigned int col_count = colCount;
  out_stream->write(reinterpret_cast<char*>(&col_count), sizeof(unsigned int));
  auto it = names.begin();
  for (; it != names.end(); ++it) {
    out_stream->write(it->c_str(), it->size() + 1);
  }
}

void bin_recorder::setStreamImpl(
    const std::shared_ptr<std::ostream>& aStreamPtr) {
  if (colCount != 0) {
    *this << close;
    if ((aStreamPtr) && (*aStreamPtr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex);
      out_stream = aStreamPtr;
      colCount = static_cast<unsigned int>(names.size());
      writeNames();
    }
  } else {
    if ((aStreamPtr) && (*aStreamPtr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex);
      out_stream = aStreamPtr;
    }
  }
}

bool bin_extractor::readRow() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((in_stream) && (*in_stream) && (colCount > 0)) {
    for (unsigned int i = 0; i < colCount; ++i) {
      double tmp = 0;
      in_stream->read(reinterpret_cast<char*>(&tmp), sizeof(double));
      if (!(*in_stream)) {
        return false;
      }
      values_rm.push(tmp);
    }
  }
  return true;
}

bool bin_extractor::readNames() {
  if ((!in_stream) || (!(*in_stream))) {
    return false;
  }
  unsigned int col_count = 0;
  in_stream->read(reinterpret_cast<char*>(&col_count), sizeof(unsigned int));
  colCount = col_count;
  std::array<char, 128> temp = {};
  for (unsigned int i = 0; i < colCount; ++i) {
    char* temp_ptr = temp.data();
    while ((temp_ptr < temp.data() + 128) && (in_stream->read(temp_ptr, 1)) &&
           (*temp_ptr != '\0')) {
      ++temp_ptr;
    }
    if ((temp_ptr >= temp.data() + 128) || (!(*in_stream))) {
      return false;
    }
    names.emplace_back(temp.data());
  }
  return true;
}

void bin_extractor::setStreamImpl(
    const std::shared_ptr<std::istream>& aStreamPtr) {
  if (colCount != 0) {
    *this >> close;
  }
  if ((aStreamPtr) && (*aStreamPtr)) {
    std::unique_lock<std::mutex> lock_here(access_mutex);
    in_stream = aStreamPtr;
    readNames();
  }
}

}  // namespace ReaK::recorder
