
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

#include <ReaK/core/recorders/ascii_recorder.hpp>

#include <mutex>

namespace ReaK::recorder {

void ascii_recorder::writeRow() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((out_stream) && (*out_stream) && (rowCount > 0) && (colCount > 0)) {
    (*out_stream) << std::endl;
    (*out_stream) << values_rm.front();
    values_rm.pop();
    for (unsigned int i = 1; i < colCount; ++i) {
      (*out_stream) << delimiter << values_rm.front();
      values_rm.pop();
    }
    --rowCount;
  }
}

void ascii_recorder::writeNames() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((!out_stream) || (!(*out_stream))) {
    return;
  }
  (*out_stream) << "%";
  auto it = names.begin();
  for (; it != names.end(); ++it) {
    (*out_stream) << delimiter << (*it);
  }
  out_stream->flush();
}

void ascii_recorder::setStreamImpl(
    const std::shared_ptr<std::ostream>& aStreamPtr) {
  if (colCount != 0) {
    *this << close;
    if ((aStreamPtr) && (*aStreamPtr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex);
      out_stream = aStreamPtr;
      out_stream->setf(std::ios::scientific, std::ios::floatfield);
      out_stream->precision(11);
      colCount = static_cast<unsigned int>(names.size());
      lock_here.unlock();
      writeNames();
    }
  } else {
    if ((aStreamPtr) && (*aStreamPtr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex);
      out_stream = aStreamPtr;
      out_stream->setf(std::ios::scientific, std::ios::floatfield);
      out_stream->precision(11);
    }
  }
}

bool ascii_extractor::readRow() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((in_stream) && (*in_stream) && (colCount > 0)) {
    std::string temp;
    std::getline(*in_stream, temp, '\n');
    if (!(*in_stream)) {
      return false;
    }

    std::size_t curr_start = 0;
    std::size_t curr_end = temp.find(delimiter, curr_start);
    std::string temp_name = temp.substr(curr_start, curr_end - curr_start);
    std::size_t i = 0;
    while (true) {
      if (!temp_name.empty()) {
        std::stringstream ss(temp_name);
        double tmp = 0;
        ss >> tmp;
        if (!ss) {
          return false;
        }
        values_rm.push(tmp);
        ++i;
      }
      if (curr_end >= temp.size()) {
        break;
      }
      curr_start = curr_end + delimiter.size();

      curr_end = temp.find(delimiter, curr_start);
      temp_name = temp.substr(curr_start, curr_end - curr_start);
    }
    if (i != colCount) {
      return false;
    }
  }
  return !((in_stream) && !(*in_stream));
}

bool ascii_extractor::readNames() {
  std::unique_lock<std::mutex> lock_here(access_mutex);
  if ((!in_stream) || (!(*in_stream))) {
    return false;
  }
  std::string temp;
  std::getline(*in_stream, temp, '\n');
  std::size_t curr_start = 0;
  std::size_t curr_end = temp.find(delimiter, curr_start);
  std::string temp_name = temp.substr(curr_start, curr_end - curr_start);
  while (true) {
    if ((!temp_name.empty()) && (temp_name[0] != '%')) {
      names.push_back(temp_name);
    }
    if (curr_end >= temp.size()) {
      break;
    }
    curr_start = curr_end + delimiter.size();

    curr_end = temp.find(delimiter, curr_start);
    temp_name = temp.substr(curr_start, curr_end - curr_start);
  }
  colCount = static_cast<unsigned int>(names.size());
  return true;
}

void ascii_extractor::setStreamImpl(
    const std::shared_ptr<std::istream>& aStreamPtr) {
  if (colCount != 0) {
    *this >> close;
  }
  if ((aStreamPtr) && (*aStreamPtr)) {
    std::unique_lock<std::mutex> lock_here(access_mutex);
    in_stream = aStreamPtr;
    lock_here.unlock();
    readNames();
  }
}

}  // namespace ReaK::recorder
