
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

#include "ReaK/core/recorders/ascii_recorder.h"

#include <mutex>

namespace ReaK::recorder {

void ascii_recorder::write_row() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((out_stream_) && (*out_stream_) && (row_count_ > 0) && (col_count_ > 0)) {
    (*out_stream_) << std::endl;
    (*out_stream_) << values_rm_.front();
    values_rm_.pop();
    for (unsigned int i = 1; i < col_count_; ++i) {
      (*out_stream_) << delimiter << values_rm_.front();
      values_rm_.pop();
    }
    --row_count_;
  }
}

void ascii_recorder::write_names() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((!out_stream_) || (!(*out_stream_))) {
    return;
  }
  (*out_stream_) << "%";
  auto it = names_.begin();
  for (; it != names_.end(); ++it) {
    (*out_stream_) << delimiter << (*it);
  }
  out_stream_->flush();
}

void ascii_recorder::set_stream_impl(
    const std::shared_ptr<std::ostream>& stream_ptr) {
  if (col_count_ != 0) {
    *this << close;
    if ((stream_ptr) && (*stream_ptr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex_);
      out_stream_ = stream_ptr;
      out_stream_->setf(std::ios::scientific, std::ios::floatfield);
      out_stream_->precision(11);
      col_count_ = static_cast<unsigned int>(names_.size());
      lock_here.unlock();
      write_names();
    }
  } else {
    if ((stream_ptr) && (*stream_ptr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex_);
      out_stream_ = stream_ptr;
      out_stream_->setf(std::ios::scientific, std::ios::floatfield);
      out_stream_->precision(11);
    }
  }
}

bool ascii_extractor::read_row() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((in_stream) && (*in_stream) && (col_count_ > 0)) {
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
        values_rm_.push(tmp);
        ++i;
      }
      if (curr_end >= temp.size()) {
        break;
      }
      curr_start = curr_end + delimiter.size();

      curr_end = temp.find(delimiter, curr_start);
      temp_name = temp.substr(curr_start, curr_end - curr_start);
    }
    if (i != col_count_) {
      return false;
    }
  }
  return !((in_stream) && !(*in_stream));
}

bool ascii_extractor::read_names() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
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
      names_.push_back(temp_name);
    }
    if (curr_end >= temp.size()) {
      break;
    }
    curr_start = curr_end + delimiter.size();

    curr_end = temp.find(delimiter, curr_start);
    temp_name = temp.substr(curr_start, curr_end - curr_start);
  }
  col_count_ = static_cast<unsigned int>(names_.size());
  return true;
}

void ascii_extractor::set_stream_impl(
    const std::shared_ptr<std::istream>& stream_ptr) {
  if (col_count_ != 0) {
    *this >> close;
  }
  if ((stream_ptr) && (*stream_ptr)) {
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    in_stream = stream_ptr;
    lock_here.unlock();
    read_names();
  }
}

}  // namespace ReaK::recorder
