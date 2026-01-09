
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

void bin_recorder::write_row() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((out_stream_) && (*out_stream_) && (row_count_ > 0) && (col_count_ > 0)) {
    for (unsigned int i = 0; i < col_count_; ++i) {
      double tmp(values_rm_.front());
      out_stream_->write(reinterpret_cast<char*>(&tmp), sizeof(double));
      values_rm_.pop();
    }
    --row_count_;
  }
}

void bin_recorder::write_names() {
  if ((!out_stream_) || (!(*out_stream_))) {
    return;
  }
  unsigned int col_count = col_count_;
  out_stream_->write(reinterpret_cast<char*>(&col_count), sizeof(unsigned int));
  auto it = names_.begin();
  for (; it != names_.end(); ++it) {
    out_stream_->write(it->c_str(), it->size() + 1);
  }
}

void bin_recorder::set_stream_impl(
    const std::shared_ptr<std::ostream>& stream_ptr) {
  if (col_count_ != 0) {
    *this << close;
    if ((stream_ptr) && (*stream_ptr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex_);
      out_stream_ = stream_ptr;
      col_count_ = static_cast<unsigned int>(names_.size());
      write_names();
    }
  } else {
    if ((stream_ptr) && (*stream_ptr)) {
      std::unique_lock<std::mutex> lock_here(access_mutex_);
      out_stream_ = stream_ptr;
    }
  }
}

bool bin_extractor::read_row() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((in_stream) && (*in_stream) && (col_count_ > 0)) {
    for (unsigned int i = 0; i < col_count_; ++i) {
      double tmp = 0;
      in_stream->read(reinterpret_cast<char*>(&tmp), sizeof(double));
      if (!(*in_stream)) {
        return false;
      }
      values_rm_.push(tmp);
    }
  }
  return true;
}

bool bin_extractor::read_names() {
  if ((!in_stream) || (!(*in_stream))) {
    return false;
  }
  unsigned int col_count = 0;
  in_stream->read(reinterpret_cast<char*>(&col_count), sizeof(unsigned int));
  col_count_ = col_count;
  std::array<char, 128> temp = {};
  for (unsigned int i = 0; i < col_count_; ++i) {
    char* temp_ptr = temp.data();
    while ((temp_ptr < temp.data() + 128) && (in_stream->read(temp_ptr, 1)) &&
           (*temp_ptr != '\0')) {
      ++temp_ptr;
    }
    if ((temp_ptr >= temp.data() + 128) || (!(*in_stream))) {
      return false;
    }
    names_.emplace_back(temp.data());
  }
  return true;
}

void bin_extractor::set_stream_impl(
    const std::shared_ptr<std::istream>& stream_ptr) {
  if (col_count_ != 0) {
    *this >> close;
  }
  if ((stream_ptr) && (*stream_ptr)) {
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    in_stream = stream_ptr;
    read_names();
  }
}

}  // namespace ReaK::recorder
