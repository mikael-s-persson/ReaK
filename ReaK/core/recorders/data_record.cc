
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

#include "ReaK/core/recorders/data_record.h"

#include <chrono>

#include <fstream>

namespace ReaK::recorder {

namespace ch = std::chrono;
using hrc = ch::high_resolution_clock;

void data_recorder::close_record_process() {
  while (row_count_ > 0) {
    write_row();
  }
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  col_count_ = 0;
  if (writing_thread_) {
    lock_here.unlock();
    if (writing_thread_->joinable()) {
      writing_thread_->join();
    }
    lock_here.lock();
    writing_thread_.reset();
  }
}

void data_recorder::record_process::operator()() const {
  hrc::time_point last_time = hrc::now();
  while (parent->col_count_ != 0) {
    do {
      parent->write_row();
    } while (parent->row_count_ > parent->max_buffer_size_);

    if (parent->flush_sample_rate_ == 0) {
      std::this_thread::yield();
    } else {
      double time_period = 1.0 / parent->flush_sample_rate_;
      hrc::time_point time_to_reach =
          last_time + ch::duration_cast<hrc::duration>(
                          ch::duration<double, std::ratio<1, 1>>(time_period));
      std::this_thread::sleep_until(time_to_reach);
    }
    last_time = hrc::now();
  }
}

data_recorder::~data_recorder() {
  close_record_process();
}

data_recorder& data_recorder::operator<<(double value) {
  if (col_count_ == 0) {
    throw end_of_record();
  }
  if (current_column_ >= col_count_) {
    throw out_of_bounds();
  }

  std::unique_lock<std::mutex> lock_here(access_mutex_);
  values_rm_.push(value);
  ++current_column_;
  return *this;
}

data_recorder& data_recorder::operator<<(const named_value_row& values) {
  if (col_count_ == 0) {
    return *this;
  }
  for (double value : values.values_) {
    (*this) << value;
  }
  return (*this) << end_value_row;
}

data_recorder& data_recorder::operator<<(const std::string& name) {
  if (col_count_ == 0) {
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    names_.push_back(name);
  }
  return *this;
}

data_recorder& data_recorder::operator<<(flag some_flag) {
  if (some_flag == end_name_row) {
    if (col_count_ != 0) {
      throw improper_flag();
    }
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    col_count_ = static_cast<unsigned int>(names_.size());
    for (std::size_t i = 0; i < col_count_; ++i) {
      named_indices_[names_[i]] = i;
    }
    lock_here.unlock();
    write_names();
    lock_here.lock();
    current_column_ = 0;
    row_count_ = 0;
    writing_thread_ = std::make_shared<std::thread>(record_process(this));
  } else if (some_flag == end_value_row) {
    if (col_count_ == 0) {
      throw improper_flag();
    }
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    for (; current_column_ < col_count_; ++current_column_) {
      values_rm_.push(0.0);
    }
    current_column_ = 0;
    ++row_count_;
  } else if (some_flag == flush) {
    // flush all data right away... normally would be done at the closure or pause...
    // not while doing other things because the function will not return until this is done.
    while (row_count_ > 0) {
      write_row();
    }
  } else if (some_flag == close) {
    // flush and stop thread.
    close_record_process();
  }
  return *this;
}

void data_recorder::set_file_name(const std::string& file_name) {
  auto file_out = std::make_shared<std::ofstream>(file_name.c_str());
  if (file_out->is_open()) {
    set_stream_impl(file_out);
  }
}

void data_recorder::save(serialization::oarchive& A,
                         unsigned int /*Version*/) const {
  shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(static_cast<unsigned int>(col_count_)) &
      RK_SERIAL_SAVE_WITH_NAME(flush_sample_rate_) &
      RK_SERIAL_SAVE_WITH_NAME(max_buffer_size_) &
      RK_SERIAL_SAVE_WITH_NAME(names_);
}

void data_recorder::load(serialization::iarchive& A, unsigned int /*Version*/) {
  close_record_process();
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  unsigned int col_count = 0;
  A& RK_SERIAL_LOAD_WITH_ALIAS("col_count_", col_count) &
      RK_SERIAL_LOAD_WITH_NAME(flush_sample_rate_) &
      RK_SERIAL_LOAD_WITH_NAME(max_buffer_size_) &
      RK_SERIAL_LOAD_WITH_NAME(names_);
  col_count_ = col_count;
  for (std::size_t i = 0; i < col_count_; ++i) {
    named_indices_[names_[i]] = i;
  }
  row_count_ = 0;
  current_column_ = 0;
  values_rm_ = std::queue<double>();
  lock_here.unlock();
  writing_thread_ = std::make_shared<std::thread>(record_process(this));
}

void data_extractor::close_extract_process() {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  while (!values_rm_.empty()) {
    values_rm_.pop();
  }
  col_count_ = 0;
  if (reading_thread_) {
    lock_here.unlock();
    if (reading_thread_->joinable()) {
      reading_thread_->join();
    }
    lock_here.lock();
    reading_thread_.reset();
  }
}

void data_extractor::extract_process::operator()() const {
  hrc::time_point last_time = hrc::now();
  while (parent->col_count_ != 0) {
    do {
      if (!parent->read_row()) {
        return;
      }
    } while (parent->values_rm_.size() <
             parent->min_buffer_size_ * parent->col_count_);
    if (parent->flush_sample_rate_ == 0) {
      std::this_thread::yield();
    } else {
      double time_period = 1.0 / parent->flush_sample_rate_;
      hrc::time_point time_to_reach =
          last_time + ch::duration_cast<hrc::duration>(
                          ch::duration<double, std::ratio<1, 1>>(time_period));
      std::this_thread::sleep_until(time_to_reach);
    }
    last_time = hrc::now();
  }
}

data_extractor::~data_extractor() {
  close_extract_process();
}

data_extractor& data_extractor::operator>>(double& value) {
  if (col_count_ == 0) {
    throw end_of_record();
  }
  if (current_column_ >= col_count_) {
    throw out_of_bounds();
  }
  while (values_rm_.size() < col_count_ - current_column_) {
    // This second size-test may look redundant,
    //  but values could have been extracted between the last test and this one
    //  (and failures have occurred because of this problem!)
    if (!read_row() && (values_rm_.size() < col_count_ - current_column_)) {
      throw end_of_record();
    }
  }
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  value = values_rm_.front();
  values_rm_.pop();
  ++current_column_;
  return *this;
}

data_extractor& data_extractor::operator>>(std::string& name) {
  if (current_name_col_ >= col_count_) {
    throw out_of_bounds();
  }

  std::unique_lock<std::mutex> lock_here(access_mutex_);
  name = names_[current_name_col_++];
  return *this;
}

data_extractor& data_extractor::operator>>(named_value_row& values) {
  if (col_count_ == 0) {
    return *this;
  }
  for (double& value : values.values_) {
    (*this) >> value;
  }
  return (*this) >> end_value_row;
}

data_extractor& data_extractor::operator>>(flag some_flag) {
  if (some_flag == end_value_row) {
    if (col_count_ == 0) {
      throw improper_flag();
    }
    std::unique_lock<std::mutex> lock_here(access_mutex_);
    for (; current_column_ < col_count_; ++current_column_) {
      values_rm_.pop();
    }
    current_column_ = 0;
  } else if (some_flag == advance) {
    // flush all data right away... normally would be done at the closure or pause...
    // not while doing other things because the function will not return until this is done.
    std::size_t last_count = values_rm_.size();
    do {
      last_count = values_rm_.size();
      if (!read_row()) {
        break;
      }
    } while ((values_rm_.size() > last_count) &&
             (values_rm_.size() < min_buffer_size_ * col_count_));
  } else if (some_flag == close) {
    // flush and stop thread.
    close_extract_process();
  }
  return *this;
}

void data_extractor::set_stream_wrapped_call(
    const std::shared_ptr<std::istream>& stream_ptr) {
  set_stream_impl(stream_ptr);
  for (std::size_t i = 0; i < col_count_; ++i) {
    named_indices_[names_[i]] = i;
  }
  current_column_ = 0;
  reading_thread_ = std::make_shared<std::thread>(extract_process(this));
}

void data_extractor::set_file_name(const std::string& file_name) {
  auto file_in = std::make_shared<std::ifstream>(file_name.c_str());
  if (file_in->is_open()) {
    set_stream_wrapped_call(file_in);
  }
}

void data_extractor::save(serialization::oarchive& A,
                          unsigned int /*Version*/) const {
  shared_object::save(A, shared_object::getStaticObjectType()->TypeVersion());
  A& RK_SERIAL_SAVE_WITH_NAME(static_cast<unsigned int>(col_count_)) &
      RK_SERIAL_SAVE_WITH_NAME(flush_sample_rate_) &
      RK_SERIAL_SAVE_WITH_NAME(min_buffer_size_) &
      RK_SERIAL_SAVE_WITH_NAME(names_);
}

void data_extractor::load(serialization::iarchive& A,
                          unsigned int /*Version*/) {
  close_extract_process();
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  shared_object::load(A, shared_object::getStaticObjectType()->TypeVersion());
  unsigned int col_count = 0;
  A& RK_SERIAL_LOAD_WITH_ALIAS("col_count_", col_count) &
      RK_SERIAL_LOAD_WITH_NAME(flush_sample_rate_) &
      RK_SERIAL_LOAD_WITH_NAME(min_buffer_size_) &
      RK_SERIAL_LOAD_WITH_NAME(names_);
  col_count_ = col_count;
  for (std::size_t i = 0; i < names_.size(); ++i) {
    named_indices_[names_[i]] = i;
  }
  current_column_ = 0;
  current_name_col_ = 0;
  values_rm_ = std::queue<double>();
  lock_here.unlock();
  reading_thread_ = std::make_shared<std::thread>(extract_process(this));
}

}  // namespace ReaK::recorder
