
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

#include "ReaK/core/recorders/vector_recorder.h"

namespace ReaK {

namespace recorder {

vector_recorder::vector_recorder() : vec_data(){};

vector_recorder::vector_recorder(std::vector<std::vector<double>>* a_vec_data)
    : vec_data(a_vec_data) {
  set_file_name("");
};

vector_recorder::~vector_recorder() = default;
;

void vector_recorder::write_row() {
  if (vec_data == nullptr) {
    return;
  }
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if ((row_count_ > 0) && (col_count_ > 0)) {
    std::vector<double> v_tmp(col_count_, 0.0);
    for (std::size_t i = 0; i < col_count_; ++i) {
      v_tmp[i] = values_rm_.front();
      values_rm_.pop();
    };
    vec_data->push_back(v_tmp);
    --row_count_;
  };
};

void vector_recorder::write_names(){};

void vector_recorder::set_file_name(const std::string& file_name){};

void vector_recorder::set_vec_data(
    std::vector<std::vector<double>>* a_vec_data) {
  vec_data = a_vec_data;
};

vector_extractor::vector_extractor() : vec_data(), cur_vec_index(0){};

vector_extractor::vector_extractor(
    const std::vector<std::vector<double>>* a_vec_data)
    : vec_data(a_vec_data),
      cur_vec_index(0){

      };

vector_extractor::~vector_extractor() = default;
;

void vector_extractor::add_name(const std::string& s) {
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  names_.push_back(s);
  ++col_count_;
};

bool vector_extractor::read_row() {
  if ((vec_data == nullptr) || (cur_vec_index >= vec_data->size())) {
    return false;
  }
  std::unique_lock<std::mutex> lock_here(access_mutex_);
  if (col_count_ > 0) {
    const std::vector<double>& v_tmp = (*vec_data)[cur_vec_index];
    for (std::size_t i = 0; (i < col_count_) && (i < v_tmp.size()); ++i) {
      values_rm_.push(v_tmp[i]);
    }
    ++cur_vec_index;
  };
  return true;
};

bool vector_extractor::read_names() {
  return true;
};

void vector_extractor::set_file_name(const std::string& file_name){};

void vector_extractor::set_vec_data(
    const std::vector<std::vector<double>>* a_vec_data) {
  vec_data = a_vec_data;
  cur_vec_index = 0;
};
};  // namespace recorder
};  // namespace ReaK
