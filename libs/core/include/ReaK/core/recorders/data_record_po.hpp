/**
 * \file data_record_po.hpp
 *
 * This library declares utility functions for creating and dealing with program-options related
 * to a data recording or extraction stream (see data_record.hpp). Here, "data" is meant as
 * columns of floating-point records of data, such as simulation results for example.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2014
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

#ifndef REAK_DATA_RECORD_PO_HPP
#define REAK_DATA_RECORD_PO_HPP

#include "ReaK/core/recorders/data_record.hpp"
#include "ReaK/core/recorders/data_record_options.hpp"

#include <filesystem>

#include "absl/flags/flag.h"

namespace ReaK::recorder {

/**
 * This function constructs a data-streaming options object from Abseil flags.
 * This function can either construct options for input or output.
 */
data_stream_options get_data_stream_options_from_flags(bool aForOutput = false);

}  // namespace ReaK::recorder

#endif  // REAK_DATA_RECORD_PO_HPP
