/**
 * \file data_record_po.h
 *
 * This library declares utility functions for creating and dealing with program-options related
 * to a data recording or extraction stream (see data_record.h). Here, "data" is meant as
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

#ifndef REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_PO_H_
#define REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_PO_H_

#include "ReaK/control/systems/satellite_invar_models.h"
#include "ReaK/control/systems/satellite_modeling_options.h"

#include "absl/flags/declare.h"

ABSL_DECLARE_FLAG(bool, sat_generate_mdl_files);

namespace ReaK::ctrl {

/**
 * This function constructs a satellite-modeling options object from Abseil flags.
 */
satellite_model_options get_satellite_model_options_from_flags();

/**
 * This function constructs a satellite-predictor options object from Abseil flags.
 */
satellite_predictor_options get_satellite_predictor_options_from_flags();

}  // namespace ReaK::ctrl

#endif  // REAK_CONTROL_SYSTEMS_SATELLITE_MODELING_PO_H_
