/**
 * \file planning_space_options_po.h
 *
 * This library defines functions to load planning-space options from Boost.Program-Options.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date February 2014
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

#ifndef REAK_PLANNING_PATH_PLANNING_PLANNING_SPACE_OPTIONS_PO_H_
#define REAK_PLANNING_PATH_PLANNING_PLANNING_SPACE_OPTIONS_PO_H_

#include "ReaK/planning/path_planning/planning_space_options.h"

namespace ReaK::pp {

planning_space_options get_planning_space_options_from_flags();

}  // namespace ReaK::pp

#endif  // REAK_PLANNING_PATH_PLANNING_PLANNING_SPACE_OPTIONS_PO_H_
