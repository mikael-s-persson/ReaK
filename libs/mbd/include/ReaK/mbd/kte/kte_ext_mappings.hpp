/**
 * \file kte_ext_mappings.hpp
 *
 * This library declares extended KTE-related mappings. Most importantly, the maps for storing
 * the kinetostatic frame states pertaining to end and intermediate frames in a chain of KTE models.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
 */

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

#ifndef REAK_KTE_EXT_MAPPINGS_HPP
#define REAK_KTE_EXT_MAPPINGS_HPP

#include "ReaK/math/kinetostatics/kinetostatics.hpp"
#include "ReaK/math/lin_alg/vect_alg.hpp"

#include <map>

namespace ReaK::kte {

struct velocity_coef_gen {
  double v{0.0};
  velocity_coef_gen() = default;
};

struct acceleration_coef_gen {
  double a{0.0};
  acceleration_coef_gen() = default;
};

struct force_coef_gen {
  double f{0.0};
  force_coef_gen() = default;
};

struct velocity_coef_2D {
  vect<double, 2> v;
  double omega{0.0};
  velocity_coef_2D() = default;
};

struct acceleration_coef_2D {
  vect<double, 2> a;
  double alpha{0.0};
  acceleration_coef_2D() = default;
};

struct force_coef_2D {
  vect<double, 2> f;
  double tau{0.0};
  force_coef_2D() = default;
};

struct velocity_coef_3D {
  vect<double, 3> v;
  vect<double, 3> omega;
  velocity_coef_3D() = default;
};

struct acceleration_coef_3D {
  vect<double, 3> a;
  vect<double, 3> alpha;
  acceleration_coef_3D() = default;
};

struct force_coef_3D {
  vect<double, 3> f;
  vect<double, 3> tau;
  force_coef_3D() = default;
};

using velocity_coef_map_gen = std::map<double*, velocity_coef_gen>;
using acceleration_coef_map_gen = std::map<double*, acceleration_coef_gen>;
using force_coef_map_gen = std::map<double*, force_coef_gen>;
using velocity_coef_map_2D = std::map<double*, velocity_coef_2D>;
using acceleration_coef_map_2D = std::map<double*, acceleration_coef_2D>;
using force_coef_map_2D = std::map<double*, force_coef_2D>;
using velocity_coef_map_3D = std::map<double*, velocity_coef_3D>;
using acceleration_coef_map_3D = std::map<double*, acceleration_coef_3D>;
using force_coef_map_3D = std::map<double*, force_coef_3D>;

/** This typedef defines a map of generalized coordinate state storage associated to a generalized coordinate pointer in
 * the KTE chain. */
using gen_coord_map = std::map<std::shared_ptr<gen_coord<double>>,
                               std::shared_ptr<gen_coord<double>>>;
/** This typedef defines a map of frame 2D state storage associated to a frame 2D pointer in the KTE chain. */
using frame_2D_map = std::map<std::shared_ptr<frame_2D<double>>,
                              std::shared_ptr<frame_2D<double>>>;
/** This typedef defines a map of frame 3D state storage associated to a frame 3D pointer in the KTE chain. */
using frame_3D_map = std::map<std::shared_ptr<frame_3D<double>>,
                              std::shared_ptr<frame_3D<double>>>;

/**
 * This struct is used as a storage repository to take a "flash" of all kinematics and dynamics states
 * of the kinetostatic frames (2D, 3D and generalized coordinates) at an instant.
 */
struct frame_storage {
  gen_coord_map
      gen_coord_mapping;          ///< Stores the generalized coordinate states.
  frame_2D_map frame_2D_mapping;  ///< Stores the frame 2D states.
  frame_3D_map frame_3D_mapping;  ///< Stores the frame 3D states.
  /** Default constructor, non-POD. */
  frame_storage() = default;
};

}  // namespace ReaK::kte

#endif
