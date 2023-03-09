/**
 * \file kte_map.h
 *
 * This library declares the basic class for all KTE models. Basically, kte_map is an interface
 * containing the three fundamental functions doMotion, clearForce, and doForce.
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

#ifndef REAK_MBD_KTE_KTE_MAP_H_
#define REAK_MBD_KTE_KTE_MAP_H_

#include "ReaK/core/base/named_object.h"
#include "ReaK/math/kinetostatics/kinetostatics.h"

#include "ReaK/mbd/kte/kte_ext_mappings.h"

namespace ReaK::kte {

/**
 * These flags signify which type of KTE-pass is requested to be performed, for self-contained algorithms.
 */
enum kte_pass_flag {
  nothing,  ///< Does nothing special, default mode.
  store_kinematics,  ///< Request to store the calculated kinematics value (in a doMotion pass).
  store_dynamics,  ///< Request to store the calculated force values (in a doForce pass).
  adapt_parameters  ///< Request to adapt the system parameters based on recorded kinematics and forces.
};

/**
 * This class is the base class for all the kinetostatic transmission elements (KTE).
 */
class kte_map : public virtual named_object {
 public:
  /**
   * Default constructor.
   */
  explicit kte_map(const std::string& aName) { this->setName(aName); }

  kte_map() : kte_map("") {}

  // Non-copyable.
  kte_map(const kte_map&) = delete;
  kte_map& operator=(const kte_map&) = delete;

  /**
   * Default destructor.
   */
  ~kte_map() override = default;

  /**
   * This method performs a motion calculation pass (i.e. computes the kinematics).
   * \pre All kinematics inputs have been set to the proper values.
   * \post All kinematics outputs have been set and the internal kinematics-related states are updated.
   * \param aFlag a flag for the type of pass requested, can be nothing, store_kinematics, or adapt_parameters.
   * \param aStorage a map of kinetostatic frame storage for stored kinematics and dynamics.
   * \note adaptivity is not yet supported.
   */
  virtual void doMotion(kte_pass_flag aFlag = nothing,
                        const std::shared_ptr<frame_storage>& aStorage =
                            std::shared_ptr<frame_storage>()) = 0;

  /**
   * This method performs a force calculation pass (i.e. computes the dynamics).
   * \pre All dynamics inputs have been set to the proper values.
   * \post All dynamics outputs have been set and the internal dynamics-related states are updated.
   * \param aFlag a flag for the type of pass requested, can be nothing, store_dynamics, or adapt_parameters.
   * \param aStorage a map of kinetostatic frame storage for stored kinematics and dynamics.
   * \note adaptivity is not yet supported.
   */
  virtual void doForce(kte_pass_flag aFlag = nothing,
                       const std::shared_ptr<frame_storage>& aStorage =
                           std::shared_ptr<frame_storage>()) = 0;

  /**
   * This method performs a force clearance pass (i.e. resets all force values to zero).
   *
   * \pre Any state, but normally before making a new doMotion + doForce pair of calculation passes.
   * \post All force outputs have been set to zero.
   *
   * \note This is required because of multiple writer schemes, each KTE ADDS FORCES AS OUTPUT, so not
   * starting from zero will have undefined results.
   */
  virtual void clearForce() = 0;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(kte_map, 0xC2100001, 1, "kte_map", named_object)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_KTE_MAP_H_
