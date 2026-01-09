/**
 * \file reacting_kte.h
 *
 * This library declares an interface for KTE models in which a reaction force feedback is necessary.
 * This mechanism has certain effects on the ordering of the KTE chain and thus, attention is advised
 * when using this mechanism.
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

#ifndef REAK_MBD_KTE_REACTING_KTE_H_
#define REAK_MBD_KTE_REACTING_KTE_H_

#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * This class declares the interface for a KTE model which receives a reaction force along a generalized coordinate.
 */
class reacting_kte_gen : public kte_map {
 public:
  /**
   * Default constructor.
   */
  explicit reacting_kte_gen(const std::string& aName = "") { set_name(aName); }

  /**
   * Default destructor.
   */
  ~reacting_kte_gen() override = default;

  /**
   * Applies to the reaction force along the relevant generalized coordinate.
   * \param aForce the reaction force to add to the sum of forces (could be a torque as well).
   */
  virtual void applyReactionForce(double aForce) = 0;

  RK_RTTI_MAKE_ABSTRACT_1BASE(reacting_kte_gen, 0xC2100016, 1,
                              "reacting_kte_gen", kte_map)
};

/**
 * This class declares the interface for a KTE model which receives a reaction force along a 2D frame.
 */
class reacting_kte_2D : public kte_map {
 public:
  /**
   * Default constructor.
   */
  explicit reacting_kte_2D(const std::string& aName = "") { set_name(aName); }

  /**
   * Default destructor.
   */
  ~reacting_kte_2D() override = default;

  /**
   * Applies to the reaction force along the relevant 2D frame.
   * \param aForce the reaction force vector to add to the sum of forces.
   * \param aTorque the reaction torque to add to the sum of torques.
   */
  virtual void applyReactionForce(vect<double, 2> aForce, double aTorque) = 0;

  RK_RTTI_MAKE_ABSTRACT_1BASE(reacting_kte_2D, 0xC2100017, 1, "reacting_kte_2D",
                              kte_map)
};

/**
 * This class declares the interface for a KTE model which receives a reaction force along a 3D frame.
 */
class reacting_kte_3D : public kte_map {
 public:
  /**
   * Default constructor.
   */
  explicit reacting_kte_3D(const std::string& aName = "") { set_name(aName); }

  /**
   * Default destructor.
   */
  ~reacting_kte_3D() override = default;

  /**
   * Applies to the reaction force along the relevant 3D frame.
   * \param aForce the reaction force vector to add to the sum of forces.
   * \param aTorque the reaction torque vector to add to the sum of torques.
   */
  virtual void applyReactionForce(vect<double, 3> aForce,
                                  vect<double, 3> aTorque) = 0;

  RK_RTTI_MAKE_ABSTRACT_1BASE(reacting_kte_3D, 0xC2100018, 1, "reacting_kte_3D",
                              kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_REACTING_KTE_H_
