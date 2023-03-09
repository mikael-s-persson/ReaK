/**
 * \file joint_backlash.h
 *
 * This library does not contain any working class. An attempt was made at modeling joints
 * with backlash, but this attempt failed.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2010
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

#ifndef REAK_MBD_KTE_JOINT_BACKLASH_H_
#define REAK_MBD_KTE_JOINT_BACKLASH_H_

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * NOT A WORKING MODEL.
 */
class joint_backlash_gen : public kte_map {
 private:
  std::shared_ptr<gen_coord<double>> mBase;
  std::shared_ptr<gen_coord<double>> mEnd;
  double mGapSize;

 public:
  double& GapSize() { return mGapSize; }
  double GapSize() const { return mGapSize; }

  /**
   * Default constructor.
   */
  explicit joint_backlash_gen(const std::string& aName = "")
      : kte_map(aName), mGapSize(0.0) {}

  /**
   * Parametrized constructor.
   */
  joint_backlash_gen(const std::string& aName,
                     std::shared_ptr<gen_coord<double>> aBase,
                     std::shared_ptr<gen_coord<double>> aEnd, double aGapSize)
      : kte_map(aName),
        mBase(std::move(aBase)),
        mEnd(std::move(aEnd)),
        mGapSize(aGapSize) {}

  /**
   * Default destructor.
   */
  ~joint_backlash_gen() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override;

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override;

  void clearForce() override;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    kte_map::save(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd) &
        RK_SERIAL_SAVE_WITH_NAME(mGapSize);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd) &
        RK_SERIAL_LOAD_WITH_NAME(mGapSize);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(joint_backlash_gen, 0xC2100022, 1,
                              "joint_backlash_gen", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_JOINT_BACKLASH_H_
