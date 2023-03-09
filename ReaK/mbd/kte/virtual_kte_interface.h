/**
 * \file virtual_kte_interface.h
 *
 * This library declares virtual model interface classes that make the boundary from a model
 * which corresponds to a real object or mechanism and a model which is virtual. This done
 * basically with Action / Reaction principle, i.e., the forces are reversed.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_MBD_KTE_VIRTUAL_KTE_INTERFACE_H_
#define REAK_MBD_KTE_VIRTUAL_KTE_INTERFACE_H_

#include "ReaK/mbd/kte/kte_map.h"

#include <utility>
#include "ReaK/math/kinetostatics/kinetostatics.h"

namespace ReaK::kte {

/**
 * This class implements a virtual model interface that makes the boundary from a model
 * which corresponds to a real object or mechanism and a model which is virtual. This done
 * basically with Action / Reaction principle, i.e., the forces are reversed.
 */
class virtual_kte_interface_gen : public kte_map {
 private:
  std::shared_ptr<gen_coord<double>>
      mBase;  ///< Holds the base frame of the interface.
  std::shared_ptr<gen_coord<double>>
      mEnd;  ///< Holds the end frame of the interface.

 public:
  /**
   * Sets the base frame of the interface (input frame).
   * \param aPtr The new base frame of the interface (input frame).
   */
  void setBaseFrame(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the base frame of the interface (input frame).
   * \return The base frame of the interface (input frame).
   */
  std::shared_ptr<gen_coord<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the end frame of the interface (output frame).
   * \param aPtr The new end frame of the interface (output frame).
   */
  void setEndFrame(const std::shared_ptr<gen_coord<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end frame of the interface (output frame).
   * \return The end frame of the interface (output frame).
   */
  std::shared_ptr<gen_coord<double>> EndFrame() const { return mEnd; }

  /**
   * Default constructor.
   */
  explicit virtual_kte_interface_gen(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aBase base frame of the interface.
   * \param aEnd end frame of the interface.
   */
  virtual_kte_interface_gen(const std::string& aName,
                            std::shared_ptr<gen_coord<double>> aBase,
                            std::shared_ptr<gen_coord<double>> aEnd)
      : kte_map(aName), mBase(std::move(aBase)), mEnd(std::move(aEnd)) {}

  /**
   * Default destructor.
   */
  ~virtual_kte_interface_gen() override = default;

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
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(virtual_kte_interface_gen, 0xC2100013, 1,
                              "virtual_kte_interface_gen", kte_map)
};

/**
 * This class implements a virtual model interface that makes the boundary from a model
 * which corresponds to a real object or mechanism and a model which is virtual. This done
 * basically with Action / Reaction principle, i.e., the forces are reversed.
 */
class virtual_kte_interface_2D : public kte_map {
 private:
  std::shared_ptr<frame_2D<double>>
      mBase;  ///< Holds the base frame of the interface.
  std::shared_ptr<frame_2D<double>>
      mEnd;  ///< Holds the end frame of the interface.

 public:
  /**
   * Sets the base frame of the interface (input frame).
   * \param aPtr The new base frame of the interface (input frame).
   */
  void setBaseFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the base frame of the interface (input frame).
   * \return The base frame of the interface (input frame).
   */
  std::shared_ptr<frame_2D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the end frame of the interface (output frame).
   * \param aPtr The new end frame of the interface (output frame).
   */
  void setEndFrame(const std::shared_ptr<frame_2D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end frame of the interface (output frame).
   * \return The end frame of the interface (output frame).
   */
  std::shared_ptr<frame_2D<double>> EndFrame() const { return mEnd; }

  /**
   * Default constructor.
   */
  explicit virtual_kte_interface_2D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aBase base frame of the interface.
   * \param aEnd end frame of the interface.
   */
  virtual_kte_interface_2D(const std::string& aName,
                           std::shared_ptr<frame_2D<double>> aBase,
                           std::shared_ptr<frame_2D<double>> aEnd)
      : kte_map(aName), mBase(std::move(aBase)), mEnd(std::move(aEnd)) {}

  /**
   * Default destructor.
   */
  ~virtual_kte_interface_2D() override = default;

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
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(virtual_kte_interface_2D, 0xC2100014, 1,
                              "virtual_kte_interface_2D", kte_map)
};

/**
 * This class implements a virtual model interface that makes the boundary from a model
 * which corresponds to a real object or mechanism and a model which is virtual. This done
 * basically with Action / Reaction principle, i.e., the forces are reversed.
 */
class virtual_kte_interface_3D : public kte_map {
 private:
  std::shared_ptr<frame_3D<double>>
      mBase;  ///< Holds the base frame of the interface.
  std::shared_ptr<frame_3D<double>>
      mEnd;  ///< Holds the end frame of the interface.

 public:
  /**
   * Sets the base frame of the interface (input frame).
   * \param aPtr The new base frame of the interface (input frame).
   */
  void setBaseFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mBase = aPtr;
  }
  /**
   * Returns the base frame of the interface (input frame).
   * \return The base frame of the interface (input frame).
   */
  std::shared_ptr<frame_3D<double>> BaseFrame() const { return mBase; }

  /**
   * Sets the end frame of the interface (output frame).
   * \param aPtr The new end frame of the interface (output frame).
   */
  void setEndFrame(const std::shared_ptr<frame_3D<double>>& aPtr) {
    mEnd = aPtr;
  }
  /**
   * Returns the end frame of the interface (output frame).
   * \return The end frame of the interface (output frame).
   */
  std::shared_ptr<frame_3D<double>> EndFrame() const { return mEnd; }

  /**
   * Default constructor.
   */
  explicit virtual_kte_interface_3D(const std::string& aName = "")
      : kte_map(aName) {}

  /**
   * Parametrized constructor.
   * \param aName name of the KTE model.
   * \param aBase base frame of the interface.
   * \param aEnd end frame of the interface.
   */
  virtual_kte_interface_3D(const std::string& aName,
                           std::shared_ptr<frame_3D<double>> aBase,
                           std::shared_ptr<frame_3D<double>> aEnd)
      : kte_map(aName), mBase(std::move(aBase)), mEnd(std::move(aEnd)) {}

  /**
   * Default destructor.
   */
  ~virtual_kte_interface_3D() override = default;

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
    A& RK_SERIAL_SAVE_WITH_NAME(mBase) & RK_SERIAL_SAVE_WITH_NAME(mEnd);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    kte_map::load(A, kte_map::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mBase) & RK_SERIAL_LOAD_WITH_NAME(mEnd);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(virtual_kte_interface_3D, 0xC2100015, 1,
                              "virtual_kte_interface_3D", kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_VIRTUAL_KTE_INTERFACE_H_
