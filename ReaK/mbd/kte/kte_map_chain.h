/**
 * \file kte_map_chain.h
 *
 * This library declares the kte_map_chain class which allows KTE models to be grouped into a
 * linear chain of KTE models. The principle is based on the fact that if all the KTEs
 * compute their kinematics in order, the dynamics can be computed in the exact reverse order.
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

#ifndef REAK_MBD_KTE_KTE_MAP_CHAIN_H_
#define REAK_MBD_KTE_KTE_MAP_CHAIN_H_

#include <vector>
#include "ReaK/mbd/kte/kte_map.h"

namespace ReaK::kte {

/**
 * This class is a container of linearly chained kinetostatic transmission elements (KTEs).
 */
class kte_map_chain : public kte_map {
 private:
  std::vector<std::shared_ptr<kte_map>> mKTEs;  ///< Stores the list of KTEs

 public:
  /**
   * This function returns the vector of KTEs contained in this chain.
   * \return A const-reference to the vector of KTEs contained in this chain.
   */
  const std::vector<std::shared_ptr<kte_map>>& getKTEs() const { return mKTEs; }

  /**
   * Default constructor.
   */
  explicit kte_map_chain(const std::string& aName = "") : kte_map(aName) {}

  /**
   * Default destructor.
   */
  ~kte_map_chain() override = default;

  void doMotion(kte_pass_flag aFlag = nothing,
                const std::shared_ptr<frame_storage>& aStorage =
                    std::shared_ptr<frame_storage>()) override {
    auto it = mKTEs.begin();
    for (; it != mKTEs.end(); ++it) {
      (*it)->doMotion(aFlag, aStorage);
    }
  }

  void doForce(kte_pass_flag aFlag = nothing,
               const std::shared_ptr<frame_storage>& aStorage =
                   std::shared_ptr<frame_storage>()) override {
    auto rit = mKTEs.rbegin();
    for (; rit != mKTEs.rend(); ++rit) {
      (*rit)->doForce(aFlag, aStorage);
    }
  }

  void clearForce() override {
    auto it = mKTEs.begin();
    for (; it != mKTEs.end(); ++it) {
      (*it)->clearForce();
    }
  }

  /**
   * This method appends a KTE model to the linear chain of models.
   * \pre Any state.
   * \post The chain will contain an additional KTE.
   * \param aKTE The KTE model to add to this chain.
   * \return reference to this chain (to chain the << operators).
   */
  kte_map_chain& operator<<(const std::shared_ptr<kte_map>& aKTE) {
    if (aKTE) {
      mKTEs.push_back(aKTE);
    }
    return *this;
  }

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mKTEs);
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mKTEs);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(kte_map_chain, 0xC2100002, 1, "kte_map_chain",
                              kte_map)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_KTE_MAP_CHAIN_H_
