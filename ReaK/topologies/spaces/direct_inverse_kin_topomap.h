/**
 * \file direct_inverse_kin_topomap.h
 *
 * This library provides class that define topological mapping between a joint-space of one
 * kinematic model and the joint-space of a serial manipulator, by assuming that the dependent
 * frames of the first kinematic models should mate exactly with the end-effector frames of the
 * manipulator (e.g., a chaser-target scenario).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_TOPOLOGIES_SPACES_DIRECT_INVERSE_KIN_TOPOMAP_H_
#define REAK_TOPOLOGIES_SPACES_DIRECT_INVERSE_KIN_TOPOMAP_H_


#include "ReaK/mbd/models/inverse_kinematics_model.h"

namespace ReaK::pp {

/**
 * This class template implements topological mapping between a joint-space of one
 * kinematic model and the joint-space of a serial manipulator, by assuming that the dependent
 * frames of the first kinematic models should mate exactly with the end-effector frames of the
 * manipulator (e.g., a chaser-target scenario).
 * \tparam DKMapType The direct-kinematics topological map type.
 * \tparam IKMapType The inverse-kinematics topological map type.
 */
template <typename DKMapType, typename IKMapType>
class manip_dk_ik_map : public shared_object {
 public:
  using self = manip_dk_ik_map<DKMapType, IKMapType>;

  DKMapType dk_map;
  IKMapType ik_map;

  /**
   * Parametrized Constructor.
   * \param aDKMap The direct-kinematics topological map.
   * \param aIKMap The inverse-kinematics topological map.
   */
  explicit manip_dk_ik_map(const DKMapType& aDKMap = DKMapType(),
                           const IKMapType& aIKMap = IKMapType())
      : dk_map(aDKMap), ik_map(aIKMap) {}

  /**
   * This function template performs a inverse kinematics calculation on the
   * manipulator model.
   * \tparam PointType The point-type of the input space.
   * \tparam InSpace The type of the input space (end-effector space).
   * \tparam OutSpace The type of the output space (joint-space).
   * \param pt The point in the input space, i.e. the end-effector coordinates.
   * \param space_in The input space, i.e. the end-effector space.
   * \param space_out The output space, i.e. the joint-space.
   * \return A point in the output space, i.e. the joint coordinates.
   */
  template <typename PointType, typename InSpace, typename OutSpace>
  typename topology_traits<OutSpace>::point_type map_to_space(
      const PointType& pt, const InSpace& space_in,
      const OutSpace& space_out) const {

    dk_map.apply_to_model(pt, space_in);

    // go through the dependent frames of each model to mate them together.
    std::size_t depgen_count = dk_map.model->getDependentCoordsCount();
    if (depgen_count > ik_map.model->getDependentCoordsCount()) {
      depgen_count = ik_map.model->getDependentCoordsCount();
    }
    for (std::size_t i = 0; i < depgen_count; ++i) {
      *(ik_map.model->getDependentCoord(i)->mFrame) =
          *(dk_map.model->getDependentCoord(i)->mFrame);
    }

    std::size_t dep2d_count = dk_map.model->getDependentFrames2DCount();
    if (dep2d_count > ik_map.model->getDependentFrames2DCount()) {
      dep2d_count = ik_map.model->getDependentFrames2DCount();
    }
    for (std::size_t i = 0; i < dep2d_count; ++i) {
      std::shared_ptr<frame_2D<double>> EE_frame =
          ik_map.model->getDependentFrame2D(i)->mFrame;
      *EE_frame = *(dk_map.model->getDependentFrame2D(i)->mFrame);
    };

    std::size_t dep3d_count = dk_map.model->getDependentFrames3DCount();
    if (dep3d_count > ik_map.model->getDependentFrames3DCount()) {
      dep3d_count = ik_map.model->getDependentFrames3DCount();
    }
    for (std::size_t i = 0; i < dep3d_count; ++i) {
      std::shared_ptr<frame_3D<double>> EE_frame =
          ik_map.model->getDependentFrame3D(i)->mFrame;
      *EE_frame = *(dk_map.model->getDependentFrame3D(i)->mFrame);
    };

    return ik_map.extract_from_model(space_out);
  }

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*Version*/) const override {
    shared_object::save(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_SAVE_WITH_NAME(dk_map) & RK_SERIAL_SAVE_WITH_NAME(ik_map);
  }
  void load(ReaK::serialization::iarchive& A,
            unsigned int /*Version*/) override {
    shared_object::load(A, shared_object::get_static_object_type()->version());
    A& RK_SERIAL_LOAD_WITH_NAME(dk_map) & RK_SERIAL_LOAD_WITH_NAME(ik_map);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(self, 0xC2400018, 1, "manip_dk_ik_map",
                              shared_object)
};

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_DIRECT_INVERSE_KIN_TOPOMAP_H_
