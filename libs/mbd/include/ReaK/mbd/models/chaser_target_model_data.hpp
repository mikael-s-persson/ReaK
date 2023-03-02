/**
 * \file chaser_target_model_data.hpp
 *
 * This library defines a class that is a meant to hold a collection of modeling data related to an
 * interception scenario (in 3D) where we have a chaser model which is meant to intercept with a target.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_CHASER_TARGET_MODEL_DATA_HPP
#define REAK_CHASER_TARGET_MODEL_DATA_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/named_object.hpp>

#include "inverse_kinematics_model.hpp"
#include "joint_space_limits.hpp"

namespace ReaK {

/* forward-declarations */
namespace geom {
class proxy_query_model_3D;
class colored_model_3D;
class proxy_query_pair_3D;
}  // namespace geom

namespace kte {

/**
 * This class is a meant to hold a collection of modeling data related to an
 * interception scenario (in 3D) where we have a chaser model (inverse-kinematics model)
 * which is meant to intercept with a target, represented as a direct-kinematics model,
 * as there is no obvious need for an inverse mapping on the side of the target.
 * This class exposes all the data members publicly because it is not a "functional" class
 * that exposes some sort of abstract functional interface, it is a data-collection class
 * that just acts as a convenient repository for the data of this scenario, and also for
 * the convenience of loading / saving the data.
 */
class chaser_target_data : public named_object {
 public:
  /** Holds a pointer to the base-frame (reference frame) of the chaser's coordinate system. */
  std::shared_ptr<frame_3D<double>> chaser_base_frame;
  /** Holds the chaser's inverse-kinematics model (polymorphic). */
  std::shared_ptr<inverse_kinematics_model> chaser_kin_model;
  /** Holds the chaser's joint limits, for velocity, acceleration and jerk. */
  std::shared_ptr<joint_limits_collection<double>> chaser_jt_limits;
  /** Holds the chaser's proximity-query model, i.e., for checking collisions involving the chaser. */
  std::shared_ptr<geom::proxy_query_model_3D> chaser_proxy;
  /** Holds the chaser's geometric model, i.e., for displaying the chaser's 3D representation. */
  std::shared_ptr<geom::colored_model_3D> chaser_geom_model;

  /** Holds the target's direct-kinematics model (polymorphic). */
  std::shared_ptr<direct_kinematics_model> target_kin_model;
  /** Holds the target's frame, i.e., the goal of the interception is that the chaser's end-effector frame meets this
   * target frame. */
  std::shared_ptr<frame_3D<double>> target_frame;
  /** Holds the target's proximity-query model, i.e., for checking collisions involving the target. */
  std::shared_ptr<geom::proxy_query_model_3D> target_proxy;
  /** Holds the target's geometric model, i.e., for displaying the target's 3D representation. */
  std::shared_ptr<geom::colored_model_3D> target_geom_model;

  /** Holds the chaser vs. target proximity-query pair, i.e., for checking collisions between the chaser and the target.
   */
  std::shared_ptr<geom::proxy_query_pair_3D> chaser_target_proxy;

  /** Holds the environment's geometric models, i.e., for displaying, in 3D, objects of the environment (static in the
   * chaser's base-frame). */
  std::vector<std::shared_ptr<geom::colored_model_3D>> env_geom_models;
  /** Holds the environment's proximity-query models, i.e., for checking collisions involving objects of the environment
   * (static in the chaser's base-frame). */
  std::vector<std::shared_ptr<geom::proxy_query_model_3D>> env_proxy_models;

  /** Holds the chaser vs. environment proximity-query pairs, i.e., for checking collisions between the chaser and the
   * environment. */
  std::vector<std::shared_ptr<geom::proxy_query_pair_3D>> chaser_env_proxies;
  /** Holds the target vs. environment proximity-query pairs, i.e., for checking collisions between the target and the
   * environment. */
  std::vector<std::shared_ptr<geom::proxy_query_pair_3D>> target_env_proxies;

  /**
   * Default constructor. Leaves all data members empty.
   */
  chaser_target_data();

  /**
   * This function loads the chaser models from a given filename. The file is expected to have
   * all the models laid out in the following order: chaser_base_frame, chaser_kin_model,
   * chaser_jt_limits, chaser_geom_model, and chaser_proxy.
   * \param fileName The filename of the file from which to load the chaser models.
   */
  void load_chaser(const std::string& fileName);

  /**
   * This function saves the chaser models to a given filename. The file will have
   * all the models laid out in the following order: chaser_base_frame, chaser_kin_model,
   * chaser_jt_limits, chaser_geom_model, and chaser_proxy.
   * \param fileName The filename of the file to which to save the chaser models.
   */
  void save_chaser(const std::string& fileName) const;

  /**
   * This function loads the target models from a given filename. The file is expected to have
   * all the models laid out in the following order: target_base, target_kin_model, target_frame,
   * target_geom_model, and target_proxy.
   * \param fileName The filename of the file from which to load the target models.
   */
  void load_target(const std::string& fileName);

  /**
   * This function saves the target models to a given filename. The file will have
   * all the models laid out in the following order: target_base, target_kin_model, target_frame,
   * target_geom_model, and target_proxy.
   * \param fileName The filename of the file to which to save the target models.
   */
  void save_target(const std::string& fileName) const;

  /**
   * This function loads (and adds) an environment object from a given filename. The file is
   * expected to have the models laid out in the following order: env_geom_model and env_proxy.
   * \param fileName The filename of the file from which to load the environment object.
   */
  void load_environment(const std::string& fileName);

  /**
   * This function saves an environment object to a given filename. The file will
   * have the models laid out in the following order: env_geom_model and env_proxy.
   * \param id The index of the environment object in the array of environment geometries.
   * \param fileName The filename of the file to which to save the environment object.
   */
  void save_environment(std::size_t id, const std::string& fileName) const;

  /**
   * This function clears the list of environment objects (and the related proximity models).
   */
  void clear_environment();

 private:
  void create_chaser_target_proxy();
  void create_chaser_env_proxies();
  void create_target_env_proxies();

 public:
  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(serialization::oarchive& A, unsigned int /*unused*/) const override;
  void load(serialization::iarchive& A, unsigned int /*unused*/) override;

  RK_RTTI_MAKE_CONCRETE_1BASE(chaser_target_data, 0xC210005A, 1,
                              "chaser_target_data", named_object)
};
}  // namespace kte
}  // namespace ReaK

#endif
