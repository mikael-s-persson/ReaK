/**
 * \file kte_chain_geometry.h
 *
 * This library declares a class
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2013
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

#ifndef REAK_MBD_KTE_KTE_CHAIN_GEOMETRY_H_
#define REAK_MBD_KTE_KTE_CHAIN_GEOMETRY_H_

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/geometry/shapes/colored_model.h"
#include "ReaK/geometry/shapes/shape_2D.h"
#include "ReaK/geometry/shapes/shape_3D.h"

#include <map>
#include <string>

namespace ReaK::kte {

class kte_map_chain;

/** This class defines a colored model for 2D geometries. */
class kte_chain_geometry_2D : public named_object {
 public:
  std::map<std::string, std::vector<geom::colored_geometry_2D>> mGeomList;
  std::map<std::string, std::vector<std::shared_ptr<geom::shape_2D>>>
      mProxyShapeList;

  /**
   * Default constructor.
   */
  explicit kte_chain_geometry_2D(const std::string& aName = "") {
    this->setName(aName);
  }

  /**
   * Default destructor.
   */
  ~kte_chain_geometry_2D() override = default;

  kte_chain_geometry_2D& addElement(
      const std::string& aKTEObjName, const geom::color& aColor,
      const std::shared_ptr<geom::geometry_2D>& aGeom) {
    mGeomList[aKTEObjName].push_back(geom::colored_geometry_2D(aColor, aGeom));
    return *this;
  }

  kte_chain_geometry_2D& addShape(
      const std::string& aKTEObjName,
      const std::shared_ptr<geom::shape_2D>& aShape) {
    mProxyShapeList[aKTEObjName].push_back(aShape);
    return *this;
  }

  std::pair<std::shared_ptr<geom::colored_model_2D>,
            std::shared_ptr<geom::proxy_query_model_2D>>
  attachToKTEChain(const kte::kte_map_chain& aKTEChain) const;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mGeomList) &
        RK_SERIAL_SAVE_WITH_NAME(mProxyShapeList);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mGeomList) &
        RK_SERIAL_LOAD_WITH_NAME(mProxyShapeList);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(kte_chain_geometry_2D, 0xC3100024, 1,
                              "kte_chain_geometry_2D", named_object)
};

/** This class defines a colored model for 3D geometries. */
class kte_chain_geometry_3D : public named_object {
 public:
  std::map<std::string, std::vector<geom::colored_geometry_3D>> mGeomList;
  std::map<std::string, std::vector<std::shared_ptr<geom::shape_3D>>>
      mProxyShapeList;

  /**
   * Default constructor.
   */
  explicit kte_chain_geometry_3D(const std::string& aName = "") {
    this->setName(aName);
  }

  /**
   * Default destructor.
   */
  ~kte_chain_geometry_3D() override = default;

  kte_chain_geometry_3D& addElement(
      const std::string& aKTEObjName, const geom::color& aColor,
      const std::shared_ptr<geom::geometry_3D>& aGeom) {
    mGeomList[aKTEObjName].push_back(geom::colored_geometry_3D(aColor, aGeom));
    return *this;
  }

  kte_chain_geometry_3D& addShape(
      const std::string& aKTEObjName,
      const std::shared_ptr<geom::shape_3D>& aShape) {
    mProxyShapeList[aKTEObjName].push_back(aShape);
    return *this;
  }

  std::pair<std::shared_ptr<geom::colored_model_3D>,
            std::shared_ptr<geom::proxy_query_model_3D>>
  attachToKTEChain(const kte::kte_map_chain& aKTEChain) const;

  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  void save(ReaK::serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    named_object::save(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_SAVE_WITH_NAME(mGeomList) &
        RK_SERIAL_SAVE_WITH_NAME(mProxyShapeList);
  }

  void load(ReaK::serialization::iarchive& A,
            unsigned int /*unused*/) override {
    named_object::load(A, named_object::getStaticObjectType()->TypeVersion());
    A& RK_SERIAL_LOAD_WITH_NAME(mGeomList) &
        RK_SERIAL_LOAD_WITH_NAME(mProxyShapeList);
  }

  RK_RTTI_MAKE_CONCRETE_1BASE(kte_chain_geometry_3D, 0xC3100025, 1,
                              "kte_chain_geometry_3D", named_object)
};

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_KTE_CHAIN_GEOMETRY_H_
