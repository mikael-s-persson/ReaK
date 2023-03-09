/**
 * \file oi_reader.h
 *
 * This library declares a class that can act as an intermediate between
 * an OpenInventor scene-graph (via Coin3D) and a collection of geometric entities for ReaK.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2013
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

#ifndef REAK_MBD_COIN3D_OI_READER_H_
#define REAK_MBD_COIN3D_OI_READER_H_

#include "ReaK/core/base/defs.h"

#include "ReaK/geometry/proximity/proxy_query_model.h"
#include "ReaK/geometry/shapes/color.h"
#include "ReaK/geometry/shapes/colored_model.h"

#include <map>
#include <vector>

// forward-declarations of the open-inventory node classes:
class SoSeparator;

namespace ReaK {

template <typename T>
class pose_2D;
template <typename T>
class pose_3D;

namespace geom {

class geometry_2D;
class geometry_3D;

/**
 * This class acts as an intermediate between an OpenInventor scene-graph (via Coin3D) and a
 * collection of geometric entities for ReaK.
 */
class oi_reader {
 protected:
  SoSeparator* mRoot;

 public:
  SoSeparator* getSceneGraph() const { return mRoot; }

  /**
   * This function computes the characteristic length used to scale the illustrative components (e.g., coordinate axes,
   * KTE representations, etc.) to a size that is proportionate to the geometries by using the bounding box computed
   * from the geometries currently present in the scene-graph.
   * \return The value obtained for the characteristic length, which will also be set internally as the current
   * characteristic-length for the scene-graph.
   */
  double computeCharacteristicLength();

  /**
   * Default constructor.
   */
  oi_reader();

  /**
   * Default constructor from a file-name.
   */
  explicit oi_reader(const std::string& aFileName);

  /**
   * Default constructor from an input stream.
   */
  explicit oi_reader(std::istream& aStream);

  oi_reader(const oi_reader& /*rhs*/);
  oi_reader& operator=(const oi_reader& /*rhs*/);

  oi_reader(oi_reader&& /*rhs*/) noexcept;
  oi_reader& operator=(oi_reader&& /*rhs*/) noexcept;

  /**
   * Default destructor.
   */
  virtual ~oi_reader();

  /**
   * Open a given file.
   */
  oi_reader& read_file(const std::string& aFileName);

  /**
   * Read scene-graph from a given stream.
   */
  oi_reader& read_stream(std::istream& aStream);

  /**
   * Check if the.
   */
  explicit operator bool() const { return (mRoot != nullptr); }

  friend oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel);

  friend oi_reader& operator>>(oi_reader& aSG, proxy_query_model_3D& aProxy);

  oi_reader& translate_into(colored_model_3D& aModel,
                            proxy_query_model_3D& aProxy);
};

// Re-declaration down in the geom namespace directly as some compilers give trouble with friend functions only declared
// in the class.

/**
 * This operator creates a 3D colored geometric model from an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph reader to construct the geometric model from.
 * \param aModel The 3D colored geometric model constructed by the scene-graph reader.
 * \return The Open-Inventor scene-graph reader given as the first parameter, by reference.
 */
oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel);

/**
 * This operator creates a 3D proximity-query model from an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph reader to construct the proximity-query model from.
 * \param aProxy The 3D proximity-query model constructed by the scene-graph reader.
 * \return The Open-Inventor scene-graph reader given as the first parameter, by reference.
 */
oi_reader& operator>>(oi_reader& aSG, proxy_query_model_3D& aProxy);

}  // namespace geom
}  // namespace ReaK

#endif  // REAK_MBD_COIN3D_OI_READER_H_
