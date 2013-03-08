/**
 * \file oi_reader.hpp
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

#ifndef REAK_OI_READER_HPP
#define REAK_OI_READER_HPP

#include "base/defs.hpp"

#include "color.hpp"
#include "colored_model.hpp"

#include "mbd_kte/kte_map.hpp"          // for kte_map
#include "mbd_kte/kte_map_chain.hpp"

#include <map>
#include <vector>

// forward-declarations of the open-inventory node classes:
class SoSeparator;


/** Main namespace for ReaK */
namespace ReaK {
  
template <typename T> class pose_2D;
template <typename T> class pose_3D;

/** Main namespace for ReaK.Geometry */
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
    
    SoSeparator* getSceneGraph() const { return mRoot; };
    
    /**
     * This function computes the characteristic length used to scale the illustrative components (e.g., coordinate axes, 
     * KTE representations, etc.) to a size that is proportionate to the geometries by using the bounding box computed
     * from the geometries currently present in the scene-graph.
     * \return The value obtained for the characteristic length, which will also be set internally as the current characteristic-length for the scene-graph.
     */
    double computeCharacteristicLength();
    
    /**
     * Default constructor.
     */
    oi_reader();
    
    /**
     * Default constructor.
     */
    oi_reader(const std::string& aFileName);
    
    /**
     * Default destructor.
     */
    virtual ~oi_reader();
    
    /**
     * Open a given file.
     */
    oi_reader& read_file(const std::string& aFileName);
    
    /**
     * Check if the.
     */
    operator bool() const { return (mRoot != NULL); };
    
    
    friend
    oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel);
    
};


// Re-declaration down in the geom namespace directly as some compilers give trouble with friend functions only declared in the class.

/**
 * This operator adds a 3D colored geometric model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the geometric model to.
 * \param aModel The 3D colored geometric model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_reader& operator>>(oi_reader& aSG, colored_model_3D& aModel);

};

};

#endif










