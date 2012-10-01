/**
 * \file oi_scene_graph.hpp
 *
 * This library declares a class that can act as an intermediate between a collection of 
 * geometric entities from ReaK and an OpenInventor scene-graph (via Coin3D).
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date September 2012
 */

/*
 *    Copyright 2012 Sven Mikael Persson
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

#ifndef REAK_OI_SCENE_GRAPH_HPP
#define REAK_OI_SCENE_GRAPH_HPP

#include "geometry_2D.hpp"
#include "geometry_3D.hpp"

#include <map>

// forward-declarations of the open-inventory node classes:
class SoSeparator;
class SoTransform;
class SoSensor;
class SoTimerSensor;


/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class defines a . */
class oi_scene_graph {
  protected:
    SoSeparator* mRoot;
    SoTimerSensor* mTimer;
    
    std::map< shared_ptr< pose_2D<double> >, std::pair<SoSeparator*, SoTransform*> > mAnchor2DMap;
    std::map< shared_ptr< pose_3D<double> >, std::pair<SoSeparator*, SoTransform*> > mAnchor3DMap;
    
    static void update_anchors(void*, SoSensor*);
    
  public:
    
    SoSeparator* getSceneGraph() const { return mRoot; };
    
    void enableAnchorUpdates();
    void disableAnchorUpdates();
    
    /**
     * Default constructor.
     */
    oi_scene_graph();
    
    /**
     * Default destructor.
     */
    virtual ~oi_scene_graph();
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_2D<double> >& aAnchor);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_3D<double> >& aAnchor);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_2D& aGeom2D);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_3D& aGeom3D);
    
    template <typename T>
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const std::vector<T>& aItems) {
      for(typename std::vector<T>::const_iterator it = aItems.begin(); it != aItems.end(); ++it)
        aSG << (*it);
      return aSG;
    };
    
    
};


};

};

#endif










