/**
 * \file oi_scene_graph.hpp
 *
 * This library declares a class that can act as an intermediate between a collection of 
 * geometric entities from ReaK and an OpenInventor scene-graph (via Coin3D). First, this 
 * class acts somewhat like a standard iostream (with << operators) to sequentially add anchors, 
 * geometries, complete models, or KTE chains to the current scene-graph stored within an object 
 * of this class. Second, this class has internal mechanisms to update all the anchors (coordinate 
 * frames, or transformations) relating to the geometries such that the scene-graph is in-sync 
 * with those anchors. This mechanism is enabled and disabled with the enableAnchorUpdates() 
 * and disableAnchorUpdates() functions, respectively. 
 * 
 * For the most part, this class acts as a single dependency point between ReaK objects and 
 * Open-Inventor. Under-the-hood, this class checks the dynamic types of the objects it is 
 * given and constructs the appropriate Open-Inventor sub-graph for each type it can recognize 
 * and support. This is not pretty, but it is the lesser evil between that and the alternative,
 * which is to give each relevant objects in ReaK a virtual method (from a common base-class)
 * which can create the Open-Inventor sub-graph for itself, which would effectively couple 
 * the entire ReaK library to Open-Inventor (or require some awkward dummy linking scheme).
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

#include "base/defs.hpp"

#include "color.hpp"
#include "colored_model.hpp"

#include "mbd_kte/kte_map.hpp"          // for kte_map
#include "mbd_kte/kte_map_chain.hpp"

#include "proximity/proxy_query_model.hpp"

#include <map>
#include <vector>

#include "base/thread_incl.hpp"

#ifdef RK_ENABLE_CXX11_FEATURES
#include <functional>
#else
#include <boost/function.hpp>
#endif

// forward-declarations of the open-inventory node classes:
class SoSeparator;
class SoSwitch;
class SoTransform;
class SoSensor;
class SoTimerSensor;


/** Main namespace for ReaK */
namespace ReaK {
  
template <typename T> class pose_2D;
template <typename T> class pose_3D;

/** Main namespace for ReaK.Geometry */
namespace geom {

class geometry_2D;
class geometry_3D;

/** 
 * This class acts as an intermediate between various entities of ReaK and the construction of an 
 * Open-Inventor scene-graph (through the Coin3D library). First, this class acts somewhat like a 
 * standard iostream (with << operators) to sequentially add anchors, geometries, complete models, 
 * or KTE chains to the current scene-graph stored within an object of this class. Second, this 
 * class has internal mechanisms to update all the anchors (coordinate frames, or transformations)
 * relating to the geometries such that the scene-graph is in-sync with those anchors. This 
 * mechanism is enabled and disabled with the enableAnchorUpdates() and disableAnchorUpdates()
 * functions, respectively. 
 * 
 * For the most part, this class acts as a single dependency point between ReaK objects and 
 * Open-Inventor. Under-the-hood, this class checks the dynamic types of the objects it is 
 * given and constructs the appropriate Open-Inventor sub-graph for each type it can recognize 
 * and support. This is not pretty, but it is the lesser evil between that and the alternative,
 * which is to give each relevant objects in ReaK a virtual method (from a common base-class)
 * which can create the Open-Inventor sub-graph for itself, which would effectively couple 
 * the entire ReaK library to Open-Inventor (or require some awkward dummy linking scheme).
 */
class oi_scene_graph {
  protected:
    SoSeparator* mRoot;
    SoSwitch* mRootSwitch;
    SoTimerSensor* mTimer;
    color mCurrentColor;
    double mCharacteristicLength;
    
    ReaKaux::recursive_mutex mAnchorUpdatingMutex;
    
    std::map< shared_ptr< pose_2D<double> >, std::pair<SoSeparator*, SoTransform*> > mAnchor2DMap;
    std::map< shared_ptr< pose_3D<double> >, std::pair<SoSeparator*, SoTransform*> > mAnchor3DMap;

#ifdef RK_ENABLE_CXX0X_FEATURES
    std::vector< std::function< void() > > mUpdateFuncs;
#else
    std::vector< boost::function< void() > > mUpdateFuncs;
#endif
    
    static void update_anchors(void*, SoSensor*);
    
  public:
    
    SoSeparator* getSceneGraph() const { return mRoot; };
    
    void setVisibility(bool aVisible) const;
    
    void clearAll();
    
    /**
     * Enable the internal sensor-mechanism that will update all the transformations in the scene-graph
     * to synchronize them with their associated anchors (pose_2D or pose_3D objects).
     */
    void enableAnchorUpdates();
    /**
     * Disable the internal sensor-mechanism that updates all the transformations in the scene-graph
     * to synchronize them with their associated anchors (pose_2D or pose_3D objects).
     */
    void disableAnchorUpdates();
    
    /**
     * This function returns the current characteristic length used to scale the illustrative components (e.g., coordinate axes, 
     * KTE representations, etc.) to a size that is proportionate to the geometries.
     * \return The current value for the characteristic length of subsequent components to be added.
     */
    double getCharacteristicLength() const { return mCharacteristicLength; };
    
    /**
     * This function computes the characteristic length used to scale the illustrative components (e.g., coordinate axes, 
     * KTE representations, etc.) to a size that is proportionate to the geometries by using the bounding box computed
     * from the geometries currently present in the scene-graph.
     * \return The value obtained for the characteristic length, which will also be set internally as the current characteristic-length for the scene-graph.
     */
    double computeCharacteristicLength();
    
    /**
     * This function sets the characteristic length used to scale the illustrative components (e.g., coordinate axes, 
     * KTE representations, etc.) to a size that is proportionate to the geometries. Use the function 
     * computeCharacteristicLength() to have this value calculated from the bounding boxes of the geometries currently
     * included in the scene-graph.
     * \param aCharacteristicLength The new value for the characteristic length for subsequent components to be added.
     */
    void setCharacteristicLength(double aCharacteristicLength) { mCharacteristicLength = aCharacteristicLength; };
    
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
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const color& aColor);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_2D& aGeom2D);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_3D& aGeom3D);
    
    /**
     * This operator adds all the elements of a vector to an Open-Inventor scene-graph.
     * \param aSG An Open-Inventor scene-graph to add the geometry to.
     * \param aItems The items to add to the scene-graph.
     * \return The Open-Inventor scene-graph given as the first parameter, by reference.
     * \tparam T Any type for which there is an appropriate scene-graph << operator.
     */
    template <typename T>
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const std::vector<T>& aItems) {
      for(typename std::vector<T>::const_iterator it = aItems.begin(); it != aItems.end(); ++it)
        aSG << (*it);
      return aSG;
    };
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_2D& aModel);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_3D& aModel);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map& aModel);
    
    friend
    oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map_chain& aModel);
    
};


// Re-declaration down in the geom namespace directly as some compilers give trouble with friend functions only declared in the class.

/**
 * This operator adds a 2D anchor to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the anchor to.
 * \param aAnchor The anchor to add to the scene-graph (this will be included in the transformation updates).
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_2D<double> >& aAnchor);
    
/**
 * This operator adds a 3D anchor to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the anchor to.
 * \param aAnchor The anchor to add to the scene-graph (this will be included in the transformation updates).
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const shared_ptr< pose_3D<double> >& aAnchor);
    
/**
 * This operator sets the current color in an Open-Inventor scene-graph, this color will be used for 
 * subsequent geometries that are added to the scene-graph.
 * \param aSG An Open-Inventor scene-graph.
 * \param aColor The new color to use for subsequent geometries added to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const color& aColor);
    
/**
 * This operator adds a 2D geometry to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the geometry to.
 * \param aGeom2D The geometry to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_2D& aGeom2D);
    
/**
 * This operator adds a 3D geometry to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the geometry to.
 * \param aGeom3D The geometry to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const geometry_3D& aGeom3D);

/**
 * This operator adds a 2D colored geometric model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the geometric model to.
 * \param aModel The 2D colored geometric model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_2D& aModel);
    
/**
 * This operator adds a 3D colored geometric model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the geometric model to.
 * \param aModel The 3D colored geometric model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const colored_model_3D& aModel);

/**
 * This operator adds a 2D proximity-query model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the proximity-query model to.
 * \param aModel The 2D proximity-query model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const proxy_query_model_2D& aModel);
    
/**
 * This operator adds a 3D proximity-query model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the proximity-query model to.
 * \param aModel The 3D proximity-query model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const proxy_query_model_3D& aModel);

/**
 * This operator adds a KTE model to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the KTE model to.
 * \param aModel The KTE model to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map& aModel);

/**
 * This operator adds a KTE chain to an Open-Inventor scene-graph.
 * \param aSG An Open-Inventor scene-graph to add the KTE chain to.
 * \param aModel The KTE chain to add to the scene-graph.
 * \return The Open-Inventor scene-graph given as the first parameter, by reference.
 */
oi_scene_graph& operator<<(oi_scene_graph& aSG, const kte::kte_map_chain& aModel);

};

};

#endif










