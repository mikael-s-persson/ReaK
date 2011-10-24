/**
 * \file hyperbox_topology.hpp
 * 
 * This library provides classes that define a hyper-box vector-topology. A hyper-box vector-topology is 
 * a vector-topology where the points are vector values and the boundary is a hyper-box.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2011
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

#ifndef REAK_HYPERBOX_TOPOLOGY_HPP
#define REAK_HYPERBOX_TOPOLOGY_HPP


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "lin_alg/vect_concepts.hpp"

#include "vector_topology.hpp"

#include <cmath>
#include "base/named_object.hpp"

#include "path_planning/global_rng.hpp"

namespace ReaK {

namespace pp {

/**
 * This library provides classes that define a hyper-box vector-topology. A hyper-box vector-topology is 
 * a vector-topology where the points are vector values and the boundary is a hyper-box.
 * \tparam Vector The vector-type for the topology, should model an Arithmetic concept and WritableVectorConcept.
 */
template <typename Vector>
class hyperbox_topology : public vector_topology<Vector>
{
  public:
    typedef hyperbox_topology<Vector> self;
    
    typedef Vector point_type;
    typedef Vector point_difference_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = vect_traits<Vector>::dimensions);
    
  protected:
    point_type lower_corner;
    point_type upper_corner;
    
  public:
    
    hyperbox_topology(const std::string& aName = "hyperbox_topology",
                      const point_type& aLowerCorner = point_type(),
		      const point_type& aUpperCorner = point_type()) : 
		      vector_topology<Vector>(aName),
		      lower_corner(aLowerCorner),
		      upper_corner(aUpperCorner) { };
		       
		       
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      using std::sqrt;
      point_difference_type dp = this->difference(b,a);
      return sqrt(dp * dp);
    }
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      using std::sqrt;
      return sqrt(delta * delta);
    }
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      boost::uniform_01<pp::global_rng_type&,double> var_rnd(pp::get_global_rng());
      
      point_type p = lower_corner;
      for(typename vect_traits<point_type>::size_type i = 0; i < p.size(); ++i)
        p[i] += var_rnd() * (upper_corner[i] - lower_corner[i]);
      
      return p;
    };

    /**
     * Takes a point and clips it to within this line-segment space.
     */
    point_type bound(point_type a) const {
      for(typename vect_traits<point_type>::size_type i = 0; i < a.size(); ++i) {
	if(lower_corner[i] < upper_corner[i]) {
	  if(a[i] < lower_corner[i])
	    a[i] = lower_corner[i];
	  else if(a[i] > upper_corner[i])
	    a[i] = upper_corner[i];
	} else {
	  if(a[i] > lower_corner[i])
	    a[i] = lower_corner[i];
	  else if(a[i] < upper_corner[i])
	    a[i] = upper_corner[i];
	};
      };
      return a;
    };

    /**
     * Returns the distance to the boundary of the space.
     */
    double distance_from_boundary(const point_type& a) const {
      double dist = std::numeric_limits< typename vect_traits<point_type>::value_type >::max();
      using std::fabs;
      for(typename vect_traits<point_type>::size_type i = 0; i < a.size(); ++i) {
	if(dist > fabs(a[i] - lower_corner[i]))
	  dist = fabs(a[i] - lower_corner[i]);
	if(dist > fabs(a[i] - upper_corner[i]))
	  dist = fabs(a[i] - upper_corner[i]);
      };
      return dist;
    };
    
    /**
     * Returns the difference to the closest boundary.
     */
    point_difference_type get_diff_to_boundary(const point_type& a) const {
      double dist = std::numeric_limits< typename vect_traits<point_type>::value_type >::max();
      typename vect_traits<point_type>::size_type j = a.size();
      bool at_upper = false;
      using std::fabs;
      for(typename vect_traits<point_type>::size_type i = 0; i < a.size(); ++i) {
	if(dist > fabs(a[i] - lower_corner[i])) {
	  j = i;
	  at_upper = false;
	  dist = fabs(a[i] - lower_corner[i]);
	};
	if(dist > fabs(a[i] - upper_corner[i])) {
	  j = i;
	  at_upper = true;
	  dist = fabs(a[i] - upper_corner[i]);
	};
      };
      point_difference_type dp = this->difference(a,a);
      if(j == a.size())
	return dp;
      if(at_upper)
	dp[j] = upper_corner[j] - a[j];
      else
	dp[j] = lower_corner[j] - a[j];
      return dp;
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return this->adjust(lower_corner,this->difference(upper_corner,lower_corner));
    };

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(lower_corner)
        & RK_SERIAL_SAVE_WITH_NAME(upper_corner);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(lower_corner)
        & RK_SERIAL_LOAD_WITH_NAME(upper_corner);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(self,0xC2400009,1,"hyperbox_topology",vector_topology<Vector>)
    
};


};

};

#endif








