/**
 * \file hyperball_topology.hpp
 * 
 * This library provides classes that define a hyper-ball vector-topology. A hyper-ball vector-topology is 
 * a vector-topology where the points are vector values and the boundary is a hyper-ellipsoid.
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

#ifndef REAK_HYPERBALL_TOPOLOGY_HPP
#define REAK_HYPERBALL_TOPOLOGY_HPP


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "lin_alg/vect_concepts.hpp"

#include "vector_topology.hpp"
#include "path_planning/metric_space_concept.hpp"
#include "default_random_sampler.hpp"

#include <cmath>
#include "base/named_object.hpp"

#include "base/global_rng.hpp"

namespace ReaK {

namespace pp {

/**
 * This library provides classes that define a hyper-ball vector-topology. A hyper-ball vector-topology is 
 * a vector-topology where the points are vector values and the boundary is a hyper-ellipsoid.
 * This class models the MetricSpaceConcept, the LieGroupConcept, the BoundedSpaceConcept, 
 * the SphereBoundedSpaceConcept, and the PointDistributionConcept.
 * \tparam Vector The vector-type for the topology, should model an Arithmetic concept and WritableVectorConcept.
 */
template <typename Vector >
class hyperball_topology : public vector_topology<Vector>
{
  public:
    typedef hyperball_topology<Vector> self;
    
    typedef Vector point_type;
    typedef Vector point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = vect_traits<Vector>::dimensions);
    
  protected:
    point_type center_point;
    double radius_value;
    
  public:
    
    hyperball_topology(const std::string& aName = "hyperball_topology",
                       const point_type& aOrigin = point_type(),
                       double aRadius = 1.0) : 
                       vector_topology<Vector>(aName),
                       center_point(aOrigin),
                       radius_value(aRadius) { };
    
   /*************************************************************************
    *                             MetricSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the distance between two points.
     */
    double distance(const point_type& a, const point_type& b) const 
    {
      return this->norm(this->difference(b,a));
    }
    
    /**
     * Returns the norm of the difference between two points.
     */
    double norm(const point_difference_type& delta) const {
      using std::sqrt;
      double result = sqrt(delta * delta);
      return result;
    }
    
   /*************************************************************************
    *                         for PointDistributionConcept
    * **********************************************************************/
    
    /**
     * Generates a random point in the space, uniformly distributed.
     */
    point_type random_point() const {
      using std::sqrt;
      
      point_difference_type dp = this->difference(center_point,center_point);
      if(dp.size() == 0)
        return center_point;
      
      double radial_dim_correction = double(dp.size());
      
      boost::variate_generator< global_rng_type&, boost::normal_distribution<typename vect_traits<point_difference_type>::value_type> > var_rnd(get_global_rng(), boost::normal_distribution<typename vect_traits<point_difference_type>::value_type>());
      
      for(typename vect_traits<point_difference_type>::size_type i = 0; i < dp.size(); ++i)
        dp[i] = var_rnd();
      
      double factor = std::pow(boost::uniform_01<global_rng_type&,double>(get_global_rng())(),1.0 / radial_dim_correction) * radius_value / sqrt(dp * dp);
      
      return this->adjust(center_point, factor * dp );
    };
    
   /*************************************************************************
    *                             BoundedSpaceConcept
    * **********************************************************************/

    /**
     * Takes a point and clips it to within this hyperball space.
     */
    void bring_point_in_bounds(point_type& a) const {
      a = this->adjust(a,this->get_diff_to_boundary(a));
    };

    /**
     * Returns the distance to the boundary of the space.
     */
    double distance_from_boundary(const point_type& a) const {
      using std::fabs;
      point_difference_type c2a = this->difference(a,center_point);
      return fabs(radius_value - this->norm(c2a));
    };
    
    /**
     * Returns the difference to the closest boundary.
     */
    point_difference_type get_diff_to_boundary(const point_type& a) const {
      point_difference_type c2a = this->difference(a,center_point);
      return (radius_value - this->norm(c2a)) * c2a;
    };
      
    /**
     * Tests if a given point is within the boundary of this space.
     */
    bool is_in_bounds(const point_type& a) const {
      point_difference_type c2a = this->difference(a,center_point);
      return radius_value >= this->norm(c2a);
    };

    /**
     * Returns the origin of the space (the lower-limit).
     */
    point_type origin() const {
      return center_point;
    };
    
   /*************************************************************************
    *                             SphereBoundedSpaceConcept
    * **********************************************************************/
    
    /**
     * Returns the radius of the space.
     */
    double get_radius() const {
      return radius_value;
    };

    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(center_point)
        & RK_SERIAL_SAVE_WITH_NAME(radius_value);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(center_point)
        & RK_SERIAL_LOAD_WITH_NAME(radius_value);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400008,1,"hyperball_topology",vector_topology<Vector>)
    
};

template <typename Vector>
struct is_metric_space< hyperball_topology<Vector> > : boost::mpl::true_ { };

template <typename Vector>
struct is_reversible_space< hyperball_topology<Vector> > : boost::mpl::true_ { };

template <typename Vector>
struct is_point_distribution< hyperball_topology<Vector> > : boost::mpl::true_ { };


};

};



#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include "lin_alg/vect_alg.hpp"

namespace ReaK {

namespace pp {

extern template class hyperball_topology< vect<double,2> >;
extern template class hyperball_topology< vect<double,3> >;
extern template class hyperball_topology< vect<double,4> >;
extern template class hyperball_topology< vect<double,6> >;
extern template class hyperball_topology< vect_n<double> >;


};

};

#endif


#endif








