/**
 * \file time_poisson_topology.hpp
 * 
 * This library provides class that define a time-topology with a Poisson distribution. 
 * A time-topology is a simple metric-space where the points are real values (doubles). 
 * However, because time is unlimited, this topology uses a Poisson distribution to 
 * generate random-points at discrete intervals.
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

#ifndef REAK_TIME_POISSON_TOPOLOGY_HPP
#define REAK_TIME_POISSON_TOPOLOGY_HPP


#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include <cmath>

#include "time_topology.hpp"
#include <boost/random/variate_generator.hpp>
#include "path_planning/global_rng.hpp"

#include "default_random_sampler.hpp"

namespace ReaK {

namespace pp {

/**
 * This class implements a time-topology with a Poisson distribution. A time-topology is a 
 * simple metric-space where the points are real values (doubles). 
 * However, because time is unlimited, this topology uses a Poisson distribution to 
 * generate random-points at discrete intervals. This class models the MetricSpaceConcept.
 */
class time_poisson_topology : public time_topology
{
  public:
    typedef double point_type;
    typedef double point_difference_type;
    
    typedef default_distance_metric distance_metric_type;
    typedef default_random_sampler random_sampler_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 1);
    
    double time_step;
    double mean_discrete_time; 
    
    explicit time_poisson_topology(const std::string& aName = "time_poisson_topology",
                                   double aTimeStep = 1.0,
				   double aMeanDiscreteTime = 1.0) : 
				   time_topology(aName),
				   time_step(aTimeStep),
				   mean_discrete_time(aMeanDiscreteTime) { };
    
    /**
     * Generates a random point in the space, uniformly distributed.
     * \note This function actually returns the origin of the space.
     */
    point_type random_point() const {
      boost::variate_generator< global_rng_type&, boost::lognormal_distribution< double > > var_gen(get_global_rng(),boost::lognormal_distribution< double >(2.0, 1.0));
      //boost::variate_generator< global_rng_type&, boost::lognormal_distribution< double > > var_gen(get_global_rng(),boost::lognormal_distribution< double >(std::log(mean_discrete_time) + 1.0, 1.0));
      //boost::variate_generator< global_rng_type&, boost::poisson_distribution< int, double > > var_gen(get_global_rng(),boost::poisson_distribution< int, double >(mean_discrete_time));
      return time_step * int(0.5 * mean_discrete_time * point_type(var_gen()));
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(time_step)
        & RK_SERIAL_SAVE_WITH_NAME(mean_discrete_time);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(time_step)
        & RK_SERIAL_LOAD_WITH_NAME(mean_discrete_time);
    };

    RK_RTTI_MAKE_ABSTRACT_1BASE(time_poisson_topology,0xC240000B,1,"time_poisson_topology",time_topology)

};


};

};

#endif








