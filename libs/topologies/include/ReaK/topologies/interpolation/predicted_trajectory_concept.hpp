/**
 * \file predicted_trajectory_concept.hpp
 *
 * This library defines the concept that represents a predicted trajectory. This concept
 * supplements the SpatialTrajectoryConcept with a few additional requirements which
 * characterize a predicted trajectory (see PredictedTrajectoryConcept).
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2011
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

#ifndef REAK_PREDICTED_TRAJECTORY_CONCEPT_HPP
#define REAK_PREDICTED_TRAJECTORY_CONCEPT_HPP

#include <ReaK/core/base/defs.hpp>

#include "spatial_path_concept.hpp"
#include <ReaK/topologies/spaces/temporal_space_concept.hpp>
#include "spatial_trajectory_concept.hpp"

#include <boost/concept_check.hpp>

namespace ReaK {

namespace pp {


/**
 * This concept represents a predicted trajectory, as used in ReaK::pp. This concept
 * supplements the SpatialTrajectoryConcept with a few additional requirements which
 * characterize a predicted trajectory (see PredictedTrajectoryConcept).
 *
 * Required concepts:
 *
 * The PredictedTrajectory should model the SpatialTrajectoryConcept on the given Topology.
 *
 * The Topology should model the TemporalSpaceConcept.
 *
 * Valid expressions:
 *
 * std::pair< const_waypoint_descriptor, point_type> w_p;  A waypoint is a pair of a const-waypoint descriptor and a
 *point on the topology.
 *
 * p.set_initial_point(pt);  The initial temporal point (pt) can be set to seed the predicted trajectory (p).
 *
 * p.set_initial_point(w_p);  The initial temporal waypoint (w_p) can be set to seed the predicted trajectory (p).
 *
 * \tparam PredictedTrajectory The trajectory type to be checked for compliance to this concept.
 * \tparam Topology The temporal-topology type that can contain the trajectory.
 */
template < typename PredictedTrajectory, typename Topology >
struct PredictedTrajectoryConcept : public SpatialTrajectoryConcept< PredictedTrajectory, Topology > {

  BOOST_CONCEPT_USAGE( PredictedTrajectoryConcept ) {
    this->p.set_initial_point( this->pt );
    this->p.set_initial_point( this->w_p );
  };
};
};
};


#endif
