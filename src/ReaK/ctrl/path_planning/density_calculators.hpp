/**
 * \file density_calculators.hpp
 * 
 * This library defines 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
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

#ifndef REAK_DENSITY_CALCULATORS_HPP
#define REAK_DENSITY_CALCULATORS_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include "planning_visitors.hpp"

namespace ReaK {
  
namespace pp {

  
struct prm_density_calculator {
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_explored(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_succeeded(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_failed(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void update_density(Vertex u, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const {
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    using std::exp;
    std::size_t deg_u = out_degree(u,g);
    if(deg_u == 0) {
      g[u].density = 0.0;
      return;
    };
    std::size_t max_node_degree = space_dim + 1;
    double sum = 0.0;
    OutEdgeIter ei, ei_end;
    for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei)
      sum += g[*ei].weight / sampling_radius;
    sum /= double(deg_u) * double(deg_u) / double(max_node_degree);
    g[u].density = exp( -sum * sum );
    
  };
};



struct sbastar_density_calculator {
  double compute_sample_similarity(double travel_dist, double event_radius, double sampling_radius, std::size_t space_dim) const {
    using std::exp; using std::log;
    double sig2_n = sampling_radius * sampling_radius;
    double sig2_x = event_radius * event_radius;
    return exp(-travel_dist * travel_dist / (sig2_x * 2.0) - 0.5 * double(space_dim) * ( sig2_n / sig2_x - 1.0 - log(sig2_n / sig2_x) ) );
  };
  
  double compute_sample_similarity(double travel_dist, double sampling_radius) const {
    using std::exp;
    return exp(-travel_dist * travel_dist / (sampling_radius * sampling_radius * 2.0));
  };

#if 0
  // keep the average sample similarity:
  void register_sample(double& density, std::size_t& count, double samp_sim, std::size_t space_dim) const {
    density = (count * density + samp_sim) / (count + 1);
    ++count;
  };
#endif
#if 0
  // keep the average sample similarity weighted by the sample probability (and its binomial converse).
  void register_sample(double& density, std::size_t& count, double samp_sim, std::size_t space_dim) const {
    double tmp_density = density * (1.0 - samp_sim) + samp_sim * samp_sim;
    density = (count * density + tmp_density) / (count + 1);
    ++count;
  };
#endif
#if 0
  // keep track only of the maximum sample similarity in the region of a node:
  void register_sample(double& density, std::size_t& count, double samp_sim, std::size_t space_dim) const {
    if(samp_sim > density)
      density = samp_sim;
    ++count;
  };
#endif
#if 1
  // keep the sample similarity weighted by the sample probability (and its binomial converse).
  // that is, assume the existing density to reflect the overall density and the newly computed 
  // sample similarity to reflect the density in its relatively probable region (binomial).
  void register_sample(double& density, std::size_t& count, double samp_sim, std::size_t space_dim) const {
    using std::pow;
    double correct_factor = 1.0 / pow(2.0 * M_PI, (double(space_dim) - 1.0) / 2.0);
    density = density * (1.0 - correct_factor * samp_sim) + correct_factor * samp_sim * samp_sim;
    ++count;
  };
#endif
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_explored(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { 
    double dist = get(distance_metric, space.get_super_space())(g[u].position, g[v].position, space.get_super_space());
    double samp_sim = compute_sample_similarity(dist, sampling_radius);
    register_sample(g[u].density, g[u].expansion_trials, samp_sim, space_dim);
    register_sample(g[v].density, g[v].expansion_trials, samp_sim, space_dim);
  };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_succeeded(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void travel_failed(Vertex u, Vertex v, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const { 
    double dist = get(distance_metric, space.get_super_space())(g[u].position, g[v].position, space.get_super_space());
    double samp_sim = compute_sample_similarity(0.5 * dist, 0.25 * dist, sampling_radius, space_dim);
    register_sample(g[u].constriction, g[u].collision_count, samp_sim, space_dim);
    register_sample(g[v].constriction, g[v].collision_count, samp_sim, space_dim);
  };
  
  template <typename Vertex, typename Graph, typename SpaceType>
  void update_density(Vertex u, Graph& g, const SpaceType& space, double sampling_radius, std::size_t space_dim) const {
    g[u].density = 0.0;
    g[u].expansion_trials = 0;
    typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
    OutEdgeIter ei, ei_end;
    for(boost::tie(ei,ei_end) = out_edges(u,g); ei != ei_end; ++ei) {
      double samp_sim = compute_sample_similarity(g[*ei].weight, sampling_radius);
      register_sample(g[u].density, g[u].expansion_trials, samp_sim, space_dim);
    };
  };
};

};

};

#endif

