/**
 * \file tree_organizer_concept.hpp
 *
 * This library provides a concept class that defines a tree-organizing visitor. Such a visitor is meant
 * to be used to lump insertion / deletion callbacks that maintain a tree structure.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2012
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

#ifndef REAK_TREE_ORGANIZER_CONCEPT_HPP
#define REAK_TREE_ORGANIZER_CONCEPT_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/config.hpp>
#include <boost/concept_check.hpp>

namespace ReaK {

namespace graph {


/**
 * This concept defines the requirements to fulfill in order to model a tree-organizing visitor concept
 * as used in ReaK::graph.
 *
 * Valid expressions (vertex v; vertex_iterator cv_it, cv_it_end; TreeType tree):
 *
 * vis.remove_vertex(v, tree);  The visitor can perform the removal of the vertex (v) from the tree.
 *
 * v = vis.add_vertex(tree_vp, tree);  The visitor can perform the addition a vertex (v) with property (vp) to the tree.
 *
 * C++0x / C++11 only:
 *
 * v = vis.add_vertex(std::move(tree_vp), tree);  The visitor can perform the addition a vertex (v) by moving in the
 *property (vp) to the tree.
 *
 * \tparam TreeOrganizerVisitor The visitor type to be checked for this concept.
 * \tparam TreeType The tree type on which the visitor must operate.
 */
template < typename TreeOrganizerVisitor, typename TreeType >
struct TreeOrganizerVisitorConcept {
  TreeType tree;
  TreeOrganizerVisitor vis;
  typename TreeType::vertex_property_type vp;
  typename boost::graph_traits< TreeType >::vertex_descriptor v;

  BOOST_CONCEPT_ASSERT( ( boost::IncidenceGraphConcept< TreeType > ) );

  BOOST_CONCEPT_USAGE( TreeOrganizerVisitorConcept ) {
    vis.remove_vertex( v, tree );
    vis.add_vertex( vp, tree );
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vis.add_vertex( std::move( vp ), tree );
#endif
  };
};
};
};


#endif
