/**
 * \file any_knn_synchro.hpp
 * 
 * This library defines a type-erasure base-class for K-nearest-neighbor synchronization objects. 
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

#ifndef REAK_ANY_KNN_SYNCHRO_HPP
#define REAK_ANY_KNN_SYNCHRO_HPP

#include "base/defs.hpp"
#include "base/shared_object.hpp"

#include <boost/config.hpp>
#include "graph_alg/any_graph.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Path-Planning */
namespace pp {


/**
 * This class can be used as the base for a dynamically polymorphic KNN synchronizer.
 */
class any_knn_synchro : public shared_object {
  public:
    typedef any_knn_synchro self;
    
  protected:
    
    virtual void added_vertex_impl(graph::any_graph::vertex_descriptor, graph::any_graph&) const { };
    virtual void removed_vertex_impl(graph::any_graph::vertex_descriptor, graph::any_graph&) const { };
    
  public:
    
    virtual ~any_knn_synchro() { };
    
    /**
     * Called to notify the synchronizer that a vertex was just added to the graph.
     * \tparam Vertex The vertex-descriptor type for the graph.
     * \tparam Graph The graph structure type.
     * \param v The vertex that was just added to the graph.
     * \param g The graph to which a vertex was just added.
     */
    template <typename Vertex, typename Graph>
    void added_vertex(Vertex v, Graph& g) const {
      type_erased_graph<Graph> teg(&g);
      this->added_vertex_impl(graph::any_graph::vertex_descriptor(boost::any(v)), teg);
    };
    
    /**
     * Called to notify the synchronizer that a vertex is about to be removed from the graph.
     * \tparam Vertex The vertex-descriptor type for the graph.
     * \tparam Graph The graph structure type.
     * \param v The vertex that is about to be removed from the graph.
     * \param g The graph from which a vertex is about to be removed.
     */
    template <typename Vertex, typename Graph>
    void removed_vertex(Vertex v, Graph& g) const {
      type_erased_graph<Graph> teg(&g);
      this->removed_vertex_impl(graph::any_graph::vertex_descriptor(boost::any(v)), teg);
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460017,1,"any_knn_synchro",shared_object)
    
};


/**
 * This class can be used to wrap a generic KNN synchronizer within a dynamically 
 * polymorphic KNN synchronizer. This operates on type-erasure via the ReaK::graph::any_graph class.
 * \tparam Graph The graph type.
 * \tparam KNNSynchro The KNN synchronizer object type to be encapsulated by this type-erasure class.
 */
template <typename Graph, typename KNNSynchro>
class type_erased_knn_synchro : public any_knn_synchro {
  public:
    typedef any_knn_synchro base_type;
    typedef type_erased_knn_synchro<Graph, KNNSynchro> self;
    
  protected:
    
    typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
    
    KNNSynchro synchro;
    
  protected:
    
    virtual void added_vertex_impl(graph::any_graph::vertex_descriptor tev, graph::any_graph& teg) const {
      synchro.added_vertex(boost::any_cast<Vertex>(tev), static_cast<graph::type_erased_graph<Graph>&>(teg).base());
    };
    
    virtual void removed_vertex_impl(graph::any_graph::vertex_descriptor tev, graph::any_graph& teg) const { 
      synchro.removed_vertex(boost::any_cast<Vertex>(tev), static_cast<graph::type_erased_graph<Graph>&>(teg).base());
    };
    
  public:
    
    explicit type_erased_knn_synchro(KNNSynchro aSynchro = KNNSynchro()) : any_knn_synchro(), synchro(aSynchro) { };
    
    virtual ~type_erased_knn_synchro() { };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/
    
    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      base_type::save(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(synchro);
    };
    
    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      base_type::load(A,base_type::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(synchro);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2460018,1,"type_erased_knn_synchro",base_type)
    
};




};

};


#endif


