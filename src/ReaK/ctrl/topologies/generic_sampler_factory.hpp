/**
 * \file generic_sampler_factory.hpp
 * 
 * This library provides an implementation of a generic sampler factory.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_GENERIC_SAMPLER_FACTORY_HPP
#define REAK_GENERIC_SAMPLER_FACTORY_HPP

#include "path_planning/tangent_bundle_concept.hpp"
#include "path_planning/temporal_space_concept.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "topologies/metric_space_tuple.hpp"

#include <boost/config.hpp>
#include <boost/concept_check.hpp>
#include <cmath>

namespace ReaK {

namespace pp {
  
namespace detail {

  

template < typename Sampler,
           typename SpaceType>
class generic_sampler_impl { 
  public:
    typedef generic_sampler_impl< Sampler, SpaceType > self;
    
    typedef typename topology_traits<SpaceType>::point_type point_type;
    
  public:
    
    template <typename Factory>
    static void generate_sample(const Sampler& sampler, point_type& result, const SpaceType& space, const Factory& factory) {
      result = sampler(space,factory);
    };
    
    static void generate_sample(const Sampler& sampler, point_type& result, const SpaceType& space) {
      result = sampler(space);
    };
};


  
  template <std::size_t Size, typename SpaceTuple, typename Sampler >
  struct generic_sampler_impl_tuple_impl { 
    //BOOST_STATIC_ASSERT(false);
  };
  
#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES
  
  template <std::size_t Size, typename Sampler, typename... Spaces>
  struct generic_sampler_impl_tuple_impl< Size, std::tuple<Spaces...>, Sampler > {
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, Spaces >... > type;
  };
  
  template <std::size_t Size, typename Sampler, typename... Spaces>
  struct generic_sampler_impl_tuple_impl< Size, arithmetic_tuple<Spaces...>, Sampler > {
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, Spaces >... > type;
  };
  
#else
  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 1, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 2, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 3, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 4, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 5, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 6, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<5,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 7, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<5,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<6,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 8, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<5,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<6,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<7,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 9, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<5,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<6,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<7,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<8,SpaceTuple>::type > > type;
  };

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple_impl< 10, SpaceTuple, Sampler > { 
    typedef arithmetic_tuple< generic_sampler_impl< Sampler, typename arithmetic_tuple_element<0,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<1,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<2,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<3,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<4,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<5,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<6,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<7,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<8,SpaceTuple>::type >,
                              generic_sampler_impl< Sampler, typename arithmetic_tuple_element<9,SpaceTuple>::type > > type;
  };

#endif

  template <typename SpaceTuple, typename Sampler>
  struct generic_sampler_impl_tuple : generic_sampler_impl_tuple_impl< arithmetic_tuple_size<SpaceTuple>::type::value, SpaceTuple, Sampler > { };
  


template <std::size_t Idx, typename SamplerTuple>
class gen_sampler_recursion_impl {
  public:
    
    template <typename Sampler, typename PointType, typename SpaceTuple, typename Factory>
    static void generate_sample(Sampler& sampler, PointType& result, const SpaceTuple& space, const Factory& factory) {
      gen_sampler_recursion_impl< Idx-1, SamplerTuple >::generate_sample(sampler, result, space, factory);
      
      arithmetic_tuple_element<Idx,SamplerTuple>::type::generate_sample(sampler, get<Idx>(result), get<Idx>(space), factory);
    };
    
    template <typename Sampler, typename PointType, typename SpaceTuple>
    static void generate_sample(Sampler& sampler, PointType& result, const SpaceTuple& space) {
      gen_sampler_recursion_impl< Idx-1, SamplerTuple >::generate_sample(sampler, result, space);
      
      arithmetic_tuple_element<Idx,SamplerTuple>::type::generate_sample(sampler, get<Idx>(result), get<Idx>(space));
    };
};


template <typename SamplerTuple>
class gen_sampler_recursion_impl<0, SamplerTuple> {
  public:
    
    template <typename Sampler, typename PointType, typename SpaceTuple, typename Factory>
    static void generate_sample(const Sampler& sampler, PointType& result, const SpaceTuple& space, const Factory& factory) {
      arithmetic_tuple_element<0,SamplerTuple>::type::generate_sample(sampler, get<0>(result), get<0>(space), factory);
    };
    
    template <typename Sampler, typename PointType, typename SpaceTuple>
    static void generate_sample(const Sampler& sampler, PointType& result, const SpaceTuple& space) {
      arithmetic_tuple_element<0,SamplerTuple>::type::generate_sample(sampler, get<0>(result), get<0>(space));
    };
    
};



template < typename Sampler,
           typename SpaceTuple, typename TupleDistMetric>
class generic_sampler_impl< Sampler, metric_space_tuple<SpaceTuple,TupleDistMetric> > {
  public:
    typedef generic_sampler_impl< Sampler, metric_space_tuple<SpaceTuple,TupleDistMetric> > self;
    
    typedef metric_space_tuple<SpaceTuple,TupleDistMetric> SpaceType;
    typedef typename topology_traits<SpaceType>::point_type point_type;
    
  public:
    
    template <typename Factory>
    static void generate_sample(const Sampler& sampler, point_type& result, const SpaceType& space, const Factory& factory) {
      gen_sampler_recursion_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, typename generic_sampler_impl_tuple< SpaceTuple, Sampler >::type >::generate_sample(sampler, result, space, factory);
    };
    
    static void generate_sample(const Sampler& sampler, point_type& result, const SpaceType& space) {
      gen_sampler_recursion_impl< arithmetic_tuple_size< SpaceTuple >::value - 1, typename generic_sampler_impl_tuple< SpaceTuple, Sampler >::type >::generate_sample(sampler, result, space);
    };
};


};


/**
 * This functor class implements a generic sampler in a topology.
 * \tparam Factory The sampler factory type which created this generic sampler (this is a fly-weight class to store whatever global parameters are needed by the sampler implementation).
 * \tparam Sampler The sampler implementation which is used to perform the specifics of the actual sampling.
 */
template <typename Sampler, typename SpaceType, typename Factory = void>
class generic_sampler {
  public:
    typedef generic_sampler<Sampler, SpaceType, Factory> self;
    typedef typename topology_traits<SpaceType>::point_type point_type;
    typedef SpaceType topology;
    
    typedef detail::generic_sampler_impl< Sampler, topology> sampler_impl_type;
    
  private:
    const Factory* parent;
    Sampler sampler;
    
  public:
    
    /**
     * Default constructor.
     */
    generic_sampler(Sampler aSampler = Sampler(), const Factory* aParent = NULL) : parent(aParent), sampler(aSampler) { };
    
    template <typename OtherFactory>
    point_type operator()(const topology& space, const OtherFactory& aNewParent) const {
      point_type result;
      sampler_impl_type::generate_sample(sampler, result, space, aNewParent);
      return result;
    };
    
    point_type operator()(const topology& space) const {
      point_type result;
      sampler_impl_type::generate_sample(sampler, result, space, *parent);
      return result;
    };
    
};

template <typename Sampler, typename SpaceType>
class generic_sampler<Sampler, SpaceType, void> {
  public:
    typedef generic_sampler<Sampler, SpaceType, void> self;
    typedef typename topology_traits<SpaceType>::point_type point_type;
    typedef SpaceType topology;
    
    typedef detail::generic_sampler_impl< Sampler, topology> sampler_impl_type;
    
  private:
    Sampler sampler;
    
  public:
    
    /**
     * Default constructor.
     */
    generic_sampler(Sampler aSampler = Sampler()) : sampler(aSampler) { };
    
    template <typename OtherFactory>
    point_type operator()(const topology& space, const OtherFactory& aParent) const {
      point_type result;
      sampler_impl_type::generate_sample(sampler, result, space, aParent);
      return result;
    };
    
    point_type operator()(const topology& space) const {
      point_type result;
      sampler_impl_type::generate_sample(sampler, result, space);
      return result;
    };
    
};



};

};

#endif









