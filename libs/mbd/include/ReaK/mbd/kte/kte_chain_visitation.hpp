/**
 * \file kte_chain_visitation.hpp
 *
 * This library provides a visitation template that can be used to perform a filtered visitation in a
 * kte_map_chain object.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2013
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

#ifndef REAK_KTE_CHAIN_VISITATION_HPP
#define REAK_KTE_CHAIN_VISITATION_HPP

#include "kte_map_chain.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/mpl/prior.hpp>

namespace ReaK {

namespace kte {

namespace detail {

template < typename Tuple, typename Idx, typename Visitor >
typename boost::disable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 0 > >, void >::type
  try_visit_kte( Visitor&, kte_map& ){};

template < typename Tuple, typename Idx, typename Visitor >
typename boost::enable_if< boost::mpl::greater< Idx, boost::mpl::size_t< 0 > >, void >::type
  try_visit_kte( Visitor& vis, kte_map& aModel ) {
  typedef typename boost::mpl::prior< Idx >::type IdxPrior;

  try_visit_kte< Tuple, IdxPrior, Visitor >( vis, aModel );

  typedef typename boost::tuples::element< IdxPrior::value, Tuple >::type TestType;
  if( aModel.castTo( TestType::getStaticObjectType() ) )
    vis( static_cast< TestType& >( aModel ) );
};
};

template < typename Tuple, typename Visitor >
void visit_kte_chain( Visitor& vis, const kte_map_chain& aChain );

template < typename Tuple, typename Visitor >
void visit_kte( Visitor& vis, kte_map& aModel ) {

  detail::try_visit_kte< Tuple, boost::mpl::size_t< boost::tuples::length< Tuple >::value >, Visitor >( vis, aModel );

  if( aModel.castTo( kte_map_chain::getStaticObjectType() ) ) {
    visit_kte_chain< Tuple, Visitor >( vis, static_cast< kte_map_chain& >( aModel ) );
  };
};


template < typename Tuple, typename Visitor >
void visit_kte_chain( Visitor& vis, const kte_map_chain& aChain ) {

  const std::vector< shared_ptr< kte_map > >& mdl_list = aChain.getKTEs();
  for( std::size_t i = 0; i < mdl_list.size(); ++i )
    visit_kte< Tuple, Visitor >( vis, *( mdl_list[i] ) );
};
};
};

#endif
