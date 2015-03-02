//  (C) Copyright Cryolite 2014. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef REAK_INDEX_SEQUENCE_HPP
#define REAK_INDEX_SEQUENCE_HPP

#include "defs.hpp"
#include <cstddef>

#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES

namespace ReaK {

template <std::size_t... Is>
struct index_sequence {};

namespace detail {

template <typename IndexSeq1, typename IndexSeq2>
struct concat_index_sequence;

template <std::size_t... Idx1, std::size_t... Idx2>
struct concat_index_sequence< index_sequence<Idx1...>, index_sequence<Idx2...> > {
  typedef index_sequence<Idx1..., Idx2 + sizeof...(Idx1)...> type;
};

};

template <std::size_t N>
struct make_index_sequence;

template <>
struct make_index_sequence<0> {
  typedef index_sequence<> type;
};

template <>
struct make_index_sequence<1> {
  typedef index_sequence<0> type;
};

template<std::size_t N>
struct make_index_sequence
  : detail::concat_index_sequence<
      typename make_index_sequence< (N + 1) / 2 >::type,
      typename make_index_sequence< N / 2 >::type >
{};

}; // ReaK

#endif

#endif


