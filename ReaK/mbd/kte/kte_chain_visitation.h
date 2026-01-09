/**
 * \file kte_chain_visitation.h
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

#ifndef REAK_MBD_KTE_KTE_CHAIN_VISITATION_H_
#define REAK_MBD_KTE_KTE_CHAIN_VISITATION_H_

#include "ReaK/mbd/kte/kte_map_chain.h"

#include <tuple>
#include <type_traits>

namespace ReaK::kte {

namespace detail {

template <typename Tuple, int Idx, typename Visitor>
void try_visit_kte(Visitor& vis, kte_map& aModel) {
  if constexpr (Idx > 0) {
    try_visit_kte<Tuple, Idx - 1, Visitor>(vis, aModel);

    using TestType = std::tuple_element_t<Idx - 1, Tuple>;
    if (aModel.cast_to(TestType::get_static_object_type())) {
      vis(static_cast<TestType&>(aModel));
    }
  }
}
}  // namespace detail

template <typename Tuple, typename Visitor>
void visit_kte_chain(Visitor& vis, const kte_map_chain& aChain);

template <typename Tuple, typename Visitor>
void visit_kte(Visitor& vis, kte_map& aModel) {

  detail::try_visit_kte<Tuple, std::tuple_size_v<Tuple>, Visitor>(vis, aModel);

  if (aModel.cast_to(kte_map_chain::get_static_object_type())) {
    visit_kte_chain<Tuple, Visitor>(vis, static_cast<kte_map_chain&>(aModel));
  }
}

template <typename Tuple, typename Visitor>
void visit_kte_chain(Visitor& vis, const kte_map_chain& aChain) {

  const std::vector<std::shared_ptr<kte_map>>& mdl_list = aChain.getKTEs();
  for (const auto& i : mdl_list) {
    visit_kte<Tuple, Visitor>(vis, *i);
  }
}

}  // namespace ReaK::kte

#endif  // REAK_MBD_KTE_KTE_CHAIN_VISITATION_H_
