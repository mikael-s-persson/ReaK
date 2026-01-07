/**
 * \file rate_limited_space_metamaps.h
 *
 * This library provides classes that define meta-maps between rate-limited and normal spaces.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date October 2013
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

#ifndef REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACE_METAMAPS_H_
#define REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACE_METAMAPS_H_


namespace ReaK::pp {

template <typename RateLimitedSpace>
struct get_rate_illimited_space {
  using type = RateLimitedSpace;
};

template <typename RateLimitedSpace>
using get_rate_illimited_space_t =
    typename get_rate_illimited_space<RateLimitedSpace>::type;

template <typename NormalSpace>
struct get_rate_limited_space {
  using type = NormalSpace;
};

template <typename NormalSpace>
using get_rate_limited_space_t =
    typename get_rate_limited_space<NormalSpace>::type;

}  // namespace ReaK::pp

#endif  // REAK_TOPOLOGIES_SPACES_RATE_LIMITED_SPACE_METAMAPS_H_
