
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

#include <ReaK/core/base/defs.hpp>

#ifndef BOOST_NO_CXX11_EXTERN_TEMPLATE

#include <ReaK/topologies/spaces/Ndof_spaces.hpp>
#include <ReaK/topologies/spaces/joint_space_limits.tpp>


namespace ReaK {

namespace pp {


#define RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( NDOF )                                      \
                                                                                                         \
  template Ndof_rl_space< double, NDOF, 0 >::type joint_limits_mapping< double >::make_rl_joint_space(   \
    const Ndof_space< double, NDOF, 0 >::type& ) const;                                                  \
  template Ndof_rl_space< double, NDOF, 1 >::type joint_limits_mapping< double >::make_rl_joint_space(   \
    const Ndof_space< double, NDOF, 1 >::type& ) const;                                                  \
  template Ndof_rl_space< double, NDOF, 2 >::type joint_limits_mapping< double >::make_rl_joint_space(   \
    const Ndof_space< double, NDOF, 2 >::type& ) const;                                                  \
                                                                                                         \
  template Ndof_space< double, NDOF, 0 >::type joint_limits_mapping< double >::make_normal_joint_space(  \
    const Ndof_rl_space< double, NDOF, 0 >::type& ) const;                                               \
  template Ndof_space< double, NDOF, 1 >::type joint_limits_mapping< double >::make_normal_joint_space(  \
    const Ndof_rl_space< double, NDOF, 1 >::type& ) const;                                               \
  template Ndof_space< double, NDOF, 2 >::type joint_limits_mapping< double >::make_normal_joint_space(  \
    const Ndof_rl_space< double, NDOF, 2 >::type& ) const;                                               \
                                                                                                         \
  template topology_traits< Ndof_rl_space< double, NDOF, 0 >::type >::point_type                         \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_space< double, NDOF, 0 >::type >::point_type& pt,                      \
      const Ndof_space< double, NDOF, 0 >::type&, const Ndof_rl_space< double, NDOF, 0 >::type& ) const; \
  template topology_traits< Ndof_rl_space< double, NDOF, 1 >::type >::point_type                         \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_space< double, NDOF, 1 >::type >::point_type& pt,                      \
      const Ndof_space< double, NDOF, 1 >::type&, const Ndof_rl_space< double, NDOF, 1 >::type& ) const; \
  template topology_traits< Ndof_rl_space< double, NDOF, 2 >::type >::point_type                         \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_space< double, NDOF, 2 >::type >::point_type& pt,                      \
      const Ndof_space< double, NDOF, 2 >::type&, const Ndof_rl_space< double, NDOF, 2 >::type& ) const; \
                                                                                                         \
  template topology_traits< Ndof_space< double, NDOF, 0 >::type >::point_type                            \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_rl_space< double, NDOF, 0 >::type >::point_type& pt,                   \
      const Ndof_rl_space< double, NDOF, 0 >::type&, const Ndof_space< double, NDOF, 0 >::type& ) const; \
  template topology_traits< Ndof_space< double, NDOF, 1 >::type >::point_type                            \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_rl_space< double, NDOF, 1 >::type >::point_type& pt,                   \
      const Ndof_rl_space< double, NDOF, 1 >::type&, const Ndof_space< double, NDOF, 1 >::type& ) const; \
  template topology_traits< Ndof_space< double, NDOF, 2 >::type >::point_type                            \
    joint_limits_mapping< double >::map_to_space(                                                        \
      const topology_traits< Ndof_rl_space< double, NDOF, 2 >::type >::point_type& pt,                   \
      const Ndof_rl_space< double, NDOF, 2 >::type&, const Ndof_space< double, NDOF, 2 >::type& ) const;

RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 0 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 1 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 2 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 3 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 4 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 5 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 6 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 7 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 8 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 9 )
RK_JOINT_SPACE_NDOF_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF( 10 )
};
};

#else

namespace ReaK {

namespace pp {

void dummy_Ndof_spaces_externs_3_symbol(){};
};
};

#endif
