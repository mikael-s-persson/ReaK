
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

#include "base/defs.hpp"

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "topologies/joint_space_topologies.hpp"
#include "topologies/joint_space_limits.tpp"


namespace ReaK {

namespace pp {



#define RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(NDOF) \
\
template Ndof_0th_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_0th_order_space<double,NDOF>::type&) const;\
template Ndof_1st_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_1st_order_space<double,NDOF>::type&) const;\
template Ndof_2nd_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_2nd_order_space<double,NDOF>::type&) const;\
\
template Ndof_0th_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_0th_order_rl_space<double,NDOF>::type&) const;\
template Ndof_1st_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_1st_order_rl_space<double,NDOF>::type&) const;\
template Ndof_2nd_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_2nd_order_rl_space<double,NDOF>::type&) const;\
\
template \
topology_traits< Ndof_0th_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_0th_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_0th_order_space<double,NDOF>::type& , \
                                                const Ndof_0th_order_rl_space<double,NDOF>::type& ) const;\
template \
topology_traits< Ndof_1st_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_1st_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_1st_order_space<double,NDOF>::type& , \
                                                const Ndof_1st_order_rl_space<double,NDOF>::type& ) const;\
template \
topology_traits< Ndof_2nd_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_2nd_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_2nd_order_space<double,NDOF>::type& , \
                                                const Ndof_2nd_order_rl_space<double,NDOF>::type& ) const;\
\
template \
topology_traits< Ndof_0th_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_0th_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_0th_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_0th_order_space<double,NDOF>::type& ) const;\
template \
topology_traits< Ndof_1st_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_1st_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_1st_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_1st_order_space<double,NDOF>::type& ) const;\
template \
topology_traits< Ndof_2nd_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_2nd_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_2nd_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_2nd_order_space<double,NDOF>::type& ) const;

RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(1)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(2)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(3)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(4)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(5)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(6)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(7)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(8)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(9)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC_DEF(10)





};

};

#else

namespace ReaK {

namespace pp {

void dummy_joint_space_topologies_externs_5_symbol() { };

};

};

#endif














