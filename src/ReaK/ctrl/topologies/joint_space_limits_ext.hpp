/**
 * \file joint_space_limits_ext.hpp
 * 
 * This library provides extern template declarations for classes defined in the joint_space_limits.hpp header file.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date November 2012
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

#ifndef REAK_JOINT_SPACE_LIMITS_EXT_HPP
#define REAK_JOINT_SPACE_LIMITS_EXT_HPP

#if (defined(RK_ENABLE_CXX11_FEATURES) && defined(RK_ENABLE_EXTERN_TEMPLATES))

#include "base/defs.hpp"
#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "differentiable_space.hpp"
#include "metric_space_tuple.hpp"
#include "rate_limited_spaces.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "kinetostatics/gen_coord.hpp"

namespace ReaK {

namespace pp {
  

extern template struct joint_limits_collection<double>;


#define RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(NDOF) \
\
extern template Ndof_0th_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_0th_order_space<double,NDOF>::type&) const;\
extern template Ndof_1st_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_1st_order_space<double,NDOF>::type&) const;\
extern template Ndof_2nd_order_rl_space<double,NDOF>::type joint_limits_collection<double>::make_rl_joint_space(const Ndof_2nd_order_space<double,NDOF>::type&) const;\
\
extern template Ndof_0th_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_0th_order_rl_space<double,NDOF>::type&) const;\
extern template Ndof_1st_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_1st_order_rl_space<double,NDOF>::type&) const;\
extern template Ndof_2nd_order_space<double,NDOF>::type joint_limits_collection<double>::make_normal_joint_space(const Ndof_2nd_order_rl_space<double,NDOF>::type&) const;\
\
extern template \
topology_traits< Ndof_0th_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_0th_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_0th_order_space<double,NDOF>::type& , \
                                                const Ndof_0th_order_rl_space<double,NDOF>::type& ) const;\
extern template \
topology_traits< Ndof_1st_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_1st_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_1st_order_space<double,NDOF>::type& , \
                                                const Ndof_1st_order_rl_space<double,NDOF>::type& ) const;\
extern template \
topology_traits< Ndof_2nd_order_rl_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_2nd_order_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_2nd_order_space<double,NDOF>::type& , \
                                                const Ndof_2nd_order_rl_space<double,NDOF>::type& ) const;\
\
extern template \
topology_traits< Ndof_0th_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_0th_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_0th_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_0th_order_space<double,NDOF>::type& ) const;\
extern template \
topology_traits< Ndof_1st_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_1st_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_1st_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_1st_order_space<double,NDOF>::type& ) const;\
extern template \
topology_traits< Ndof_2nd_order_space<double,NDOF>::type >::point_type \
  joint_limits_collection<double>::map_to_space(const topology_traits< Ndof_2nd_order_rl_space<double,NDOF>::type >::point_type& pt,\
                                                const Ndof_2nd_order_rl_space<double,NDOF>::type& , \
                                                const Ndof_2nd_order_space<double,NDOF>::type& ) const;

RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(1)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(2)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(3)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(4)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(5)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(6)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(7)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(8)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(9)
RK_JOINT_SPACE_LIMITS_EXT_MAKE_MEMBER_FUNC(10)



extern template metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_0th_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_1st_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type&) const;

extern template metric_space_array< se2_0th_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_1st_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se2_2nd_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type&) const;

extern template 
topology_traits< metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_0th_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_0th_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_1st_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_1st_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_2nd_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type& ) const;

extern template 
topology_traits< metric_space_array< se2_0th_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_0th_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_0th_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_1st_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_1st_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_1st_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se2_2nd_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se2_2nd_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se2_2nd_order_topology<double>::type, 1>::type& ) const;



extern template metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_0th_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_1st_order_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type joint_limits_collection<double>::make_rl_joint_space(const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type&) const;

extern template metric_space_array< se3_0th_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se3_1st_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type&) const;
extern template metric_space_array< se3_2nd_order_topology<double>::type, 1>::type joint_limits_collection<double>::make_normal_joint_space(const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type&) const;

extern template 
topology_traits< metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_0th_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_0th_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_1st_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_1st_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_2nd_order_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type& ) const;

extern template 
topology_traits< metric_space_array< se3_0th_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_0th_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_0th_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se3_1st_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_1st_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_1st_order_topology<double>::type, 1>::type& ) const;
extern template 
topology_traits< metric_space_array< se3_2nd_order_topology<double>::type, 1>::type >::point_type 
  joint_limits_collection<double>::map_to_space(const topology_traits< metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type >::point_type& pt,
                                                const metric_space_array< se3_2nd_order_rl_topology<double>::type, 1>::type& , 
                                                const metric_space_array< se3_2nd_order_topology<double>::type, 1>::type& ) const;


};

};

#endif


#endif








