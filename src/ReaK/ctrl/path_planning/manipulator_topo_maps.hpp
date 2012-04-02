/**
 * \file manipulator_topo_maps.hpp
 * 
 * This library provides classes that define topological mappings between a joint-space (generalized 
 * coordinates) and the end-effector frame of a serial manipulator. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date March 2012
 */

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

#ifndef REAK_MANIPULATOR_TOPO_MAPS_HPP
#define REAK_MANIPULATOR_TOPO_MAPS_HPP


#include "base/defs.hpp"

#include <boost/config.hpp> // For BOOST_STATIC_CONSTANT

#include "topologies/joint_space_topologies.hpp"
#include "topologies/joint_space_limits.hpp"
#include "topologies/se3_topologies.hpp"
#include "topologies/se2_topologies.hpp"

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"

#include "kinetostatics/gen_coord.hpp"

#include "mbd_kte/manipulator_model.hpp"
#include "mbd_kte/manipulator_model_helper.hpp"

#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>

namespace ReaK {

namespace pp {


namespace detail {
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type write_one_joint_coord_impl( const PointType& pt,
					   const InSpace&,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
				           const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->Coords()[gen_i++]) = get_gen_coord(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type write_one_joint_coord_impl( const PointType& pt,
					   const InSpace&,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
				           const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->Frames2D()[f2d_i++]) = get_frame_2D(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type write_one_joint_coord_impl( const PointType& pt,
					   const InSpace&,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
				           const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->Frames3D()[f3d_i++]) = get_frame_3D(pt);
  };
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const PointType& pt,
					     const InSpaceTuple& space_in,
					     std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				             const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const PointType& pt,
					     const InSpaceTuple& space_in,
					     std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				             const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>
    >,  
  void >::type write_one_joint_coord_impl( const PointType& pt,
					   const InSpace& space_in,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				           const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size<InSpace> >::type >(pt, space_in, gen_i, f2d_i, f3d_i, model);
  };
  
  
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const PointType& pt,
					     const InSpaceTuple& space_in,
					     std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl< typename boost::mpl::prior<Idx>::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
    write_one_joint_coord_impl(get<Idx::type::value>(pt),get<Idx::type::value>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const PointType& pt,
					     const InSpaceTuple& space_in,
					     std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_one_joint_coord_impl(get<0>(pt),get<0>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type write_joint_coordinates_impl( const PointType& pt,
				     const InSpaceTuple& space_in,
				     const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size<InSpaceTuple> >::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type write_joint_coordinates_impl( const PointType& pt,
				     const InSpaceTuple& space_in,
				     const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_joint_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
					  const InSpace&,
					  std::size_t& gen_i, std::size_t&, std::size_t&,
				          const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_gen_coord(pt,*(model->Coords()[gen_i++]));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
					  const InSpace&,
					  std::size_t&, std::size_t& f2d_i, std::size_t&,
				          const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_frame_2D(pt,*(model->Frames2D()[f2d_i++]));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
					  const InSpace&,
					  std::size_t&, std::size_t&, std::size_t& f3d_i,
				          const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_frame_3D(pt,*(model->Frames3D()[f3d_i++]));
  };
  
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( PointType& pt,
					    const InSpaceTuple& space_in,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				            const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( PointType& pt,
					    const InSpaceTuple& space_in,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				            const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>
    >,  
  void >::type read_one_joint_coord_impl( PointType& pt,
					  const InSpace& space_in,
					  std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				          const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_joint_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size<InSpace> >::type >(pt, space_in, gen_i, f2d_i, f3d_i, model);
  };
  
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( PointType& pt,
					    const InSpaceTuple& space_in,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_joint_coordinates_impl< typename boost::mpl::prior<Idx>::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
    read_one_joint_coord_impl(get<Idx::type::value>(pt),get<Idx::type::value>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( PointType& pt,
					    const InSpaceTuple& space_in,
					    std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_one_joint_coord_impl(get<0>(pt),get<0>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type read_joint_coordinates_impl( PointType& pt,
				    const InSpaceTuple& space_in,
				    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_joint_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size< InSpaceTuple > >::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type read_joint_coordinates_impl( PointType& pt,
				    const InSpaceTuple& space_in,
				    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_one_joint_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
					       const InSpace&,
					       std::size_t& gen_i, std::size_t&, std::size_t&,
				               const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->DependentCoords()[gen_i++]->mFrame) = get_gen_coord(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
					       const InSpace&,
					       std::size_t&, std::size_t& f2d_i, std::size_t&,
				               const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->DependentFrames2D()[f2d_i++]->mFrame) = get_frame_2D(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
					       const InSpace&,
					       std::size_t&, std::size_t&, std::size_t& f3d_i,
				               const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->DependentFrames3D()[f3d_i++]->mFrame) = get_frame_3D(pt);
  };
  
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_dependent_coordinates_impl( const PointType& pt,
					         const InSpaceTuple& space_in,
					         std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                 const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_dependent_coordinates_impl( const PointType& pt,
					         const InSpaceTuple& space_in,
					         std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                 const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>
    >,  
  void >::type write_one_dependent_coord_impl( const PointType& pt,
					       const InSpace& space_in,
					       std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				               const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_dependent_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size<InSpace> >::type >(pt, space_in, gen_i, f2d_i, f3d_i, model);
  };
  
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_dependent_coordinates_impl( const PointType& pt,
					         const InSpaceTuple& space_in,
					         std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                 const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_dependent_coordinates_impl< typename boost::mpl::prior<Idx>::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
    write_one_dependent_coord_impl(get<Idx::type::value>(pt),get<Idx::type::value>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_dependent_coordinates_impl( const PointType& pt,
					         const InSpaceTuple& space_in,
					         std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                 const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_one_dependent_coord_impl(get<0>(pt),get<0>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type write_dependent_coordinates_impl( const PointType& pt,
				         const InSpaceTuple& space_in,
				         const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_dependent_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size< InSpaceTuple > >::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type write_dependent_coordinates_impl( const PointType& pt,
				                 const InSpaceTuple& space_in,
				                 const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_dependent_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type read_one_dependent_coord_impl( PointType& pt,
					      const InSpace&,
					      std::size_t& gen_i, std::size_t&, std::size_t&,
				              const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_gen_coord(pt,*(model->DependentCoords()[gen_i++]->mFrame));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type read_one_dependent_coord_impl( PointType& pt,
					      const InSpace&,
					      std::size_t&, std::size_t& f2d_i, std::size_t&,
				              const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_frame_2D(pt,*(model->DependentFrames2D()[f2d_i++]->mFrame));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type read_one_dependent_coord_impl( PointType& pt,
					      const InSpace&,
					      std::size_t&, std::size_t&, std::size_t& f3d_i,
				              const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_frame_3D(pt,*(model->DependentFrames3D()[f3d_i++]->mFrame));
  };
  
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_dependent_coordinates_impl( PointType& pt,
					        const InSpaceTuple& space_in,
					        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_dependent_coordinates_impl( PointType& pt,
					        const InSpaceTuple& space_in,
					        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>
    >,  
  void >::type read_one_dependent_coord_impl( PointType& pt,
					      const InSpace& space_in,
					      std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				              const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_dependent_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size<InSpace> >::type >(pt, space_in, gen_i, f2d_i, f3d_i, model);
  };
  
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_dependent_coordinates_impl( PointType& pt,
					        const InSpaceTuple& space_in,
					        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_dependent_coordinates_impl< typename boost::mpl::prior<Idx>::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
    read_one_dependent_coord_impl(get<Idx::type::value>(pt),get<Idx::type::value>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_dependent_coordinates_impl( PointType& pt,
					        const InSpaceTuple& space_in,
					        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				                const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_one_dependent_coord_impl(get<0>(pt),get<0>(space_in),gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type read_dependent_coordinates_impl( PointType& pt,
				        const InSpaceTuple& space_in,
				        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_dependent_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size< InSpaceTuple > >::type >(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>
    >,  
  void >::type read_dependent_coordinates_impl( PointType& pt,
				        const InSpaceTuple& space_in,
				        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_one_dependent_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
  
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_normal_joint_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,  
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    ik_calc.lower_bounds[gen_i] = get<0>(space_out).origin() - get<0>(space_out).get_radius();
    ik_calc.upper_bounds[gen_i] = get<0>(space_out).origin() + get<0>(space_out).get_radius();
    ik_calc.lower_bounds[offset + gen_i] = -std::numeric_limits<double>::infinity();
    ik_calc.upper_bounds[offset + gen_i] = std::numeric_limits<double>::infinity();
    
    ++gen_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_normal_joint_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,  
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    ik_calc.lower_bounds[gen_i] = get<0>(space_out).origin() - get<0>(space_out).get_radius();
    ik_calc.upper_bounds[gen_i] = get<0>(space_out).origin() + get<0>(space_out).get_radius();
    ik_calc.lower_bounds[offset + gen_i] = get<1>(space_out).origin() - get<1>(space_out).get_radius();
    ik_calc.upper_bounds[offset + gen_i] = get<1>(space_out).origin() + get<1>(space_out).get_radius();
    
    ++gen_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se2_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i] = p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 1] = p_lo[1];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i] = p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 1] = p_up[1];
    
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 2] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 3] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = std::numeric_limits< double >::infinity();
    
    ++f2d_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se2_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i] = p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 1] = p_lo[1];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i] = p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 1] = p_up[1];
    
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 2] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 3] = std::numeric_limits< double >::infinity();
    
    vect<double,2> v_o = get<1>(get<0>(space_out)).origin(); 
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i] = v_o[0] - get<1>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = v_o[1] - get<1>(get<0>(space_out)).get_radius();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i] = v_o[0] + get<1>(get<0>(space_out)).get_radius();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = v_o[1] + get<1>(get<0>(space_out)).get_radius();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = get<1>(get<1>(space_out)).origin() - get<1>(get<1>(space_out)).get_radius();
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = get<1>(get<1>(space_out)).origin() + get<1>(get<1>(space_out)).get_radius();
    
    ++f2d_i;
  };
  
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se3_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i] = p_lo[0];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 1] = p_lo[1];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 2] = p_lo[2];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i] = p_up[0];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 1] = p_up[1];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 2] = p_up[2];
    
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 1] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 1] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 2] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    
    ++f3d_i;
    
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se3_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i] = p_lo[0];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 1] = p_lo[1];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 2] = p_lo[2];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i] = p_up[0];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 1] = p_up[1];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 2] = p_up[2];
    
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i] = v_o[0] - v_r;
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 1] = v_o[1] - v_r;
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 2] = v_o[2] - v_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i] = v_o[0] + v_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 1] = v_o[1] + v_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 2] = v_o[2] + v_r;
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 3] = w_o[0] - w_r;
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 4] = w_o[1] - w_r;
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 5] = w_o[2] - w_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 3] = w_o[0] + w_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 4] = w_o[1] + w_r;
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 5] = w_o[2] + w_r;
    
    ++f3d_i;
    
  };
  
  
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  
  template <typename OutSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpace>,
      is_se2_space<OutSpace>,
      is_se3_space<OutSpace>
    >,  
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_bounds_impl< typename boost::mpl::prior< arithmetic_tuple_size<OutSpace> >::type >(space_out, offset, gen_i, f2d_i, f3d_i, ik_calc, model);
  };
  
  
  
  template <typename Idx, typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_bounds_impl< typename boost::mpl::prior<Idx>::type >(space_out,offset,gen_i,f2d_i,f3d_i,ik_calc,model);
    
    write_one_joint_bound_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,ik_calc,model);
  };
  
  template <typename Idx, typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_one_joint_bound_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,ik_calc,model);
  };
  
  template <typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpaceTuple>,
      is_se2_space<OutSpaceTuple>,
      is_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_bounds_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > >::type >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,ik_calc,model);
  };
  
  template <typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpaceTuple>,
      is_se2_space<OutSpaceTuple>,
      is_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_joint_bound_impl(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,ik_calc,model);
  };
  
  
  
  
  
  
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_normal_joint_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,  
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
			                   const shared_ptr< kte::manipulator_kinematics_model >&,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    centers[offset + gen_i] = 0.0;
    diag_gains[gen_i] = get<0>(space_out).get_radius(); diag_gains[gen_i] = 1.0 / (diag_gains[gen_i] * diag_gains[gen_i]);
    diag_gains[offset + gen_i] = 1.0;
    ++gen_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_normal_joint_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,  
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    centers[offset + gen_i] = get<1>(space_out).origin();
    diag_gains[gen_i] = get<0>(space_out).get_radius(); diag_gains[gen_i] = 1.0 / (diag_gains[gen_i] * diag_gains[gen_i]);
    diag_gains[offset + gen_i] = get<1>(space_out).get_radius(); diag_gains[offset + gen_i] = 1.0 / (diag_gains[offset + gen_i] * diag_gains[offset + gen_i]);
    ++gen_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se2_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    
    centers[offset + model->Coords().size() + 3 * f2d_i] = 0.0;
    centers[offset + model->Coords().size() + 3 * f2d_i + 1] = 0.0;
    
    centers[offset + model->Coords().size() + 3 * f2d_i + 2] = 0.0;
    
    diag_gains[model->Coords().size() + 4 * f2d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    
    diag_gains[model->Coords().size() + 4 * f2d_i + 2] = 1.0; 
    diag_gains[model->Coords().size() + 4 * f2d_i + 3] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i] = 1.0; 
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 1] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 2] = 1.0; 
    
    ++f2d_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se2_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,2> v_o = get<1>(get<0>(space_out)).origin(); 
    
    centers[offset + model->Coords().size() + 3 * f2d_i] = v_o[0];
    centers[offset + model->Coords().size() + 3 * f2d_i + 1] = v_o[1];
    
    centers[offset + model->Coords().size() + 3 * f2d_i + 2] = get<1>(get<1>(space_out)).origin();
    
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    
    diag_gains[model->Coords().size() + 4 * f2d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    
    diag_gains[model->Coords().size() + 4 * f2d_i + 2] = 1.0; 
    diag_gains[model->Coords().size() + 4 * f2d_i + 3] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i] = 1.0 / (get<1>(get<0>(space_out)).get_radius() * get<1>(get<0>(space_out)).get_radius()); 
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 1] = 1.0 / (get<1>(get<0>(space_out)).get_radius() * get<1>(get<0>(space_out)).get_radius()); 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 2] = 1.0 / (get<1>(get<1>(space_out)).get_radius() * get<1>(get<1>(space_out)).get_radius()); 
    
    ++f2d_i;
  };
  
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se3_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    
    diag_gains[base_offset + 7 * f3d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[base_offset + 7 * f3d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    diag_gains[base_offset + 7 * f3d_i + 2] = 4.0 / ((p_up[2] - p_lo[2]) * (p_up[2] - p_lo[2])); 
    
    diag_gains[base_offset + 7 * f3d_i + 3] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 4] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 5] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 6] = 1.0; 
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    centers[offset + base_offset + 6 * f3d_i] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 1] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 2] = 0.0;
    
    centers[offset + base_offset + 6 * f3d_i + 3] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 4] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 5] = 0.0;
    
    diag_gains[offset + base_offset + 6 * f3d_i] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 1] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 2] = 1.0; 
    
    diag_gains[offset + base_offset + 6 * f3d_i + 3] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 4] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 5] = 1.0; 
    
    ++f3d_i;
  };
  
  template <typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_se3_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    
    diag_gains[base_offset + 7 * f3d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[base_offset + 7 * f3d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    diag_gains[base_offset + 7 * f3d_i + 2] = 4.0 / ((p_up[2] - p_lo[2]) * (p_up[2] - p_lo[2])); 
    
    diag_gains[base_offset + 7 * f3d_i + 3] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 4] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 5] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 6] = 1.0; 
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    
    centers[offset + base_offset + 6 * f3d_i] = v_o[0];
    centers[offset + base_offset + 6 * f3d_i + 1] = v_o[1];
    centers[offset + base_offset + 6 * f3d_i + 2] = v_o[2];
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    
    centers[offset + base_offset + 6 * f3d_i + 3] = w_o[0];
    centers[offset + base_offset + 6 * f3d_i + 4] = w_o[1];
    centers[offset + base_offset + 6 * f3d_i + 5] = w_o[2];
    
    v_r = 1.0 / (v_r * v_r);
    diag_gains[offset + base_offset + 6 * f3d_i] = v_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 1] = v_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 2] = v_r; 
    
    w_r = 1.0 / (w_r * w_r);
    diag_gains[offset + base_offset + 6 * f3d_i + 3] = w_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 4] = w_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 5] = w_r; 
    
    ++f3d_i;
  };
  
  
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains);
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains);
  
  
  template <typename OutSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpace>,
      is_se2_space<OutSpace>,
      is_se3_space<OutSpace>
    >,  
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    write_joints_quad_info_impl< typename boost::mpl::prior< arithmetic_tuple_size<OutSpace> >::type >(space_out, offset, gen_i, f2d_i, f3d_i, model, centers, diag_gains);
  };
  
  
  
  template <typename Idx, typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_joints_quad_info_impl< typename boost::mpl::prior<Idx>::type >(space_out,offset,gen_i,f2d_i,f3d_i,model,centers,diag_gains);
    
    write_one_joint_quad_info_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,model,centers,diag_gains);
  };
  
  template <typename Idx, typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_one_joint_quad_info_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpaceTuple>,
      is_se2_space<OutSpaceTuple>,
      is_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joints_quad_info_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > >::type >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<OutSpaceTuple>,
      is_se2_space<OutSpaceTuple>,
      is_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_joint_quad_info_impl(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,model,centers,diag_gains);
  };
  
  
  
  
  
  
  
  
  
  
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_joint_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,  
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    ik_calc.lower_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() - get<0>(space_out).get_radius());
    ik_calc.upper_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() + get<0>(space_out).get_radius());
    ik_calc.lower_bounds[offset + gen_i] = -j_limits.gen_speed_limits[gen_i];
    ik_calc.upper_bounds[offset + gen_i] = j_limits.gen_speed_limits[gen_i];
    
    ++gen_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_joint_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,  
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    ik_calc.lower_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() - get<0>(space_out).get_radius());
    ik_calc.upper_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() + get<0>(space_out).get_radius());
    ik_calc.lower_bounds[offset + gen_i] = j_limits.gen_accel_limits[gen_i] * (get<1>(space_out).origin() - get<1>(space_out).get_radius());
    ik_calc.upper_bounds[offset + gen_i] = j_limits.gen_accel_limits[gen_i] * (get<1>(space_out).origin() + get<1>(space_out).get_radius());
    ++gen_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se2_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 1] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_lo[1];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 1] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_up[1];
    
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 2] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 3] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i] = -j_limits.frame2D_speed_limits[2 * f2d_i];
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = -j_limits.frame2D_speed_limits[2 * f2d_i];
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i] = j_limits.frame2D_speed_limits[2 * f2d_i];
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = j_limits.frame2D_speed_limits[2 * f2d_i];
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = -j_limits.frame2D_speed_limits[2 * f2d_i + 1];
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = j_limits.frame2D_speed_limits[2 * f2d_i + 1];
    
    ++f2d_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se2_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    vect<double,2> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,2> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 1] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_lo[1];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 1] = j_limits.frame2D_speed_limits[2 * f2d_i] * p_up[1];
    
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 4 * f2d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 2] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 4 * f2d_i + 3] = std::numeric_limits< double >::infinity();
    
    vect<double,2> v_o = get<1>(get<0>(space_out)).origin(); 
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i] = j_limits.frame2D_accel_limits[2 * f2d_i] * (v_o[0] - get<1>(get<0>(space_out)).get_radius());
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = j_limits.frame2D_accel_limits[2 * f2d_i] * (v_o[1] - get<1>(get<0>(space_out)).get_radius());
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i] = j_limits.frame2D_accel_limits[2 * f2d_i] * (v_o[0] + get<1>(get<0>(space_out)).get_radius());
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 1] = j_limits.frame2D_accel_limits[2 * f2d_i] * (v_o[1] + get<1>(get<0>(space_out)).get_radius());
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = j_limits.frame2D_accel_limits[2 * f2d_i + 1] * (get<1>(get<1>(space_out)).origin() - get<1>(get<1>(space_out)).get_radius());
    ik_calc.upper_bounds[offset + model->Coords().size() + 3 * f2d_i + 2] = j_limits.frame2D_accel_limits[2 * f2d_i + 1] * (get<1>(get<1>(space_out)).origin() + get<1>(get<1>(space_out)).get_radius());
    
    ++f2d_i;
  };
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se3_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[0];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[1];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[2];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[0];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[1];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[2];
    
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 1] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 2] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i];
    
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 3] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 4] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 5] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 3] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 4] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 5] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    
    ++f3d_i;
    
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se3_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[0];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[1];
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[2];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[0];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[1];
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[2];
    
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[base_offset + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[base_offset + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[0] - v_r);
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[1] - v_r);
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[2] - v_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[0] + v_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[1] + v_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[2] + v_r);
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[0] - w_r);
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[1] - w_r);
    ik_calc.lower_bounds[offset + base_offset + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[2] - w_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[0] + w_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[1] + w_r);
    ik_calc.upper_bounds[offset + base_offset + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[2] + w_r);
    
    ++f3d_i;
    
  };
  
  
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model);
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpace>,
      is_rate_limited_se2_space<OutSpace>,
      is_rate_limited_se3_space<OutSpace>
    >, 
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_bounds_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, offset, gen_i, f2d_i, f3d_i, j_limits, ik_calc, model);
  };
  
  
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_bounds_impl< typename boost::mpl::prior<Idx>::type >(space_out,offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model);
    
    write_one_joint_bound_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model);
  };
  
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					kte::manip_clik_calculator& ik_calc,
					const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_one_joint_bound_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model);
  };
  
  template <typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpaceTuple>,
      is_rate_limited_se2_space<OutSpaceTuple>,
      is_rate_limited_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        const RateLimitMap& j_limits,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_bounds_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > >::type >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,j_limits,ik_calc,model);
  };
  
  template <typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpaceTuple>,
      is_rate_limited_se2_space<OutSpaceTuple>,
      is_rate_limited_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        const RateLimitMap& j_limits,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_joint_bound_impl(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,j_limits,ik_calc,model);
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_joint_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<1>
      >
    >,  
  void >::type write_one_joint_quad_info_impl( const OutSpace&,
					       std::size_t offset,
					       std::size_t& gen_i, std::size_t&, std::size_t&,
					       const RateLimitMap& j_limits,
			                       const shared_ptr< kte::manipulator_kinematics_model >&,
				               vect_n<double>& centers, vect_n<double>& diag_gains) {
    centers[offset + gen_i] = 0.0;
    diag_gains[gen_i] = j_limits.gen_speed_limits[gen_i]; diag_gains[gen_i] = 1.0 / (diag_gains[gen_i] * diag_gains[gen_i]);
    diag_gains[offset + gen_i] = 1.0;
    ++gen_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_joint_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< OutSpace, time_topology >,
	boost::mpl::size_t<0>
      >
    >,  
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >&,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    centers[offset + gen_i] = j_limits.gen_accel_limits[gen_i] * get<1>(space_out).origin();
    diag_gains[gen_i] = j_limits.gen_speed_limits[gen_i]; diag_gains[gen_i] = 1.0 / (diag_gains[gen_i] * diag_gains[gen_i]);
    diag_gains[offset + gen_i] = j_limits.gen_accel_limits[gen_i]; diag_gains[offset + gen_i] = 1.0 / (diag_gains[offset + gen_i] * diag_gains[offset + gen_i]);
    ++gen_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se2_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    std::size_t base_offset = model->Coords().size();
    centers[offset + base_offset + 3 * f2d_i] = 0.0;
    centers[offset + base_offset + 3 * f2d_i + 1] = 0.0;
    
    centers[offset + base_offset + 3 * f2d_i + 2] = 0.0;
    
    diag_gains[base_offset + 4 * f2d_i] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i] * j_limits.frame2D_speed_limits[2 * f2d_i]); 
    diag_gains[base_offset + 4 * f2d_i + 1] = diag_gains[base_offset + 4 * f2d_i]; 
    
    diag_gains[base_offset + 4 * f2d_i + 2] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    diag_gains[base_offset + 4 * f2d_i + 3] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    
    diag_gains[offset + base_offset + 3 * f2d_i] = 1.0; 
    diag_gains[offset + base_offset + 3 * f2d_i + 1] = 1.0; 
    
    diag_gains[offset + base_offset + 3 * f2d_i + 2] = 1.0; 
    
    ++f2d_i;
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se2_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    std::size_t base_offset = model->Coords().size();
    
    vect<double,2> v_o = get<1>(get<0>(space_out)).origin(); 
    centers[offset + base_offset + 3 * f2d_i] = j_limits.frame2D_accel_limits[2 * f2d_i] * v_o[0];
    centers[offset + base_offset + 3 * f2d_i + 1] = j_limits.frame2D_accel_limits[2 * f2d_i] * v_o[1];
    
    centers[offset + base_offset + 3 * f2d_i + 2] = j_limits.frame2D_accel_limits[2 * f2d_i + 1] * get<1>(get<1>(space_out)).origin();
    
    diag_gains[base_offset + 4 * f2d_i] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i] * j_limits.frame2D_speed_limits[2 * f2d_i]); 
    diag_gains[base_offset + 4 * f2d_i + 1] = diag_gains[model->Coords().size() + 4 * f2d_i]; 
    
    diag_gains[base_offset + 4 * f2d_i + 2] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    diag_gains[base_offset + 4 * f2d_i + 3] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    
    diag_gains[offset + base_offset + 3 * f2d_i] = 1.0 / (j_limits.frame2D_accel_limits[2 * f2d_i] * j_limits.frame2D_accel_limits[2 * f2d_i]); 
    diag_gains[offset + base_offset + 3 * f2d_i + 1] = diag_gains[offset + base_offset + 3 * f2d_i]; 
    
    diag_gains[offset + base_offset + 3 * f2d_i + 2] = 1.0 / (j_limits.frame2D_accel_limits[2 * f2d_i + 1] * j_limits.frame2D_accel_limits[2 * f2d_i + 1]); 
    
    ++f2d_i;
  };
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se3_space<OutSpace>,
      boost::mpl::less<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<1>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace&,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    
    diag_gains[base_offset + 7 * f3d_i] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i] * j_limits.frame3D_speed_limits[2 * f3d_i]); 
    diag_gains[base_offset + 7 * f3d_i + 1] = diag_gains[base_offset + 7 * f3d_i]; 
    diag_gains[base_offset + 7 * f3d_i + 2] = diag_gains[base_offset + 7 * f3d_i]; 
    
    diag_gains[base_offset + 7 * f3d_i + 3] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i + 1] * j_limits.frame3D_speed_limits[2 * f3d_i + 1]); 
    diag_gains[base_offset + 7 * f3d_i + 4] = diag_gains[base_offset + 7 * f3d_i + 3]; 
    diag_gains[base_offset + 7 * f3d_i + 5] = diag_gains[base_offset + 7 * f3d_i + 3]; 
    diag_gains[base_offset + 7 * f3d_i + 6] = diag_gains[base_offset + 7 * f3d_i + 3]; 
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    centers[offset + base_offset + 6 * f3d_i] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 1] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 2] = 0.0;
    
    centers[offset + base_offset + 6 * f3d_i + 3] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 4] = 0.0;
    centers[offset + base_offset + 6 * f3d_i + 5] = 0.0;
    
    diag_gains[offset + base_offset + 6 * f3d_i] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 1] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 2] = 1.0; 
    
    diag_gains[offset + base_offset + 6 * f3d_i + 3] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 4] = 1.0; 
    diag_gains[offset + base_offset + 6 * f3d_i + 5] = 1.0; 
    
    ++f3d_i;
    
  };
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::and_<
      is_rate_limited_se3_space<OutSpace>,
      boost::mpl::greater<
        max_derivation_order< arithmetic_tuple_element<0,OutSpace>, time_topology >,
	boost::mpl::size_t<0>
      >
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    std::size_t base_offset = model->Coords().size() + 4 * model->Frames2D().size();
    
    diag_gains[base_offset + 7 * f3d_i] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i] * j_limits.frame3D_speed_limits[2 * f3d_i]); 
    diag_gains[base_offset + 7 * f3d_i + 1] = diag_gains[base_offset + 7 * f3d_i]; 
    diag_gains[base_offset + 7 * f3d_i + 2] = diag_gains[base_offset + 7 * f3d_i]; 
    
    diag_gains[base_offset + 7 * f3d_i + 3] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 4] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 5] = 1.0; 
    diag_gains[base_offset + 7 * f3d_i + 6] = 1.0; 
    
    base_offset = model->Coords().size() + 3 * model->Frames2D().size();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    centers[offset + base_offset + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[0];
    centers[offset + base_offset + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[1];
    centers[offset + base_offset + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[2];
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    centers[offset + base_offset + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[0];
    centers[offset + base_offset + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[1];
    centers[offset + base_offset + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[2];
    
    v_r = 1.0 / (j_limits.frame3D_accel_limits[2 * f3d_i] * j_limits.frame3D_accel_limits[2 * f3d_i]);
    diag_gains[offset + base_offset + 6 * f3d_i] = v_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 1] = v_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 2] = v_r; 
    
    w_r = 1.0 / (j_limits.frame3D_accel_limits[2 * f3d_i + 1] * j_limits.frame3D_accel_limits[2 * f3d_i + 1]);
    diag_gains[offset + base_offset + 6 * f3d_i + 3] = w_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 4] = w_r; 
    diag_gains[offset + base_offset + 6 * f3d_i + 5] = w_r; 
    
    ++f3d_i;
    
  };
  
  
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains);
  
  //declaration only.
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains);
  
  
  template <typename OutSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpace>,
      is_rate_limited_se2_space<OutSpace>,
      is_rate_limited_se3_space<OutSpace>
    >, 
  void >::type write_one_joint_quad_info_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					   const RateLimitMap& j_limits,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    write_joints_quad_info_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpace > >::type >(space_out, offset, gen_i, f2d_i, f3d_i, j_limits, model, centers, diag_gains);
  };
  
  
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_joints_quad_info_impl< typename boost::mpl::prior<Idx>::type >(space_out,offset,gen_i,f2d_i,f3d_i,j_limits,model,centers,diag_gains);
    
    write_one_joint_quad_info_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,model,centers,diag_gains);
  };
  
  template <typename Idx, typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
					std::size_t offset,
				        std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
					const RateLimitMap& j_limits,
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_one_joint_quad_info_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpaceTuple>,
      is_rate_limited_se2_space<OutSpaceTuple>,
      is_rate_limited_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
			        const RateLimitMap& j_limits,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joints_quad_info_impl< typename boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > >::type >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,j_limits,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_rate_limited_joint_space<OutSpaceTuple>,
      is_rate_limited_se2_space<OutSpaceTuple>,
      is_rate_limited_se3_space<OutSpaceTuple>
    >,  
  void >::type write_joints_quad_info_impl( const OutSpaceTuple& space_out,
			        const RateLimitMap& j_limits,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_joint_quad_info_impl(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,j_limits,model,centers,diag_gains);
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
  
};





/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates (both 
 * generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_direct_kin_map : public shared_object {
  public:
    
    typedef manip_direct_kin_map self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    
    manip_direct_kin_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >()) :
                         model(aModel) { };
    
    /**
     * This function template performs a forward kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (joint-space).
     * \tparam OutSpace The type of the output space (end-effector space).
     * \param pt The point in the input space, i.e. the joint coordinates.
     * \param space_in The input space, i.e. the joint-space.
     * \param space_out The output space, i.e. the end-effector space.
     * \return A point in the output space, i.e. the end-effector coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_joint_coordinates_impl(pt, space_in, model);
    
      model->doMotion();
    
      detail::read_dependent_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400012,1,"manip_direct_kin_map",shared_object)
    
    
};



/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates 
 * (both generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename RateLimitMap>
class manip_rl_direct_kin_map : public shared_object {
  public:
    
    typedef manip_rl_direct_kin_map<RateLimitMap> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    manip_rl_direct_kin_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
			    const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >()) : 
                            model(aModel),
                            joint_limits_map(aJointLimitMap) { };
    
    /**
     * This function template performs a forward kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (rate-limited joint-space).
     * \tparam OutSpace The type of the output space (end-effector space).
     * \param pt The point in the input space, i.e. the rate-limited joint coordinates.
     * \param space_in The input space, i.e. the rate-limited joint-space.
     * \param space_out The output space, i.e. the end-effector space.
     * \return A point in the output space, i.e. the end-effector coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      typedef typename get_rate_illimited_space< InSpace >::type NormalJointSpace;
      NormalJointSpace normal_j_space;
      typename topology_traits<NormalJointSpace>::point_type pt_inter = joint_limits_map->map_to_space(pt, space_in, normal_j_space);
      detail::write_joint_coordinates_impl(pt_inter, normal_j_space, model);
    
      model->doMotion();
      
      detail::read_dependent_coordinates_impl(result, space_out, model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400013,1,"manip_rl_direct_kin_map",shared_object)
    
    
};




/**
 * This class is a factory class for creating quadratic cost evaluators for the cost 
 * associated to a joint-state of the manipulator.
 */
class clik_quad_cost_factory : public shared_object {
  public:
    typedef clik_quad_cost_factory self;
    
    /** This holds the posture that is the most desirable for the manipulator, that is, a set of 
     * joint coordinates (in normal joint-space) that is well conditioned w.r.t. some performance
     * index like manipulability or condition number of the Jacobian. This is used in the inverse
     * kinematics mapping to optimize the condition of the posture of the robot, that is, if there 
     * is any freedom in the matter (if the manipulator has some redundancy to exploit) otherwise 
     * it has no effect. */
    vect_n<double> preferred_posture;
    
    /**
     * Parametrized constructor.
     * \param aPreferredPosture The preferred posture for the manipulator.
     */
    clik_quad_cost_factory(const vect_n<double>& aPreferredPosture = vect_n<double>()) : 
                           preferred_posture(aPreferredPosture) { };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace>
    shared_ptr< optim::quadratic_cost_evaluator > create_evaluator(
      const JointSpace& aJSpace,
      const shared_ptr< kte::manipulator_kinematics_model >& aModel
    ) const {
      
      vect_n<double> centers(aModel->getJointPositionsCount() + aModel->getJointVelocitiesCount());
      vect_n<double> diag_gains(aModel->getJointPositionsCount() + aModel->getJointVelocitiesCount());
      if(preferred_posture.size() == aModel->getJointPositionsCount())
        std::copy(preferred_posture.begin(), preferred_posture.end(), centers.begin());
      else
	std::fill(centers.begin(), centers.begin() + aModel->getJointPositionsCount(), 0.0);
      
      detail::write_joints_quad_info_impl(
        aJSpace, 
        aModel,
        centers,
        diag_gains
      );
      
      if(preferred_posture.size() != aModel->getJointPositionsCount())
	std::fill(diag_gains.begin(), diag_gains.begin() + aModel->getJointPositionsCount(), 0.0);
      
      return shared_ptr< optim::quadratic_cost_evaluator >( 
        new optim::quadratic_cost_evaluator(centers, mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(diag_gains))));
    };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aJLimits The joint space rate-limits associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace, typename RateLimitMap>
    shared_ptr< optim::quadratic_cost_evaluator > create_evaluator(
      const JointSpace& aJSpace,
      const RateLimitMap& aJLimits,
      const shared_ptr< kte::manipulator_kinematics_model >& aModel
    ) const {
      
      vect_n<double> centers(aModel->getJointPositionsCount() + aModel->getJointVelocitiesCount());
      vect_n<double> diag_gains(aModel->getJointPositionsCount() + aModel->getJointVelocitiesCount());
      if(preferred_posture.size() == aModel->getJointPositionsCount())
        std::copy(preferred_posture.begin(), preferred_posture.end(), centers.begin());
      else
	std::fill(centers.begin(), centers.begin() + aModel->getJointPositionsCount(), 0.0);
      
      detail::write_joints_quad_info_impl(
        aJSpace, 
	aJLimits,
        aModel,
        centers,
        diag_gains
      );
      
      if(preferred_posture.size() != aModel->getJointPositionsCount())
	std::fill(diag_gains.begin(), diag_gains.begin() + aModel->getJointPositionsCount(), 0.0);
      
      return shared_ptr< optim::quadratic_cost_evaluator >( 
        new optim::quadratic_cost_evaluator(centers, mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(diag_gains))));
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(preferred_posture);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(preferred_posture);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400017,1,"clik_quad_cost_factory",shared_object)
    
};



/**
 * This cost evaluator attempts to keep two joints from being too straight, i.e., tries 
 * to keep them at a +/- 90 degrees angle. The cost function itself is a sin-squared 
 * function applied on each joint in its internal list.
 */
class clik_bent_joints_cost_eval : public ReaK::optim::cost_evaluator {
  public:
    std::vector<int> joint_ids; ///< Holds the indices of the joints to keep bent.
    
    /**
     * Parametrized Constructor.
     * \param aJointIDs The indices of the joints to keep bent.
     */
    clik_bent_joints_cost_eval(const std::vector<int>& aJointIDs) : joint_ids(aJointIDs) { };
    
    virtual double compute_cost(const vect_n<double>& x) const {
      double sum = 0.0;
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
	double s = std::sin( x[joint_ids[i]] );
	sum += 1.0 - s * s;
      };
      return sum;
    };
    
    virtual vect_n<double> compute_cost_grad(const vect_n<double>& x) const {
      vect_n<double> result(x.size(),0.0);
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
	result[ joint_ids[i] ] = -2.0 * std::sin( x[joint_ids[i]] ) * std::cos( x[joint_ids[i]] );
      };
      return result;
    };
    
    virtual void compute_cost_hessian(mat<double,mat_structure::symmetric>& H, const vect_n<double>& x, double, const vect_n<double>&) const {
      H.set_row_count(x.size());
      for(std::size_t i = 0; i < joint_ids.size(); ++i) {
	double s = std::sin( x[joint_ids[i]] ); s *= s;
	double c = std::cos( x[joint_ids[i]] ); c *= c;
	H(joint_ids[i], joint_ids[i]) = 2.0 * (s - c);
      };
    };
  
};




/**
 * This class is a factory class for creating quadratic cost evaluators for the cost 
 * associated to a joint-state of the manipulator.
 */
class clik_bent_joints_cost_factory : public shared_object {
  public:
    typedef clik_bent_joints_cost_factory self;
    
    /** Holds the indices of the joints to keep bent. */
    std::vector<int> joint_ids;
    
    /**
     * Parametrized constructor.
     * \param aJointIDs The indices of the joints to keep bent.
     */
    clik_bent_joints_cost_factory(const std::vector<int>& aJointIDs) : joint_ids(aJointIDs) { };
    
    clik_bent_joints_cost_factory(int J1 = -1, int J2 = -1, int J3 = -1, int J4 = -1) {
      if(J1 >= 0)
	joint_ids.push_back(J1);
      if(J2 >= 0)
	joint_ids.push_back(J2);
      if(J3 >= 0)
	joint_ids.push_back(J3);
      if(J4 >= 0)
	joint_ids.push_back(J4);
    };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace>
    shared_ptr< clik_bent_joints_cost_eval > create_evaluator(
      const JointSpace&,
      const shared_ptr< kte::manipulator_kinematics_model >&
    ) const {
      return shared_ptr< clik_bent_joints_cost_eval >( 
        new clik_bent_joints_cost_eval(joint_ids));
    };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aJLimits The joint space rate-limits associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace, typename RateLimitMap>
    shared_ptr< clik_bent_joints_cost_eval > create_evaluator(
      const JointSpace& ,
      const RateLimitMap& ,
      const shared_ptr< kte::manipulator_kinematics_model >& 
    ) const {
      return shared_ptr< clik_bent_joints_cost_eval >( 
        new clik_bent_joints_cost_eval(joint_ids));
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(joint_ids);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(joint_ids);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400019,1,"clik_bent_joints_cost_factory",shared_object)
    
};




/**
 * This class is a factory class for creating mixed cost evaluators for the cost 
 * associated to a joint-state of the manipulator. The costs will be calculated as 
 * an addition of a quadratic function for the velocities and a user-given cost 
 * evaluator.
 */
template <typename CostEvalFactory>
class clik_mixed_cost_factory : public shared_object {
  public:
    typedef clik_mixed_cost_factory<CostEvalFactory> self;
    
    /** Points to a user-given cost evaluator factory */
    CostEvalFactory user_factory;
    
    /**
     * Parametrized constructor.
     * \param aUserFactory The user-given cost evaluator factory for the manipulator.
     */
    clik_mixed_cost_factory(const CostEvalFactory& aUserFactory = CostEvalFactory()) : 
                            user_factory(aUserFactory) { };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace>
    shared_ptr< optim::added_cost_evaluator > create_evaluator(
      const JointSpace& aJSpace,
      const shared_ptr< kte::manipulator_kinematics_model >& aModel
    ) const {
      
      clik_quad_cost_factory quad_fact;
      
      return shared_ptr< optim::added_cost_evaluator >(
	new optim::added_cost_evaluator(
	  quad_fact.create_evaluator(aJSpace, aModel),
	  user_factory.create_evaluator(aJSpace, aModel)
	)
      );
    };
    
    /**
     * This function creates a quadratic cost evaluator for a manipulator model and 
     * joint space.
     * \param aJSpace The joint space associated to the joints of the manipulator.
     * \param aJLimits The joint space rate-limits associated to the joints of the manipulator.
     * \param aModel The kinematics model of the manipulator.
     * \return The quadratic cost evaluator for the manipulator model.
     */
    template <typename JointSpace, typename RateLimitMap>
    shared_ptr< optim::quadratic_cost_evaluator > create_evaluator(
      const JointSpace& aJSpace,
      const RateLimitMap& aJLimits,
      const shared_ptr< kte::manipulator_kinematics_model >& aModel
    ) const {
      
      clik_quad_cost_factory quad_fact;
      
      return shared_ptr< optim::added_cost_evaluator >(
	new optim::added_cost_evaluator(
	  quad_fact.create_evaluator(aJSpace, aJLimits, aModel),
	  user_factory.create_evaluator(aJSpace, aJLimits, aModel)
	)
      );
    };
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(user_factory);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(user_factory);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400018,1,"clik_mixed_cost_factory",shared_object)
    
};


/**
 * This class template is a factory class to create inverse kinematics calculators.
 * \tparam JointSpace The joint-space type.
 * \tparam CostEvalFactory The factory type for creating cost evaluators.
 */
template <typename JointSpace, typename CostEvalFactory>
class manip_clik_calc_factory : public shared_object {
  public:
    typedef manip_clik_calc_factory<JointSpace,CostEvalFactory> self;
    
    typedef JointSpace joint_space_type;
    typedef CostEvalFactory cost_eval_factory_type;
    
    /** Points to the joint-space. */
    shared_ptr< JointSpace > j_space;
    /** Holds the cost evaluator factory. */
    CostEvalFactory cost_eval_factory;
    
    /** Holds the maximum radius use in the optimization steps. */
    double radius;
    /** Holds the initial mu-value used in the optimization process. */
    double mu;
    /** Holds the maximum number of optimization steps. */
    double max_iter; 
    /** Holds the absolute tolerance on the error on the end-effector state. */
    double tol; 
    /** Holds the eta value for the optimization process (close to 0). */
    double eta; 
    /** Holds the tau value for the optimization process (close to 1). */
    double tau;
    
    shared_ptr< optim::cost_evaluator > cost_eval;
    
    /**
     * Parametrized constructor.
     * \param aJSpace A pointer to the joint-space.
     * \param aCostEvalFactory The cost-evaluator factory object.
     * \param aRadius The maximum radius use in the optimization steps.
     * \param aMu The initial mu-value used in the optimization process.
     * \param aMaxIter The maximum number of optimization steps.
     * \param aTol The absolute tolerance on the error on the end-effector state.
     * \param aEta The eta value for the optimization process (close to 0).
     * \param aTau The tau value for the optimization process (close to 1).
     */
    manip_clik_calc_factory(
      const shared_ptr< JointSpace >& aJSpace = shared_ptr< JointSpace >(),
      const CostEvalFactory& aCostEvalFactory = CostEvalFactory(),
      double aRadius = 5.0,
      double aMu = 0.1,
      double aMaxIter = 100,
      double aTol = 8e-4,
      double aEta = 1e-2,
      double aTau = 0.95) : 
      j_space(aJSpace), cost_eval_factory(aCostEvalFactory),
      radius(aRadius), mu(aMu), max_iter(aMaxIter), tol(aTol), eta(aEta), tau(aTau), cost_eval() { };
    
    /**
     * This function creates a CLIK calculator for the given manipulator model.
     * \param aModel The kinematics model for the manipulator.
     * \return The CLIK calculator for this given manipulator model.
     */
    kte::manip_clik_calculator create_calculator(const shared_ptr< kte::manipulator_kinematics_model >& aModel) {
      
      cost_eval = cost_eval_factory.create_evaluator(*j_space,aModel);
      
      kte::manip_clik_calculator ik_calc(
        aModel.get(), cost_eval,
        radius, mu, max_iter, tol, eta, tau
      );
      
      detail::write_joint_bounds_impl(
        *j_space, 
        ik_calc,
        aModel
      );
      
      return ik_calc;
    };
    
    /**
     * This function creates a CLIK calculator for the given manipulator model.
     * \param aJLimits The joint rate-limits associated to the manipulator model and joint-space.
     * \param aModel The kinematics model for the manipulator.
     * \return The CLIK calculator for this given manipulator model.
     */
    template <typename RateLimitMap>
    kte::manip_clik_calculator create_calculator(const RateLimitMap& aJLimits, const shared_ptr< kte::manipulator_kinematics_model >& aModel) {
      
      cost_eval = cost_eval_factory.create_evaluator(*j_space,aJLimits,aModel);
      
      kte::manip_clik_calculator ik_calc(
        aModel.get(), cost_eval,
        radius, mu, max_iter, tol, eta, tau
      );
    
      detail::write_joint_bounds_impl(
        *j_space, 
	aJLimits,
        ik_calc,
        aModel
      );
      
      return ik_calc;
    };
    
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(j_space)
        & RK_SERIAL_SAVE_WITH_NAME(cost_eval_factory)
        & RK_SERIAL_SAVE_WITH_NAME(radius)
        & RK_SERIAL_SAVE_WITH_NAME(mu)
        & RK_SERIAL_SAVE_WITH_NAME(max_iter)
        & RK_SERIAL_SAVE_WITH_NAME(tol)
        & RK_SERIAL_SAVE_WITH_NAME(eta)
        & RK_SERIAL_SAVE_WITH_NAME(tau);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(j_space)
        & RK_SERIAL_LOAD_WITH_NAME(cost_eval_factory)
        & RK_SERIAL_LOAD_WITH_NAME(radius)
        & RK_SERIAL_LOAD_WITH_NAME(mu)
        & RK_SERIAL_LOAD_WITH_NAME(max_iter)
        & RK_SERIAL_LOAD_WITH_NAME(tol)
        & RK_SERIAL_LOAD_WITH_NAME(eta)
        & RK_SERIAL_LOAD_WITH_NAME(tau);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400016,1,"manip_clik_calc_factory",shared_object)
    
};




/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam CLIKCalcFactory The factory type for creating a CLIK calculator.
 */
template <typename CLIKCalcFactory>
class manip_inverse_kin_map : public shared_object {
  public:
    
    typedef manip_inverse_kin_map<CLIKCalcFactory> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    /** This holds the inverse kinematics calculator factory. */
    CLIKCalcFactory clik_calc_factory;
    
    mutable kte::manip_clik_calculator ik_calc;
    
    manip_inverse_kin_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
                          const CLIKCalcFactory& aCLIKCalcFactory = CLIKCalcFactory()) :
                          model(aModel),
                          clik_calc_factory(aCLIKCalcFactory), ik_calc(aModel.get()) {
      ik_calc = clik_calc_factory.create_calculator(model);
    };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the joint-space.
     * \return A point in the output space, i.e. the joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
       
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      ik_calc.solveInverseKinematics();
      
      detail::read_joint_coordinates_impl(result,space_out,model);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc_factory);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc_factory);
      ik_calc = clik_calc_factory.create_calculator(model);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400014,1,"manip_inverse_kin_map",shared_object)
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam CLIKCalcFactory The factory type for creating a CLIK calculator.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename CLIKCalcFactory, typename RateLimitMap>
class manip_rl_inverse_kin_map : public shared_object {
  public:
    
    typedef manip_rl_inverse_kin_map<CLIKCalcFactory,RateLimitMap> self;
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model;
    /** This holds the inverse kinematics calculator factory. */
    CLIKCalcFactory clik_calc_factory;
    mutable kte::manip_clik_calculator ik_calc;
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    manip_rl_inverse_kin_map(const shared_ptr< kte::manipulator_kinematics_model >& aModel = shared_ptr< kte::manipulator_kinematics_model >(),
			     const shared_ptr< RateLimitMap >& aJointLimitMap = shared_ptr< RateLimitMap >(),
                             const CLIKCalcFactory& aCLIKCalcFactory = CLIKCalcFactory()) :
                             model(aModel),
                             clik_calc_factory(aCLIKCalcFactory), ik_calc(aModel.get()),
                             joint_limits_map(aJointLimitMap) {
      ik_calc = clik_calc_factory.create_calculator(*joint_limits_map,model);
    };
    
    /**
     * This function template performs a inverse kinematics calculation on the 
     * manipulator model.
     * \tparam PointType The point-type of the input space.
     * \tparam InSpace The type of the input space (end-effector space).
     * \tparam OutSpace The type of the output space (rate-limited joint-space).
     * \param pt The point in the input space, i.e. the end-effector coordinates.
     * \param space_in The input space, i.e. the end-effector space.
     * \param space_out The output space, i.e. the rate-limited joint-space.
     * \return A point in the output space, i.e. the rate-limited joint coordinates.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::write_dependent_coordinates_impl(pt,space_in,model);
      
      ik_calc.solveInverseKinematics();
    
      typedef typename RateLimitMap::normal_space_type NormalJointSpace;
      typename topology_traits<NormalJointSpace>::point_type result_inter;
      detail::read_joint_coordinates_impl< typename boost::mpl::prior< arithmetic_tuple_size< NormalJointSpace > >::type >(result_inter, model);
      result = joint_limits_map->map_to_space(result_inter, NormalJointSpace(), space_out);
      
      return result;
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(model)
        & RK_SERIAL_SAVE_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_SAVE_WITH_NAME(joint_limits_map);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(model)
        & RK_SERIAL_LOAD_WITH_NAME(clik_calc_factory)
        & RK_SERIAL_LOAD_WITH_NAME(joint_limits_map);
      ik_calc = clik_calc_factory.create_calculator(*joint_limits_map,model);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2400015,1,"manip_rl_inverse_kin_map",shared_object)
    
    
};



};


};

#endif








