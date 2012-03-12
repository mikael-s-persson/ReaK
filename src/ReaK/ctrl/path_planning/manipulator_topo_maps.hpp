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
    write_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size<InSpace> > >(pt, space_in, gen_i, f2d_i, f3d_i, model);
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
    write_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
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
  void write_joint_coordinates_impl( const PointType& pt,
				     const InSpaceTuple& space_in,
				     const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size<InSpaceTuple> > >(pt,space_in,gen_i,f2d_i,f3d_i,model);
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
    set_frame_3D(pt,*(model->Frames2D()[f2d_i++]));
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
    read_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size<InSpace> > >(pt, space_in, gen_i, f2d_i, f3d_i, model);
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
    read_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
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
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  void read_joint_coordinates_impl( PointType& pt,
				    const InSpaceTuple& space_in,
				    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_joint_coordinates_impl< Idx >(pt,space_in,gen_i,f2d_i,f3d_i,model);
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
    write_dependent_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size<InSpace> > >(pt, space_in, gen_i, f2d_i, f3d_i, model);
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
    write_dependent_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
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
  
  template <typename Idx, typename PointType, typename InSpaceTuple>
  void write_dependent_coordinates_impl( const PointType& pt,
				         const InSpaceTuple& space_in,
				         const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_dependent_coordinates_impl< Idx >(pt,space_in,gen_i,f2d_i,f3d_i,model);
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
    set_frame_3D(pt,*(model->DependentFrames2D()[f2d_i++]->mFrame));
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
    read_dependent_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size<InSpace> > >(pt, space_in, gen_i, f2d_i, f3d_i, model);
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
    read_dependent_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,gen_i,f2d_i,f3d_i,model);
    
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
  void read_dependent_coordinates_impl( const PointType& pt,
				        const InSpaceTuple& space_in,
				        const shared_ptr< kte::manipulator_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_dependent_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size< InSpaceTuple > > >(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
  
  
  template <typename InSpace, typename OutSpace>
  void compute_serial_DK_impl( typename topology_traits<OutSpace>::point_type& result,
			       const typename topology_traits<InSpace>::point_type& pt,
			       const InSpace& space_in,
			       const OutSpace& space_out,
			       const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl(pt, space_in, model);
    
    model->doMotion();
    
    read_dependent_coordinates_impl(result, space_out, model);
  };
  
  
  template <typename InSpace, typename OutSpace, typename RateLimitMap>
  void compute_serial_DK_impl( typename topology_traits<OutSpace>::point_type& result,
			       const typename topology_traits<InSpace>::point_type& pt,
			       const InSpace& space_in,
			       const OutSpace& space_out,
			       const RateLimitMap& j_limits,
			       const shared_ptr< kte::manipulator_kinematics_model >& model) {
    typedef typename get_rate_illimited_space< InSpace >::type NormalJointSpace;
    NormalJointSpace normal_j_space;
    typename topology_traits<NormalJointSpace>::point_type pt_inter = j_limits.map_to_space(pt, space_in, normal_j_space);
    write_joint_coordinates_impl(pt_inter, normal_j_space, model);
    
    model->doMotion();
    
    read_dependent_coordinates_impl(result, space_out, model);
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[gen_i] = get<0>(space_out).origin() - get<0>(space_out).get_radius();
    ik_calc.upper_bounds[gen_i] = get<0>(space_out).origin() + get<0>(space_out).get_radius();
    ik_calc.lower_bounds[offset + gen_i] = -std::numeric_limits<double>::infinity();
    ik_calc.upper_bounds[offset + gen_i] = std::numeric_limits<double>::infinity();
    
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[gen_i] = get<0>(space_out).origin() - get<0>(space_out).get_radius();
    ik_calc.upper_bounds[gen_i] = get<0>(space_out).origin() + get<0>(space_out).get_radius();
    ik_calc.lower_bounds[offset + gen_i] = get<1>(space_out).origin() - get<1>(space_out).get_radius();
    ik_calc.upper_bounds[offset + gen_i] = get<1>(space_out).origin() + get<1>(space_out).get_radius();
    
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
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
    
    centers[offset + model->Coords().size() + 3 * f2d_i] = v_o[0];
    centers[offset + model->Coords().size() + 3 * f2d_i + 1] = v_o[1];
    
    centers[offset + model->Coords().size() + 3 * f2d_i + 2] = get<1>(get<1>(space_out)).origin();
    
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t&, std::size_t& f3d_i,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i] = p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 1] = p_lo[1];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 2] = p_lo[2];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i] = p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 1] = p_up[1];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 2] = p_up[2];
    
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    
    centers[offset + model->Coords().size() + 6 * f3d_i] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 1] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 2] = 0.0;
    
    centers[offset + model->Coords().size() + 6 * f3d_i + 3] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 4] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 5] = 0.0;
    
    diag_gains[model->Coords().size() + 7 * f3d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 2] = 4.0 / ((p_up[2] - p_lo[2]) * (p_up[2] - p_lo[2])); 
    
    diag_gains[model->Coords().size() + 7 * f3d_i + 3] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 4] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 5] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 6] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 6 * f3d_i] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 1] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 2] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 3] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 4] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 5] = 1.0; 
    
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i] = p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 1] = p_lo[1];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 2] = p_lo[2];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i] = p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 1] = p_up[1];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 2] = p_up[2];
    
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i] = v_o[0] - v_r;
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = v_o[1] - v_r;
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = v_o[2] - v_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i] = v_o[0] + v_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = v_o[1] + v_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = v_o[2] + v_r;
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = w_o[0] - w_r;
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = w_o[1] - w_r;
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = w_o[2] - w_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = w_o[0] + w_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = w_o[1] + w_r;
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = w_o[2] + w_r;
    
    centers[offset + model->Coords().size() + 6 * f3d_i] = v_o[0];
    centers[offset + model->Coords().size() + 6 * f3d_i + 1] = v_o[1];
    centers[offset + model->Coords().size() + 6 * f3d_i + 2] = v_o[2];
    
    centers[offset + model->Coords().size() + 6 * f3d_i + 3] = w_o[0];
    centers[offset + model->Coords().size() + 6 * f3d_i + 4] = w_o[1];
    centers[offset + model->Coords().size() + 6 * f3d_i + 5] = w_o[2];
    
    diag_gains[model->Coords().size() + 7 * f3d_i] = 4.0 / ((p_up[0] - p_lo[0]) * (p_up[0] - p_lo[0])); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 1] = 4.0 / ((p_up[1] - p_lo[1]) * (p_up[1] - p_lo[1])); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 2] = 4.0 / ((p_up[2] - p_lo[2]) * (p_up[2] - p_lo[2])); 
    
    diag_gains[model->Coords().size() + 7 * f3d_i + 3] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 4] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 5] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 6] = 1.0; 
    
    v_r = 1.0 / (v_r * v_r);
    diag_gains[offset + model->Coords().size() + 6 * f3d_i] = v_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 1] = v_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 2] = v_r; 
    
    w_r = 1.0 / (w_r * w_r);
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 3] = w_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 4] = w_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 5] = w_r; 
    
    ++f3d_i;
    
  };
  
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size<OutSpace> > >(space_out, offset, gen_i, f2d_i, f3d_i, ik_calc, model, centers, diag_gains);
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
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior<Idx> >(space_out,offset,gen_i,f2d_i,f3d_i,ik_calc,model,centers,diag_gains);
    
    write_one_joint_bound_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,ik_calc,model,centers,diag_gains);
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
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_one_joint_bound_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,ik_calc,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple>
  void write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > > >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,ik_calc,model,centers,diag_gains);
  };
  
  
  
  template <typename InSpace, typename OutSpace>
  void compute_serial_IK_impl( typename topology_traits<OutSpace>::point_type& result,
			       const typename topology_traits<InSpace>::point_type& pt,
			       const InSpace& space_in,
			       const OutSpace& space_out,
			       const shared_ptr< kte::manipulator_kinematics_model >& model,
                               const vect_n<double>& preferred_posture) {
    
    write_dependent_coordinates_impl(pt,space_in,model);
    
    kte::manip_clik_calculator ik_calc = kte::manip_clik_calculator(model.get());
    
    vect_n<double> centers(model->getJointPositionsCount() + model->getJointVelocitiesCount());
    vect_n<double> diag_gains(model->getJointPositionsCount() + model->getJointVelocitiesCount());
    std::copy(preferred_posture.begin(), preferred_posture.end(), centers.begin());
    
    write_joint_bounds_impl(
      space_out, 
      ik_calc,
      model,
      centers,
      diag_gains
    );
    
    kte::quadratic_cost_evaluator e(centers, mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(diag_gains)));
    ik_calc.cost_evaluator = &e;
    
    ik_calc.solveInverseKinematics();
    
    read_joint_coordinates_impl(result,space_out,model);
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() - get<0>(space_out).get_radius());
    ik_calc.upper_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() + get<0>(space_out).get_radius());
    ik_calc.lower_bounds[offset + gen_i] = -j_limits.gen_speed_limits[gen_i];
    ik_calc.upper_bounds[offset + gen_i] = j_limits.gen_speed_limits[gen_i];
    
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t& gen_i, std::size_t&, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() - get<0>(space_out).get_radius());
    ik_calc.upper_bounds[gen_i] = j_limits.gen_speed_limits[gen_i] * (get<0>(space_out).origin() + get<0>(space_out).get_radius());
    ik_calc.lower_bounds[offset + gen_i] = j_limits.gen_accel_limits[gen_i] * (get<1>(space_out).origin() - get<1>(space_out).get_radius());
    ik_calc.upper_bounds[offset + gen_i] = j_limits.gen_accel_limits[gen_i] * (get<1>(space_out).origin() + get<1>(space_out).get_radius());
    
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
  void >::type write_one_joint_bound_impl( const OutSpace& space_out,
					   std::size_t offset,
					   std::size_t&, std::size_t& f2d_i, std::size_t&,
					   const RateLimitMap& j_limits,
				           kte::manip_clik_calculator& ik_calc,
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
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
    
    centers[offset + model->Coords().size() + 3 * f2d_i] = 0.0;
    centers[offset + model->Coords().size() + 3 * f2d_i + 1] = 0.0;
    
    centers[offset + model->Coords().size() + 3 * f2d_i + 2] = 0.0;
    
    diag_gains[model->Coords().size() + 4 * f2d_i] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i] * j_limits.frame2D_speed_limits[2 * f2d_i]); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 1] = diag_gains[model->Coords().size() + 4 * f2d_i]; 
    
    diag_gains[model->Coords().size() + 4 * f2d_i + 2] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 3] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i] = 1.0; 
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 1] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 2] = 1.0; 
    
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
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
    
    centers[offset + model->Coords().size() + 3 * f2d_i] = j_limits.frame2D_accel_limits[2 * f2d_i] * v_o[0];
    centers[offset + model->Coords().size() + 3 * f2d_i + 1] = j_limits.frame2D_accel_limits[2 * f2d_i] * v_o[1];
    
    centers[offset + model->Coords().size() + 3 * f2d_i + 2] = j_limits.frame2D_accel_limits[2 * f2d_i + 1] * get<1>(get<1>(space_out)).origin();
    
    diag_gains[model->Coords().size() + 4 * f2d_i] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i] * j_limits.frame2D_speed_limits[2 * f2d_i]); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 1] = diag_gains[model->Coords().size() + 4 * f2d_i]; 
    
    diag_gains[model->Coords().size() + 4 * f2d_i + 2] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    diag_gains[model->Coords().size() + 4 * f2d_i + 3] = 1.0 / (j_limits.frame2D_speed_limits[2 * f2d_i + 1] * j_limits.frame2D_speed_limits[2 * f2d_i + 1]); 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i] = 1.0 / (j_limits.frame2D_accel_limits[2 * f2d_i] * j_limits.frame2D_accel_limits[2 * f2d_i]); 
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 1] = diag_gains[offset + model->Coords().size() + 3 * f2d_i]; 
    
    diag_gains[offset + model->Coords().size() + 3 * f2d_i + 2] = 1.0 / (j_limits.frame2D_accel_limits[2 * f2d_i + 1] * j_limits.frame2D_accel_limits[2 * f2d_i + 1]); 
    
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[1];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[2];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[1];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[2];
    
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = -j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i];
    
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = -j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = j_limits.frame3D_speed_limits[2 * f3d_i + 1];
    
    centers[offset + model->Coords().size() + 6 * f3d_i] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 1] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 2] = 0.0;
    
    centers[offset + model->Coords().size() + 6 * f3d_i + 3] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 4] = 0.0;
    centers[offset + model->Coords().size() + 6 * f3d_i + 5] = 0.0;
    
    diag_gains[model->Coords().size() + 7 * f3d_i] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i] * j_limits.frame3D_speed_limits[2 * f3d_i]); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 1] = diag_gains[model->Coords().size() + 7 * f3d_i]; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 2] = diag_gains[model->Coords().size() + 7 * f3d_i]; 
    
    diag_gains[model->Coords().size() + 7 * f3d_i + 3] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i + 1] * j_limits.frame3D_speed_limits[2 * f3d_i + 1]); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 4] = diag_gains[model->Coords().size() + 7 * f3d_i + 3]; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 5] = diag_gains[model->Coords().size() + 7 * f3d_i + 3]; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 6] = diag_gains[model->Coords().size() + 7 * f3d_i + 3]; 
    
    diag_gains[offset + model->Coords().size() + 6 * f3d_i] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 1] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 2] = 1.0; 
    
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 3] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 4] = 1.0; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 5] = 1.0; 
    
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    
    vect<double,3> p_lo = get<0>(get<0>(space_out)).get_lower_bound();
    vect<double,3> p_up = get<0>(get<0>(space_out)).get_upper_bound();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[0];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[1];
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_lo[2];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[0];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 1] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[1];
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 2] = j_limits.frame3D_speed_limits[2 * f3d_i] * p_up[2];
    
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 3] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 4] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 5] = -std::numeric_limits< double >::infinity();
    ik_calc.lower_bounds[model->Coords().size() + 7 * f3d_i + 6] = -std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 3] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 4] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 5] = std::numeric_limits< double >::infinity();
    ik_calc.upper_bounds[model->Coords().size() + 7 * f3d_i + 6] = std::numeric_limits< double >::infinity();
    
    vect<double,3> v_o = get<1>(get<0>(space_out)).origin(); 
    double v_r = get<1>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[0] - v_r);
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[1] - v_r);
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[2] - v_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[0] + v_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[1] + v_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * (v_o[2] + v_r);
    
    vect<double,3> w_o = get<1>(get<1>(space_out)).origin();
    double w_r = get<1>(get<1>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[0] - w_r);
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[1] - w_r);
    ik_calc.lower_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[2] - w_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[0] + w_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[1] + w_r);
    ik_calc.upper_bounds[offset + model->Coords().size() + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * (w_o[2] + w_r);
    
    centers[offset + model->Coords().size() + 6 * f3d_i] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[0];
    centers[offset + model->Coords().size() + 6 * f3d_i + 1] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[1];
    centers[offset + model->Coords().size() + 6 * f3d_i + 2] = j_limits.frame3D_accel_limits[2 * f3d_i] * v_o[2];
    
    centers[offset + model->Coords().size() + 6 * f3d_i + 3] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[0];
    centers[offset + model->Coords().size() + 6 * f3d_i + 4] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[1];
    centers[offset + model->Coords().size() + 6 * f3d_i + 5] = j_limits.frame3D_accel_limits[2 * f3d_i + 1] * w_o[2];
    
    diag_gains[model->Coords().size() + 7 * f3d_i] = 1.0 / (j_limits.frame3D_speed_limits[2 * f3d_i] * j_limits.frame3D_speed_limits[2 * f3d_i]); 
    diag_gains[model->Coords().size() + 7 * f3d_i + 1] = diag_gains[model->Coords().size() + 7 * f3d_i]; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 2] = diag_gains[model->Coords().size() + 7 * f3d_i]; 
    
    diag_gains[model->Coords().size() + 7 * f3d_i + 3] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 4] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 5] = 1.0; 
    diag_gains[model->Coords().size() + 7 * f3d_i + 6] = 1.0; 
    
    v_r = 1.0 / (j_limits.frame3D_accel_limits[2 * f3d_i] * j_limits.frame3D_accel_limits[2 * f3d_i]);
    diag_gains[offset + model->Coords().size() + 6 * f3d_i] = v_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 1] = v_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 2] = v_r; 
    
    w_r = 1.0 / (j_limits.frame3D_accel_limits[2 * f3d_i + 1] * j_limits.frame3D_accel_limits[2 * f3d_i + 1]);
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 3] = w_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 4] = w_r; 
    diag_gains[offset + model->Coords().size() + 6 * f3d_i + 5] = w_r; 
    
    ++f3d_i;
    
  };
  
  
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
			                   const shared_ptr< kte::manipulator_kinematics_model >& model,
				           vect_n<double>& centers, vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpace > > >(space_out, offset, gen_i, f2d_i, f3d_i, j_limits, ik_calc, model, centers, diag_gains);
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
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior<Idx> >(space_out,offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model,centers,diag_gains);
    
    write_one_joint_bound_impl(get<Idx::type::value>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model,centers,diag_gains);
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
					const shared_ptr< kte::manipulator_kinematics_model >& model,
					vect_n<double>& centers,
					vect_n<double>& diag_gains) {
    write_one_joint_bound_impl(get<0>(space_out),offset,gen_i,f2d_i,f3d_i,j_limits,ik_calc,model,centers,diag_gains);
  };
  
  template <typename OutSpaceTuple, typename RateLimitMap>
  void write_joint_bounds_impl( const OutSpaceTuple& space_out,
			        const RateLimitMap& j_limits,
			        kte::manip_clik_calculator& ik_calc,
			        const shared_ptr< kte::manipulator_kinematics_model >& model,
			        vect_n<double>& centers,
			        vect_n<double>& diag_gains) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpaceTuple > > >(space_out,model->getJointPositionsCount(),gen_i,f2d_i,f3d_i,j_limits,ik_calc,model,centers,diag_gains);
  };
  
  
  
  
  
  template <typename InSpace, typename OutSpace, typename RateLimitMap>
  void compute_serial_IK_impl( typename topology_traits<OutSpace>::point_type& result,
			       const typename topology_traits<InSpace>::point_type& pt,
			       const InSpace& space_in,
			       const OutSpace& space_out,
			       const RateLimitMap& j_limits,
			       const shared_ptr< kte::manipulator_kinematics_model >& model,
                               const vect_n<double>& preferred_posture) {
    
    write_dependent_coordinates_impl(pt,space_in,model);
    
    kte::manip_clik_calculator ik_calc = kte::manip_clik_calculator(model.get());
    
    vect_n<double> centers(model->getJointPositionsCount() + model->getJointVelocitiesCount());
    vect_n<double> diag_gains(model->getJointPositionsCount() + model->getJointVelocitiesCount());
    std::copy(preferred_posture.begin(), preferred_posture.end(), centers.begin());
    
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpace > > >(
      space_out, 
      j_limits,
      ik_calc,
      centers,
      diag_gains
    );
    
    kte::quadratic_cost_evaluator e(centers, mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(diag_gains)));
    ik_calc.cost_evaluator = &e;
    
    ik_calc.solveInverseKinematics();
    
    typedef typename RateLimitMap::normal_space_type NormalJointSpace;
    typename topology_traits<NormalJointSpace>::point_type result_inter;
    read_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size< NormalJointSpace > > >(result_inter, model);
    result = j_limits.map_to_space(result_inter, NormalJointSpace(), space_out);
  };
  
  
  
  
};





/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates (both 
 * generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_direct_kinematics_map {
  public:
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    
    /**
     * This function template performs a forward kinematics calculation on the 
     * manipulator model.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::compute_serial_DK_impl(result, pt, space_in, space_out, model);
      
      return result;
    };
    
};



/**
 * This class implements the forward kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates 
 * (both generalized and frames), and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename RateLimitMap>
class manip_rl_direct_kinematics_map {
  public:
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    /**
     * This function template performs a forward kinematics calculation on the 
     * manipulator model.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::compute_serial_DK_impl(result, pt, space_in, space_out, *joint_limits_map, model);
      
      return result;
    };
    
};



/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of joint coordinates, and that 
 * it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 */
class manip_inverse_kinematics_map {
  public:
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model; 
    /** This holds the posture that is the most desirable for the manipulator, that is, a set of 
     * joint coordinates (in normal joint-space) that is well conditioned w.r.t. some performance
     * index like manipulability or condition number of the Jacobian. This is used in the inverse
     * kinematics mapping to optimize the condition of the posture of the robot, that is, if there 
     * is any freedom in the matter (if the manipulator has some redundancy to exploit) otherwise 
     * it has no effect. */
    vect_n<double> preferred_posture;
    
    /**
     * This function template performs a forward or inverse kinematics calculation on the manipulator model.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::compute_serial_IK_impl(result, pt, space_in, space_out, model, preferred_posture);
      
      return result;
    };
    
};


/**
 * This class implements the inverse kinematics mappings associated to a given manipulator kinematics 
 * model. This class assumes that the manipulator model has a number of rate-limited joint coordinates, 
 * and that it has dependent coordinate frames (gen, 2D or 3D) as end-effectors.
 * \tparam RateLimitMap The type of the mapping between rate-limited joint-spaces and normal joint-spaces.
 */
template <typename RateLimitMap>
class manip_rl_inverse_kinematics_map {
  public:
    
    /** This data member points to a manipulator kinematics model to use for the mappings performed. */
    shared_ptr< kte::manipulator_kinematics_model > model;
    /** This holds the posture that is the most desirable for the manipulator, that is, a set of 
     * joint coordinates (in normal joint-space) that is well conditioned w.r.t. some performance
     * index like manipulability or condition number of the Jacobian. This is used in the inverse
     * kinematics mapping to optimize the condition of the posture of the robot, that is, if there 
     * is any freedom in the matter (if the manipulator has some redundancy to exploit) otherwise 
     * it has no effect. */
    vect_n<double> preferred_posture;
    /** This data member holds a mapping between the rate-limited joint space and the normal joint-space. */
    shared_ptr< RateLimitMap > joint_limits_map;
    
    /**
     * This function template performs a forward kinematics calculation on the manipulator model.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::compute_serial_IK_impl(result, pt, space_in, space_out, *joint_limits_map, model, preferred_posture);
      
      return result;
    };
    
};



};


};

#endif








