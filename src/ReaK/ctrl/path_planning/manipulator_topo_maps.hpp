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

#include "lin_alg/arithmetic_tuple.hpp"
#include "lin_alg/vect_alg.hpp"

#include "kinetostatics/gen_coord.hpp"

#include "mbd_kte/manipulator_model.hpp"
#include "mbd_kte/manipulator_model_helper.hpp"

#include <boost/mpl/less.hpp>

namespace ReaK {

namespace pp {


namespace detail {
  
  
  template <typename Idx, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const typename topology_traits<InSpace>::point_type& pt,
				             const InSpace& space_in,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,model);
    
    *(model->Coords()[Idx::type::value]) = get_gen_coord(get<Idx::type::value>(pt));
    
  };
  
  template <typename Idx, typename InSpace>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const typename topology_traits<InSpace>::point_type& pt,
				             const InSpace& space_in,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->Coords()[0]) = get_gen_coord(get<0>(pt));
  };
  
  
  template <typename Idx, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( typename topology_traits<InSpace>::point_type& pt,
				            const InSpace& space_in,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,model);
    
    set_gen_coord(get<Idx::type::value>(pt),*(model->Coords()[Idx::type::value]));
  };
  
  template <typename Idx, typename InSpace>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( typename topology_traits<InSpace>::point_type& pt,
				            const InSpace& space_in,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_gen_coord(get<0>(pt),*(model->Coords()[0]));
  };
  
  
  
  template <typename Idx, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const typename topology_traits<InSpace>::point_type& pt,
				             const InSpace& space_in,
					     const RateLimitMap& j_limits,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,model);
    
    *(model->Coords()[Idx::type::value]) = get_gen_coord(get<Idx::type::value>(pt));
    
  };
  
  template <typename Idx, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_coordinates_impl( const typename topology_traits<InSpace>::point_type& pt,
				             const InSpace& space_in,
					     const RateLimitMap& j_limits,
				             const shared_ptr< kte::manipulator_kinematics_model >& model) {
    *(model->Coords()[0]) = get_gen_coord(get<0>(pt));
  };
  
  
  template <typename Idx, typename InSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( typename topology_traits<InSpace>::point_type& pt,
				            const InSpace& space_in,
					    const RateLimitMap& j_limits,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    read_joint_coordinates_impl< boost::mpl::prior<Idx> >(pt,space_in,model);
    
    set_gen_coord(get<Idx::type::value>(pt),*(model->Coords()[Idx::type::value]));
  };
  
  template <typename Idx, typename InSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type read_joint_coordinates_impl( typename topology_traits<InSpace>::point_type& pt,
				            const InSpace& space_in,
					    const RateLimitMap& j_limits,
				            const shared_ptr< kte::manipulator_kinematics_model >& model) {
    set_gen_coord(get<0>(pt),*(model->Coords()[0]));
  };
  // TODO: Find a way to deal with joints that are not expressed in the correct units (reach-time joint-spaces!).
  
  
  
  
  template <typename InSpace, typename OutSpace>
  void compute_DK_3dEE_impl( typename topology_traits<OutSpace>::point_type& result,
			     const typename topology_traits<InSpace>::point_type& pt,
			     const InSpace& space_in,
			     const OutSpace&,
			     const shared_ptr< kte::manipulator_kinematics_model >& model) {
    write_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size< space_in > > >(pt, space_in, model);
    
    model->doMotion();
    
    set_frame_3D(result, *(model->DependentFrames3D()[0]->mFrame));
  };
  
  
  template <typename InSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< hyperbox_topology< vect<double,3> > >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< quaternion_topology<double> >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type& result,
    const typename topology_traits<
      InSpaceTuple
    >::point_type& pt,
    const InSpaceTuple& space_in,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< hyperbox_topology< vect<double,3> > >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< quaternion_topology<double> >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_out,
	  const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_DK_3dEE_impl(result, pt, space_in, space_out,model);
  };
  
  template <typename InSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    hyperbox_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> > 
	  >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    quaternion_topology<double>,
	    ang_velocity_3D_topology<double>
	  >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type& result,
    const typename topology_traits<
      InSpaceTuple
    >::point_type& pt,
    const InSpaceTuple& space_in,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        hyperbox_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> > 
	      >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        quaternion_topology<double>,
	        ang_velocity_3D_topology<double>
	      >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_out,
	  const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_DK_3dEE_impl(result, pt, space_in, space_out, model);
  };
  
  template <typename InSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    hyperbox_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> > 
	  >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    quaternion_topology<double>,
            ang_velocity_3D_topology<double>,
	    ang_accel_3D_topology<double>
	  >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type& result,
    const typename topology_traits<
      InSpaceTuple
    >::point_type& pt,
    const InSpaceTuple& space_in,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        hyperbox_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> > 
	      >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        quaternion_topology<double>,
	        ang_velocity_3D_topology<double>,
	        ang_accel_3D_topology<double>
	      >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_out,
	  const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_DK_3dEE_impl(result, pt, space_in, space_out, model);
  };
  
  
  
  
  
  
  
  
  
  template <typename Idx, typename OutSpace>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpace& space_out,
					std::size_t offset,
				        kte::manip_clik_calculator& ik_calc,
				        vect_n<double>& centers,
				        vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior<Idx> >(space_in,offset,ik_calc);
    
    ik_calc.lower_bounds[Idx::type::value] = get<0>(get<Idx::type::value>(space_out)).origin() - get<0>(get<Idx::type::value>(space_out)).get_radius();
    ik_calc.upper_bounds[Idx::type::value] = get<0>(get<Idx::type::value>(space_out)).origin() + get<0>(get<Idx::type::value>(space_out)).get_radius();
    ik_calc.lower_bounds[offset + Idx::type::value] = get<1>(get<Idx::type::value>(space_out)).origin() - get<1>(get<Idx::type::value>(space_out)).get_radius();
    ik_calc.upper_bounds[offset + Idx::type::value] = get<1>(get<Idx::type::value>(space_out)).origin() + get<1>(get<Idx::type::value>(space_out)).get_radius();
    
    centers[offset + Idx::type::value] = get<1>(get<Idx::type::value>(space_out)).origin();
    diag_gains[Idx::type::value] = get<0>(get<Idx::type::value>(space_out)).get_radius(); diag_gains[Idx::type::value] = 1.0 / (diag_gains[Idx::type::value] * diag_gains[Idx::type::value]);
    diag_gains[offset + Idx::type::value] = get<1>(get<Idx::type::value>(space_out)).get_radius(); diag_gains[offset + Idx::type::value] = 1.0 / (diag_gains[offset + Idx::type::value] * diag_gains[offset + Idx::type::value]);
  };
  
  template <typename Idx, typename OutSpace>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpace& space_out,
					std::size_t offset,
				        kte::manip_clik_calculator& ik_calc,
				        vect_n<double>& centers,
				        vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[0] = get<0>(get<0>(space_out)).origin() - get<0>(get<0>(space_out)).get_radius();
    ik_calc.upper_bounds[0] = get<0>(get<0>(space_out)).origin() + get<0>(get<0>(space_out)).get_radius();
    ik_calc.lower_bounds[offset] = get<1>(get<0>(space_out)).origin() - get<1>(get<0>(space_out)).get_radius();
    ik_calc.upper_bounds[offset] = get<1>(get<0>(space_out)).origin() + get<1>(get<0>(space_out)).get_radius();
    
    centers[offset] = get<1>(get<0>(space_out)).origin();
    diag_gains[0] = get<0>(get<0>(space_out)).get_radius(); diag_gains[0] = 1.0 / (diag_gains[0] * diag_gains[0]);
    diag_gains[offset] = get<1>(get<0>(space_out)).get_radius(); diag_gains[offset] = 1.0 / (diag_gains[offset] * diag_gains[offset]);
  };
  
  
  template <typename InSpace, typename OutSpace>
  void compute_IK_3dEE_impl( typename topology_traits<OutSpace>::point_type& result,
			     const typename topology_traits<InSpace>::point_type& pt,
			     const InSpace& space_in,
			     const OutSpace& space_out,
			     const shared_ptr< kte::manipulator_kinematics_model >& model,
			     const vect_n<double>& ideal_posture) {
    
    *(model->DependentFrames3D()[0]->mFrame) = get_frame_3D(pt);
    
    kte::manip_clik_calculator ik_calc = kte::manip_clik_calculator(model.get());
    
    vect_n<double> centers(2 * arithmetic_tuple_size< OutSpace >::type::value);
    vect_n<double> diag_gains(2 * arithmetic_tuple_size< OutSpace >::type::value);
    std::copy(ideal_posture.begin(), ideal_posture.end(), centers.begin());
    
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< space_out > > >(
      space_out, 
      arithmetic_tuple_size< OutSpace >::type::value,
      ik_calc,
      centers,
      diag_gains
    );
    
    kte::quadratic_cost_evaluator e(centers, mat<double,mat_structure::symmetric>(mat<double,mat_structure::diagonal>(diag_gains)));
    ik_calc.cost_evaluator = &e;
    
    ik_calc.solveInverseKinematics();
    
    read_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpace > > >(result, space_out, model);
  };
  
  
  
  
  
  
  template <typename Idx, typename OutSpace, typename RateLimitMap>
  typename boost::disable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpace& space_out,
					std::size_t offset,
					const RateLimitMap& j_limits,
				        kte::manip_clik_calculator& ik_calc,
				        vect_n<double>& centers,
				        vect_n<double>& diag_gains) {
    write_joint_bounds_impl< boost::mpl::prior<Idx> >(space_in,offset,ik_calc);
    
    ik_calc.lower_bounds[Idx::type::value] = j_limits.speed_limits[Idx::type::value] * (get<0>(get<Idx::type::value>(space_out)).origin() - get<0>(get<Idx::type::value>(space_out)).get_radius());
    ik_calc.upper_bounds[Idx::type::value] = j_limits.speed_limits[Idx::type::value] * (get<0>(get<Idx::type::value>(space_out)).origin() + get<0>(get<Idx::type::value>(space_out)).get_radius());
    ik_calc.lower_bounds[offset + Idx::type::value] = j_limits.accel_limits[Idx::type::value] * (get<1>(get<Idx::type::value>(space_out)).origin() - get<1>(get<Idx::type::value>(space_out)).get_radius());
    ik_calc.upper_bounds[offset + Idx::type::value] = j_limits.accel_limits[Idx::type::value] * (get<1>(get<Idx::type::value>(space_out)).origin() + get<1>(get<Idx::type::value>(space_out)).get_radius());
    
    centers[offset + Idx::type::value] = j_limits.accel_limits[Idx::type::value] * get<1>(get<Idx::type::value>(space_out)).origin();
    diag_gains[Idx::type::value] = j_limits.speed_limits[Idx::type::value]; diag_gains[Idx::type::value] = 1.0 / (diag_gains[Idx::type::value] * diag_gains[Idx::type::value]);
    diag_gains[offset + Idx::type::value] = j_limits.accel_limits[Idx::type::value]; diag_gains[offset + Idx::type::value] = 1.0 / (diag_gains[offset + Idx::type::value] * diag_gains[offset + Idx::type::value]);
  };
  
  template <typename Idx, typename OutSpace, typename RateLimitMap>
  typename boost::enable_if< 
    boost::mpl::less< 
      Idx, 
      boost::mpl::size_t<1> 
    >,
  void >::type write_joint_bounds_impl( const OutSpace& space_out,
					std::size_t offset,
					const RateLimitMap& j_limits,
				        kte::manip_clik_calculator& ik_calc,
				        vect_n<double>& centers,
				        vect_n<double>& diag_gains) {
    ik_calc.lower_bounds[0] = j_limits.speed_limits[0] * (get<0>(get<0>(space_out)).origin() - get<0>(get<0>(space_out)).get_radius());
    ik_calc.upper_bounds[0] = j_limits.speed_limits[0] * (get<0>(get<0>(space_out)).origin() + get<0>(get<0>(space_out)).get_radius());
    ik_calc.lower_bounds[offset] = j_limits.accel_limits[0] * (get<1>(get<0>(space_out)).origin() - get<1>(get<0>(space_out)).get_radius());
    ik_calc.upper_bounds[offset] = j_limits.accel_limits[0] * (get<1>(get<0>(space_out)).origin() + get<1>(get<0>(space_out)).get_radius());
    
    centers[offset] = j_limits.accel_limits[0] * get<1>(get<0>(space_out)).origin();
    diag_gains[0] = j_limits.speed_limits[0]; diag_gains[0] = 1.0 / (diag_gains[0] * diag_gains[0]);
    diag_gains[offset] = j_limits.accel_limits[0]; diag_gains[offset] = 1.0 / (diag_gains[offset] * diag_gains[offset]);
  };
  
  
  template <typename InSpace, typename OutSpace, typename RateLimitMap>
  void compute_IK_3dEE_impl( typename topology_traits<OutSpace>::point_type& result,
			     const typename topology_traits<InSpace>::point_type& pt,
			     const InSpace& space_in,
			     const OutSpace& space_out,
			     const RateLimitMap& j_limits,
			     const shared_ptr< kte::manipulator_kinematics_model >& model,
			     const vect_n<double>& ideal_posture) {
    
    *(model->DependentFrames3D()[0]->mFrame) = get_frame_3D(pt);
    
    kte::manip_clik_calculator ik_calc = kte::manip_clik_calculator(model.get());
    
    vect_n<double> centers(2 * arithmetic_tuple_size< OutSpace >::type::value);
    vect_n<double> diag_gains(2 * arithmetic_tuple_size< OutSpace >::type::value);
    std::copy(ideal_posture.begin(), ideal_posture.end(), centers.begin());
    
    write_joint_bounds_impl< boost::mpl::prior< arithmetic_tuple_size< OutSpace > > >(
      space_out, 
      arithmetic_tuple_size< OutSpace >::type::value,
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
    read_joint_coordinates_impl< boost::mpl::prior< arithmetic_tuple_size< NormalJointSpace > > >(result_inter, space_out, model);
    result = j_limits.map_to_space(result_inter, NormalJointSpace(), space_out);
  };
  
  
  
  template <typename OutSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits<
      OutSpaceTuple
    >::point_type& result,
    const typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< hyperbox_topology< vect<double,3> > >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< quaternion_topology<double> >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type& pt,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< hyperbox_topology< vect<double,3> > >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< quaternion_topology<double> >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_in,
    const OutSpaceTuple& space_out,
    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_IK_3dEE_impl(result, pt, space_in, space_out,model);
  };
  
  template <typename OutSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits<
      OutSpaceTuple
    >::point_type& result,
    const typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    hyperbox_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> > 
	  >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    quaternion_topology<double>,
	    ang_velocity_3D_topology<double>
	  >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type & pt,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        hyperbox_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> > 
	      >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        quaternion_topology<double>,
	        ang_velocity_3D_topology<double>
	      >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_in,
    const OutSpaceTuple& space_out,
    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_IK_3dEE_impl(result, pt, space_in, space_out, model);
  };
  
  template <typename OutSpaceTuple, typename DistanceMetric>
  void compute_kinematics_impl(
    typename topology_traits<
      OutSpaceTuple
    >::point_type& result,
    const typename topology_traits< 
      metric_space_tuple< arithmetic_tuple<
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    hyperbox_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> >,
	    hyperball_topology< vect<double,3> > 
	  >, 
	  DistanceMetric 
        >,
        differentiable_space< 
          time_topology, 
	  arithmetic_tuple< 
	    quaternion_topology<double>,
            ang_velocity_3D_topology<double>,
	    ang_accel_3D_topology<double>
	  >, 
	  DistanceMetric 
        > >,
        DistanceMetric 
      >
    >::point_type& pt,
    const metric_space_tuple< arithmetic_tuple<
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        hyperbox_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> >,
	        hyperball_topology< vect<double,3> > 
	      >, 
	      DistanceMetric 
            >,
            differentiable_space< 
              time_topology, 
	      arithmetic_tuple< 
	        quaternion_topology<double>,
	        ang_velocity_3D_topology<double>,
	        ang_accel_3D_topology<double>
	      >, 
	      DistanceMetric 
            > >,
            DistanceMetric 
          >& space_in,
    const OutSpaceTuple& space_out,
    const shared_ptr< kte::manipulator_kinematics_model >& model) {
    compute_IK_3dEE_impl(result, pt, space_in, space_out, model);
  };
  
  
  
  
  
  
};



class manipulator_kinematic_map {
  public:
    
    shared_ptr< kte::manipulator_kinematics_model > model;
    
    
    
    /**
     * This function template performs a forward kinematics calculation on the manipulator model.
     */
    template <typename PointType, typename InSpace, typename OutSpace>
    typename topology_traits< OutSpace >::point_type
    map_to_space(const PointType& pt, const InSpace& space_in, const OutSpace& space_out) const {
      typename topology_traits< OutSpace >::point_type result;
      
      detail::compute_kinematics_impl(result, pt, space_in, space_out, model);
      
    };
    
    
    
    
  
  
  
  
};



};


};

#endif








