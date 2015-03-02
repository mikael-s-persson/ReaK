/*
 * \file inverse_kinematics_topomap_detail.hpp
 * 
 * DETAILS OF:
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

#ifndef REAK_INVERSE_KINEMATICS_TOPOMAP_DETAIL_HPP
#define REAK_INVERSE_KINEMATICS_TOPOMAP_DETAIL_HPP

#include <ReaK/math/lin_alg/arithmetic_tuple.hpp>

#include "joint_space_topologies.hpp"
#include "joint_space_limits.hpp"
#include "se3_topologies.hpp"
#include "se2_topologies.hpp"
#include "Ndof_spaces.hpp"

namespace ReaK {

namespace pp {

namespace detail { namespace {
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace&,
                                          std::size_t& gen_i, std::size_t&, std::size_t&,
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
    set_gen_coord(pt,*(model->getCoord(gen_i++)));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace&,
                                          std::size_t&, std::size_t& f2d_i, std::size_t&,
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
    set_frame_2D(pt,*(model->getFrame2D(f2d_i++)));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace&,
                                          std::size_t&, std::size_t&, std::size_t& f3d_i,
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
    set_frame_3D(pt,*(model->getFrame3D(f3d_i++)));
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_Ndof_space<InSpace>,
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace&,
                                          std::size_t& gen_i, std::size_t&, std::size_t&,
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
    for(std::size_t i = 0; i < get<0>(pt).size(); ++i) 
      set_gen_coord(pt, i, *(model->getCoord(gen_i++)));
  };
  
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>,
      is_Ndof_space<InSpace>
    >,  
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace& space_in,
                                          std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                          const shared_ptr< kte::direct_kinematics_model >& model);
  
  struct joint_coordinates_reader_impl {
    shared_ptr< kte::direct_kinematics_model > p_model;
    std::size_t* p_gen_i;
    std::size_t* p_f2d_i;
    std::size_t* p_f3d_i;
    joint_coordinates_reader_impl(const shared_ptr< kte::direct_kinematics_model >& model, 
                                  std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) :
                                  p_model(model), p_gen_i(&gen_i), p_f2d_i(&f2d_i), p_f3d_i(&f3d_i) {};
    template <typename Point, typename Space>
    void operator()(Point& pt, const Space& space) const {
      read_one_joint_coord_impl(pt, space, *p_gen_i, *p_f2d_i, *p_f3d_i, p_model);
    };
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>,
      is_Ndof_space<InSpaceTuple>
    >,  
  void >::type read_joint_coordinates_impl( PointType& pt,
                                            const InSpaceTuple& space_in,
                                            const shared_ptr< kte::direct_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    tuple_for_each(pt, space_in, joint_coordinates_reader_impl(model, gen_i, f2d_i, f3d_i));
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>,
      is_Ndof_space<InSpaceTuple>
    >,  
  void >::type read_joint_coordinates_impl( PointType& pt,
                                    const InSpaceTuple& space_in,
                                    const shared_ptr< kte::direct_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    read_one_joint_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>,
      is_Ndof_space<InSpace>
    >,  
  void >::type read_one_joint_coord_impl( PointType& pt,
                                          const InSpace& space_in,
                                          std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
    tuple_for_each(pt, space_in, joint_coordinates_reader_impl(model, gen_i, f2d_i, f3d_i));
  };
  
  
  
  
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_normal_joint_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace&,
                                               std::size_t& gen_i, std::size_t&, std::size_t&,
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
    *(model->getDependentCoord(gen_i++)->mFrame) = get_gen_coord(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se2_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace&,
                                               std::size_t&, std::size_t& f2d_i, std::size_t&,
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
    *(model->getDependentFrame2D(f2d_i++)->mFrame) = get_frame_2D(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_se3_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace&,
                                               std::size_t&, std::size_t&, std::size_t& f3d_i,
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
    *(model->getDependentFrame3D(f3d_i++)->mFrame) = get_frame_3D(pt);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::enable_if< 
    is_Ndof_space<InSpace>,
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace&,
                                               std::size_t& gen_i, std::size_t&, std::size_t&,
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
    for(std::size_t i = 0; i < get<0>(pt).size(); ++i)
      *(model->getDependentCoord(gen_i++)->mFrame) = get_gen_coord(pt, i);
  };
  
  
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>,
      is_Ndof_space<InSpace>
    >,  
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace& space_in,
                                               std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                               const shared_ptr< kte::direct_kinematics_model >& model);
  
  struct dep_coordinates_writer_impl {
    shared_ptr< kte::direct_kinematics_model > p_model;
    std::size_t* p_gen_i;
    std::size_t* p_f2d_i;
    std::size_t* p_f3d_i;
    dep_coordinates_writer_impl(const shared_ptr< kte::direct_kinematics_model >& model, 
                                std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i) :
                                p_model(model), p_gen_i(&gen_i), p_f2d_i(&f2d_i), p_f3d_i(&f3d_i) {};
    template <typename Point, typename Space>
    void operator()(const Point& pt, const Space& space) const {
      write_one_dependent_coord_impl(pt, space, *p_gen_i, *p_f2d_i, *p_f3d_i, p_model);
    };
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>,
      is_Ndof_space<InSpaceTuple>
    >,  
  void >::type write_dependent_coordinates_impl( const PointType& pt,
                                         const InSpaceTuple& space_in,
                                         const shared_ptr< kte::direct_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    tuple_for_each(pt, space_in, dep_coordinates_writer_impl(model, gen_i, f2d_i, f3d_i));
  };
  
  template <typename PointType, typename InSpaceTuple>
  typename boost::enable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpaceTuple>,
      is_se2_space<InSpaceTuple>,
      is_se3_space<InSpaceTuple>,
      is_Ndof_space<InSpaceTuple>
    >,  
  void >::type write_dependent_coordinates_impl( const PointType& pt,
                                                 const InSpaceTuple& space_in,
                                                 const shared_ptr< kte::direct_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_dependent_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  template <typename PointType, typename InSpace>
  typename boost::disable_if< 
    boost::mpl::or_<
      is_normal_joint_space<InSpace>,
      is_se2_space<InSpace>,
      is_se3_space<InSpace>,
      is_Ndof_space<InSpace>
    >,  
  void >::type write_one_dependent_coord_impl( const PointType& pt,
                                               const InSpace& space_in,
                                               std::size_t& gen_i, std::size_t& f2d_i, std::size_t& f3d_i,
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
    tuple_for_each(pt, space_in, dep_coordinates_writer_impl(model, gen_i, f2d_i, f3d_i));
  };
  
  
  
}; };

};

};


#endif




