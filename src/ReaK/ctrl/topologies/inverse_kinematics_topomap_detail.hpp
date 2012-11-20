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


#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>


namespace ReaK {

namespace pp {


namespace detail {
  
  
  
  
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
                                            const shared_ptr< kte::direct_kinematics_model >& model);
  
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
                                            const shared_ptr< kte::direct_kinematics_model >& model);
  
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
                                          const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                            const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                            const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                    const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                    const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                                 const shared_ptr< kte::direct_kinematics_model >& model);
  
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
                                                 const shared_ptr< kte::direct_kinematics_model >& model);
  
  
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
                                               const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                                 const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                                 const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                         const shared_ptr< kte::direct_kinematics_model >& model) {
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
                                                 const shared_ptr< kte::direct_kinematics_model >& model) {
    std::size_t gen_i = 0;
    std::size_t f2d_i = 0;
    std::size_t f3d_i = 0;
    write_one_dependent_coord_impl(pt,space_in,gen_i,f2d_i,f3d_i,model);
  };
  
  
  
  
};

};

};


#endif




