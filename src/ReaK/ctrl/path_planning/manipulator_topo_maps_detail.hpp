/*
 * \file manipulator_topo_maps_details.hpp
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

};

};







