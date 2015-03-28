/**
 * \file manip_SCARA_arm.hpp
 *
 * This library declares a class to represent a kte-based model of a 3R manipulator in 2D or 3D.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
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

#ifndef REAK_MANIP_SCARA_ARM_HPP
#define REAK_MANIP_SCARA_ARM_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/mbd/kte/kte_map_chain.hpp>

#include "inverse_kinematics_model.hpp"

namespace ReaK {

namespace kte {


/**
 * This class that models a 3D manipulator with 2 revolute joints and 1 prismatic joint
 * arranged in a SCARA configuration (i.e., shoulder-elbow joints and prismatic wrist,
 * all aligned along the z-axis). This class is only a kinematics model.
 */
class manip_SCARA_kinematics : public inverse_kinematics_model {
private:
  shared_ptr< frame_3D< double > > m_base_frame;
  std::vector< shared_ptr< gen_coord< double > > > m_joints;
  shared_ptr< joint_dependent_frame_3D > m_EE;
  double m_link1_length;
  double m_link2_length;
  double m_link3_height;

  shared_ptr< kte_map_chain > m_chain;

public:
  BOOST_STATIC_CONSTANT( std::size_t, degrees_of_freedom = 3 );

  /**
   * Default constructor.
   */
  manip_SCARA_kinematics( const std::string& aName = "",
                          const shared_ptr< frame_3D< double > >& aBaseFrame = shared_ptr< frame_3D< double > >(),
                          double aLink1Length = 1.0, double aLink2Length = 1.0, double aLink3Height = 0.0 );

  virtual ~manip_SCARA_kinematics(){};

  virtual std::size_t getJointPositionsCount() const { return 3; };

  virtual std::size_t getJointVelocitiesCount() const { return 3; };

  virtual std::size_t getJointAccelerationsCount() const { return 3; };

  virtual std::size_t getDependentPositionsCount() const { return 7; };

  virtual std::size_t getDependentVelocitiesCount() const { return 6; };

  virtual std::size_t getDependentAccelerationsCount() const { return 6; };

  virtual std::size_t getCoordsCount() const { return 3; };

  virtual shared_ptr< gen_coord< double > > getCoord( std::size_t i ) const { return m_joints[i]; };

  virtual std::size_t getDependentFrames3DCount() const { return 1; };

  virtual shared_ptr< joint_dependent_frame_3D > getDependentFrame3D( std::size_t i ) const { return m_EE; };

  virtual void doDirectMotion();

  virtual void doInverseMotion();

  virtual void getJacobianMatrix( mat< double, mat_structure::rectangular >& Jac ) const;

  virtual void getJacobianMatrixAndDerivative( mat< double, mat_structure::rectangular >& Jac,
                                               mat< double, mat_structure::rectangular >& JacDot ) const;

  virtual vect_n< double > getJointPositionLowerBounds() const {
    return vect_n< double >( -std::numeric_limits< double >::infinity(), -std::numeric_limits< double >::infinity(),
                             -std::numeric_limits< double >::infinity() );
  };

  virtual void setJointPositionLowerBounds( const vect_n< double >& aJointLowerBounds ){};

  virtual vect_n< double > getJointPositionUpperBounds() const {
    return vect_n< double >( std::numeric_limits< double >::infinity(), std::numeric_limits< double >::infinity(),
                             std::numeric_limits< double >::infinity() );
  };

  virtual void setJointPositionUpperBounds( const vect_n< double >& aJointUpperBounds ){};

  virtual vect_n< double > getJointPositions() const;

  virtual void setJointPositions( const vect_n< double >& aJointPositions );

  virtual vect_n< double > getJointVelocities() const;

  virtual void setJointVelocities( const vect_n< double >& aJointVelocities );

  virtual vect_n< double > getJointAccelerations() const;

  virtual void setJointAccelerations( const vect_n< double >& aJointAccelerations );

  virtual vect_n< double > getDependentPositions() const;

  virtual vect_n< double > getDependentVelocities() const;

  virtual vect_n< double > getDependentAccelerations() const;

  virtual void setDependentPositions( const vect_n< double >& aDepPositions );

  virtual void setDependentVelocities( const vect_n< double >& aDepVelocities );

  virtual void setDependentAccelerations( const vect_n< double >& aDepAccelerations );

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const;

  virtual void RK_CALL load( serialization::iarchive& A, unsigned int );

  RK_RTTI_MAKE_CONCRETE_1BASE( manip_SCARA_kinematics, 0xC2100054, 1, "manip_SCARA_kinematics",
                               inverse_kinematics_model )
};
};
};

#endif
