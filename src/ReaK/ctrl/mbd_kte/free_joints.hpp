/**
 * \file free_joints.hpp
 * 
 * This library provides classes to handle free-joints. A free joint is simple a 
 * frame transformation that is applied to an input frame to obtain the output 
 * frame. This way, one can represent a general motion of a body.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef FREE_JOINTS_HPP
#define FREE_JOINTS_HPP


#include "reacting_kte.hpp"
#include "math/motion_jacobians.hpp"

#include <boost/shared_ptr.hpp>

namespace ReaK {

namespace kte {


/**
 * This class implements a free joint in 2D space. A 2D frame is used to represent the
 * joint's motion between a base coordinate frame to an end coordinate frame.
 */
class free_joint_2D : public reacting_kte_2D {
  protected:
    boost::shared_ptr< frame_2D<double> > mCoord; ///< The coordinate frame representing the joint's motion.
    boost::shared_ptr< frame_2D<double> > mBase; ///< The coordinate frame at the base of the joint.
    boost::shared_ptr< frame_2D<double> > mEnd; ///< The coordinate frame just after the joint transformations are applied.

    boost::shared_ptr< jacobian_2D_2D<double> > mJacobian; ///< The Jacobian frame produced by this joint.

  public:

    /**
     * Default constructor.
     */
    free_joint_2D(const std::string& aName = "") : reacting_kte_2D(aName), mCoord(), mBase(), mEnd(), mJacobian() { };

    /**
     * Parametrized constructor.
     * \param aName the name of this KTE model.
     * \param aCoord the coordinate frame associated with the displacement of this joint.
     * \param aBase the coordinate frame at the base of the joint.
     * \param aEnd the coordinate frame just after the joint transformations are applied.
     * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the Jacobian frame's calculation.
     */
    free_joint_2D(const std::string& aName,
		       boost::shared_ptr< frame_2D<double> > aCoord,
		       boost::shared_ptr< frame_2D<double> > aBase,
		       boost::shared_ptr< frame_2D<double> > aEnd,
                       boost::shared_ptr< jacobian_2D_2D<double> > aJacobian = boost::shared_ptr< jacobian_2D_2D<double> >()) :
		       reacting_kte_2D(aName),
		       mCoord(aCoord),
		       mBase(aBase),
		       mEnd(aEnd),
                       mJacobian(aJacobian) { };

    /**
     * Default destructor.
     */
    virtual ~free_joint_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void applyReactionForce(vect<double,2> aForce, double aTorque);

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int Version) const {
      reacting_kte_2D::save(A,reacting_kte_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCoord)
        & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
	& RK_SERIAL_SAVE_WITH_NAME(mJacobian);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int Version) {
      reacting_kte_2D::load(A,reacting_kte_2D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCoord)
        & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
	& RK_SERIAL_LOAD_WITH_NAME(mJacobian);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(free_joint_2D,0xC2100041,1,"free_joint_2D",reacting_kte_2D)

};

/**
 * This class implements a free joint in 3D space. A 3D frame is used to represent the
 * joint's motion between a base coordinate frame to an end coordinate frame.
 */
class free_joint_3D : public reacting_kte_3D {
  protected:
    boost::shared_ptr< frame_3D<double> > mCoord; ///< The coordinate frame representing the joint's motion.
    boost::shared_ptr< frame_3D<double> > mBase; ///< The coordinate frame at the base of the joint.
    boost::shared_ptr< frame_3D<double> > mEnd; ///< The coordinate frame just after the joint transformations are applied.

    boost::shared_ptr< jacobian_3D_3D<double> > mJacobian; ///< The Jacobian frame produced by this joint.

  public:

    /**
     * Default constructor.
     */
    free_joint_3D(const std::string& aName = "") : reacting_kte_3D(aName), mCoord(), mBase(), mEnd(), mJacobian() { };

    /**
     * Parametrized constructor.
     * \param aName the name of this KTE model.
     * \param aCoord the coordinate frame associated with the motion of this joint.
     * \param aBase the coordinate frame at the base of the joint.
     * \param aEnd the coordinate frame just after the joint transformations are applied.
     * \param aJacobian a pointer to contain the Jacobian frame produced by this joint, default value will disable the Jacobian frame's calculation.
     */
    free_joint_3D(const std::string& aName,
		  boost::shared_ptr< frame_3D<double> > aCoord,
		  boost::shared_ptr< frame_3D<double> > aBase,
		  boost::shared_ptr< frame_3D<double> > aEnd,
                  boost::shared_ptr< jacobian_3D_3D<double> > aJacobian = boost::shared_ptr< jacobian_3D_3D<double> >()) :
		  reacting_kte_3D(aName),
		  mCoord(aCoord),
		  mBase(aBase),
		  mEnd(aEnd),
                  mJacobian(aJacobian){ };

    /**
     * Default destructor.
     */
    virtual ~free_joint_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, boost::shared_ptr<frame_storage> aStorage = boost::shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void applyReactionForce(vect<double,3> aForce, vect<double,3> aTorque);

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      reacting_kte_3D::save(A,reacting_kte_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mCoord)
        & RK_SERIAL_SAVE_WITH_NAME(mBase)
        & RK_SERIAL_SAVE_WITH_NAME(mEnd)
	& RK_SERIAL_SAVE_WITH_NAME(mJacobian);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      reacting_kte_3D::load(A,reacting_kte_3D::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mCoord)
        & RK_SERIAL_LOAD_WITH_NAME(mBase)
        & RK_SERIAL_LOAD_WITH_NAME(mEnd)
	& RK_SERIAL_LOAD_WITH_NAME(mJacobian);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(free_joint_3D,0xC2100042,1,"free_joint_3D",reacting_kte_3D)

};




};

};






#endif










