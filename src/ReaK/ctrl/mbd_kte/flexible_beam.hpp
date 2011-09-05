/**
 * \file flexible_beam.hpp
 *
 * This library declares KTE models for a flexible beam in 2D and 3D. The flexible beam is
 * implemented as a linear iso-tropic stiffness for both elongation and torsion. This
 * constitutes a simplified Euler-Bernoulli beam, akin to beam elements in finite element methods.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2010
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

#ifndef REAK_FLEXIBLE_BEAM_HPP
#define REAK_FLEXIBLE_BEAM_HPP

#include "kte_map.hpp"
#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements the 2D model of a flexible beam. The beam is made of three 2D frames, two inputs and
 * one output (kinematically speaking). The differencial displacement between the two anchors is used to compute
 * restitution forces and torques at the anchors. The average or common-mode or net motion of the two anchors is
 * used to compute the kinematics of the center of the beam, the object-frame. Furthermore, the net forces applied
 * on the object-frame is rigidly transmitted and evenly split to the anchors, in addition to the restitution forces.
 */
class flexible_beam_2D : public kte_map {
  private:
    shared_pointer< frame_2D<double> >::type mAnchor1; ///< Holds the first end of the beam.
    shared_pointer< frame_2D<double> >::type mAnchor2; ///< Holds the second end of the beam.
    shared_pointer< frame_2D<double> >::type mObjectFrame; ///< Holds the center frame of the beam (i.e. its bulk).
    double mRestLength; ///< The undeformed length of the beam.
    double mStiffness; ///< The linear stiffness of the beam, stress-strain relation, iso-tropically.
    double mTorsionStiffness; ///< The angular or torsion stiffness of the beam, iso-tropically.

  public:

    /** Get read-write access to the mRestLength. */
    double& RestLength() { return mRestLength; };
    /** Get read-only access to the mRestLength. */
    double RestLength() const { return mRestLength; };

    /** Get read-write access to the mStiffness. */
    double& Stiffness() { return mStiffness; };
    /** Get read-only access to the mStiffness. */
    double Stiffness() const { return mStiffness; };

    /** Get read-write access to the mTorsionStiffness. */
    double& TorsionStiffness() { return mTorsionStiffness; };
    /** Get read-only access to the mTorsionStiffness. */
    double TorsionStiffness() const { return mTorsionStiffness; };


    /**
     * Default constructor.
     */
    flexible_beam_2D(const std::string& aName = "") : kte_map(aName),
                                                      mAnchor1(),
                                                      mAnchor2(),
                                                      mObjectFrame(),
                                                      mRestLength(0.0),
                                                      mStiffness(0.0),
                                                      mTorsionStiffness(0.0){ };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aAnchor1 the first end of the beam (kinematic input).
     * \param aAnchor2 the second end of the beam (kinematic input).
     * \param aObjectFrame the center frame of the beam (i.e. its bulk).
     * \param aRestLength the undeformed length of the beam.
     * \param aStiffness the linear stiffness of the beam, stress-strain relation, iso-tropically.
     * \param aTorsionStiffness the angular or torsion stiffness of the beam, iso-tropically.
     */
    flexible_beam_2D(const std::string& aName,
                     const shared_pointer< frame_2D<double> >::type& aAnchor1,
                     const shared_pointer< frame_2D<double> >::type& aAnchor2,
                     const shared_pointer< frame_2D<double> >::type& aObjectFrame,
                     double aRestLength,
                     double aStiffness,
                     double aTorsionStiffness) :
                     kte_map(aName),
                     mAnchor1(aAnchor1),
                     mAnchor2(aAnchor2),
                     mObjectFrame(aObjectFrame),
                     mRestLength(aRestLength),
                     mStiffness(aStiffness),
                     mTorsionStiffness(aTorsionStiffness){ };

    /**
     * Default destructor.
     */
    virtual ~flexible_beam_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
        & RK_SERIAL_SAVE_WITH_NAME(mObjectFrame)
        & RK_SERIAL_SAVE_WITH_NAME(mRestLength)
        & RK_SERIAL_SAVE_WITH_NAME(mStiffness)
        & RK_SERIAL_SAVE_WITH_NAME(mTorsionStiffness);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
        & RK_SERIAL_LOAD_WITH_NAME(mObjectFrame)
        & RK_SERIAL_LOAD_WITH_NAME(mRestLength)
        & RK_SERIAL_LOAD_WITH_NAME(mStiffness)
        & RK_SERIAL_LOAD_WITH_NAME(mTorsionStiffness);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(flexible_beam_2D,0xC210001D,1,"flexible_beam_2D",kte_map)

};

/**
 * This class implements the 3D model of a flexible beam. The beam is made of three 3D frames, two inputs and
 * one output (kinematically speaking). The differencial displacement between the two anchors is used to compute
 * restitution forces and torques at the anchors. The average or common-mode or net motion of the two anchors is
 * used to compute the kinematics of the center of the beam, the object-frame. Furthermore, the net forces applied
 * on the object-frame is rigidly transmitted and evenly split to the anchors, in addition to the restitution forces.
 */
class flexible_beam_3D : public kte_map {
  private:
    shared_pointer< frame_3D<double> >::type mAnchor1; ///< Holds the first end of the beam.
    shared_pointer< frame_3D<double> >::type mAnchor2; ///< Holds the second end of the beam.
    shared_pointer< frame_3D<double> >::type mObjectFrame; ///< Holds the center frame of the beam (i.e. its bulk).
    double mRestLength; ///< The undeformed length of the beam.
    double mStiffness; ///< The linear stiffness of the beam, stress-strain relation, iso-tropically.
    double mTorsionStiffness; ///< The angular or torsion stiffness of the beam, iso-tropically.

  public:

    /** Get read-write access to the mRestLength. */
    double& RestLength() { return mRestLength; };
    /** Get read-only access to the mRestLength. */
    double RestLength() const { return mRestLength; };

    /** Get read-write access to the mStiffness. */
    double& Stiffness() { return mStiffness; };
    /** Get read-only access to the mStiffness. */
    double Stiffness() const { return mStiffness; };

    /** Get read-write access to the mTorsionStiffness. */
    double& TorsionStiffness() { return mTorsionStiffness; };
    /** Get read-only access to the mTorsionStiffness. */
    double TorsionStiffness() const { return mTorsionStiffness; };

    /**
     * Default constructor.
     */
    flexible_beam_3D(const std::string& aName = "") : kte_map(aName),
                                                      mAnchor1(),
                                                      mAnchor2(),
                                                      mObjectFrame(),
                                                      mRestLength(0.0),
                                                      mStiffness(0.0),
                                                      mTorsionStiffness(0.0){ };

    /**
     * Parametrized constructor.
     * \param aName the name of the KTE model.
     * \param aAnchor1 the first end of the beam (kinematic input).
     * \param aAnchor2 the second end of the beam (kinematic input).
     * \param aObjectFrame the center frame of the beam (i.e. its bulk).
     * \param aRestLength the undeformed length of the beam.
     * \param aStiffness the linear stiffness of the beam, stress-strain relation, iso-tropically.
     * \param aTorsionStiffness the angular or torsion stiffness of the beam, iso-tropically.
     */
    flexible_beam_3D(const std::string& aName,
                     shared_pointer< frame_3D<double> >::type& aAnchor1,
                     shared_pointer< frame_3D<double> >::type& aAnchor2,
                     shared_pointer< frame_3D<double> >::type& aObjectFrame,
                     double aRestLength,
                     double aStiffness,
                     double aTorsionStiffness) :
                     kte_map(aName),
                     mAnchor1(aAnchor1),
                     mAnchor2(aAnchor2),
                     mObjectFrame(aObjectFrame),
                     mRestLength(aRestLength),
                     mStiffness(aStiffness),
                     mTorsionStiffness(aTorsionStiffness){ };

    /**
     * Default destructor.
     */
    virtual ~flexible_beam_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
        & RK_SERIAL_SAVE_WITH_NAME(mObjectFrame)
        & RK_SERIAL_SAVE_WITH_NAME(mRestLength)
        & RK_SERIAL_SAVE_WITH_NAME(mStiffness)
        & RK_SERIAL_SAVE_WITH_NAME(mTorsionStiffness);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
        & RK_SERIAL_LOAD_WITH_NAME(mObjectFrame)
        & RK_SERIAL_LOAD_WITH_NAME(mRestLength)
        & RK_SERIAL_LOAD_WITH_NAME(mStiffness)
        & RK_SERIAL_LOAD_WITH_NAME(mTorsionStiffness);

    };

    RK_RTTI_MAKE_CONCRETE_1BASE(flexible_beam_3D,0xC210001D,1,"flexible_beam_3D",kte_map)
};


};

};



#endif








