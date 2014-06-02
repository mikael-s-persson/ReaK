/**
 * \file torsion_damper.hpp
 *
 * This library declares the KTE models for linear torsion dampers.
 * Here damper classes are available for 2D and 3D point-to-point angular difference. The model of the
 * damper is a basic linear, constant
 * damping coefficient that multiplies the angular velocity to obtain the torque that will oppose the relative
 * angular velocity between the anchors.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date May 2010
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

#ifndef REAK_TORSION_DAMPER_HPP
#define REAK_TORSION_DAMPER_HPP

#include "kte_map.hpp"

#include <ReaK/core/kinetostatics/kinetostatics.hpp>

namespace ReaK {

namespace kte {


/** This class defines a torsion damper acting between two 2D frames. */
class torsion_damper_2D : public kte_map {
  private:
    shared_ptr< frame_2D<double> > mAnchor1; ///< Holds the first 2D frame.
    shared_ptr< frame_2D<double> > mAnchor2; ///< Holds the second 2D frame.
    double mDamping; ///< The damping coefficient (in Nms/rad).

  public:
    
    /**
     * Sets the first anchor frame of the torsion damper.
     * \param aPtr The new first anchor frame of the torsion damper.
     */
    void setAnchor1(const shared_ptr< frame_2D<double> >& aPtr) { mAnchor1 = aPtr; };
    /**
     * Returns the first anchor frame of the torsion damper.
     * \return The first anchor frame of the torsion damper.
     */
    shared_ptr< frame_2D<double> > Anchor1() const { return mAnchor1; };
    
    /**
     * Sets the second anchor frame of the torsion damper.
     * \param aPtr The new second anchor frame of the torsion damper.
     */
    void setAnchor2(const shared_ptr< frame_2D<double> >& aPtr) { mAnchor2 = aPtr; };
    /**
     * Returns the second anchor frame of the torsion damper.
     * \return The second anchor frame of the torsion damper.
     */
    shared_ptr< frame_2D<double> > Anchor2() const { return mAnchor2; };
    
    /**
     * Sets the damping factor of the torsion damper.
     * \param aValue The new damping factor of the torsion damper.
     */
    void setDamping(double aValue) { mDamping = aValue; };
    /**
     * Returns the damping factor of the torsion damper.
     * \return The damping factor of the torsion damper.
     */
    double Damping() const { return mDamping; };
    
    /**
     * Default constructor.
     */
    torsion_damper_2D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mDamping(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the damper.
     * \param aAnchor2 second attach point of the damper.
     * \param aDamping damping coefficient (in Nms/rad).
     */
    torsion_damper_2D(const std::string& aName,
                      const shared_ptr< frame_2D<double> >& aAnchor1,
                      const shared_ptr< frame_2D<double> >& aAnchor2,
                      double aDamping) :
                      kte_map(aName),
                      mAnchor1(aAnchor1),
                      mAnchor2(aAnchor2),
                      mDamping(aDamping) { };

    /**
     * Default destructor.
     */
    virtual ~torsion_damper_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
        & RK_SERIAL_SAVE_WITH_NAME(mDamping);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
        & RK_SERIAL_LOAD_WITH_NAME(mDamping);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(torsion_damper_2D,0xC210002E,1,"torsion_damper_2D",kte_map)

};

/** This class defines a torsion damper acting between two 3D frames. */
class torsion_damper_3D : public kte_map {
  private:
    shared_ptr< frame_3D<double> > mAnchor1; ///< Holds the first 3D frame.
    shared_ptr< frame_3D<double> > mAnchor2; ///< Holds the second 3D frame.
    double mDamping; ///< The damping coefficient (in Nms/rad).

  public:
    
    /**
     * Sets the first anchor frame of the torsion damper.
     * \param aPtr The new first anchor frame of the torsion damper.
     */
    void setAnchor1(const shared_ptr< frame_3D<double> >& aPtr) { mAnchor1 = aPtr; };
    /**
     * Returns the first anchor frame of the torsion damper.
     * \return The first anchor frame of the torsion damper.
     */
    shared_ptr< frame_3D<double> > Anchor1() const { return mAnchor1; };
    
    /**
     * Sets the second anchor frame of the torsion damper.
     * \param aPtr The new second anchor frame of the torsion damper.
     */
    void setAnchor2(const shared_ptr< frame_3D<double> >& aPtr) { mAnchor2 = aPtr; };
    /**
     * Returns the second anchor frame of the torsion damper.
     * \return The second anchor frame of the torsion damper.
     */
    shared_ptr< frame_3D<double> > Anchor2() const { return mAnchor2; };
    
    /**
     * Sets the damping factor of the torsion damper.
     * \param aValue The new damping factor of the torsion damper.
     */
    void setDamping(double aValue) { mDamping = aValue; };
    /**
     * Returns the damping factor of the torsion damper.
     * \return The damping factor of the torsion damper.
     */
    double Damping() const { return mDamping; };
    
    /**
     * Default constructor.
     */
    torsion_damper_3D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mDamping(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the damper.
     * \param aAnchor2 second attach point of the damper.
     * \param aDamping damping coefficient (in Nms/rad).
     */
    torsion_damper_3D(const std::string& aName,
                      const shared_ptr< frame_3D<double> >& aAnchor1,
                      const shared_ptr< frame_3D<double> >& aAnchor2,
                      double aDamping) :
                      kte_map(aName),
                      mAnchor1(aAnchor1),
                      mAnchor2(aAnchor2),
                      mDamping(aDamping) { };

    /**
     * Default destructor.
     */
    virtual ~torsion_damper_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
        & RK_SERIAL_SAVE_WITH_NAME(mDamping);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
        & RK_SERIAL_LOAD_WITH_NAME(mDamping);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(torsion_damper_3D,0xC210002F,1,"torsion_damper_3D",kte_map)
};


};

};

#endif





