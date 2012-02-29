/**
 * \file spring.hpp
 *
 * This library declares the KTE models for linear springs.
 * Here spring classes are available for 2D and 3D point-to-point springs as well as a
 * spring between two generalized coordinates. The model of the spring is a basic linear, constant
 * stiffness with optional staturation value based on Hooke's Law to obtain the force that will oppose the relative
 * distance between the anchors, offsetted by the rest-length of the spring.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
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

#ifndef REAK_SPRING_HPP
#define REAK_SPRING_HPP

#include "kte_map.hpp"

#include "kinetostatics/kinetostatics.hpp"

namespace ReaK {

namespace kte {


/**
 * This class implements a spring model on a generalized coordinate system.The model of the spring is a basic linear, constant
 * stiffness with optional staturation value based on Hooke's Law to obtain the force that will oppose the relative
 * distance between the anchors, offsetted by the rest-length of the spring.
 */
class spring_gen : public kte_map {
  private:
    shared_ptr< gen_coord<double> > mAnchor1; ///< Holds the first generalized coordinate.
    shared_ptr< gen_coord<double> > mAnchor2; ///< Holds the second generalized coordinate.
    double mRestLength; ///< Holds the rest-length of the spring.
    double mStiffness; ///< Holds the stiffness of the spring.
    double mSaturation; ///< Holds the saturation force, or maximum force the spring can exert, if 0 there is no saturation.

  public:

    /** Get read-write access to mRestLength. */
    double& RestLength() { return mRestLength; };
    /** Get read-only access to mRestLength. */
    double RestLength() const { return mRestLength; };

    /** Get read-write access to mStiffness. */
    double& Stiffness() { return mStiffness; };
    /** Get read-only access to mStiffness. */
    double Stiffness() const { return mStiffness; };

    /** Get read-write access to mSaturation. */
    double& Saturation() { return mSaturation; };
    /** Get read-only access to mSaturation. */
    double Saturation() const { return mSaturation; };

    /**
     * Default constructor.
     */
    spring_gen(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mRestLength(0.0), mStiffness(0.0), mSaturation(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the spring.
     * \param aAnchor2 second attach point of the spring.
     * \param aRestLength rest-length of the spring (m or rad).
     * \param aStiffness stiffness coefficient (in N/m for a linear generalized coord. or Nm/rad for an angular generalized coord.).
     * \param aSaturation saturation force of the spring, default 0 will disable saturation.
     */
    spring_gen(const std::string& aName,
	       const shared_ptr< gen_coord<double> >& aAnchor1,
	       const shared_ptr< gen_coord<double> >& aAnchor2,
	       double aRestLength,
	       double aStiffness,
               double aSaturation = 0.0) :
	       kte_map(aName),
	       mAnchor1(aAnchor1),
	       mAnchor2(aAnchor2),
	       mRestLength(aRestLength),
	       mStiffness(aStiffness),
	       mSaturation(aSaturation) { };

    /**
     * Default destructor.
     */
    virtual ~spring_gen() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();


    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mRestLength)
	& RK_SERIAL_SAVE_WITH_NAME(mStiffness)
	& RK_SERIAL_SAVE_WITH_NAME(mSaturation);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int Version) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mRestLength)
	& RK_SERIAL_LOAD_WITH_NAME(mStiffness)
	& RK_SERIAL_LOAD_WITH_NAME(mSaturation);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(spring_gen,0xC210000D,1,"spring_gen",kte_map)

};

/**
 * This class implements a spring model in 2D space. The model of the spring is a basic linear, constant
 * stiffness with optional staturation value based on Hooke's Law to obtain the force that will oppose the relative
 * distance between the anchors, offsetted by the rest-length of the spring.
 */
class spring_2D : public kte_map {
  private:
    shared_ptr< frame_2D<double> > mAnchor1; ///< Holds the first 2D frame.
    shared_ptr< frame_2D<double> > mAnchor2; ///< Holds the second 2D frame.
    double mRestLength; ///< Holds the rest-length of the spring.
    double mStiffness; ///< Holds the stiffness of the spring.
    double mSaturation; ///< Holds the saturation force, or maximum force the spring can exert, if 0 there is no saturation.

  public:

    /** Get read-write access to mRestLength. */
    double& RestLength() { return mRestLength; };
    /** Get read-only access to mRestLength. */
    double RestLength() const { return mRestLength; };

    /** Get read-write access to mStiffness. */
    double& Stiffness() { return mStiffness; };
    /** Get read-only access to mStiffness. */
    double Stiffness() const { return mStiffness; };

    /** Get read-write access to mSaturation. */
    double& Saturation() { return mSaturation; };
    /** Get read-only access to mSaturation. */
    double Saturation() const { return mSaturation; };

    /**
     * Default constructor.
     */
    spring_2D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mRestLength(0.0), mStiffness(0.0), mSaturation(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the spring.
     * \param aAnchor2 second attach point of the spring.
     * \param aRestLength rest-length of the spring (m).
     * \param aStiffness stiffness coefficient (in N/m).
     * \param aSaturation saturation force of the spring, default 0 will disable saturation.
     */
    spring_2D(const std::string& aName,
	      const shared_ptr< frame_2D<double> >& aAnchor1,
	      const shared_ptr< frame_2D<double> >& aAnchor2,
	      double aRestLength,
	      double aStiffness,
              double aSaturation = 0.0) :
	      kte_map(aName),
	      mAnchor1(aAnchor1),
	      mAnchor2(aAnchor2),
	      mRestLength(aRestLength),
	      mStiffness(aStiffness),
	      mSaturation(aSaturation) { };

    /**
     * Default destructor.
     */
    virtual ~spring_2D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mRestLength)
	& RK_SERIAL_SAVE_WITH_NAME(mStiffness)
	& RK_SERIAL_SAVE_WITH_NAME(mSaturation);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int Version) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mRestLength)
	& RK_SERIAL_LOAD_WITH_NAME(mStiffness)
	& RK_SERIAL_LOAD_WITH_NAME(mSaturation);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(spring_2D,0xC210000E,1,"spring_2D",kte_map)

};

/**
 * This class implements a spring model in 3D space. The model of the spring is a basic linear, constant
 * stiffness with optional staturation value based on Hooke's Law to obtain the force that will oppose the relative
 * distance between the anchors, offsetted by the rest-length of the spring.
 */
class spring_3D : public kte_map {
  private:
    shared_ptr< frame_3D<double> > mAnchor1; ///< Holds the first 3D frame.
    shared_ptr< frame_3D<double> > mAnchor2; ///< Holds the second 3D frame.
    double mRestLength; ///< Holds the rest-length of the spring.
    double mStiffness; ///< Holds the stiffness of the spring.
    double mSaturation; ///< Holds the saturation force, or maximum force the spring can exert, if 0 there is no saturation.

  public:

    /** Get read-write access to mRestLength. */
    double& RestLength() { return mRestLength; };
    /** Get read-only access to mRestLength. */
    double RestLength() const { return mRestLength; };

    /** Get read-write access to mStiffness. */
    double& Stiffness() { return mStiffness; };
    /** Get read-only access to mStiffness. */
    double Stiffness() const { return mStiffness; };

    /** Get read-write access to mSaturation. */
    double& Saturation() { return mSaturation; };
    /** Get read-only access to mSaturation. */
    double Saturation() const { return mSaturation; };

    /**
     * Default constructor.
     */
    spring_3D(const std::string& aName = "") : kte_map(aName), mAnchor1(), mAnchor2(), mRestLength(0.0), mStiffness(0.0), mSaturation(0.0) { };

    /**
     * Parametrized constructor.
     * \param aName name of the KTE model.
     * \param aAnchor1 first attach point of the spring.
     * \param aAnchor2 second attach point of the spring.
     * \param aRestLength rest-length of the spring (m).
     * \param aStiffness stiffness coefficient (in N/m).
     * \param aSaturation saturation force of the spring, default 0 will disable saturation.
     */
    spring_3D(const std::string& aName,
	      const shared_ptr< frame_3D<double> >& aAnchor1,
	      const shared_ptr< frame_3D<double> >& aAnchor2,
	      double aRestLength,
	      double aStiffness,
              double aSaturation = 0.0) :
	      kte_map(aName),
	      mAnchor1(aAnchor1),
	      mAnchor2(aAnchor2),
	      mRestLength(aRestLength),
	      mStiffness(aStiffness),
	      mSaturation(aSaturation) { };

    /**
     * Default destructor.
     */
    virtual ~spring_3D() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_ptr<frame_storage>& aStorage = shared_ptr<frame_storage>());

    virtual void clearForce();

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      kte_map::save(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mAnchor1)
        & RK_SERIAL_SAVE_WITH_NAME(mAnchor2)
	& RK_SERIAL_SAVE_WITH_NAME(mRestLength)
	& RK_SERIAL_SAVE_WITH_NAME(mStiffness)
	& RK_SERIAL_SAVE_WITH_NAME(mSaturation);
    };

    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int Version) {
      kte_map::load(A,kte_map::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mAnchor1)
        & RK_SERIAL_LOAD_WITH_NAME(mAnchor2)
	& RK_SERIAL_LOAD_WITH_NAME(mRestLength)
	& RK_SERIAL_LOAD_WITH_NAME(mStiffness)
	& RK_SERIAL_LOAD_WITH_NAME(mSaturation);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(spring_3D,0xC210000F,1,"spring_3D",kte_map)
};


};

};

#endif










