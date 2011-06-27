
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

#include "spring.hpp"

namespace ReaK {

namespace kte {





void spring_gen::doMotion(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->gen_coord_mapping[mAnchor1]))
      aStorage->gen_coord_mapping[mAnchor1] = boost::shared_ptr< ReaK::gen_coord<double> >(new ReaK::gen_coord<double>((*mAnchor1)),ReaK::scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->gen_coord_mapping[mAnchor2]))
      aStorage->gen_coord_mapping[mAnchor2] = boost::shared_ptr< ReaK::gen_coord<double> >(new ReaK::gen_coord<double>((*mAnchor2)),ReaK::scoped_deleter());
    else
      (*(aStorage->gen_coord_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void spring_gen::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;
  using std::fabs;
  
  if(mAnchor1->q > mAnchor2->q) {
    double force_mag = (mAnchor1->q - mAnchor2->q - mRestLength) * mStiffness; //Positive if tension, negative if compression.
    if((mSaturation > 0) && (fabs(force_mag) > mSaturation)) {
      if(force_mag > 0) {
        mAnchor1->f -= mSaturation;
        mAnchor2->f += mSaturation;
      } else {
        mAnchor1->f += mSaturation;
        mAnchor2->f -= mSaturation;
      };
    } else {
      mAnchor1->f -= force_mag;
      mAnchor2->f += force_mag;
    };
  } else {
    double force_mag = (mAnchor2->q - mAnchor1->q - mRestLength) * mStiffness; //Positive if tension, negative if compression.
    if((mSaturation > 0) && (fabs(force_mag) > mSaturation)) {
      if(force_mag > 0) {
        mAnchor1->f += mSaturation;
        mAnchor2->f -= mSaturation;
      } else {
        mAnchor1->f -= mSaturation;
        mAnchor2->f += mSaturation;
      };
    } else {
      mAnchor1->f += force_mag;
      mAnchor2->f -= force_mag;
    };
  };
};


void spring_gen::clearForce() {
  if(mAnchor1) {
    mAnchor1->f = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->f = 0.0;
  };
};





void spring_2D::doMotion(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_2D_mapping[mAnchor1]))
      aStorage->frame_2D_mapping[mAnchor1] = boost::shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_2D_mapping[mAnchor2]))
      aStorage->frame_2D_mapping[mAnchor2] = boost::shared_ptr< frame_2D<double> >(new frame_2D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_2D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void spring_2D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;

  using std::fabs;

  vect<double,2> diff = mAnchor1->Position - mAnchor2->Position;
  double diff_mag = norm(diff);
  if(diff_mag > 1E-7) {
    double force_mag = (diff_mag - mRestLength) * mStiffness;

    if((mSaturation > 0) && (fabs(force_mag) > mSaturation)) {
      diff *= mSaturation / diff_mag; //Positive if tension, negative if compression.
      if(force_mag > 0) {
        mAnchor1->Force -= diff * mAnchor1->Rotation;
        mAnchor2->Force += diff * mAnchor2->Rotation;
      } else {
        mAnchor1->Force += diff * mAnchor1->Rotation;
        mAnchor2->Force -= diff * mAnchor2->Rotation;
      };
    } else {
      diff *= force_mag / diff_mag; //Positive if tension, negative if compression.
      mAnchor1->Force -= diff * mAnchor1->Rotation;
      mAnchor2->Force += diff * mAnchor2->Rotation;
    };
  };

};


void spring_2D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,2>();
    mAnchor1->Torque = 0.0;
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,2>();
    mAnchor2->Torque = 0.0;
  };
};






void spring_3D::doMotion(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) {
  if((!mAnchor1) || (!mAnchor2))
    return;

  if((aFlag == store_kinematics) && (aStorage)) {
    if(!(aStorage->frame_3D_mapping[mAnchor1]))
      aStorage->frame_3D_mapping[mAnchor1] = boost::shared_ptr< ReaK::frame_3D<double> >(new frame_3D<double>((*mAnchor1)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor1])) = (*mAnchor1);
    if(!(aStorage->frame_3D_mapping[mAnchor2]))
      aStorage->frame_3D_mapping[mAnchor2] = boost::shared_ptr< ReaK::frame_3D<double> >(new frame_3D<double>((*mAnchor2)),scoped_deleter());
    else
      (*(aStorage->frame_3D_mapping[mAnchor2])) = (*mAnchor2);
  };
};

void spring_3D::doForce(kte_pass_flag aFlag, boost::shared_ptr<frame_storage> aStorage) { RK_UNUSED(aFlag); RK_UNUSED(aStorage);
  if((!mAnchor1) || (!mAnchor2))
    return;
  
  using std::fabs;

  vect<double,3> diff = mAnchor1->Position - mAnchor2->Position;
  double diff_mag = norm(diff);
  if(diff_mag > 1E-7) {

    double force_mag = (diff_mag - mRestLength) * mStiffness;

    if((mSaturation > 0) && (fabs(force_mag) > mSaturation)) {
      diff *= mSaturation / diff_mag; //Positive if tension, negative if compression.
      if(force_mag > 0) {
        mAnchor1->Force -= invert(mAnchor1->Quat) * diff;
        mAnchor2->Force += invert(mAnchor2->Quat) * diff;
      } else {
        mAnchor1->Force += invert(mAnchor1->Quat) * diff;
        mAnchor2->Force -= invert(mAnchor2->Quat) * diff;
      };
    } else {
      diff *= force_mag / diff_mag; //Positive if tension, negative if compression.
      mAnchor1->Force -= invert(mAnchor1->Quat) * diff;
      mAnchor2->Force += invert(mAnchor2->Quat) * diff;
    };

  };

};


void spring_3D::clearForce() {
  if(mAnchor1) {
    mAnchor1->Force = vect<double,3>();
    mAnchor1->Torque = vect<double,3>();
  };
  if(mAnchor2) {
    mAnchor2->Force = vect<double,3>();
    mAnchor2->Torque = vect<double,3>();
  };
};



};


};







