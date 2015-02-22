/**
 * \file prox_fundamentals_3D.hpp
 *
 * This library declares fundamental proximity query methods for 3D.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date April 2012
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

#ifndef REAK_PROX_FUNDAMENTALS_3D_HPP
#define REAK_PROX_FUNDAMENTALS_3D_HPP

#include "proximity_finder_3D.hpp"

#include <ReaK/geometry/shapes/box.hpp>
#include <ReaK/geometry/shapes/cylinder.hpp>
#include <ReaK/geometry/shapes/capped_cylinder.hpp>

#include <ReaK/math/lin_alg/mat_alg.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D findProximityBoxToPoint(const box& aBox, const pose_3D<double>& aBoxGblPose, const vect<double,3>& aPoint);


proximity_record_3D findProximityBoxToLine(const box& aBox, const pose_3D<double>& aBoxGblPose, const vect<double,3>& aCenter, const vect<double,3>& aTangent, double aHalfLength);



struct slack_minimize_func {
  double operator()(const vect_n<double>& aX) {
    return aX[0];
  };
};

struct slack_minimize_grad {
  vect_n<double> operator()(const vect_n<double>& aX) {
    return vect_n<double>(1.0, 0.0, 0.0, 0.0);
    //return vect_n<double>(aX[0], 0.0, 0.0, 0.0);
  };
};

struct slack_minimize_hess {
  template <typename Matrix>
  void operator()(Matrix& H, const vect_n<double>&, double, const vect_n<double>&) const {
    H = mat<double,mat_structure::nil>(4,4);
    //H(0,0) = 1.0;
  };
};


struct cylinder_slacking_func {
  const cylinder* mCylinder;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 3;
  
  cylinder_slacking_func(const cylinder& aCylinder, const pose_3D<double>& aGblPose) : 
                         mCylinder(&aCylinder), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    using std::sqrt;
    vect_n<double> result(3);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mGblPose->transformFromGlobal(pt);
    result[0] = pt_rel[2] + 0.5 * aX[0] * mCylinder->getLength(); // lower-bound.
    result[1] = 0.5 * aX[0] * mCylinder->getLength() - pt_rel[2]; // upper-bound.
    result[2] = aX[0] * mCylinder->getRadius() - sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1]);
    //result[2] = aX[0] * aX[0] * mCylinder->getRadius() * mCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1];
    return result;
  };
};

struct cylinder_slacking_jac {
  const cylinder* mCylinder;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 3;
  
  cylinder_slacking_jac(const cylinder& aCylinder, const pose_3D<double>& aGblPose) : 
                        mCylinder(&aCylinder), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    using std::sqrt;
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mGblPose->transformFromGlobal(pt);
    
    vect<double,3> x_rel = mGblPose->rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
    //aJac.resize(3,4);
    
    //result[0] = pt_rel[2] + 0.5 * aX[0] * mCylinder->getLength(); // lower-bound.
    aJac(0,0) = 0.5 * mCylinder->getLength();
    aJac(0,1) = x_rel[2];
    aJac(0,2) = y_rel[2];
    aJac(0,3) = z_rel[2];
    
    //result[1] = 0.5 * aX[0] * mCylinder->getLength() - pt_rel[2]; // upper-bound.
    aJac(1,0) = 0.5 * mCylinder->getLength();
    aJac(1,1) = -x_rel[2];
    aJac(1,2) = -y_rel[2];
    aJac(1,3) = -z_rel[2];
    
    //result[2] = aX[0] * aX[0] * mCylinder->getRadius() * mCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1];
    aJac(2,0) = mCylinder->getRadius();
    double tmp = sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1]);
    tmp = (tmp < 1e-6 ? 1e-6 : tmp);
    aJac(2,1) = -(pt_rel[0] * x_rel[0] + pt_rel[1] * x_rel[1]) / tmp;
    aJac(2,2) = -(pt_rel[0] * y_rel[0] + pt_rel[1] * y_rel[1]) / tmp;
    aJac(2,3) = -(pt_rel[0] * z_rel[0] + pt_rel[1] * z_rel[1]) / tmp;
    //aJac(2,0) = 2.0 * aX[0] * mCylinder->getRadius() * mCylinder->getRadius();
    //aJac(2,1) = -2.0 * (pt_rel[0] * x_rel[0] + pt_rel[1] * x_rel[1]);
    //aJac(2,2) = -2.0 * (pt_rel[0] * y_rel[0] + pt_rel[1] * y_rel[1]);
    //aJac(2,3) = -2.0 * (pt_rel[0] * z_rel[0] + pt_rel[1] * z_rel[1]);
    
  };
};

struct ccylinder_slacking_func {
  const capped_cylinder* mCCylinder;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 1;
  
  ccylinder_slacking_func(const capped_cylinder& aCCylinder, const pose_3D<double>& aGblPose) :
                          mCCylinder(&aCCylinder), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    using std::sqrt;
    vect_n<double> result(1);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mGblPose->transformFromGlobal(pt);
    
    if(pt_rel[2] > 0.5 * aX[0] * mCCylinder->getLength()) {
      double tmp = pt_rel[2] - 0.5 * aX[0] * mCCylinder->getLength();
      result[0] = aX[0] * mCCylinder->getRadius() - sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1] + tmp * tmp);
    } else if(pt_rel[2] < -0.5 * aX[0] * mCCylinder->getLength()) {
      double tmp = pt_rel[2] + 0.5 * aX[0] * mCCylinder->getLength();
      result[0] = aX[0] * mCCylinder->getRadius() - sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1] + tmp * tmp);
    } else {
      result[0] = aX[0] * mCCylinder->getRadius() - sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1]); 
    };
    
    return result;
  };
};

struct ccylinder_slacking_jac {
  const capped_cylinder* mCCylinder;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 1;
  
  ccylinder_slacking_jac(const capped_cylinder& aCCylinder, const pose_3D<double>& aGblPose) :
                         mCCylinder(&aCCylinder), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    using std::sqrt;
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mGblPose->transformFromGlobal(pt);
    
    vect<double,3> x_rel = mGblPose->rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
    //aJac.resize(1,4);
    aJac(0,0) = mCCylinder->getRadius();
    
    double tmp = 0.0;
    double dist = 0.0;
    if(pt_rel[2] > 0.5 * aX[0] * mCCylinder->getLength()) {
      tmp = pt_rel[2] - 0.5 * aX[0] * mCCylinder->getLength();
      dist = sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1] + tmp * tmp);
      dist = (dist < 1e-6 ? 1e-6 : dist);
      aJac(0,0) -= 0.5 * tmp * mCCylinder->getLength() / dist;
    } else if(pt_rel[2] < -0.5 * aX[0] * mCCylinder->getLength()) {
      tmp = pt_rel[2] + 0.5 * aX[0] * mCCylinder->getLength();
      dist = sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1] + tmp * tmp);
      dist = (dist < 1e-6 ? 1e-6 : dist);
      aJac(0,0) += 0.5 * tmp * mCCylinder->getLength() / dist;
    } else {
      dist = sqrt(pt_rel[0] * pt_rel[0] + pt_rel[1] * pt_rel[1]);
      dist = (dist < 1e-6 ? 1e-6 : dist);
    };
    aJac(0,1) = -(pt_rel[0] * x_rel[0] + pt_rel[1] * x_rel[1] + tmp * x_rel[2]) / dist;
    aJac(0,2) = -(pt_rel[0] * y_rel[0] + pt_rel[1] * y_rel[1] + tmp * y_rel[2]) / dist;
    aJac(0,3) = -(pt_rel[0] * z_rel[0] + pt_rel[1] * z_rel[1] + tmp * z_rel[2]) / dist;
    
  };
};


struct box_slacking_func {
  const box* mBox;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 6;
  
  box_slacking_func(const box& aBox, const pose_3D<double>& aGblPose) : 
                    mBox(&aBox), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    vect_n<double> result(6);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mGblPose->transformFromGlobal(pt);
    result[0] = pt_rel[0] + 0.5 * aX[0] * mBox->getDimensions()[0]; // lower-bound.
    result[1] = pt_rel[1] + 0.5 * aX[0] * mBox->getDimensions()[1]; // lower-bound.
    result[2] = pt_rel[2] + 0.5 * aX[0] * mBox->getDimensions()[2]; // lower-bound.
    result[3] = 0.5 * aX[0] * mBox->getDimensions()[0] - pt_rel[0]; // upper-bound.
    result[4] = 0.5 * aX[0] * mBox->getDimensions()[1] - pt_rel[1]; // upper-bound.
    result[5] = 0.5 * aX[0] * mBox->getDimensions()[2] - pt_rel[2]; // upper-bound.
    return result;
  };
};

struct box_slacking_jac {
  const box* mBox;
  const pose_3D<double>* mGblPose;
  
  static const std::size_t size = 6;
  
  box_slacking_jac(const box& aBox, const pose_3D<double>& aGblPose) : 
                   mBox(&aBox), mGblPose(&aGblPose) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    vect<double,3> x_rel = mGblPose->rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mGblPose->rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
    //aJac.resize(6,4);
    aJac(0,0) = 0.5 * mBox->getDimensions()[0]; aJac(0,1) =  x_rel[0]; aJac(0,2) =  y_rel[0]; aJac(0,3) =  z_rel[0];
    aJac(1,0) = 0.5 * mBox->getDimensions()[1]; aJac(1,1) =  x_rel[1]; aJac(1,2) =  y_rel[1]; aJac(1,3) =  z_rel[1];
    aJac(2,0) = 0.5 * mBox->getDimensions()[2]; aJac(2,1) =  x_rel[2]; aJac(2,2) =  y_rel[2]; aJac(2,3) =  z_rel[2];
    aJac(3,0) = 0.5 * mBox->getDimensions()[0]; aJac(3,1) = -x_rel[0]; aJac(3,2) = -y_rel[0]; aJac(3,3) = -z_rel[0];
    aJac(4,0) = 0.5 * mBox->getDimensions()[1]; aJac(4,1) = -x_rel[1]; aJac(4,2) = -y_rel[1]; aJac(4,3) = -z_rel[1];
    aJac(5,0) = 0.5 * mBox->getDimensions()[2]; aJac(5,1) = -x_rel[2]; aJac(5,2) = -y_rel[2]; aJac(5,3) = -z_rel[2];
  };
};

template <typename BoundFunc1, typename BoundFunc2>
struct dual_slacking_func {
  BoundFunc1 mBound1;
  BoundFunc2 mBound2;
  
  dual_slacking_func(BoundFunc1 aBound1, BoundFunc2 aBound2) : mBound1(aBound1), mBound2(aBound2) { };
  
  vect_n<double> operator()(const vect_n<double>& aX) {
    vect_n<double> r1 = mBound1(aX);
    vect_n<double> r2 = mBound2(aX);
    vect_n<double> result(r1.size() + r2.size() + 1);
    result[0] = aX[0]; // slack-bound.
    std::copy(r2.begin(), r2.end(), std::copy(r1.begin(), r1.end(), result.begin() + 1));
    return result;
  };
  
};

template <typename BoundJac1, typename BoundJac2>
struct dual_slacking_jac {
  BoundJac1 mBound1;
  BoundJac2 mBound2;
  
  dual_slacking_jac(BoundJac1 aBound1, BoundJac2 aBound2) : mBound1(aBound1), mBound2(aBound2) { };
  
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    //aJac.resize(aH.size(), aX.size());
    aJac(0,0) = 1.0; aJac(0,1) = 0.0; aJac(0,2) = 0.0; aJac(0,3) = 0.0; // slack-bound.
    mat_sub_block<Matrix> aJac_s1(aJac, BoundJac1::size, 4, 1, 0);
    mat_sub_block<Matrix> aJac_s2(aJac, BoundJac2::size, 4, BoundJac1::size + 1, 0);
    vect_n<double> h1(BoundJac1::size);
    vect_n<double> h2(BoundJac2::size);
    std::copy(aH.begin(), aH.begin() + BoundJac1::size, h1.begin());
    std::copy(aH.begin() + BoundJac1::size, aH.begin() + BoundJac1::size + BoundJac2::size, h2.begin());
    
    mBound1(aJac_s1, aX, h1);
    mBound2(aJac_s2, aX, h2);
  };
  
};





struct cylinder_boundary_func {
  const cylinder* mCylinder;
  const pose_3D<double>* mGblPose;
  
  cylinder_boundary_func(const cylinder& aCylinder, const pose_3D<double>& aGblPose) : 
                        mCylinder(&aCylinder), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  vect<double,3> operator()(vect<double,3> v);
};

struct cylinder_boundary_jac {
  const cylinder* mCylinder;
  const pose_3D<double>* mGblPose;
  
  cylinder_boundary_jac(const cylinder& aCylinder, const pose_3D<double>& aGblPose) : 
                        mCylinder(&aCylinder), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  mat<double,mat_structure::square> operator()(vect<double,3> v);
};

struct ccylinder_boundary_func {
  const capped_cylinder* mCCylinder;
  const pose_3D<double>* mGblPose;
  
  ccylinder_boundary_func(const capped_cylinder& aCCylinder, const pose_3D<double>& aGblPose) :
                          mCCylinder(&aCCylinder), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  vect<double,3> operator()(vect<double,3> v);
};

struct ccylinder_boundary_jac {
  const capped_cylinder* mCCylinder;
  const pose_3D<double>* mGblPose;
  
  ccylinder_boundary_jac(const capped_cylinder& aCCylinder, const pose_3D<double>& aGblPose) :
                         mCCylinder(&aCCylinder), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  mat<double,mat_structure::square> operator()(vect<double,3> v);
};


struct box_boundary_func {
  const box* mBox;
  const pose_3D<double>* mGblPose;
  
  box_boundary_func(const box& aBox, const pose_3D<double>& aGblPose) : 
                    mBox(&aBox), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  vect<double,3> operator()(vect<double,3> v);
};

struct box_boundary_jac {
  const box* mBox;
  const pose_3D<double>* mGblPose;
  
  box_boundary_jac(const box& aBox, const pose_3D<double>& aGblPose) : 
                   mBox(&aBox), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  mat<double,mat_structure::square> operator()(vect<double,3> v);
};


struct support_func_base {
  virtual ~support_func_base() { };
  
  virtual vect<double,3> operator()(vect<double,3> v) const = 0;
};

struct cylinder_support_func : support_func_base {
  const cylinder* mCylinder;
  const pose_3D<double>* mGblPose;
  
  cylinder_support_func(const cylinder& aCylinder, const pose_3D<double>& aGblPose) : 
                        mCylinder(&aCylinder), mGblPose(&aGblPose) { };
  
  // aX is the query-direction (aX[0],aX[1],aX[2]).
  vect<double,3> operator()(vect<double,3> v) const;
};

struct ccylinder_support_func : support_func_base {
  const capped_cylinder* mCCylinder;
  const pose_3D<double>* mGblPose;
  
  ccylinder_support_func(const capped_cylinder& aCCylinder, const pose_3D<double>& aGblPose) :
                         mCCylinder(&aCCylinder), mGblPose(&aGblPose) { };
  
  // v is the query-direction.
  vect<double,3> operator()(vect<double,3> v) const;
};

struct box_support_func : support_func_base {
  const box* mBox;
  const pose_3D<double>* mGblPose;
  
  box_support_func(const box& aBox, const pose_3D<double>& aGblPose) : 
                   mBox(&aBox), mGblPose(&aGblPose) { };
  
  // v is the query-direction.
  vect<double,3> operator()(vect<double,3> v) const;
};



proximity_record_3D findProximityByGJKEPA(const support_func_base& aShape1, const support_func_base& aShape2);




};

};

#endif










