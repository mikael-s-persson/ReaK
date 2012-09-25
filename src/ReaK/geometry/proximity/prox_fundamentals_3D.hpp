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

#include "shapes/box.hpp"
#include "shapes/cylinder.hpp"
#include "shapes/capped_cylinder.hpp"

#include "lin_alg/mat_alg.hpp"

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


proximity_record_3D findProximityBoxToPoint(const shared_ptr< box >& aBox, const vect<double,3>& aPoint);


proximity_record_3D findProximityBoxToLine(const shared_ptr< box >& aBox, const vect<double,3>& aCenter, const vect<double,3>& aTangent, double aHalfLength);



struct slack_minimize_func {
  double operator()(const vect_n<double>& aX) {
    return 0.5 * aX[0] * aX[0];
  };
};

struct slack_minimize_grad {
  vect_n<double> operator()(const vect_n<double>& aX) {
    return vect_n<double>(aX[0], 0.0, 0.0, 0.0);
  };
};

struct slack_minimize_hess {
  template <typename Matrix>
  void operator()(Matrix& H, const vect_n<double>&, double, const vect_n<double>&) const {
    H = mat<double,mat_structure::nil>(4,4);
    H(0,0) = 1.0;
  };
};


struct cylinder_boundary_func {
  shared_ptr< cylinder > mCylinder;
  
  static const std::size_t size = 3;
  
  cylinder_boundary_func(const shared_ptr< cylinder >& aCylinder) : mCylinder(aCylinder) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    vect_n<double> result(3);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mCylinder->getPose().transformFromGlobal(pt);
    result[0] = pt_rel[2] + 0.5 * aX[0] * mCylinder->getLength(); // lower-bound.
    result[1] = 0.5 * aX[0] * mCylinder->getLength() - pt_rel[2]; // upper-bound.
    result[2] = aX[0] * aX[0] * mCylinder->getRadius() * mCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1];
    return result;
  };
};

struct cylinder_boundary_jac {
  shared_ptr< cylinder > mCylinder;
  
  static const std::size_t size = 3;
  
  cylinder_boundary_jac(const shared_ptr< cylinder >& aCylinder) : mCylinder(aCylinder) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mCylinder->getPose().transformFromGlobal(pt);
    
    vect<double,3> x_rel = mCylinder->getPose().rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mCylinder->getPose().rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mCylinder->getPose().rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
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
    aJac(2,0) = 2.0 * aX[0] * mCylinder->getRadius() * mCylinder->getRadius();
    aJac(2,1) = -2.0 * (pt_rel[0] * x_rel[0] + pt_rel[1] * x_rel[1]);
    aJac(2,2) = -2.0 * (pt_rel[0] * y_rel[0] + pt_rel[1] * y_rel[1]);
    aJac(2,3) = -2.0 * (pt_rel[0] * z_rel[0] + pt_rel[1] * z_rel[1]);
    
  };
};

struct ccylinder_boundary_func {
  shared_ptr< capped_cylinder > mCCylinder;
  
  static const std::size_t size = 1;
  
  ccylinder_boundary_func(const shared_ptr< capped_cylinder >& aCCylinder) : mCCylinder(aCCylinder) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    vect_n<double> result(1);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mCCylinder->getPose().transformFromGlobal(pt);
    
    if(pt_rel[2] > 0.5 * aX[0] * mCCylinder->getLength()) {
      double tmp = pt_rel[2] - 0.5 * aX[0] * mCCylinder->getLength();
      result[0] = aX[0] * aX[0] * mCCylinder->getRadius() * mCCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1] - tmp * tmp;
    } else if(pt_rel[2] < -0.5 * aX[0] * mCCylinder->getLength()) {
      double tmp = pt_rel[2] + 0.5 * aX[0] * mCCylinder->getLength();
      result[0] = aX[0] * aX[0] * mCCylinder->getRadius() * mCCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1] - tmp * tmp;
    } else {
      result[0] = aX[0] * aX[0] * mCCylinder->getRadius() * mCCylinder->getRadius() - pt_rel[0] * pt_rel[0] - pt_rel[1] * pt_rel[1]; 
    };
    
    return result;
  };
};

struct ccylinder_boundary_jac {
  shared_ptr< capped_cylinder > mCCylinder;
  
  static const std::size_t size = 1;
  
  ccylinder_boundary_jac(const shared_ptr< capped_cylinder >& aCCylinder) : mCCylinder(aCCylinder) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mCCylinder->getPose().transformFromGlobal(pt);
    
    vect<double,3> x_rel = mCCylinder->getPose().rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mCCylinder->getPose().rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mCCylinder->getPose().rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
    //aJac.resize(1,4);
    aJac(0,0) = 2.0 * aX[0] * mCCylinder->getRadius() * mCCylinder->getRadius();
    
    double tmp = 0.0;
    if(pt_rel[2] > 0.5 * aX[0] * mCCylinder->getLength()) {
      tmp = pt_rel[2] - 0.5 * aX[0] * mCCylinder->getLength();
      aJac(0,0) -= tmp * mCCylinder->getLength();
    } else if(pt_rel[2] < -0.5 * aX[0] * mCCylinder->getLength()) {
      tmp = pt_rel[2] + 0.5 * aX[0] * mCCylinder->getLength();
      aJac(0,0) += tmp * mCCylinder->getLength();
    };
    aJac(0,1) = -2.0 * (pt_rel[0] * x_rel[0] + pt_rel[1] * x_rel[1] + tmp * x_rel[2]);
    aJac(0,2) = -2.0 * (pt_rel[0] * y_rel[0] + pt_rel[1] * y_rel[1] + tmp * y_rel[2]);
    aJac(0,3) = -2.0 * (pt_rel[0] * z_rel[0] + pt_rel[1] * z_rel[1] + tmp * z_rel[2]);
    
  };
};


struct box_boundary_func {
  shared_ptr< box > mBox;
  
  static const std::size_t size = 6;
  
  box_boundary_func(const shared_ptr< box >& aBox) : mBox(aBox) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  vect_n<double> operator()(const vect_n<double>& aX) {
    vect_n<double> result(6);
    vect<double,3> pt(aX[1],aX[2],aX[3]);
    vect<double,3> pt_rel = mBox->getPose().transformFromGlobal(pt);
    result[0] = pt_rel[0] + 0.5 * aX[0] * mBox->getDimensions()[0]; // lower-bound.
    result[1] = pt_rel[1] + 0.5 * aX[0] * mBox->getDimensions()[1]; // lower-bound.
    result[2] = pt_rel[2] + 0.5 * aX[0] * mBox->getDimensions()[2]; // lower-bound.
    result[3] = 0.5 * aX[0] * mBox->getDimensions()[0] - pt_rel[0]; // upper-bound.
    result[4] = 0.5 * aX[0] * mBox->getDimensions()[1] - pt_rel[1]; // upper-bound.
    result[5] = 0.5 * aX[0] * mBox->getDimensions()[2] - pt_rel[2]; // upper-bound.
    return result;
  };
};

struct box_boundary_jac {
  shared_ptr< box > mBox;
  
  static const std::size_t size = 6;
  
  box_boundary_jac(const shared_ptr< box >& aBox) : mBox(aBox) { };
  
  // aX is the (slack (aX[0]), query-point (aX[1],aX[2],aX[3])).
  template <typename Matrix>
  void operator()(Matrix& aJac, const vect_n<double>& aX, const vect_n<double>& aH) {
    vect<double,3> x_rel = mBox->getPose().rotateFromGlobal(vect<double,3>(1.0,0.0,0.0));
    vect<double,3> y_rel = mBox->getPose().rotateFromGlobal(vect<double,3>(0.0,1.0,0.0));
    vect<double,3> z_rel = mBox->getPose().rotateFromGlobal(vect<double,3>(0.0,0.0,1.0));
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
struct dual_boundary_func {
  BoundFunc1 mBound1;
  BoundFunc2 mBound2;
  
  dual_boundary_func(BoundFunc1 aBound1, BoundFunc2 aBound2) : mBound1(aBound1), mBound2(aBound2) { };
  
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
struct dual_boundary_jac {
  BoundJac1 mBound1;
  BoundJac2 mBound2;
  
  dual_boundary_jac(BoundJac1 aBound1, BoundJac2 aBound2) : mBound1(aBound1), mBound2(aBound2) { };
  
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


};

};

#endif










