
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

#include "manip_clik_calculator.hpp"

#include "optimization/nl_interior_points_methods.hpp"
#include "optimization/function_types.hpp"

namespace ReaK {

namespace kte {
  

namespace detail {
  
  
struct clik_ineq_function {
  const manip_clik_calculator* parent;
  
  clik_ineq_function(const manip_clik_calculator* aParent) : parent(aParent) { };
  
  vect_n<double> operator()(const vect_n<double>& x) const {
    std::size_t l_size = 0;
    for(std::size_t i = 0; i < parent->lower_bounds.size(); ++i) {
      if(parent->lower_bounds[i] != -std::numeric_limits<double>::infinity())
        ++l_size;
    };
    std::size_t u_size = 0;
    for(std::size_t i = 0; i < parent->upper_bounds.size(); ++i) {
      if(parent->upper_bounds[i] != std::numeric_limits<double>::infinity())
        ++u_size;
    };
    vect_n<double> result(l_size + u_size);
    std::size_t j = 0;
    for(std::size_t i = 0; (i < parent->lower_bounds.size()) && (i < x.size()); ++i) {
      if(parent->lower_bounds[i] != -std::numeric_limits<double>::infinity()) {
        result[j] = x[i] - parent->lower_bounds[i];
        ++j;
      };
    };
    for(std::size_t i = 0; (i < parent->upper_bounds.size()) && (i < x.size()); ++i) {
      if(parent->upper_bounds[i] != std::numeric_limits<double>::infinity()) {
        result[j] = parent->upper_bounds[i] - x[i];
        ++j;
      };
    };
    return result;
  };
  
};

struct clik_ineq_jac_filler {
  const manip_clik_calculator* parent;
  
  clik_ineq_jac_filler(const manip_clik_calculator* aParent) : parent(aParent) { };
  
  template <typename Matrix>
  void operator()(Matrix& J, const vect_n<double>& x, const vect_n<double>& h) const {
    J = mat<double,mat_structure::nil>(h.size(), x.size());
    std::size_t j = 0;
    for(std::size_t i = 0; (i < parent->lower_bounds.size()) && (i < x.size()); ++i) {
      if(parent->lower_bounds[i] != -std::numeric_limits<double>::infinity()) {
        J(j,i) = 1.0;
        ++j;
      };
    };
    for(std::size_t i = 0; (i < parent->upper_bounds.size()) && (i < x.size()); ++i) {
      if(parent->upper_bounds[i] != std::numeric_limits<double>::infinity()) {
        J(j,i) = -1.0;
        ++j;
      };
    };
  };
  
};

struct clik_eq_function {
  const manip_clik_calculator* parent;
  
  clik_eq_function(const manip_clik_calculator* aParent) : parent(aParent) { };
  
  vect_n<double> operator()(const vect_n<double>& x) const {
    
    direct_kinematics_model* pmdl = parent->model;
    
    const std::vector< shared_ptr< joint_dependent_gen_coord > >& dep_gen_coords = pmdl->getDependentCoord();
    const std::vector< shared_ptr< joint_dependent_frame_2D > >& dep_frames_2D = pmdl->getDependentFrame2D();
    const std::vector< shared_ptr< joint_dependent_frame_3D > >& dep_frames_3D = pmdl->DependentFrames3D();
    
    if( ( pmdl->getDependentCoordsCount() != parent->desired_gen_coords.size() ) ||
        ( pmdl->getDependentFrames2DCount() != parent->desired_frame_2D.size() ) ||
        ( pmdl->getDependentFrames3DCount() != parent->desired_frame_3D.size() ) ) 
      throw std::range_error("Improper inverse-kinematics problem, the number of desired frames does not match the number of end-effector frames!");
    
    manip_kin_mdl_joint_io(pmdl).setJointPositions(&x[0]);
    manip_kin_mdl_joint_io(pmdl).setJointVelocities(&x[pmdl->getJointPositionsCount()]);
    
    pmdl->doMotion();
    
    vect_n<double> result(pmdl->getDependentVelocitiesCount() * 2
                        + pmdl->Frames2D().size() + pmdl->Frames3D().size());
    
    
    // enforce the desired 'end-effector' frames.
    std::size_t j = 0;
    std::size_t k = pmdl->getDependentVelocitiesCount();
    
    for(std::size_t i = 0; i < pmdl->getDependentCoordsCount(); ++i) {
      result[j] = pmdl->getDependentCoord(i)->mFrame->q - parent->desired_gen_coords[i].q; ++j;
      result[k] = pmdl->getDependentCoord(i)->mFrame->q_dot - parent->desired_gen_coords[i].q_dot; ++k;
    };
    
    for(std::size_t i = 0; i < pmdl->getDependentFrames2DCount(); ++i) {
      frame_2D<double> err = parent->desired_frame_2D[i].getFrameRelativeTo(pmdl->getDependentFrame2D(i)->mFrame);
      result[j] = -err.Position[0]; ++j;
      result[j] = -err.Position[1]; ++j;
      result[j] = -err.Rotation.getAngle(); ++j;
      result[k] = -err.Velocity[0]; ++k;
      result[k] = -err.Velocity[1]; ++k;
      result[k] = -err.AngVelocity; ++k;
    };
    
    for(std::size_t i = 0; i < pmdl->getDependentFrames3DCount(); ++i) {
      frame_3D<double> err = parent->desired_frame_3D[i].getFrameRelativeTo(pmdl->getDependentFrame3D(i)->mFrame);
      result[j] = -err.Position[0]; ++j;
      result[j] = -err.Position[1]; ++j;
      result[j] = -err.Position[2]; ++j;
      axis_angle<double> aa = axis_angle<double>(err.Quat);
      vect<double,3> v = aa.angle() * aa.axis(); 
      result[j] = -v[0]; ++j;
      result[j] = -v[1]; ++j;
      result[j] = -v[2]; ++j;
      result[k] = -err.Velocity[0]; ++k;
      result[k] = -err.Velocity[1]; ++k;
      result[k] = -err.Velocity[2]; ++k;
      result[k] = -err.AngVelocity[0]; ++k;
      result[k] = -err.AngVelocity[1]; ++k;
      result[k] = -err.AngVelocity[2]; ++k;
    };
    
    // enforce the normality of the rotation representation.
    j = pmdl->getCoordsCount();
    for(std::size_t i = 0; i < pmdl->getFrames2DCount(); ++i) {
      j += 2;
      result[k] = 1.0 - x[j] * x[j] - x[j+1] * x[j+1]; ++k;
      j += 2;
    };
    
    for(std::size_t i = 0; i < pmdl->getFrames3DCount(); ++i) {
      j += 3;
      result[k] = 1.0 - x[j] * x[j] - x[j+1] * x[j+1] - x[j+2] * x[j+2] - x[j+3] * x[j+3]; ++k;
      j += 4;
    };
    
    return result;
  };
  
};

struct clik_eq_jac_filler {
  const manip_clik_calculator* parent;
  mat<double,mat_structure::rectangular> jac_tmp;
  mat<double,mat_structure::rectangular> jacdot_tmp;
  
  clik_eq_jac_filler(const manip_clik_calculator* aParent) : parent(aParent) { };
  
  template <typename Matrix>
  void operator()(Matrix& J, const vect_n<double>& x, const vect_n<double>& h) const {
    
    direct_kinematics_model* pmdl = parent->model;
    
    if( ( pmdl->getDependentCoordsCount() != parent->desired_gen_coords.size() ) ||
        ( pmdl->getDependentFrames2DCount() != parent->desired_frame_2D.size() ) ||
        ( pmdl->getDependentFrames3DCount() != parent->desired_frame_3D.size() ) ) 
      throw std::range_error("Improper inverse-kinematics problem, the number of desired frames does not match the number of end-effector frames!");
    
    manip_kin_mdl_joint_io(pmdl).setJointPositions(&x[0]);
    manip_kin_mdl_joint_io(pmdl).setJointVelocities(&x[pmdl->getJointPositionsCount()]);
    
    pmdl->doDirectMotion();
    
    J = mat<double,mat_structure::nil>(h.size(), x.size());
    
    manip_kin_mdl_jac_calculator(pmdl).getJacobianMatrixAndDerivative(jac_tmp, jacdot_tmp);
    
    std::size_t n1 = pmdl->getDependentVelocitiesCount();
    std::size_t m1 = (x.size() + pmdl->getFrames2DCount() + pmdl->Frames3D().size()) / 2;
    for(std::size_t i = 0; i < (x.size() - pmdl->getFrames2DCount() - pmdl->Frames3D().size()) / 2; ++i) {
      for(std::size_t j = 0; j < n1; ++j)
        J(j, i) = jac_tmp(j, i);
      for(std::size_t j = 0; j < n1; ++j)
        J(j + n1, i) = jacdot_tmp(j, i);
      for(std::size_t j = 0; j < n1; ++j)
        J(j + n1, i + m1) = jac_tmp(j, i);
    };
    
    // j is the index to the last position element of x.
    std::size_t j = (x.size() + pmdl->getFrames2DCount() + pmdl->getFrames3DCount()) / 2 - 1;
    // k is the index to the last valid column of J.
    std::size_t k = (x.size() - pmdl->getFrames2DCount() - pmdl->getFrames3DCount()) / 2 - 1;
    // l is the index to the last normality-constraint row of J.
    std::size_t l = h.size() - 1;
    
    for(std::size_t i = 0; i < pmdl->getFrames3DCount(); ++i) {
      mat<double, mat_structure::rectangular> H_inv(3,4);  // NOTE : Check this again, shouldn't there be a factor of 2 or 0.5 ???
      H_inv(0,0) = -x[j - 2]; H_inv(1,0) = -x[j - 1]; H_inv(2,0) = -x[j];
      H_inv(0,1) =  x[j - 3]; H_inv(1,1) = -x[j]; H_inv(2,1) = x[j - 1];
      H_inv(0,2) =  x[j]; H_inv(1,2) =  x[j - 3]; H_inv(2,2) = -x[j - 2];
      H_inv(0,3) =  -x[j - 1]; H_inv(1,3) =  x[j - 2]; H_inv(2,3) = x[j - 3];
      
      // apply the transformation from quat_dot to omega:
      sub(J)(range(0,n1 * 2 - 1),
             range(j-3,j)) = sub(J)(range(0,n1 * 2 - 1),
                                    range(k-2,k)) * H_inv;
      // fill in the normality-constraint jacobians:
      J(l,j-3) = -2.0 * x[j-3];
      J(l,j-2) = -2.0 * x[j-2];
      J(l,j-1) = -2.0 * x[j-1];
      J(l,j) = -2.0 * x[j];
      k -= 3;
      j -= 4;
      --l;
      // copy the position row:
      for(std::size_t r = 0; r < 3; ++r)
        for(std::size_t s = 0; s < n1 * 2; ++s)
          J(s,j-r) = J(s,k-r);
      k -= 3;
      j -= 3;
      
    };
    
    for(std::size_t i = 0; i < pmdl->getFrames2DCount(); ++i) {
      // apply the transformation from quat_dot to omega:
      for(std::size_t s = 0; s < n1 * 2; ++s) {
        J(s,j-1) = -J(s,k) * x[j];
        J(s,j) = J(s,k) * x[j-1];
      };
      // fill in the normality-constraint jacobians:
      J(l,j-1) = -2.0 * x[j-1];
      J(l,j) = -2.0 * x[j];
      k -= 1;
      j -= 2;
      --l;
      // copy the position row:
      for(std::size_t r = 0; r < 2; ++r)
        for(std::size_t s = 0; s < n1 * 2; ++s)
          J(s,j-r) = J(s,k-r);
      k -= 2;
      i -= 2;
      
    };
    
  };
  
};

  
};



void manip_clik_calculator::readDesiredFromModel() {
  
  desired_gen_coords.resize( model->getDependentCoordsCount() );
  for(std::size_t i = 0; i < desired_gen_coords.size(); ++i)
    desired_gen_coords[i] = *(model->getDependentCoord(i)->mFrame);
  
  desired_frame_2D.resize( model->getDependentFrames2DCount() );
  for(std::size_t i = 0; i < desired_frame_2D.size(); ++i)
    desired_frame_2D[i] = *(model->getDependentFrame2D(i)->mFrame);
  
  desired_frame_3D.resize( model->getDependentFrames3DCount() );
  for(std::size_t i = 0; i < desired_frame_3D.size(); ++i)
    desired_frame_3D[i] = *(model->getDependentFrame3D(i)->mFrame);
  
};


vect_n<double> manip_clik_calculator::readJointStatesFromModel() {
  vect_n<double> x(model->getJointPositionsCount() + model->getJointVelocitiesCount());
  
  manip_kin_mdl_joint_io(model).getJointPositions(&x[0]);
  manip_kin_mdl_joint_io(model).getJointVelocities(&x[model->getJointPositionsCount()]);
  
  return x;
};

void manip_clik_calculator::writeJointStatesToModel(const vect_n<double>& x) const {
  manip_kin_mdl_joint_io(model).setJointPositions(&x[0]);
  manip_kin_mdl_joint_io(model).setJointVelocities(&x[model->getJointPositionsCount()]);
};

void manip_clik_calculator::runOptimizer( vect_n<double>& x ) {
  shared_ptr<const optim::cost_evaluator> tmp_cost_eval = cost_eval;
  if(!cost_eval) {
    tmp_cost_eval = shared_ptr<const optim::cost_evaluator>(
      new optim::quadratic_cost_evaluator(
        vect_n<double>(x.size(), double(0.0)),
        mat<double,mat_structure::symmetric>(
          ( mat<double,mat_structure::nil>(model->getJointPositionsCount(), model->getJointPositionsCount() + model->getJointVelocitiesCount()) |
          ( mat<double,mat_structure::nil>(model->getJointVelocitiesCount(), model->getJointPositionsCount()) & mat<double,mat_structure::identity>(model->getJointVelocitiesCount()) ) )
        ) 
      )
    );
  };
  
  typedef optim::nlip_newton_tr_factory<optim::oop_cost_function,
                                        optim::oop_cost_grad,
                                        optim::oop_cost_hess,
                                        double,
                                        detail::clik_eq_function,
                                        detail::clik_eq_jac_filler,
                                        detail::clik_ineq_function,
                                        detail::clik_ineq_jac_filler> optim_factory_type;
//  typedef optim::nlip_quasi_newton_tr_factory<optim::oop_cost_function,
//                                              optim::oop_cost_grad,
//                                              double,
//                                              eq_function,
//                                              eq_jac_filler,
//                                              ineq_function,
//                                              ineq_jac_filler> optim_factory_type;
  
  optim_factory_type optimizer = optim_factory_type(
    optim::oop_cost_function(tmp_cost_eval),
    optim::oop_cost_grad(tmp_cost_eval),
    optim::oop_cost_hess(tmp_cost_eval),
    max_radius, mu, max_iter,
    eq_function(this), eq_jac_filler(this),
    ineq_function(this), ineq_jac_filler(this),
    tol, eta, tau);
  
  optimizer( x );
};

void manip_clik_calculator::solveInverseKinematics() {
  
  readDesiredFromModel();
  
  vect_n<double> x = readJointStatesFromModel();
  
  runOptimizer(x);
  
  writeJointStatesToModel(x);
};



};

};










