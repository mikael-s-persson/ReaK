/**
 *\file py_kinetostatics.cpp
 *
 * This source file defines export functions for the python bindings on kinetostatics classes 
 * of the ReaK platform.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date June 2012
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



#include "base/defs.hpp"

#include "gen_coord.hpp"
#include "frame_2D.hpp"
#include "frame_3D.hpp"
#include "rotations.hpp"
#include "quat_alg.hpp"

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "base/py_fixes.hpp"

#include <sstream>



namespace ReaK {
  
rot_mat_2D<double> transpose(const rot_mat_2D<double>&);
rot_mat_2D<double> invert(const rot_mat_2D<double>&);
double trace(const rot_mat_2D<double>&);
double determinant(const rot_mat_2D<double>&);

trans_mat_2D<double> invert(const trans_mat_2D<double>&);
double trace(const trans_mat_2D<double>&);
double determinant(const trans_mat_2D<double>&);

rot_mat_3D<double> transpose(const rot_mat_3D<double>&);
rot_mat_3D<double> invert(const rot_mat_3D<double>&);
double trace(const rot_mat_3D<double>&);
double determinant(const rot_mat_3D<double>&);

quaternion<double> transpose(const quaternion<double>&);
quaternion<double> invert(const quaternion<double>&);
double trace(const quaternion<double>&);
double determinant(const quaternion<double>&);

euler_angles_TB<double> transpose(const euler_angles_TB<double>&);
euler_angles_TB<double> invert(const euler_angles_TB<double>&);
double trace(const euler_angles_TB<double>&);
double determinant(const euler_angles_TB<double>&);

axis_angle<double> transpose(const axis_angle<double>&);
axis_angle<double> invert(const axis_angle<double>&);
double trace(const axis_angle<double>&);
double determinant(const axis_angle<double>&);

trans_mat_3D<double> invert(const trans_mat_3D<double>&);
double trace(const trans_mat_3D<double>&);
double determinant(const trans_mat_3D<double>&);
  
};


namespace PyReaK {

template <typename Vector>
double vect_getitem(const Vector& v, std::size_t i) {
  return v[i]; 
};

template <typename Vector>
void vect_setitem(Vector& v, std::size_t i, double d) {
  v[i] = d;
};

template <typename T>
std::string obj_to_string(const T& a) {
  std::stringstream ss;
  ss << a;
  return ss.str();
};

double ea_get_yaw(const ReaK::euler_angles_TB<double>& e) { return e.yaw(); };
void ea_set_yaw(ReaK::euler_angles_TB<double>& e, double y) { e.yaw() = y; };
double ea_get_pitch(const ReaK::euler_angles_TB<double>& e) { return e.pitch(); };
void ea_set_pitch(ReaK::euler_angles_TB<double>& e, double y) { e.pitch() = y; };
double ea_get_roll(const ReaK::euler_angles_TB<double>& e) { return e.roll(); };
void ea_set_roll(ReaK::euler_angles_TB<double>& e, double y) { e.roll() = y; };

ReaK::vect<double,3> aa_get_axis(const ReaK::axis_angle<double>& a) { return a.axis(); };
void aa_set_axis(ReaK::axis_angle<double>& a, const ReaK::vect<double,3>& y) { a.axis() = y; };
double aa_get_angle(const ReaK::axis_angle<double>& a) { return a.angle(); };
void aa_set_angle(ReaK::axis_angle<double>& a, double y) { a.angle() = y; };


void export_kinetostatics() {

  using namespace boost::python;
  
  
  
  /********************************************************************************
   *        Generalized coordinates  (gen_coord.hpp)
   * *****************************************************************************/
  
  
  class_< ReaK::gen_coord<double>,
          bases< ReaK::shared_object >,
	  ReaK::shared_ptr< ReaK::gen_coord<double> >
        >("GenCoord")
    .def(init<double,double,double,double>())
    .def_readwrite("q",&ReaK::gen_coord<double>::q)
    .def_readwrite("q_dot",&ReaK::gen_coord<double>::q_dot)
    .def_readwrite("q_ddot",&ReaK::gen_coord<double>::q_ddot)
    .def_readwrite("f",&ReaK::gen_coord<double>::f)
    .def("add_pos", &ReaK::gen_coord<double>::add_Q, return_internal_reference<>())
    .def("add_vel",&ReaK::gen_coord<double>::add_Q_dot, return_internal_reference<>())
    .def("add_acc",&ReaK::gen_coord<double>::add_Q_ddot, return_internal_reference<>())
    .def("add_force",&ReaK::gen_coord<double>::add_F, return_internal_reference<>())
    .def(self + self)
    .def(self - self)
    .def(self += self)
    .def(self -= self)
    .def("__str__",obj_to_string< ReaK::gen_coord<double> >);
  
  def("gen_coord_pos",ReaK::gen_coord_pos<double>);
  def("gen_coord_vel",ReaK::gen_coord_vel<double>);
  def("gen_coord_acc",ReaK::gen_coord_acc<double>);
  def("gen_coord_force",ReaK::gen_coord_force<double>);
  
  def("create_gen_coord",ReaK::rk_create< ReaK::gen_coord<double> >);
  def("create_gen_coord",ReaK::rk_create< ReaK::gen_coord<double>, const double&, const double&, const double&, const double& >);
  
  
  /********************************************************************************
   *        2D Rotation representations  (rotations_2D.hpp)
   * *****************************************************************************/
  
  
  class_< ReaK::rot_mat_2D<double>,
          bases< ReaK::serialization::serializable >
        >("Rotation2D")
    .def(init<double>())
    .def(init< ReaK::vect<double,2> >())
    .add_property("angle",&ReaK::rot_mat_2D<double>::getAngle,&ReaK::rot_mat_2D<double>::setAngle)
    .def("__call__",&ReaK::rot_mat_2D<double>::operator())
    .def(self * self)
    .def(self *= self)
    .def(self * other< ReaK::vect<double,2> >())
    .def(other< ReaK::vect<double,2> >() * self)
    .def("transpose",static_cast<ReaK::rot_mat_2D<double>(*)(const ReaK::rot_mat_2D<double>&) >(&ReaK::transpose))
    .def("invert",static_cast<ReaK::rot_mat_2D<double>(*)(const ReaK::rot_mat_2D<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::rot_mat_2D<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::rot_mat_2D<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::rot_mat_2D<double> >);
    
  
  class_< ReaK::trans_mat_2D<double>,
          bases< ReaK::serialization::serializable >
        >("Transform2D")
    .def(init< double, ReaK::vect<double,2> >())
    .def(init< const ReaK::rot_mat_2D<double>&, ReaK::vect<double,2> >())
    .add_property("rot_mat",&ReaK::trans_mat_2D<double>::getRotMat,&ReaK::trans_mat_2D<double>::setRotMat)
    .add_property("angle",&ReaK::trans_mat_2D<double>::getAngle,&ReaK::trans_mat_2D<double>::setAngle)
    .add_property("translation",&ReaK::trans_mat_2D<double>::getTranslation,&ReaK::trans_mat_2D<double>::setTranslation)
    .def("__call__",&ReaK::trans_mat_2D<double>::operator(),return_value_policy<copy_const_reference>())
    .def(self * self)
    .def(self *= self)
    .def(self * other< ReaK::vect<double,2> >())
    .def(self * other< ReaK::vect<double,3> >())
    .def("rotate",&ReaK::trans_mat_2D<double>::rotate)
    .def("invert",static_cast<ReaK::trans_mat_2D<double>(*)(const ReaK::trans_mat_2D<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::trans_mat_2D<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::trans_mat_2D<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::trans_mat_2D<double> >);
  
  
  /********************************************************************************
   *        2D Pose and Frame  (pose_2D.hpp and frame_2D.hpp)
   * *****************************************************************************/
  
  
  class_< ReaK::pose_2D<double>,
          bases< ReaK::shared_object >,
	  ReaK::shared_ptr< ReaK::pose_2D<double> >
        >("Pose2D")
    .def(init< ReaK::shared_ptr< ReaK::pose_2D<double> >, ReaK::vect<double,2>, ReaK::rot_mat_2D<double> >())
    .def_readwrite("parent",&ReaK::pose_2D<double>::Parent)
    .def_readwrite("position",&ReaK::pose_2D<double>::Position)
    .def_readwrite("rotation",&ReaK::pose_2D<double>::Rotation)
    .def(self * self)
    .def(self *= self)
    .def(~self)
    .def("get_global_pose",&ReaK::pose_2D<double>::getGlobalPose)
    .def("is_parent_pose",&ReaK::pose_2D<double>::isParentPose)
    .def("get_pose_relative_to",&ReaK::pose_2D<double>::getPoseRelativeTo)
    .def("rotate_to_parent",&ReaK::pose_2D<double>::rotateToParent)
    .def("rotate_to_global",&ReaK::pose_2D<double>::rotateToGlobal)
    .def("rotate_from_parent",&ReaK::pose_2D<double>::rotateFromParent)
    .def("rotate_from_global",&ReaK::pose_2D<double>::rotateFromGlobal)
    .def("transform_to_parent",&ReaK::pose_2D<double>::transformToParent)
    .def("transform_to_global",&ReaK::pose_2D<double>::transformToGlobal)
    .def("transform_from_parent",&ReaK::pose_2D<double>::transformFromParent)
    .def("transform_from_global",&ReaK::pose_2D<double>::transformFromGlobal)
    .def("add_pose_before",&ReaK::pose_2D<double>::addBefore, return_internal_reference<>())
    .def("add_pose_after",&ReaK::pose_2D<double>::addAfter, return_internal_reference<>())
    .def("translate_pose_local",&ReaK::pose_2D<double>::translateLocal, return_internal_reference<>())
    .def("translate_pose_global",&ReaK::pose_2D<double>::translateGlobal, return_internal_reference<>())
    .def("rotate_pose",&ReaK::pose_2D<double>::rotate, return_internal_reference<>())
    .def("__str__",obj_to_string< ReaK::pose_2D<double> >);
    
  class_< ReaK::weak_ptr< ReaK::pose_2D<double> > >("WeakPtrPose2D");
  implicitly_convertible< std::shared_ptr< ReaK::pose_2D<double> >, 
                          std::weak_ptr< ReaK::pose_2D<double> > >();
  
  def("create_pose_2d",ReaK::rk_create< ReaK::pose_2D<double> >);
  def("create_pose_2d",ReaK::rk_create< ReaK::pose_2D<double>, const ReaK::shared_ptr< ReaK::pose_2D<double> >&, const ReaK::vect<double,2>&, const ReaK::rot_mat_2D<double>& >);
  
  
  class_< ReaK::frame_2D<double>,
          bases< ReaK::pose_2D<double> >,
	  ReaK::shared_ptr< ReaK::frame_2D<double> >
        >("Frame2D")
    .def(init< ReaK::shared_ptr< ReaK::pose_2D<double> >, 
	       ReaK::vect<double,2>, ReaK::rot_mat_2D<double>,
	       ReaK::vect<double,2>, double,
	       ReaK::vect<double,2>, double,
	       ReaK::vect<double,2>, double >())
    .def_readwrite("velocity",&ReaK::frame_2D<double>::Velocity)
    .def_readwrite("ang_velocity",&ReaK::frame_2D<double>::AngVelocity)
    .def_readwrite("acceleration",&ReaK::frame_2D<double>::Acceleration)
    .def_readwrite("ang_acceleration",&ReaK::frame_2D<double>::AngAcceleration)
    .def_readwrite("force",&ReaK::frame_2D<double>::Force)
    .def_readwrite("torque",&ReaK::frame_2D<double>::Torque)
    .def(self * self)
    .def(self *= self)
    .def(self * other< ReaK::frame_2D<double> >())
    .def(other< ReaK::frame_2D<double> >() * self)
    .def(self *= other< ReaK::frame_2D<double> >())
    .def(~self)
    .def("get_global_frame",&ReaK::frame_2D<double>::getGlobalFrame)
    .def("get_frame_relative_to",&ReaK::frame_2D<double>::getFrameRelativeTo)
    .def("add_frame_before",&ReaK::frame_2D<double>::addBefore, return_internal_reference<>())
    .def("add_frame_after",&ReaK::frame_2D<double>::addAfter, return_internal_reference<>())
    .def("__str__",obj_to_string< ReaK::frame_2D<double> >);
    
  class_< ReaK::weak_ptr< ReaK::frame_2D<double> > >("WeakPtrFrame2D");
  implicitly_convertible< std::shared_ptr< ReaK::frame_2D<double> >, 
                          std::weak_ptr< ReaK::frame_2D<double> > >();
  
  def("create_frame_2d",ReaK::rk_create< ReaK::frame_2D<double> >);
  def("create_frame_2d",ReaK::rk_create< ReaK::frame_2D<double>, 
       const ReaK::shared_ptr< ReaK::pose_2D<double> >&, 
       const ReaK::vect<double,2>&, const ReaK::rot_mat_2D<double>&, 
       const ReaK::vect<double,2>&, const double&, 
       const ReaK::vect<double,2>&, const double&, 
       const ReaK::vect<double,2>&, const double& >);
  
  
  /********************************************************************************
   *        3D Rotation representations  (rotations_3D.hpp)
   * *****************************************************************************/
  
  
  class_< ReaK::rot_mat_3D<double>,
          bases< ReaK::serialization::serializable >
        >("Rotation3D")
    .def("__call__",&ReaK::rot_mat_3D<double>::operator())
    .def(self * self)
    .def(self *= self)
    .def(self * other< ReaK::vect<double,3> >())
    .def(other< ReaK::vect<double,3> >() * self)
    .def("transpose",static_cast<ReaK::rot_mat_3D<double>(*)(const ReaK::rot_mat_3D<double>&) >(&ReaK::transpose))
    .def("invert",static_cast<ReaK::rot_mat_3D<double>(*)(const ReaK::rot_mat_3D<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::rot_mat_3D<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::rot_mat_3D<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::rot_mat_3D<double> >);
  
  
  class_< ReaK::quaternion<double>,
          bases< ReaK::serialization::serializable >
        >("QuaternionRot")
    .def(init< ReaK::vect<double,4> >())
    .def(init< ReaK::rot_mat_3D<double> >())
    .def("__getitem__",&ReaK::quaternion<double>::operator[],return_value_policy<copy_const_reference>())
    .def("get_rotation_matrix",&ReaK::quaternion<double>::getRotMat)
    .def(self * self)
    .def(self * other< ReaK::rot_mat_3D<double> >())
    .def(other< ReaK::rot_mat_3D<double> >() * self)
    .def(self *= self)
    .def(self *= other< ReaK::rot_mat_3D<double> >())
    .def(self *= other< ReaK::axis_angle<double> >())
    .def(self *= other< ReaK::euler_angles_TB<double> >())
    .def(self * other< ReaK::vect<double,3> >())
    .def("transpose",static_cast<ReaK::quaternion<double>(*)(const ReaK::quaternion<double>&) >(&ReaK::transpose))
    .def("invert",static_cast<ReaK::quaternion<double>(*)(const ReaK::quaternion<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::quaternion<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::quaternion<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::quaternion<double> >);
  
  
  class_< ReaK::euler_angles_TB<double>,
          bases< ReaK::serialization::serializable >
        >("EulerAngles")
    .def(init< double, double, double >())
    .def(init< ReaK::rot_mat_3D<double> >())
    .def(init< ReaK::quaternion<double> >())
    .add_property("yaw",&ea_get_yaw,&ea_set_yaw)
    .add_property("pitch",&ea_get_pitch,&ea_set_pitch)
    .add_property("roll",&ea_get_roll,&ea_set_roll)
    .def("get_rotation_matrix",&ReaK::euler_angles_TB<double>::getRotMat)
    .def("get_quaternion",&ReaK::euler_angles_TB<double>::getQuaternion)
    .def(self * self)
    .def(self * other< ReaK::rot_mat_3D<double> >())
    .def(self * other< ReaK::quaternion<double> >())
    .def(other< ReaK::rot_mat_3D<double> >() * self)
    .def(other< ReaK::quaternion<double> >() * self)
    .def(self *= self)
    .def(self *= other< ReaK::rot_mat_3D<double> >())
    .def(self *= other< ReaK::axis_angle<double> >())
    .def(self *= other< ReaK::quaternion<double> >())
    .def(self * other< ReaK::vect<double,3> >())
    .def("transpose",static_cast<ReaK::euler_angles_TB<double>(*)(const ReaK::euler_angles_TB<double>&) >(&ReaK::transpose))
    .def("invert",static_cast<ReaK::euler_angles_TB<double>(*)(const ReaK::euler_angles_TB<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::euler_angles_TB<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::euler_angles_TB<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::euler_angles_TB<double> >);
  
  
  class_< ReaK::axis_angle<double>,
          bases< ReaK::serialization::serializable >
        >("AxisAngle")
    .def(init< double, ReaK::vect<double,3> >())
    .def(init< ReaK::rot_mat_3D<double> >())
    .def(init< ReaK::quaternion<double> >())
    .def(init< ReaK::euler_angles_TB<double> >())
    .add_property("angle",&aa_get_angle,&aa_set_angle)
    .add_property("axis",&aa_get_axis,&aa_set_axis)
    .def("get_rotation_matrix",&ReaK::axis_angle<double>::getRotMat)
    .def("get_quaternion",&ReaK::axis_angle<double>::getQuaternion)
    .def("get_quaternion",&ReaK::axis_angle<double>::getEulerAnglesTB)
    .def(self * self)
    .def(self * other< ReaK::rot_mat_3D<double> >())
    .def(self * other< ReaK::euler_angles_TB<double> >())
    .def(self * other< ReaK::quaternion<double> >())
    .def(other< ReaK::rot_mat_3D<double> >() * self)
    .def(other< ReaK::euler_angles_TB<double> >() * self)
    .def(other< ReaK::quaternion<double> >() * self)
    .def(self *= self)
    .def(self *= other< ReaK::rot_mat_3D<double> >())
    .def(self *= other< ReaK::euler_angles_TB<double> >())
    .def(self *= other< ReaK::quaternion<double> >())
    .def(self * other< ReaK::vect<double,3> >())
    .def("transpose",static_cast<ReaK::axis_angle<double>(*)(const ReaK::axis_angle<double>&) >(&ReaK::transpose))
    .def("invert",static_cast<ReaK::axis_angle<double>(*)(const ReaK::axis_angle<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::axis_angle<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::axis_angle<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::axis_angle<double> >);
  
  
  class_< ReaK::trans_mat_3D<double>,
          bases< ReaK::serialization::serializable >
        >("Transform3D")
    .def(init< ReaK::quaternion<double>, ReaK::vect<double,3> >())
    .def(init< ReaK::rot_mat_3D<double>, ReaK::vect<double,3> >())
    .def(init< ReaK::euler_angles_TB<double>, ReaK::vect<double,3> >())
    .def(init< ReaK::axis_angle<double>, ReaK::vect<double,3> >())
    .add_property("rot_mat",&ReaK::trans_mat_3D<double>::getRotMat,&ReaK::trans_mat_3D<double>::setRotMat)
    .add_property("quaternion",&ReaK::trans_mat_3D<double>::getQuaternion,&ReaK::trans_mat_3D<double>::setQuaternion)
    .add_property("euler_angles",&ReaK::trans_mat_3D<double>::getEulerAnglesTB,&ReaK::trans_mat_3D<double>::setEulerAnglesTB)
    .add_property("axis_angle",&ReaK::trans_mat_3D<double>::getAxisAngle,&ReaK::trans_mat_3D<double>::setAxisAngle)
    .add_property("translation",&ReaK::trans_mat_3D<double>::getTranslation,&ReaK::trans_mat_3D<double>::setTranslation)
    .def("__call__",&ReaK::trans_mat_3D<double>::operator(),return_value_policy<copy_const_reference>())
    .def(self * self)
    .def(self * other< ReaK::rot_mat_3D<double> >())
    .def(self * other< ReaK::euler_angles_TB<double> >())
    .def(self * other< ReaK::axis_angle<double> >())
    .def(self * other< ReaK::quaternion<double> >())
    .def(other< ReaK::rot_mat_3D<double> >() * self)
    .def(other< ReaK::euler_angles_TB<double> >() * self)
    .def(other< ReaK::axis_angle<double> >() * self)
    .def(other< ReaK::quaternion<double> >() * self)
    .def(self *= self)
    .def(self *= other< ReaK::rot_mat_3D<double> >())
    .def(self *= other< ReaK::euler_angles_TB<double> >())
    .def(self *= other< ReaK::axis_angle<double> >())
    .def(self *= other< ReaK::quaternion<double> >())
    .def(self * other< ReaK::vect<double,3> >())
    .def(self * other< ReaK::vect<double,4> >())
    .def("rotate",&ReaK::trans_mat_3D<double>::rotate)
    .def("invert",static_cast<ReaK::trans_mat_3D<double>(*)(const ReaK::trans_mat_3D<double>&) >(&ReaK::invert))
    .def("trace",static_cast<double(*)(const ReaK::trans_mat_3D<double>&) >(&ReaK::trace))
    .def("determinant",static_cast<double(*)(const ReaK::trans_mat_3D<double>&) >(&ReaK::determinant))
    .def("__str__",obj_to_string< ReaK::trans_mat_3D<double> >);
  
  
  
  /********************************************************************************
   *        3D Pose and Frame  (pose_3D.hpp and frame_3D.hpp)
   * *****************************************************************************/
  
  
  class_< ReaK::pose_3D<double>,
          bases< ReaK::shared_object >,
	  ReaK::shared_ptr< ReaK::pose_3D<double> >
        >("Pose3D")
    .def(init< ReaK::shared_ptr< ReaK::pose_3D<double> >, ReaK::vect<double,3>, ReaK::quaternion<double> >())
    .def_readwrite("parent",&ReaK::pose_3D<double>::Parent)
    .def_readwrite("position",&ReaK::pose_3D<double>::Position)
    .def_readwrite("quaternion",&ReaK::pose_3D<double>::Quat)
    .def(self * self)
    .def(self *= self)
    .def(~self)
    .def("get_global_pose",&ReaK::pose_3D<double>::getGlobalPose)
    .def("is_parent_pose",&ReaK::pose_3D<double>::isParentPose)
    .def("get_pose_relative_to",&ReaK::pose_3D<double>::getPoseRelativeTo)
    .def("rotate_to_parent",&ReaK::pose_3D<double>::rotateToParent)
    .def("rotate_to_global",&ReaK::pose_3D<double>::rotateToGlobal)
    .def("rotate_from_parent",&ReaK::pose_3D<double>::rotateFromParent)
    .def("rotate_from_global",&ReaK::pose_3D<double>::rotateFromGlobal)
    .def("transform_to_parent",&ReaK::pose_3D<double>::transformToParent)
    .def("transform_to_global",&ReaK::pose_3D<double>::transformToGlobal)
    .def("transform_from_parent",&ReaK::pose_3D<double>::transformFromParent)
    .def("transform_from_global",&ReaK::pose_3D<double>::transformFromGlobal)
    .def("add_pose_before",&ReaK::pose_3D<double>::addBefore, return_internal_reference<>())
    .def("add_pose_after",&ReaK::pose_3D<double>::addAfter, return_internal_reference<>())
    .def("translate_pose_local",&ReaK::pose_3D<double>::translateLocal, return_internal_reference<>())
    .def("translate_pose_global",&ReaK::pose_3D<double>::translateGlobal, return_internal_reference<>())
    .def("rotate_pose_local",&ReaK::pose_3D<double>::rotateLocal, return_internal_reference<>())
    .def("rotate_pose_global",&ReaK::pose_3D<double>::rotateGlobal, return_internal_reference<>())
    .def("__str__",obj_to_string< ReaK::pose_3D<double> >);
    
  class_< ReaK::weak_ptr< ReaK::pose_3D<double> > >("WeakPtrPose3D");
  implicitly_convertible< std::shared_ptr< ReaK::pose_3D<double> >, 
                          std::weak_ptr< ReaK::pose_3D<double> > >();
  
  def("create_pose_3d",ReaK::rk_create< ReaK::pose_3D<double> >);
  def("create_pose_3d",ReaK::rk_create< ReaK::pose_3D<double>, 
        const ReaK::shared_ptr< ReaK::pose_3D<double> >&, 
        const ReaK::vect<double,3>&, const ReaK::quaternion<double>& >);
  
  
  class_< ReaK::frame_3D<double>,
          bases< ReaK::pose_3D<double> >,
	  ReaK::shared_ptr< ReaK::frame_3D<double> >
        >("Frame3D")
    .def(init< ReaK::shared_ptr< ReaK::pose_3D<double> >, 
	       ReaK::vect<double,3>, ReaK::quaternion<double>,
	       ReaK::vect<double,3>, ReaK::vect<double,3>,
	       ReaK::vect<double,3>, ReaK::vect<double,3>,
	       ReaK::vect<double,3>, ReaK::vect<double,3> >())
    .def_readwrite("velocity",&ReaK::frame_3D<double>::Velocity)
    .def_readwrite("ang_velocity",&ReaK::frame_3D<double>::AngVelocity)
    .def_readwrite("acceleration",&ReaK::frame_3D<double>::Acceleration)
    .def_readwrite("ang_acceleration",&ReaK::frame_3D<double>::AngAcceleration)
    .def_readwrite("force",&ReaK::frame_3D<double>::Force)
    .def_readwrite("torque",&ReaK::frame_3D<double>::Torque)
    .def(self * self)
    .def(self *= self)
    .def(self * other< ReaK::frame_3D<double> >())
    .def(other< ReaK::frame_3D<double> >() * self)
    .def(self *= other< ReaK::frame_3D<double> >())
    .def(~self)
    .def("get_global_frame",&ReaK::frame_3D<double>::getGlobalFrame)
    .def("get_frame_relative_to",&ReaK::frame_3D<double>::getFrameRelativeTo)
    .def("add_frame_before",static_cast<ReaK::frame_3D<double>&(ReaK::frame_3D<double>::*)(const ReaK::pose_3D<double>&)>(&ReaK::frame_3D<double>::addBefore), return_internal_reference<>())
    .def("add_frame_after",static_cast<ReaK::frame_3D<double>&(ReaK::frame_3D<double>::*)(const ReaK::pose_3D<double>&)>(&ReaK::frame_3D<double>::addAfter), return_internal_reference<>())
    .def("add_frame_before",static_cast<ReaK::frame_3D<double>&(ReaK::frame_3D<double>::*)(const ReaK::frame_3D<double>&)>(&ReaK::frame_3D<double>::addBefore), return_internal_reference<>())
    .def("add_frame_after",static_cast<ReaK::frame_3D<double>&(ReaK::frame_3D<double>::*)(const ReaK::frame_3D<double>&)>(&ReaK::frame_3D<double>::addAfter), return_internal_reference<>())
    .def("__str__",obj_to_string< ReaK::frame_3D<double> >);
    
  class_< ReaK::weak_ptr< ReaK::frame_3D<double> > >("WeakPtrFrame3D");
  implicitly_convertible< std::shared_ptr< ReaK::frame_3D<double> >, 
                          std::weak_ptr< ReaK::frame_3D<double> > >();
  
  def("create_frame_3d",ReaK::rk_create< ReaK::frame_3D<double> >);
  def("create_frame_3d",ReaK::rk_create< ReaK::frame_3D<double>, 
       const ReaK::shared_ptr< ReaK::pose_3D<double> >&, 
       const ReaK::vect<double,3>&, const ReaK::quaternion<double>&, 
       const ReaK::vect<double,3>&, const ReaK::vect<double,3>&, 
       const ReaK::vect<double,3>&, const ReaK::vect<double,3>&, 
       const ReaK::vect<double,3>&, const ReaK::vect<double,3>& >);
  
  
  
#ifdef RK_ENABLE_CXX0X_FEATURES
  implicitly_convertible< std::shared_ptr< ReaK::gen_coord<double> >, 
                          std::shared_ptr< ReaK::shared_object > >();
  implicitly_convertible< std::shared_ptr< ReaK::pose_2D<double> >, 
                          std::shared_ptr< ReaK::shared_object > >();
  implicitly_convertible< std::shared_ptr< ReaK::frame_2D<double> >, 
                          std::shared_ptr< ReaK::pose_2D<double> > >();
  implicitly_convertible< std::shared_ptr< ReaK::pose_3D<double> >, 
                          std::shared_ptr< ReaK::shared_object > >();
  implicitly_convertible< std::shared_ptr< ReaK::frame_3D<double> >, 
                          std::shared_ptr< ReaK::pose_3D<double> > >();
#endif
  
};


};
















