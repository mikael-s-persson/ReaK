
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

#ifndef MOTION_JACOBIANS_HPP
#define MOTION_JACOBIANS_HPP

#include "kinetostatics.hpp"

namespace ReaK {
  

template <typename T> 
class jacobian_gen_gen : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_gen_gen<T> self;
  
  value_type qd_qd;
  value_type qd_qdd;
  
  jacobian_gen_gen(const value_type& aQdQd = value_type(), 
		   const value_type& aQdQdd = value_type()) : qd_qd(aQdQd), qd_qdd(aQdQdd) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(qd_qd)
      & RK_SERIAL_SAVE_WITH_NAME(qd_qdd);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(qd_qd)
      & RK_SERIAL_LOAD_WITH_NAME(qd_qdd);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000022,1,"jacobian_gen_gen",shared_object)
};

template <typename T> 
class jacobian_gen_2D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_gen_2D<T> self;
  
  boost::weak_ptr< frame_2D<value_type> > Parent;
  vect<value_type,2> qd_vel;
  value_type qd_avel;
  vect<value_type,2> qd_acc;
  value_type qd_aacc;
  
  jacobian_gen_2D(boost::weak_ptr< frame_2D<value_type> > aParent = boost::weak_ptr< frame_2D<value_type> >(),
                  const vect<value_type,2>& aQdVel = vect<value_type,2>(),
		  const value_type& aQdAVel = value_type(), 
		  const vect<value_type,2>& aQdAcc = vect<value_type,2>(),
		  const value_type& aQdAAcc = value_type()) : 
		  Parent(aParent),
		  qd_vel(aQdVel),
		  qd_avel(aQdAVel),
		  qd_acc(aQdAcc),
		  qd_aacc(aQdAAcc) { };
  
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(qd_vel)
      & RK_SERIAL_SAVE_WITH_NAME(qd_avel)
      & RK_SERIAL_SAVE_WITH_NAME(qd_acc)
      & RK_SERIAL_SAVE_WITH_NAME(qd_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_2D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(qd_vel)
      & RK_SERIAL_LOAD_WITH_NAME(qd_avel)
      & RK_SERIAL_LOAD_WITH_NAME(qd_acc)
      & RK_SERIAL_LOAD_WITH_NAME(qd_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000023,1,"jacobian_gen_2D",shared_object)
};

template <typename T> 
class jacobian_gen_3D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_gen_3D<T> self;
  
  boost::weak_ptr< frame_3D<value_type> > Parent;
  vect<value_type,3> qd_vel;
  vect<value_type,3> qd_avel;
  vect<value_type,3> qd_acc;
  vect<value_type,3> qd_aacc;
  
  jacobian_gen_3D(boost::weak_ptr< frame_3D<value_type> > aParent = boost::weak_ptr< frame_3D<value_type> >(),
                  const vect<value_type,3>& aQdVel = vect<value_type,3>(),
		  const vect<value_type,3>& aQdAVel = vect<value_type,3>(), 
		  const vect<value_type,3>& aQdAcc = vect<value_type,3>(),
		  const vect<value_type,3>& aQdAAcc = vect<value_type,3>()) : 
		  Parent(aParent),
		  qd_vel(aQdVel),
		  qd_avel(aQdAVel),
		  qd_acc(aQdAcc),
		  qd_aacc(aQdAAcc) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(qd_vel)
      & RK_SERIAL_SAVE_WITH_NAME(qd_avel)
      & RK_SERIAL_SAVE_WITH_NAME(qd_acc)
      & RK_SERIAL_SAVE_WITH_NAME(qd_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_3D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(qd_vel)
      & RK_SERIAL_LOAD_WITH_NAME(qd_avel)
      & RK_SERIAL_LOAD_WITH_NAME(qd_acc)
      & RK_SERIAL_LOAD_WITH_NAME(qd_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000024,1,"jacobian_gen_3D",shared_object)
};





template <typename T> 
class jacobian_2D_gen : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_2D_gen<T> self;
  
  vect<value_type,2> vel_qd;
  value_type avel_qd;
  vect<value_type,2> vel_qdd;
  value_type avel_qdd;
  
  jacobian_2D_gen(const vect<value_type,2>& aVelQd = vect<value_type,2>(),
		  const value_type& aAVelQd = value_type(), 
		  const vect<value_type,2>& aVelQdd = vect<value_type,2>(),
		  const value_type& aAVelQdd = value_type()) : 
		  vel_qd(aVelQd),
		  avel_qd(aAVelQd),
		  vel_qdd(aVelQdd),
		  avel_qdd(aAVelQdd) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(vel_qd)
      & RK_SERIAL_SAVE_WITH_NAME(avel_qd)
      & RK_SERIAL_SAVE_WITH_NAME(vel_qdd)
      & RK_SERIAL_SAVE_WITH_NAME(avel_qdd);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(vel_qd)
      & RK_SERIAL_LOAD_WITH_NAME(avel_qd)
      & RK_SERIAL_LOAD_WITH_NAME(vel_qdd)
      & RK_SERIAL_LOAD_WITH_NAME(avel_qdd);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000025,1,"jacobian_2D_gen",shared_object)
};

template <typename T> 
class jacobian_2D_2D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_2D_2D<T> self;
  
  boost::weak_ptr< frame_2D<value_type> > Parent;
  vect<vect<value_type,2>,2> vel_vel;
  vect<value_type,2> vel_avel;
  vect<value_type,2> avel_vel;
  value_type avel_avel;
  vect<vect<value_type,2>,2> vel_acc;
  vect<value_type,2> vel_aacc;
  vect<value_type,2> avel_acc;
  value_type avel_aacc;
  
  jacobian_2D_2D(boost::weak_ptr< frame_2D<value_type> > aParent = boost::weak_ptr< frame_2D<value_type> >(),
                 const vect<vect<value_type,2>,2>& aVelVel = vect<vect<value_type,2>,2>(),
		 const vect<value_type,2>& aVelAVel = vect<value_type,2>(),
		 const vect<value_type,2>& aAVelVel = vect<value_type,2>(),
		 const value_type& aAVelAVel = value_type(), 
		 const vect<vect<value_type,2>,2>& aVelAcc = vect<vect<value_type,2>,2>(),
		 const vect<value_type,2>& aVelAAcc = vect<value_type,2>(),
		 const vect<value_type,2>& aAVelAcc = vect<value_type,2>(),
		 const value_type& aAVelAAcc = value_type()) : 
		 Parent(aParent),
		 vel_vel(aVelVel),
		 vel_avel(aVelAVel),
		 avel_vel(aAVelVel),
		 avel_avel(aAVelAVel),
		 vel_acc(aVelAcc),
		 vel_aacc(aVelAAcc),
		 avel_acc(aAVelAcc),
		 avel_aacc(aAVelAAcc) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(vel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(vel_aacc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_2D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(vel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(vel_aacc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000026,1,"jacobian_2D_2D",shared_object)
};

template <typename T> 
class jacobian_2D_3D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_2D_3D<T> self;
  
  boost::weak_ptr< frame_3D<value_type> > Parent;
  vect<vect<value_type,3>,2> vel_vel;
  vect<vect<value_type,3>,2> vel_avel;
  vect<value_type,3> avel_vel;
  vect<value_type,3> avel_avel;
  vect<vect<value_type,3>,2> vel_acc;
  vect<vect<value_type,3>,2> vel_aacc;
  vect<value_type,3> avel_acc;
  vect<value_type,3> avel_aacc;
  
  jacobian_2D_3D(boost::weak_ptr< frame_3D<value_type> > aParent = boost::weak_ptr< frame_3D<value_type> >(),
                 const vect<vect<value_type,3>,2>& aVelVel = vect<vect<value_type,3>,2>(),
		 const vect<vect<value_type,3>,2>& aVelAVel = vect<vect<value_type,3>,2>(),
		 const vect<value_type,3>& aAVelVel = vect<value_type,3>(),
		 const vect<value_type,3>& aAVelAVel = vect<value_type,3>(), 
		 const vect<vect<value_type,3>,2>& aVelAcc = vect<vect<value_type,3>,2>(),
		 const vect<vect<value_type,3>,2>& aVelAAcc = vect<vect<value_type,3>,2>(),
		 const vect<value_type,3>& aAVelAcc = vect<value_type,3>(),
		 const vect<value_type,3>& aAVelAAcc = vect<value_type,3>()) : 
		 Parent(aParent),
		 vel_vel(aVelVel),
		 vel_avel(aVelAVel),
		 avel_vel(aAVelVel),
		 avel_avel(aAVelAVel),
		 vel_acc(aVelAcc),
		 vel_aacc(aVelAAcc),
		 avel_acc(aAVelAcc),
		 avel_aacc(aAVelAAcc) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(vel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(vel_aacc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_3D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(vel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(vel_aacc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000027,1,"jacobian_2D_3D",shared_object)
};





template <typename T> 
class jacobian_3D_gen : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_3D_gen<T> self;
  
  vect<value_type,3> vel_qd;
  vect<value_type,3> avel_qd;
  vect<value_type,3> vel_qdd;
  vect<value_type,3> avel_qdd;
  
  jacobian_3D_gen(const vect<value_type,3>& aVelQd = vect<value_type,3>(),
		  const vect<value_type,3>& aAVelQd = vect<value_type,3>(), 
		  const vect<value_type,3>& aVelQdd = vect<value_type,3>(),
		  const vect<value_type,3>& aAVelQdd = vect<value_type,3>()) : 
		  vel_qd(aVelQd),
		  avel_qd(aAVelQd),
		  vel_qdd(aVelQdd),
		  avel_qdd(aAVelQdd) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    A & RK_SERIAL_SAVE_WITH_NAME(vel_qd)
      & RK_SERIAL_SAVE_WITH_NAME(avel_qd)
      & RK_SERIAL_SAVE_WITH_NAME(vel_qdd)
      & RK_SERIAL_SAVE_WITH_NAME(avel_qdd);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    A & RK_SERIAL_LOAD_WITH_NAME(vel_qd)
      & RK_SERIAL_LOAD_WITH_NAME(avel_qd)
      & RK_SERIAL_LOAD_WITH_NAME(vel_qdd)
      & RK_SERIAL_LOAD_WITH_NAME(avel_qdd);
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000028,1,"jacobian_3D_gen",shared_object)
};

template <typename T> 
class jacobian_3D_2D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_3D_2D<T> self;
  
  boost::weak_ptr< frame_2D<value_type> > Parent;
  vect<vect<value_type,2>,3> vel_vel;
  vect<value_type,3> vel_avel;
  vect<vect<value_type,2>,3> avel_vel;
  vect<value_type,3> avel_avel;
  vect<vect<value_type,2>,3> vel_acc;
  vect<value_type,3> vel_aacc;
  vect<vect<value_type,2>,3> avel_acc;
  vect<value_type,3> avel_aacc;
  
  jacobian_3D_2D(boost::weak_ptr< frame_2D<value_type> > aParent = boost::weak_ptr< frame_2D<value_type> >(),
                 const vect<vect<value_type,2>,3>& aVelVel = vect<vect<value_type,2>,3>(),
		 const vect<value_type,3>& aVelAVel = vect<value_type,3>(),
		 const vect<vect<value_type,2>,3>& aAVelVel = vect<vect<value_type,2>,3>(),
		 const vect<value_type,3>& aAVelAVel = vect<value_type,3>(), 
		 const vect<vect<value_type,2>,3>& aVelAcc = vect<vect<value_type,2>,3>(),
		 const vect<value_type,3>& aVelAAcc = vect<value_type,3>(),
		 const vect<vect<value_type,2>,3>& aAVelAcc = vect<vect<value_type,2>,3>(),
		 const vect<value_type,3>& aAVelAAcc = vect<value_type,3>()) : 
		 Parent(aParent),
		 vel_vel(aVelVel),
		 vel_avel(aVelAVel),
		 avel_vel(aAVelVel),
		 avel_avel(aAVelAVel),
		 vel_acc(aVelAcc),
		 vel_aacc(aVelAAcc),
		 avel_acc(aAVelAcc),
		 avel_aacc(aAVelAAcc) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(vel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(vel_aacc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_2D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(vel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(vel_aacc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x00000029,1,"jacobian_3D_2D",shared_object)
};

template <typename T> 
class jacobian_3D_3D : public shared_object {
public: 
  typedef T value_type;
  typedef jacobian_3D_3D<T> self;
  
  boost::weak_ptr< frame_3D<value_type> > Parent;
  vect<vect<value_type,3>,3> vel_vel;
  vect<vect<value_type,3>,3> vel_avel;
  vect<vect<value_type,3>,3> avel_vel;
  vect<vect<value_type,3>,3> avel_avel;
  vect<vect<value_type,3>,3> vel_acc;
  vect<vect<value_type,3>,3> vel_aacc;
  vect<vect<value_type,3>,3> avel_acc;
  vect<vect<value_type,3>,3> avel_aacc;
  
  jacobian_3D_3D(boost::weak_ptr< frame_3D<value_type> > aParent = boost::weak_ptr< frame_3D<value_type> >(),
                 const vect<vect<value_type,3>,3>& aVelVel = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aVelAVel = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aAVelVel = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aAVelAVel = vect<vect<value_type,3>,3>(), 
		 const vect<vect<value_type,3>,3>& aVelAcc = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aVelAAcc = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aAVelAcc = vect<vect<value_type,3>,3>(),
		 const vect<vect<value_type,3>,3>& aAVelAAcc = vect<vect<value_type,3>,3>()) : 
		 Parent(aParent),
		 vel_vel(aVelVel),
		 vel_avel(aVelAVel),
		 avel_vel(aAVelVel),
		 avel_avel(aAVelAVel),
		 vel_acc(aVelAcc),
		 vel_aacc(aVelAAcc),
		 avel_acc(aAVelAcc),
		 avel_aacc(aAVelAAcc) { };
  
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

  virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
    if(Parent.expired())
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",boost::shared_ptr<serialization::serializable>());
    else
      A & RK_SERIAL_SAVE_WITH_ALIAS("Parent",Parent.lock());
    A & RK_SERIAL_SAVE_WITH_NAME(vel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_vel)
      & RK_SERIAL_SAVE_WITH_NAME(avel_avel)
      & RK_SERIAL_SAVE_WITH_NAME(vel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(vel_aacc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_acc)
      & RK_SERIAL_SAVE_WITH_NAME(avel_aacc);
  };
  virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
    boost::shared_ptr< frame_3D<value_type> > tmp;
    A & RK_SERIAL_LOAD_WITH_ALIAS("Parent",tmp)
      & RK_SERIAL_LOAD_WITH_NAME(vel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_vel)
      & RK_SERIAL_LOAD_WITH_NAME(avel_avel)
      & RK_SERIAL_LOAD_WITH_NAME(vel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(vel_aacc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_acc)
      & RK_SERIAL_LOAD_WITH_NAME(avel_aacc);
    Parent = tmp;
  };

  RK_RTTI_MAKE_CONCRETE_1BASE(self,0x0000002A,1,"jacobian_3D_3D",shared_object)
};





  
};

#endif




