
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

#ifndef RIGID_BODY_STATES_HPP
#define RIGID_BODY_STATES_HPP

#include "base/shared_object.hpp"
#include "ctrl_sys/state_vector_concept.hpp"

#include "math/vect_alg.hpp"
#include "math/rotations.hpp"
#include "math/quat_alg.hpp"
#include "math/vect_index_iterator.hpp"

namespace ReaK {

namespace ctrl {
  

template <typename T>
class rigid_body2D_state : public shared_object {
  public:
    typedef rigid_body2D_state<T> self;
    typedef T value_type;
    typedef self state_type;
    typedef vect_n<value_type> state_difference_type;
    typedef typename vect_traits< state_difference_type >::size_type size_type;
    typedef typename vect_traits< state_difference_type >::difference_type difference_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef typename vect_traits< state_difference_type >::allocator_type allocator_type; //just in case it is cast to variable-length vector.
  
    typedef vect_index_iter<self> iterator;
    typedef vect_index_const_iter<self> const_iterator;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 7);
    
    vect< value_type, 2> position;
    rot_mat_2D<value_type> rotation;
    vect< value_type, 2> momentum;
    value_type angular_momentum;
    
    rigid_body2D_state(const vect<value_type,2>& aPosition = vect<value_type,2>(),
		       const rot_mat_2D<value_type>& aRotation = rot_mat_2D<value_type>(),
		       const vect<value_type,2>& aMomentum = vect<value_type,2>(),
		       const value_type& aAngMomentum = value_type(0.0)) :
		       position(aPosition),
		       rotation(aRotation),
		       momentum(aMomentum),
		       angular_momentum(aAngMomentum) { };
		       
    value_type operator[](size_type i) const {
      if(i < 2)
	return position[i];
      else if(i < 4)
	return rotation[i - 2];
      else if(i < 6)
	return momentum[i - 4];
      else
	return angular_momentum;
    };
    
    size_type size() const { return 7; };
    
    iterator begin() { return iterator(*this,0); };
    const_iterator begin() const { return const_iterator(*this,0); };
    iterator end() { return iterator(*this); };
    const_iterator end() const { return const_iterator(*this); };
    
    friend
    state_difference_type diff(const self& a, const self& b) {
      return state_difference_type(a.position[0] - b.position[0],
				   a.position[1] - b.position[1],
				   (invert(b.rotation) * a.rotation).getAngle(),
				   a.momentum[0] - b.momentum[0],
				   a.momentum[1] - b.momentum[1],
				   a.angular_momentum - b.angular_momentum);
    };
    
    friend
    self add(const self& a, const state_difference_type& da) {
      return self(a.position + vect<value_type,2>(da[0], da[1]),
		  a.rotation * rot_mat_2D<value_type>(da[2]),
		  a.momentum + vect<value_type,2>(da[3],da[4]),
		  a.angular_momentum + da[5]);
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(position)
        & RK_SERIAL_SAVE_WITH_NAME(rotation)
        & RK_SERIAL_SAVE_WITH_NAME(momentum)
        & RK_SERIAL_SAVE_WITH_NAME(angular_momentum);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(position)
        & RK_SERIAL_LOAD_WITH_NAME(rotation)
        & RK_SERIAL_LOAD_WITH_NAME(momentum)
        & RK_SERIAL_LOAD_WITH_NAME(angular_momentum);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2310010,1,"rigid_body2D_state",shared_object)
    
};


template <typename T>
struct is_readable_vector< rigid_body2D_state<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< rigid_body2D_state<T> > type;
};



template <typename T>
class rigid_body3D_state : public shared_object {
  public:
    typedef rigid_body3D_state<T> self;
    typedef T value_type;
    typedef self state_type;
    typedef vect_n<value_type> state_difference_type;
    typedef typename vect_traits< state_difference_type >::size_type size_type;
    typedef typename vect_traits< state_difference_type >::difference_type difference_type;
    
    typedef T& reference;
    typedef const T& const_reference;
    typedef T* pointer;
    typedef const T* const_pointer;
    typedef typename vect_traits< state_difference_type >::allocator_type allocator_type; //just in case it is cast to variable-length vector.
  
    typedef vect_index_iter<self> iterator;
    typedef vect_index_const_iter<self> const_iterator;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 13);
    
    vect< value_type, 3> position;
    quaternion<value_type> rotation;
    vect< value_type, 3> momentum;
    vect< value_type, 3> angular_momentum;
    
    rigid_body2D_state(const vect<value_type,3>& aPosition = vect<value_type,3>(),
		       const quaternion<value_type>& aRotation = quaternion<value_type>(),
		       const vect<value_type,3>& aMomentum = vect<value_type,3>(),
		       const vect<value_type,3>& aAngMomentum = vect<value_type,3>()) :
		       position(aPosition),
		       rotation(aRotation),
		       momentum(aMomentum),
		       angular_momentum(aAngMomentum) { };
		       
    value_type operator[](size_type i) const {
      if(i < 3)
	return position[i];
      else if(i < 7)
	return rotation[i - 3];
      else if(i < 10)
	return momentum[i - 7];
      else
	return angular_momentum[i - 10];
    };
    
    size_type size() const { return 7; };
    
    iterator begin() { return iterator(*this,0); };
    const_iterator begin() const { return const_iterator(*this,0); };
    iterator end() { return iterator(*this); };
    const_iterator end() const { return const_iterator(*this); };
    
    
    
    friend
    state_difference_type diff(const self& a, const self& b) {
      quaternion<value_type> q_diff = invert(b.rotation) * a.rotation;
      quat<value_type> ang_diff = log(quat<value_type>(q_diff[0],q_diff[1],q_diff[2],q_diff[3]));
      vect<value_type,3> ang_mom_diff = a.angular_momentum - invert(q_diff) * b.angular_momentum;
      return state_difference_type(a.position[0] - b.position[0],
				   a.position[1] - b.position[1],
				   a.position[2] - b.position[2],
				   2.0 * ang_diff[1],
				   2.0 * ang_diff[2],
				   2.0 * ang_diff[3],
				   a.momentum[0] - b.momentum[0],
				   a.momentum[1] - b.momentum[1],
				   a.momentum[2] - b.momentum[2],
				   ang_mom_diff[0],
				   ang_mom_diff[1],
				   ang_mom_diff[2]);
    };
    
    friend
    self add(const self& a, const state_difference_type& da) {
      quaternion<value_type> q_diff(exp(quat<value_type>(0.0, 0.5 * da[3], 0.5 * da[4], 0.5 * da[5])));
      return self(a.position + vect<value_type,3>(da[0], da[1], da[2]),
		  a.rotation * q_diff,
		  a.momentum + vect<value_type,3>(da[6],da[7],da[8]),
		  invert(q_diff) * a.angular_momentum + vect<value_type,3>(da[9],da[10],da[11]));
    };
    
    
/*******************************************************************************
                   ReaK's RTTI and Serialization interfaces
*******************************************************************************/

    virtual void RK_CALL save(ReaK::serialization::oarchive& A, unsigned int) const {
      shared_object::save(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(position)
        & RK_SERIAL_SAVE_WITH_NAME(rotation)
        & RK_SERIAL_SAVE_WITH_NAME(momentum)
        & RK_SERIAL_SAVE_WITH_NAME(angular_momentum);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& A, unsigned int) {
      shared_object::load(A,shared_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(position)
        & RK_SERIAL_LOAD_WITH_NAME(rotation)
        & RK_SERIAL_LOAD_WITH_NAME(momentum)
        & RK_SERIAL_LOAD_WITH_NAME(angular_momentum);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2310011,1,"rigid_body3D_state",shared_object)
    
};


template <typename T>
struct is_readable_vector< rigid_body3D_state<T> > {
  BOOST_STATIC_CONSTANT( bool, value = true );
  typedef is_readable_vector< rigid_body3D_state<T> > type;
};




};

};


#endif



















