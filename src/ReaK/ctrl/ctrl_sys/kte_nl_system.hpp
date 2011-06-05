
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

#ifndef KTE_NL_SYSTEM_HPP
#define KTE_NL_SYSTEM_HPP



#include "base/named_object.hpp"

#include "math/mat_alg.hpp"
#include "math/vect_alg.hpp"

#include "math/kinetostatics.hpp"

#include "mbd_kte/kte_map_chain.hpp"
#include "mbd_kte/kte_system_input.hpp"
#include "mbd_kte/kte_system_output.hpp"
#include "mbd_kte/mass_matrix_calculator.hpp"

#include "math/mat_cholesky.hpp"

namespace ReaK {

namespace ctrl {


class kte_nl_system : public named_object {
  public:
    
    std::vector< boost::shared_ptr< gen_coord<double> > > dofs_gen;
    std::vector< boost::shared_ptr< frame_2D<double> > >  dofs_2D;
    std::vector< boost::shared_ptr< frame_3D<double> > >  dofs_3D;
    
    std::vector< boost::shared_ptr< kte::system_input > > inputs;
    std::vector< boost::shared_ptr< kte::system_output > > outputs;
    
    boost::shared_ptr< kte::kte_map_chain > chain;
    boost::shared_ptr< kte::mass_matrix_calc > mass_calc;
    
    typedef kte_nl_system self;
    typedef double value_type;
    typedef std::size_t size_type;
    
    typedef vect_n<double> point_type;
    typedef vect_n<double> point_difference_type;
    typedef vect_n<double> point_derivative_type;
    typedef self topology;
  
    typedef double time_type;
    typedef double time_difference_type;
  
    typedef vect_n<double> input_type;
    typedef vect_n<double> output_type;
    
    BOOST_STATIC_CONSTANT(std::size_t, dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, input_dimensions = 0);
    BOOST_STATIC_CONSTANT(std::size_t, output_dimensions = 0);
        
    kte_nl_system(const std::string& aName = "") {
      setName(aName);
    };

    kte_nl_system(const self& rhs) : named_object(), 
                                     dofs_gen(rhs.dofs_gen), 
                                     dofs_2D(rhs.dofs_2D), 
                                     dofs_3D(rhs.dofs_3D), 
                                     inputs(rhs.inputs), 
                                     outputs(rhs.outputs),
                                     chain(rhs.chain), 
                                     mass_calc(rhs.mass_calc) {
      setName(rhs.getName());
    };
      
    friend void swap(self& lhs, self& rhs) throw() {
      using std::swap;
      swap(lhs.dofs_gen,rhs.dofs_gen);
      swap(lhs.dofs_2D, rhs.dofs_2D);
      swap(lhs.dofs_3D, rhs.dofs_3D);
      swap(lhs.inputs, rhs.inputs);
      swap(lhs.chain, rhs.chain);
      swap(lhs.mass_calc, rhs.mass_calc);
    };

    self& operator =(self rhs) {
      swap(*this,rhs);
      return *this;
    };
    
    size_type get_state_count() const { 
      return 2 * dofs_gen.size() + 7 * dofs_2D.size() + 13 * dofs_3D.size(); 
    };
    size_type get_input_count() const { 
      size_type sum = 0;
      for(unsigned int i = 0; i < inputs.size(); ++i)
	sum += inputs[i]->getInputCount();
      return sum; 
    };
    size_type get_output_count() const {
      size_type sum = 0;
      for(unsigned int i = 0; i < outputs.size(); ++i)
	sum += outputs[i]->getOutputCount();
      return sum;
    };

    void apply_states_and_inputs(const point_type& p, const input_type& u) const {
      if(p.size() != get_state_count()) {
	std::cout << "Size of obtained state vector is: " << p.size() << std::endl;
	std::cout << "But expected the size of the state vector to be: " << get_state_count() << std::endl;
	throw std::range_error("State vector dimension mismatch!");
      };
      if(u.size() != get_input_count())
	throw std::range_error("Input vector dimension mismatch!");
      size_type i = 0;
      for(size_type j = 0; j < dofs_gen.size(); ++j) {
	dofs_gen[j]->q = p[i++]; 
	dofs_gen[j]->q_dot = p[i++]; 
	dofs_gen[j]->q_ddot = 0.0;
      };
      for(size_type j = 0; j < dofs_2D.size(); ++j) {
	dofs_2D[j]->Position[0] = p[i++];
	dofs_2D[j]->Position[1] = p[i++];
	vect<double,2> tmp; tmp[0] = p[i++]; tmp[1] = p[i++];
	dofs_2D[j]->Rotation = frame_2D<double>::rotation_type(tmp);
	dofs_2D[j]->Velocity[0] = p[i++];
	dofs_2D[j]->Velocity[1] = p[i++];
	dofs_2D[j]->AngVelocity = p[i++];
	dofs_2D[j]->Acceleration = vect<double,2>();
	dofs_2D[j]->AngAcceleration = 0.0;
      };
      for(size_type j = 0; j < dofs_3D.size(); ++j) {
	dofs_3D[j]->Position[0] = p[i++];
	dofs_3D[j]->Position[1] = p[i++];
	dofs_3D[j]->Position[2] = p[i++];
	vect<double,4> tmp; tmp[0] = p[i++]; tmp[1] = p[i++]; tmp[2] = p[i++]; tmp[3] = p[i++];
	dofs_3D[j]->Quat = quaternion<double>(tmp);
	dofs_3D[j]->Velocity[0] = p[i++];
	dofs_3D[j]->Velocity[1] = p[i++];
	dofs_3D[j]->Velocity[2] = p[i++];
	dofs_3D[j]->AngVelocity[0] = p[i++];
	dofs_3D[j]->AngVelocity[1] = p[i++];
	dofs_3D[j]->AngVelocity[2] = p[i++];
	dofs_3D[j]->Acceleration = vect<double,3>();
	dofs_3D[j]->AngAcceleration = vect<double,3>();
      };
      
      i = 0;
      for(size_type j = 0; j < inputs.size(); ++j)
	for(size_type k = 0; k < inputs[j]->getInputCount(); k++)
	  inputs[j]->getInput(k) = u[i++];
    };
    
    point_derivative_type get_state_derivative(const point_type& p, const input_type& u, const time_type& t = 0) const {
      apply_states_and_inputs(p,u);
      
      point_derivative_type pd(get_state_count());
      
      if(chain && mass_calc) {
        chain->doMotion();
        chain->clearForce();
        chain->doForce();
      
	mat<double, mat_structure::symmetric> M(dofs_gen.size() + 3 * dofs_2D.size() + 6 * dofs_3D.size());
	vect_n<double> f(dofs_gen.size() + 3 * dofs_2D.size() + 6 * dofs_3D.size());
	
	for(size_type i = 0; i < dofs_gen.size(); ++i)
	  f[i] = dofs_gen[i]->f;
	
	size_type base_i = dofs_gen.size();
	for(size_type i = 0; i < dofs_2D.size(); ++i) {
	  f[base_i + 3*i] = dofs_2D[i]->Force[0];
	  f[base_i + 3*i + 1] = dofs_2D[i]->Force[1];
	  f[base_i + 3*i + 2] = dofs_2D[i]->Torque;
	};
	
	base_i = dofs_gen.size() + 3 * dofs_2D.size();
	for(size_type i = 0; i < dofs_3D.size(); ++i) {
	  f[base_i + 6*i] = dofs_3D[i]->Force[0];
	  f[base_i + 6*i + 1] = dofs_3D[i]->Force[1];
	  f[base_i + 6*i + 2] = dofs_3D[i]->Force[2];
	  f[base_i + 6*i + 3] = dofs_3D[i]->Torque[0];
	  f[base_i + 6*i + 4] = dofs_3D[i]->Torque[1];
	  f[base_i + 6*i + 5] = dofs_3D[i]->Torque[2];
	};
	
	mass_calc->getMassMatrix(M);
	mat_vect_adaptor< vect_n<double>, mat_alignment::column_major > f_adapt(f);
	linsolve_Cholesky(M, f_adapt);
	
	for(size_type i = 0; i < dofs_gen.size(); ++i) {
	  pd[2*i] = dofs_gen[i]->q_dot;
	  pd[2*i + 1] = f[i];
	};
	
	base_i = dofs_gen.size();
	for(size_type i = 0; i < dofs_2D.size(); ++i) {
	  pd[2*base_i + 7*i] = dofs_2D[i]->Velocity[0];     
	  pd[2*base_i + 7*i + 1] = dofs_2D[i]->Velocity[1]; 
	  pd[2*base_i + 7*i + 2] = -p[2*base_i + 7*i + 3] * dofs_2D[i]->AngVelocity; 
	  pd[2*base_i + 7*i + 3] =  p[2*base_i + 7*i + 2] * dofs_2D[i]->AngVelocity; 
	  pd[2*base_i + 7*i + 4] = f[base_i + 3*i];       
	  pd[2*base_i + 7*i + 5] = f[base_i + 3*i + 1];   
	  pd[2*base_i + 7*i + 6] = f[base_i + 3*i + 2];   
	};
	
	base_i = dofs_gen.size() + 3 * dofs_2D.size();
	for(size_type i = 0; i < dofs_3D.size(); ++i) {
	  pd[2*base_i + dofs_2D.size() + 13*i] = dofs_3D[i]->Velocity[0];
	  pd[2*base_i + dofs_2D.size() + 13*i + 1] = dofs_3D[i]->Velocity[1];
	  pd[2*base_i + dofs_2D.size() + 13*i + 2] = dofs_3D[i]->Velocity[2];
	  dofs_3D[i]->UpdateQuatDot();
	  pd[2*base_i + dofs_2D.size() + 13*i + 3] = dofs_3D[i]->QuatDot[0];
	  pd[2*base_i + dofs_2D.size() + 13*i + 4] = dofs_3D[i]->QuatDot[1];
	  pd[2*base_i + dofs_2D.size() + 13*i + 5] = dofs_3D[i]->QuatDot[2];
	  pd[2*base_i + dofs_2D.size() + 13*i + 6] = dofs_3D[i]->QuatDot[3];
	  pd[2*base_i + dofs_2D.size() + 13*i + 7] = f[base_i + 6*i];
	  pd[2*base_i + dofs_2D.size() + 13*i + 8] = f[base_i + 6*i + 1];
	  pd[2*base_i + dofs_2D.size() + 13*i + 9] = f[base_i + 6*i + 2];
	  pd[2*base_i + dofs_2D.size() + 13*i + 10] = f[base_i + 6*i + 3];
	  pd[2*base_i + dofs_2D.size() + 13*i + 11] = f[base_i + 6*i + 4];
	  pd[2*base_i + dofs_2D.size() + 13*i + 12] = f[base_i + 6*i + 5];
	};
      } else {
	for(size_type i = 0; i < dofs_gen.size(); ++i) {
	  pd[2*i] = dofs_gen[i]->q_dot;
	  pd[2*i + 1] = 0.0;
	};
	
	size_type base_i = dofs_gen.size();
	for(size_type i = 0; i < dofs_2D.size(); ++i) {
	  pd[2*base_i + 7*i] = dofs_2D[i]->Velocity[0];     
	  pd[2*base_i + 7*i + 1] = dofs_2D[i]->Velocity[1]; 
	  pd[2*base_i + 7*i + 2] = -p[2*base_i + 7*i + 3] * dofs_2D[i]->AngVelocity; 
	  pd[2*base_i + 7*i + 3] =  p[2*base_i + 7*i + 2] * dofs_2D[i]->AngVelocity; 
	  pd[2*base_i + 7*i + 4] = 0.0;       
	  pd[2*base_i + 7*i + 5] = 0.0;   
	  pd[2*base_i + 7*i + 6] = 0.0;   
	};
	
	base_i = dofs_gen.size() + 3 * dofs_2D.size();
	for(size_type i = 0; i < dofs_3D.size(); ++i) {
	  pd[2*base_i + dofs_2D.size() + 13*i] = dofs_3D[i]->Velocity[0];
	  pd[2*base_i + dofs_2D.size() + 13*i + 1] = dofs_3D[i]->Velocity[1];
	  pd[2*base_i + dofs_2D.size() + 13*i + 2] = dofs_3D[i]->Velocity[2];
	  dofs_3D[i]->UpdateQuatDot();
	  pd[2*base_i + dofs_2D.size() + 13*i + 3] = dofs_3D[i]->QuatDot[0];
	  pd[2*base_i + dofs_2D.size() + 13*i + 4] = dofs_3D[i]->QuatDot[1];
	  pd[2*base_i + dofs_2D.size() + 13*i + 5] = dofs_3D[i]->QuatDot[2];
	  pd[2*base_i + dofs_2D.size() + 13*i + 6] = dofs_3D[i]->QuatDot[3];
	  pd[2*base_i + dofs_2D.size() + 13*i + 7] = 0.0;
	  pd[2*base_i + dofs_2D.size() + 13*i + 8] = 0.0;
	  pd[2*base_i + dofs_2D.size() + 13*i + 9] = 0.0;
	  pd[2*base_i + dofs_2D.size() + 13*i + 10] = 0.0;
	  pd[2*base_i + dofs_2D.size() + 13*i + 11] = 0.0;
	  pd[2*base_i + dofs_2D.size() + 13*i + 12] = 0.0;
	};
      };
      
      return pd;
    };
    
    output_type get_output(const point_type& p, const input_type& u, const time_type& t = 0) const {
      apply_states_and_inputs(p,u);
      
      chain->doMotion();
      chain->clearForce();
      chain->doForce();
      
      output_type y(get_output_count());
      
      size_type i = 0;
      for(size_type j = 0; j < outputs.size(); ++j)
	for(size_type k = 0; k < outputs[j]->getOutputCount(); k++)
	  y[i++] = outputs[j]->getOutput(k);
	
      return y;
    };
    
    point_type adjust(const point_type& p, const point_difference_type& dp) const {
      return p + dp;
    };
    
    
    
    virtual void RK_CALL save(ReaK::serialization::oarchive& aA, unsigned int) const {
      ReaK::named_object::save(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_SAVE_WITH_NAME(dofs_gen)
         & RK_SERIAL_SAVE_WITH_NAME(dofs_2D)
         & RK_SERIAL_SAVE_WITH_NAME(dofs_3D)
         & RK_SERIAL_SAVE_WITH_NAME(inputs)
         & RK_SERIAL_SAVE_WITH_NAME(outputs)
         & RK_SERIAL_SAVE_WITH_NAME(chain)
         & RK_SERIAL_SAVE_WITH_NAME(mass_calc);
    };
    virtual void RK_CALL load(ReaK::serialization::iarchive& aA, unsigned int) {
      ReaK::named_object::load(aA,ReaK::named_object::getStaticObjectType()->TypeVersion());
      aA & RK_SERIAL_LOAD_WITH_NAME(dofs_gen)
         & RK_SERIAL_LOAD_WITH_NAME(dofs_2D)
         & RK_SERIAL_LOAD_WITH_NAME(dofs_3D)
         & RK_SERIAL_LOAD_WITH_NAME(inputs)
         & RK_SERIAL_LOAD_WITH_NAME(outputs)
         & RK_SERIAL_LOAD_WITH_NAME(chain)
         & RK_SERIAL_LOAD_WITH_NAME(mass_calc);
    };
    
    RK_RTTI_MAKE_CONCRETE_1BASE(self,0xC2300002,1,"kte_nl_system",named_object)
  
};




  
};

};

#endif














