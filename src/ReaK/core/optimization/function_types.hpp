
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

#ifndef FUNCTION_TYPES_HPP
#define FUNCTION_TYPES_HPP

namespace ReaK {

namespace optim {


class NoFunctionAvailable { };

template <class Input, class Output, class HostClass>
class MemberFunctionPtr {
  public:
    typedef Output ( HostClass::*FuncPtr ) ( const Input& );
  private:
    HostClass* Obj;
    FuncPtr Func;
  public:
    MemberFunctionPtr ( HostClass* aObj, FuncPtr aFunc ) : Obj ( aObj ), Func ( aFunc ) { };
    Output operator() ( const Input& aInput ) {
      if ( Obj )
        return ( Obj->*Func ) ( aInput );
      else
        return Output();
    };
};

template <class Input, class Output, class HostClass>
inline MemberFunctionPtr<Input,Output,HostClass> make_function_ptr(HostClass* aObj, Output (HostClass::*aFunc) (const Input&) ) {
  return MemberFunctionPtr<Input,Output,HostClass>(aObj,aFunc);
};






/**
 * This class is the function-object interface for a cost function computation.
 */
template <class T>
class cost_function_1D
{
  public:
    /**
     * This function computes the cost for a given parameter.
     *
     * \param aParameter scalar parameter to the cost function
     * \return the scalar value of the cost function.
     */
    virtual T computeCost ( const T& aParameter ) const = 0;

    /** Default destructor. */
    virtual ~cost_function_1D() { };

};

/**
 * This structure holds a pair of value and cost values for a given optimization estimate, i.e. (x, f(x)).
 */
template <class T>
struct value_cost_pair
{
  T value, cost;
  value_cost_pair ( T aValue, T aCost ) : value ( aValue ), cost ( aCost ) { };
};

/**
 * This function template creates a value cost pair with requiring the explicit template instantiation.
 */
template <class T>
inline value_cost_pair<T> make_value_cost_pair ( T value, T cost ) {
  return value_cost_pair<T> ( value,cost );
};



};


};

#endif







