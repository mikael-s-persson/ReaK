/**
 * \file expected.hpp
 * 
 * This library declares the "expected<T>" class template which is based on Andrei Alexandrescu's 
 * presentation of the topic at the "C++ and Beyond 2012" conference. This particular implementation
 * is the work of Mikael Persson, but the entire idea is public domain.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date December 2012
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

#ifndef REAK_EXPECTED_HPP
#define REAK_EXPECTED_HPP

#include "defs.hpp"

#ifdef RK_ENABLE_CXX11_FEATURES

#include <exception>
#include <stdexcept>

namespace ReaK {

/**
 * This class template represents a value-expectation. This can be used to store 
 * either a value of type T or an exception associated to a failure to produce the 
 * value in question (e.g., a return value), i.e., this class allows one to fail to
 * meet the expectation of a value.
 * This class template is based on Andrei Alexandrescu's presentation of "Expected<T>"
 * at the "C++ and Beyond 2012" conference. This particular implementation
 * is the work of Mikael Persson, but the entire idea is public domain.
 * This class template can be used as follows:
 * \code
 * expected<double> divide(int n, int d) {
 *   if( d == 0 )
 *     return expected<double>::failed(std::invalid_argument("The denominator of the fraction is zero."));
 *   return expected<double>(double(n) / d);
 * };
 * \endcode
 */
template <typename T>
class expected {
  private:
    union { 
      T value;
      std::exception_ptr except;
    };
    bool hasValue;
    
    struct construct_failure { };
    
    expected(construct_failure, std::exception_ptr&& aExcept) : except(std::move(aExcept)), hasValue(false) { };
    
  public:
    template <typename U>
    friend class expected;
    
    /**
     * Constructor from value (implicit conversion).
     * \param aValue The value that fulfills the value-expectation.
     */
    expected(const T& aValue = T()) : value(aValue), hasValue(true) { };
    
    /**
     * Constructor from moved value (implicit move-conversion).
     * \param aValue The value that fulfills the value-expectation (by rvalue-ref).
     */
    expected(T&& aValue) : value(std::move(aValue)), hasValue(true) { };
    
    /**
     * Copy-constructor.
     * \tparam U A value-type which is copyable to a value of type T.
     * \param rhs The right-hand-side of the copy operation.
     */
    template <typename U>
    expected(const expected<U>& rhs) : hasValue(rhs.hasValue) {
      if(hasValue)
        new(&value) T(rhs.value);
      else
        new(&except) std::exception_ptr(rhs.except);
    };
    
    /**
     * Move-constructor.
     * \tparam U A value-type which is moveable to a value of type T.
     * \param rhs The right-hand-side of the move operation.
     */
    template <typename U>
    expected(expected<U>&& rhs) : hasValue(rhs.hasValue) {
      if(hasValue)
        new(&value) T(std::move(rhs.value));
      else
        new(&except) std::exception_ptr(std::move(rhs.except));
    };
    
    /**
     * Standard swap function for the value-expectation objects.
     * \param lhs The left-hand-side of the swap.
     * \param rhs The right-hand-side of the swap.
     */
    friend void swap(expected<T>& lhs, expected<T>& rhs) {
      switch( (int(lhs.hasValue) << 1) + int(rhs.hasValue) ) {
        case 0:
          lhs.except.swap(rhs.except);
          break;
        case 1:
          std::exception_ptr tmp(std::move(lhs.except));
          new(&lhs.value) T(std::move(rhs.value));
          new(&rhs.except) std::exception_ptr(std::move(tmp));
          rhs.hasValue = false; lhs.hasValue = true;
          break;
        case 2:
          std::exception_ptr tmp(std::move(rhs.except));
          new(&rhs.value) T(std::move(lhs.value));
          new(&lhs.except) std::exception_ptr(std::move(tmp));
          lhs.hasValue = false; rhs.hasValue = true;
          break;
        case 3:
          using std::swap;
          swap(lhs.value, rhs.value);
          break;
      };
    };
    
    /**
     * Copy- and Move-assignment operator using the dual copy-and-swap / move-and-swap idiom.
     */
    expected<T>& operator=(expected<T> rhs) {
      swap(*this, rhs);
      return *this;
    };
    
    /**
     * Destructor.
     */
    ~expected() {
      using std::exception_ptr;
      if(hasValue) value.~T();
      else except.~exception_ptr();
    };
    
    /**
     * This function checks if this value-expectation object stores a valid object or not (failure).
     * \return True if the expectation was met with a valid value; false if it failed.
     */
    bool valid() const {
      return hasValue;
    };
    
    /**
     * This function returns a reference to the internally-stored value.
     * \note If valid, one should consider copying out or moving out the value instead of 
     *       constantly using it with this get() function.
     * \throw std::exception_ptr If this value-expectation object contains a failure, it will re-throw that failure as an exception.
     * \return A reference to the value stored in this value-expectation (if met).
     */
    T& get() {
      if(!hasValue)
        std::rethrow_exception(except);
      return value;
    };
    
    /**
     * This function returns a const-reference to the internally-stored value.
     * \note If valid, one should consider copying out or moving out the value instead of 
     *       constantly using it with this get() function.
     * \throw std::exception_ptr If this value-expectation object contains a failure, it will re-throw that failure as an exception.
     * \return A const-reference to the value stored in this value-expectation (if met).
     */
    const T& get() const {
      if(!hasValue)
        std::rethrow_exception(except);
      return value;
    };
    
    /**
     * This function checks if this value-expectation object contains a particular 
     * type of exception. This function is slow in the sense that it must re-throw 
     * the exception-pointer, then catch it again (because std::exception_ptr is opaque).
     */
    template <typename E>
    bool has_exception() const {
      try {
        if(!hasValue)
          std::rethrow_exception(except);
      } catch(const E& e) { RK_UNUSED(e);
        return true;
      } catch(...) { };
      return false;
    };
    
    
    /**
     * This static function is used to create a value-expectation object with a 
     * failure to meet the expectation and provide, instead, an exception object,
     * via a standard exception pointer.
     */
    static expected<T> failed(std::exception_ptr aExceptPtr) {
      return expected<T>(construct_failure(), std::move(aExceptPtr));
    };
    
    /**
     * This static function is used to create a value-expectation object with a 
     * failure to meet the expectation and provide, instead, an exception object,
     * via the current exception (the one obtained by std::current_exception()).
     */
    static expected<T> failed() {
      return expected<T>(construct_failure(), std::current_exception());
    };
    
    /**
     * This static function is used to create a value-expectation object with a 
     * failure to meet the expectation and provide, instead, an exception object.
     */
    template <typename E>
    static expected<T> failed(const E& aExcept) {
      if( typeid(aExcept) != typeid(E) )
        throw std::invalid_argument("Exception raised on failure of value-expectation cannot be carried through without slicing!");
      return failed(std::make_exception_ptr(aExcept));
    };
    
    /**
     * This static function takes a functor (function object, function pointer, 
     * lambda expression, etc.) and tries to execute it, catching any exception that
     * would arise from it and storing it as a failure on the value-expectation that 
     * it returns.
     * \tparam F Some callable type (functor) which can be called with no arguments and returns something convertible to type T.
     * \param fun The callable object that is to be tried.
     * \return A value-expectation object which contains either the result of calling fun() or the exception that the fun() call produced.
     */
    template <typename F>
    static expected<T> try_code(F fun) {
      try {
        return expected<T>(fun());
      } catch(...) {
        return failed();
      };
    };
    
};


};


#else


namespace ReaK {

/**
 * 
 */
template <typename T>
class expected {
  public:
    
    
};

};


#endif

#endif


