/**
 * \file backtrace_exception.hpp
 * 
 * This library declares the "backtrace_except<E>" class template that can be used to bundle an 
 * exception with file-line back-trace of where it occurred and in what contexts it was caught and 
 * re-thrown.
 * 
 * \author Mikael Persson (mikael.s.persson@gmail.com)
 * \date June 2013
 */

/*
 *    Copyright 2013 Sven Mikael Persson
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

#ifndef REAK_BACKTRACE_EXCEPTION_HPP
#define REAK_BACKTRACE_EXCEPTION_HPP

#include "defs.hpp"


#include <string>
#include <sstream>
#include <exception>


namespace ReaK {


/**
 * This class template
 */
template <typename BaseException>
class backtrace_except : public BaseException {
  private:
    std::string backtrace_str;
    
    void create_backtrace_string(const char* aFileName, const char* aFuncName, int aLineNum, const char* aMessage) {
      std::stringstream ss;
      ss << "From '" << aFileName << "' in function '" << aFuncName << "' at line " << aLineNum << ":\n" << aMessage;
      backtrace_str = ss.str();
    };
    
    void create_compound_backtrace_string(const char* aFileName, const char* aFuncName, int aLineNum, const char* aNewMessage, const char* aOldMessage) {
      std::stringstream ss;
      ss << "From '" << aFileName << "' in function '" << aFuncName << "' at line " << aLineNum << ":\n" 
         << aNewMessage
         << "\nCaused by the following underlying issue:\n" 
         << aOldMessage;
      backtrace_str = ss.str();
    };
    
  public:
    
    enum transformed_exception_t { transformed_exception };
    
#ifdef RK_ENABLE_CXX11_FEATURES
    template <typename... Args>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, Args&&... args) :
                     BaseException(std::forward<Args>(args)...) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    
    template <typename... Args>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, Args&&... args) :
                     BaseException(std::forward<Args>(args)...) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
#else
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum) :
                     BaseException() {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    template <typename Arg1>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1) :
                     BaseException(arg1) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    template <typename Arg1, typename Arg2>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2) :
                     BaseException(arg1, arg2) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    template <typename Arg1, typename Arg2, typename Arg3>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3) :
                     BaseException(arg1, arg2, arg3) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4) :
                     BaseException(arg1, arg2, arg3, arg4) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
    backtrace_except(const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4, const Arg5& arg5) :
                     BaseException(arg1, arg2, arg3, arg4, arg5) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what());
    };
    
    
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum) :
                     BaseException() {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
    template <typename Arg1>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1) :
                     BaseException(arg1) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
    template <typename Arg1, typename Arg2>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2) :
                     BaseException(arg1, arg2) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
    template <typename Arg1, typename Arg2, typename Arg3>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3) :
                     BaseException(arg1, arg2, arg3) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4) :
                     BaseException(arg1, arg2, arg3, arg4) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
    template <typename Arg1, typename Arg2, typename Arg3, typename Arg4, typename Arg5>
    backtrace_except(transformed_exception_t, const std::exception& e, 
                     const char* aFileName, const char* aFuncName, int aLineNum, 
                     const Arg1& arg1, const Arg2& arg2, const Arg3& arg3, const Arg4& arg4, const Arg5& arg5) :
                     BaseException(arg1, arg2, arg3, arg4, arg5) {
      create_compound_backtrace_string(aFileName, aFuncName, aLineNum, BaseException::what(), e.what());
    };
#endif
    
    backtrace_except(const std::exception& e, const char* aFileName, const char* aFuncName, int aLineNum) :
                     BaseException(static_cast< const BaseException& >(e)) {
      create_backtrace_string(aFileName, aFuncName, aLineNum, e.what());
    };
    
    
    virtual ~backtrace_except() 
#ifdef RK_ENABLE_CXX11_FEATURES
      noexcept
#else
      throw()
#endif
      { };
    
    virtual const char* what() const 
#ifdef RK_ENABLE_CXX11_FEATURES
      noexcept
#else
      throw()
#endif
    {
      return backtrace_str.c_str();
    };
    
};


template <typename BaseException>
backtrace_except<BaseException> make_backtrace_exception(const BaseException& e, const char* aFileName, const char* aFuncName, int aLineNum) {
  return backtrace_except<BaseException>(e, aFileName, aFuncName, aLineNum);
};


};


#ifdef RK_ENABLE_EXCEPTION_BACKTRACING


#define RK_THROW_WITH_BACKTRACE_0(EXCEPTION) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__)
#define RK_THROW_WITH_BACKTRACE_1(EXCEPTION, ARG1) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1)
#define RK_THROW_WITH_BACKTRACE_2(EXCEPTION, ARG1, ARG2) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2)
#define RK_THROW_WITH_BACKTRACE_3(EXCEPTION, ARG1, ARG2, ARG3) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3)
#define RK_THROW_WITH_BACKTRACE_4(EXCEPTION, ARG1, ARG2, ARG3, ARG4) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3, ARG4)
#define RK_THROW_WITH_BACKTRACE_5(EXCEPTION, ARG1, ARG2, ARG3, ARG4, ARG5) throw ReaK::backtrace_except< EXCEPTION >(__FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3, ARG4, ARG5)

#define RK_RETHROW_WITH_BACKTRACE(EXCEPT_NAME) throw ReaK::make_backtrace_exception(EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__)

#define RK_THROW_NEW_WITH_BACKTRACE_0(EXCEPT_NAME, EXCEPTION) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__)
#define RK_THROW_NEW_WITH_BACKTRACE_1(EXCEPT_NAME, EXCEPTION, ARG1) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1)
#define RK_THROW_NEW_WITH_BACKTRACE_2(EXCEPT_NAME, EXCEPTION, ARG1, ARG2) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2)
#define RK_THROW_NEW_WITH_BACKTRACE_3(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3)
#define RK_THROW_NEW_WITH_BACKTRACE_4(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3, ARG4) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3, ARG4)
#define RK_THROW_NEW_WITH_BACKTRACE_5(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3, ARG4, ARG5) throw ReaK::backtrace_except< EXCEPTION >(ReaK::backtrace_except< EXCEPTION >::transformed_exception, EXCEPT_NAME, __FILE__, RK_CURRENT_FUNCTION, __LINE__, ARG1, ARG2, ARG3, ARG4, ARG5)


#else


#define RK_THROW_WITH_BACKTRACE_0(EXCEPTION) throw EXCEPTION()
#define RK_THROW_WITH_BACKTRACE_1(EXCEPTION, ARG1) throw EXCEPTION(ARG1)
#define RK_THROW_WITH_BACKTRACE_2(EXCEPTION, ARG1, ARG2) throw EXCEPTION(ARG1, ARG2)
#define RK_THROW_WITH_BACKTRACE_3(EXCEPTION, ARG1, ARG2, ARG3) throw EXCEPTION(ARG1, ARG2, ARG3)
#define RK_THROW_WITH_BACKTRACE_4(EXCEPTION, ARG1, ARG2, ARG3, ARG4) throw EXCEPTION(ARG1, ARG2, ARG3, ARG4)
#define RK_THROW_WITH_BACKTRACE_5(EXCEPTION, ARG1, ARG2, ARG3, ARG4, ARG5) throw EXCEPTION(ARG1, ARG2, ARG3, ARG4, ARG5)

#define RK_RETHROW_WITH_BACKTRACE(EXCEPT_NAME) throw EXCEPT_NAME

#define RK_THROW_NEW_WITH_BACKTRACE_0(EXCEPT_NAME, EXCEPTION) throw EXCEPTION()
#define RK_THROW_NEW_WITH_BACKTRACE_1(EXCEPT_NAME, EXCEPTION, ARG1) throw EXCEPTION(ARG1)
#define RK_THROW_NEW_WITH_BACKTRACE_2(EXCEPT_NAME, EXCEPTION, ARG1, ARG2) throw EXCEPTION(ARG1, ARG2)
#define RK_THROW_NEW_WITH_BACKTRACE_3(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3) throw EXCEPTION(ARG1, ARG2, ARG3)
#define RK_THROW_NEW_WITH_BACKTRACE_4(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3, ARG4) throw EXCEPTION(ARG1, ARG2, ARG3, ARG4)
#define RK_THROW_NEW_WITH_BACKTRACE_5(EXCEPT_NAME, EXCEPTION, ARG1, ARG2, ARG3, ARG4, ARG5) throw EXCEPTION(ARG1, ARG2, ARG3, ARG4, ARG5)


#endif


#endif


