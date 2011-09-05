/**
 * \file mat_esn_expressions.hpp
 * 
 * This library is a work in progress. As of now, it does nothing useful. The idea of this 
 * library is to use Einstein Summation Notation (ESN) to construct matrix expression templates.
 */

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

#ifndef REAK_MAT_ESN_EXPRESSIONS_HPP
#define REAK_MAT_ESN_EXPRESSIONS_HPP

#include "mat_traits.hpp"
#include "vect_traits.hpp"

namespace ReaK {


enum esn_i_tag { esn_i };
enum esn_j_tag { esn_j };
enum esn_k_tag { esn_k };
enum esn_l_tag { esn_l };
enum esn_m_tag { esn_m };
enum esn_n_tag { esn_n };
enum esn_o_tag { esn_o };
enum esn_p_tag { esn_p };
enum esn_q_tag { esn_q };
enum esn_r_tag { esn_r };
enum esn_s_tag { esn_s };
enum esn_t_tag { esn_t };
enum esn_u_tag { esn_u };


namespace detail {
  
  enum op_class {
    term,
    split,
    selector,
    permutor
  };
  
  template <typename T,typename U>
  struct op_addition {
    BOOST_STATIC_CONSTANT(unsigned int, arity = 2);
    BOOST_STATIC_CONSTANT(op_class, classification = split);
#ifndef RK_ENABLE_CXX0X_FEATURES
    static T evaluate(const T& t, const U& u) { return t + u; };
#else
    static auto evaluate(T&& t, U&& u) -> decltype(std::move(t) + std::move(u)) { return std::move(t) + std::move(u); };
    static auto evaluate(const T& t, const U& u) -> decltype(t + u) { return t + u; };
#endif
  };
  
  template <typename T,typename U>
  struct op_subtraction {
    BOOST_STATIC_CONSTANT(unsigned int, arity = 2);
    BOOST_STATIC_CONSTANT(op_class, classification = split);
#ifndef RK_ENABLE_CXX0X_FEATURES
    static T evaluate(const T& t, const U& u) { return t - u; };
#else
    static auto evaluate(T&& t, U&& u) -> decltype(std::move(t) - std::move(u)) { return std::move(t) - std::move(u); };
    static auto evaluate(const T& t, const U& u) -> decltype(t - u) { return t - u; };
#endif
  };
  
  template <typename T,typename U>
  struct op_negation {
    BOOST_STATIC_CONSTANT(unsigned int, arity = 1);
    BOOST_STATIC_CONSTANT(op_class, classification = split);
#ifndef RK_ENABLE_CXX0X_FEATURES
    static T evaluate(const T& t) { return -t; };
#else
    static auto evaluate(T&& t) -> decltype(-std::move(t)) { return -std::move(t); };
    static auto evaluate(const T& t) -> decltype(-t) { return -t; };
#endif
  };
  
  template <typename T,typename U>
  struct op_multiplication {
    BOOST_STATIC_CONSTANT(unsigned int, arity = 2);
    BOOST_STATIC_CONSTANT(op_class, classification = term);
#ifndef RK_ENABLE_CXX0X_FEATURES
    static T evaluate(const T& t, const U& u) { return t * u; };
#else
    static auto evaluate(T&& t, U&& u) -> decltype(std::move(t) * std::move(u)) { return std::move(t) * std::move(u); };
    static auto evaluate(const T& t, const U& u) -> decltype(t * u) { return t * u; };
#endif
  };
  
  template <typename T,typename U>
  struct op_division {
    BOOST_STATIC_CONSTANT(unsigned int, arity = 2);
    BOOST_STATIC_CONSTANT(op_class, classification = term);
#ifndef RK_ENABLE_CXX0X_FEATURES
    static T evaluate(const T& t, const U& u) { return t / u; };
#else
    static auto evaluate(T&& t, U&& u) -> decltype(std::move(t) / std::move(u)) { return std::move(t) / std::move(u); };
    static auto evaluate(const T& t, const U& u) -> decltype(t / u) { return t / u; };
#endif
  };
  
  
  
  
  
  
  
};







};

#endif












