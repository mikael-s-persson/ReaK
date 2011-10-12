/**
 * \file arithmetic_tuple.hpp
 * 
 * This library provides the arithmetic tuple class. This class is basically just a wrapper 
 * of either the std::tuple class or the boost::tuples::tuple class, either way, it provides 
 * a meta-programming interface that is equivalent to the wrapped class and with the addition 
 * of the support for all the basic arithmetic operators, requiring, of course, that these 
 * arithmetic operators are also available on all the types contained in the tuple.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2011
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

#ifndef REAK_ARITHMETIC_TUPLE_HPP
#define REAK_ARITHMETIC_TUPLE_HPP

#include "base/defs.hpp"  //For the RK_ENABLE_CXX0X_FEATURES compilation flag.

#ifdef RK_ENABLE_CXX0X_FEATURES
#include <tuple>
#include <type_traits>
#else
#include <boost/tuple/tuple.hpp>
#endif

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/prior.hpp>

#include "serialization/archiver.hpp"
#include "rtti/so_register_type.hpp"
#include "rtti/typed_primitives.hpp"

namespace ReaK {

  
  
namespace detail {
  
/*****************************************************************************************
                             Implementation details
*****************************************************************************************/


  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_addassign_impl( Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_addassign_impl( Tuple& lhs, const Tuple& rhs) {
    tuple_addassign_impl< boost::mpl::prior<Idx>,Tuple >(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) += get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx,
      boost::mpl::size_t<0>
    >,
  void >::type tuple_subassign_impl( Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx,
      boost::mpl::size_t<0>
    >,
  void >::type tuple_subassign_impl( Tuple& lhs, const Tuple& rhs) {
    tuple_subassign_impl< boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) -= get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_mulassign_impl( Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_mulassign_impl( Tuple& lhs, const Tuple& rhs) {
    tuple_mulassign_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) *= get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_divassign_impl( Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_divassign_impl( Tuple& lhs, const Tuple& rhs) {
    tuple_divassign_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) /= get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_smulassign_impl( Tuple& lhs, const Scalar& rhs) { };

  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_smulassign_impl( Tuple& lhs, const Scalar& rhs) {
    tuple_smulassign_impl<boost::mpl::prior<Idx>,Tuple,Scalar>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) *= rhs;
  };
  
  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_sdivassign_impl( Tuple& lhs, const Scalar& rhs) { };

  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_sdivassign_impl( Tuple& lhs, const Scalar& rhs) {
    tuple_sdivassign_impl<boost::mpl::prior<Idx>,Tuple,Scalar>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(lhs) /= rhs;
  };
  
  
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_add_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_add_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
    tuple_add_impl<boost::mpl::prior<Idx>,Tuple>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) + get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_sub_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_sub_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
    tuple_sub_impl<boost::mpl::prior<Idx>,Tuple>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) - get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_neg_impl( Tuple& result, const Tuple& lhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_neg_impl( Tuple& result, const Tuple& lhs) {
    tuple_neg_impl<boost::mpl::prior<Idx>,Tuple>(result,lhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = -get<boost::mpl::prior<Idx>::type::value>(lhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_mul_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_mul_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
    tuple_mul_impl<boost::mpl::prior<Idx>,Tuple>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) * get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_div_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_div_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
    tuple_div_impl<boost::mpl::prior<Idx>,Tuple>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) / get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_muls_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) { };

  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_muls_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) {
    tuple_muls_impl<boost::mpl::prior<Idx>,Tuple,Scalar>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) * rhs;
  };
  
  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_smul_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_smul_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) {
    tuple_smul_impl<boost::mpl::prior<Idx>,Tuple,Scalar>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = lhs * get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_divs_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) { };

  template <typename Idx, typename Tuple, typename Scalar>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_divs_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) {
    tuple_divs_impl<boost::mpl::prior<Idx>,Tuple,Scalar>(result,lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) / rhs;
  };
  
  
  
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_save_impl( serialization::oarchive& lhs, const Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_save_impl( serialization::oarchive& lhs, const Tuple& rhs) {
    tuple_save_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    lhs << get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_save_nvp_impl( serialization::oarchive& lhs, const std::pair< std::string, const Tuple& >& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_save_nvp_impl( serialization::oarchive& lhs, const std::pair< std::string, const Tuple& >& rhs) {
    tuple_save_nvp_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    std::stringstream ss(rhs.first); ss << "_q" << boost::mpl::prior<Idx>::type::value;
    lhs & serialization::make_save_nvp(ss.str(),get<boost::mpl::prior<Idx>::type::value>(rhs.second));
  };
  
  
  
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_load_impl( serialization::iarchive& lhs, Tuple& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_load_impl( serialization::iarchive& lhs, Tuple& rhs) {
    tuple_load_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    lhs >> get<boost::mpl::prior<Idx>::type::value>(rhs);
  };
  
  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::equal_to< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_load_nvp_impl( serialization::iarchive& lhs, const std::pair< std::string, Tuple& >& rhs) { };

  template <typename Idx, typename Tuple>
  inline 
  typename boost::enable_if< 
    boost::mpl::greater< 
      Idx, 
      boost::mpl::size_t<0> 
    >,
  void >::type tuple_load_nvp_impl( serialization::iarchive& lhs, const std::pair< std::string, Tuple& >& rhs) {
    tuple_load_nvp_impl<boost::mpl::prior<Idx>,Tuple>(lhs,rhs);
#ifdef RK_ENABLE_CXX0X_FEATURES
    using std::get;
#else
    using boost::tuples::get;
#endif
    std::stringstream ss(rhs.first); ss << "_q" << boost::mpl::prior<Idx>::type::value;
    lhs & serialization::make_load_nvp(ss.str(),get<boost::mpl::prior<Idx>::type::value>(rhs.second));
  };
  
  
/*****************************************************************************************
                           END OF Implementation details
*****************************************************************************************/

};



/**
 * This meta-function computes a bool integral constant if the given type is an arithmetic-tuple.
 * \tparam Tuple The type to be tested as being an arithmetic-tuple or not.
 */
template <typename Tuple>
struct is_arithmetic_tuple : 
  boost::mpl::false_ { };

/**
 * This meta-function computes an integral constant describing the size (or length) of an arithmetic-tuple.
 * \tparam Tuple The arithmetic-tuple type.
 */
template <typename Tuple>
struct arithmetic_tuple_size : 
  boost::mpl::size_t< 0 > { };

#ifdef RK_ENABLE_CXX0X_FEATURES

/**
 * This class template is a simple wrapper of a tuple with the addition of arithmetic operators. 
 * This class is basically just a wrapper of the std::tuple class, and it provides 
 * a meta-programming interface that is equivalent to std::tuple and with the addition 
 * of the support for all the basic arithmetic operators, requiring, of course, that these 
 * arithmetic operators are also available on all the types contained in the tuple.
 * \tparam T The types contained in the arithmetic-tuple.
 */
template <typename... T>
class arithmetic_tuple : public std::tuple< T... > {
  private:
    typedef std::tuple< T... > arithmetic_tuple_base_class;
  public:
    
#if 1
    constexpr arithmetic_tuple() : arithmetic_tuple_base_class() { };
    
    explicit arithmetic_tuple(const T&... t) : arithmetic_tuple_base_class(t...) { };
    
    template <typename... U>
    explicit arithmetic_tuple(U&&... u) noexcept : arithmetic_tuple_base_class(std::forward<U>(u)...) { };
    
    arithmetic_tuple(const arithmetic_tuple< T... >&) = default;
    arithmetic_tuple(arithmetic_tuple< T... >&&) = default;
    
    //TODO: missing other standard-specified constructors (with other tuple types, and std::pair).
    
#else
    using arithmetic_tuple_base_class::arithmetic_tuple_base_class;
#endif
};

/* Specialization, see general template docs. */
template <typename... T>
struct is_arithmetic_tuple< arithmetic_tuple< T... > > : 
  boost::mpl::true_ { };

/* Specialization, see general template docs. */
template <typename... T>
struct arithmetic_tuple_size< arithmetic_tuple< T... > > : 
  boost::mpl::size_t< sizeof... T > { };
  
/**
 * This function template can be used to create an arithmetic-tuple.
 * \tparam T The types contained in the arithmetic-tuple.
 * \param t The values that make up the arithmetic-tuple.
 * \return An arithmetic-tuple.
 */
template <typename... T>
inline 
arithmetic_tuple< T... > make_arithmetic_tuple(T&&... t) {
  return arithmetic_tuple< T... >(std::forward<T>(t)...);
};

#else

/**
 * This class template is a simple wrapper of a tuple with the addition of arithmetic operators.
 * This class is basically just a wrapper of the boost::tuples::tuple class, and it provides 
 * a meta-programming interface that is equivalent to boost::tuples::tuple and with the addition 
 * of the support for all the basic arithmetic operators, requiring, of course, that these 
 * arithmetic operators are also available on all the types contained in the tuple.
 * \tparam Tn The types contained in the arithmetic-tuple.
 */
template <typename T1, typename T2 = void, typename T3 = void, typename T4 = void, typename T5 = void, 
          typename T6 = void, typename T7 = void, typename T8 = void, typename T9 = void, typename T10 = void>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
                     const T4& t4 = T4(), const T5& t5 = T5(), const T6& t6 = T6(),
		     const T7& t7 = T7(), const T8& t8 = T8(), const T9& t9 = T9(), 
		     const T10& t10 = T10()) : arithmetic_tuple_base_class(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4 , typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10>
struct is_arithmetic_tuple< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > > : 
  boost::mpl::true_ { };

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4 , typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10>
struct arithmetic_tuple_size< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > > : 
  boost::mpl::size_t< boost::tuples::length<arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >::arithmetic_tuple_base_class>::value > { };
  


/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4(), const T5& t5 = T5(), const T6& t6 = T6(),
			      const T7& t7 = T7(), const T8& t8 = T8(), const T9& t9 = T9()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4,t5,t6,t7,t8,t9) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4(), const T5& t5 = T5(), const T6& t6 = T6(),
			      const T7& t7 = T7(), const T8& t8 = T8()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4,t5,t6,t7,t8) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5, T6, T7 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4(), const T5& t5 = T5(), const T6& t6 = T6(),
			      const T7& t7 = T7()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4,t5,t6,t7) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5, T6 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5, T6 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4(), const T5& t5 = T5(), const T6& t6 = T6()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4,t5,t6) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4, T5 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4, T5 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4(), const T5& t5 = T5()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4,t5) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3, T4 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3, T4 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3(),
			      const T4& t4 = T4()) : 
			      arithmetic_tuple_base_class(t1,t2,t3,t4) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2, T3 > {
  private:
    typedef boost::tuples::tuple< T1, T2, T3 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2(), const T3& t3 = T3()) : 
                              arithmetic_tuple_base_class(t1,t2,t3) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1, typename T2>
class arithmetic_tuple : public boost::tuples::tuple< T1, T2 > {
  private:
    typedef boost::tuples::tuple< T1, T2 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1(), const T2& t2 = T2()) : 
		              arithmetic_tuple_base_class(t1,t2) { };
    
    
    
};

/* Specialization, see general template docs. */
template <typename T1>
class arithmetic_tuple : public boost::tuples::tuple< T1 > {
  private:
    typedef boost::tuples::tuple< T1 > arithmetic_tuple_base_class;
  public:
    
    
    explicit arithmetic_tuple(const T1& t1 = T1()) : 
		              arithmetic_tuple_base_class(t1) { };
    
    
    
};



/* Specialization, see general template docs. */
template <typename T1>
inline 
arithmetic_tuple< T1 > make_arithmetic_tuple(const T1& t1) {
  return arithmetic_tuple< T1 >(t1);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2>
inline 
arithmetic_tuple< T1, T2 > make_arithmetic_tuple(const T1& t1, const T2& t2) {
  return arithmetic_tuple< T1, T2 >(t1,t2);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3>
inline 
arithmetic_tuple< T1, T2, T3 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3) {
  return arithmetic_tuple< T1, T2, T3 >(t1,t2,t3);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4>
inline 
arithmetic_tuple< T1, T2, T3, T4 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
							 const T4& t4) {
  return arithmetic_tuple< T1, T2, T3, T4 >(t1,t2,t3,t4);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
							     const T4& t4, const T5& t5) {
  return arithmetic_tuple< T1, T2, T3, T4, T5 >(t1,t2,t3,t4,t5);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5, T6 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
							         const T4& t4, const T5& t5, const T6& t6) {
  return arithmetic_tuple< T1, T2, T3, T4, T5, T6 >(t1,t2,t3,t4,t5,t6);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
								     const T4& t4, const T5& t5, const T6& t6,
								     const T7& t7) {
  return arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7 >(t1,t2,t3,t4,t5,t6,t7);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
									 const T4& t4, const T5& t5, const T6& t6,
									 const T7& t7, const T8& t8) {
  return arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8 >(t1,t2,t3,t4,t5,t6,t7,t8);
};

/* Specialization, see general template docs. */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
									     const T4& t4, const T5& t5, const T6& t6,
									     const T7& t7, const T8& t8, const T9& t9) {
  return arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 >(t1,t2,t3,t4,t5,t6,t7,t8,t9);
};

/**
 * This function template can be used to create an arithmetic-tuple.
 * \tparam Tn The types contained in the arithmetic-tuple.
 * \param tn The values that make up the arithmetic-tuple.
 * \return An arithmetic-tuple.
 */
template <typename T1, typename T2, typename T3, typename T4, typename T5,
          typename T6, typename T7, typename T8, typename T9, typename T10>
inline 
arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > make_arithmetic_tuple(const T1& t1, const T2& t2, const T3& t3, 
										  const T4& t4, const T5& t5, const T6& t6,
										  const T7& t7, const T8& t8, const T9& t9, const T10& t10) {
  return arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10);
};




#endif

  
/**
 * This function template is an overload of the addition operator on arithmetic-tuples. 
 * This function performs the addition of each element of the tuple, will only compile if 
 * all elements of the tuple support the addition operator.
 * \param lhs Left-hand side of the addition.
 * \param rhs Right-hand side of the addition.
 * \return Added tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) + get<0>(rhs), get<1>(lhs) + get<1>(rhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator +(const Tuple& lhs, const Tuple& rhs) {
  Tuple result;
  detail::tuple_add_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the subtraction operator on arithmetic-tuples. 
 * This function performs the subtraction of each element of the tuple, will only compile if 
 * all elements of the tuple support the subtraction operator.
 * \param lhs Left-hand side of the subtraction.
 * \param rhs Right-hand side of the subtraction.
 * \return Subtracted tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) - get<0>(rhs), get<1>(lhs) - get<1>(rhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator -(const Tuple& lhs, const Tuple& rhs) {
  Tuple result;
  detail::tuple_sub_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the negation operator on arithmetic-tuples. 
 * This function performs the negation of each element of the tuple, will only compile if 
 * all elements of the tuple support the negation and assignment operators.
 * \param lhs Left-hand side of the negation.
 * \return Negated tuple, equivalent to make_arithmetic_tuple(-get<0>(lhs), -get<1>(lhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator -(const Tuple& lhs) const {
  Tuple result;
  detail::tuple_neg_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs);
  return result;
};

/**
 * This function template is an overload of the add-and-assign operator on arithmetic-tuples. 
 * This function performs the add-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the add-and-assign operator.
 * \param lhs Left-hand side of the add-and-assign.
 * \param rhs Right-hand side of the add-and-assign.
 * \return Added-and-assigned tuple reference, equivalent to get<0>(lhs) += get<0>(rhs); get<1>(lhs) += get<1>(rhs); ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator +=(Tuple& lhs, const Tuple& rhs) {
  detail::tuple_addassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
  return lhs;
};
    
/**
 * This function template is an overload of the sub-and-assign operator on arithmetic-tuples. 
 * This function performs the sub-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the sub-and-assign operator.
 * \param lhs Left-hand side of the sub-and-assign.
 * \param rhs Right-hand side of the sub-and-assign.
 * \return Subtracted-and-assigned tuple reference, equivalent to get<0>(lhs) -= get<0>(rhs); get<1>(lhs) -= get<1>(rhs); ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator -=(Tuple& lhs, const Tuple& rhs) {
  detail::tuple_subassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
  return lhs;
};
    
/**
 * This function template is an overload of the multiplication operator on arithmetic-tuples. 
 * This function performs the multiplication of each element of the tuple, will only compile if 
 * all elements of the tuple support the multiplication operator.
 * \param lhs Left-hand side of the multiplication.
 * \param rhs Right-hand side of the multiplication.
 * \return Multiplied tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) * get<0>(rhs), get<1>(lhs) * get<1>(rhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator *(const Tuple& lhs, const Tuple& rhs) {
  Tuple result;
  detail::tuple_mul_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the division operator on arithmetic-tuples. 
 * This function performs the division of each element of the tuple, will only compile if 
 * all elements of the tuple support the division operator.
 * \param lhs Left-hand side of the division.
 * \param rhs Right-hand side of the division.
 * \return Divided tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) / get<0>(rhs), get<1>(lhs) / get<1>(rhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator /(const Tuple& lhs, const Tuple& rhs) {
  Tuple result;
  detail::tuple_div_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the scalar-multiplication operator on arithmetic-tuples. 
 * This function performs the scalar-multiplication of each element of the tuple, will only compile if 
 * all elements of the tuple support the scalar-multiplication operator.
 * \param lhs Left-hand side of the scalar-multiplication.
 * \param rhs Right-hand side of the scalar-multiplication (the scalar).
 * \return Multiplied tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) * rhs, get<1>(lhs) * rhs, ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same type again).
 */  
template <typename Tuple, typename Scalar>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator *(const Tuple& lhs, const Scalar& rhs) {
  Tuple result;
  detail::tuple_muls_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the scalar-multiplication operator on arithmetic-tuples. 
 * This function performs the scalar-multiplication of each element of the tuple, will only compile if 
 * all elements of the tuple support the scalar-multiplication operator.
 * \param lhs Left-hand side of the scalar-multiplication (the scalar).
 * \param rhs Right-hand side of the scalar-multiplication.
 * \return Multiplied tuple, equivalent to make_arithmetic_tuple(lhs * get<0>(rhs), lhs * get<1>(rhs), ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same type again).
 */  
template <typename Tuple, typename Scalar>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator *(const Scalar& lhs, const Tuple& rhs) {
  Tuple result;
  detail::tuple_smul_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the scalar-division operator on arithmetic-tuples. 
 * This function performs the scalar-division of each element of the tuple, will only compile if 
 * all elements of the tuple support the scalar-division operator.
 * \param lhs Left-hand side of the scalar-division.
 * \param rhs Right-hand side of the scalar-division (the scalar).
 * \return Divided tuple, equivalent to make_arithmetic_tuple(get<0>(lhs) / rhs, get<1>(lhs) / rhs, ... etc.).
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same type again).
 */  
template <typename Tuple, typename Scalar>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple >::type operator /(const Tuple& lhs, const Scalar& rhs) {
  Tuple result;
  detail::tuple_divs_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
  return result;
};
    
/**
 * This function template is an overload of the multiply-and-assign operator on arithmetic-tuples. 
 * This function performs the multiply-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the multiply-and-assign operator.
 * \param lhs Left-hand side of the multiply-and-assign.
 * \param rhs Right-hand side of the multiply-and-assign.
 * \return Multiplied tuple reference, equivalent to get<0>(lhs) *= get<0>(rhs); get<1>(lhs) *= get<1>(rhs); ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator *=(Tuple& lhs, const Tuple& rhs) {
  detail::tuple_mulassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
  return lhs;
};
    
/**
 * This function template is an overload of the divide-and-assign operator on arithmetic-tuples. 
 * This function performs the divide-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the divide-and-assign operator.
 * \param lhs Left-hand side of the divide-and-assign.
 * \param rhs Right-hand side of the divide-and-assign.
 * \return Divided tuple reference, equivalent to get<0>(lhs) /= get<0>(rhs); get<1>(lhs) /= get<1>(rhs); ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator /=(Tuple& lhs, const Tuple& rhs) {
  detail::tuple_divassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
  return lhs;
};
    
/**
 * This function template is an overload of the scalar-multiply-and-assign operator on arithmetic-tuples. 
 * This function performs the scalar-multiply-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the scalar-multiply-and-assign operator.
 * \param lhs Left-hand side of the scalar-multiply-and-assign.
 * \param rhs Right-hand side of the scalar-multiply-and-assign (the scalar).
 * \return Multiplied tuple reference, equivalent to get<0>(lhs) *= rhs; get<1>(lhs) *= rhs; ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same type again).
 */  
template <typename Tuple, typename Scalar>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator *=(Tuple& lhs, const Scalar& rhs) {
  detail::tuple_smulassign_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(lhs, rhs);
  return lhs;
};
    
/**
 * This function template is an overload of the scalar-divide-and-assign operator on arithmetic-tuples. 
 * This function performs the scalar-divide-and-assign of each element of the tuple, will only compile if 
 * all elements of the tuple support the scalar-divide-and-assign operator.
 * \param lhs Left-hand side of the scalar-divide-and-assign.
 * \param rhs Right-hand side of the scalar-divide-and-assign (the scalar).
 * \return Divided tuple reference, equivalent to get<0>(lhs) /= rhs; get<1>(lhs) /= rhs; ... etc.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 * \tparam Scalar Any type whose multiplication with all types in the Tuple are possible and closed (results in the same type again).
 */  
template <typename Tuple, typename Scalar>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
Tuple& >::type operator /=(Tuple& lhs, const Scalar& rhs) {
  detail::tuple_sdivassign_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(lhs, rhs);
  return lhs;
};


namespace serialization {


/**
 * This function template is an overload of the unnamed archive-output operator for arithmetic-tuples. 
 * This function performs the unnamed archive-output of each element of the tuple, will only compile if 
 * all elements of the tuple support the unnamed archive-output operator.
 * \param out The output archive.
 * \param rhs The arithmetic-tuple object to output on the archive.
 * \return The output archive.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
oarchive& >::type operator <<(oarchive& out, const Tuple& rhs) {
  detail::tuple_save_impl<arithmetic_tuple_size<Tuple>, Tuple>(out,rhs);
  return out;
};

/**
 * This function template is an overload of the named archive-output operator for arithmetic-tuples. 
 * This function performs the named archive-output of each element of the tuple, will only compile if 
 * all elements of the tuple support the named archive-output operator.
 * \param out The output archive.
 * \param rhs The arithmetic-tuple object to output on the archive.
 * \return The output archive.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
oarchive& >::type operator &(oarchive& out, const std::pair<std::string, const Tuple& >& rhs) {
  detail::tuple_save_nvp_impl<arithmetic_tuple_size<Tuple>, Tuple>(out,rhs);
  return out;
};

/**
 * This function template is an overload of the unnamed archive-input operator for arithmetic-tuples. 
 * This function performs the unnamed archive-input of each element of the tuple, will only compile if 
 * all elements of the tuple support the unnamed archive-input operator.
 * \param in The input archive.
 * \param rhs The arithmetic-tuple object to input from the archive.
 * \return The input archive.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
iarchive& >::type operator >>(iarchive& in, Tuple& rhs) {
  detail::tuple_load_impl<arithmetic_tuple_size<Tuple>, Tuple>(in,rhs);
  return in;
};

/**
 * This function template is an overload of the named archive-input operator for arithmetic-tuples. 
 * This function performs the named archive-input of each element of the tuple, will only compile if 
 * all elements of the tuple support the named archive-input operator.
 * \param in The input archive.
 * \param rhs The arithmetic-tuple object to input from the archive.
 * \return The input archive.
 * \tparam Tuple An arithmetic-tuple type (see is_arithmetic_tuple).
 */  
template <typename Tuple>
typename boost::enable_if< is_arithmetic_tuple<Tuple>,
iarchive& >::type operator &(iarchive& in, const std::pair<std::string, Tuple& >& rhs) {
  detail::tuple_load_nvp_impl<arithmetic_tuple_size<Tuple>, Tuple>(in,rhs);
  return in;
};


};



namespace rtti {

#ifdef RK_ENABLE_CXX0X_FEATURES
  
template <typename... T>
struct get_type_id< arithmetic_tuple< T... > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000002C);
  static std::string type_name() { return "arithmetic_tuple"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const arithmetic_tuple< T... >& save_type;
  typedef arithmetic_tuple< T... >& load_type;
};

template <typename Tail, typename... T>
struct get_type_info< arithmetic_tuple< T... >, Tail > {
  typedef detail::type_id< arithmetic_tuple< T... > , typename get_type_info_seq< T... >::template with_tail<Tail>::type > type;
  static std::string type_name() { return get_type_id< arithmetic_tuple< T... > >::type_name() + "<" + get_type_info_seq< T... >::type_name() + ">" + "," + Tail::type_name(); };
};

#else

template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10>
struct get_type_id< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0x0000002C);
  static std::string type_name() { return "arithmetic_tuple"; };
  static construct_ptr CreatePtr() { return NULL; };
  
  typedef const arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >& save_type;
  typedef arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >& load_type;
};

template <typename T1, typename T2, typename T3, typename T4, typename T5, 
          typename T6, typename T7, typename T8, typename T9, typename T10, 
	  typename Tail>
struct get_type_info< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >, Tail > {
  typedef detail::type_id< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > , typename get_type_info_seq< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >::template with_tail<Tail>::type > type;
  static std::string type_name() { return get_type_id< arithmetic_tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 > >::type_name() + "<" + get_type_info_seq< T1, T2, T3, T4, T5, T6, T7, T8, T9, T10 >::type_name() + ">" + "," + Tail::type_name(); };
};

#endif

};



};





#ifdef RK_ENABLE_CXX0X_FEATURES

namespace std {
  
/* Specialization, see general template docs. */
  template <typename... T>
  class tuple_size< ReaK::arithmetic_tuple<T...> > : 
    public integral_constant< size_t, sizeof... (T) > { };
    
/* Specialization, see general template docs. */
  template <typename... T>
  class tuple_size< const ReaK::arithmetic_tuple<T...> > : 
    public integral_constant< size_t, sizeof... (T) > { };
    
/* Specialization, see general template docs. */
  template <typename... T>
  class tuple_size< volatile ReaK::arithmetic_tuple<T...> > : 
    public integral_constant< size_t, sizeof... (T) > { };
    
/* Specialization, see general template docs. */
  template <typename... T>
  class tuple_size< const volatile ReaK::arithmetic_tuple<T...> > : 
    public integral_constant< size_t, sizeof... (T) > { };
  
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class tuple_element< Idx, ReaK::arithmetic_tuple<T...> > {
    public:
      typedef typename tuple_element< Idx, typename ReaK::arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class tuple_element< Idx, const ReaK::arithmetic_tuple<T...> > {
    public:
      typedef typename tuple_element< Idx, const typename ReaK::arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class tuple_element< Idx, volatile ReaK::arithmetic_tuple<T...> > {
    public:
      typedef typename tuple_element< Idx, volatile typename ReaK::arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class tuple_element< Idx, const volatile ReaK::arithmetic_tuple<T...> > {
    public:
      typedef typename tuple_element< Idx, const volatile typename ReaK::arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
};


namespace ReaK {
  
  template <int Idx, typename Tuple>
  class arithmetic_tuple_element {
    public:
      typedef typename std::tuple_element< Idx, Tuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class arithmetic_tuple_element< Idx, arithmetic_tuple<T...> > {
    public:
      typedef typename std::tuple_element< Idx, typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class arithmetic_tuple_element< Idx, const arithmetic_tuple<T...> > {
    public:
      typedef typename std::tuple_element< Idx, const typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class arithmetic_tuple_element< Idx, volatile arithmetic_tuple<T...> > {
    public:
      typedef typename std::tuple_element< Idx, volatile typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename... T>
  class arithmetic_tuple_element< Idx, const volatile arithmetic_tuple<T...> > {
    public:
      typedef typename std::tuple_element< Idx, const volatile typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
  };
  
};
  
#else

namespace boost {
  
namespace tuples {
  
/* Specialization, see general template docs. */
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  struct length< ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > : 
    ReaK::arithmetic_tuple_size< ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { };
    
/* Specialization, see general template docs. */
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  struct length< const ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > : 
    ReaK::arithmetic_tuple_size< ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { };
    
/* Specialization, see general template docs. */
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class length< volatile ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > : 
    ReaK::arithmetic_tuple_size< ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { };
    
/* Specialization, see general template docs. */
  template <typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class length< const volatile ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > : 
    ReaK::arithmetic_tuple_size< ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > { };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class element< Idx, ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename element< Idx, typename ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class element< Idx, const ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename element< Idx, const typename ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class>::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class element< Idx, volatile ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename element< Idx, volatile typename ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class>::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class element< Idx, const volatile ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename element< Idx, const volatile typename ReaK::arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class>::type type;
  };
  
};
  
};

namespace ReaK {
  
  template <int Idx, typename Tuple>
  class arithmetic_tuple_element {
    public:
      typedef typename boost::tuples::element< Idx, Tuple >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class arithmetic_tuple_element< Idx, arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename boost::tuples::element< Idx, typename arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class arithmetic_tuple_element< Idx, const arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename boost::tuples::element< Idx, const typename arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class arithmetic_tuple_element< Idx, volatile arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename boost::tuples::element< Idx, volatile typename arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class >::type type;
  };
  
/* Specialization, see general template docs. */
  template <int Idx, typename T1, typename T2, typename T3, typename T4, typename T5, 
            typename T6, typename T7, typename T8, typename T9, typename T10>
  class arithmetic_tuple_element< Idx, const volatile arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10> > {
    public:
      typedef typename boost::tuples::element< Idx, const volatile typename arithmetic_tuple<T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::arithmetic_tuple_base_class >::type type;
  };
  
};
  
  
#endif





#endif








