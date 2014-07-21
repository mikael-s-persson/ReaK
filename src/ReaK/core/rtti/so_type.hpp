/**
 * \file so_type.hpp
 * 
 * This library defines a set of template classes used to create the type descriptors for 
 * the types registered to the ReaK::rtti system.
 * 
 * \author Mikael Persson <mikael.s.persson@gmail.com>
 * \date april 2011
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

#ifndef REAK_SO_TYPE_HPP
#define REAK_SO_TYPE_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/shared_object_base.hpp>

#include <boost/config.hpp>
#include <boost/cstdint.hpp>

#include <string>
#include <sstream>
#include <set>
#include <boost/utility/enable_if.hpp>
#include <boost/function.hpp>

#include <boost/mpl/integral_c.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/type_traits.hpp>


/** Main namespace for ReaK */
namespace ReaK {


class shared_object; //forward-declaration.

/** Main namespace for ReaK's Serialization */
namespace serialization {

  class serializable; //forward-declaration
  
};

/** Main namespace for ReaK's Run-time Type Identification (RTTI) */
namespace rtti {
  
typedef ReaK::shared_ptr<shared_object> shared_object_shared_pointer;

typedef shared_object_shared_pointer (RK_CALL *construct_ptr)();

//this is really the only thing that the class needs to define
template <typename T>
struct get_type_id {
  BOOST_STATIC_CONSTANT(unsigned int, ID = T::rk_rtti_ID);
  static std::string type_name() { return T::rk_rtti_TypeName(); };
  static construct_ptr CreatePtr() { return T::rk_rtti_CreatePtr(); };
  
  typedef const serialization::serializable& save_type;
  typedef serialization::serializable& load_type;
};


template <unsigned int U>
struct get_type_id< boost::mpl::integral_c<unsigned int,U> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = U);
  static std::string type_name() { std::stringstream ss; ss << ID; return ss.str(); };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::true_ > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 2);
  static std::string type_name() { return "true"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <>
struct get_type_id< boost::mpl::false_ > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 1);
  static std::string type_name() { return "false"; };
  static construct_ptr CreatePtr() { return NULL; };
};

template <int I>
struct get_type_id< boost::mpl::int_<I> > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = I);
  static std::string type_name() { std::stringstream ss; ss << I; return ss.str(); };
  static construct_ptr CreatePtr() { return NULL; };
};


  
namespace detail {
  
  struct null_type_id { 
    BOOST_STATIC_CONSTANT(unsigned int, ID = 0);
    
  };
  
  template <typename T, typename Tail = null_type_id>
  struct type_id {
    typedef Tail tail;
    BOOST_STATIC_CONSTANT(unsigned int, ID = ::ReaK::rtti::get_type_id<T>::ID);
  };
  
//   template <typename T>
//   struct type_id<T, null_type_id> { 
//     BOOST_STATIC_CONSTANT(unsigned int, ID = get_type_id<T>::ID);
//   };
  
  template <typename T>
  struct type_id_count {
    BOOST_STATIC_CONSTANT(unsigned int, value = type_id_count< typename T::tail >::value + 1);
  };
  
  template <>
  struct type_id_count<null_type_id> {
    BOOST_STATIC_CONSTANT(unsigned int, value = 1);
  };
  
  template <typename T>
  struct get_type_id {
    unsigned int at(unsigned int i) const {
      if(i == 0)
        return T::ID;
      else
        return get_type_id< typename T::tail >().at(--i);
    };
  };
  
  template <>
  struct get_type_id<null_type_id> {
    unsigned int at(unsigned int) const { return 0; };
  };
  
};


struct null_type_info {
  typedef detail::null_type_id type;
  static std::string type_name() { return std::string(); };
};

template <>
struct get_type_id< null_type_info > {
  BOOST_STATIC_CONSTANT(unsigned int, ID = 0);
  static std::string type_name() { return ""; };
  static construct_ptr CreatePtr() { return NULL; };
};

#ifndef BOOST_NO_CXX11_VARIADIC_TEMPLATES

template <typename... T> struct get_type_info_seq;

template <typename T1, typename... T>
struct get_type_info_seq< T1, T... > {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef detail::type_id< T1,
            typename get_type_info_seq< T... >::template with_tail<Tail>::type > type;
  };
  
  static std::string type_name() { return get_type_id<T1>::type_name() + ","
                                        + get_type_info_seq< T... >::type_name(); };
};

template <>
struct get_type_info_seq<> {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef typename Tail::type type;
  };
  
  static std::string type_name() { return ""; };
};

#else

template <typename T1, typename T2 = void, typename T3 = void, typename T4 = void, typename T5 = void, 
          typename T6 = void, typename T7 = void, typename T8 = void, typename T9 = void, typename T10 = void>
struct get_type_info_seq {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef detail::type_id< T1, 
            typename get_type_info_seq< T2, T3, T4, T5, T6, T7, T8, T9, T10 >::template with_tail<Tail>::type > type;
  };
  
  static std::string type_name() { return get_type_id<T1>::type_name() + ","
                                        + get_type_info_seq< T2, T3, T4, T5, T6, T7, T8, T9, T10 >::type_name(); };
};

template <typename T1>
struct get_type_info_seq<T1,void,void,void,void,void,void,void,void,void> {
  template <typename Tail = null_type_info>
  struct with_tail {
    typedef detail::type_id< T1, typename Tail::type > type;
  };
  
  static std::string type_name() { return get_type_id<T1>::type_name(); };
};

#endif


template <typename T, typename Tail = null_type_info>
struct get_type_info {
  typedef detail::type_id<T, typename Tail::type> type;
  static std::string type_name() { return get_type_id<T>::type_name() + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename> class U, typename T, typename Tail>
struct get_type_info< U<T>, Tail > {
  typedef detail::type_id< U<T>, typename get_type_info<T, Tail>::type > type;
  static std::string type_name() { return get_type_id< U<T> >::type_name() + "<" + get_type_id<T>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};


template <template <typename,typename> class U, typename T1, 
                                                typename T2, typename Tail>
struct get_type_info< U<T1,
                        T2>, Tail > {
  typedef detail::type_id< U<T1,
                             T2>, typename get_type_info<T1, 
                                           get_type_info<T2, Tail> >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2> >::type_name() + "<" + get_type_id<T1>::type_name() + "," 
                                                                                  + get_type_id<T2>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3, Tail> > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4, Tail> > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5, Tail> > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, 
                                       typename T6, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5,
                        T6>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5,
                             T6>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5,
                                           get_type_info<T6, Tail> > > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5,
                                                         T6> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ","
                                                                                  + get_type_id<T6>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, 
                                       typename T6, 
                                       typename T7, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5,
                        T6,
                        T7>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5,
                             T6,
                             T7>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5,
                                           get_type_info<T6,
                                           get_type_info<T7, Tail> > > > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5,
                                                         T6,
                                                         T7> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ","
                                                                                  + get_type_id<T6>::type_name() + ","
                                                                                  + get_type_id<T7>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, 
                                       typename T6, 
                                       typename T7, 
                                       typename T8, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5,
                        T6,
                        T7,
                        T8>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5,
                             T6,
                             T7,
                             T8>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5,
                                           get_type_info<T6,
                                           get_type_info<T7,
                                           get_type_info<T8, Tail> > > > > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5,
                                                         T6,
                                                         T7,
                                                         T8> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ","
                                                                                  + get_type_id<T6>::type_name() + ","
                                                                                  + get_type_id<T7>::type_name() + ","
                                                                                  + get_type_id<T8>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, 
                                       typename T6, 
                                       typename T7, 
                                       typename T8, 
                                       typename T9, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5,
                        T6,
                        T7,
                        T8,
                        T9>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5,
                             T6,
                             T7,
                             T8,
                             T9>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5,
                                           get_type_info<T6,
                                           get_type_info<T7,
                                           get_type_info<T8,
                                           get_type_info<T9, Tail> > > > > > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5,
                                                         T6,
                                                         T7,
                                                         T8,
                                                         T9> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ","
                                                                                  + get_type_id<T6>::type_name() + ","
                                                                                  + get_type_id<T7>::type_name() + ","
                                                                                  + get_type_id<T8>::type_name() + ","
                                                                                  + get_type_id<T9>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <template <typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename,
                    typename> class U, typename T1, 
                                       typename T2, 
                                       typename T3, 
                                       typename T4, 
                                       typename T5, 
                                       typename T6, 
                                       typename T7, 
                                       typename T8, 
                                       typename T9, 
                                       typename T10, typename Tail>
struct get_type_info< U<T1,
                        T2,
                        T3,
                        T4,
                        T5,
                        T6,
                        T7,
                        T8,
                        T9,
                        T10>, Tail > {
  typedef detail::type_id< U<T1,
                             T2,
                             T3,
                             T4,
                             T5,
                             T6,
                             T7,
                             T8,
                             T9,
                             T10>, typename get_type_info<T1, 
                                           get_type_info<T2,
                                           get_type_info<T3,
                                           get_type_info<T4,
                                           get_type_info<T5,
                                           get_type_info<T6,
                                           get_type_info<T7,
                                           get_type_info<T8,
                                           get_type_info<T9,
                                           get_type_info<T10, Tail> > > > > > > > > >::type > type;
  static std::string type_name() { return get_type_id< U<T1,
                                                         T2,
                                                         T3,
                                                         T4,
                                                         T5,
                                                         T6,
                                                         T7,
                                                         T8,
                                                         T9,
                                                         T10> >::type_name() + "<" + get_type_id<T1>::type_name() + ","
                                                                                  + get_type_id<T2>::type_name() + ","
                                                                                  + get_type_id<T3>::type_name() + ","
                                                                                  + get_type_id<T4>::type_name() + ","
                                                                                  + get_type_id<T5>::type_name() + ","
                                                                                  + get_type_id<T6>::type_name() + ","
                                                                                  + get_type_id<T7>::type_name() + ","
                                                                                  + get_type_id<T8>::type_name() + ","
                                                                                  + get_type_id<T9>::type_name() + ","
                                                                                  + get_type_id<T10>::type_name() + ">" + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};


template <unsigned int U, typename Tail>
struct get_type_info< boost::mpl::integral_c<unsigned int,U>, Tail > {
  typedef detail::type_id<boost::mpl::integral_c<unsigned int,U>, typename Tail::type> type;
  static std::string type_name() { return get_type_id< boost::mpl::integral_c<unsigned int,U> >::type_name() + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <typename Tail>
struct get_type_info< boost::mpl::true_, Tail > {
  typedef detail::type_id<boost::mpl::true_, typename Tail::type> type;
  static std::string type_name() { return get_type_id<boost::mpl::true_>::type_name() + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <typename Tail>
struct get_type_info< boost::mpl::false_, Tail > {
  typedef detail::type_id<boost::mpl::false_, typename Tail::type> type;
  static std::string type_name() { return get_type_id<boost::mpl::false_>::type_name() + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};

template <int I, typename Tail>
struct get_type_info< boost::mpl::int_<I>, Tail > {
  typedef detail::type_id<boost::mpl::int_<I>, typename Tail::type> type;
  static std::string type_name() { return get_type_id< boost::mpl::int_<I> >::type_name() + (boost::is_same< Tail, null_type_info >::value ? "" : "," + Tail::type_name()); };
};





/**
 * This class is used to identify and register in the platform a new object
 * type. When a plugin is loaded, at initialization of the plugin, it is responsible
 * for registering all the new types it contains.
 */
class so_type : public shared_object_base {
private:
  so_type(const so_type&);
  so_type& operator =(const so_type&);
public:
  typedef ReaK::weak_ptr<so_type> weak_pointer;
  typedef ReaK::shared_ptr<so_type> shared_pointer;
protected:
  
  so_type() : shared_object_base() { };

  ///This function finds a TypeID in the descendants (recusively) of this.
  virtual weak_pointer RK_CALL findDescendant_impl(const unsigned int* aTypeID ) const = 0;

  ///This function gets the number of direct descendants of this.
  virtual unsigned int RK_CALL getDescendantCount_impl() const = 0;

  ///This function gets a Type record by index in the direct descendants of this.
  virtual shared_pointer RK_CALL getDescendant_impl(unsigned int aIndex) const = 0;

  ///This function checks if a typeID is parent to this.
  virtual weak_pointer RK_CALL findAncestor_impl(const unsigned int* aTypeID ) const = 0;
  
  static bool compare_equal(const unsigned int* pid1, const unsigned int* pid2);

public:
 
  virtual void RK_CALL destroy() { delete this; };

  virtual ~so_type() { };

  ///This function adds a Descendant of this.
  virtual shared_pointer RK_CALL addDescendant(const shared_pointer& aObj ) = 0;
  
  virtual weak_pointer RK_CALL addAncestor(shared_pointer& aThis, const weak_pointer& aObj) = 0;
    
  ///This function finds a TypeID in the descendants (recusively) of this.
  weak_pointer RK_CALL findDescendant(const unsigned int* aTypeID ) const {
    return this->findDescendant_impl(aTypeID);
  };

  ///This function gets the number of direct descendants of this.
  unsigned int RK_CALL getDirectDescendantCount() const {
    return this->getDescendantCount_impl();
  };

  ///This function gets a type record by index in the direct descendants of this.
  shared_pointer RK_CALL getDirectDescendant(unsigned int aIndex) const {
    return this->getDescendant_impl(aIndex);
  };

  ///This function checks if a typeID is parent to this.
  weak_pointer RK_CALL findAncestor(const unsigned int* aTypeID ) const {
    return this->findAncestor_impl(aTypeID);
  };
  
  ///This function finds a TypeID in the descendants (recusively) of this.
  weak_pointer RK_CALL findDescendant(const shared_pointer& aTypeID ) const {
    return this->findDescendant_impl(aTypeID->TypeID_begin());
  };

  ///This function checks if a typeID is parent to this.
  weak_pointer RK_CALL findAncestor(shared_pointer& aThis, const shared_pointer& aTypeID ) const {
    return this->findAncestor_impl(aTypeID->TypeID_begin());
  };

  ///This function inserts this into a global repo.
  virtual void RK_CALL insertToRepo(const shared_pointer& aThis, shared_pointer& aRepo) = 0;

  virtual const unsigned int* RK_CALL TypeID_begin() const = 0;
  
  virtual unsigned int RK_CALL TypeVersion() const = 0;

  virtual const std::string& RK_CALL TypeName() const = 0;

  virtual shared_object_shared_pointer RK_CALL CreateObject() const = 0;
  
  virtual bool isConcrete() const = 0;
  
  friend bool operator ==(const so_type& t1, const so_type& t2) {
    return compare_equal(t1.TypeID_begin(),t2.TypeID_begin());
  };
  
  friend bool operator !=(const so_type& t1, const so_type& t2) {
    return !compare_equal(t1.TypeID_begin(),t2.TypeID_begin());
  };

};


class so_type_impl : public so_type {
protected:
  static bool compare_shared(const shared_pointer& t1, const shared_pointer& t2);
  static bool compare_weak(const weak_pointer& t1, const weak_pointer& t2);
  
  typedef bool (*compare_shared_t)(const shared_pointer&,const shared_pointer&);
  typedef bool (*compare_weak_t)(const weak_pointer&,const weak_pointer&);
  
  std::set< shared_pointer, compare_shared_t > mDescendants;
  std::set< weak_pointer, compare_weak_t > mAncestors;
  
  so_type_impl();
  
  ///This function finds a TypeID in the descendants (recusively) of this.
  virtual weak_pointer RK_CALL findDescendant_impl(const unsigned int* aTypeID ) const;

  ///This function gets the number of direct descendants of this.
  virtual unsigned int RK_CALL getDescendantCount_impl() const;

  ///This function gets a Type record by index in the direct descendants of this.
  virtual shared_pointer RK_CALL getDescendant_impl(unsigned int aIndex) const;

  ///This function checks if a typeID is parent to this.
  virtual weak_pointer RK_CALL findAncestor_impl(const unsigned int* aTypeID ) const;
  
public:
  
  ///This function adds a Descendant of this.
  virtual shared_pointer RK_CALL addDescendant(const shared_pointer& aObj);

  virtual weak_pointer RK_CALL addAncestor(shared_pointer& aThis, const weak_pointer& aObj);
    
  ///This function inserts this into a global repo.
  virtual void RK_CALL insertToRepo(const shared_pointer& aThis, shared_pointer& aRepo);

};


template <typename T, unsigned int Version = 1>
class so_type_descriptor : public so_type_impl {
private:
  BOOST_STATIC_CONSTANT(unsigned int, mTypeIDLength = detail::type_id_count< typename get_type_info<T>::type >::value);
  BOOST_STATIC_CONSTANT(unsigned int, mTypeVersion = Version);
  unsigned int mTypeID[mTypeIDLength];
  std::string mTypeName;
  construct_ptr mConstruct;
public:
  so_type_descriptor() : so_type_impl() {
    for(unsigned int i = 0; i < mTypeIDLength; ++i) {
      mTypeID[i] = detail::get_type_id< typename get_type_info<T>::type >().at(i);
    };
    mTypeName = get_type_info<T>::type_name();
    mConstruct = get_type_id<T>::CreatePtr();
  };

  virtual ~so_type_descriptor() {
    RK_NOTICE(6,"type " << mTypeName << "\t version " << mTypeVersion << " just deleted!");
  };
  
  virtual const unsigned int* RK_CALL TypeID_begin() const { return mTypeID; };

  virtual unsigned int RK_CALL TypeVersion() const { return mTypeVersion; };

  virtual const std::string& RK_CALL TypeName() const { return mTypeName; };

  virtual shared_object_shared_pointer RK_CALL CreateObject() const {
    if(mConstruct)
      return mConstruct();
    else
      return shared_object_shared_pointer();
  };
  
  virtual bool isConcrete() const {
    return (mConstruct);
  };
  
};


// template <typename T, unsigned int Version>
// unsigned int so_type_descriptor<T,Version>::mTypeID[so_type_descriptor<T,Version>::mTypeIDLength];
// 
// template <typename T, unsigned int Version>
// std::string so_type_descriptor<T,Version>::mTypeName;
// 
// template <typename T, unsigned int Version>
// construct_ptr so_type_descriptor<T,Version>::mConstruct;



namespace detail {
  
  class dummy_so_type : public so_type_impl {
  private:
    const unsigned int* mTypeID;
    std::string mTypeName;
        
  public:
    dummy_so_type(const unsigned int* aTypeID) : mTypeID(aTypeID), mTypeName("Root") { };
    virtual ~dummy_so_type() { };
    
    virtual const unsigned int* RK_CALL TypeID_begin() const { return mTypeID; };
    virtual unsigned int RK_CALL TypeVersion() const { return 0; };
    virtual const std::string& RK_CALL TypeName() const { return mTypeName; };
    virtual shared_object_shared_pointer RK_CALL CreateObject() const { return shared_object_shared_pointer(); };
    virtual bool isConcrete() const { return false; };
    
  };
  
  
  
};





};

};

#endif


#include "typed_primitives.hpp"
#include "typed_containers.hpp"










