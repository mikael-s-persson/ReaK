/**
 * \file color.hpp
 *
 * This library declares a class to represent a color.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date October 2012
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

#ifndef REAK_COLOR_HPP
#define REAK_COLOR_HPP

#include <ReaK/core/base/defs.hpp>
#include <ReaK/core/base/serializable.hpp>

/** Main namespace for ReaK */
namespace ReaK {

/** Main namespace for ReaK.Geometry */
namespace geom {


/** This class represents a color. */
struct color : public serializable {
public:
  double R;
  double G;
  double B;

  /**
   * Default constructor.
   * \param aR The red component (0.0 to 1.0).
   * \param aR The green component (0.0 to 1.0).
   * \param aR The blue component (0.0 to 1.0).
   */
  color( double aR = 1.0, double aG = 1.0, double aB = 1.0 ) : R( aR ), G( aG ), B( aB ){};

  /**
   * Default destructor.
   */
  virtual ~color(){};


  /*******************************************************************************
                     ReaK's RTTI and Serialization interfaces
  *******************************************************************************/

  virtual void RK_CALL save( serialization::oarchive& A, unsigned int ) const {
    A& std::pair< std::string, typename ReaK::rtti::get_type_id< double >::save_type >( "red", R )
      & std::pair< std::string, typename ReaK::rtti::get_type_id< double >::save_type >( "green", G )
      & std::pair< std::string, typename ReaK::rtti::get_type_id< double >::save_type >( "blue", B );
  };
  virtual void RK_CALL load( serialization::iarchive& A, unsigned int ) {
    A& std::pair< std::string, typename ReaK::rtti::get_type_id< double >::load_type >( "red", R )
      & std::pair< std::string, typename ReaK::rtti::get_type_id< double >::load_type >( "green", G )
      & std::pair< std::string, typename ReaK::rtti::get_type_id< double >::load_type >( "blue", B );
  };

  RK_RTTI_REGISTER_CLASS_1BASE( color, 1, serializable )
};
};


namespace rtti {

template <>
struct get_type_id< geom::color > {
  BOOST_STATIC_CONSTANT( unsigned int, ID = 0x00000031 );
#ifdef RK_RTTI_USE_CONSTEXPR_STRINGS
  BOOST_STATIC_CONSTEXPR auto type_name = RK_LSA( "color" );
#else
  static const char* type_name() BOOST_NOEXCEPT { return "color"; };
#endif
  static construct_ptr CreatePtr() BOOST_NOEXCEPT { return nullptr; };

  typedef const serializable& save_type;
  typedef serializable& load_type;
};
};
};

#endif
