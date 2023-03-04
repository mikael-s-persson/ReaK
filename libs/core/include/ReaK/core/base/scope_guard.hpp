/**
 * \file scope_guard.hpp
 *
 * This library declares the "scope_guard<T>" class template which is based on Andrei Alexandrescu's
 * presentation of the topic at the "C++ and Beyond 2012" conference, called ScopeGuard11. This
 * particular implementation is the work of Mikael Persson, but the entire idea is public domain.
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

#ifndef REAK_SCOPE_GUARD_HPP
#define REAK_SCOPE_GUARD_HPP

namespace ReaK {

/**
 * This class template implements the idea of a scope-guard, which is a method
 * to attach a functor (or lambda) to the scope-exit, which can be used for clean
 * up in case of an exception being thrown, or dimissed if all went OK.
 * This particular implementation is the work of Mikael Persson, but the design was
 * introduced by Andrei Alexandrescu as "ScopeGuard11", which is public domain.
 * The following class is easier to use with either the make_scope_guard function (for
 * functor type deduction) or using the MACRO RK_SCOPE_EXIT:
 * \code
 * try {
 *    ...
 *   auto sg_close_file = ReaK::make_scope_guard([&]() {
 *     my_file.close();
 *   });
 *    ...
 *   sg_close_file.dismiss();
 *    ...
 * } catch(...) { };
 * \endcode
 * which is also equivalent to:
 * \code
 * try {
 *    ...
 *   auto sg_close_file = RK_SCOPE_EXIT {
 *     my_file.close();
 *   };
 *    ...
 *   sg_close_file.dismiss();
 *    ...
 * } catch(...) { };
 * \endcode
 * \tparam Functor A callable type (functor) that can be called with no argument (and discarded result) at scope exit.
 */
template <typename Functor>
class scope_guard {
 private:
  Functor f;
  bool isActive;

 public:
  /**
   * Construct the scope-guard with a functor to be called on the scope exit.
   * \param aF The functor to be called with no argument at scope exit.
   */
  explicit scope_guard(Functor aF) : f(std::move(aF)), isActive(true) {}
  /**
   * Destructor.
   */
  ~scope_guard() {
    if (isActive) {
      f();
    }
  }
  /**
   * Dismiss the scope-guard such that the functor won't be called.
   */
  void dismiss() { isActive = false; }

  scope_guard() = delete;
  scope_guard(const scope_guard&) = delete;
  scope_guard& operator=(const scope_guard&) = delete;

  /**
   * Standard move-constructor.
   */
  scope_guard(scope_guard&& rhs) noexcept
      : f(std::move(rhs.f)), isActive(rhs.isActive) {
    rhs.dismiss();
  }
  /**
   * Standard move-assignment operator.
   */
  scope_guard& operator=(scope_guard&& rhs) noexcept {
    if (isActive) {
      f();
    }
    f = std::move(rhs.f);
    isActive = rhs.isActive;
    rhs.dismiss();
  }
};

/**
 * This function template can be used to deduce the type of a functor and create a
 * scope-guard object out of that. This function template can also be used with a
 * lambda expression, as:
 * \code
 * try {
 *    ...
 *   auto sg_close_file = ReaK::make_scope_guard([&]() {
 *     my_file.close();
 *   });
 *    ...
 *   sg_close_file.dismiss();
 *    ...
 * } catch(...) { };
 * \endcode
 * which is also equivalent to:
 * \code
 * try {
 *    ...
 *   auto sg_close_file = RK_SCOPE_EXIT {
 *     my_file.close();
 *   };
 *    ...
 *   sg_close_file.dismiss();
 *    ...
 * } catch(...) { };
 * \endcode
 */
template <typename Functor>
scope_guard<Functor> make_scope_guard(Functor f) {
  return scope_guard<Functor>(std::move(f));
}

namespace detail {

struct scope_guard_on_exit {};

template <typename Functor>
scope_guard<Functor> operator+(scope_guard_on_exit /*unused*/, Functor&& f) {
  return scope_guard<Functor>(std::forward<Functor>(f));
}
}  // namespace detail

#define RK_SCOPE_EXIT ReaK::detail::scope_guard_on_exit() + [&]()

#define RK_SCOPE_EXIT_ROUTINE(X) \
  auto X = ReaK::detail::scope_guard_on_exit() + [&]()

#define RK_SCOPE_EXIT_DISMISS(X) X.dismiss()
}  // namespace ReaK

#endif
