/**
 * \file rpc_function.hpp
 *
 * This library declares a class template to represent a published and subscribed-to
 * remote procedure call (rpc) function.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date June 2014
 */

/*
 *    Copyright 2014 Sven Mikael Persson
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

#ifndef REAK_RPC_FUNCTION_HPP
#define REAK_RPC_FUNCTION_HPP

#include "detail/remote_function.hpp"
#include "detail/rpc_function_helpers.hpp"

#include <functional>
#include <variant>
#include "ReaK/core/rpc/rpc_exceptions.hpp"

namespace ReaK::rpc {

template <typename F>
class function;

template <typename R, typename... Args>
class function<R(Args...)> : public detail::remote_function {
 private:
  using FunctorType = std::function<R(Args...)>;
  std::variant<FunctorType, std::string> func_id;

 public:
  std::size_t get_params_hash() const override {
    static std::size_t h_val =
        detail::get_fnv_1a_hash(detail::concat_arg_types<0, Args...>());
    return h_val;
  }

  std::string get_host() const override {
    if (func_id.index() == 0) {
      return "localhost";
    }
    return std::get<std::string>(func_id);
  }

  using self = function<R(Args...)>;
  using input_tuple = typename detail::get_func_input_tuple<Args...>::type;

  function() : detail::remote_function(), func_id(std::string("")) {}

  function(const std::string& aName, const std::string& aHost)
      : detail::remote_function(aName), func_id(aHost) {}

  function(const std::string& aName, FunctorType aFuncObj)
      : detail::remote_function(aName), func_id(std::move(aFuncObj)) {
    detail::remote_function::publish();
  }

  ~function() override {
    if (func_id.index() == 0) {
      detail::remote_function::unpublish();
    }
  }

  void execute(serialization::iarchive& inputs,
               serialization::oarchive& outputs) override {
    if (func_id.index() != 0) {
      throw std::bad_function_call();
    }

    input_tuple input_data;
    try {
      // generically deserialize the inputs:
      detail::tuple_loader_impl<std::tuple_size<input_tuple>::value>::apply(
          inputs, input_data);
    } catch (std::exception& e) {
      throw marshalling_error(this, e.what());
    }

    // call the func_obj with the unrolled parameter tuple.
    detail::generic_return_type<R> result_val =
        detail::input_tuple_unroller<std::tuple_size<input_tuple>::value>::
            apply(std::get<FunctorType>(func_id), input_data);

    try {
      // code to populate the output stream.
      result_val.save_to(outputs);
      detail::save_output_from_input<0, input_tuple, Args...>(outputs,
                                                              input_data);
    } catch (std::exception& e) {
      throw marshalling_error(this, e.what());
    }
  }

  template <typename... Args2>
  R operator()(Args2&&... args) {
    if (func_id.index() == 0) {
      return std::get<FunctorType>(func_id)(std::forward<Args2>(args)...);
    }
    try {
      // Put args into the call preparations
      detail::call_preparations pre_data = this->prepare_for_call();
      detail::input_saver<Args...>::apply(*pre_data.p_aro,
                                          std::forward<Args2>(args)...);

      // Do the remote call
      detail::call_results res = this->do_remote_call(std::move(pre_data));

      // Retrieve args from p_ari
      detail::generic_return_type<R> result_val;
      result_val.load_from(*res.p_ari);
      detail::output_loader<Args...>::apply(*res.p_ari,
                                            std::forward<Args2>(args)...);

      return result_val.take_value();
    } catch (std::exception& e) {
      throw marshalling_error(this, e.what());
    }
  }

  void publish(const std::string& aName, std::function<R(Args...)> aFuncObj) {
    if ((func_id.index() == 0) && (aName != this->name)) {
      detail::remote_function::unpublish();
    }
    this->name = aName;
    func_id = std::move(aFuncObj);
    detail::remote_function::publish();
  }

  void from_remote(const std::string& aName, const std::string& aHost) {
    if (func_id.index() == 0) {
      detail::remote_function::unpublish();
    }
    this->name = aName;
    func_id = aHost;
  }
};

}  // namespace ReaK::rpc

#endif  // REAK_RPC_FUNCTION_HPP
