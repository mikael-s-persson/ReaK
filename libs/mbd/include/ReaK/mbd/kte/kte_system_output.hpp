/**
 * \file kte_system_output.hpp
 *
 * This library defines the base class for system outputs to a KTE model. A system output
 * is simply a vector of values which serve as an output to a KTE model. This model is useful
 * when using a KTE model into a state-space system definition.
 *
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date May 2011
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

#ifndef REAK_KTE_SYSTEM_OUTPUT_HPP
#define REAK_KTE_SYSTEM_OUTPUT_HPP

#include "ReaK/core/base/named_object.hpp"

#include <vector>

namespace ReaK::kte {

/**
 * This class is a base class for system outputs of a KTE model. A system output
 * is simply a vector of values which serve as an output to a KTE model. This model is useful
 * when using a KTE model into a state-space system definition.
 */
class system_output : public virtual named_object {
 public:
  /**
   * Constructs a system input class with the given name.
   */
  explicit system_output(const std::string& aName = "") {
    this->setName(aName);
  }

  /**
   * Destructor.
   */
  ~system_output() override = default;
  ;

  /**
   * Returns the number of output variables provided by this system output.
   * \return the number of output variables provided by this system output.
   */
  virtual unsigned int getOutputCount() const = 0;

  /**
   * Returns the output variable at index i, with read-only access.
   * \param i The index of the output variable.
   * \return The variable at index i.
   */
  virtual double getOutput(unsigned int i) const = 0;

  void save(serialization::oarchive& A,
            unsigned int /*unused*/) const override {
    ReaK::named_object::save(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  void load(serialization::iarchive& A, unsigned int /*unused*/) override {
    ReaK::named_object::load(
        A, named_object::getStaticObjectType()->TypeVersion());
  }

  RK_RTTI_MAKE_ABSTRACT_1BASE(system_output, 0xC2100034, 1, "system_output",
                              named_object)
};

}  // namespace ReaK::kte

#endif
