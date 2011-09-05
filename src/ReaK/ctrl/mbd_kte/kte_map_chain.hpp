/**
 * \file kte_map_chain.hpp
 *
 * This library declares the kte_map_chain class which allows KTE models to be grouped into a
 * linear chain of KTE models. The principle is based on the fact that if all the KTEs
 * compute their kinematics in order, the dynamics can be computed in the exact reverse order.
 *
 * \author Mikael Persson, <mikael.s.persson@gmail.com>
 * \date March 2010
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

#ifndef REAK_KTE_MAP_CHAIN_HPP
#define REAK_KTE_MAP_CHAIN_HPP

#include "kte_map.hpp"
#include <vector>

namespace ReaK {

namespace kte {



/**
 * This class is a container of linearly chained kinetostatic transmission elements (KTEs).
 */
class kte_map_chain : public kte_map {
  private:
    std::vector< shared_pointer<kte_map>::type > mKTEs; ///< Stores the list of KTEs

  public:

    /**
     * Default constructor.
     */
    kte_map_chain(const std::string& aName = "") : kte_map(aName), mKTEs() { };

    /**
     * Default destructor.
     */
    virtual ~kte_map_chain() { };

    virtual void doMotion(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type()) {
      std::vector< shared_pointer<kte_map>::type >::iterator it = mKTEs.begin();
      for(;it != mKTEs.end();++it) {
        (*it)->doMotion(aFlag,aStorage);
      };
    };

    virtual void doForce(kte_pass_flag aFlag = nothing, const shared_pointer<frame_storage>::type& aStorage = shared_pointer<frame_storage>::type()) {
      std::vector< shared_pointer<kte_map>::type >::reverse_iterator rit = mKTEs.rbegin();
      for(;rit != mKTEs.rend();++rit) {
        (*rit)->doForce(aFlag,aStorage);
      };
    };

    virtual void clearForce() {
      std::vector< shared_pointer<kte_map>::type >::iterator it = mKTEs.begin();
      for(;it != mKTEs.end();++it)
        (*it)->clearForce();
    };

    /**
     * This method appends a KTE model to the linear chain of models.
     * \pre Any state.
     * \post The chain will contain an additional KTE.
     * \param aKTE The KTE model to add to this chain.
     * \return reference to this chain (to chain the << operators).
     */
    kte_map_chain& operator <<(const shared_pointer<kte_map>::type& aKTE) {
      if(aKTE)
        mKTEs.push_back(aKTE);
      return *this;
    };

    virtual void RK_CALL save(serialization::oarchive& A, unsigned int) const {
      ReaK::named_object::save(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_SAVE_WITH_NAME(mKTEs);
    };

    virtual void RK_CALL load(serialization::iarchive& A, unsigned int) {
      ReaK::named_object::load(A,named_object::getStaticObjectType()->TypeVersion());
      A & RK_SERIAL_LOAD_WITH_NAME(mKTEs);
    };

    RK_RTTI_MAKE_CONCRETE_1BASE(kte_map_chain,0xC2100002,1,"kte_map_chain",kte_map)

};


};

};

#endif





