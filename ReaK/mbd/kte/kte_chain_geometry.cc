
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

#include "ReaK/mbd/kte/kte_chain_geometry.h"

#include "ReaK/mbd/kte/inertia.h"
#include "ReaK/mbd/kte/kte_chain_visitation.h"
#include "ReaK/mbd/kte/kte_map_chain.h"
#include "ReaK/mbd/kte/prismatic_joint.h"
#include "ReaK/mbd/kte/revolute_joint.h"
#include "ReaK/mbd/kte/rigid_link.h"

#include <map>
#include <string>

namespace ReaK::kte {

namespace detail {

using kte_geom_visit_2D_types = std::tuple<rigid_link_2D, prismatic_joint_2D,
                                           revolute_joint_2D, inertia_2D>;
using kte_geom_visit_3D_types = std::tuple<rigid_link_3D, prismatic_joint_3D,
                                           revolute_joint_3D, inertia_3D>;

struct kte_geom_2D_visitor {
  const kte_chain_geometry_2D* parent;

  explicit kte_geom_2D_visitor(const kte_chain_geometry_2D* aParent)
      : parent(aParent) {}

  template <typename KTEType>
  void connect_to_base_or_end(const KTEType& aObj) {
    const std::string& obj_name = aObj.get_name();

    auto it = parent->mGeomList.find(obj_name + "_base");
    if ((it != parent->mGeomList.end()) && (!it->second.empty())) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.BaseFrame());
        }
      }
    }

    auto it_prox = parent->mProxyShapeList.find(obj_name + "_base");
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (!it_prox->second.empty())) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.BaseFrame());
        }
      }
    }

    it = parent->mGeomList.find(obj_name + "_end");
    if ((it != parent->mGeomList.end()) && (!it->second.empty())) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.EndFrame());
        }
      }
    }

    it_prox = parent->mProxyShapeList.find(obj_name + "_end");
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (!it_prox->second.empty())) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.EndFrame());
        }
      }
    }
  }

  void operator()(const rigid_link_2D& aObj) { connect_to_base_or_end(aObj); }

  void operator()(const prismatic_joint_2D& aObj) {
    connect_to_base_or_end(aObj);
  }

  void operator()(const revolute_joint_2D& aObj) {
    connect_to_base_or_end(aObj);
  }

  void operator()(const inertia_2D& aObj) const {
    const std::string& obj_name = aObj.get_name();

    auto it = parent->mGeomList.find(obj_name);
    if ((it != parent->mGeomList.end()) &&
        (static_cast<unsigned int>(!it->second.empty()) != 0U)) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.CenterOfMass()->mFrame);
        }
      }
    }

    auto it_prox = parent->mProxyShapeList.find(obj_name);
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (static_cast<unsigned int>(!it_prox->second.empty()) != 0U)) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.CenterOfMass()->mFrame);
        }
      }
    }
  }
};

struct kte_geom_3D_visitor {
  const kte_chain_geometry_3D* parent;

  explicit kte_geom_3D_visitor(const kte_chain_geometry_3D* aParent)
      : parent(aParent){};

  template <typename KTEType>
  void connect_to_base_or_end(const KTEType& aObj) {
    const std::string& obj_name = aObj.get_name();

    auto it = parent->mGeomList.find(obj_name + "_base");
    if ((it != parent->mGeomList.end()) && (!it->second.empty())) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.BaseFrame());
        }
      }
    }

    auto it_prox = parent->mProxyShapeList.find(obj_name + "_base");
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (!it_prox->second.empty())) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.BaseFrame());
        }
      }
    }

    it = parent->mGeomList.find(obj_name + "_end");
    if ((it != parent->mGeomList.end()) && (!it->second.empty())) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.EndFrame());
        }
      }
    }

    it_prox = parent->mProxyShapeList.find(obj_name + "_end");
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (!it_prox->second.empty())) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.EndFrame());
        }
      }
    }
  }

  void operator()(const rigid_link_3D& aObj) { connect_to_base_or_end(aObj); }

  void operator()(const prismatic_joint_3D& aObj) {
    connect_to_base_or_end(aObj);
  }

  void operator()(const revolute_joint_3D& aObj) {
    connect_to_base_or_end(aObj);
  }

  void operator()(const inertia_3D& aObj) const {
    const std::string& obj_name = aObj.get_name();

    auto it = parent->mGeomList.find(obj_name);
    if ((it != parent->mGeomList.end()) &&
        (static_cast<unsigned int>(!it->second.empty()) != 0U)) {
      for (const auto& i : it->second) {
        if (i.mGeom) {
          i.mGeom->setAnchor(aObj.CenterOfMass()->mFrame);
        }
      }
    }

    auto it_prox = parent->mProxyShapeList.find(obj_name);
    if ((it_prox != parent->mProxyShapeList.end()) &&
        (static_cast<unsigned int>(!it_prox->second.empty()) != 0U)) {
      for (const auto& i : it_prox->second) {
        if (i) {
          i->setAnchor(aObj.CenterOfMass()->mFrame);
        }
      }
    }
  }
};

}  // namespace detail

std::pair<std::shared_ptr<geom::colored_model_2D>,
          std::shared_ptr<geom::proxy_query_model_2D>>
kte_chain_geometry_2D::attachToKTEChain(const kte_map_chain& aKTEChain) const {
  // first, visit the KTE chain to associated each geometry to its KTE-related anchor.
  auto vis = detail::kte_geom_2D_visitor(this);
  visit_kte_chain<detail::kte_geom_visit_2D_types>(vis, aKTEChain);

  // then, construct the colored_model_2D from the resulting (attached) geometries:
  auto result_geom =
      std::make_shared<geom::colored_model_2D>(aKTEChain.get_name() + "_geom");
  for (const auto& it : mGeomList) {
    for (const auto& i : it.second) {
      if (i.mGeom) {
        result_geom->addAnchor(i.mGeom->getAnchor());
        result_geom->addElement(i.mColor, i.mGeom);
      }
    }
  }

  // then, construct the proxy_query_model_2D from the resulting (attached) geometries:
  auto result_prox = std::make_shared<geom::proxy_query_model_2D>(
      aKTEChain.get_name() + "_proxy");
  for (const auto& it : mProxyShapeList) {
    for (const auto& i : it.second) {
      if (i) {
        result_prox->addShape(i);
      }
    }
  }

  return {result_geom, result_prox};
}

std::pair<std::shared_ptr<geom::colored_model_3D>,
          std::shared_ptr<geom::proxy_query_model_3D>>
kte_chain_geometry_3D::attachToKTEChain(const kte_map_chain& aKTEChain) const {
  auto vis = detail::kte_geom_3D_visitor(this);
  visit_kte_chain<detail::kte_geom_visit_3D_types>(vis, aKTEChain);

  // then, construct the colored_model_3D from the resulting (attached) geometries:
  auto result_geom =
      std::make_shared<geom::colored_model_3D>(aKTEChain.get_name() + "_geom");
  for (const auto& it : mGeomList) {
    for (const auto& i : it.second) {
      if (i.mGeom) {
        result_geom->addAnchor(i.mGeom->getAnchor());
        result_geom->addElement(i.mColor, i.mGeom);
      }
    }
  }

  // then, construct the proxy_query_model_3D from the resulting (attached) geometries:
  auto result_prox = std::make_shared<geom::proxy_query_model_3D>(
      aKTEChain.get_name() + "_proxy");
  for (const auto& it : mProxyShapeList) {
    for (const auto& i : it.second) {
      if (i) {
        result_prox->addShape(i);
      }
    }
  }

  return {result_geom, result_prox};
}

}  // namespace ReaK::kte
