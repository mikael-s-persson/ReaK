
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

#include <ReaK/core/serialization/scheme_builder.hpp>

#include <ReaK/core/base/shared_object.hpp>
#include <ReaK/core/rtti/rtti.hpp>

#include <ReaK/core/serialization/type_schemes.hpp>

#include <string>

namespace ReaK::serialization {

std::map<std::string, std::shared_ptr<type_scheme>>& get_global_schemes() {
  static std::map<std::string, std::shared_ptr<type_scheme>> gbl_inst;
  return gbl_inst;
};

scheme_builder::scheme_builder() {
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("RootScheme", nullptr, 0));
}

scheme_builder::~scheme_builder() {
  while (!field_stack.empty()) {
    field_stack
        .pop();  // normally, there should be only one object here at the destruction.
  }
}

oarchive& scheme_builder::saveToNewArchive_impl(
    const serializable_shared_pointer& Item, const std::string& FileName) {
  return scheme_builder::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("item", Item));
}

oarchive& scheme_builder::saveToNewArchiveNamed_impl(
    const std::pair<std::string, const serializable_shared_pointer&>& Item,
    const std::string& FileName) {
  return scheme_builder::save_serializable_ptr(Item);
}

oarchive& scheme_builder::save_serializable_ptr(
    const serializable_shared_pointer& Item) {
  return scheme_builder::save_serializable_ptr(
      std::pair<std::string, const serializable_shared_pointer&>("item_ptr",
                                                                 Item));
}

oarchive& scheme_builder::save_serializable_ptr(
    const std::pair<std::string, const serializable_shared_pointer&>& Item) {
  if (!Item.second) {
    return *this;
  }

  auto it = mObjRegMap.find(Item.second);

  if (it != mObjRegMap.end()) {
    return *this;  // this object was already processed, in which case we just don't do anything.
  }

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  rtti::so_type* so_type_ptr = Item.second->getObjectType();

  // See if we need to add the std::shared_ptr< Type > type-scheme for this object type.
  constexpr auto tname =
      rtti::get_type_id<serializable_shared_pointer>::type_name;
  std::string ptr_type_name = std::string(tname);
  ptr_type_name += "<" + so_type_ptr->TypeName() + ">";
  auto itm = scheme_map.find(ptr_type_name);
  if (itm == scheme_map.end()) {
    constexpr auto iname = rtti::get_type_id<unsigned int>::type_name;
    auto objID_itm = scheme_map.find(std::string(iname));
    std::shared_ptr<type_scheme> objID_sch;
    if (objID_itm == scheme_map.end()) {
      objID_sch = std::make_shared<primitive_scheme<std::size_t>>();
      scheme_map[std::string(iname)] = objID_sch;
    } else {
      objID_sch = objID_itm->second;
    }

    auto ser_sch_ptr = std::make_shared<serializable_ptr_scheme>(
        so_type_ptr->TypeName(), so_type_ptr->TypeID_begin(),
        so_type_ptr->TypeVersion(), objID_sch);
    scheme_map[ptr_type_name] = ser_sch_ptr;
  }

  // Create a dummy field stack entry, and then save the pointee object (to make sure the type-scheme for the pointee
  // type is registered).
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));

  save_serializable(
      std::pair<std::string, const serializable&>("Dummy", *(Item.second)));

  field_stack.pop();

  return *this;
}

oarchive& scheme_builder::save_serializable(const serializable& Item) {
  return scheme_builder::save_serializable(
      std::pair<std::string, const serializable&>("item", Item));
}

oarchive& scheme_builder::save_serializable(
    const std::pair<std::string, const serializable&>& Item) {

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  rtti::so_type* so_type_ptr = Item.second.getObjectType();

  auto itm = scheme_map.find(so_type_ptr->TypeName());
  std::shared_ptr<type_scheme> sch_ptr;
  if (itm == scheme_map.end()) {
    auto ser_sch_ptr = std::make_shared<serializable_obj_scheme>(
        so_type_ptr->TypeName(), so_type_ptr->TypeID_begin(),
        so_type_ptr->TypeVersion());
    scheme_map[so_type_ptr->TypeName()] = ser_sch_ptr;
    sch_ptr = ser_sch_ptr;

    field_stack.push(ser_sch_ptr);

    Item.second.save(*this, so_type_ptr->TypeVersion());

    field_stack.pop();
  } else {
    sch_ptr = itm->second;
  }

  field_stack.top()->add_field(Item.first, sch_ptr);
  return *this;
}

template <typename T>
void scheme_builder::save_primitive(const std::string& aName) {

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  constexpr auto tname = rtti::get_type_id<T>::type_name;
  auto itm = scheme_map.find(std::string(tname));
  std::shared_ptr<type_scheme> sch_ptr;
  if (itm == scheme_map.end()) {
    sch_ptr = std::make_shared<primitive_scheme<T>>();
    scheme_map[std::string(tname)] = sch_ptr;
  } else {
    sch_ptr = itm->second;
  }

  field_stack.top()->add_field(aName, sch_ptr);
}

oarchive& scheme_builder::save_char(char i) {
  return scheme_builder::save_char(std::pair<std::string, char>("i", i));
}

oarchive& scheme_builder::save_char(const std::pair<std::string, char>& i) {
  save_primitive<char>(i.first);
  return *this;
}

oarchive& scheme_builder::save_unsigned_char(unsigned char u) {
  return scheme_builder::save_unsigned_char(
      std::pair<std::string, unsigned char>("u", u));
}

oarchive& scheme_builder::save_unsigned_char(
    const std::pair<std::string, unsigned char>& u) {
  save_primitive<unsigned char>(u.first);
  return *this;
}

oarchive& scheme_builder::save_int(std::ptrdiff_t i) {
  return scheme_builder::save_int(
      std::pair<std::string, std::ptrdiff_t>("i", i));
}

oarchive& scheme_builder::save_int(
    const std::pair<std::string, std::ptrdiff_t>& i) {
  save_primitive<std::ptrdiff_t>(i.first);
  return *this;
}

oarchive& scheme_builder::save_unsigned_int(std::size_t u) {
  return scheme_builder::save_unsigned_int(
      std::pair<std::string, std::size_t>("u", u));
}

oarchive& scheme_builder::save_unsigned_int(
    const std::pair<std::string, std::size_t>& u) {
  save_primitive<std::size_t>(u.first);
  return *this;
}

oarchive& scheme_builder::save_float(float f) {
  return scheme_builder::save_float(std::pair<std::string, float>("f", f));
}

oarchive& scheme_builder::save_float(const std::pair<std::string, float>& f) {
  save_primitive<float>(f.first);
  return *this;
}

oarchive& scheme_builder::save_double(double d) {
  return scheme_builder::save_double(std::pair<std::string, double>("d", d));
}

oarchive& scheme_builder::save_double(const std::pair<std::string, double>& d) {
  save_primitive<double>(d.first);
  return *this;
}

oarchive& scheme_builder::save_bool(bool b) {
  return scheme_builder::save_bool(std::pair<std::string, bool>("b", b));
}

oarchive& scheme_builder::save_bool(const std::pair<std::string, bool>& b) {
  save_primitive<bool>(b.first);
  return *this;
}

oarchive& scheme_builder::save_string(const std::string& s) {
  return scheme_builder::save_string(
      std::pair<std::string, const std::string&>("str", s));
}

oarchive& scheme_builder::save_string(
    const std::pair<std::string, const std::string&>& s) {
  save_primitive<std::string>(s.first);
  return *this;
}

void scheme_builder::signal_polymorphic_field(const std::string& aBaseTypeName,
                                              const unsigned int* aTypeID,
                                              const std::string& aFieldName) {

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  // check if the scheme already exists
  constexpr auto tname =
      rtti::get_type_id<serializable_shared_pointer>::type_name;
  std::string ptr_type_name = std::string(tname);
  ptr_type_name += "<" + aBaseTypeName + ">";
  auto itm = scheme_map.find(ptr_type_name);
  std::shared_ptr<type_scheme> sch_ptr;
  if (itm == scheme_map.end()) {
    constexpr auto iname = rtti::get_type_id<unsigned int>::type_name;
    auto objID_itm = scheme_map.find(std::string(iname));
    std::shared_ptr<type_scheme> objID_sch;
    if (objID_itm == scheme_map.end()) {
      objID_sch = std::make_shared<primitive_scheme<std::size_t>>();
      scheme_map[std::string(iname)] = objID_sch;
    } else {
      objID_sch = objID_itm->second;
    }

    rtti::so_type* so_type_sptr =
        rtti::getRKSharedObjTypeRepo().findType(aTypeID);
    if (so_type_sptr != nullptr) {
      auto ser_sch_ptr = std::make_shared<serializable_ptr_scheme>(
          so_type_sptr->TypeName(), so_type_sptr->TypeID_begin(),
          so_type_sptr->TypeVersion(), objID_sch);
      scheme_map[ptr_type_name] = ser_sch_ptr;
      sch_ptr = ser_sch_ptr;

      // create a dummy type-scheme to accumulate the "fields" of all the descendants.
      field_stack.push(
          std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));

      for (unsigned int i = 0; i < so_type_sptr->getDirectDescendantCount();
           ++i) {
        rtti::so_type* tmp = so_type_sptr->getDirectDescendant(i);
        signal_polymorphic_field(tmp->TypeName(), tmp->TypeID_begin(), "Dummy");
        save_serializable_ptr(
            tmp->CreateObject());  // if creation fails, the std::shared_ptr saver will not save it (abstract class).
      }

      field_stack.pop();  // discard the dummy type-scheme.
    }
  } else {
    sch_ptr = itm->second;
  }

  // add the field to the current list (top of stack).
  field_stack.top()->add_field(aFieldName, sch_ptr);
}

void scheme_builder::start_repeated_field(const std::string& aTypeName) {
  field_stack.top()->pop_last_field();  // pop the count field.
  value_name_stack.push(
      std::pair<std::string, std::string>("values", aTypeName));
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));
}

void scheme_builder::start_repeated_field(const std::string& aTypeName,
                                          const std::string& aName) {
  field_stack.top()->pop_last_field();  // pop the count field.
  value_name_stack.push(std::pair<std::string, std::string>(aName, aTypeName));
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));
}

void scheme_builder::finish_repeated_field() {

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  // look for the current top of the value_name_stack to find the recorded type scheme.
  auto itm = scheme_map.find(value_name_stack.top().second);
  std::shared_ptr<type_scheme> val_sch_ptr;
  if (itm != scheme_map.end()) {
    val_sch_ptr = itm->second;
  }

  field_stack.pop();  // pop the dummy type scheme.

  auto sch_ptr = std::make_shared<vector_type_scheme>(
      "std::vector<" + value_name_stack.top().second + ">", val_sch_ptr);

  field_stack.top()->add_field(value_name_stack.top().first, sch_ptr);
  scheme_map[sch_ptr->get_type_name()] = sch_ptr;

  value_name_stack.pop();
}

void scheme_builder::start_repeated_pair(const std::string& aTypeName1,
                                         const std::string& aTypeName2) {
  field_stack.top()->pop_last_field();  // pop the count field.
  value_name_stack.push(std::pair<std::string, std::string>("map", aTypeName1));
  value_name_stack.push(std::pair<std::string, std::string>("map", aTypeName2));
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));
}

void scheme_builder::start_repeated_pair(const std::string& aTypeName1,
                                         const std::string& aTypeName2,
                                         const std::string& aName) {
  field_stack.top()->pop_last_field();  // pop the count field.
  value_name_stack.push(std::pair<std::string, std::string>(aName, aTypeName1));
  value_name_stack.push(std::pair<std::string, std::string>(aName, aTypeName2));
  field_stack.push(
      std::make_shared<serializable_obj_scheme>("DummyType", nullptr, 0));
}

void scheme_builder::finish_repeated_pair() {

  std::map<std::string, std::shared_ptr<type_scheme>>& scheme_map =
      get_global_schemes();

  // look for the current top of the value_name_stack to find the recorded type scheme.
  auto val_itm = scheme_map.find(value_name_stack.top().second);
  std::shared_ptr<type_scheme> val_sch_ptr;
  if (val_itm != scheme_map.end()) {
    val_sch_ptr = val_itm->second;
  }

  std::string value_type_name = value_name_stack.top().second;
  value_name_stack.pop();

  auto key_itm = scheme_map.find(value_name_stack.top().second);
  std::shared_ptr<type_scheme> key_sch_ptr;
  if (key_itm != scheme_map.end()) {
    key_sch_ptr = key_itm->second;
  }

  field_stack.pop();  // pop the dummy type scheme.

  auto sch_ptr = std::make_shared<map_type_scheme>(
      "std::map<" + value_name_stack.top().second + "," + value_type_name + ">",
      key_sch_ptr, val_sch_ptr);

  field_stack.top()->add_field(value_name_stack.top().first, sch_ptr);
  scheme_map[sch_ptr->get_type_name()] = sch_ptr;

  value_name_stack.pop();
}

}  // namespace ReaK::serialization
