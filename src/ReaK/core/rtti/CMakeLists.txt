
set(RTTI_SOURCES 
  "${SRCROOT}${RKRTTIDIR}/so_type.cpp"
  "${SRCROOT}${RKRTTIDIR}/so_type_repo.cpp"
)

set(RTTI_HEADERS 
  "${RKRTTIDIR}/so_register_type.hpp"
  "${RKRTTIDIR}/so_type.hpp"
  "${RKRTTIDIR}/so_type_repo.hpp"
  "${RKRTTIDIR}/typed_containers.hpp"
  "${RKRTTIDIR}/typed_object.hpp"
  "${RKRTTIDIR}/typed_primitives.hpp"
)

add_library(reakobj_rtti OBJECT ${RTTI_SOURCES})

setup_headers("${RTTI_HEADERS}" "${RKRTTIDIR}")




