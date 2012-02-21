
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

add_library(reak_rtti STATIC ${RTTI_SOURCES})
setup_custom_target(reak_rtti "${SRCROOT}${RKRTTIDIR}")

setup_headers("${RTTI_HEADERS}" "${RKRTTIDIR}")

#add_executable(test_math "${SRCROOT}${RKMATHDIR}/test.cpp")
#setup_custom_target(test_math "${SRCROOT}${RKMATHDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")
#target_link_libraries(test_math ${Boost_LIBRARIES})




