
set(SERIALIZATION_SOURCES 
  "${SRCROOT}${RKSERIALIZATIONDIR}/xml_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/bin_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/protobuf_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/objtree_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/scheme_builder.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/archiver_factory.cpp"
)

set(SERIALIZATION_HEADERS 
  "${RKSERIALIZATIONDIR}/archiver.hpp"
  "${RKSERIALIZATIONDIR}/xml_archiver.hpp"
  "${RKSERIALIZATIONDIR}/bin_archiver.hpp"
  "${RKSERIALIZATIONDIR}/protobuf_archiver.hpp"
  "${RKSERIALIZATIONDIR}/objtree_archiver.hpp"
  "${RKSERIALIZATIONDIR}/scheme_builder.hpp"
  "${RKSERIALIZATIONDIR}/type_schemes.hpp"
  "${RKSERIALIZATIONDIR}/archiver_factory.hpp"
)

add_library(reakobj_serialization OBJECT ${SERIALIZATION_SOURCES})

setup_headers("${SERIALIZATION_HEADERS}" "${RKSERIALIZATIONDIR}")

add_executable(unit_test_serialization "${SRCROOT}${RKSERIALIZATIONDIR}/unit_test_serialization.cpp")
setup_custom_test_program(unit_test_serialization "${SRCROOT}${RKSERIALIZATIONDIR}")
target_link_libraries(unit_test_serialization reak_core)




