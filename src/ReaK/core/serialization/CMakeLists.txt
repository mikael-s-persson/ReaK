
set(SERIALIZATION_SOURCES 
  "${SRCROOT}${RKSERIALIZATIONDIR}/xml_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/bin_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/protobuf_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/scheme_builder.cpp"
)

set(SERIALIZATION_HEADERS 
  "${RKSERIALIZATIONDIR}/archiver.hpp"
  "${RKSERIALIZATIONDIR}/xml_archiver.hpp"
  "${RKSERIALIZATIONDIR}/bin_archiver.hpp"
  "${RKSERIALIZATIONDIR}/protobuf_archiver.hpp"
  "${RKSERIALIZATIONDIR}/scheme_builder.hpp"
  "${RKSERIALIZATIONDIR}/type_schemes.hpp"
)

add_library(reak_serialization STATIC ${SERIALIZATION_SOURCES})
setup_custom_target(reak_serialization "${SRCROOT}${RKSERIALIZATIONDIR}")
target_link_libraries(reak_serialization ${EXTRA_SYSTEM_LIBS})

setup_headers("${SERIALIZATION_HEADERS}" "${RKSERIALIZATIONDIR}")

#add_executable(test_math "${SRCROOT}${RKMATHDIR}/test.cpp")
#setup_custom_target(test_math "${SRCROOT}${RKMATHDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

#target_link_libraries(test_math ${Boost_LIBRARIES} ${EXTRA_SYSTEM_LIBS})




