
set(SERIALIZATION_SOURCES 
  "${SRCROOT}${RKSERIALIZATIONDIR}/xml_archiver.cpp"
  "${SRCROOT}${RKSERIALIZATIONDIR}/bin_archiver.cpp"
)

set(SERIALIZATION_HEADERS 
  "${RKSERIALIZATIONDIR}/archiver.hpp"
  "${RKSERIALIZATIONDIR}/xml_archiver.hpp"
  "${RKSERIALIZATIONDIR}/bin_archiver.hpp"
)

add_library(reak_serialization STATIC ${SERIALIZATION_SOURCES})
setup_custom_target(reak_serialization "${SRCROOT}${RKSERIALIZATIONDIR}")

setup_headers("${SERIALIZATION_HEADERS}" "${RKSERIALIZATIONDIR}")

#add_executable(test_math "${SRCROOT}${RKMATHDIR}/test.cpp")
#setup_custom_target(test_math "${SRCROOT}${RKMATHDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

#target_link_libraries(test_math ${BOOST_THREAD_LIB})




