
set(RECORDERS_SOURCES 
  "${SRCROOT}${RKRECORDERSDIR}/data_record.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/ssv_recorder.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/tsv_recorder.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/bin_recorder.cpp"
)

set(RECORDERS_HEADERS 
  "${RKRECORDERSDIR}/data_record.hpp"
  "${RKRECORDERSDIR}/ssv_recorder.hpp"
  "${RKRECORDERSDIR}/tsv_recorder.hpp"
  "${RKRECORDERSDIR}/bin_recorder.hpp"
)



add_library(reak_recorders STATIC ${RECORDERS_SOURCES})
setup_custom_target(reak_recorders "${SRCROOT}${RKRECORDERSDIR}")
target_link_libraries(reak_recorders ${Boost_LIBRARIES})


setup_headers("${RECORDERS_HEADERS}" "${RKRECORDERSDIR}")

add_executable(test_rec "${SRCROOT}${RKRECORDERSDIR}/test_rec.cpp")
setup_custom_target(test_rec "${SRCROOT}${RKRECORDERSDIR}")

include_directories(BEFORE ${BOOST_INCLUDE_DIRS})
include_directories(AFTER "${SRCROOT}${RKCOREDIR}")

target_link_libraries(test_rec ${Boost_LIBRARIES})


target_link_libraries(test_rec reak_recorders reak_serialization reak_rtti)




