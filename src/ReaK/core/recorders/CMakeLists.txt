
set(RECORDERS_SOURCES 
  "${SRCROOT}${RKRECORDERSDIR}/data_record.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/ssv_recorder.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/tsv_recorder.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/bin_recorder.cpp"
  "${SRCROOT}${RKRECORDERSDIR}/tcp_recorder.cpp"
)

set(RECORDERS_HEADERS 
  "${RKRECORDERSDIR}/data_record.hpp"
  "${RKRECORDERSDIR}/ssv_recorder.hpp"
  "${RKRECORDERSDIR}/tsv_recorder.hpp"
  "${RKRECORDERSDIR}/bin_recorder.hpp"
  "${RKRECORDERSDIR}/tcp_recorder.hpp"
)



add_library(reakobj_recorders OBJECT ${RECORDERS_SOURCES})

setup_headers("${RECORDERS_HEADERS}" "${RKRECORDERSDIR}")

add_executable(unit_test_recorders "${SRCROOT}${RKRECORDERSDIR}/unit_test_recorders.cpp")
setup_custom_test_program(unit_test_recorders "${SRCROOT}${RKRECORDERSDIR}")
target_link_libraries(unit_test_recorders reak_core)






