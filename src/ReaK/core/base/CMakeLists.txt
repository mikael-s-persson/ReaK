
#set(BASE_SOURCES "")

set(BASE_HEADERS 
  "${RKBASEDIR}/defs.hpp"
  "${RKBASEDIR}/exec_time_profiler.hpp"
  "${RKBASEDIR}/named_object.hpp"
  "${RKBASEDIR}/serializable.hpp"
  "${RKBASEDIR}/shared_object.hpp"
  "${RKBASEDIR}/shared_object_base.hpp"
  "${RKBASEDIR}/scope_guard.hpp"
  "${RKBASEDIR}/expected.hpp"
  "${RKBASEDIR}/thread_incl.hpp"
)


#add_library(reak_base STATIC ${BASE_SOURCES})
#setup_custom_target(reak_base "${SRCROOT}${RKBASEDIR}")
setup_headers("${BASE_HEADERS}" "${RKBASEDIR}")

add_executable(test_base "${SRCROOT}${RKBASEDIR}/test.cpp")
setup_custom_target(test_base "${SRCROOT}${RKBASEDIR}")
target_link_libraries(test_base ${EXTRA_SYSTEM_LIBS})

include_directories(AFTER "${SRCROOT}${RKCOREDIR}")




