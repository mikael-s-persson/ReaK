
#set(BASE_SOURCES "")

set(BASE_HEADERS 
  "${RKBASEDIR}/backtrace_exception.hpp"
  "${RKBASEDIR}/chrono_incl.hpp"
  "${RKBASEDIR}/defs.hpp"
  "${RKBASEDIR}/exec_time_profiler.hpp"
  "${RKBASEDIR}/expected.hpp"
  "${RKBASEDIR}/function_incl.hpp"
  "${RKBASEDIR}/global_rng.hpp"
  "${RKBASEDIR}/misc_math.hpp"
  "${RKBASEDIR}/named_object.hpp"
  "${RKBASEDIR}/py_fixes.hpp"
  "${RKBASEDIR}/scope_guard.hpp"
  "${RKBASEDIR}/serializable.hpp"
  "${RKBASEDIR}/shared_object.hpp"
  "${RKBASEDIR}/shared_object_base.hpp"
  "${RKBASEDIR}/thread_incl.hpp"
)


#add_library(reak_base STATIC ${BASE_SOURCES})
#setup_custom_target(reak_base "${SRCROOT}${RKBASEDIR}")
setup_headers("${BASE_HEADERS}" "${RKBASEDIR}")

add_executable(test_base "${SRCROOT}${RKBASEDIR}/test.cpp")
setup_custom_target(test_base "${SRCROOT}${RKBASEDIR}")
target_link_libraries(test_base reak_core)


