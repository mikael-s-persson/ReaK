
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
  "${RKBASEDIR}/typed_object.hpp"
)


#add_library(reak_base STATIC ${BASE_SOURCES})
#setup_custom_target(reak_base)
setup_headers("${BASE_HEADERS}" "${RKBASEDIR}")

add_executable(test_base "${SRCROOT}${RKBASEDIR}/test.cpp")
setup_custom_target(test_base)
target_link_libraries(test_base reak_core)

add_executable(unit_test_cnst_string "${SRCROOT}${RKBASEDIR}/unit_test_cnst_string.cpp")
setup_custom_test_program(unit_test_cnst_string)
target_link_libraries(unit_test_cnst_string reak_core)

add_executable(unit_test_endian_conversions "${SRCROOT}${RKBASEDIR}/unit_test_endian_conversions.cpp")
setup_custom_test_program(unit_test_endian_conversions)
target_link_libraries(unit_test_endian_conversions reak_core)


